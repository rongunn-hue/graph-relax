"""Controlled comparison of existing RouteSelector vs capacity-aware selection."""

from __future__ import annotations

from dataclasses import dataclass
from itertools import combinations, permutations
from typing import Callable, Iterable, Sequence

from .capacity_aware_selector import CapacityAwareSelector, fixed_capacity_aware_key
from .route_committer import RouteCommitter
from .route_selector import RouteSelector
from .router import PointRouter, RouteCandidate
from .universe import Line, Point, PointState, PointUniverse, TerminalPair, WhiteSpaceViewer


ConnectionName = str


@dataclass(frozen=True)
class ConnectionSpec:
    name: ConnectionName
    start: Point
    goal: Point


@dataclass(frozen=True)
class FixtureSpec:
    width: int
    height: int
    obstacle_lines: tuple[tuple[str, tuple[Point, ...]], ...]
    connections: tuple[ConnectionSpec, ConnectionSpec, ConnectionSpec]
    max_paths: int
    examined_index: int


@dataclass(frozen=True)
class FrozenFixture:
    spec: FixtureSpec
    a_candidates: tuple[tuple[Point, ...], ...]
    baseline_a_points: tuple[Point, ...]
    capacity_a_points: tuple[Point, ...]
    baseline_key: tuple[object, ...]
    capacity_key: tuple[object, ...]


@dataclass(frozen=True)
class SequentialRunResult:
    committed_count: int
    selected: tuple[tuple[str, tuple[Point, ...] | None], ...]
    candidate_counts: tuple[tuple[str, int], ...]
    final_white_region_count: int


def build_universe(spec: FixtureSpec) -> PointUniverse:
    universe = PointUniverse(spec.width, spec.height, line_separation=0, square_clearance=0)
    for line_id, points in spec.obstacle_lines:
        universe.add_line(Line(line_id, "OBSTACLE", points))
    return universe


def white_points(universe: PointUniverse) -> tuple[Point, ...]:
    return tuple(sorted(universe.layer(PointState.WHITE), key=lambda p: (p[1], p[0])))


def fixture_catalog() -> Iterable[tuple[int, int, tuple[tuple[str, tuple[Point, ...]], ...]]]:
    yield (5, 3, ())
    yield (7, 3, (("WALL", ((3, 0), (3, 2))),))
    yield (5, 4, ())
    yield (6, 4, (("WALL", ((3, 0), (3, 3))),))


def curated_connection_triples() -> Iterable[tuple[ConnectionSpec, ConnectionSpec, ConnectionSpec]]:
    yield (
        ConnectionSpec("A", (0, 0), (0, 2)),
        ConnectionSpec("B", (2, 0), (2, 2)),
        ConnectionSpec("C", (4, 1), (6, 1)),
    )
    yield (
        ConnectionSpec("A", (0, 0), (2, 0)),
        ConnectionSpec("B", (0, 1), (2, 1)),
        ConnectionSpec("C", (4, 1), (6, 1)),
    )
    yield (
        ConnectionSpec("A", (0, 2), (2, 2)),
        ConnectionSpec("B", (0, 1), (2, 1)),
        ConnectionSpec("C", (4, 1), (6, 1)),
    )
    yield (
        ConnectionSpec("A", (0, 1), (2, 1)),
        ConnectionSpec("B", (2, 0), (2, 2)),
        ConnectionSpec("C", (4, 1), (6, 1)),
    )


def deterministic_fixture_specs(max_paths: int = 6) -> Iterable[FixtureSpec]:
    examined = 0
    for width, height, obstacles in fixture_catalog():
        for connections in curated_connection_triples():
            if all(0 <= point[0] < width and 0 <= point[1] < height for connection in connections for point in (connection.start, connection.goal)):
                yield FixtureSpec(width, height, obstacles, connections, max_paths, examined)
                examined += 1
        probe = FixtureSpec(width, height, obstacles, (), max_paths, examined)
        universe = build_universe(probe)
        points = white_points(universe)
        ordered_pairs = tuple((a, b) for a, b in permutations(points, 2) if a != b)
        emitted_for_fixture = 0
        for a_start, a_goal in ordered_pairs:
            for b_start, b_goal in ordered_pairs:
                if len({a_start, a_goal, b_start, b_goal}) < 4:
                    continue
                for c_start, c_goal in ordered_pairs:
                    if len({a_start, a_goal, b_start, b_goal, c_start, c_goal}) < 6:
                        continue
                    connections = (
                        ConnectionSpec("A", a_start, a_goal),
                        ConnectionSpec("B", b_start, b_goal),
                        ConnectionSpec("C", c_start, c_goal),
                    )
                    yield FixtureSpec(width, height, obstacles, connections, max_paths, examined)
                    examined += 1
                    emitted_for_fixture += 1
                    if emitted_for_fixture >= 250:
                        break
                if emitted_for_fixture >= 250:
                    break
            if emitted_for_fixture >= 250:
                break


def generated_candidates(universe: PointUniverse, connection: ConnectionSpec, max_paths: int) -> tuple[RouteCandidate, ...]:
    try:
        return tuple(PointRouter(universe).candidate_paths(connection.start, connection.goal, max_paths=max_paths))
    except ValueError:
        return ()


def pending_pairs(connections: Sequence[ConnectionSpec]) -> tuple[TerminalPair, ...]:
    return tuple(TerminalPair(connection.name, connection.start, connection.goal) for connection in connections)


def candidate_avoids_future_terminals(candidate: RouteCandidate, future: Sequence[ConnectionSpec]) -> bool:
    future_terminals = {connection.start for connection in future} | {connection.goal for connection in future}
    return not (set(candidate.points) & future_terminals)


def discover_fixture() -> FrozenFixture:
    selector = RouteSelector()
    capacity_selector = CapacityAwareSelector()
    for spec in deterministic_fixture_specs():
        universe = build_universe(spec)
        a, b, c = spec.connections
        a_candidates = generated_candidates(universe, a, spec.max_paths)
        if len(a_candidates) < 2:
            continue
        if not all(candidate_avoids_future_terminals(candidate, (b, c)) for candidate in a_candidates):
            continue
        try:
            baseline = selector.select(a_candidates)
            capacity_selection = capacity_selector.select(
                a_candidates,
                viewer=WhiteSpaceViewer(universe),
                pending_pairs=pending_pairs((b, c)),
            )
        except ValueError:
            continue
        candidate_infos = []
        for candidate in a_candidates:
            try:
                multi = WhiteSpaceViewer(universe).hypothetical_multi_terminal_capacity_consequence(
                    candidate.points,
                    pending_pairs((b, c)),
                )
                aggregate = WhiteSpaceViewer(universe).aggregate_terminal_capacity_consequence(multi)
            except ValueError:
                continue
            candidate_infos.append((candidate, fixed_capacity_aware_key(candidate, aggregate)))
        if len(candidate_infos) < 2:
            continue
        baseline_info = next((item for item in candidate_infos if item[0] is baseline), None)
        if baseline_info is None:
            continue
        best_info = min(candidate_infos, key=lambda item: item[1])
        if best_info[0] is baseline:
            continue
        if not best_info[1] < baseline_info[1]:
            continue
        return FrozenFixture(
            spec=spec,
            a_candidates=tuple(candidate.points for candidate in a_candidates),
            baseline_a_points=baseline.points,
            capacity_a_points=capacity_selection.candidate.points,
            baseline_key=baseline_info[1],
            capacity_key=best_info[1],
        )
    raise AssertionError("no deterministic fixture satisfied the capacity-aware comparison preconditions")


def run_sequential(
    spec: FixtureSpec,
    *,
    select_current: Callable[
        [PointUniverse, tuple[RouteCandidate, ...], tuple[ConnectionSpec, ...]],
        RouteCandidate,
    ],
) -> SequentialRunResult:
    universe = build_universe(spec)
    selected: list[tuple[str, tuple[Point, ...] | None]] = []
    candidate_counts: list[tuple[str, int]] = []
    committed_count = 0
    connections = spec.connections
    for index, connection in enumerate(connections):
        candidates = generated_candidates(universe, connection, spec.max_paths)
        candidate_counts.append((connection.name, len(candidates)))
        if not candidates:
            selected.append((connection.name, None))
            continue
        pending = connections[index + 1 :]
        chosen = select_current(universe, candidates, pending)
        RouteCommitter(universe).commit(chosen, line_id=connection.name, net_id=connection.name)
        selected.append((connection.name, chosen.points))
        committed_count += 1
    return SequentialRunResult(
        committed_count=committed_count,
        selected=tuple(selected),
        candidate_counts=tuple(candidate_counts),
        final_white_region_count=WhiteSpaceViewer(universe).white_region_count(),
    )


def baseline_select(
    _universe: PointUniverse,
    candidates: tuple[RouteCandidate, ...],
    _pending: tuple[ConnectionSpec, ...],
) -> RouteCandidate:
    return RouteSelector().select(candidates)


def capacity_aware_select(
    universe: PointUniverse,
    candidates: tuple[RouteCandidate, ...],
    pending: tuple[ConnectionSpec, ...],
) -> RouteCandidate:
    if not pending:
        return RouteSelector().select(candidates)
    return CapacityAwareSelector().select(
        candidates,
        viewer=WhiteSpaceViewer(universe),
        pending_pairs=pending_pairs(pending),
    ).candidate


def test_capacity_aware_selection_controlled_sequential_comparison() -> None:
    fixture = discover_fixture()
    spec = fixture.spec
    a, b, c = spec.connections

    assert (a.name, b.name, c.name) == ("A", "B", "C")
    assert len(fixture.a_candidates) >= 2
    assert fixture.baseline_a_points in fixture.a_candidates
    assert fixture.capacity_a_points in fixture.a_candidates
    assert fixture.capacity_key < fixture.baseline_key
    assert fixture.baseline_a_points != fixture.capacity_a_points

    baseline_universe = build_universe(spec)
    capacity_universe = build_universe(spec)
    assert [
        candidate.points
        for candidate in generated_candidates(baseline_universe, a, spec.max_paths)
    ] == [
        candidate.points
        for candidate in generated_candidates(capacity_universe, a, spec.max_paths)
    ] == list(fixture.a_candidates)

    baseline = run_sequential(spec, select_current=baseline_select)
    capacity_aware = run_sequential(spec, select_current=capacity_aware_select)

    assert 0 <= baseline.committed_count <= 3
    assert 0 <= capacity_aware.committed_count <= 3
    assert tuple(name for name, _ in baseline.selected) == ("A", "B", "C")
    assert tuple(name for name, _ in capacity_aware.selected) == ("A", "B", "C")
    assert baseline.selected[0][1] == fixture.baseline_a_points
    assert capacity_aware.selected[0][1] == fixture.capacity_a_points
    assert len(baseline.candidate_counts) == 3
    assert len(capacity_aware.candidate_counts) == 3
