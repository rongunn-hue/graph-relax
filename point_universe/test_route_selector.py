"""Source-only tests for deterministic route candidate selection."""

from __future__ import annotations

import pytest

from .route_selector import RouteSelector
from .router import PointRouter, RouteCandidate
from .universe import HypotheticalWhiteConsequence, Point, PointState, PointUniverse, Square


def obstacle_universe() -> PointUniverse:
    universe = PointUniverse(7, 7, line_separation=0, square_clearance=0)
    universe.add_square(Square("OBSTACLE", 2, 2, 2, 2))
    return universe


def snapshot(universe: PointUniverse) -> tuple[
    tuple[tuple[PointState, ...], ...],
    tuple[str, ...],
    tuple[str, ...],
    tuple[tuple[Point, tuple[str, ...]], ...],
]:
    matrix = tuple(
        tuple(universe.point_state((x, y)) for x in range(universe.width))
        for y in range(universe.height)
    )
    return (
        matrix,
        tuple(sorted(universe.squares)),
        tuple(sorted(universe.lines)),
        tuple(sorted((point, tuple(sorted(lines))) for point, lines in universe.point_lines.items())),
    )


def consequence(delta_R: int) -> HypotheticalWhiteConsequence:
    return HypotheticalWhiteConsequence(
        R_before=1,
        R_after=1 + delta_R,
        delta_R=delta_R,
        removed_white_points=0,
        resulting_region_sizes=(),
    )


def candidate(points: tuple[Point, ...], delta_R: int) -> RouteCandidate:
    return RouteCandidate(points=points, consequence=consequence(delta_R))


def test_selector_prefers_lower_delta_r_over_shorter_path_in_established_router_case() -> None:
    universe = obstacle_universe()
    before = snapshot(universe)
    candidates = PointRouter(universe).candidate_paths((0, 3), (6, 3), max_paths=2)
    selector = RouteSelector()

    assert len(candidates) == 2
    assert len(candidates[0].points) == 31
    assert candidates[0].consequence.delta_R == 0
    assert len(candidates[1].points) == 29
    assert candidates[1].consequence.delta_R == 1

    selected = selector.select(candidates)
    assert selected is candidates[0]
    assert selector.select(tuple(reversed(candidates))) is candidates[0]
    assert snapshot(universe) == before


def test_selector_uses_fewer_points_to_break_delta_r_tie() -> None:
    longer = candidate(((0, 0), (1, 0), (2, 0)), delta_R=0)
    shorter = candidate(((0, 0), (1, 0)), delta_R=0)
    selector = RouteSelector()

    assert selector.select((longer, shorter)) is shorter
    assert selector.select((shorter, longer)) is shorter


def test_selector_uses_point_tuple_to_break_length_tie() -> None:
    lexicographically_larger = candidate(((0, 0), (1, 0)), delta_R=0)
    lexicographically_smaller = candidate(((0, 0), (0, 1)), delta_R=0)
    selector = RouteSelector()

    assert selector.select((lexicographically_larger, lexicographically_smaller)) is lexicographically_smaller
    assert selector.select((lexicographically_smaller, lexicographically_larger)) is lexicographically_smaller


def test_selector_rejects_empty_input() -> None:
    with pytest.raises(ValueError, match="empty"):
        RouteSelector().select(())
