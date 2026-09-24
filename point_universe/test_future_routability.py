"""Source-only tests for future-routability perception."""

from __future__ import annotations

from .future_routability import future_routability_after_candidate
from .router import PointRouter, RouteCandidate
from .universe import Line, Point, PointState, PointUniverse


def future_universe() -> PointUniverse:
    universe = PointUniverse(7, 7, line_separation=0, square_clearance=0)
    universe.add_line(Line("WALL_TOP", "WALL", ((0, 0), (6, 0))))
    universe.add_line(Line("WALL_BOTTOM", "WALL", ((0, 6), (6, 6))))
    universe.add_line(Line("WALL_LEFT_Y1", "WALL", ((0, 1), (2, 1))))
    universe.add_line(Line("WALL_LEFT_Y2", "WALL", ((0, 2), (2, 2))))
    universe.add_line(Line("WALL_LEFT_Y4", "WALL", ((0, 4), (2, 4))))
    universe.add_line(Line("WALL_LEFT_Y5", "WALL", ((0, 5), (2, 5))))
    universe.add_line(Line("WALL_RIGHT", "WALL", ((6, 1), (6, 5))))
    return universe


def authoritative_snapshot(universe: PointUniverse) -> tuple[
    tuple[tuple[PointState, ...], ...],
    tuple[str, ...],
    tuple[tuple[str, tuple[Point, ...]], ...],
    tuple[tuple[Point, tuple[str, ...]], ...],
    tuple[tuple[str, tuple[Point, ...]], ...],
]:
    matrix = tuple(
        tuple(universe.point_state((x, y)) for x in range(universe.width))
        for y in range(universe.height)
    )
    return (
        matrix,
        tuple(sorted(universe.lines)),
        tuple(sorted((line_id, line.points) for line_id, line in universe.lines.items())),
        tuple(sorted((point, tuple(sorted(line_ids))) for point, line_ids in universe.point_lines.items())),
        tuple(sorted((line_id, tuple(sorted(points))) for line_id, points in universe.line_points.items())),
    )


def assert_candidate_legal_in_disposable_copy(candidate: RouteCandidate) -> None:
    disposable = future_universe()
    disposable.add_line(Line("A", "NET_A", candidate.points))
    assert "A" in disposable.lines
    assert disposable.lines["A"].points == candidate.points


def disposable_b_count_after(candidate: RouteCandidate, b_start: Point, b_goal: Point, b_max_paths: int) -> int:
    disposable = future_universe()
    disposable.add_line(Line("A", "NET_A", candidate.points))
    assert disposable.point_state(b_start) & PointState.WHITE
    assert disposable.point_state(b_goal) & PointState.WHITE
    return len(PointRouter(disposable).candidate_paths(b_start, b_goal, max_paths=b_max_paths))


def classify_generated_candidates(
    universe: PointUniverse,
    candidates: tuple[RouteCandidate, ...],
    b_start: Point,
    b_goal: Point,
    b_max_paths: int,
) -> tuple[RouteCandidate, object, RouteCandidate, object]:
    blocking: RouteCandidate | None = None
    blocking_result: object | None = None
    preserving: RouteCandidate | None = None
    preserving_result: object | None = None
    for candidate in candidates:
        if b_start in candidate.points or b_goal in candidate.points:
            continue
        result = future_routability_after_candidate(
            universe,
            candidate,
            B_start=b_start,
            B_goal=b_goal,
            B_max_paths=b_max_paths,
        )
        if not result.B_routable and result.B_candidate_count == 0 and blocking is None:
            blocking = candidate
            blocking_result = result
        if result.B_routable and result.B_candidate_count >= 1 and preserving is None:
            preserving = candidate
            preserving_result = result
        if blocking is not None and preserving is not None:
            return blocking, blocking_result, preserving, preserving_result
    raise AssertionError("generated A candidates did not include both blocking and preserving examples")


def test_future_routability_distinguishes_blocking_and_preserving_a_candidates_without_mutation() -> None:
    universe = future_universe()
    snapshot = authoritative_snapshot(universe)
    a_start = (3, 1)
    a_goal = (3, 5)
    a_max_paths = 16
    b_start = (1, 3)
    b_goal = (4, 3)
    b_max_paths = 2
    a_candidates = tuple(PointRouter(universe).candidate_paths(a_start, a_goal, max_paths=a_max_paths))
    all_candidates_before = tuple((candidate.points, candidate.consequence) for candidate in a_candidates)
    blocking, blocking_result, preserving, preserving_result = classify_generated_candidates(
        universe,
        a_candidates,
        b_start,
        b_goal,
        b_max_paths,
    )
    blocking_points_before = blocking.points
    preserving_points_before = preserving.points
    blocking_consequence_before = blocking.consequence
    preserving_consequence_before = preserving.consequence

    assert len(a_candidates) >= 2
    assert blocking in a_candidates
    assert preserving in a_candidates
    assert b_start not in blocking.points
    assert b_goal not in blocking.points
    assert b_start not in preserving.points
    assert b_goal not in preserving.points
    assert universe.point_state(b_start) & PointState.WHITE
    assert universe.point_state(b_goal) & PointState.WHITE
    assert_candidate_legal_in_disposable_copy(blocking)
    assert_candidate_legal_in_disposable_copy(preserving)
    assert blocking.points != preserving.points
    assert authoritative_snapshot(universe) == snapshot

    assert not blocking_result.B_routable
    assert blocking_result.B_candidate_count == 0
    assert blocking_result.B_candidates == ()
    assert preserving_result.B_routable
    assert preserving_result.B_candidate_count >= 1
    assert blocking_result.B_candidate_count != preserving_result.B_candidate_count

    assert blocking_result.B_candidate_count == disposable_b_count_after(blocking, b_start, b_goal, b_max_paths)
    assert preserving_result.B_candidate_count == disposable_b_count_after(preserving, b_start, b_goal, b_max_paths)

    assert blocking.points == blocking_points_before
    assert preserving.points == preserving_points_before
    assert blocking.consequence == blocking_consequence_before
    assert preserving.consequence == preserving_consequence_before
    assert tuple((candidate.points, candidate.consequence) for candidate in a_candidates) == all_candidates_before
    assert "HYPOTHETICAL_A" not in universe.lines
    assert "A" not in universe.lines
