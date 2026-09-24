"""Source-only tests for the first point-universe router."""

from __future__ import annotations

import pytest

from .router import PointRouter
from .universe import Line, Point, PointState, PointUniverse, Square, WhiteSpaceViewer


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


def assert_legal_simple_white_path(universe: PointUniverse, path: tuple[Point, ...], start: Point, goal: Point) -> None:
    assert path[0] == start
    assert path[-1] == goal
    assert len(set(path)) == len(path)
    for point in path:
        assert universe.point_state(point) & PointState.WHITE
    for a, b in zip(path, path[1:]):
        assert abs(a[0] - b[0]) + abs(a[1] - b[1]) == 1


def test_router_generates_multiple_legal_deterministic_candidates_without_mutation() -> None:
    universe = obstacle_universe()
    before = snapshot(universe)
    start = (0, 3)
    goal = (6, 3)
    assert universe.point_state(start) & PointState.WHITE
    assert universe.point_state(goal) & PointState.WHITE
    router = PointRouter(universe)

    first = router.candidate_paths(start, goal, max_paths=2)
    second = router.candidate_paths(start, goal, max_paths=2)

    assert len(first) == 2
    assert len(second) == 2
    assert [candidate.points for candidate in first] == [candidate.points for candidate in second]
    assert first[0].points != first[1].points

    for candidate in first:
        assert_legal_simple_white_path(universe, candidate.points, start, goal)
        expected = WhiteSpaceViewer(universe).hypothetical_white_consequence(candidate.points)
        assert candidate.consequence == expected

    one_candidate = router.candidate_paths(start, goal, max_paths=1)
    assert len(one_candidate) == 1
    assert one_candidate[0].points == first[0].points

    assert snapshot(universe) == before


def test_blocked_connection_returns_no_candidate() -> None:
    universe = PointUniverse(3, 3, line_separation=0, square_clearance=0)
    universe.add_line(Line("SEPARATOR", "N", ((1, 0), (1, 2))))
    router = PointRouter(universe)
    assert router.candidate_paths((0, 1), (2, 1), max_paths=2) == []


def test_start_equals_goal_returns_single_point_candidate() -> None:
    universe = obstacle_universe()
    point = (0, 0)
    router = PointRouter(universe)
    candidates = router.candidate_paths(point, point, max_paths=3)
    assert len(candidates) == 1
    assert candidates[0].points == (point,)
    assert candidates[0].consequence == WhiteSpaceViewer(universe).hypothetical_white_consequence((point,))


def test_router_rejects_invalid_inputs() -> None:
    universe = obstacle_universe()
    router = PointRouter(universe)

    with pytest.raises(ValueError, match="max_paths"):
        router.candidate_paths((0, 0), (1, 0), max_paths=0)

    with pytest.raises(ValueError, match="max_paths"):
        router.candidate_paths((0, 0), (1, 0), max_paths=-1)

    with pytest.raises(ValueError, match="start point outside"):
        router.candidate_paths((-1, 0), (1, 0), max_paths=1)

    with pytest.raises(ValueError, match="goal point outside"):
        router.candidate_paths((0, 0), (7, 0), max_paths=1)

    with pytest.raises(ValueError, match="start point is not WHITE"):
        router.candidate_paths((2, 2), (0, 0), max_paths=1)

    with pytest.raises(ValueError, match="goal point is not WHITE"):
        router.candidate_paths((0, 0), (2, 2), max_paths=1)
