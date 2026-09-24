"""Source-only tests for terminal-specific minimum WHITE vertex cuts."""

from __future__ import annotations

from itertools import combinations

import pytest

from .universe import Line, Point, PointState, PointUniverse, WhiteSpaceViewer, four_neighbors


def snapshot(universe: PointUniverse) -> tuple[
    tuple[tuple[PointState, ...], ...],
    tuple[str, ...],
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
        tuple(sorted((point, tuple(sorted(line_ids))) for point, line_ids in universe.point_lines.items())),
        tuple(sorted((line_id, tuple(sorted(points))) for line_id, points in universe.line_points.items())),
    )


def white_points(universe: PointUniverse) -> set[Point]:
    return set(universe.layer(PointState.WHITE))


def connected_after_removing(universe: PointUniverse, start: Point, goal: Point, removed: set[Point]) -> bool:
    remaining = white_points(universe) - set(removed)
    if start not in remaining or goal not in remaining:
        return False
    stack = [start]
    seen = {start}
    while stack:
        point = stack.pop()
        if point == goal:
            return True
        for neighbor in four_neighbors(point):
            if neighbor in remaining and neighbor not in seen:
                seen.add(neighbor)
                stack.append(neighbor)
    return False


def branched_corridor_universe() -> PointUniverse:
    universe = PointUniverse(7, 5, line_separation=0, square_clearance=0)
    universe.add_line(Line("WALL_TOP", "WALL", ((0, 0), (6, 0))))
    universe.add_line(Line("WALL_BOTTOM", "WALL", ((0, 4), (6, 4))))
    universe.add_line(Line("WALL_LEFT", "WALL", ((0, 0), (0, 4))))
    universe.add_line(Line("WALL_RIGHT", "WALL", ((6, 0), (6, 4))))
    universe.add_line(Line("BRANCH_LEFT_WALL", "WALL", ((1, 1), (1, 2))))
    universe.add_line(Line("BRANCH_RIGHT_WALL", "WALL", ((3, 1), (3, 2))))
    universe.add_line(Line("UPPER_RIGHT_WALL", "WALL", ((3, 2), (5, 2))))
    return universe


def test_terminal_minimum_vertex_cuts_find_all_two_point_cuts_in_two_passage_graph() -> None:
    universe = PointUniverse(3, 2, line_separation=0, square_clearance=0)
    viewer = WhiteSpaceViewer(universe)
    start = (0, 0)
    goal = (2, 1)
    expected = (
        ((1, 0), (0, 1)),
        ((1, 0), (1, 1)),
        ((2, 0), (1, 1)),
    )

    assert viewer.same_white_region(start, goal)
    assert viewer.terminal_separating_cut_points(start, goal) == ()
    for point in white_points(universe) - {start, goal}:
        assert connected_after_removing(universe, start, goal, {point})

    cuts = viewer.terminal_minimum_vertex_cuts(start, goal)

    assert cuts == expected
    for cut in cuts:
        assert len(cut) == 2
        assert start not in cut
        assert goal not in cut
        assert all(universe.point_state(point) & PointState.WHITE for point in cut)
        assert not connected_after_removing(universe, start, goal, set(cut))
        for point in cut:
            assert connected_after_removing(universe, start, goal, {point})


def test_terminal_minimum_vertex_cuts_agree_with_singleton_cut_points() -> None:
    universe = branched_corridor_universe()
    viewer = WhiteSpaceViewer(universe)
    start = (1, 3)
    goal = (5, 3)
    cut_points = viewer.terminal_separating_cut_points(start, goal)
    expected = tuple((point,) for point in cut_points)

    assert cut_points == ((2, 3), (3, 3), (4, 3))
    assert viewer.terminal_minimum_vertex_cuts(start, goal) == expected


def test_terminal_minimum_vertex_cuts_find_all_three_point_cuts_in_open_grid() -> None:
    universe = PointUniverse(3, 3, line_separation=0, square_clearance=0)
    viewer = WhiteSpaceViewer(universe)
    start = (0, 1)
    goal = (2, 1)
    expected = (
        ((0, 0), (1, 1), (0, 2)),
        ((0, 0), (1, 1), (1, 2)),
        ((0, 0), (1, 1), (2, 2)),
        ((1, 0), (1, 1), (0, 2)),
        ((1, 0), (1, 1), (1, 2)),
        ((1, 0), (1, 1), (2, 2)),
        ((2, 0), (1, 1), (0, 2)),
        ((2, 0), (1, 1), (1, 2)),
        ((2, 0), (1, 1), (2, 2)),
    )

    assert viewer.same_white_region(start, goal)
    assert viewer.terminal_separating_cut_points(start, goal) == ()
    for size in (1, 2):
        for removed in combinations(white_points(universe) - {start, goal}, size):
            assert connected_after_removing(universe, start, goal, set(removed))

    cuts = viewer.terminal_minimum_vertex_cuts(start, goal)

    assert cuts == expected
    for cut in cuts:
        assert len(cut) == 3
        assert not connected_after_removing(universe, start, goal, set(cut))


def test_terminal_minimum_vertex_cuts_residual_network_handles_bidirectional_adjacency() -> None:
    universe = PointUniverse(2, 3, line_separation=0, square_clearance=0)
    viewer = WhiteSpaceViewer(universe)
    start = (0, 0)
    goal = (1, 2)
    expected = (
        ((1, 0), (0, 1)),
        ((0, 1), (1, 1)),
        ((1, 1), (0, 2)),
    )

    assert viewer.same_white_region(start, goal)
    assert viewer.terminal_separating_cut_points(start, goal) == ()
    for point in white_points(universe) - {start, goal}:
        assert connected_after_removing(universe, start, goal, {point})

    cuts = viewer.terminal_minimum_vertex_cuts(start, goal)

    assert cuts == expected
    for cut in cuts:
        assert len(cut) == 2
        assert not connected_after_removing(universe, start, goal, set(cut))
        for point in cut:
            assert connected_after_removing(universe, start, goal, {point})


def test_terminal_minimum_vertex_cuts_do_not_mutate_authoritative_universe() -> None:
    universe = PointUniverse(3, 2, line_separation=0, square_clearance=0)
    before = snapshot(universe)

    WhiteSpaceViewer(universe).terminal_minimum_vertex_cuts((0, 0), (2, 1))

    assert snapshot(universe) == before


def test_terminal_minimum_vertex_cuts_input_semantics() -> None:
    universe = branched_corridor_universe()
    viewer = WhiteSpaceViewer(universe)

    with pytest.raises(ValueError, match="start point is not WHITE"):
        viewer.terminal_minimum_vertex_cuts((0, 0), (5, 3))

    with pytest.raises(ValueError, match="goal point is not WHITE"):
        viewer.terminal_minimum_vertex_cuts((1, 3), (0, 0))

    with pytest.raises(ValueError, match="start point outside"):
        viewer.terminal_minimum_vertex_cuts((-1, 3), (5, 3))

    with pytest.raises(ValueError, match="goal point outside"):
        viewer.terminal_minimum_vertex_cuts((1, 3), (7, 3))

    assert viewer.terminal_minimum_vertex_cuts((1, 3), (1, 3)) == ()

    disconnected = PointUniverse(3, 3, line_separation=0, square_clearance=0)
    disconnected.add_line(Line("WALL", "WALL", ((1, 0), (1, 2))))
    with pytest.raises(ValueError, match="start and goal must be in the same WHITE region"):
        WhiteSpaceViewer(disconnected).terminal_minimum_vertex_cuts((0, 1), (2, 1))

    adjacent = PointUniverse(2, 1, line_separation=0, square_clearance=0)
    assert WhiteSpaceViewer(adjacent).terminal_minimum_vertex_cuts((0, 0), (1, 0)) == ()
