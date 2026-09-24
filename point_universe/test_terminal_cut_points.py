"""Source-only tests for terminal-specific WHITE cut-point perception."""

from __future__ import annotations

import pytest

from .universe import Line, Point, PointState, PointUniverse, WhiteSpaceViewer


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


def connected_without(universe: PointUniverse, start: Point, goal: Point, removed: Point) -> bool:
    if start == removed or goal == removed:
        return False
    remaining = set(universe.layer(PointState.WHITE))
    remaining.discard(removed)
    if start not in remaining or goal not in remaining:
        return False
    stack = [start]
    seen = {start}
    while stack:
        x, y = stack.pop()
        for neighbor in ((x, y - 1), (x - 1, y), (x + 1, y), (x, y + 1)):
            if neighbor in remaining and neighbor not in seen:
                seen.add(neighbor)
                stack.append(neighbor)
    return goal in seen


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


def test_terminal_cut_points_identify_chokepoint_between_terminals() -> None:
    universe = branched_corridor_universe()
    viewer = WhiteSpaceViewer(universe)
    start = (1, 3)
    goal = (5, 3)
    chokepoint = (3, 3)
    unrelated_articulation = (2, 2)
    branch_leaf = (2, 1)

    assert universe.point_state(start) & PointState.WHITE
    assert universe.point_state(goal) & PointState.WHITE
    assert universe.point_state(chokepoint) & PointState.WHITE
    assert universe.point_state(unrelated_articulation) & PointState.WHITE
    assert universe.point_state(branch_leaf) & PointState.WHITE
    assert viewer.same_white_region(start, goal)
    assert viewer.same_white_region(start, branch_leaf)
    assert not connected_without(universe, start, branch_leaf, unrelated_articulation)
    assert connected_without(universe, start, goal, unrelated_articulation)

    cut_points = viewer.terminal_separating_cut_points(start, goal)

    assert chokepoint in cut_points
    assert unrelated_articulation not in cut_points
    assert start not in cut_points
    assert goal not in cut_points
    assert cut_points == ((2, 3), (3, 3), (4, 3))


def test_terminal_cut_points_empty_when_redundant_passages_exist() -> None:
    universe = PointUniverse(3, 3, line_separation=0, square_clearance=0)
    viewer = WhiteSpaceViewer(universe)
    assert viewer.terminal_separating_cut_points((0, 1), (2, 1)) == ()


def test_terminal_cut_points_do_not_mutate_authoritative_universe() -> None:
    universe = branched_corridor_universe()
    before = snapshot(universe)

    WhiteSpaceViewer(universe).terminal_separating_cut_points((1, 3), (5, 3))

    assert snapshot(universe) == before


def test_terminal_cut_points_input_semantics() -> None:
    universe = branched_corridor_universe()
    viewer = WhiteSpaceViewer(universe)

    with pytest.raises(ValueError, match="start point is not WHITE"):
        viewer.terminal_separating_cut_points((0, 0), (5, 3))

    with pytest.raises(ValueError, match="goal point is not WHITE"):
        viewer.terminal_separating_cut_points((1, 3), (0, 0))

    with pytest.raises(ValueError, match="start point outside"):
        viewer.terminal_separating_cut_points((-1, 3), (5, 3))

    with pytest.raises(ValueError, match="goal point outside"):
        viewer.terminal_separating_cut_points((1, 3), (7, 3))

    assert viewer.terminal_separating_cut_points((1, 3), (1, 3)) == ()
