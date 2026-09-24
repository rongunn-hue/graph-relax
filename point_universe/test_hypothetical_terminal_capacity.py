"""Source-only tests for hypothetical terminal capacity consequences."""

from __future__ import annotations

from itertools import combinations

import pytest

from .universe import Line, Point, PointState, PointUniverse, WhiteSpaceViewer, four_neighbors


def capacity_universe() -> PointUniverse:
    universe = PointUniverse(5, 3, line_separation=0, square_clearance=0)
    universe.add_line(Line("WALL", "WALL", ((3, 0), (3, 2))))
    return universe


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


def connected_in_points(white: set[Point], start: Point, goal: Point) -> bool:
    if start not in white or goal not in white:
        return False
    stack = [start]
    seen = {start}
    while stack:
        point = stack.pop()
        if point == goal:
            return True
        for neighbor in four_neighbors(point):
            if neighbor in white and neighbor not in seen:
                seen.add(neighbor)
                stack.append(neighbor)
    return False


def expected_minimum_cuts(white: set[Point], start: Point, goal: Point) -> tuple[tuple[Point, ...], ...]:
    nonterminals = tuple(sorted(white - {start, goal}, key=lambda p: (p[1], p[0])))
    for size in range(1, len(nonterminals) + 1):
        cuts: list[tuple[Point, ...]] = []
        for candidate in combinations(nonterminals, size):
            if not connected_in_points(white - set(candidate), start, goal):
                cuts.append(tuple(sorted(candidate, key=lambda p: (p[1], p[0]))))
        if cuts:
            return tuple(sorted(cuts, key=lambda cut: tuple((point[1], point[0]) for point in cut)))
    return ()


B_START = (0, 1)
B_GOAL = (2, 1)
MINIMUM_CUTS_BEFORE = (
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
UNCHANGED_CANDIDATE = ((4, 1),)
REDUCED_CANDIDATE = ((1, 1),)
DISCONNECTING_CANDIDATE = ((0, 0), (1, 0), (1, 1), (1, 2), (0, 2))
REDUCED_CUTS_AFTER = (
    ((0, 0), (0, 2)),
    ((0, 0), (1, 2)),
    ((0, 0), (2, 2)),
    ((1, 0), (0, 2)),
    ((1, 0), (1, 2)),
    ((1, 0), (2, 2)),
    ((2, 0), (0, 2)),
    ((2, 0), (1, 2)),
    ((2, 0), (2, 2)),
)


def assert_capacity_oracle(
    universe: PointUniverse,
    candidate: tuple[Point, ...],
    expected_cuts: tuple[tuple[Point, ...], ...],
) -> None:
    white_after = white_points(universe) - set(candidate)
    assert connected_in_points(white_after, B_START, B_GOAL)
    assert expected_minimum_cuts(white_after, B_START, B_GOAL) == expected_cuts
    capacity = len(expected_cuts[0])
    nonterminals = white_after - {B_START, B_GOAL}
    for size in range(1, capacity):
        for removed in combinations(nonterminals, size):
            assert connected_in_points(white_after - set(removed), B_START, B_GOAL)
    for cut in expected_cuts:
        assert not connected_in_points(white_after - set(cut), B_START, B_GOAL)


def test_hypothetical_terminal_capacity_distinguishes_unchanged_reduced_and_disconnected_cases() -> None:
    universe = capacity_universe()
    viewer = WhiteSpaceViewer(universe)
    before_snapshot = snapshot(universe)
    candidates_before = (UNCHANGED_CANDIDATE, REDUCED_CANDIDATE, DISCONNECTING_CANDIDATE)

    assert viewer.terminal_minimum_vertex_cuts(B_START, B_GOAL) == MINIMUM_CUTS_BEFORE
    assert expected_minimum_cuts(white_points(universe), B_START, B_GOAL) == MINIMUM_CUTS_BEFORE

    unchanged = viewer.hypothetical_terminal_capacity_consequence(
        UNCHANGED_CANDIDATE,
        B_start=B_START,
        B_goal=B_GOAL,
    )
    reduced = viewer.hypothetical_terminal_capacity_consequence(
        REDUCED_CANDIDATE,
        B_start=B_START,
        B_goal=B_GOAL,
    )
    disconnected = viewer.hypothetical_terminal_capacity_consequence(
        DISCONNECTING_CANDIDATE,
        B_start=B_START,
        B_goal=B_GOAL,
    )

    assert unchanged.capacity_before == 3
    assert unchanged.capacity_after == 3
    assert unchanged.delta_capacity == 0
    assert unchanged.B_routable_after is True
    assert unchanged.finite_cut_after is True
    assert unchanged.minimum_cuts_before == MINIMUM_CUTS_BEFORE
    assert unchanged.minimum_cuts_after == MINIMUM_CUTS_BEFORE
    assert_capacity_oracle(universe, UNCHANGED_CANDIDATE, MINIMUM_CUTS_BEFORE)

    assert reduced.capacity_before == 3
    assert reduced.capacity_after == 2
    assert reduced.delta_capacity == -1
    assert reduced.B_routable_after is True
    assert reduced.finite_cut_after is True
    assert reduced.minimum_cuts_before == MINIMUM_CUTS_BEFORE
    assert reduced.minimum_cuts_after == REDUCED_CUTS_AFTER
    assert_capacity_oracle(universe, REDUCED_CANDIDATE, REDUCED_CUTS_AFTER)

    assert disconnected.capacity_before == 3
    assert disconnected.capacity_after == 0
    assert disconnected.delta_capacity == -3
    assert disconnected.B_routable_after is False
    assert disconnected.finite_cut_after is False
    assert disconnected.minimum_cuts_before == MINIMUM_CUTS_BEFORE
    assert disconnected.minimum_cuts_after == ()
    white_after_disconnection = white_points(universe) - set(DISCONNECTING_CANDIDATE)
    assert B_START in white_after_disconnection
    assert B_GOAL in white_after_disconnection
    assert not connected_in_points(white_after_disconnection, B_START, B_GOAL)

    assert snapshot(universe) == before_snapshot
    assert candidates_before == (UNCHANGED_CANDIDATE, REDUCED_CANDIDATE, DISCONNECTING_CANDIDATE)


def test_hypothetical_terminal_capacity_rejects_invalid_inputs() -> None:
    universe = capacity_universe()
    viewer = WhiteSpaceViewer(universe)

    with pytest.raises(ValueError, match="B_start point outside"):
        viewer.hypothetical_terminal_capacity_consequence(UNCHANGED_CANDIDATE, B_start=(-1, 1), B_goal=B_GOAL)

    with pytest.raises(ValueError, match="B_goal point outside"):
        viewer.hypothetical_terminal_capacity_consequence(UNCHANGED_CANDIDATE, B_start=B_START, B_goal=(5, 1))

    with pytest.raises(ValueError, match="B_start point is not WHITE"):
        viewer.hypothetical_terminal_capacity_consequence(UNCHANGED_CANDIDATE, B_start=(3, 1), B_goal=B_GOAL)

    with pytest.raises(ValueError, match="B_goal point is not WHITE"):
        viewer.hypothetical_terminal_capacity_consequence(UNCHANGED_CANDIDATE, B_start=B_START, B_goal=(3, 1))

    with pytest.raises(ValueError, match="B_start and B_goal must be distinct"):
        viewer.hypothetical_terminal_capacity_consequence(UNCHANGED_CANDIDATE, B_start=B_START, B_goal=B_START)

    with pytest.raises(ValueError, match="B_start and B_goal must be in the same WHITE region"):
        viewer.hypothetical_terminal_capacity_consequence(UNCHANGED_CANDIDATE, B_start=B_START, B_goal=(4, 1))

    adjacent = PointUniverse(3, 1, line_separation=0, square_clearance=0)
    with pytest.raises(ValueError, match="finite non-terminal terminal cut is required"):
        WhiteSpaceViewer(adjacent).hypothetical_terminal_capacity_consequence(
            ((2, 0),),
            B_start=(0, 0),
            B_goal=(1, 0),
        )

    with pytest.raises(ValueError, match="hypothetical candidate is empty"):
        viewer.hypothetical_terminal_capacity_consequence((), B_start=B_START, B_goal=B_GOAL)

    with pytest.raises(ValueError, match="hypothetical candidate repeats point"):
        viewer.hypothetical_terminal_capacity_consequence(((4, 0), (4, 0)), B_start=B_START, B_goal=B_GOAL)

    with pytest.raises(ValueError, match="hypothetical candidate point is not WHITE"):
        viewer.hypothetical_terminal_capacity_consequence(((3, 1),), B_start=B_START, B_goal=B_GOAL)

    with pytest.raises(ValueError, match="hypothetical candidate point outside universe"):
        viewer.hypothetical_terminal_capacity_consequence(((5, 1),), B_start=B_START, B_goal=B_GOAL)

    with pytest.raises(ValueError, match="hypothetical candidate points must be four-neighbor adjacent"):
        viewer.hypothetical_terminal_capacity_consequence(((4, 0), (4, 2)), B_start=B_START, B_goal=B_GOAL)

    with pytest.raises(ValueError, match="hypothetical candidate may not include B_start"):
        viewer.hypothetical_terminal_capacity_consequence((B_START,), B_start=B_START, B_goal=B_GOAL)

    with pytest.raises(ValueError, match="hypothetical candidate may not include B_goal"):
        viewer.hypothetical_terminal_capacity_consequence((B_GOAL,), B_start=B_START, B_goal=B_GOAL)
