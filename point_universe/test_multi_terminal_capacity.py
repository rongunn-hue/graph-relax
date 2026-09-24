"""Source-only tests for multi-terminal hypothetical capacity perception."""

from __future__ import annotations

from itertools import combinations

import pytest

from .universe import (
    Line,
    Point,
    PointState,
    PointUniverse,
    TerminalPair,
    WhiteSpaceViewer,
    four_neighbors,
)


def multi_capacity_universe() -> PointUniverse:
    universe = PointUniverse(7, 3, line_separation=0, square_clearance=0)
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


CANDIDATE_A = ((0, 0), (1, 0), (1, 1), (1, 2), (0, 2))
PAIR_B = TerminalPair("B", (2, 0), (2, 2))
PAIR_C = TerminalPair("C", (4, 1), (6, 1))
PAIR_D = TerminalPair("D", (0, 1), (2, 1))
PRIMARY_PAIRS = (PAIR_B, PAIR_C, PAIR_D)

B_CUTS_BEFORE = (
    ((1, 0), (2, 1)),
    ((2, 1), (1, 2)),
)
B_CUTS_AFTER = (((2, 1),),)
C_CUTS = (
    ((4, 0), (5, 1), (4, 2)),
    ((4, 0), (5, 1), (5, 2)),
    ((4, 0), (5, 1), (6, 2)),
    ((5, 0), (5, 1), (4, 2)),
    ((5, 0), (5, 1), (5, 2)),
    ((5, 0), (5, 1), (6, 2)),
    ((6, 0), (5, 1), (4, 2)),
    ((6, 0), (5, 1), (5, 2)),
    ((6, 0), (5, 1), (6, 2)),
)
D_CUTS_BEFORE = (
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
EXPECTED_BY_ID = {
    "B": (2, 1, -1, True, True, B_CUTS_BEFORE, B_CUTS_AFTER),
    "C": (3, 3, 0, True, True, C_CUTS, C_CUTS),
    "D": (3, 0, -3, False, False, D_CUTS_BEFORE, ()),
}


def assert_pair_oracle(
    universe: PointUniverse,
    pair: TerminalPair,
    expected_before: tuple[tuple[Point, ...], ...],
    expected_after: tuple[tuple[Point, ...], ...],
    connected_after: bool,
) -> None:
    before = white_points(universe)
    after = before - set(CANDIDATE_A)
    assert pair.start in before
    assert pair.goal in before
    assert connected_in_points(before, pair.start, pair.goal)
    assert expected_minimum_cuts(before, pair.start, pair.goal) == expected_before
    assert pair.start in after
    assert pair.goal in after
    assert connected_in_points(after, pair.start, pair.goal) is connected_after
    if connected_after:
        assert expected_minimum_cuts(after, pair.start, pair.goal) == expected_after
        capacity_after = len(expected_after[0])
        for size in range(1, capacity_after):
            for removed in combinations(after - {pair.start, pair.goal}, size):
                assert connected_in_points(after - set(removed), pair.start, pair.goal)
        for cut in expected_after:
            assert not connected_in_points(after - set(cut), pair.start, pair.goal)
    else:
        assert expected_after == ()


def assert_consequence(pair_id: str, consequence) -> None:
    before, after, delta, routable, finite, cuts_before, cuts_after = EXPECTED_BY_ID[pair_id]
    assert consequence.capacity_before == before
    assert consequence.capacity_after == after
    assert consequence.delta_capacity == delta
    assert consequence.B_routable_after is routable
    assert consequence.finite_cut_after is finite
    assert consequence.minimum_cuts_before == cuts_before
    assert consequence.minimum_cuts_after == cuts_after


def test_multi_terminal_capacity_reports_multiple_pair_consequences_without_mutation() -> None:
    universe = multi_capacity_universe()
    viewer = WhiteSpaceViewer(universe)
    before_snapshot = snapshot(universe)
    candidate_before = CANDIDATE_A
    pairs_before = PRIMARY_PAIRS

    result = viewer.hypothetical_multi_terminal_capacity_consequence(CANDIDATE_A, PRIMARY_PAIRS)

    assert result.candidate_points == CANDIDATE_A
    assert tuple(pair_id for pair_id, _ in result.pair_consequences) == ("B", "C", "D")
    for pair_id, consequence in result.pair_consequences:
        assert_consequence(pair_id, consequence)

    assert_pair_oracle(universe, PAIR_B, B_CUTS_BEFORE, B_CUTS_AFTER, True)
    assert_pair_oracle(universe, PAIR_C, C_CUTS, C_CUTS, True)
    assert_pair_oracle(universe, PAIR_D, D_CUTS_BEFORE, (), False)
    assert snapshot(universe) == before_snapshot
    assert CANDIDATE_A == candidate_before
    assert PRIMARY_PAIRS == pairs_before


def test_multi_terminal_capacity_preserves_pair_order_and_matches_single_pair_results() -> None:
    universe = multi_capacity_universe()
    viewer = WhiteSpaceViewer(universe)
    reordered = (PAIR_D, PAIR_B, PAIR_C)

    first = viewer.hypothetical_multi_terminal_capacity_consequence(CANDIDATE_A, PRIMARY_PAIRS)
    second = viewer.hypothetical_multi_terminal_capacity_consequence(CANDIDATE_A, reordered)

    assert tuple(pair_id for pair_id, _ in first.pair_consequences) == ("B", "C", "D")
    assert tuple(pair_id for pair_id, _ in second.pair_consequences) == ("D", "B", "C")
    by_id = dict(first.pair_consequences)
    for pair_id, consequence in second.pair_consequences:
        assert consequence == by_id[pair_id]

    for pair in PRIMARY_PAIRS:
        single = viewer.hypothetical_terminal_capacity_consequence(
            CANDIDATE_A,
            B_start=pair.start,
            B_goal=pair.goal,
        )
        assert by_id[pair.pair_id] == single


def test_multi_terminal_capacity_rejects_invalid_pairs_and_candidates_without_mutation() -> None:
    universe = multi_capacity_universe()
    viewer = WhiteSpaceViewer(universe)
    before = snapshot(universe)

    with pytest.raises(ValueError, match="at least one TerminalPair is required"):
        viewer.hypothetical_multi_terminal_capacity_consequence(CANDIDATE_A, ())

    with pytest.raises(ValueError, match="pair_id must be non-empty"):
        viewer.hypothetical_multi_terminal_capacity_consequence(CANDIDATE_A, (TerminalPair("", (2, 0), (2, 2)),))

    with pytest.raises(ValueError, match="pair_id must be a string"):
        viewer.hypothetical_multi_terminal_capacity_consequence(CANDIDATE_A, (TerminalPair(3, (2, 0), (2, 2)),))

    with pytest.raises(ValueError, match="duplicate pair_id"):
        viewer.hypothetical_multi_terminal_capacity_consequence(CANDIDATE_A, (PAIR_B, TerminalPair("B", (4, 1), (6, 1))))

    with pytest.raises(ValueError, match="start point outside"):
        viewer.hypothetical_multi_terminal_capacity_consequence(CANDIDATE_A, (TerminalPair("X", (-1, 0), (2, 2)),))

    with pytest.raises(ValueError, match="goal point outside"):
        viewer.hypothetical_multi_terminal_capacity_consequence(CANDIDATE_A, (TerminalPair("X", (2, 0), (7, 2)),))

    with pytest.raises(ValueError, match="start point is not WHITE"):
        viewer.hypothetical_multi_terminal_capacity_consequence(CANDIDATE_A, (TerminalPair("X", (3, 0), (2, 2)),))

    with pytest.raises(ValueError, match="goal point is not WHITE"):
        viewer.hypothetical_multi_terminal_capacity_consequence(CANDIDATE_A, (TerminalPair("X", (2, 0), (3, 0)),))

    with pytest.raises(ValueError, match="terminal pair endpoints must be distinct"):
        viewer.hypothetical_multi_terminal_capacity_consequence(CANDIDATE_A, (TerminalPair("X", (2, 0), (2, 0)),))

    with pytest.raises(ValueError, match="terminal pair endpoints must be in the same WHITE region"):
        viewer.hypothetical_multi_terminal_capacity_consequence(CANDIDATE_A, (TerminalPair("X", (2, 0), (4, 1)),))

    adjacent = PointUniverse(3, 1, line_separation=0, square_clearance=0)
    with pytest.raises(ValueError, match="a finite non-terminal terminal cut is required"):
        WhiteSpaceViewer(adjacent).hypothetical_multi_terminal_capacity_consequence(
            ((2, 0),),
            (TerminalPair("X", (0, 0), (1, 0)),),
        )

    with pytest.raises(ValueError, match="hypothetical candidate may not include terminal pair start"):
        viewer.hypothetical_multi_terminal_capacity_consequence(CANDIDATE_A, (TerminalPair("X", (1, 0), (2, 2)),))

    with pytest.raises(ValueError, match="hypothetical candidate may not include terminal pair goal"):
        viewer.hypothetical_multi_terminal_capacity_consequence(CANDIDATE_A, (TerminalPair("X", (2, 2), (1, 0)),))

    with pytest.raises(ValueError, match="hypothetical candidate is empty"):
        viewer.hypothetical_multi_terminal_capacity_consequence((), PRIMARY_PAIRS)

    with pytest.raises(ValueError, match="hypothetical candidate repeats point"):
        viewer.hypothetical_multi_terminal_capacity_consequence(((4, 0), (4, 0)), PRIMARY_PAIRS)

    with pytest.raises(ValueError, match="hypothetical candidate point is not WHITE"):
        viewer.hypothetical_multi_terminal_capacity_consequence(((3, 1),), PRIMARY_PAIRS)

    with pytest.raises(ValueError, match="hypothetical candidate point outside universe"):
        viewer.hypothetical_multi_terminal_capacity_consequence(((7, 1),), PRIMARY_PAIRS)

    with pytest.raises(ValueError, match="hypothetical candidate points must be four-neighbor adjacent"):
        viewer.hypothetical_multi_terminal_capacity_consequence(((4, 0), (4, 2)), PRIMARY_PAIRS)

    assert snapshot(universe) == before
