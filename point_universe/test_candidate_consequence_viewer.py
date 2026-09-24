"""Tests for WhiteSpaceViewer hypothetical white-space consequences."""

from __future__ import annotations

import pytest

from .universe import Line, PointState, PointUniverse, Square, WhiteSpaceViewer


Point = tuple[int, int]


def obstacle_universe() -> PointUniverse:
    universe = PointUniverse(9, 9, line_separation=0, square_clearance=0)
    universe.add_square(Square("OBSTACLE", 3, 3, 2, 2))
    return universe


def matrix_snapshot(universe: PointUniverse) -> tuple[tuple[PointState, ...], ...]:
    return tuple(
        tuple(universe.point_state((x, y)) for x in range(universe.width))
        for y in range(universe.height)
    )


def authoritative_snapshot(universe: PointUniverse) -> tuple[
    tuple[tuple[PointState, ...], ...],
    tuple[str, ...],
    tuple[str, ...],
    tuple[tuple[Point, tuple[str, ...]], ...],
]:
    return (
        matrix_snapshot(universe),
        tuple(sorted(universe.squares)),
        tuple(sorted(universe.lines)),
        tuple(sorted((point, tuple(sorted(line_ids))) for point, line_ids in universe.point_lines.items())),
    )


def committed_region_data(candidate: tuple[Point, ...]) -> tuple[int, tuple[int, ...]]:
    universe = obstacle_universe()
    universe.add_line(Line("CANDIDATE", "N", candidate))
    viewer = WhiteSpaceViewer(universe)
    sizes = tuple(sorted(viewer.white_region_size(region_id) for region_id in range(1, viewer.white_region_count() + 1)))
    return viewer.white_region_count(), sizes


def test_hypothetical_consequence_reproduces_established_boundary_and_interior_cases() -> None:
    universe = obstacle_universe()
    snapshot = authoritative_snapshot(universe)
    viewer = WhiteSpaceViewer(universe)

    boundary_candidate = ((2, 3), (2, 4), (2, 5))
    interior_candidate = ((2, 3), (1, 3), (1, 4), (1, 5), (2, 5))

    boundary = viewer.hypothetical_white_consequence(boundary_candidate)
    interior = viewer.hypothetical_white_consequence(interior_candidate)

    assert boundary.R_before == 1
    assert boundary.R_after == 1
    assert boundary.delta_R == 0
    assert boundary.removed_white_points == 3

    assert interior.R_before == 1
    assert interior.R_after == 2
    assert interior.delta_R == 1
    assert interior.removed_white_points == 5
    assert 1 in interior.resulting_region_sizes

    committed_boundary_count, committed_boundary_sizes = committed_region_data(boundary_candidate)
    committed_interior_count, committed_interior_sizes = committed_region_data(interior_candidate)

    assert boundary.R_after == committed_boundary_count
    assert boundary.resulting_region_sizes == committed_boundary_sizes
    assert interior.R_after == committed_interior_count
    assert interior.resulting_region_sizes == committed_interior_sizes

    assert authoritative_snapshot(universe) == snapshot


def test_hypothetical_consequence_rejects_invalid_candidates() -> None:
    universe = obstacle_universe()
    viewer = WhiteSpaceViewer(universe)

    with pytest.raises(ValueError, match="empty"):
        viewer.hypothetical_white_consequence(())

    with pytest.raises(ValueError, match="outside universe"):
        viewer.hypothetical_white_consequence(((2, 3), (9, 3)))

    with pytest.raises(ValueError, match="not WHITE"):
        viewer.hypothetical_white_consequence(((3, 3),))

    with pytest.raises(ValueError, match="four-neighbor"):
        viewer.hypothetical_white_consequence(((2, 3), (2, 5)))

    with pytest.raises(ValueError, match="repeats"):
        viewer.hypothetical_white_consequence(((2, 3), (2, 4), (2, 3)))
