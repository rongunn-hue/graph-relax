"""Source-only tests for committing one selected route candidate."""

from __future__ import annotations

import pytest

from .route_committer import RouteCommitter
from .route_selector import RouteSelector
from .router import PointRouter, RouteCandidate
from .universe import HypotheticalWhiteConsequence, Point, PointState, PointUniverse, Square, WhiteSpaceViewer


def obstacle_universe() -> PointUniverse:
    universe = PointUniverse(7, 7, line_separation=0, square_clearance=0)
    universe.add_square(Square("OBSTACLE", 2, 2, 2, 2))
    return universe


def matrix_snapshot(universe: PointUniverse) -> tuple[tuple[PointState, ...], ...]:
    return tuple(
        tuple(universe.point_state((x, y)) for x in range(universe.width))
        for y in range(universe.height)
    )


def region_sizes(viewer: WhiteSpaceViewer) -> tuple[int, ...]:
    return tuple(sorted(viewer.white_region_size(region_id) for region_id in range(1, viewer.white_region_count() + 1)))


def selected_candidate_case() -> tuple[PointUniverse, list[RouteCandidate], RouteCandidate]:
    universe = obstacle_universe()
    candidates = PointRouter(universe).candidate_paths((0, 3), (6, 3), max_paths=2)
    selected = RouteSelector().select(candidates)
    return universe, candidates, selected


def test_selected_candidate_commits_exact_path_and_matches_hypothetical_consequence() -> None:
    universe, candidates, selected = selected_candidate_case()
    nonselected = next(candidate for candidate in candidates if candidate is not selected)
    selected_points_before = selected.points
    selected_consequence_before = selected.consequence
    nonselected_only = set(nonselected.points) - set(selected.points)
    white_nonselected_only_before = {
        point for point in nonselected_only if universe.point_state(point) & PointState.WHITE
    }

    assert len(selected.points) == 31
    assert selected.consequence.R_before == 1
    assert selected.consequence.R_after == 1
    assert selected.consequence.delta_R == 0

    line = RouteCommitter(universe).commit(selected, line_id="SELECTED", net_id="N")

    assert universe.lines["SELECTED"] is line
    assert line.points == selected.points
    assert universe.line_points["SELECTED"] == set(selected.points)

    after_viewer = WhiteSpaceViewer(universe)
    actual_r_after = after_viewer.white_region_count()
    actual_delta_r = actual_r_after - selected.consequence.R_before
    assert actual_r_after == selected.consequence.R_after
    assert actual_r_after == 1
    assert actual_delta_r == selected.consequence.delta_R
    assert actual_delta_r == 0
    assert region_sizes(after_viewer) == selected.consequence.resulting_region_sizes

    for point in selected.points:
        state = universe.point_state(point)
        assert state & PointState.LINE
        assert not (state & PointState.WHITE)

    for point in white_nonselected_only_before:
        assert "SELECTED" not in universe.point_lines.get(point, set())

    assert selected.points == selected_points_before
    assert selected.consequence == selected_consequence_before


def test_invalid_commit_is_rejected_without_partial_mutation() -> None:
    universe = obstacle_universe()
    before_matrix = matrix_snapshot(universe)
    invalid = RouteCandidate(
        points=((0, 2), (1, 2), (2, 2), (3, 2)),
        consequence=HypotheticalWhiteConsequence(
            R_before=1,
            R_after=1,
            delta_R=0,
            removed_white_points=4,
            resulting_region_sizes=(),
        ),
    )

    with pytest.raises(ValueError):
        RouteCommitter(universe).commit(invalid, line_id="INVALID", net_id="N")

    assert "INVALID" not in universe.lines
    assert "INVALID" not in universe.line_points
    assert matrix_snapshot(universe) == before_matrix


def test_duplicate_line_id_is_rejected_without_replacing_original_line() -> None:
    universe = PointUniverse(7, 7, line_separation=0, square_clearance=0)
    first = RouteCandidate(
        points=((0, 0), (1, 0)),
        consequence=HypotheticalWhiteConsequence(1, 1, 0, 2, ()),
    )
    duplicate = RouteCandidate(
        points=((0, 1), (1, 1)),
        consequence=HypotheticalWhiteConsequence(1, 1, 0, 2, ()),
    )
    committer = RouteCommitter(universe)

    original_line = committer.commit(first, line_id="DUPLICATE", net_id="N")
    original_points = original_line.points

    with pytest.raises(ValueError):
        committer.commit(duplicate, line_id="DUPLICATE", net_id="N")

    assert universe.lines["DUPLICATE"] is original_line
    assert universe.lines["DUPLICATE"].points == original_points
    assert universe.line_points["DUPLICATE"] == set(original_points)
