"""Source-only tests for two sequential point-universe connections."""

from __future__ import annotations

from .route_committer import RouteCommitter
from .route_selector import RouteSelector
from .router import PointRouter, RouteCandidate
from .universe import HypotheticalWhiteConsequence, Point, PointState, PointUniverse, Square, WhiteSpaceViewer


def corridor_universe() -> PointUniverse:
    universe = PointUniverse(7, 5, line_separation=0, square_clearance=0)
    universe.add_square(Square("CENTER", 2, 1, 2, 2))
    return universe


def region_sizes(viewer: WhiteSpaceViewer) -> tuple[int, ...]:
    return tuple(sorted(viewer.white_region_size(region_id) for region_id in range(1, viewer.white_region_count() + 1)))


def assert_committed_points_are_line_not_white(universe: PointUniverse, points: tuple[Point, ...]) -> None:
    for point in points:
        state = universe.point_state(point)
        assert state & PointState.LINE
        assert not (state & PointState.WHITE)


def test_second_connection_routes_against_universe_mutated_by_first_connection() -> None:
    universe = corridor_universe()
    a_start = (0, 2)
    a_goal = (6, 2)
    b_start = (0, 4)
    b_goal = (6, 4)
    for point in (a_start, a_goal, b_start, b_goal):
        assert universe.point_state(point) & PointState.WHITE

    a_candidates = PointRouter(universe).candidate_paths(a_start, a_goal, max_paths=2)
    assert len(a_candidates) >= 2
    for candidate in a_candidates:
        assert candidate.points[0] == a_start
        assert candidate.points[-1] == a_goal
        assert candidate.consequence == WhiteSpaceViewer(universe).hypothetical_white_consequence(candidate.points)

    selected_a = RouteSelector().select(a_candidates)
    selected_a_number = a_candidates.index(selected_a) + 1
    assert selected_a_number in (1, 2)
    a_line = RouteCommitter(universe).commit(selected_a, line_id="A", net_id="NET_A")

    assert "A" in universe.lines
    assert a_line.points == selected_a.points
    assert_committed_points_are_line_not_white(universe, selected_a.points)

    post_a_viewer = WhiteSpaceViewer(universe)
    post_a_count = post_a_viewer.white_region_count()
    post_a_delta = post_a_count - selected_a.consequence.R_before
    assert post_a_count == selected_a.consequence.R_after
    assert post_a_delta == selected_a.consequence.delta_R
    post_a_sizes = region_sizes(post_a_viewer)

    a_line_object_before_b = universe.lines["A"]
    a_points_before_b = universe.lines["A"].points
    a_line_points_before_b = set(universe.line_points["A"])

    assert "A" in universe.lines
    b_router = PointRouter(universe)
    b_candidates = b_router.candidate_paths(b_start, b_goal, max_paths=2)
    assert len(b_candidates) >= 1
    for candidate in b_candidates:
        assert candidate.consequence.R_before == post_a_count
        for point in candidate.points:
            assert "A" not in universe.point_lines.get(point, set())

    selected_b = RouteSelector().select(b_candidates)
    b_line = RouteCommitter(universe).commit(selected_b, line_id="B", net_id="NET_B")

    assert "B" in universe.lines
    assert b_line.points == selected_b.points
    assert_committed_points_are_line_not_white(universe, selected_b.points)
    assert universe.lines["A"] is a_line_object_before_b
    assert universe.lines["A"].points == a_points_before_b
    assert universe.line_points["A"] == a_line_points_before_b

    final_viewer = WhiteSpaceViewer(universe)
    final_count = final_viewer.white_region_count()
    final_delta = final_count - selected_b.consequence.R_before
    assert final_count == selected_b.consequence.R_after
    assert final_delta == selected_b.consequence.delta_R
    assert region_sizes(final_viewer) == selected_b.consequence.resulting_region_sizes

    # Retain observable values in assertions so the test represents the
    # complete sequential mechanism without adding a new routing algorithm.
    assert len(a_candidates) == 2
    assert len(b_candidates) >= 1
    assert post_a_sizes == selected_a.consequence.resulting_region_sizes


def test_blocked_second_connection_reports_no_candidate_without_rerouting_first_connection() -> None:
    universe = PointUniverse(3, 3, line_separation=0, square_clearance=0)
    a_candidate = RouteCandidate(
        points=((1, 0), (1, 1), (1, 2)),
        consequence=HypotheticalWhiteConsequence(
            R_before=1,
            R_after=2,
            delta_R=1,
            removed_white_points=3,
            resulting_region_sizes=(3, 3),
        ),
    )

    a_line = RouteCommitter(universe).commit(a_candidate, line_id="A", net_id="NET_A")
    a_points = a_line.points
    a_line_points = set(universe.line_points["A"])

    b_candidates = PointRouter(universe).candidate_paths((0, 1), (2, 1), max_paths=2)

    assert b_candidates == []
    assert universe.lines["A"] is a_line
    assert universe.lines["A"].points == a_points
    assert universe.line_points["A"] == a_line_points
