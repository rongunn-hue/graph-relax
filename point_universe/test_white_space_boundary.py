"""Unexecuted tests for white-space boundary perception."""

from __future__ import annotations

from .universe import Line, PointState, PointUniverse, Square, WhiteSpaceViewer


def assert_boundary_path_laws(viewer: WhiteSpaceViewer, region_id: int, path: list[tuple[int, int]]) -> None:
    boundary = viewer.white_region_boundary(region_id)
    assert path
    for point in path:
        assert viewer.universe.point_state(point) & PointState.WHITE
        assert viewer.white_region_at(*point) == region_id
        assert point in boundary
    for a, b in zip(path, path[1:]):
        assert abs(a[0] - b[0]) + abs(a[1] - b[1]) == 1


def square_obstacle_universe() -> PointUniverse:
    universe = PointUniverse(9, 9, line_separation=0, square_clearance=0)
    universe.add_square(Square("OBSTACLE", 3, 3, 2, 2))
    return universe


def has_square_neighbor(universe: PointUniverse, point: tuple[int, int]) -> bool:
    x, y = point
    for neighbor in ((x, y - 1), (x, y + 1), (x - 1, y), (x + 1, y)):
        if universe.in_bounds(neighbor) and universe.point_state(neighbor) & PointState.SQUARE:
            return True
    return False


def test_boundary_points_are_derived_from_matrix_white_state() -> None:
    universe = PointUniverse(5, 5, line_separation=0, square_clearance=0)
    viewer = WhiteSpaceViewer(universe)
    region_id = viewer.white_region_at(0, 0)
    boundary = viewer.white_region_boundary(region_id)
    assert (0, 0) in boundary
    assert (2, 2) not in boundary
    assert universe.point_state((0, 0)) & PointState.WHITE
    assert universe.point_state((2, 2)) & PointState.WHITE


def test_four_neighbor_boundary_connectivity_and_path_validity() -> None:
    universe = PointUniverse(5, 5, line_separation=0, square_clearance=0)
    viewer = WhiteSpaceViewer(universe)
    region_id = viewer.white_region_at(0, 0)
    a = (0, 0)
    b = (4, 0)
    assert viewer.same_boundary_component(region_id, a, b)
    path = viewer.boundary_path(region_id, a, b)
    assert path[0] == a
    assert path[-1] == b
    assert_boundary_path_laws(viewer, region_id, path)


def test_diagonal_only_contact_does_not_establish_boundary_connectivity() -> None:
    universe = PointUniverse(3, 3, line_separation=0, square_clearance=0)
    universe.add_line(Line("DIAGONAL_SEPARATOR", "CUT", ((1, 0), (2, 0), (2, 2), (0, 2), (0, 1))))
    viewer = WhiteSpaceViewer(universe)
    a = (0, 0)
    b = (1, 1)
    assert universe.point_state(a) & PointState.WHITE
    assert universe.point_state(b) & PointState.WHITE
    assert abs(a[0] - b[0]) == 1
    assert abs(a[1] - b[1]) == 1
    region_a = viewer.white_region_at(*a)
    region_b = viewer.white_region_at(*b)
    assert region_a != region_b
    assert a in viewer.white_region_boundary(region_a)
    assert b in viewer.white_region_boundary(region_b)
    assert not viewer.same_boundary_component(region_a, a, b)
    assert viewer.boundary_path(region_a, a, b) is None


def test_no_path_between_separate_boundary_components() -> None:
    universe = square_obstacle_universe()
    viewer = WhiteSpaceViewer(universe)
    region_id = viewer.white_region_at(0, 0)
    outer_boundary_point = (0, 0)
    object_boundary_point = (2, 3)
    assert viewer.white_region_at(*outer_boundary_point) == region_id
    assert viewer.white_region_at(*object_boundary_point) == region_id
    assert has_square_neighbor(universe, object_boundary_point)
    assert len(viewer.boundary_components(region_id)) >= 2
    assert not viewer.same_boundary_component(region_id, outer_boundary_point, object_boundary_point)
    assert viewer.boundary_path(region_id, outer_boundary_point, object_boundary_point) is None


def test_boundary_following_line_preserves_white_region_count() -> None:
    before = square_obstacle_universe()
    before_viewer = WhiteSpaceViewer(before)
    a = (2, 3)
    b = (2, 5)
    region_id = before_viewer.white_region_at(*a)
    assert before.point_state(a) & PointState.WHITE
    assert before.point_state(b) & PointState.WHITE
    assert before_viewer.white_region_at(*b) == region_id
    assert a in before_viewer.white_region_boundary(region_id)
    assert b in before_viewer.white_region_boundary(region_id)
    assert before_viewer.same_boundary_component(region_id, a, b)
    path = before_viewer.boundary_path(region_id, a, b)
    assert path[0] == a
    assert path[-1] == b
    assert_boundary_path_laws(before_viewer, region_id, path)
    assert any(has_square_neighbor(before, point) for point in path)

    after = square_obstacle_universe()
    after.add_line(Line("BOUNDARY_CONNECTION", "N", tuple(path)))
    after_viewer = WhiteSpaceViewer(after)

    r_before = before_viewer.white_region_count()
    r_after = after_viewer.white_region_count()
    assert r_after - r_before == 0


def test_interior_cut_increases_white_region_count() -> None:
    before = square_obstacle_universe()
    before_viewer = WhiteSpaceViewer(before)
    r_before = before_viewer.white_region_count()

    after = square_obstacle_universe()
    after.add_line(Line("INTERIOR_CUT", "N", ((0, 1), (8, 1))))
    after_viewer = WhiteSpaceViewer(after)
    r_after = after_viewer.white_region_count()

    assert r_after - r_before > 0
