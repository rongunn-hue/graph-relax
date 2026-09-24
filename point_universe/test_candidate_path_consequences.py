"""Source-only tests for pre-commit candidate path consequences.

These tests intentionally do not add routing policy.  They define two explicit
legal candidate paths for the same connection and measure the white-space
topology that would result if each candidate were committed in isolation.
"""

from __future__ import annotations

from .universe import Line, PointState, PointUniverse, Square, WhiteSpaceViewer


Point = tuple[int, int]


def obstacle_universe() -> PointUniverse:
    universe = PointUniverse(9, 9, line_separation=0, square_clearance=0)
    universe.add_square(Square("OBSTACLE", 3, 3, 2, 2))
    return universe


def point_matrix_snapshot(universe: PointUniverse) -> tuple[tuple[PointState, ...], ...]:
    return tuple(
        tuple(universe.point_state((x, y)) for x in range(universe.width))
        for y in range(universe.height)
    )


def expanded_path(points: tuple[Point, ...]) -> list[Point]:
    expanded: list[Point] = []
    for start, end in zip(points, points[1:]):
        x1, y1 = start
        x2, y2 = end
        if x1 != x2 and y1 != y2:
            raise AssertionError("candidate segment is not orthogonal")
        if x1 == x2:
            step = 1 if y2 >= y1 else -1
            segment = [(x1, y) for y in range(y1, y2 + step, step)]
        else:
            step = 1 if x2 >= x1 else -1
            segment = [(x, y1) for x in range(x1, x2 + step, step)]
        if expanded:
            assert expanded[-1] == segment[0]
            expanded.extend(segment[1:])
        else:
            expanded.extend(segment)
    return expanded


def is_orthogonal_path(points: tuple[Point, ...]) -> bool:
    return all(a[0] == b[0] or a[1] == b[1] for a, b in zip(points, points[1:]))


def has_square_neighbor(universe: PointUniverse, point: Point) -> bool:
    x, y = point
    for neighbor in ((x, y - 1), (x, y + 1), (x - 1, y), (x + 1, y)):
        if universe.in_bounds(neighbor) and universe.point_state(neighbor) & PointState.SQUARE:
            return True
    return False


def is_outer_edge_only(universe: PointUniverse, path: list[Point]) -> bool:
    return all(
        x == 0 or y == 0 or x == universe.width - 1 or y == universe.height - 1
        for x, y in path
    )


def assert_path_is_uncommitted_white(universe: PointUniverse, path: list[Point]) -> None:
    for point in path:
        assert universe.point_state(point) & PointState.WHITE


def hypothetical_region_delta(candidate: tuple[Point, ...]) -> tuple[int, int, int]:
    before = obstacle_universe()
    r_before = WhiteSpaceViewer(before).white_region_count()

    after = obstacle_universe()
    after.add_line(Line("CANDIDATE", "N", candidate))
    r_after = WhiteSpaceViewer(after).white_region_count()

    return r_before, r_after, r_after - r_before


def test_same_connection_candidate_paths_have_different_white_space_consequences() -> None:
    start = obstacle_universe()
    start_snapshot = point_matrix_snapshot(start)
    viewer = WhiteSpaceViewer(start)

    a = (2, 3)
    b = (2, 5)
    boundary_candidate = (a, (2, 4), b)
    interior_candidate = (a, (1, 3), (1, 5), b)

    assert boundary_candidate[0] == a
    assert boundary_candidate[-1] == b
    assert interior_candidate[0] == a
    assert interior_candidate[-1] == b
    assert boundary_candidate != interior_candidate
    assert is_orthogonal_path(boundary_candidate)
    assert is_orthogonal_path(interior_candidate)

    boundary_points = expanded_path(boundary_candidate)
    interior_points = expanded_path(interior_candidate)
    assert_path_is_uncommitted_white(start, boundary_points)
    assert_path_is_uncommitted_white(start, interior_points)

    region_id = viewer.white_region_at(*a)
    assert viewer.white_region_at(*b) == region_id
    region_boundary = viewer.white_region_boundary(region_id)
    assert all(point in region_boundary for point in boundary_points)
    assert any(has_square_neighbor(start, point) for point in boundary_points)
    assert not is_outer_edge_only(start, boundary_points)

    boundary_before, boundary_after, boundary_delta = hypothetical_region_delta(boundary_candidate)
    interior_before, interior_after, interior_delta = hypothetical_region_delta(interior_candidate)

    assert boundary_before == 1
    assert boundary_after == 1
    assert boundary_delta == 0
    assert interior_before == 1
    assert interior_after > interior_before
    assert interior_delta > 0

    assert point_matrix_snapshot(start) == start_snapshot
