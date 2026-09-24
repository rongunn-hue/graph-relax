"""Source-only tests for the experiment-only structure-guided router."""

from __future__ import annotations

from pathlib import Path

from .router import PointRouter
from .structure_guided_router import HorizontalRunRouter
from .universe import Line, Point, PointState, PointUniverse, four_neighbors


def open_space_fixture() -> tuple[PointUniverse, Point, Point]:
    return PointUniverse(9, 7, line_separation=0, square_clearance=0), (1, 3), (7, 3)


def single_detour_fixture() -> tuple[PointUniverse, Point, Point]:
    universe = PointUniverse(9, 7, line_separation=0, square_clearance=0)
    universe.add_line(Line("WALL", "WALL", ((4, 1), (4, 2), (4, 3), (4, 4), (4, 5))))
    return universe, (1, 3), (7, 3)


def corridor_chamber_fixture() -> tuple[PointUniverse, Point, Point]:
    universe = PointUniverse(11, 9, line_separation=0, square_clearance=0)
    universe.add_line(Line("TOP_LEFT", "WALL", ((2, 2), (3, 2), (4, 2))))
    universe.add_line(Line("TOP_RIGHT", "WALL", ((6, 2), (7, 2), (8, 2))))
    universe.add_line(Line("BOTTOM_LEFT", "WALL", ((2, 6), (3, 6), (4, 6))))
    universe.add_line(Line("BOTTOM_RIGHT", "WALL", ((6, 6), (7, 6), (8, 6))))
    universe.add_line(Line("LEFT_INNER", "WALL", ((4, 3), (4, 4), (4, 5))))
    universe.add_line(Line("RIGHT_INNER", "WALL", ((6, 3), (6, 4), (6, 5))))
    return universe, (1, 4), (9, 4)


FIXTURES = (
    ("open", open_space_fixture),
    ("detour", single_detour_fixture),
    ("corridor", corridor_chamber_fixture),
)


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


def assert_legal_path(universe: PointUniverse, path: tuple[Point, ...], start: Point, goal: Point) -> None:
    assert path[0] == start
    assert path[-1] == goal
    assert len(path) == len(set(path))
    for point in path:
        assert universe.point_state(point) & PointState.WHITE
    for a, b in zip(path, path[1:]):
        assert b in four_neighbors(a)


def assert_cells_are_maximal_horizontal_runs(universe: PointUniverse, cells) -> None:
    white = set(universe.layer(PointState.WHITE))
    covered = set()
    expected_id = 1
    for cell in cells:
        assert cell.cell_id == expected_id
        expected_id += 1
        assert cell.points == tuple((x, cell.y) for x in range(cell.x_start, cell.x_end + 1))
        assert all(point in white for point in cell.points)
        assert (cell.x_start - 1, cell.y) not in white
        assert (cell.x_end + 1, cell.y) not in white
        covered.update(cell.points)
    assert covered == white


def assert_adjacency_matches_vertical_overlap(cells, adjacency) -> None:
    by_id = {cell.cell_id: cell for cell in cells}
    observed = {(a, b): overlap for a, b, overlap in adjacency}
    for a, b, overlap in adjacency:
        left = by_id[a]
        right = by_id[b]
        assert abs(left.y - right.y) == 1
        expected = tuple(range(max(left.x_start, right.x_start), min(left.x_end, right.x_end) + 1))
        assert overlap == expected
    for a in cells:
        for b in cells:
            if a.cell_id >= b.cell_id or abs(a.y - b.y) != 1:
                continue
            start = max(a.x_start, b.x_start)
            end = min(a.x_end, b.x_end)
            if start <= end:
                assert observed[(a.cell_id, b.cell_id)] == tuple(range(start, end + 1))


def test_structure_guided_router_locates_one_legal_path_in_fixed_fixtures_without_mutation() -> None:
    for _name, fixture in FIXTURES:
        universe, start, goal = fixture()
        before = snapshot(universe)
        baseline = PointRouter(universe).candidate_paths(start, goal, max_paths=1)
        result = HorizontalRunRouter(universe).locate_path(start, goal)

        assert baseline
        assert baseline[0].points[0] == start
        assert baseline[0].points[-1] == goal
        assert result.route_found is True
        assert result.complete_point_paths_enumerated == 1
        assert result.start_cell_id is not None
        assert result.goal_cell_id is not None
        assert result.cell_path[0] == result.start_cell_id
        assert result.cell_path[-1] == result.goal_cell_id
        assert result.dequeued_cell_count >= 1
        assert_legal_path(universe, result.point_path, start, goal)
        assert_cells_are_maximal_horizontal_runs(universe, result.cells)
        assert_adjacency_matches_vertical_overlap(result.cells, result.adjacency)
        assert snapshot(universe) == before


def test_structure_guided_router_uses_exactly_three_fixed_fixtures() -> None:
    assert tuple(name for name, _fixture in FIXTURES) == ("open", "detour", "corridor")


def test_structure_guided_router_source_has_no_disallowed_dependencies_or_search() -> None:
    source = Path(__file__).with_name("structure_guided_router.py").read_text()

    assert "PointRouter" not in source
    assert "hypothetical_terminal_capacity_consequence" not in source
    assert "hypothetical_multi_terminal_capacity_consequence" not in source
    assert "aggregate_terminal_capacity_consequence" not in source
    assert "RouteSelector" not in source
    assert "CapacityAwareSelector" not in source
    assert "candidate_paths" not in source
    assert "fixture" not in source
