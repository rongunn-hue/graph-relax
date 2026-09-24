"""Experiment-only horizontal-run cell-decomposition router.

This router locates one path through WHITE space by constructing a compact
cell graph from maximal horizontal WHITE runs. It does not use the existing
point-path enumerator, score alternatives, or mutate the authoritative universe.
"""

from __future__ import annotations

from collections import deque
from dataclasses import dataclass
from typing import Dict, FrozenSet, Iterable, Optional, Tuple

from .universe import Point, PointState, PointUniverse


@dataclass(frozen=True)
class HorizontalCell:
    cell_id: int
    y: int
    x_start: int
    x_end: int
    points: Tuple[Point, ...]


@dataclass(frozen=True)
class StructureGuidedRoute:
    route_found: bool
    white_point_count: int
    cells: Tuple[HorizontalCell, ...]
    adjacency: Tuple[Tuple[int, int, Tuple[int, ...]], ...]
    start_cell_id: Optional[int]
    goal_cell_id: Optional[int]
    dequeued_cell_count: int
    cell_path: Tuple[int, ...]
    point_path: Tuple[Point, ...]
    complete_point_paths_enumerated: int


class HorizontalRunRouter:
    def __init__(self, universe: PointUniverse):
        self.universe = universe

    def locate_path(self, start: Point, goal: Point) -> StructureGuidedRoute:
        self._require_white(start, "start")
        self._require_white(goal, "goal")
        white_points = set(self.universe.layer(PointState.WHITE))
        cells = self._build_cells(white_points)
        point_to_cell = {
            point: cell.cell_id
            for cell in cells
            for point in cell.points
        }
        adjacency_map = self._build_adjacency(cells)
        adjacency = tuple(
            (min(a, b), max(a, b), overlap)
            for (a, b), overlap in sorted(adjacency_map.items())
            if a < b
        )
        start_cell = point_to_cell[start]
        goal_cell = point_to_cell[goal]
        cell_path, dequeued = self._bfs_cell_path(start_cell, goal_cell, adjacency_map)
        if not cell_path:
            return StructureGuidedRoute(
                route_found=False,
                white_point_count=len(white_points),
                cells=cells,
                adjacency=adjacency,
                start_cell_id=start_cell,
                goal_cell_id=goal_cell,
                dequeued_cell_count=dequeued,
                cell_path=(),
                point_path=(),
                complete_point_paths_enumerated=0,
            )
        point_path = self._reconstruct_point_path(start, goal, cell_path, {cell.cell_id: cell for cell in cells}, adjacency_map)
        return StructureGuidedRoute(
            route_found=True,
            white_point_count=len(white_points),
            cells=cells,
            adjacency=adjacency,
            start_cell_id=start_cell,
            goal_cell_id=goal_cell,
            dequeued_cell_count=dequeued,
            cell_path=cell_path,
            point_path=point_path,
            complete_point_paths_enumerated=1,
        )

    def _require_white(self, point: Point, label: str) -> None:
        if not self.universe.in_bounds(point):
            raise ValueError(f"{label} point outside universe: {point}")
        if not (self.universe.point_state(point) & PointState.WHITE):
            raise ValueError(f"{label} point is not WHITE: {point}")

    def _build_cells(self, white_points: set[Point]) -> Tuple[HorizontalCell, ...]:
        cells = []
        cell_id = 1
        for y in range(self.universe.height):
            x = 0
            while x < self.universe.width:
                if (x, y) not in white_points:
                    x += 1
                    continue
                x_start = x
                while x + 1 < self.universe.width and (x + 1, y) in white_points:
                    x += 1
                x_end = x
                points = tuple((px, y) for px in range(x_start, x_end + 1))
                cells.append(HorizontalCell(cell_id, y, x_start, x_end, points))
                cell_id += 1
                x += 1
        return tuple(cells)

    def _build_adjacency(self, cells: Tuple[HorizontalCell, ...]) -> Dict[Tuple[int, int], Tuple[int, ...]]:
        by_row: Dict[int, Tuple[HorizontalCell, ...]] = {}
        for cell in cells:
            by_row.setdefault(cell.y, ())
            by_row[cell.y] = by_row[cell.y] + (cell,)
        adjacency: Dict[Tuple[int, int], Tuple[int, ...]] = {}
        for y in range(self.universe.height - 1):
            for lower in by_row.get(y, ()):
                for upper in by_row.get(y + 1, ()):
                    start = max(lower.x_start, upper.x_start)
                    end = min(lower.x_end, upper.x_end)
                    if start <= end:
                        overlap = tuple(range(start, end + 1))
                        adjacency[(lower.cell_id, upper.cell_id)] = overlap
                        adjacency[(upper.cell_id, lower.cell_id)] = overlap
        return adjacency

    def _bfs_cell_path(
        self,
        start_cell: int,
        goal_cell: int,
        adjacency: Dict[Tuple[int, int], Tuple[int, ...]],
    ) -> Tuple[Tuple[int, ...], int]:
        queue = deque([start_cell])
        parent: Dict[int, Optional[int]] = {start_cell: None}
        dequeued = 0
        while queue:
            cell = queue.popleft()
            dequeued += 1
            if cell == goal_cell:
                break
            neighbors = sorted(b for a, b in adjacency if a == cell)
            for neighbor in neighbors:
                if neighbor not in parent:
                    parent[neighbor] = cell
                    queue.append(neighbor)
        if goal_cell not in parent:
            return (), dequeued
        path = []
        cursor: Optional[int] = goal_cell
        while cursor is not None:
            path.append(cursor)
            cursor = parent[cursor]
        path.reverse()
        return tuple(path), dequeued

    def _reconstruct_point_path(
        self,
        start: Point,
        goal: Point,
        cell_path: Tuple[int, ...],
        cells: Dict[int, HorizontalCell],
        adjacency: Dict[Tuple[int, int], Tuple[int, ...]],
    ) -> Tuple[Point, ...]:
        route = [start]
        current = start
        for current_cell_id, next_cell_id in zip(cell_path, cell_path[1:]):
            overlap = adjacency[(current_cell_id, next_cell_id)]
            portal_x = min(overlap, key=lambda x: (abs(x - current[0]), x))
            current = self._append_horizontal(route, current, portal_x)
            next_y = cells[next_cell_id].y
            step = (current[0], next_y)
            if step != current:
                route.append(step)
                current = step
        current = self._append_horizontal(route, current, goal[0])
        if current != goal:
            raise ValueError("deterministic reconstruction did not end at goal")
        if len(set(route)) != len(route):
            raise ValueError("deterministic reconstruction repeated a point")
        return tuple(route)

    def _append_horizontal(self, route: list[Point], current: Point, target_x: int) -> Point:
        x, y = current
        if x == target_x:
            return current
        step = 1 if target_x > x else -1
        while x != target_x:
            x += step
            point = (x, y)
            if point != route[-1]:
                route.append(point)
        return (target_x, y)
