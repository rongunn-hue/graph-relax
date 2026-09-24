"""Small point-universe candidate path enumerator.

This module intentionally does not choose, rank, or commit routes.  It only
enumerates bounded simple WHITE-point paths and attaches WhiteSpaceViewer's
hypothetical consequence perception to each complete candidate.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import List, Tuple

from .universe import HypotheticalWhiteConsequence, Point, PointState, PointUniverse, WhiteSpaceViewer


NEIGHBOR_DELTAS: Tuple[Point, Point, Point, Point] = (
    (0, -1),   # UP
    (-1, 0),   # LEFT
    (1, 0),    # RIGHT
    (0, 1),    # DOWN
)


@dataclass(frozen=True)
class RouteCandidate:
    points: Tuple[Point, ...]
    consequence: HypotheticalWhiteConsequence


class PointRouter:
    def __init__(self, universe: PointUniverse):
        self.universe = universe

    def candidate_paths(self, start: Point, goal: Point, *, max_paths: int) -> List[RouteCandidate]:
        self._validate_request(start, goal, max_paths)
        viewer = WhiteSpaceViewer(self.universe)
        if start == goal:
            path = (start,)
            return [RouteCandidate(path, viewer.hypothetical_white_consequence(path))]

        candidates: List[RouteCandidate] = []
        visited = {start}
        stack: List[Point] = [start]

        def search(point: Point) -> None:
            if len(candidates) >= max_paths:
                return
            if point == goal:
                path = tuple(stack)
                candidates.append(RouteCandidate(path, viewer.hypothetical_white_consequence(path)))
                return
            for neighbor in self._neighbors(point):
                if neighbor in visited:
                    continue
                if not self._is_white(neighbor):
                    continue
                visited.add(neighbor)
                stack.append(neighbor)
                search(neighbor)
                stack.pop()
                visited.remove(neighbor)
                if len(candidates) >= max_paths:
                    return

        search(start)
        return candidates

    def _validate_request(self, start: Point, goal: Point, max_paths: int) -> None:
        if max_paths <= 0:
            raise ValueError("max_paths must be positive")
        self._require_white_point(start, "start")
        self._require_white_point(goal, "goal")

    def _require_white_point(self, point: Point, label: str) -> None:
        if not self.universe.in_bounds(point):
            raise ValueError(f"{label} point outside universe: {point}")
        if not self._is_white(point):
            raise ValueError(f"{label} point is not WHITE: {point}")

    def _is_white(self, point: Point) -> bool:
        return bool(self.universe.point_state(point) & PointState.WHITE)

    def _neighbors(self, point: Point) -> Tuple[Point, ...]:
        x, y = point
        return tuple((x + dx, y + dy) for dx, dy in NEIGHBOR_DELTAS if self.universe.in_bounds((x + dx, y + dy)))
