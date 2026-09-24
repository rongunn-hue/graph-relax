"""Commit an already-selected point-universe route candidate."""

from __future__ import annotations

from .router import RouteCandidate
from .universe import Line, PointUniverse


class RouteCommitter:
    def __init__(self, universe: PointUniverse):
        self.universe = universe

    def commit(self, candidate: RouteCandidate, *, line_id: str, net_id: str) -> Line:
        line = Line(line_id, net_id, candidate.points)
        self.universe.add_line(line)
        return line
