"""Perceive whether a future connection remains routable after a candidate.

This module does not select routes.  It evaluates one already-generated
candidate in an independent hypothetical PointUniverse copy, commits that
candidate in the copy using normal PointUniverse line legality, and asks the
existing PointRouter whether a future connection has candidates there.
"""

from __future__ import annotations

from copy import deepcopy
from dataclasses import dataclass
from typing import Tuple

from .router import PointRouter, RouteCandidate
from .universe import Line, Point, PointUniverse


@dataclass(frozen=True)
class FutureRoutabilityResult:
    B_routable: bool
    B_candidate_count: int
    B_candidates: Tuple[RouteCandidate, ...]


def future_routability_after_candidate(
    universe: PointUniverse,
    candidate: RouteCandidate,
    *,
    B_start: Point,
    B_goal: Point,
    B_max_paths: int,
    hypothetical_line_id: str = "HYPOTHETICAL_A",
    hypothetical_net_id: str = "HYPOTHETICAL_NET_A",
) -> FutureRoutabilityResult:
    hypothetical = deepcopy(universe)
    hypothetical.add_line(Line(hypothetical_line_id, hypothetical_net_id, candidate.points))
    b_candidates = tuple(PointRouter(hypothetical).candidate_paths(B_start, B_goal, max_paths=B_max_paths))
    return FutureRoutabilityResult(
        B_routable=bool(b_candidates),
        B_candidate_count=len(b_candidates),
        B_candidates=b_candidates,
    )
