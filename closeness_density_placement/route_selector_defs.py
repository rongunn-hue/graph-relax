"""Deterministic selection among already-generated route candidates."""

from __future__ import annotations

from collections.abc import Sequence

from router_defs import RouteCandidate


class RouteSelector:
    def select(self, candidates: Sequence[RouteCandidate]) -> RouteCandidate:
        if not candidates:
            raise ValueError("cannot select from empty candidate collection")
        return min(
            candidates,
            key=lambda candidate: (
                candidate.consequence.delta_R,
                len(candidate.points),
                candidate.points,
            ),
        )
