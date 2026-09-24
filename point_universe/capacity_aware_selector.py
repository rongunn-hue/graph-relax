"""Experiment-only capacity-aware selection among existing route candidates.

This module does not generate routes, commit routes, or modify established
selection policy.  It only selects from candidates supplied by PointRouter for
the current controlled experiment.
"""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass
from typing import Tuple

from .router import RouteCandidate
from .universe import (
    AggregateTerminalCapacityConsequence,
    MultiTerminalCapacityConsequence,
    Point,
    TerminalPair,
    WhiteSpaceViewer,
)


@dataclass(frozen=True)
class CapacityAwareSelection:
    candidate: RouteCandidate
    multi_consequence: MultiTerminalCapacityConsequence | None
    aggregate: AggregateTerminalCapacityConsequence | None
    key: Tuple[object, ...]


class CapacityAwareSelector:
    """Select candidates by the fixed experiment-only capacity-aware key."""

    def select(
        self,
        candidates: Sequence[RouteCandidate],
        *,
        viewer: WhiteSpaceViewer,
        pending_pairs: Sequence[TerminalPair],
    ) -> CapacityAwareSelection:
        if not candidates:
            raise ValueError("cannot select from empty candidate collection")
        if not pending_pairs:
            selected = min(
                candidates,
                key=lambda candidate: (
                    candidate.consequence.delta_R,
                    len(candidate.points),
                    candidate.points,
                ),
            )
            return CapacityAwareSelection(
                candidate=selected,
                multi_consequence=None,
                aggregate=None,
                key=(
                    selected.consequence.delta_R,
                    len(selected.points),
                    selected.points,
                ),
            )

        selections = []
        for candidate in candidates:
            multi = viewer.hypothetical_multi_terminal_capacity_consequence(candidate.points, pending_pairs)
            aggregate = viewer.aggregate_terminal_capacity_consequence(multi)
            key = (
                aggregate.disconnected_after_count,
                -aggregate.total_delta_capacity,
                candidate.consequence.delta_R,
                len(candidate.points),
                candidate.points,
            )
            selections.append(
                CapacityAwareSelection(
                    candidate=candidate,
                    multi_consequence=multi,
                    aggregate=aggregate,
                    key=key,
                )
            )
        return min(selections, key=lambda selection: selection.key)


def fixed_capacity_aware_key(
    candidate: RouteCandidate,
    aggregate: AggregateTerminalCapacityConsequence,
) -> Tuple[object, ...]:
    return (
        aggregate.disconnected_after_count,
        -aggregate.total_delta_capacity,
        candidate.consequence.delta_R,
        len(candidate.points),
        candidate.points,
    )
