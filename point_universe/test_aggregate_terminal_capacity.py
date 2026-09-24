"""Source-only tests for aggregate terminal-capacity consequence representation."""

from __future__ import annotations

from dataclasses import replace

import pytest

from .test_multi_terminal_capacity import (
    CANDIDATE_A,
    PAIR_B,
    PAIR_C,
    PAIR_D,
    PRIMARY_PAIRS,
    multi_capacity_universe,
    snapshot,
)
from .universe import (
    AggregateTerminalCapacityConsequence,
    HypotheticalTerminalCapacityConsequence,
    MultiTerminalCapacityConsequence,
    TerminalPair,
    WhiteSpaceViewer,
)


def valid_consequence(
    *,
    before: int = 2,
    after: int = 2,
    delta: int = 0,
    routable: bool = True,
    finite: bool = True,
    before_cuts: tuple[tuple[tuple[int, int], ...], ...] = (((1, 0),),),
    after_cuts: tuple[tuple[tuple[int, int], ...], ...] = (((1, 0),),),
) -> HypotheticalTerminalCapacityConsequence:
    return HypotheticalTerminalCapacityConsequence(
        capacity_before=before,
        capacity_after=after,
        delta_capacity=delta,
        B_routable_after=routable,
        finite_cut_after=finite,
        minimum_cuts_before=before_cuts,
        minimum_cuts_after=after_cuts,
    )


def test_aggregate_terminal_capacity_summarizes_established_multi_pair_case_without_mutation() -> None:
    universe = multi_capacity_universe()
    viewer = WhiteSpaceViewer(universe)
    before = snapshot(universe)
    multi = viewer.hypothetical_multi_terminal_capacity_consequence(CANDIDATE_A, PRIMARY_PAIRS)
    multi_before = multi

    aggregate = viewer.aggregate_terminal_capacity_consequence(multi)

    assert aggregate == AggregateTerminalCapacityConsequence(
        pair_count=3,
        routable_after_count=2,
        disconnected_after_count=1,
        unchanged_capacity_pair_ids=("C",),
        reduced_positive_capacity_pair_ids=("B",),
        disconnected_pair_ids=("D",),
        total_capacity_before=8,
        total_capacity_after=4,
        total_delta_capacity=-4,
    )
    assert aggregate.pair_count == aggregate.routable_after_count + aggregate.disconnected_after_count
    assert aggregate.pair_count == (
        len(aggregate.unchanged_capacity_pair_ids)
        + len(aggregate.reduced_positive_capacity_pair_ids)
        + len(aggregate.disconnected_pair_ids)
    )
    assert aggregate.total_delta_capacity == aggregate.total_capacity_after - aggregate.total_capacity_before
    assert snapshot(universe) == before
    assert multi == multi_before


def test_aggregate_terminal_capacity_preserves_pair_order_classes() -> None:
    universe = multi_capacity_universe()
    viewer = WhiteSpaceViewer(universe)
    reordered = viewer.hypothetical_multi_terminal_capacity_consequence(CANDIDATE_A, (PAIR_D, PAIR_B, PAIR_C))

    aggregate = viewer.aggregate_terminal_capacity_consequence(reordered)

    assert tuple(pair_id for pair_id, _ in reordered.pair_consequences) == ("D", "B", "C")
    assert aggregate == AggregateTerminalCapacityConsequence(
        pair_count=3,
        routable_after_count=2,
        disconnected_after_count=1,
        unchanged_capacity_pair_ids=("C",),
        reduced_positive_capacity_pair_ids=("B",),
        disconnected_pair_ids=("D",),
        total_capacity_before=8,
        total_capacity_after=4,
        total_delta_capacity=-4,
    )

    synthetic = MultiTerminalCapacityConsequence(
        candidate_points=((9, 9),),
        pair_consequences=(
            ("X2", valid_consequence()),
            ("X1", valid_consequence()),
        ),
    )
    synthetic_aggregate = viewer.aggregate_terminal_capacity_consequence(synthetic)

    assert synthetic_aggregate.unchanged_capacity_pair_ids == ("X2", "X1")
    assert synthetic_aggregate.pair_count == 2
    assert synthetic_aggregate.routable_after_count == 2
    assert synthetic_aggregate.disconnected_after_count == 0
    assert synthetic_aggregate.total_capacity_before == 4
    assert synthetic_aggregate.total_capacity_after == 4
    assert synthetic_aggregate.total_delta_capacity == 0


def test_aggregate_terminal_capacity_rejects_malformed_inputs() -> None:
    viewer = WhiteSpaceViewer(multi_capacity_universe())
    good = valid_consequence()

    malformed_cases = (
        ("wrong top-level type", object(), "consequence must be a MultiTerminalCapacityConsequence"),
        (
            "empty pair_consequences",
            MultiTerminalCapacityConsequence(candidate_points=(), pair_consequences=()),
            "pair_consequences must be non-empty",
        ),
        (
            "empty pair_id",
            MultiTerminalCapacityConsequence(candidate_points=(), pair_consequences=(("", good),)),
            "pair_id must be non-empty",
        ),
        (
            "non-string pair_id",
            MultiTerminalCapacityConsequence(candidate_points=(), pair_consequences=((3, good),)),
            "pair_id must be a string",
        ),
        (
            "duplicate pair_id",
            MultiTerminalCapacityConsequence(candidate_points=(), pair_consequences=(("X", good), ("X", good))),
            "duplicate pair_id",
        ),
        (
            "wrong individual consequence type",
            MultiTerminalCapacityConsequence(candidate_points=(), pair_consequences=(("X", object()),)),
            "pair consequence must be a HypotheticalTerminalCapacityConsequence",
        ),
        (
            "capacity_before == 0",
            MultiTerminalCapacityConsequence(candidate_points=(), pair_consequences=(("X", replace(good, capacity_before=0)),)),
            "capacity_before must be a positive integer",
        ),
        (
            "capacity_before < 0",
            MultiTerminalCapacityConsequence(candidate_points=(), pair_consequences=(("X", replace(good, capacity_before=-1)),)),
            "capacity_before must be a positive integer",
        ),
        (
            "capacity_after < 0",
            MultiTerminalCapacityConsequence(candidate_points=(), pair_consequences=(("X", replace(good, capacity_after=-1)),)),
            "capacity_after must be a non-negative integer",
        ),
        (
            "capacity_after is None",
            MultiTerminalCapacityConsequence(candidate_points=(), pair_consequences=(("X", replace(good, capacity_after=None)),)),
            "capacity_after must be an integer",
        ),
        (
            "delta_capacity is None",
            MultiTerminalCapacityConsequence(candidate_points=(), pair_consequences=(("X", replace(good, delta_capacity=None)),)),
            "delta_capacity must be an integer",
        ),
        (
            "delta_capacity inconsistent",
            MultiTerminalCapacityConsequence(candidate_points=(), pair_consequences=(("X", replace(good, delta_capacity=1)),)),
            "delta_capacity must equal capacity_after - capacity_before",
        ),
        (
            "disconnected but capacity_after != 0",
            MultiTerminalCapacityConsequence(
                candidate_points=(),
                pair_consequences=(("X", valid_consequence(after=1, delta=-1, routable=False, finite=False, after_cuts=())),),
            ),
            "disconnected consequence must have capacity_after == 0",
        ),
        (
            "disconnected but finite_cut_after == True",
            MultiTerminalCapacityConsequence(
                candidate_points=(),
                pair_consequences=(("X", valid_consequence(after=0, delta=-2, routable=False, finite=True, after_cuts=())),),
            ),
            "disconnected consequence must have finite_cut_after == False",
        ),
        (
            "disconnected but minimum_cuts_after non-empty",
            MultiTerminalCapacityConsequence(
                candidate_points=(),
                pair_consequences=(("X", valid_consequence(after=0, delta=-2, routable=False, finite=False)),),
            ),
            "disconnected consequence must have empty minimum_cuts_after",
        ),
        (
            "routable but capacity_after == 0",
            MultiTerminalCapacityConsequence(
                candidate_points=(),
                pair_consequences=(("X", valid_consequence(after=0, delta=-2)),),
            ),
            "routable consequence must have positive capacity_after",
        ),
        (
            "routable but finite_cut_after == False",
            MultiTerminalCapacityConsequence(
                candidate_points=(),
                pair_consequences=(("X", valid_consequence(finite=False)),),
            ),
            "routable consequence must have finite_cut_after == True",
        ),
        (
            "routable but minimum_cuts_after == ()",
            MultiTerminalCapacityConsequence(
                candidate_points=(),
                pair_consequences=(("X", valid_consequence(after_cuts=())),),
            ),
            "routable consequence must have non-empty minimum_cuts_after",
        ),
        (
            "empty minimum_cuts_before",
            MultiTerminalCapacityConsequence(
                candidate_points=(),
                pair_consequences=(("X", valid_consequence(before_cuts=())),),
            ),
            "minimum_cuts_before must be non-empty",
        ),
        (
            "routable delta_capacity > 0",
            MultiTerminalCapacityConsequence(
                candidate_points=(),
                pair_consequences=(("X", valid_consequence(before=2, after=3, delta=1)),),
            ),
            "routable capacity increase is not an accepted aggregate class",
        ),
    )

    for _label, malformed, message in malformed_cases:
        with pytest.raises(ValueError, match=message):
            viewer.aggregate_terminal_capacity_consequence(malformed)


def test_aggregate_terminal_capacity_is_pure_with_respect_to_universe_geometry() -> None:
    universe = multi_capacity_universe()
    viewer = WhiteSpaceViewer(universe)
    before = snapshot(universe)
    multi = viewer.hypothetical_multi_terminal_capacity_consequence(CANDIDATE_A, PRIMARY_PAIRS)
    multi_before = multi

    viewer.aggregate_terminal_capacity_consequence(multi)

    assert snapshot(universe) == before
    assert multi == multi_before
    assert universe.lines == dict(universe.lines)
