"""Unexecuted unit-test source for corrected point-universe laws."""

from __future__ import annotations

import pytest

from .universe import Line, LineEndpoint, PointState, PointUniverse, Square, WhiteSpaceViewer


def test_minimum_2x2_square_and_one_attachment_per_side() -> None:
    square = Square("S", 2, 2, 2, 2)
    square.add_attachment("S_TOP", "top", 1)
    square.add_attachment("S_RIGHT", "right", 1)
    square.add_attachment("S_BOTTOM", "bottom", 1)
    square.add_attachment("S_LEFT", "left", 1)
    universe = PointUniverse(8, 8)
    universe.add_square(square)
    assert universe.point_state((3, 2)) & PointState.ATTACHMENT
    assert universe.point_state((4, 3)) & PointState.ATTACHMENT
    assert universe.point_state((3, 4)) & PointState.ATTACHMENT
    assert universe.point_state((2, 3)) & PointState.ATTACHMENT


def test_vertex_and_illegal_attachment_offsets_rejected() -> None:
    square = Square("S", 0, 0, 2, 2)
    with pytest.raises(ValueError):
        square.add_attachment("VERTEX", "top", 0)
    with pytest.raises(ValueError):
        square.add_attachment("OTHER_VERTEX", "right", 2)


def test_90_degree_rotation_of_non_square_swaps_footprint_and_attachment_boundary() -> None:
    square = Square("R", 5, 5, 2, 4, rotation=90)
    square.add_attachment("R_LEFT", "left", 2)
    universe = PointUniverse(20, 20)
    universe.add_square(square)
    assert (9, 7) in universe.layer(PointState.SQUARE)
    assert (7, 5) == universe.attachment_points["R_LEFT"]
    assert universe.point_state((7, 5)) & PointState.ATTACHMENT


def test_line_through_square_interior_rejected() -> None:
    square = Square("S", 4, 4, 2, 2)
    universe = PointUniverse(12, 12)
    universe.add_square(square)
    with pytest.raises(ValueError):
        universe.add_line(Line("L", "N", ((0, 5), (10, 5))))


def test_line_endpoint_not_coincident_with_declared_attachment_rejected() -> None:
    square = Square("S", 4, 4, 2, 2)
    square.add_attachment("S_LEFT", "left", 1)
    universe = PointUniverse(12, 12)
    universe.add_square(square)
    with pytest.raises(ValueError):
        universe.add_line(Line(
            "L",
            "N",
            ((0, 5), (5, 5)),
            endpoints=(LineEndpoint(), LineEndpoint(attachment_id="S_LEFT")),
        ))


def test_legal_ordinary_bend_accepted() -> None:
    universe = PointUniverse(12, 12, line_separation=0)
    universe.add_line(Line("BEND", "N", ((1, 1), (6, 1), (6, 6))))
    assert universe.point_state((6, 1)) & PointState.LINE


def test_illegal_self_crossing_rejected() -> None:
    universe = PointUniverse(12, 12, line_separation=0)
    with pytest.raises(ValueError):
        universe.add_line(Line("SELF_CROSS", "N", ((1, 1), (7, 1), (7, 6), (3, 6), (3, 0))))


def test_illegal_self_retrace_rejected() -> None:
    universe = PointUniverse(12, 12, line_separation=0)
    with pytest.raises(ValueError):
        universe.add_line(Line("SELF_RETRACE", "N", ((1, 1), (7, 1), (7, 4), (4, 4), (4, 1))))


def test_unrelated_line_overlap_rejected() -> None:
    universe = PointUniverse(12, 12, line_separation=0)
    universe.add_line(Line("A", "NA", ((1, 1), (5, 1))))
    with pytest.raises(ValueError):
        universe.add_line(Line("B", "NB", ((3, 1), (7, 1))))


def test_line_spacing_violation_rejected() -> None:
    universe = PointUniverse(12, 12, line_separation=1)
    universe.add_line(Line("A", "NA", ((1, 1), (5, 1))))
    with pytest.raises(ValueError):
        universe.add_line(Line("B", "NB", ((1, 2), (5, 2))))


def test_legal_connected_line_junction_accepted() -> None:
    universe = PointUniverse(12, 12, line_separation=0)
    universe.add_line(Line("A", "N", ((1, 5), (5, 5))))
    universe.add_line(Line(
        "B",
        "N",
        ((5, 5), (5, 8)),
        endpoints=(LineEndpoint(line_id="A"), LineEndpoint()),
    ))
    assert universe.connectivity.point_to_lines[(5, 5)] == {"A", "B"}


def test_more_than_two_lines_at_one_point_rejected() -> None:
    universe = PointUniverse(12, 12, line_separation=0)
    universe.add_line(Line("A", "N", ((1, 5), (5, 5))))
    universe.add_line(Line("B", "N", ((5, 5), (5, 8)), endpoints=(LineEndpoint(line_id="A"), LineEndpoint())))
    with pytest.raises(ValueError):
        universe.add_line(Line("C", "N", ((5, 5), (8, 5)), endpoints=(LineEndpoint(line_id="A"), LineEndpoint())))


def test_connectivity_attachment_endpoint_not_line_endpoint() -> None:
    square = Square("S", 4, 4, 2, 2)
    square.add_attachment("S_LEFT", "left", 1)
    universe = PointUniverse(12, 12, line_separation=0)
    universe.add_square(square)
    universe.add_line(Line(
        "ATTACHED",
        "N",
        ((1, 5), (4, 5)),
        endpoints=(LineEndpoint(), LineEndpoint(attachment_id="S_LEFT")),
    ))
    assert universe.connectivity.line_to_attachment[("ATTACHED", 1)] == "S_LEFT"
    assert ("ATTACHED", 1) not in universe.connectivity.line_to_line


def test_connectivity_line_endpoint_not_attachment_endpoint() -> None:
    universe = PointUniverse(12, 12, line_separation=0)
    universe.add_line(Line("A", "N", ((1, 5), (5, 5))))
    universe.add_line(Line(
        "B",
        "N",
        ((5, 5), (5, 8)),
        endpoints=(LineEndpoint(line_id="A"), LineEndpoint()),
    ))
    assert universe.connectivity.line_to_line[("B", 0)] == "A"
    assert ("B", 0) not in universe.connectivity.line_to_attachment


def test_white_unusable_classification_comes_from_matrix_state() -> None:
    universe = PointUniverse(7, 7, square_clearance=1)
    square = Square("S", 2, 2, 2, 2)
    universe.add_square(square)
    assert universe.point_state((2, 2)) & PointState.SQUARE
    assert universe.point_state((1, 2)) & PointState.UNUSABLE
    assert universe.point_state((0, 0)) & PointState.WHITE


def test_viewer_uses_matrix_white_space_state() -> None:
    universe = PointUniverse(5, 3, line_separation=0)
    universe.add_line(Line("CUT", "N", ((2, 0), (2, 2))))
    viewer = WhiteSpaceViewer(universe)
    assert viewer.white_region_count() == 2
    assert not viewer.same_white_region((0, 1), (4, 1))


def test_11x21_split_specification_remains_present() -> None:
    universe = PointUniverse(width=11, height=21, line_separation=0, square_clearance=0)
    assert WhiteSpaceViewer(universe).white_region_size(1) == 231
    universe.add_line(Line("MIDDLE_ROW", "CUT", ((0, 10), (10, 10))))
    sizes = sorted(WhiteSpaceViewer(universe).white_region_size(region_id) for region_id in (1, 2))
    assert sizes == [110, 110]
