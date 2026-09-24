"""Unexecuted unit tests for point-universe white-space perception."""

from __future__ import annotations

from .universe import Line, PointState, PointUniverse, WhiteSpaceViewer


def test_horizontal_line_divides_11_by_21_white_region() -> None:
    universe = PointUniverse(width=11, height=21, line_separation=0, square_clearance=0)

    initial = WhiteSpaceViewer(universe)
    assert initial.white_region_count() == 1
    assert initial.white_region_size(1) == 231

    universe.add_line(Line(
        "MIDDLE_ROW",
        "CUT",
        points=((0, 10), (10, 10)),
    ))

    divided = WhiteSpaceViewer(universe)
    sizes = sorted(divided.white_region_size(region_id) for region_id in range(1, divided.white_region_count() + 1))

    assert divided.white_region_count() == 2
    assert sizes == [110, 110]
    assert divided.same_white_region((0, 0), (10, 9))
    assert divided.same_white_region((0, 11), (10, 20))
    assert not divided.same_white_region((0, 0), (0, 11))
    assert not (universe.point_state((0, 10)) & PointState.WHITE)
