"""Small deterministic graph input for the point-universe prototype.

This file defines data only. It does not execute the example.
"""

from __future__ import annotations

from .universe import Line, LineEndpoint, PointUniverse, Square


def build_small_example() -> PointUniverse:
    universe = PointUniverse(width=34, height=22, line_separation=1, square_clearance=1)

    amp = Square("AMP", x=3, y=4, width=4, height=4)
    amp.add_attachment("AMP_TOP", "top", 2)
    amp.add_attachment("AMP_RIGHT", "right", 2)
    amp.add_attachment("AMP_BOTTOM", "bottom", 2)

    sensor = Square("SENSOR", x=14, y=3, width=2, height=2)
    sensor.add_attachment("SENSOR_LEFT", "left", 1)
    sensor.add_attachment("SENSOR_RIGHT", "right", 1)

    filter_sq = Square("FILTER", x=23, y=4, width=3, height=3)
    filter_sq.add_attachment("FILTER_LEFT", "left", 1)
    filter_sq.add_attachment("FILTER_BOTTOM", "bottom", 2)

    bias = Square("BIAS", x=6, y=14, width=2, height=2)
    bias.add_attachment("BIAS_TOP", "top", 1)
    bias.add_attachment("BIAS_RIGHT", "right", 1)

    output = Square("OUTPUT", x=18, y=14, width=4, height=2)
    output.add_attachment("OUTPUT_LEFT", "left", 1)
    output.add_attachment("OUTPUT_RIGHT", "right", 1)

    trim = Square("TRIM", x=28, y=12, width=2, height=4, rotation=90)
    trim.add_attachment("TRIM_LEFT", "left", 2)
    trim.add_attachment("TRIM_TOP", "top", 1)
    trim.add_attachment("TRIM_RIGHT", "right", 2)

    for square in (amp, sensor, filter_sq, bias, output, trim):
        universe.add_square(square)

    universe.add_line(Line(
        "L_AMP_SENSOR",
        "NET_A",
        points=((7, 6), (11, 6), (11, 4), (14, 4)),
        endpoints=(LineEndpoint(attachment_id="AMP_RIGHT"), LineEndpoint(attachment_id="SENSOR_LEFT")),
    ))
    universe.add_line(Line(
        "L_SENSOR_FILTER",
        "NET_B",
        points=((16, 4), (20, 4), (20, 5), (23, 5)),
        endpoints=(LineEndpoint(attachment_id="SENSOR_RIGHT"), LineEndpoint(attachment_id="FILTER_LEFT")),
    ))
    universe.add_line(Line(
        "L_AMP_BIAS",
        "NET_C",
        points=((5, 8), (5, 12), (7, 12), (7, 14)),
        endpoints=(LineEndpoint(attachment_id="AMP_BOTTOM"), LineEndpoint(attachment_id="BIAS_TOP")),
    ))
    universe.add_line(Line(
        "L_BIAS_OUTPUT",
        "NET_D",
        points=((8, 15), (12, 15), (18, 15)),
        endpoints=(LineEndpoint(attachment_id="BIAS_RIGHT"), LineEndpoint(attachment_id="OUTPUT_LEFT")),
    ))
    universe.add_line(Line(
        "L_FILTER_JOINED_SEGMENTS",
        "NET_E",
        points=((25, 7), (25, 10), (30, 10), (30, 12)),
        endpoints=(LineEndpoint(attachment_id="FILTER_BOTTOM"), LineEndpoint(attachment_id="TRIM_LEFT")),
    ))
    universe.add_line(Line(
        "L_OUTPUT_TRIM",
        "NET_F",
        points=((22, 15), (26, 15), (30, 15), (30, 14)),
        endpoints=(LineEndpoint(attachment_id="OUTPUT_RIGHT"), LineEndpoint(attachment_id="TRIM_RIGHT")),
    ))

    return universe
