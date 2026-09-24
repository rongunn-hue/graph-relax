#!/usr/bin/env python3
"""Real schematic symbol artwork, reusing the existing production symbol library's DATA (not its drawing
functions, which assume one fixed aspect ratio per device and a single uniform scale factor - the opposite of
this Universe's rule that every device is an equal-sided square). Same normalized geometry, a transform of my
own that fits it non-uniformly into whatever square the device actually has, so pins land exactly on the real
routed attachment points regardless of aspect ratio.

Physical pin layout is NOT optional for these three device types - the symbol's internal wiring (which line goes
to which triangle vertex) is fixed to a PHYSICAL pin position (e.g. LM1458 pin 2 is always the left side's 2nd
pin from the top), not to whichever net it happens to carry. `forced_pinout()` below returns the forced
side/offset assignment for these types; `symbol_side_len()` returns the real package's required square side (an
8-pin DIP needs 4 pins on one side, which needs side >= 5, not the generic pin-count formula's 3).
"""
import sys
sys.path.insert(0, "/home/rgunn/picaxe/schematic_router")
from symbol_library import (
    resistor_symbol_points, capacitor_symbol_segments, diode_symbol_segments, zener_symbol_segments,
    SYMBOL_PIN_ORDER, SYMBOL_2PIN_TYPES,
    OPAMP_DUAL_PINOUT, DUAL_OPAMP_NORM_ORIGIN, DUAL_OPAMP_NORM_SIZE, DUAL_OPAMP_BODY, DUAL_OPAMP_NOTCH_BEZIER,
    DUAL_OPAMP_TERMINALS, DUAL_OPAMP_PIN_PATHS, DUAL_OPAMP_TRIANGLE_A, DUAL_OPAMP_TRIANGLE_B,
    OPAMP_DIP8_PINOUT, OPAMP_DIP8_NORM_ORIGIN, OPAMP_DIP8_NORM_SIZE, OPAMP_DIP8_BODY, OPAMP_DIP8_NOTCH_BEZIER,
    OPAMP_DIP8_TERMINALS, OPAMP_DIP8_LEAD_PATHS, OPAMP_DIP8_TRIANGLE,
    TRANSISTOR_PINOUT, NPN_NORM_ORIGIN, NPN_NORM_SIZE, NPN_TERMINALS, NPN_BODY, NPN_BASE_BAR, NPN_LEAD_PATHS,
    NPN_ARROWHEAD,
)

DEVICE_TYPES_NEEDING_FORCED_PINOUT = set(OPAMP_DUAL_PINOUT) | set(OPAMP_DIP8_PINOUT) | set(TRANSISTOR_PINOUT)


# Square side and offsets chosen to make the DRAWN wire and the symbol's OWN lead artwork meet at the same
# point, not just close - the Universe's own attachment offsets are evenly spaced integers (1..L-1), but the
# real symbol data's pin positions are NOT evenly spaced fractions of the package height, so the generic evenly-
# spaced formula leaves a small but real, visible gap between where the wire ends and where the symbol's lead
# is drawn. Found by searching for the smallest square side whose integer offsets land closest to the real
# fractions (see the search that produced these): DIP8 op-amps (both LM1458 and LM741 are genuinely the same
# physical 8-pin DIP size) match almost exactly at side 8; the NPN transistor matches almost exactly at side 6.
DIP8_SIDE = 8
DIP8_OFFSETS = {"1": 1, "2": 3, "3": 5, "4": 7, "8": 1, "7": 3, "6": 5, "5": 7}   # top-to-bottom, both sides
NPN_SIDE = 6
NPN_BASE_OFFSET = 3       # left side, matches the real symbol's near-exact-middle base lead
NPN_CE_OFFSET = 4         # top/bottom side, matches the real symbol's collector/emitter x-position


def symbol_side_len(dev_type, generic_side_len):
    """The real package's required square side, when this device gets real symbol artwork - chosen so the
    Universe's own integer attachment offsets land where the real symbol's pin leads actually are (see above),
    not just "big enough to fit 4 pins a side". Anything else keeps the generic formula."""
    if dev_type in OPAMP_DUAL_PINOUT or dev_type in OPAMP_DIP8_PINOUT:
        return DIP8_SIDE
    if dev_type in TRANSISTOR_PINOUT:
        return NPN_SIDE
    return generic_side_len


def forced_pinout(dev_type, L):
    """{pin: (side, offset)} - the real package's fixed physical layout, only for the three device types that
    need one. L is the device's actual (already symbol_side_len-sized) square side."""
    if dev_type in OPAMP_DUAL_PINOUT:
        po = OPAMP_DUAL_PINOUT[dev_type]
    elif dev_type in OPAMP_DIP8_PINOUT:
        po = OPAMP_DIP8_PINOUT[dev_type]
    else:
        po = None
    if po is not None:
        out = {}
        for pin in po["left_order"]:
            out[pin] = ("left", DIP8_OFFSETS[pin])
        for pin in po["right_order"]:
            out[pin] = ("right", DIP8_OFFSETS[pin])
        return out
    if dev_type in TRANSISTOR_PINOUT:
        tp = TRANSISTOR_PINOUT[dev_type]
        return {tp["base"]: ("left", NPN_BASE_OFFSET), tp["collector"]: ("top", NPN_CE_OFFSET),
                tp["emitter"]: ("bottom", NPN_CE_OFFSET)}
    return {}


def _fit_xform(x0, y0, x1, y1, origin, size):
    """Non-uniform (independent x/y) scale+translate fitting `size` (normalized) into the real box [x0,x1]x
    [y0,y1] - unlike the library's own xforms, which use one scale derived from width alone and require the box
    to already be in the exact normalized aspect ratio. This is what lets a real (non-square-proportioned)
    symbol sit correctly inside this Universe's equal-sided square."""
    ox, oy = origin
    sx = (x1 - x0) / size[0]
    sy = (y1 - y0) / size[1]
    def T(nx, ny):
        return (x0 + (nx - ox) * sx, y0 + (ny - oy) * sy)
    return T


def _piecewise_map(breakpoints):
    """1-D piecewise-linear map through `breakpoints` ([(norm, real), ...], >=2, sorted by norm). Interior values
    interpolate between their bracketing pair; values outside the range extrapolate along the nearest segment's
    slope, so body/notch/triangle geometry beyond the outermost pin still transforms smoothly."""
    pts = sorted(breakpoints)
    def f(v):
        if v <= pts[0][0]:
            (n0, r0), (n1, r1) = pts[0], pts[1]
        elif v >= pts[-1][0]:
            (n0, r0), (n1, r1) = pts[-2], pts[-1]
        else:
            for i in range(len(pts) - 1):
                if pts[i][0] <= v <= pts[i + 1][0]:
                    (n0, r0), (n1, r1) = pts[i], pts[i + 1]
                    break
        return r0 + (v - n0) / (n1 - n0) * (r1 - r0)
    return f


def _fit_xform_y(x0, x1, origin_x, size_x, y_breakpoints):
    """Same job as _fit_xform, but the y-axis is a piecewise-linear map calibrated to the REAL per-pin
    attachment y-positions (y_breakpoints: [(norm_y, real_y), ...], one pair per pin actually wired) instead of
    one uniform scale - the real datasheet's pin spacing isn't evenly divisible by any reasonably-sized integer
    square side, so a uniform scale leaves a small but real gap at some pins even after resizing the square (see
    DUAL_OPAMP/LM1458: pin-to-pin spacing is 22/27/25, not equal). Because every point of the WHOLE symbol -
    body, notch, triangles, every lead waypoint - passes through this one T, and because each lead segment is
    already axis-aligned in normalized space (constant nx or constant ny), mapping x and y independently can
    never introduce a diagonal: a vertical segment (nx fixed) stays at one x; a horizontal segment (ny fixed)
    stays at one y. This is the "change the scale, not the points" fix - no per-lead shifting afterward."""
    sx = (x1 - x0) / size_x
    fy = _piecewise_map(y_breakpoints)
    def T(nx, ny):
        return (x0 + (nx - origin_x) * sx, fy(ny))
    return T


def _y_breakpoints(terminals, real_points):
    """One (norm_y, real_y) pair per distinct norm_y - left/right sides share the same norm_y per row (see
    DUAL_OPAMP_TERMINALS/OPAMP_DIP8_TERMINALS), and the router keeps both sides' real Y aligned per row too, so
    duplicates are averaged (should already agree) rather than left in, which would zero out a segment's span."""
    by_norm = {}
    for pin, (rx, ry) in real_points.items():
        if pin not in terminals:
            continue
        ny = terminals[pin][1]
        by_norm.setdefault(ny, []).append(ry)
    return [(ny, sum(rys) / len(rys)) for ny, rys in by_norm.items()]


def _poly(svg, T, norm_pts, **attrs):
    pts = " ".join(f"{px:.1f},{py:.1f}" for px, py in (T(nx, ny) for nx, ny in norm_pts))
    a = " ".join(f'{k.replace("_", "-")}="{v}"' for k, v in attrs.items())
    svg.append(f'<polyline points="{pts}" fill="none" stroke="black" stroke-width="1" class="symbol" {a}/>')


def _body_path(svg, T, body, notch):
    bl, bt, br, bb = body
    p0 = T(bl, bt)
    p_ns = T(*notch["start"])
    p_c1a, p_c1b, p_mid = T(*notch["c1a"]), T(*notch["c1b"]), T(*notch["mid"])
    p_c2a, p_c2b, p_ne = T(*notch["c2a"]), T(*notch["c2b"]), T(*notch["end"])
    p_tr, p_br, p_bl = T(br, bt), T(br, bb), T(bl, bb)
    svg.append(
        f'<path d="M {p0[0]:.1f},{p0[1]:.1f} L {p_ns[0]:.1f},{p_ns[1]:.1f} '
        f'C {p_c1a[0]:.1f},{p_c1a[1]:.1f} {p_c1b[0]:.1f},{p_c1b[1]:.1f} {p_mid[0]:.1f},{p_mid[1]:.1f} '
        f'C {p_c2a[0]:.1f},{p_c2a[1]:.1f} {p_c2b[0]:.1f},{p_c2b[1]:.1f} {p_ne[0]:.1f},{p_ne[1]:.1f} '
        f'L {p_tr[0]:.1f},{p_tr[1]:.1f} L {p_br[0]:.1f},{p_br[1]:.1f} '
        f'L {p_bl[0]:.1f},{p_bl[1]:.1f} Z" fill="none" stroke="black" stroke-width="1" class="symbol"/>'
    )


def _poly_anchored(svg, T, norm_pts, real_first, **attrs):
    """same as _poly, but shifted as a whole by (real_first - the library's own transformed terminal), so the
    lead starts exactly at the real wire-attachment point with NO gap - and, critically, every segment stays
    perfectly axis-aligned (a uniform translation can't introduce a slant; replacing only the outer point and
    leaving the rest fixed can and did). The library's own terminal fraction is close but not exact (real
    datasheet pin spacing isn't evenly divisible by any reasonably-sized integer square side), so a whole-path
    shift is the only way to make the two agree exactly without bending a lead off-axis."""
    t0 = T(*norm_pts[0])
    dx, dy = real_first[0] - t0[0], real_first[1] - t0[1]
    pts = [(px + dx, py + dy) for px, py in (T(nx, ny) for nx, ny in norm_pts)]
    spts = " ".join(f"{px:.1f},{py:.1f}" for px, py in pts)
    a = " ".join(f'{k.replace("_", "-")}="{v}"' for k, v in attrs.items())
    svg.append(f'<polyline points="{spts}" fill="none" stroke="black" stroke-width="1" class="symbol" {a}/>')


def dual_opamp_svg(x0, y0, x1, y1, ref, pinout, real_points=None):
    real_points = real_points or {}
    y_bps = _y_breakpoints(DUAL_OPAMP_TERMINALS, real_points)
    if len(y_bps) >= 2:
        T = _fit_xform_y(x0, x1, DUAL_OPAMP_NORM_ORIGIN[0], DUAL_OPAMP_NORM_SIZE[0], y_bps)
    else:
        T = _fit_xform(x0, y0, x1, y1, DUAL_OPAMP_NORM_ORIGIN, DUAL_OPAMP_NORM_SIZE)
    svg = []
    _body_path(svg, T, DUAL_OPAMP_BODY, DUAL_OPAMP_NOTCH_BEZIER)
    p0 = T(DUAL_OPAMP_BODY[0], DUAL_OPAMP_BODY[1])
    svg.append(f'<text x="{p0[0]+4:.1f}" y="{p0[1]-8:.1f}" font-size="8" font-family="monospace" text-anchor="start">{ref}</text>')
    for pin_name in pinout["left_order"]:
        _poly(svg, T, DUAL_OPAMP_PIN_PATHS[pin_name])
        tx, ty = real_points.get(pin_name, T(*DUAL_OPAMP_TERMINALS[pin_name]))
        svg.append(f'<text x="{tx-3:.1f}" y="{ty+3:.1f}" font-size="7" font-family="monospace" text-anchor="end" class="symbol">{pin_name}</text>')
    for pin_name in pinout["right_order"]:
        _poly(svg, T, DUAL_OPAMP_PIN_PATHS[pin_name])
        tx, ty = real_points.get(pin_name, T(*DUAL_OPAMP_TERMINALS[pin_name]))
        svg.append(f'<text x="{tx+3:.1f}" y="{ty+3:.1f}" font-size="7" font-family="monospace" text-anchor="start" class="symbol">{pin_name}</text>')
    for section, spec in (("A", DUAL_OPAMP_TRIANGLE_A), ("B", DUAL_OPAMP_TRIANGLE_B)):
        _poly(svg, T, [spec["base_l"], spec["apex"], spec["base_r"], spec["base_l"]])
        lx, ly = T(*spec["label"])
        svg.append(f'<text x="{lx:.1f}" y="{ly:.1f}" font-size="9" font-family="monospace" text-anchor="middle" class="symbol">{section}</text>')
    return svg


def opamp8_svg(x0, y0, x1, y1, ref, pinout, real_points=None):
    real_points = real_points or {}
    y_bps = _y_breakpoints(OPAMP_DIP8_TERMINALS, real_points)
    if len(y_bps) >= 2:
        T = _fit_xform_y(x0, x1, OPAMP_DIP8_NORM_ORIGIN[0], OPAMP_DIP8_NORM_SIZE[0], y_bps)
    else:
        T = _fit_xform(x0, y0, x1, y1, OPAMP_DIP8_NORM_ORIGIN, OPAMP_DIP8_NORM_SIZE)
    svg = []
    _body_path(svg, T, OPAMP_DIP8_BODY, OPAMP_DIP8_NOTCH_BEZIER)
    p0 = T(OPAMP_DIP8_BODY[0], OPAMP_DIP8_BODY[1])
    svg.append(f'<text x="{p0[0]+4:.1f}" y="{p0[1]-8:.1f}" font-size="8" font-family="monospace" text-anchor="start">{ref}</text>')
    for pin_name in pinout["left_order"]:
        _poly(svg, T, OPAMP_DIP8_LEAD_PATHS[pin_name])
        tx, ty = real_points.get(pin_name, T(*OPAMP_DIP8_TERMINALS[pin_name]))
        svg.append(f'<text x="{tx-3:.1f}" y="{ty+3:.1f}" font-size="7" font-family="monospace" text-anchor="end" class="symbol">{pin_name}</text>')
    for pin_name in pinout["right_order"]:
        _poly(svg, T, OPAMP_DIP8_LEAD_PATHS[pin_name])
        tx, ty = real_points.get(pin_name, T(*OPAMP_DIP8_TERMINALS[pin_name]))
        svg.append(f'<text x="{tx+3:.1f}" y="{ty+3:.1f}" font-size="7" font-family="monospace" text-anchor="start" class="symbol">{pin_name}</text>')
    tri = OPAMP_DIP8_TRIANGLE
    _poly(svg, T, [tri["base_top"], tri["base_bottom"], tri["apex"], tri["base_top"]])
    return svg


def npn_svg(x0, y0, x1, y1, ref, real_points=None):
    T = _fit_xform(x0, y0, x1, y1, NPN_NORM_ORIGIN, NPN_NORM_SIZE)
    real_points = real_points or {}
    svg = []
    p0 = T(NPN_NORM_ORIGIN[0], NPN_NORM_ORIGIN[1])
    svg.append(f'<text x="{p0[0]+4:.1f}" y="{p0[1]-4:.1f}" font-size="8" font-family="monospace" text-anchor="start">{ref}</text>')
    cx, cy = T(NPN_BODY["cx"], NPN_BODY["cy"])
    rx = NPN_BODY["r"] * (T(1, 0)[0] - T(0, 0)[0])
    ry = NPN_BODY["r"] * (T(0, 1)[1] - T(0, 0)[1])
    svg.append(f'<ellipse cx="{cx:.1f}" cy="{cy:.1f}" rx="{rx:.1f}" ry="{ry:.1f}" fill="none" stroke="black" stroke-width="1" class="symbol"/>')
    _poly(svg, T, NPN_BASE_BAR)
    for role in ("base", "collector", "emitter"):
        rp = real_points.get(role, T(*NPN_TERMINALS[role]))
        _poly_anchored(svg, T, NPN_LEAD_PATHS[role], rp)
    tri = " ".join(f"{px:.1f},{py:.1f}" for px, py in (T(nx, ny) for nx, ny in NPN_ARROWHEAD))
    svg.append(f'<polygon points="{tri}" fill="black" stroke="none" class="symbol"/>')
    return svg


def two_pin_svg(dev_type, p1, p2, ref, label_pt):
    svg = []
    if dev_type == "GENERIC_RESISTOR":
        segs = [resistor_symbol_points(p1, p2)]
    elif dev_type == "GENERIC_CAPACITOR":
        segs = capacitor_symbol_segments(p1, p2)
    elif dev_type == "ZENER_2PIN":
        segs = zener_symbol_segments(p1, p2)
    else:  # 1N4148 - p1/p2 must be (anode, cathode)
        segs = diode_symbol_segments(p1, p2)
    for seg in segs:
        pts = " ".join(f"{px:.1f},{py:.1f}" for px, py in seg)
        svg.append(f'<polyline points="{pts}" fill="none" stroke="black" stroke-width="1" class="symbol"/>')
    lx, ly = label_pt
    svg.append(f'<text x="{lx:.1f}" y="{ly:.1f}" font-size="8" font-family="monospace" text-anchor="middle">{ref}</text>')
    return svg
