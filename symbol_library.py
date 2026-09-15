import math

# Every device's schematic artwork -- geometry constants, coordinate
# transforms, and the draw functions that turn them into SVG -- lives in
# this file, separate from schematic_router.py's placement/routing/
# verification code. render_svg queries SYMBOL_LIBRARY (built at the
# bottom of this file) to decide whether a device type has real artwork;
# if not, it draws a plain labeled box and reports the device type as
# needing a real symbol here. Adding a new device to an EXISTING family
# (another part number in TRANSISTOR_PINOUT, OPAMP_DUAL_PINOUT,
# OPAMP_DIP8_PINOUT, SW_DPST_PINOUT, or another simple 2-pin type
# following the SYMBOL_2PIN_TYPES pattern) is a data-only edit here --
# schematic_router.py's dispatch, build_wires' terminal allocation, and
# verify()/electrical_equivalence_check() never need to change for that.
# A genuinely new pin-ARRANGEMENT shape still needs a new xform + draw
# function here, plus one registration touch in schematic_router.py's
# render_svg/build_wires/alloc()/device_box_pts() to wire the new family
# in -- this file only removes the need to touch those for data that
# fits a pattern already wired in.

SYMBOL_2PIN_TYPES = {'GENERIC_RESISTOR', 'GENERIC_CAPACITOR', '1N4148', 'ZENER_2PIN', 'GENERIC_LED', 'GENERIC_METER'}

# Real pin ORDER (not just count) for every symbol-drawn 2-pin type, used to
# look up each named pin's exact routed coordinate directly from
# pin_points instead of an order-agnostic touch-point search. Resistors and
# capacitors are symmetric (order doesn't matter, either pin can be "p1"),
# but a diode is polarized -- the anode MUST map to the triangle's base and
# the cathode to its bar side, so this table's order is load-bearing for
# those two types, not just documentation.
SYMBOL_PIN_ORDER = {
    'GENERIC_RESISTOR': ('1', '2'),
    'GENERIC_CAPACITOR': ('1', '2'),
    '1N4148': ('A', 'K'),
    'ZENER_2PIN': ('A', 'K'),
    'GENERIC_LED': ('A', 'K'),
    'GENERIC_METER': ('POS', 'NEG'),
}

# Single-section op-amp (LM741), real physical DIP-8 package -- same
# architecture as OPAMP_DUAL_PINOUT/dual_opamp_dip8_svg below: only the
# physical left/right pin ORDER is needed here (for forced_side/
# forced_rank in build_wires), not a role mapping -- which physical pin
# connects to which triangle vertex is entirely fixed in
# OPAMP_DIP8_LEAD_PATHS, never derived from a role dict at render time.
# Verified against the actual baseline netlists (mono741/HRNG): U1.2=IN-,
# U1.3=IN+, U1.4=V-, U1.6=OUT, U1.7=V+ on every LM741 instance checked --
# matches the reference file's own drawn connections exactly (upper input
# pin row 2, lower input pin row 3, output pin row 3 on the right = pin
# 6). Pins 1/5 (offset null) and 8 (NC) carry no net in either circuit;
# they still get a real package-pin stub (no internal connector), same as
# the dual op-amp's V+/V- pins.
OPAMP_DIP8_PINOUT = {
    'LM741': {
        'left_order': ['1', '2', '3', '4'],
        'right_order': ['8', '7', '6', '5'],
    },
}

# Calculated vector single-section op-amp DIP-8 symbol -- same fixed-
# template + shared-transform architecture as the dual-op-amp and the NPN
# transistor: nothing here is ever reshaped by live routing data. Only the
# terminal points (the box's own left/right pin-row edges) are electrical;
# everything else is class="symbol" artwork. Coordinates copied verbatim
# from the authoritative reference file
# /mnt/ZYXEL/DLINK/dump/generic_8pin_opamp_traced.svg (117x168 canvas,
# body 20..103 x 17..150, small centered notch, 4 pins per side at
# y=34/67/100/134, right-pointing triangle apex at (90,84)).
OPAMP_DIP8_NORM_ORIGIN = (7, 17)   # leftmost pin x, topmost body y
OPAMP_DIP8_NORM_SIZE = (109, 133)  # 116-7, 150-17 -- bounds body AND pins, not just one

OPAMP_DIP8_TERMINALS = {  # electrical terminal, normalized -- ONLY these 8 points are electrical
    '1': (7, 34), '2': (7, 67), '3': (7, 100), '4': (7, 134),
    '8': (116, 34), '7': (116, 67), '6': (116, 100), '5': (116, 134),
}
OPAMP_DIP8_BODY = (20, 17, 103, 150)  # left, top, right, bottom, normalized
OPAMP_DIP8_NOTCH_BEZIER = {  # small centered notch, two cubic beziers, normalized
    'start': (49, 17), 'c1a': (49, 24), 'c1b': (53, 28), 'mid': (60, 28),
    'c2a': (67, 28), 'c2b': (71, 24), 'end': (71, 17),
}
OPAMP_DIP8_TRIANGLE = {'base_top': (47, 64), 'base_bottom': (47, 102), 'apex': (90, 84)}

# Terminal all the way to the body (triangle interior base point for the
# two inputs, apex for the output, plain package-wall stub for every
# other pin), one flat waypoint list per physical pin, copied verbatim
# from the reference file's own <path> elements. Graphical only.
OPAMP_DIP8_LEAD_PATHS = {
    '1': [(7, 34), (20, 34)],
    '2': [(7, 67), (20, 67), (34, 67), (34, 73), (47, 73)],
    '3': [(7, 100), (20, 100), (35, 100), (35, 95), (47, 95)],
    '4': [(7, 134), (20, 134)],
    '8': [(116, 34), (103, 34)],
    '7': [(116, 67), (103, 67)],
    '6': [(116, 100), (103, 100), (90, 100), (90, 84)],
    '5': [(116, 134), (103, 134)],
}


def opamp8_xform(bx0, by0, bx1, by1):
    # ONE transform, shared by every caller (drawing AND electrical
    # terminal placement) -- scale derived from width alone because
    # device_box_pts() always sizes this device's box in the EXACT
    # OPAMP_DIP8_NORM_SIZE proportion (109:133).
    ox, oy = OPAMP_DIP8_NORM_ORIGIN
    scale = (bx1 - bx0) / OPAMP_DIP8_NORM_SIZE[0]
    def T(nx, ny):
        return (bx0 + (nx-ox)*scale, by0 + (ny-oy)*scale)
    return T


OPAMP_DIP8_TRIANGLE_PINS = {'2', '3', '6'}  # IN-, IN+, OUT -- the only pins
# this glyph structurally requires to carry a real net. Pins 1/5 (offset
# null) and 8 (NC) are standard LM741 pins that no real circuit wires --
# they will NEVER appear in pin_points, by design, so gating on ALL 8
# pins here (a real bug this pass) made this function always return None
# and silently fall back to the generic box for every LM741. 4/7 (V-/V+)
# are expected to be wired in practice but aren't required structurally
# either -- the fixed lead artwork draws from OPAMP_DIP8_LEAD_PATHS
# regardless, with or without a live net at that pin.


def opamp8_symbol_svg(x0, y0, x1, y1, ref, pinout, pin_points):
    # Renders one single-section op-amp DIP-8 package as CALCULATED
    # VECTOR GEOMETRY -- same pattern as dual_opamp_dip8_svg/npn_symbol_
    # svg: one normalized shape, transformed into this device's real box.
    # No raster image, nothing reshaped by live routing data.
    p = {n: pin_points.get((ref, n)) for n in OPAMP_DIP8_TRIANGLE_PINS}
    if any(v is None for v in p.values()):
        return None  # caller falls back to the generic box

    T = opamp8_xform(x0, y0, x1, y1)
    svg = []
    svg.append(f'<rect x="{x0:.1f}" y="{y0:.1f}" width="{x1-x0:.1f}" height="{y1-y0:.1f}" fill="none" stroke="none"/>')
    svg.append(f'<text x="{x0+4:.1f}" y="{y0-4:.1f}" font-size="8" font-family="monospace" text-anchor="start">{ref}</text>')

    def poly(norm_pts):
        pts = ' '.join(f'{px:.1f},{py:.1f}' for px, py in (T(nx, ny) for nx, ny in norm_pts))
        svg.append(f'<polyline points="{pts}" fill="none" stroke="black" stroke-width="1" class="symbol"/>')

    # Package body outline with the small centered notch -- one closed
    # path, two cubic beziers for the bump (an affine scale+translate
    # carries a cubic bezier's control points through unchanged).
    bl, bt, br, bb = OPAMP_DIP8_BODY
    nb = OPAMP_DIP8_NOTCH_BEZIER
    p0 = T(bl, bt)
    p_ns = T(*nb['start'])
    p_c1a, p_c1b, p_mid = T(*nb['c1a']), T(*nb['c1b']), T(*nb['mid'])
    p_c2a, p_c2b, p_ne = T(*nb['c2a']), T(*nb['c2b']), T(*nb['end'])
    p_tr, p_br, p_bl = T(br, bt), T(br, bb), T(bl, bb)
    svg.append(
        f'<path d="M {p0[0]:.1f},{p0[1]:.1f} L {p_ns[0]:.1f},{p_ns[1]:.1f} '
        f'C {p_c1a[0]:.1f},{p_c1a[1]:.1f} {p_c1b[0]:.1f},{p_c1b[1]:.1f} {p_mid[0]:.1f},{p_mid[1]:.1f} '
        f'C {p_c2a[0]:.1f},{p_c2a[1]:.1f} {p_c2b[0]:.1f},{p_c2b[1]:.1f} {p_ne[0]:.1f},{p_ne[1]:.1f} '
        f'L {p_tr[0]:.1f},{p_tr[1]:.1f} L {p_br[0]:.1f},{p_br[1]:.1f} '
        f'L {p_bl[0]:.1f},{p_bl[1]:.1f} Z" fill="none" stroke="black" stroke-width="1" class="symbol"/>'
    )

    # One combined path per pin (terminal -> ... -> triangle, or terminal
    # -> wall for unused/power pins) + pin-number labels.
    for pin_name in pinout['left_order']:
        poly(OPAMP_DIP8_LEAD_PATHS[pin_name])
        tx, ty = T(*OPAMP_DIP8_TERMINALS[pin_name])
        svg.append(f'<text x="{tx-3:.1f}" y="{ty+3:.1f}" font-size="7" font-family="monospace" text-anchor="end" class="symbol">{pin_name}</text>')
    for pin_name in pinout['right_order']:
        poly(OPAMP_DIP8_LEAD_PATHS[pin_name])
        tx, ty = T(*OPAMP_DIP8_TERMINALS[pin_name])
        svg.append(f'<text x="{tx+3:.1f}" y="{ty+3:.1f}" font-size="7" font-family="monospace" text-anchor="start" class="symbol">{pin_name}</text>')

    # Right-pointing triangle.
    tri = OPAMP_DIP8_TRIANGLE
    poly([tri['base_top'], tri['base_bottom'], tri['apex'], tri['base_top']])

    return svg

# NPN transistor physical pin roles, per device type. P2N2222A verified
# two ways, not assumed from the symbol name alone: (1) the real ON
# Semiconductor/Motorola P2N2222A datasheet -- this specific TO-92 part
# marking is C-B-E pin1/2/3, NOT the E-B-C of the plastic PN2222A (a real,
# documented gotcha, easy to get backwards); (2) HRNG's own netlist --
# Q1.1->RNG_RAW (a collector output), Q1.2->Q1_BASE (explicitly named),
# Q1.3->0V (grounded emitter) -- independently agrees with C-B-E.
TRANSISTOR_PINOUT = {
    'P2N2222A': {'base': '2', 'collector': '1', 'emitter': '3'},
}

# 3-terminal potentiometer: end_a/end_b (the resistive element's two
# ends, forced onto opposite sides like a plain resistor -- see the
# forced_side block in build_wires) and wiper (the middle tap, side left
# free, same "leads land where routed" principle as every other symbol).
POTENTIOMETER_PINOUT = {
    'GENERIC_POTENTIOMETER': {'end_a': '1', 'wiper': '2', 'end_b': '3'},
}

# Calculated vector NPN transistor symbol -- same architecture as the
# dual-op-amp DIP-8 symbol: ONE fixed normalized template (never reshaped
# by where a real wire routes to -- see the long correction history on
# this in project memory if touching this again), transformed into the
# real box by a single uniform scale+translate, used identically by
# build_wires' alloc() (the real electrical terminal) and by the artwork
# renderer, so the drawn lead and the wire it really connects to can
# never drift apart. Coordinates below are copied verbatim from the
# authoritative reference file
# /mnt/ZYXEL/DLINK/dump/transistor_npn_traced_v2.svg (480x488 canvas,
# body circle center 242,254 radius 140).
NPN_NORM_ORIGIN = (31, 44)   # normalized point that maps to the box's (x0,y0)
# normalized width,height -- MUST bound the circle too, not just the three
# terminal points: the circle (center 242,254, r=140) bulges out to
# x=382 at its equator, well past the collector/emitter terminals' own
# x=266 -- a real, confirmed defect this pass sized this box from the
# terminals alone (235 wide), which let the circle stick out past the
# box's own right edge into space a neighboring wire was routed through.
# User: "the edge cuts through the circle... Routing error." Width is
# therefore (circle_right - origin_x) = (242+140) - 31 = 351, not
# (266-31)=235; height was already correctly bounded by the collector/
# emitter terminals (44/465), which exceed the circle's own y-extent
# (114/394), so it needed no change.
NPN_NORM_SIZE = (351, 421)

NPN_TERMINALS = {  # electrical terminal, normalized -- ONLY these 3 points are electrical
    'base': (31, 254), 'collector': (266, 44), 'emitter': (266, 465),
}
NPN_BODY = {'cx': 242, 'cy': 254, 'r': 140}
NPN_BASE_BAR = [(171, 172), (171, 336)]  # vertical base plate, normalized

# Terminal all the way to the body, one flat waypoint list per lead,
# copied verbatim from the reference's own <path> elements. Graphical
# only: never electrical, never routed, never added to nets or junction
# geometry.
NPN_LEAD_PATHS = {
    'base': [(31, 254), (171, 254)],
    'collector': [(266, 44), (266, 137), (171, 208)],
    'emitter': [(266, 465), (266, 372), (171, 302)],
}
# Filled arrowhead, fixed, aligned with the emitter lead and pointing
# outward (away from the base) -- exact vertices from the reference, not
# recomputed from a generic perpendicular-offset formula.
NPN_ARROWHEAD = [(239.40, 352.40), (209.01, 349.88), (227.99, 324.12)]


def npn_xform(bx0, by0, bx1, by1):
    # ONE transform, shared by every caller (drawing AND electrical
    # terminal placement): scale is derived from width alone because
    # device_box_pts() always sizes this device's box in the EXACT
    # NPN_NORM_SIZE proportion (235:421) -- so box_w/235 == box_h/421 by
    # construction.
    ox, oy = NPN_NORM_ORIGIN
    scale = (bx1 - bx0) / NPN_NORM_SIZE[0]
    def T(nx, ny):
        return (bx0 + (nx-ox)*scale, by0 + (ny-oy)*scale)
    return T


def npn_symbol_svg(x0, y0, x1, y1, ref, pinout, pin_points):
    # Renders one NPN transistor as CALCULATED VECTOR GEOMETRY -- same
    # pattern as dual_opamp_dip8_svg: one normalized shape, transformed
    # into this device's real box by npn_xform. No raster image, and
    # nothing here is reshaped by live routing data -- the base bar,
    # both leads, the circle, and the arrowhead are 100% fixed
    # proportions; only the box's own size (from device_box_pts) scales
    # the whole thing uniformly.
    #
    # ELECTRICAL geometry: none drawn here. The router already drew the
    # real external wire for base/collector/emitter elsewhere in
    # render_svg, terminating exactly at pin_points[(ref, pin)] -- which
    # build_wires placed using this SAME npn_xform, at NPN_TERMINALS.
    # That is the only electrically meaningful geometry this device has.
    #
    # Everything drawn below is tagged class="symbol" and is pure
    # artwork -- electrical_equivalence_check already skips any
    # class="symbol" polyline/polygon when reconstructing connectivity.
    p = {role: pin_points.get((ref, pinout[role])) for role in ('base', 'collector', 'emitter')}
    if any(v is None for v in p.values()):
        return None  # caller falls back to the generic box

    T = npn_xform(x0, y0, x1, y1)
    svg = []
    svg.append(f'<rect x="{x0:.1f}" y="{y0:.1f}" width="{x1-x0:.1f}" height="{y1-y0:.1f}" fill="none" stroke="none"/>')
    svg.append(f'<text x="{x0+4:.1f}" y="{y0-4:.1f}" font-size="8" font-family="monospace" text-anchor="start">{ref}</text>')

    def poly(norm_pts):
        pts = ' '.join(f'{px:.1f},{py:.1f}' for px, py in (T(nx, ny) for nx, ny in norm_pts))
        svg.append(f'<polyline points="{pts}" fill="none" stroke="black" stroke-width="1" class="symbol"/>')

    cx, cy = T(NPN_BODY['cx'], NPN_BODY['cy'])
    r = NPN_BODY['r'] * (T(1, 0)[0] - T(0, 0)[0])
    svg.append(f'<circle cx="{cx:.1f}" cy="{cy:.1f}" r="{r:.1f}" fill="none" stroke="black" stroke-width="1" class="symbol"/>')

    poly(NPN_BASE_BAR)
    for role in ('base', 'collector', 'emitter'):
        poly(NPN_LEAD_PATHS[role])

    tri = ' '.join(f'{px:.1f},{py:.1f}' for px, py in (T(nx, ny) for nx, ny in NPN_ARROWHEAD))
    svg.append(f'<polygon points="{tri}" fill="black" stroke="none" class="symbol"/>')

    return svg


# Calculated vector DPST switch symbol -- same fixed-template + shared-
# transform architecture as the transistor/dual-op-amp: nothing here is
# ever reshaped by live routing data. Two independent break-contact poles
# (pins 3/4 = top pole, pins 1/2 = bottom pole) plus the mechanical-gang
# linkage between them, standard multi-pole-switch convention. Coordinates
# copied verbatim from the reference file exported via `kicad-cli sym
# export svg` (Switch.kicad_sym "SW_DPST", 13.1x14.7mm canvas).
SW_DPST_PINOUT = {
    'GENERIC_SWITCH': {'left_order': ['3', '1'], 'right_order': ['4', '2']},
}
SW_DPST_NORM_ORIGIN = (1.4732, 3.2639)
SW_DPST_NORM_SIZE = (10.16, 6.604)  # 11.6332-1.4732, 9.8679-3.2639 -- bounds every element, not just terminals

SW_DPST_TERMINALS = {
    '3': (1.4732, 4.7879), '1': (1.4732, 9.8679),
    '4': (11.6332, 4.7879), '2': (11.6332, 9.8679),
}
SW_DPST_LEAD_PATHS = {  # terminal -> contact-dot-adjacent point; dot drawn separately
    '3': [(1.4732, 4.7879), (4.0132, 4.7879)],
    '1': [(1.4732, 9.8679), (4.0132, 9.8679)],
    '4': [(11.6332, 4.7879), (9.0932, 4.7879)],
    '2': [(11.6332, 9.8679), (9.0932, 9.8679)],
}
SW_DPST_DOTS = {  # contact points, radius 0.508, normalized
    '3': (4.5212, 4.7879), '1': (4.5212, 9.8679),
    '4': (8.5852, 4.7879), '2': (8.5852, 9.8679),
}
SW_DPST_DOT_R = 0.508
SW_DPST_BLADES = [  # open-contact blades, one per pole, fixed
    [(5.0292, 4.5339), (7.8232, 3.2639)],   # top pole (3-4)
    [(5.0292, 9.6139), (7.8232, 8.3439)],   # bottom pole (1-2)
]
SW_DPST_GANG_DASHES = [  # mechanical-linkage dashes between the two poles
    [(6.5532, 4.1529), (6.5532, 4.7879)],
    [(6.5532, 5.4229), (6.5532, 6.0579)],
    [(6.5532, 6.6929), (6.5532, 7.3279)],
    [(6.5532, 7.9629), (6.5532, 8.5979)],
]


def sw_dpst_xform(bx0, by0, bx1, by1):
    ox, oy = SW_DPST_NORM_ORIGIN
    scale = (bx1 - bx0) / SW_DPST_NORM_SIZE[0]
    def T(nx, ny):
        return (bx0 + (nx-ox)*scale, by0 + (ny-oy)*scale)
    return T


def sw_dpst_symbol_svg(x0, y0, x1, y1, ref, pinout, pin_points):
    p = {n: pin_points.get((ref, n)) for n in (pinout['left_order'] + pinout['right_order'])}
    if any(v is None for v in p.values()):
        return None  # caller falls back to the generic box

    T = sw_dpst_xform(x0, y0, x1, y1)
    svg = []
    svg.append(f'<rect x="{x0:.1f}" y="{y0:.1f}" width="{x1-x0:.1f}" height="{y1-y0:.1f}" fill="none" stroke="none"/>')
    svg.append(f'<text x="{x0+4:.1f}" y="{y0-4:.1f}" font-size="8" font-family="monospace" text-anchor="start">{ref}</text>')

    def poly(norm_pts):
        pts = ' '.join(f'{px:.1f},{py:.1f}' for px, py in (T(nx, ny) for nx, ny in norm_pts))
        svg.append(f'<polyline points="{pts}" fill="none" stroke="black" stroke-width="1" class="symbol"/>')

    for pin_name in pinout['left_order']:
        poly(SW_DPST_LEAD_PATHS[pin_name])
        tx, ty = T(*SW_DPST_TERMINALS[pin_name])
        svg.append(f'<text x="{tx-3:.1f}" y="{ty+3:.1f}" font-size="7" font-family="monospace" text-anchor="end" class="symbol">{pin_name}</text>')
    for pin_name in pinout['right_order']:
        poly(SW_DPST_LEAD_PATHS[pin_name])
        tx, ty = T(*SW_DPST_TERMINALS[pin_name])
        svg.append(f'<text x="{tx+3:.1f}" y="{ty+3:.1f}" font-size="7" font-family="monospace" text-anchor="start" class="symbol">{pin_name}</text>')

    scale = (x1 - x0) / SW_DPST_NORM_SIZE[0]
    for pin_name, (dcx, dcy) in SW_DPST_DOTS.items():
        cx, cy = T(dcx, dcy)
        svg.append(f'<circle cx="{cx:.1f}" cy="{cy:.1f}" r="{SW_DPST_DOT_R*scale:.1f}" fill="none" stroke="black" stroke-width="1" class="symbol"/>')
    for blade in SW_DPST_BLADES:
        poly(blade)
    for dash in SW_DPST_GANG_DASHES:
        pts = ' '.join(f'{px:.1f},{py:.1f}' for px, py in (T(nx, ny) for nx, ny in dash))
        svg.append(f'<polyline points="{pts}" fill="none" stroke="black" stroke-width="1" stroke-dasharray="2,1" class="symbol"/>')

    return svg


# Dual-section packages (two op-amp triangles inside ONE physical DIP body,
# not two independently-placed units). Pin numbers verified against TI's
# actual LM1458/LM1558 datasheet (Figure 2, Dual-In-Line Package) AND
# independently against HRNG's own pre-existing netlist (U1A_OUT/MINUS/PLUS,
# U1B_..., matching exactly, pin for pin) -- both agree: left side top to
# bottom is 1,2,3,4 = OUT_A,IN-A,IN+A,V-; right side top to bottom is
# 8,7,6,5 = V+,OUT_B,IN-B,IN+B. Every 8-pin dual op-amp (LM358, TL072,
# NE5532, OPA2604, ...) shares this same JEDEC-standard pinout.
OPAMP_DUAL_PINOUT = {
    'LM1458': {
        'left_order': ['1', '2', '3', '4'],
        'right_order': ['8', '7', '6', '5'],
        'A': {'out': '1', 'in_minus': '2', 'in_plus': '3'},
        'B': {'out': '7', 'in_minus': '6', 'in_plus': '5'},
        'vplus': '8', 'vminus': '4',
    },
}



def _lerp(p1, p2, t):
    return (p1[0] + (p2[0]-p1[0])*t, p1[1] + (p2[1]-p1[1])*t)


def resistor_symbol_points(p1, p2):
    # Real IEEE zigzag, drawn directly between the two pin points the
    # router already computed -- since the pin-side fix above (see
    # SYMBOL_2PIN_TYPES / fixed_side in build_wires) guarantees these two
    # points are diametrically opposite on the box, the zigzag reads as a
    # straight axial component with zero extra routing: the existing wire
    # polylines already terminate exactly at p1/p2, so this glyph is the
    # only thing drawn INSIDE the box, seamlessly continuing each lead.
    dx, dy = p2[0]-p1[0], p2[1]-p1[1]
    length = math.hypot(dx, dy)
    if length < 1e-6:
        return [p1, p2]
    ux, uy = dx/length, dy/length
    vx, vy = -uy, ux
    body_lo, body_hi = 0.28, 0.72
    n = 6
    amp = 6.0
    pts = [p1, _lerp(p1, p2, body_lo)]
    for i in range(1, n):
        t = body_lo + (body_hi - body_lo) * i / n
        cx, cy = _lerp(p1, p2, t)
        s = amp if i % 2 == 1 else -amp
        pts.append((cx + vx*s, cy + vy*s))
    pts.append(_lerp(p1, p2, body_hi))
    pts.append(p2)
    return pts


def potentiometer_symbol_segments(p_end_a, p_end_b, p_wiper):
    # Same IEEE zigzag as resistor_symbol_points, between the two end
    # terminals -- reused directly. The wiper arrow points from the
    # zigzag's own midpoint OUTWARD toward wherever the real wiper pin
    # actually landed (not a fixed perpendicular assumption), same
    # "leads land exactly where routed" principle as every other symbol
    # here. Real KiCad R_Potentiometer_US proportions: filled arrowhead,
    # tip at the zigzag's middle peak, base out toward the wiper lead.
    segs = [resistor_symbol_points(p_end_a, p_end_b)]
    mx, my = _lerp(p_end_a, p_end_b, 0.5)
    dx, dy = p_wiper[0]-mx, p_wiper[1]-my
    length = math.hypot(dx, dy)
    if length < 1e-6:
        return segs
    ux, uy = dx/length, dy/length
    vx, vy = -uy, ux
    tip = (mx, my)
    base_c = (mx + ux*7.5, my + uy*7.5)
    half_w = 4.0
    b1 = (base_c[0]+vx*half_w, base_c[1]+vy*half_w)
    b2 = (base_c[0]-vx*half_w, base_c[1]-vy*half_w)
    segs.append([b1, tip, b2, b1])
    segs.append([base_c, p_wiper])
    return segs


def capacitor_symbol_segments(p1, p2):
    # Two parallel plates centered between the pin points, same reasoning
    # as resistor_symbol_points -- p1/p2 are real routed-wire endpoints,
    # not independently chosen.
    dx, dy = p2[0]-p1[0], p2[1]-p1[1]
    length = math.hypot(dx, dy)
    if length < 1e-6:
        return [[p1, p2]]
    ux, uy = dx/length, dy/length
    vx, vy = -uy, ux
    half_gap = min(4.0, length*0.12)
    plate_half = 7.0
    mid = _lerp(p1, p2, 0.5)
    c1 = (mid[0]-ux*half_gap, mid[1]-uy*half_gap)
    c2 = (mid[0]+ux*half_gap, mid[1]+uy*half_gap)
    plate1 = [(c1[0]+vx*plate_half, c1[1]+vy*plate_half), (c1[0]-vx*plate_half, c1[1]-vy*plate_half)]
    plate2 = [(c2[0]+vx*plate_half, c2[1]+vy*plate_half), (c2[0]-vx*plate_half, c2[1]-vy*plate_half)]
    return [[p1, c1], plate1, [c2, p2], plate2]


def diode_symbol_segments(p_anode, p_cathode):
    # Standard triangle+bar diode glyph, current-flow direction anode ->
    # cathode -- POLARIZED, unlike resistor/capacitor, so p_anode/p_cathode
    # must be the real named pins (looked up from pin_points via
    # SYMBOL_PIN_ORDER's ('A','K') order in render_svg), never an
    # order-agnostic touch-point pair. Same "leads land exactly on the
    # router's own already-computed pin coordinates" principle as the
    # other two symbols.
    dx, dy = p_cathode[0]-p_anode[0], p_cathode[1]-p_anode[1]
    length = math.hypot(dx, dy)
    if length < 1e-6:
        return [[p_anode, p_cathode]]
    ux, uy = dx/length, dy/length
    vx, vy = -uy, ux
    half_w = 6.0
    base_c = _lerp(p_anode, p_cathode, 0.40)
    tip = _lerp(p_anode, p_cathode, 0.60)
    bar_c = tip
    base1 = (base_c[0]+vx*half_w, base_c[1]+vy*half_w)
    base2 = (base_c[0]-vx*half_w, base_c[1]-vy*half_w)
    bar1 = (bar_c[0]+vx*half_w, bar_c[1]+vy*half_w)
    bar2 = (bar_c[0]-vx*half_w, bar_c[1]-vy*half_w)
    return [[p_anode, base_c], [base1, base2, tip, base1], [bar1, bar2], [bar_c, p_cathode]]


def zener_symbol_segments(p_anode, p_cathode):
    # Same triangle+bar as diode_symbol_segments (reused, not
    # reimplemented) plus the Zener's characteristic two-footed "Z" bend
    # on the cathode bar -- geometry measured directly off the user's own
    # reference (/mnt/ZYXEL/DLINK/dump/Zener.svg): each end of the bar
    # kicks out diagonally, rotationally symmetric about the bar's center
    # (one foot bends toward the anode side AND further outward on one
    # end, the other bends toward the cathode side AND further outward on
    # the other end) -- NOT a single one-sided foot (an earlier version of
    # this, matched to KiCad's own "D_Zener" export instead, only bent one
    # end; the user's reference is the authoritative shape). Measured
    # ratios off that reference (foot's along-axis kick vs. its outward
    # kick, both relative to the bar's own half-length): ~0.4 along-axis,
    # ~0.33 outward.
    segs = diode_symbol_segments(p_anode, p_cathode)
    dx, dy = p_cathode[0]-p_anode[0], p_cathode[1]-p_anode[1]
    length = math.hypot(dx, dy)
    if length < 1e-6:
        return segs
    ux, uy = dx/length, dy/length
    vx, vy = -uy, ux
    half_w = 6.0
    bar_c = _lerp(p_anode, p_cathode, 0.60)  # same point diode_symbol_segments calls bar_c/tip
    bar1 = (bar_c[0]+vx*half_w, bar_c[1]+vy*half_w)
    bar2 = (bar_c[0]-vx*half_w, bar_c[1]-vy*half_w)
    tab_u, tab_v = 0.4*half_w, 0.33*half_w
    foot1 = (bar1[0]+vx*tab_v+ux*tab_u, bar1[1]+vy*tab_v+uy*tab_u)
    foot2 = (bar2[0]-vx*tab_v-ux*tab_u, bar2[1]-vy*tab_v-uy*tab_u)
    return segs + [[bar1, foot1], [bar2, foot2]]


def led_symbol_segments(p_anode, p_cathode):
    # Same triangle+bar as diode_symbol_segments (an LED is electrically a
    # diode; current flows anode->cathode identically) -- reused directly,
    # not reimplemented, so any future diode-glyph tweak carries over
    # automatically. The circular enclosure and two emission arrows are
    # what distinguish an LED -- geometry measured directly off the
    # user's own reference (/mnt/ZYXEL/DLINK/dump/LED.svg), not KiCad's:
    # its own proportions land the enclosure circle at exactly r=12.0 --
    # the SAME absolute radius already used for the meter's body -- and
    # centered halfway between the triangle base and the cathode bar
    # (_lerp .50, between diode_symbol_segments' own .40/.60 split, which
    # this reference's triangle/bar positions independently confirm). The
    # two arrows sit at ~59 degrees off the anode->cathode axis, not 45 --
    # an earlier version of this, matched to KiCad's own "LED" export
    # instead, used 45; the user's reference is the authoritative shape.
    # The second arrow is the first shifted along its own perpendicular,
    # not along the diode's v-axis (measured ratio, not assumed).
    segs = diode_symbol_segments(p_anode, p_cathode)
    dx, dy = p_cathode[0]-p_anode[0], p_cathode[1]-p_anode[1]
    length = math.hypot(dx, dy)
    circle_c = _lerp(p_anode, p_cathode, 0.50)
    circle_r = 12.0
    if length < 1e-6:
        return segs, (circle_c, circle_r)
    ux, uy = dx/length, dy/length
    vx, vy = -uy, ux
    bar_c = _lerp(p_anode, p_cathode, 0.60)  # same point diode_symbol_segments calls bar_c/tip
    dir_x, dir_y = ux*0.514 + vx*-0.857, uy*0.514 + vy*-0.857
    perp_x, perp_y = -dir_y, dir_x
    tail0 = (bar_c[0] + ux*-6.28 + vx*-4.69, bar_c[1] + uy*-6.28 + vy*-4.69)
    for shift in (0.0, 2.42):
        tail = (tail0[0] + perp_x*shift, tail0[1] + perp_y*shift)
        tip = (tail[0] + dir_x*6.6, tail[1] + dir_y*6.6)
        back = (tip[0] - dir_x*1.7, tip[1] - dir_y*1.7)
        h1 = (back[0] + perp_x*0.8, back[1] + perp_y*0.8)
        h2 = (back[0] - perp_x*0.8, back[1] - perp_y*0.8)
        segs.append([tail, tip])
        segs.append([h1, tip, h2])
    return segs, (circle_c, circle_r)


def meter_symbol_segments(p_pos, p_neg):
    # Real KiCad Ammeter_DC proportions (Device.kicad_sym): a circle
    # (radius 2.54, local units) with the leads entering top/bottom, plus
    # a small "+" polarity tick next to the positive lead -- this
    # circuit's own M1.POS/M1.NEG pin names are exactly this polarity, so
    # the tick is placed at whichever real point is p_pos, not assumed.
    # The circle's own "A" label is drawn separately by the caller (a
    # <circle> isn't a polyline, so this symbol needs its own render_svg
    # hook, same as the transistor's circle).
    dx, dy = p_neg[0]-p_pos[0], p_neg[1]-p_pos[1]
    length = math.hypot(dx, dy)
    if length < 1e-6:
        return [[p_pos, p_neg]], (p_pos, p_pos, 1.0)
    ux, uy = dx/length, dy/length
    vx, vy = -uy, ux
    r = 12.0
    center = ((p_pos[0]+p_neg[0])/2, (p_pos[1]+p_neg[1])/2)
    edge_pos = (center[0]-ux*r, center[1]-uy*r)
    edge_neg = (center[0]+ux*r, center[1]+uy*r)
    segs = [[p_pos, edge_pos], [p_neg, edge_neg]]
    # "+" tick, offset to the side of the positive lead, clear of the
    # circle itself.
    tick_c = (p_pos[0] - ux*4.0 + vx*6.0, p_pos[1] - uy*4.0 + vy*6.0)
    segs.append([(tick_c[0]-ux*2.0, tick_c[1]-uy*2.0), (tick_c[0]+ux*2.0, tick_c[1]+uy*2.0)])
    segs.append([(tick_c[0]-vx*2.0, tick_c[1]-vy*2.0), (tick_c[0]+vx*2.0, tick_c[1]+vy*2.0)])
    return segs, (center, edge_pos, r)


# Calculated vector dual-op-amp DIP-8 symbol -- same architecture as
# resistor/capacitor/diode: fixed geometry in one normalized coordinate
# system, transformed into real coordinates by a single uniform scale +
# translate, never a raster image. Coordinates below are copied verbatim
# from the authoritative reference file
# /mnt/ZYXEL/DLINK/dump/dual_opamp_reference_traced.svg (body
# 30..102 x 18..123, notch a two-bezier bump, terminals at x=23/109).
#
# The transform is ONE formula, used for every element -- package body,
# notch, pin stubs, internal artwork, AND the real electrical terminals
# (build_wires' alloc() calls dual_opamp_xform for the latter) -- so
# artwork and electrical terminals can never drift apart.
DUAL_OPAMP_NORM_ORIGIN = (23, 18)   # normalized point that maps to the box's (x0,y0)
DUAL_OPAMP_NORM_SIZE = (86, 105)    # normalized width,height (109-23, 123-18)

DUAL_OPAMP_TERMINALS = {  # electrical terminal, normalized -- ONLY these 8 points are electrical
    '1': (23, 34), '2': (23, 56), '3': (23, 83), '4': (23, 108),
    '8': (109, 34), '7': (109, 56), '6': (109, 83), '5': (109, 108),
}
DUAL_OPAMP_BODY = (30, 18, 102, 123)  # left, top, right, bottom, normalized
DUAL_OPAMP_NOTCH_BEZIER = {  # rounded top-notch bump, two cubic beziers, normalized
    'start': (56, 18), 'c1a': (56, 24), 'c1b': (60, 27), 'mid': (66, 27),
    'c2a': (72, 27), 'c2b': (76, 24), 'end': (76, 18),
}

DUAL_OPAMP_TRIANGLE_A = {'apex': (56, 43), 'base_l': (42, 67), 'base_r': (70, 67), 'label': (56, 61)}
DUAL_OPAMP_TRIANGLE_B = {'apex': (76, 70), 'base_l': (61, 95), 'base_r': (91, 95), 'label': (76, 89)}

# ONE combined path per physical pin, outer electrical terminal all the
# way to its triangle (inputs land on an interior base point, never a
# vertex, per the reference's own contract note), or terminal to the
# package wall only for the two supply pins -- copied verbatim from the
# reference file's own per-pin stub + internal-connector elements.
# Graphical only: never electrical, never routed, never added to nets or
# junction geometry, regardless of how closely a bend approaches a
# triangle vertex.
DUAL_OPAMP_PIN_PATHS = {
    '1': [(23, 34), (30, 34), (56, 34), (56, 43)],
    '2': [(23, 56), (30, 56), (36, 56), (36, 74), (52, 74), (52, 67)],
    '3': [(23, 83), (30, 83), (60, 83), (60, 67)],
    '4': [(23, 108), (30, 108)],
    '8': [(109, 34), (102, 34)],
    '7': [(109, 56), (102, 56), (76, 56), (76, 70)],
    '6': [(109, 83), (102, 83), (96, 83), (96, 103), (81, 103), (81, 95)],
    '5': [(109, 108), (102, 108), (68, 108), (68, 95)],
}


def dual_opamp_xform(bx0, by0, bx1, by1):
    # ONE transform, shared by every caller (drawing AND electrical
    # terminal placement): scale is derived from width alone because
    # device_box_pts() always sizes this device's box in the EXACT
    # DUAL_OPAMP_NORM_SIZE proportion (86:105) -- so box_w/86 ==
    # box_h/105 by construction.
    ox, oy = DUAL_OPAMP_NORM_ORIGIN
    scale = (bx1 - bx0) / DUAL_OPAMP_NORM_SIZE[0]
    def T(nx, ny):
        return (bx0 + (nx-ox)*scale, by0 + (ny-oy)*scale)
    return T


def dual_opamp_dip8_svg(x0, y0, x1, y1, ref, pinout, pin_points):
    # Renders one dual-op-amp DIP-8 package as CALCULATED VECTOR GEOMETRY
    # -- same pattern as resistor_symbol_points/diode_symbol_segments: one
    # normalized shape, transformed into this device's real box by
    # dual_opamp_xform. No raster image anywhere.
    #
    # ELECTRICAL geometry: none drawn here. The router already drew the
    # real external wire for every pin elsewhere in render_svg, terminating
    # exactly at pin_points[(ref, pin_number)] -- which build_wires placed
    # using this SAME dual_opamp_xform, at DUAL_OPAMP_TERMINALS. That is
    # the only electrically meaningful geometry this device has.
    #
    # Everything drawn below (body, notch, pin numbers, both triangles,
    # A/B labels, per-pin connector paths) is tagged class="symbol" and is
    # pure artwork -- electrical_equivalence_check already skips any
    # class="symbol" polyline when reconstructing connectivity, so none of
    # it can be routed, added to nets, or treated as junction/wire
    # geometry.
    p = {n: pin_points.get((ref, n)) for n in (pinout['left_order'] + pinout['right_order'])}
    if any(v is None for v in p.values()):
        return None  # caller falls back to the generic box

    T = dual_opamp_xform(x0, y0, x1, y1)
    svg = []

    def poly(norm_pts):
        pts = ' '.join(f'{px:.1f},{py:.1f}' for px, py in (T(nx, ny) for nx, ny in norm_pts))
        svg.append(f'<polyline points="{pts}" fill="none" stroke="black" stroke-width="1" class="symbol"/>')

    # Package body outline with the rounded top notch -- one closed path,
    # two cubic beziers for the bump (an affine scale+translate carries a
    # cubic bezier's control points through unchanged, so each point below
    # is simply transformed individually).
    bl, bt, br, bb = DUAL_OPAMP_BODY
    nb = DUAL_OPAMP_NOTCH_BEZIER
    p0 = T(bl, bt)
    p_ns = T(*nb['start'])
    p_c1a, p_c1b, p_mid = T(*nb['c1a']), T(*nb['c1b']), T(*nb['mid'])
    p_c2a, p_c2b, p_ne = T(*nb['c2a']), T(*nb['c2b']), T(*nb['end'])
    p_tr, p_br, p_bl = T(br, bt), T(br, bb), T(bl, bb)
    svg.append(
        f'<path d="M {p0[0]:.1f},{p0[1]:.1f} L {p_ns[0]:.1f},{p_ns[1]:.1f} '
        f'C {p_c1a[0]:.1f},{p_c1a[1]:.1f} {p_c1b[0]:.1f},{p_c1b[1]:.1f} {p_mid[0]:.1f},{p_mid[1]:.1f} '
        f'C {p_c2a[0]:.1f},{p_c2a[1]:.1f} {p_c2b[0]:.1f},{p_c2b[1]:.1f} {p_ne[0]:.1f},{p_ne[1]:.1f} '
        f'L {p_tr[0]:.1f},{p_tr[1]:.1f} L {p_br[0]:.1f},{p_br[1]:.1f} '
        f'L {p_bl[0]:.1f},{p_bl[1]:.1f} Z" fill="none" stroke="black" stroke-width="1" class="symbol"/>'
    )
    svg.append(f'<text x="{p0[0]+4:.1f}" y="{p0[1]-8:.1f}" font-size="8" font-family="monospace" text-anchor="start">{ref}</text>')

    # One combined path per pin (terminal -> ... -> triangle, or terminal
    # -> wall for the two supply pins) + pin-number labels.
    for pin_name in pinout['left_order']:
        poly(DUAL_OPAMP_PIN_PATHS[pin_name])
        tx, ty = T(*DUAL_OPAMP_TERMINALS[pin_name])
        svg.append(f'<text x="{tx-3:.1f}" y="{ty+3:.1f}" font-size="7" font-family="monospace" text-anchor="end" class="symbol">{pin_name}</text>')
    for pin_name in pinout['right_order']:
        poly(DUAL_OPAMP_PIN_PATHS[pin_name])
        tx, ty = T(*DUAL_OPAMP_TERMINALS[pin_name])
        svg.append(f'<text x="{tx+3:.1f}" y="{ty+3:.1f}" font-size="7" font-family="monospace" text-anchor="start" class="symbol">{pin_name}</text>')

    # Both triangles + their A/B labels.
    for section, spec in (('A', DUAL_OPAMP_TRIANGLE_A), ('B', DUAL_OPAMP_TRIANGLE_B)):
        poly([spec['base_l'], spec['apex'], spec['base_r'], spec['base_l']])
        lx, ly = T(*spec['label'])
        svg.append(f'<text x="{lx:.1f}" y="{ly:.1f}" font-size="9" font-family="monospace" text-anchor="middle" class="symbol">{section}</text>')

    return svg

# Every device type with real, library-provided artwork. render_svg
# queries this set once per device, BEFORE its existing dispatch chain
# runs, purely to detect and report a missing symbol -- it does not
# change what that chain draws.
SYMBOL_LIBRARY = (set(SYMBOL_2PIN_TYPES) | set(POTENTIOMETER_PINOUT) | set(TRANSISTOR_PINOUT)
                  | set(OPAMP_DIP8_PINOUT) | set(SW_DPST_PINOUT) | set(OPAMP_DUAL_PINOUT))
