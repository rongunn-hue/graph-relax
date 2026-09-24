#!/usr/bin/env python3
"""Route the signal wires through the finished closeness+density+group placement, inside a real Point Universe.

Fresh code (not built on any of the earlier abandoned routing work) - only the Universe's own rules are load-
bearing here: orthogonal lines, no shared points between different nets, tunnels as the one legal way to cross
another net's wire (kept in universe_defs.py because the user said it would be useful).

Steps:
  1. take the finished device layout from grouped_placer (positions + sizes only, nothing else reused)
  2. give every device real Attachments: wired (signal-net) pins point toward the net's other terminal(s), in
     the same angular order as the wire actually leaves the device, so nothing has to cross right at the body;
     power-net pins (never wired, drawn as stubs by convention) fill the leftover slots
  3. commit every Square into one PointUniverse
  4. route every signal net (2-terminal = one wire; 3-terminal = a wire plus a T-join), shortest legal orthogonal
     path first-and-only pass, no tunnels yet - report whatever doesn't route rather than forcing it
"""
import copy, json, math, shutil, subprocess, sys, heapq, itertools
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from universe_defs import Attachment, HypotheticalWhiteConsequence, Line, LineEndpoint, PointState, PointUniverse, Square, WhiteSpaceViewer
from router_defs import RouteCandidate
from route_selector_defs import RouteSelector
import closeness_placer as CP
import grouped_placer as GP
import topology as TOPO
import symbols as SYM

MARGIN = 6
DIR_ORDERS = (((0, -1), (-1, 0), (1, 0), (0, 1)), ((1, 0), (0, 1), (-1, 0), (0, -1)))
BENDS = (0, 1, 4, 10)
BIG = (1 << 60, 1 << 60)


def crossable(universe, q, d, net):
    """Can a tunnel of `net` heading d pass over point q? (a straight 90-degree stretch of another net's wire) -
    same rule as the old hrng_universe_route_tunnels.py, reused at the user's request ("use our old router code")."""
    ls = universe.point_lines.get(q)
    if not ls or len(ls) != 1:
        return False
    (lid,) = ls
    line = universe.lines[lid]
    if line.tunnel or line.net_id == net:
        return False
    pts = universe.line_points[lid]
    ax, ay = abs(d[1]), abs(d[0])
    if not ((q[0] + ax, q[1] + ay) in pts and (q[0] - ax, q[1] - ay) in pts):
        return False
    return (q[0] + abs(d[0]), q[1] + abs(d[1])) not in pts and (q[0] - abs(d[0]), q[1] - abs(d[1])) not in pts


def find_route(universe, white, start, goal_pts, goal_att, net, bend, dirs):
    heap = [((0, 0), 0, start, None)]
    best = {(start, None): (0, 0)}
    prev = {}
    tie = itertools.count(1)
    while heap:
        cost, _, p, hd = heapq.heappop(heap)
        if cost > best.get((p, hd), BIG):
            continue
        if p in goal_pts or p in goal_att:
            moves, cur = [], (p, hd)
            while cur in prev:
                pp, phd, kind, crossed = prev[cur]
                moves.append((pp, cur[0], kind, crossed))
                cur = (pp, phd)
            moves.reverse()
            return moves
        tn, ln = cost
        for d in dirs:
            if hd is not None and d == (-hd[0], -hd[1]):
                continue
            q = (p[0] + d[0], p[1] + d[1])
            if not universe.in_bounds(q):
                continue
            turn = bend if hd is not None and d != hd else 0
            if q in white or q in goal_pts or q in goal_att:
                nc, key = (tn, ln + 1 + turn), (q, d)
                if nc < best.get(key, BIG):
                    best[key] = nc; prev[key] = (p, hd, "step", None)
                    heapq.heappush(heap, (nc, next(tie), q, d))
            elif crossable(universe, q, d, net):
                crossed, r = [], q
                while crossable(universe, r, d, net):
                    crossed.append(r); r = (r[0] + d[0], r[1] + d[1])
                if universe.in_bounds(r) and (r in white or r in goal_att):
                    nc, key = (tn + 1, ln + len(crossed) + 1 + turn), (r, d)
                    if nc < best.get(key, BIG):
                        best[key] = nc; prev[key] = (p, hd, "jump", crossed)
                        heapq.heappush(heap, (nc, next(tie), r, d))
    return None


def path_points_and_pieces(start, moves):
    pts, pieces, cur = [start], [], [start]
    for p, q, kind, crossed in moves:
        if kind == "step":
            cur.append(q); pts.append(q)
        else:
            if len(cur) >= 2:
                pieces.append(("wire", cur))
            pieces.append(("tunnel", (p, q)))
            pts.extend(crossed + [q]); cur = [q]
    if len(cur) >= 2:
        pieces.append(("wire", cur))
    return tuple(pts), pieces


def candidate_routes(universe, start, goal_pts, goal_att, net):
    white = universe.layer(PointState.WHITE)
    viewer = WhiteSpaceViewer(universe)
    seen, out = {}, []
    for bend in BENDS:
        for dirs in DIR_ORDERS:
            moves = find_route(universe, white, start, goal_pts, goal_att, net, bend, dirs)
            if not moves:
                continue
            pts, pieces = path_points_and_pieces(start, moves)
            if pts in seen or len(pts) < 2:
                continue
            seen[pts] = pieces
            if any(kind == "jump" for _, _, kind, _ in moves):
                removed = [p for p in pts if p in white]
                after = set(viewer._region_at) - set(removed)
                sizes = viewer._white_component_sizes(after)
                rb = viewer.white_region_count()
                cons = HypotheticalWhiteConsequence(rb, len(sizes), len(sizes) - rb, len(removed), tuple(sorted(sizes)))
            elif len(pts) > 2:
                cons = viewer.hypothetical_white_consequence(pts[1:-1])
            else:
                sizes = viewer._white_component_sizes(set(viewer._region_at))
                rb = viewer.white_region_count()
                cons = HypotheticalWhiteConsequence(rb, rb, 0, 0, tuple(sorted(sizes)))
            out.append(RouteCandidate(pts, cons))
    return out, seen


def build_lines(lid, net, pieces, start_decl, end_decl):
    lines, prev_id = [], None
    for k, (kind, geom) in enumerate(pieces):
        pid = f"{lid}.{k}"
        first, last = k == 0, k == len(pieces) - 1
        s = start_decl if first else LineEndpoint(line_id=prev_id)
        e = end_decl if last else LineEndpoint()
        if kind == "wire":
            lines.append(Line(pid, net, corners(tuple(geom)), endpoints=(s, e)))
        else:
            lines.append(Line(pid, net, tuple(geom), endpoints=(s, e), tunnel=True))
        prev_id = pid
    return lines


def perimeter_slots(L):
    r = range(1, L)
    return ([("top", o) for o in r] + [("right", o) for o in r]
            + [("bottom", o) for o in reversed(r)] + [("left", o) for o in reversed(r)])


def local_point(L, side, off):
    return {"top": (off, 0), "right": (L, off), "bottom": (off, L), "left": (0, off)}[side]


def build_attachments(comps, nbrs_unused, nets, squares_xyl, forced_slots=None):
    """squares_xyl: ref -> (x, y, L) from the placement. `forced_slots`: {(ref,pin): (side,offset)} - pins pinned
    to a specific slot before the rest are assigned (used to swap a blocked pin onto a never-wired power/NC slot,
    since which physical slot a pin sits at is our own choice, not an electrical fact)."""
    forced_slots = forced_slots or {}
    centers = {r: (x + L / 2, y + L / 2) for r, (x, y, L) in squares_xyl.items()}
    pins_of_net = {}
    net_of_pin = {}
    for n in nets:
        terms = [CP.split_term(t) for t in n["terminals"]]
        pins_of_net[n["name"]] = terms
        for r, p in terms:
            net_of_pin[(r, p)] = n["name"]

    squares = {}
    for r in sorted(comps):
        x, y, L = squares_xyl[r]
        sq = Square(r, x, y, L, L)
        cx, cy = centers[r]
        pins = [str(p) for p in comps[r]["pins"]]
        signal = []
        power = []
        for p in pins:
            net = net_of_pin.get((r, p))
            if net is None:
                power.append(p); continue           # NC pin: no net at all, treat like a power pin (leftover slot)
            others = [(rr, pp) for rr, pp in pins_of_net[net] if rr != r]
            is_power = net in ("+5V", "+9V", "-9V", "0V")
            if others and not is_power:
                ox = sum(centers[rr][0] for rr, pp in others) / len(others)
                oy = sum(centers[rr][1] for rr, pp in others) / len(others)
                ang = math.atan2(oy - cy, ox - cx)
                signal.append((ang, p))
            else:
                power.append(p)
        signal.sort()
        slots = perimeter_slots(L)
        slot_ang = [(math.atan2(local_point(L, s, o)[1] - L / 2, local_point(L, s, o)[0] - L / 2), (s, o))
                    for s, o in slots]
        avail = sorted(slot_ang)
        used = set()
        assign = {}
        for p in pins:
            if (r, p) in forced_slots:
                so = forced_slots[(r, p)]
                assign[p] = so; used.add(so)
        # A 2-pin symbol glyph (resistor zigzag, capacitor plates, diode arrow...) is drawn as a straight shape
        # directly between its two pins, and only looks right if those two pins are on OPPOSITE sides of the
        # square - adjacent sides (e.g. left+top) draw a diagonal that reads as an X, not the real symbol. Force
        # opposite sides for every such device, choosing the axis from the signal target direction when there is
        # one, defaulting to left/right for a device with no signal net at all (a power-only floater).
        dev_type = comps[r]["device"]
        if dev_type in SYM.SYMBOL_2PIN_TYPES and len(pins) == 2 and not any(p in [pp for pp in assign] for p in pins):
            pin_a, pin_b = SYM.SYMBOL_PIN_ORDER[dev_type]
            ref_ang = signal[0][0] if signal else 0.0
            horiz = abs(math.cos(ref_ang)) >= abs(math.sin(ref_ang))
            if horiz:
                side_a, side_b = ("left", "right") if math.cos(ref_ang) <= 0 else ("right", "left")
            else:
                side_a, side_b = ("top", "bottom") if math.sin(ref_ang) <= 0 else ("bottom", "top")
            # whichever of the two pins actually carries the signal (if any) gets the side facing its target
            first_pin = signal[0][1] if signal else pin_a
            order = [first_pin] + [p for p in (pin_a, pin_b) if p != first_pin]
            assign[order[0]] = (side_a, 1); used.add((side_a, 1))
            assign[order[1]] = (side_b, 1); used.add((side_b, 1))
        signal = [(ang, p) for ang, p in signal if p not in assign]
        power = [p for p in power if p not in assign]
        for ang, p in signal:
            best_i = min(range(len(avail)), key=lambda i: abs(((avail[i][0] - ang + math.pi) % (2 * math.pi)) - math.pi)
                         if avail[i][1] not in used else 1e9)
            s2, o2 = avail[best_i][1]
            used.add((s2, o2))
            assign[p] = (s2, o2)
        rest = [so for a, so in avail if so not in used]
        for p, so in zip(sorted(power), rest):
            assign[p] = so
        for p in pins:
            s, o = assign[p]
            sq.add_attachment(f"{r}.{p}", s, o)
        squares[r] = sq
    return squares


def corners(path):
    out = [path[0]]
    for a, b, c in zip(path, path[1:], path[2:]):
        if (b[0] - a[0], b[1] - a[1]) != (c[0] - b[0], c[1] - b[1]):
            out.append(b)
    out.append(path[-1])
    return tuple(out)


def expand_points(corner_pts):
    """corner points -> every point along the line - a T-join is legal anywhere on a line's interior, not just
    at its bend points, so the join-target set has to be the full expansion, not the bare corner list."""
    out = [corner_pts[0]]
    for a, b in zip(corner_pts, corner_pts[1:]):
        dx = (b[0] > a[0]) - (b[0] < a[0]); dy = (b[1] > a[1]) - (b[1] < a[1])
        x, y = a
        while (x, y) != b:
            x += dx; y += dy; out.append((x, y))
    return out


def route(universe, signal_nets):
    """Same shape as the old hrng_universe_route_tunnels.py main() loop: netlist order, RouteSelector's own law
    (delta_R, length, points) among candidates, a trial on a deepcopy before ever committing for real, tunnels
    used only when candidate_routes finds no tunnel-free path at all."""
    import time
    t0 = time.time()
    selector = RouteSelector()
    committed, blocked, tunnel_uses = [], [], []
    # user: route the device that has a blockage first - U1/U2 now have all 8 pins forced onto just their left
    # and right sides (the real DIP8 package layout), which crowds those two edges far more than before; give
    # their nets first pick of the surrounding space instead of leaving them to fight over whatever's left.
    PRIORITY_DEVICES = {"U1", "U2"}
    def priority(n):
        refs = {CP.split_term(t)[0] for t in n["terminals"]}
        return 0 if refs & PRIORITY_DEVICES else 1
    signal_nets = sorted(signal_nets, key=priority)
    for n in signal_nets:
        name = n["name"]
        print(f"  routing {name} ... [{time.time()-t0:.1f}s]", flush=True)
        terms = [f"{r}.{p}" for r, p in (CP.split_term(t) for t in n["terminals"])]
        join_owner = {}
        wire1_ok = False
        for i in range(len(terms) - 1):
            a = terms[0] if i == 0 else terms[2]
            b = terms[1] if i == 0 else None
            lid = f"{name}#{i}"
            if i == 1 and not wire1_ok:
                blocked.append((lid, a, "T on wire 1 (wire 1 not routed)")); continue
            start = universe.attachment_points[a]
            if i == 0:
                goal_pts, goal_att = set(), {universe.attachment_points[b]}
            else:
                goal_pts, goal_att = set(join_owner), set()
            cands, pieces_by_path = candidate_routes(universe, start, goal_pts, goal_att, name)
            if not cands:
                blocked.append((lid, a, b or "T on wire 1")); continue
            dRs = sorted({c.consequence.delta_R for c in cands})
            lens = sorted({len(c.points) for c in cands})
            print(f"    {lid}: {len(cands)} candidates, delta_R values {dRs}, lengths {lens}", flush=True)
            chosen = selector.select(cands)
            pieces = pieces_by_path[chosen.points]
            start_decl = LineEndpoint(attachment_id=a)
            end_decl = LineEndpoint(attachment_id=b) if i == 0 else LineEndpoint(line_id=join_owner[chosen.points[-1]])
            lines = build_lines(lid, name, pieces, start_decl, end_decl)
            trial = copy.deepcopy(universe)
            try:
                for ln in lines:
                    trial.add_line(ln)
            except ValueError as e:
                blocked.append((lid, a, str(e))); continue
            for ln in lines:
                universe.add_line(ln)
            if i == 0:
                wire1_ok = True
                for ln in lines:
                    if not ln.tunnel:
                        for pt in expand_points(list(ln.points))[1:-1]:
                            join_owner[pt] = ln.line_id
            ntun = sum(1 for ln in lines if ln.tunnel)
            if ntun:
                tunnel_uses.append((lid, ntun))
            committed.extend(lines)
    return committed, blocked, tunnel_uses


def power_or_nc_pins(comps, nets, ref):
    """pins on `ref` that are never wired at all (power rails or NC) - safe to swap a blocked pin onto, since
    nothing already-committed depends on where these sit."""
    wired = set()
    for n in nets:
        if n["name"] in ("+5V", "+9V", "-9V", "0V"):
            continue
        for t in n["terminals"]:
            r, p = CP.split_term(t)
            if r == ref:
                wired.add(p)
    return [str(p) for p in comps[ref]["pins"] if str(p) not in wired]


def run_once(comps, nbrs, squares0, W, H, nets, signal_nets, forced_slots=None):
    xyl = {r: (s.x, s.y, s.width) for r, s in squares0.items()}
    squares = build_attachments(comps, nbrs, nets, xyl, forced_slots)
    universe = PointUniverse(W + MARGIN, H + MARGIN, line_separation=0, square_clearance=0)
    universe.defer_refresh = True
    for r in sorted(squares):
        universe.add_square(squares[r])
    universe.defer_refresh = False
    universe._refresh_empty_space()
    committed, blocked, tunnel_uses = route(universe, signal_nets)
    return universe, squares, committed, blocked, tunnel_uses


def main():
    sem = json.loads((CP.BASE / "semantic.json").read_text())
    comps, nbrs, pin_fn = TOPO.build_weighted_closeness(sem, CP.split_term)
    squares0, W, H, nudges, group_of = GP.place_grouped(comps, nbrs, sem)
    print(f"derived groups: {sorted(set(group_of.values()))}")
    nets = sem["electrical"]["nets"]
    signal_nets = [n for n in nets if n["name"] not in ("+5V", "+9V", "-9V", "0V")]

    forced = {}
    for r in comps:
        dev_type = comps[r]["device"]
        if dev_type in SYM.DEVICE_TYPES_NEEDING_FORCED_PINOUT:
            L = squares0[r].width
            for pin, slot in SYM.forced_pinout(dev_type, L).items():
                forced[(r, pin)] = slot
    print(f"forced real physical pinout for symbol-drawn devices: {sorted(set(r for r, p in forced))}")
    universe, squares, committed, blocked, tunnel_uses = run_once(comps, nbrs, squares0, W, H, nets, signal_nets, forced)
    print(f"first pass: {len(blocked)} blocked: {[(lid, a, b) for lid, a, b in blocked]}")

    # a blocked wire's start pin (`a`) may just be sitting in an unlucky slot - the slot a pin occupies is our own
    # choice, not an electrical fact, so try swapping it onto a same-device pin that's never wired at all (a power
    # or NC pin - safe, since nothing committed depends on where THAT one sits) and see if it actually routes.
    tried = set()
    for lid, a, b in list(blocked):
        ref, pin = a.split(".", 1)
        if (ref, pin) in forced or (ref, pin) in tried:
            continue
        tried.add((ref, pin))
        for swap_pin in power_or_nc_pins(comps, nets, ref):
            if (ref, swap_pin) in forced:
                continue
            trial_forced = dict(forced)
            xyl = {r: (s.x, s.y, s.width) for r, s in squares0.items()}
            base = build_attachments(comps, nbrs, nets, xyl, forced)
            cur_slot = next((a2.side, a2.offset) for a2 in base[ref].attachments if a2.attachment_id == f"{ref}.{pin}")
            swap_slot = next((a2.side, a2.offset) for a2 in base[ref].attachments if a2.attachment_id == f"{ref}.{swap_pin}")
            trial_forced[(ref, pin)] = swap_slot
            trial_forced[(ref, swap_pin)] = cur_slot
            print(f"  trying {ref}.{pin} <-> {ref}.{swap_pin} ...")
            u2, sq2, c2, b2, t2 = run_once(comps, nbrs, squares0, W, H, nets, signal_nets, trial_forced)
            if len(b2) < len(blocked):
                print(f"  IMPROVED: {len(blocked)} -> {len(b2)} blocked; keeping the swap")
                forced = trial_forced
                universe, squares, committed, blocked, tunnel_uses = u2, sq2, c2, b2, t2
                break
            else:
                print(f"  no improvement ({len(b2)} blocked), reverting")

    ntunnels = sum(n for _, n in tunnel_uses)
    wires_needed = sum(2 if len(n["terminals"]) == 3 else 1 for n in signal_nets)
    print(f"placement: window {W} x {H}; {len(squares)} devices, {len(universe.attachments)} attachments")
    print(f"routing: {len(committed)} line pieces committed, {len(blocked)} blocked (of {wires_needed} wires needed "
          f"for {len(signal_nets)} signal nets); tunnels used: {ntunnels} ({tunnel_uses}); pin swaps used: {forced}")
    for lid, a, b in blocked:
        print(f"  BLOCKED {lid}: {a} -> {b}")

    print()
    print("=== independent electrical check (geometry only, ignores routing labels) ===")
    import electrical_check as EC
    EC.check(universe, squares, sem, circuit_path="/home/rgunn/Chat-Projects/deterministic_schematic/circuit_examples/hrng.circuit")
    print()

    S = 8          # EXPERIMENT (2026-09-23): was 14 - lowered alongside grouped_placer.GAP to target ~1300x800px
    MARGIN_Q = 6   # equal blank margin (quanta) kept on all four sides of the actual schematic, not the header
    HEADER = 40    # pixel strip above the margin for the title, outside the centered/margined schematic area

    # true content bounds, computed directly (not by cropping an image after the fact) - devices, every wire
    # point, and the blocked-wire markers, so the margin is exact and equal on every side by construction
    xs = [sq.x for sq in squares.values()] + [sq.x + sq.width for sq in squares.values()]
    ys = [sq.y for sq in squares.values()] + [sq.y + sq.height for sq in squares.values()]
    for line in committed:
        xs += [p[0] for p in line.points]; ys += [p[1] for p in line.points]
    for lid, a, b in blocked:
        pa = universe.attachment_points.get(a)
        if pa:
            xs.append(pa[0]); ys.append(pa[1])
    cx0, cx1 = min(xs) - MARGIN_Q, max(xs) + MARGIN_Q
    cy0, cy1 = min(ys) - MARGIN_Q, max(ys) + MARGIN_Q

    pal = ["#d62728", "#1f77b4", "#2ca02c", "#9467bd", "#ff7f0e", "#8c564b", "#e377c2", "#17becf", "#bcbd22", "#7f7f7f"]
    netcolor = {n["name"]: pal[i % len(pal)] for i, n in enumerate(sorted(signal_nets, key=lambda n: n["name"]))}
    Wp, Hp = (cx1 - cx0) * S, (cy1 - cy0) * S + HEADER
    svg = [f'<svg xmlns="http://www.w3.org/2000/svg" width="{Wp}" height="{Hp}" viewBox="0 0 {Wp} {Hp}">',
           f'<rect width="{Wp}" height="{Hp}" fill="white"/>',
           f'<text x="10" y="26" font-size="18" font-family="DejaVu Sans, sans-serif">HRNG - grouped placement, wires routed '
           f'({len(committed)} pieces committed, {len(blocked)} blocked, {ntunnels} tunnels)</text>']
    X = lambda x: (x - cx0) * S; Y = lambda y: (y - cy0) * S + HEADER
    for ref, sq in sorted(squares.items()):
        dev_type = comps[ref]["device"]
        bx0, by0, bx1, by1 = X(sq.x), Y(sq.y), X(sq.x + sq.width), Y(sq.y + sq.height)
        if dev_type in SYM.OPAMP_DUAL_PINOUT:
            po = SYM.OPAMP_DUAL_PINOUT[dev_type]
            rp = {p: (X(universe.attachment_points[f"{ref}.{p}"][0]), Y(universe.attachment_points[f"{ref}.{p}"][1]))
                  for p in po["left_order"] + po["right_order"]}
            svg.extend(SYM.dual_opamp_svg(bx0, by0, bx1, by1, ref, po, rp))
        elif dev_type in SYM.OPAMP_DIP8_PINOUT:
            po = SYM.OPAMP_DIP8_PINOUT[dev_type]
            rp = {p: (X(universe.attachment_points[f"{ref}.{p}"][0]), Y(universe.attachment_points[f"{ref}.{p}"][1]))
                  for p in po["left_order"] + po["right_order"]}
            svg.extend(SYM.opamp8_svg(bx0, by0, bx1, by1, ref, po, rp))
        elif dev_type in SYM.TRANSISTOR_PINOUT:
            tp = SYM.TRANSISTOR_PINOUT[dev_type]
            rp = {role: (X(universe.attachment_points[f"{ref}.{pin}"][0]), Y(universe.attachment_points[f"{ref}.{pin}"][1]))
                  for role, pin in tp.items()}
            svg.extend(SYM.npn_svg(bx0, by0, bx1, by1, ref, rp))
        elif dev_type in SYM.SYMBOL_2PIN_TYPES:
            pin_a, pin_b = SYM.SYMBOL_PIN_ORDER[dev_type]
            p1 = universe.attachment_points.get(f"{ref}.{pin_a}")
            p2 = universe.attachment_points.get(f"{ref}.{pin_b}")
            if p1 and p2:
                p1x, p2x = (X(p1[0]), Y(p1[1])), (X(p2[0]), Y(p2[1]))
                svg.extend(SYM.two_pin_svg(dev_type, p1x, p2x, ref, (X(sq.x + sq.width / 2), Y(sq.y) - 3)))
            else:
                svg.append(f'<rect x="{bx0}" y="{by0}" width="{bx1-bx0}" height="{by1-by0}" fill="#f4f4f4" stroke="black" stroke-width="1.3"/>')
                svg.append(f'<text x="{X(sq.x+sq.width/2)}" y="{Y(sq.y)-3}" font-size="11" font-family="DejaVu Sans, sans-serif" text-anchor="middle">{ref}</text>')
        else:
            svg.append(f'<rect x="{bx0}" y="{by0}" width="{bx1-bx0}" height="{by1-by0}" fill="#f4f4f4" stroke="black" stroke-width="1.3"/>')
            svg.append(f'<text x="{X(sq.x+sq.width/2)}" y="{Y(sq.y)-3}" font-size="11" font-family="DejaVu Sans, sans-serif" text-anchor="middle">{ref}</text>')
    for line in committed:
        d = " ".join(f"{X(x)},{Y(y)}" for x, y in line.points)
        dash = ' stroke-dasharray="1 3"' if line.tunnel else ""
        w = 4 if line.tunnel else 2.2
        svg.append(f'<polyline points="{d}" fill="none" stroke="{netcolor[line.net_id]}" stroke-width="{w}"{dash}/>')
    # No junction dots: a dot marks where two DIFFERENT nets connect. That never happens in this Universe - a
    # tunnel is exactly what lets one net's wire cross another's without connecting, so the only way two nets'
    # drawn wires ever touch a point is a non-connecting tunnel crossing, which must NOT get a dot (dot = joined).
    # Same-net T-joins (two pieces of one net meeting) are the same net continuing - unambiguous by colour alone.
    for lid, a, b in blocked:
        pa = universe.attachment_points.get(a)
        if pa:
            svg.append(f'<circle cx="{X(pa[0])}" cy="{Y(pa[1])}" r="6" fill="none" stroke="red" stroke-width="2.5"/>')
    svg.append("</svg>")
    out_png = HERE / "hrng_routed.png"; sv = out_png.with_suffix(".svg"); sv.write_text("\n".join(svg))
    # no -trim here: the SVG canvas is already sized to content + an exact equal margin on every side, by
    # construction (see MARGIN_Q above) - trimming afterward would just undo that and re-introduce asymmetry
    subprocess.run(["convert", "-background", "white", str(sv), "-depth", "8", str(out_png)], check=True)
    for f in (out_png, sv):
        shutil.copy2(f, Path("/mnt/ZYXEL/DLINK/dump") / f.name); shutil.copy2(f, Path("/home/rgunn/picaxe") / f.name)
    print("saved", out_png)


if __name__ == "__main__":
    main()
