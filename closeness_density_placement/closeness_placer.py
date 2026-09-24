#!/usr/bin/env python3
"""Constructive placer: closeness (device-to-device, from the .circuit file's own nets) + density (one tunable
target area per device -> one target distance between connected devices), on the Point Universe's integer grid.

Definitions used, not derived from any of the earlier (abandoned) routing/placement work:
  * closeness: two devices are close if a net connects them. A net with k terminals contributes weight 1/(k-1)
    between every pair of devices on it (so a 2-terminal net ties its two devices together fully; a 22-terminal
    rail like 0V says almost nothing about any single pair). Every net counts, power included, straight from the
    .circuit file - nothing is excluded by convention here.
  * density: DENSITY_TARGET_AREA (quanta^2, ONE tunable number) -> TARGET_DISTANCE = sqrt(DENSITY_TARGET_AREA),
    the same for every connected pair (density is uniform, so is the target distance).
  * device footprint: Square, straight from the Universe's own definition (universe_defs.py, byte-identical to
    the untouched original, tunnel definition restored) - equal sides, side length the smallest that gives the
    device's pins 1 slot each (4*(L-1) >= pin count).
  * placement ("orbit" rule, replacing an earlier least-squares fit that could get stuck collapsing devices onto
    a straight line - see the R3/R4/C2/R5-around-U1 case): the most-connected device sits at the centre. Every
    other device is placed at exactly the target distance from its single STRONGEST already-placed closeness tie
    (its "anchor") - not fit against every tie at once - at a random angle around that anchor, rejecting angles
    too close to siblings already orbiting the same anchor. Ties for strongest anchor go to whichever candidate
    is itself more globally connected. This is what makes several devices that share one op-amp end up at close
    to the same distance from it, spread around it, instead of stacked in a line at growing distances.
"""
import json, math, random, shutil, subprocess, sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from universe_defs import Square

BASE = Path("/home/rgunn/Chat-Projects/graph-relax/output/hrng_baseline")

DENSITY_TARGET_AREA = 256          # quanta^2 per device - THE tunable knob; placeholder default (see report)
TARGET_DISTANCE = round(math.sqrt(DENSITY_TARGET_AREA))
MIN_GAP = 2                        # minimum clear quanta between two device footprints (density's hard floor)
MIN_ORBIT_ANGLE = math.radians(20) # siblings orbiting the same anchor stay at least this far apart, angularly
ORBIT_TRIES = 200                  # bounded random draws per device before settling for the best one found
SEED = 20260921                    # fixed seed: the randomness is reproducible, not different every run


def split_term(t):
    ref, pin = t.rsplit(".", 1)
    return ref, pin


def side_len(n_pins):
    L = 2
    while 4 * (L - 1) < n_pins:
        L += 1
    return L


def build_closeness(sem):
    comps = {c["ref"]: c for c in sem["electrical"]["components"]}
    weight = {}
    for net in sem["electrical"]["nets"]:
        refs = sorted({split_term(t)[0] for t in net["terminals"]})
        k = len(net["terminals"])
        if k < 2:
            continue
        w = 1.0 / (k - 1)
        for i in range(len(refs)):
            for j in range(i + 1, len(refs)):
                key = (refs[i], refs[j])
                weight[key] = weight.get(key, 0.0) + w
    nbrs = {r: {} for r in comps}
    for (a, b), w in weight.items():
        nbrs[a][b] = w; nbrs[b][a] = w
    return comps, nbrs


def place(comps, nbrs):
    rng = random.Random(SEED)
    L = {r: side_len(len(comps[r]["pins"])) for r in comps}
    N = len(comps)
    total_w = {r: sum(nbrs[r].values()) for r in comps}
    order0 = sorted(comps, key=lambda r: (-total_w[r], r))
    start = order0[0]
    pos = {start: (0, 0)}
    placed = [start]
    remaining = set(comps) - {start}
    orbits = {r: [] for r in comps}          # anchor -> [angle, angle, ...] of its already-placed orbiters

    def attach_weight(r):
        return sum(nbrs[r].get(p, 0.0) for p in placed)

    def overlaps(r, x, y, skip=None):
        l = L[r]
        for p in placed:
            if p == skip:
                continue
            lp = L[p]; px, py = pos[p]
            if x - MIN_GAP < px + lp and px - MIN_GAP < x + l and y - MIN_GAP < py + lp and py - MIN_GAP < y + l:
                return True
        return False

    def min_gap_to_siblings(anchor, angle):
        existing = orbits[anchor]
        if not existing:
            return math.pi
        return min(abs(((angle - a + math.pi) % (2 * math.pi)) - math.pi) for a in existing)

    nudges = 0
    while remaining:
        r = max(remaining, key=lambda r: (attach_weight(r), r))
        anchor = max((p for p in placed if p in nbrs[r]),
                     key=lambda p: (nbrs[r][p], total_w[p], p))
        ax, ay = pos[anchor]
        best = None
        for _ in range(ORBIT_TRIES):
            a = rng.uniform(0, 2 * math.pi)
            gap = min_gap_to_siblings(anchor, a)
            if best is None or gap > best[0]:
                best = (gap, a)
            if gap >= MIN_ORBIT_ANGLE:
                break
        angle = best[1]
        x = round(ax + TARGET_DISTANCE * math.cos(angle))
        y = round(ay + TARGET_DISTANCE * math.sin(angle))
        step = 0
        while overlaps(r, x, y) and step < 60:
            step += 1; nudges += 1
            x = round(ax + (TARGET_DISTANCE + step) * math.cos(angle))
            y = round(ay + (TARGET_DISTANCE + step) * math.sin(angle))
        pos[r] = (x, y)
        orbits[anchor].append(angle)
        placed.append(r); remaining.discard(r)

    minx = min(pos[r][0] - L[r] // 2 for r in comps); miny = min(pos[r][1] - L[r] // 2 for r in comps)
    squares = {}
    for r in comps:
        x = pos[r][0] - L[r] // 2 - minx + 4
        y = pos[r][1] - L[r] // 2 - miny + 4
        squares[r] = Square(r, x, y, L[r], L[r])
    W = max(s.x + s.width for s in squares.values()) + 5
    H = max(s.y + s.height for s in squares.values()) + 5
    return squares, W, H, nudges


def report(comps, nbrs, squares):
    strong = [(a, b, w) for a in nbrs for b, w in nbrs[a].items() if a < b and w >= 0.999]
    def center(r):
        s = squares[r]; return (s.x + s.width / 2, s.y + s.height / 2)
    ds = []
    for a, b, w in strong:
        ax, ay = center(a); bx, by = center(b)
        ds.append(math.hypot(ax - bx, ay - by))
    if ds:
        print(f"2-terminal-net pairs (weight 1.0, the strongest closeness ties): {len(ds)}; "
              f"actual distance vs target {TARGET_DISTANCE}: mean {sum(ds)/len(ds):.1f}, "
              f"min {min(ds):.1f}, max {max(ds):.1f}")
    for hub in ("U1",):
        kids = [r for r in comps if r != hub and nbrs[r].get(hub, 0) >= 0.4]
        dists = [math.hypot(*(a - b for a, b in zip(center(r), center(hub)))) for r in kids]
        if dists:
            print(f"devices strongly tied to {hub} ({len(kids)}): distance from it - "
                  f"mean {sum(dists)/len(dists):.1f}, min {min(dists):.1f}, max {max(dists):.1f}")


EDGE_DRAW_THRESHOLD = 0.4   # only ties from nets with <=3 terminals show up in the picture; weaker (rail) ties
                            # still count for placement, they just aren't drawn - drawing all of them (300+ edges
                            # for 35 devices, mostly from 0V/+9V/etc.) makes every picture unreadable regardless
                            # of placement quality


def render(squares, nbrs, W, H, out_png, S=14):
    MARGIN_Q = 6   # equal blank margin (quanta) on all four sides of the schematic itself, not the header
    HEADER = 44    # pixel strip above the margined schematic, for the title only
    xs = [s.x for s in squares.values()] + [s.x + s.width for s in squares.values()]
    ys = [s.y for s in squares.values()] + [s.y + s.height for s in squares.values()]
    cx0, cx1 = min(xs) - MARGIN_Q, max(xs) + MARGIN_Q
    cy0, cy1 = min(ys) - MARGIN_Q, max(ys) + MARGIN_Q
    Wp, Hp = (cx1 - cx0) * S, (cy1 - cy0) * S + HEADER
    svg = [f'<svg xmlns="http://www.w3.org/2000/svg" width="{Wp}" height="{Hp}" viewBox="0 0 {Wp} {Hp}">',
           f'<rect width="{Wp}" height="{Hp}" fill="white"/>',
           f'<text x="10" y="26" font-size="18" font-family="DejaVu Sans, sans-serif">HRNG - closeness + density placement, '
           f'orbit rule (preview: straight closeness-graph edges, not routed wires; density target {DENSITY_TARGET_AREA} q^2/device, '
           f'target distance {TARGET_DISTANCE}q)</text>']
    X = lambda x: (x - cx0) * S; Y = lambda y: (y - cy0) * S + HEADER
    seen = set()
    drawn = skipped = 0
    for a in nbrs:
        for b, w in nbrs[a].items():
            if a >= b or (a, b) in seen:
                continue
            seen.add((a, b))
            if w < EDGE_DRAW_THRESHOLD:
                skipped += 1
                continue
            drawn += 1
            sa, sb = squares[a], squares[b]
            ax, ay = sa.x + sa.width / 2, sa.y + sa.height / 2
            bx, by = sb.x + sb.width / 2, sb.y + sb.height / 2
            op = min(1.0, 0.15 + 0.7 * w)
            svg.append(f'<line x1="{X(ax)}" y1="{Y(ay)}" x2="{X(bx)}" y2="{Y(by)}" stroke="#1f77b4" '
                       f'stroke-width="{1 + 2*w:.2f}" stroke-opacity="{op:.2f}"/>')
    for ref, s in sorted(squares.items()):
        svg.append(f'<rect x="{X(s.x)}" y="{Y(s.y)}" width="{s.width*S}" height="{s.height*S}" fill="#f4f4f4" stroke="black" stroke-width="1.3"/>')
        svg.append(f'<text x="{X(s.x + s.width/2)}" y="{Y(s.y + s.height/2) + 4}" font-size="11" font-family="DejaVu Sans, sans-serif" text-anchor="middle">{ref}</text>')
    svg.append("</svg>")
    print(f"drawing {drawn} closeness ties (weight >= {EDGE_DRAW_THRESHOLD}), leaving out {skipped} weaker ones "
          f"(placement still used all of them, just not drawn)")
    out_png = Path(out_png); sv = out_png.with_suffix(".svg"); sv.write_text("\n".join(svg))
    # no -trim: canvas is already sized to content + an exact equal margin on every side, by construction
    subprocess.run(["convert", "-background", "white", str(sv), "-depth", "8", str(out_png)], check=True)
    for f in (out_png, sv):
        shutil.copy2(f, Path("/mnt/ZYXEL/DLINK/dump") / f.name); shutil.copy2(f, Path("/home/rgunn/picaxe") / f.name)
    return out_png


def main():
    sem = json.loads((BASE / "semantic.json").read_text())
    comps, nbrs = build_closeness(sem)
    squares, W, H, nudges = place(comps, nbrs)
    print(f"density target {DENSITY_TARGET_AREA} q^2/device -> target distance {TARGET_DISTANCE}q; "
          f"{len(comps)} devices; window {W} x {H} (aspect {W/H:.2f}); collision nudges: {nudges}")
    report(comps, nbrs, squares)
    out = render(squares, nbrs, W, H, HERE / "hrng_closeness_density.png")
    print("saved", out, "(copied to ~/picaxe and dump)")


if __name__ == "__main__":
    main()
