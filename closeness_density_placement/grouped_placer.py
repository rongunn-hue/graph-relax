#!/usr/bin/env python3
"""Two-level version of the orbit placer, tried per the user's idea: one density figure used two ways.

  * within a functional group: devices orbit each other at the density-implied spacing directly (tight - same
    number as the plain device-level placer, TARGET_DISTANCE).
  * between groups: NOT device distance at all. Each group's own footprint (the box its packed members actually
    occupy, from the within-group pass) becomes its size; groups then orbit each other at a distance built from
    THEIR sizes, the same way two devices would be size-spaced, plus the same TARGET_DISTANCE as the gap term.
    A group's size varies a lot by how many devices it holds, unlike individual parts (mostly the same size),
    which was the objection to sizing plain device-to-device distance by part size.

Groups are DERIVED, not hand-assigned (see topology.py): community detection on a closeness graph whose
feedback-loop nets (an op-amp's IN-/OUT pins, from the .circuit file's own pin-function metadata) are weighted
up. Replaces an earlier hardcoded 4-box mapping that mis-sorted R2 (an exact 50/50 tie between C1 and U1) into
"Amplification" just because a source label said so - see topology.py's docstring for the full reasoning.
Power-only devices (no signal net at all) are excluded from grouping and placed by the same single-strongest-tie
orbit rule as before, unchanged, as a separate pass after the derived groups are down.
"""
import json, math, random, sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from universe_defs import Square
import closeness_placer as CP
import topology as TOPO
GAP = 8                              # EXPERIMENT (2026-09-23): was CP.TARGET_DISTANCE (=16) - halved to compress
                                     # the layout toward a readable ~1300x800px target; original density-derived
                                     # value below is the tunable: clear space added to real size, everywhere -
                                     # not a flat replacement distance. target_distance(a,b) = (size_a+size_b)/2 + GAP,
                                     # used for device-device, device-anchor (power-only), and group-group alike, so
                                     # spacing is correct by construction (big IC vs small R get different spacing)
                                     # and the overlap-avoidance nudge stops being load-bearing.


def orbit_place(items, weight, size, target_distance, seed, start=None, preplaced=None):
    """generic version of closeness_placer.place()'s growth rule: most-connected item at the centre, everything
    else placed at target_distance(a, b) from its single strongest already-placed tie, at a random angle,
    spread away from siblings already orbiting the same anchor, and checked against EVERY already-placed item
    (not just its own anchor) using the same size-aware target_distance - so overlap avoidance comes from the
    same distance rule everything else uses, not a separate patch. Works for devices, groups, or (via `preplaced`)
    devices being added around an already-finished layout.
    `preplaced`: dict of ref -> fixed (x, y) for items that are already down and usable as anchors, but are not
    themselves being placed (e.g. the four functional groups' devices, when placing the power-only leftovers)."""
    rng = random.Random(seed)
    preplaced = preplaced or {}
    total_w = {r: sum(weight[r].values()) for r in items}
    pos = dict(preplaced)
    placed = list(preplaced)
    if start is None and not placed:
        start = sorted(items, key=lambda r: (-total_w[r], r))[0]
    if start is not None:
        pos[start] = (0.0, 0.0)
        placed.append(start)
    remaining = set(items) - set(placed)
    orbits = {r: [] for r in items}

    def attach_weight(r):
        return sum(weight[r].get(p, 0.0) for p in placed)

    def overlaps(r, x, y):
        for p in placed:
            d = target_distance(r, p) - CP.MIN_GAP
            if math.hypot(pos[p][0] - x, pos[p][1] - y) < max(d, 1.0):
                return True
        return False

    def min_gap_to_siblings(anchor, angle):
        existing = orbits.setdefault(anchor, [])
        return math.pi if not existing else min(
            abs(((angle - a + math.pi) % (2 * math.pi)) - math.pi) for a in existing)

    nudges = 0
    while remaining:
        r = max(remaining, key=lambda r: (attach_weight(r), r))
        anchor = max((p for p in placed if p in weight[r]), key=lambda p: (weight[r][p], total_w.get(p, 0.0), p))
        ax, ay = pos[anchor]
        T = target_distance(r, anchor)
        best = None
        for _ in range(CP.ORBIT_TRIES):
            a = rng.uniform(0, 2 * math.pi)
            gap = min_gap_to_siblings(anchor, a)
            if best is None or gap > best[0]:
                best = (gap, a)
            if gap >= CP.MIN_ORBIT_ANGLE:
                break
        angle = best[1]
        x = ax + T * math.cos(angle); y = ay + T * math.sin(angle)
        step = 0
        while overlaps(r, x, y) and step < 60:
            step += 1; nudges += 1
            x = ax + (T + step) * math.cos(angle); y = ay + (T + step) * math.sin(angle)
        pos[r] = (x, y)
        orbits[anchor].append(angle)
        placed.append(r); remaining.discard(r)
    return pos, nudges


def place_grouped(comps, nbrs, sem):
    import symbols as SYM
    L = {r: SYM.symbol_side_len(comps[r]["device"], CP.side_len(len(comps[r]["pins"]))) for r in comps}
    power_only = TOPO.power_only_devices(sem, comps)
    GROUPS, _ = TOPO.derive_groups(comps, nbrs, power_only)
    group_of = {r: g for g, refs in GROUPS.items() for r in refs}
    lone = sorted(power_only)

    # ---- level 1: pack each group's own devices
    local_pos, group_size, total_nudges = {}, {}, 0
    for g, refs in GROUPS.items():
        w = {r: {p: nbrs[r].get(p, 0.0) for p in refs if p != r and p in nbrs[r]} for r in refs}
        pos, nudges = orbit_place(refs, w, L, lambda a, b: (L[a] + L[b]) / 2 + GAP, seed=CP.SEED, start=None)
        total_nudges += nudges
        for r in refs:
            local_pos[r] = pos[r]
        xs = [pos[r][0] for r in refs]; ys = [pos[r][1] for r in refs]
        span = max(max(xs) - min(xs), max(ys) - min(ys), 0) + max(L[r] for r in refs)
        group_size[g] = span

    # ---- level 2: place the four groups relative to each other, sized by what they actually hold
    gnbrs = {g: {} for g in GROUPS}
    for g1, refs1 in GROUPS.items():
        for g2, refs2 in GROUPS.items():
            if g1 >= g2:
                continue
            w = sum(nbrs[a].get(b, 0.0) for a in refs1 for b in refs2)
            if w > 0:
                gnbrs[g1][g2] = w; gnbrs[g2][g1] = w

    def group_target_distance(a, b):
        return (group_size[a] + group_size[b]) / 2 + GAP

    gpos, gnudges = orbit_place(list(GROUPS), gnbrs, group_size, group_target_distance, seed=CP.SEED)

    # ---- combine: device position = its group's placed centre + its local offset within the group
    pos = {r: (gpos[group_of[r]][0] + local_pos[r][0], gpos[group_of[r]][1] + local_pos[r][1]) for r in group_of}

    # ---- power-only devices: ONE tight cluster, held well clear of the real circuit - not near their IC.
    # First tried "near" placement per the .circuit file's own SOFT_NEARNESS hints, but that put them close
    # enough to the IC that real signal wires had to route around them (seen on C16/C17). Since these devices
    # are never actually wired (power stubs only, by convention), physical proximity to "their" IC buys nothing
    # electrically - being OUT OF THE ROUTER'S WAY matters more, so they all go in one packed group, far enough
    # from the real circuit's bounding box that no signal wire's corridor runs through it.
    FLOAT_GAP = 2          # packed tight against each other
    CLEAR_MARGIN = 3 * GAP  # pushed well clear of the real circuit, not just past the nearest device

    def occupied(x, y, l, extra=0):
        for r2, (px, py) in pos.items():
            l2 = L[r2]
            if x - extra < px + l2 / 2 and px - l2 / 2 - extra < x + l and \
               y - extra < py + l2 / 2 and py - l2 / 2 - extra < y + l:
                return True
        return False

    floaters = sorted(lone)
    ncols = 5
    max_y = max(py + L[r] / 2 for r, (px, py) in pos.items())
    min_x = min(px - L[r] / 2 for r, (px, py) in pos.items())
    ay = max_y + CLEAR_MARGIN
    for i, r in enumerate(floaters):
        l = L[r]
        col, row = i % ncols, i // ncols
        x = min_x + col * (l + FLOAT_GAP) + l / 2
        y = ay + row * (l + FLOAT_GAP) + l / 2
        while occupied(x, y, l, extra=FLOAT_GAP):
            y += l + FLOAT_GAP
        pos[r] = (x, y)

    # ---- squares on the integer grid
    minx = min(pos[r][0] - L[r] / 2 for r in comps); miny = min(pos[r][1] - L[r] / 2 for r in comps)
    squares = {}
    for r in comps:
        x = round(pos[r][0] - L[r] / 2 - minx) + 4
        y = round(pos[r][1] - L[r] / 2 - miny) + 4
        squares[r] = Square(r, x, y, L[r], L[r])
    W = max(s.x + s.width for s in squares.values()) + 5
    H = max(s.y + s.height for s in squares.values()) + 5
    return squares, W, H, total_nudges + gnudges, group_of


def main():
    sem = json.loads((CP.BASE / "semantic.json").read_text())
    comps, nbrs, pin_fn = TOPO.build_weighted_closeness(sem, CP.split_term)
    squares, W, H, nudges, group_of = place_grouped(comps, nbrs, sem)
    print(f"derived groups: {sorted(set(group_of.values()))}")
    print(f"grouped placement: window {W} x {H} (aspect {W/H:.2f}); collision nudges: {nudges}")
    CP.report(comps, nbrs, squares)
    out = CP.render(squares, nbrs, W, H, HERE / "hrng_grouped_density.png")
    print("saved", out, "(copied to ~/picaxe and dump)")


if __name__ == "__main__":
    main()
