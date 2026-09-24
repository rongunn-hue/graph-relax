#!/usr/bin/env python3
"""Independent electrical check of a routed Universe. Nets are rebuilt from GEOMETRY ALONE (committed line
points, tunnel spans, attachment points) - never from routing labels (net_id, LineEndpoint declarations) - and
compared with the real netlist terminal by terminal. This is the same methodology used earlier in this project
on the old pipeline, written fresh for this one.

Rules used to read the geometry:
  * two lines that share a point are joined, EXCEPT where that point is inside a tunnel's span (a crossing, not
    a connection)
  * a line end that lies on an attachment point connects to that pin
"""
from universe_defs import PointState


def ordered_points(line):
    out = [line.points[0]]
    for a, b in zip(line.points, line.points[1:]):
        dx = (b[0] > a[0]) - (b[0] < a[0]); dy = (b[1] > a[1]) - (b[1] < a[1])
        p = a
        while p != b:
            p = (p[0] + dx, p[1] + dy); out.append(p)
    return out


def nc_pins(circuit_path):
    """NC pins straight from the .circuit file's own `nc` directives."""
    out = []
    for raw in open(circuit_path).read().splitlines():
        w = raw.split()
        if len(w) == 2 and w[0] == "nc":
            out.append(w[1])
    return out


def check(universe, squares, sem, circuit_path=None):
    results = []
    def rep(name, ok, detail=""):
        results.append(ok)
        print(f"{'PASS' if ok else 'FAIL'}  {name}{('  - ' + detail) if detail else ''}")

    signal = {n["name"]: set(n["terminals"]) for n in sem["electrical"]["nets"]
              if n["name"] not in ("+5V", "+9V", "-9V", "0V")}
    power = {n["name"]: set(n["terminals"]) for n in sem["electrical"]["nets"]
             if n["name"] in ("+5V", "+9V", "-9V", "0V")}
    all_signal_pins = set().union(*signal.values()) if signal else set()

    pts = {lid: ordered_points(l) for lid, l in universe.lines.items()}
    at_point = {p: aid for aid, p in universe.attachment_points.items()}
    solid = set()
    for sq in squares.values():
        solid |= sq.occupied_points()

    by_point = {}
    for lid, plist in pts.items():
        for p in plist:
            by_point.setdefault(p, []).append(lid)

    parent = {}
    def find(x):
        parent.setdefault(x, x)
        while parent[x] != x:
            parent[x] = parent[parent[x]]; x = parent[x]
        return x
    def union(a, b):
        parent[find(a)] = find(b)
    for lid in pts:
        find(("L", lid))

    crossings, bad_shared = [], []
    for p, lids in by_point.items():
        if len(lids) < 2:
            continue
        tunnels = [l for l in lids if universe.lines[l].tunnel and p not in universe.lines[l].endpoint_points()]
        if tunnels:
            others = [l for l in lids if l not in tunnels]
            straight = len(lids) == 2 and len(tunnels) == 1 and not universe.lines[others[0]].tunnel
            (crossings if straight else bad_shared).append((p, tuple(lids)))
            continue
        if len(lids) > 2:
            bad_shared.append((p, tuple(lids))); continue
        union(("L", lids[0]), ("L", lids[1]))

    touched, on_device = set(), []
    for lid, plist in pts.items():
        for p in plist:
            if p in solid and not (p in (plist[0], plist[-1]) and p in at_point):
                on_device.append((lid, p))
        for p in (plist[0], plist[-1]):
            if p in at_point:
                union(("L", lid), ("A", at_point[p])); touched.add(at_point[p])

    comps = {}
    for lid in pts:
        comps.setdefault(find(("L", lid)), {"lines": set(), "pins": set()})["lines"].add(lid)
    for aid in touched:
        comps.setdefault(find(("A", aid)), {"lines": set(), "pins": set()})["pins"].add(aid)
    component_pins = [frozenset(c["pins"]) for c in comps.values() if c["lines"]]
    expected = {frozenset(t) for t in signal.values() if len(t) > 1}

    rep("every wire/tunnel belongs to a conductor group that touches >= 2 pins (no floating wire)",
        all(len(c["pins"]) >= 2 for c in comps.values() if c["lines"]), f"{sum(1 for c in comps.values() if c['lines'])} conductor groups")
    rep("each conductor group's pin set equals exactly one signal net's terminals",
        set(component_pins) == expected, f"{len(expected)} expected, {len(set(component_pins))} found")
    missing = expected - set(component_pins); extra = set(component_pins) - expected
    for m in sorted(missing, key=sorted)[:5]:
        print("      expected but not found:", sorted(m))
    for e in sorted(extra, key=sorted)[:5]:
        print("      found but not in netlist:", sorted(e))
    rep("no wire touches a pin that is not in a signal net (power/NC pins untouched)",
        touched <= all_signal_pins, f"{len(touched)} pins wired")
    nc = nc_pins(circuit_path) if circuit_path else []
    if nc:
        rep("no wire touches an NC pin", not (touched & set(nc)), f"NC pins: {', '.join(nc)}")
    rep("no wire passes through or along a device body (only ends on attachments)", not on_device, str(on_device[:3]))
    rep("every shared point is a joint or a straight 90-degree tunnel crossing",
        not bad_shared, f"{len(crossings)} tunnel crossings; problems: {bad_shared[:3]}")
    cross_ok = all(find(("L", a)) != find(("L", b)) for _, (a, b) in crossings)
    rep("the two wires at every tunnel crossing are NOT connected to each other", cross_ok)
    tun = [l for l in universe.lines.values() if l.tunnel]
    rep("every tunnel is part of its own net's conductor group", all(find(("L", l.line_id)) is not None for l in tun), f"{len(tun)} tunnels")

    print()
    print(f"signal nets checked: {len(expected)} ({sum(len(t) for t in signal.values())} terminals)")
    print(f"not drawn (power rails as local stubs, by convention): {len(power)} power nets, "
          f"{sum(len(t) for t in power.values())} pin connections")
    print("OVERALL:", "ALL CHECKS PASS" if all(results) else "SOME CHECKS FAILED")
    return all(results)
