#!/usr/bin/env python3
"""Derive functional grouping from the circuit's own data instead of hand-copying block labels - the fix for the
R2 misplacement (it was hand-assigned to "Amplification" by a label in the source file, even though its actual
closeness is an exact 50/50 tie between C1 and U1; nothing about that tie said which side it should join).

Two real signals, both used together:
  1. Pin-function metadata already in the .circuit file's source_metadata (e.g. U1 pin 2 = "IN- A", pin 3 =
     "IN+ A"). A net touching an op-amp's IN- or OUT pin is structurally a FEEDBACK connection (closes the loop
     around that op-amp); a net touching IN+ is an INPUT/bias connection (ties to whatever drives that op-amp,
     not to the op-amp's own feedback network). Feedback-net ties get a weight boost; nothing else does.
  2. Community detection (networkx's greedy modularity maximization) run on the resulting weighted closeness
     graph - not a hand-picked box per device. A device like R2, still tied evenly between two neighbours after
     the boost, is decided by which side's WHOLE structure it fits into best, not by which label someone wrote
     in a source file.
Also reports articulation points (real graph-theoretic bridge devices) as a sanity check on the result.
"""
import json, re
import networkx as nx
from networkx.algorithms.community import greedy_modularity_communities

FEEDBACK_MARKERS = ("IN-", "OUT")     # closes the loop around an op-amp; IN+ deliberately excluded
FEEDBACK_BOOST = 3.0

CIRCUIT_PATH = "/home/rgunn/Chat-Projects/deterministic_schematic/circuit_examples/hrng.circuit"


def load_pin_functions():
    """{(ref, pin): function name} straight from the .circuit file's own source_metadata, e.g. ("U1","2") ->
    "IN- A". Devices with plain numbered pins (resistors, caps...) simply have no entry."""
    text = open(CIRCUIT_PATH).read()
    for raw in text.splitlines():
        if raw.startswith("source_metadata "):
            meta = json.loads(raw.split(" ", 1)[1].strip()[1:-1])
            break
    else:
        return {}
    out = {}
    for p in meta["normalized"]["parts"]:
        for pin, info in p.get("terminal_metadata", {}).items():
            name = info.get("name", pin)
            if name != pin:                      # only real function names, not a name equal to the bare number
                out[(p["ref"], pin)] = name
    return out


def build_weighted_closeness(sem, split_term):
    """Same base closeness as before (every net, weight 1/(k-1) per pair), plus a feedback-net boost from pin
    function metadata. Returns (comps, nbrs) exactly like closeness_placer.build_closeness, weights refined."""
    comps = {c["ref"]: c for c in sem["electrical"]["components"]}
    pin_fn = load_pin_functions()
    weight = {}
    for net in sem["electrical"]["nets"]:
        terms = [split_term(t) for t in net["terminals"]]
        refs = sorted({r for r, p in terms})
        k = len(net["terminals"])
        if k < 2:
            continue
        is_feedback = any(any(m in pin_fn.get((r, p), "") for m in FEEDBACK_MARKERS) for r, p in terms)
        w = (FEEDBACK_BOOST if is_feedback else 1.0) / (k - 1)
        for i in range(len(refs)):
            for j in range(i + 1, len(refs)):
                key = (refs[i], refs[j])
                weight[key] = weight.get(key, 0.0) + w
    nbrs = {r: {} for r in comps}
    for (a, b), w in weight.items():
        nbrs[a][b] = w; nbrs[b][a] = w
    return comps, nbrs, pin_fn


def derive_groups(comps, nbrs, exclude):
    """Community detection on the weighted closeness graph, excluding power-only devices (`exclude`). Returns
    {group_name: [refs]}, named after each community's highest-weighted-degree member for readability."""
    G = nx.Graph()
    members = [r for r in comps if r not in exclude]
    G.add_nodes_from(members)
    for a in members:
        for b, w in nbrs[a].items():
            if b in G and a < b:
                G.add_edge(a, b, weight=w)
    communities = list(greedy_modularity_communities(G, weight="weight"))
    groups = {}
    for com in communities:
        total_w = {r: sum(nbrs[r].get(p, 0) for p in com if p != r) for r in com}
        leader = max(com, key=lambda r: (total_w[r], r))
        groups[leader] = sorted(com)
    return groups, G


def articulation_report(G):
    if G.number_of_nodes() < 3 or not nx.is_connected(G):
        return []
    return sorted(nx.articulation_points(G))


def power_only_devices(sem, comps):
    """devices with no signal net at all - every net they touch is a power rail - excluded from grouping."""
    def is_power_only(r):
        nets = [n["name"] for n in sem["electrical"]["nets"] for t in n["terminals"] if t.startswith(r + ".")]
        return all(n in ("+5V", "+9V", "-9V", "0V") for n in nets)
    return {r for r in comps if is_power_only(r)}


def main():
    from pathlib import Path
    import closeness_placer as CP
    sem = json.loads((CP.BASE / "semantic.json").read_text())
    comps, nbrs, pin_fn = build_weighted_closeness(sem, CP.split_term)
    power_only = power_only_devices(sem, comps)
    groups, G = derive_groups(comps, nbrs, power_only)
    print(f"{len(groups)} derived groups ({len(power_only)} power-only devices excluded, placed separately as before):")
    for leader, members in sorted(groups.items()):
        print(f"  [{leader}] {members}")
    arts = articulation_report(G)
    print(f"\narticulation points (structural bridges) in the signal closeness graph: {arts}")
    print("\nR2's group:", next(g for g, m in groups.items() if "R2" in m))
    print("R2's boosted ties:", {k: round(v, 2) for k, v in nbrs["R2"].items() if v > 0})


if __name__ == "__main__":
    main()


def placement_hints():
    """{ref: target_ref} straight from the .circuit file's own top-level `placement '{...}'` directives
    (SOFT_NEARNESS facts, e.g. C13 -> U1). Devices with no such directive are absent from the returned dict."""
    out = {}
    for raw in open(CIRCUIT_PATH).read().splitlines():
        if not raw.startswith("placement "):
            continue
        body = raw[len("placement "):].strip()
        if body[0] == body[-1] == "'":
            body = body[1:-1]
        info = json.loads(body)
        if info.get("kind") == "SOFT_NEARNESS" and info.get("near") and info.get("ref"):
            out[info["ref"]] = info["near"]
    return out
