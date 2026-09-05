#!/usr/bin/env python3
"""One frozen crowding snapshot, one global least-squares displacement."""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import numpy as np

import graph_first_experiment_common_metric as cm
import graph_first_experiment_common_metric_guard as guard
import graph_first_experiment_metric_derived as derived
import graph_first_experiment_planar_region_inventory as inventory
import graph_first_experiment_polygon_area_expansion as polygon
import placement_geometry as geometry
from render_iamp_common_metric_direct_clean import clean_svg, bounds, graph_elements

TOL = 1e-7
CATEGORIES = ("polygon", "node_node", "node_edge", "unrelated_edge_edge", "shared_endpoint_overlay")


def edge(a, b):
    return tuple(sorted((a, b)))


def segment_parameter(point, first, second):
    dx, dy = second[0] - first[0], second[1] - first[1]
    length2 = dx * dx + dy * dy
    if length2 <= TOL * TOL:
        return 0.0
    return max(0.0, min(1.0, ((point[0] - first[0]) * dx + (point[1] - first[1]) * dy) / length2))


def polygon_relations(faces, coordinates):
    """Deficient nearest-boundary pairs, deduplicated by actual vertex pair."""
    strongest = {}
    duplicates = []
    for face in sorted(faces, key=lambda item: item["face_id"]):
        if not (face["is_bounded"] and face["simple_boundary"]):
            continue
        walk = face["ordered_boundary_walk"]
        vertices = sorted(set(walk))
        D = polygon.metrics(walk, coordinates)["average_nearest_boundary_node_distance"]
        seen = set()
        for vertex in vertices:
            other = min((q for q in vertices if q != vertex),
                        key=lambda q: (math.dist(coordinates[vertex], coordinates[q]), q))
            pair = tuple(sorted((vertex, other)))
            if pair in seen:
                continue
            seen.add(pair)
            d = math.dist(coordinates[pair[0]], coordinates[pair[1]])
            if d >= D - TOL:
                continue
            item = {"type": "POLYGON", "category": "polygon", "pair": pair,
                    "face_id": face["face_id"], "source_faces": [face["face_id"]],
                    "affected": list(pair), "d": d, "D": D, "ratio": d / D}
            previous = strongest.get(pair)
            if previous is None:
                strongest[pair] = item
            else:
                keep, drop = (item, previous) if (item["ratio"], item["face_id"]) < (previous["ratio"], previous["face_id"]) else (previous, item)
                keep["source_faces"] = sorted(set(keep["source_faces"] + drop["source_faces"]))
                strongest[pair] = keep
                duplicates.append({"pair": list(pair), "kept_face": keep["face_id"],
                                   "discarded_face": drop["face_id"],
                                   "reason": "same polygon boundary-node geometric pair; strongest normalized demand retained"})
    return sorted(strongest.values(), key=lambda x: (x["ratio"], x["pair"])), duplicates


def active_relations(faces, coordinates, edges):
    detected = cm.corrected_map(faces, coordinates, edges)
    polygons, duplicates = polygon_relations(faces, coordinates)
    relations = {key: [] for key in CATEGORIES}
    relations["polygon"] = polygons
    for key in CATEGORIES[1:]:
        relations[key] = [dict(item, category=key) for item in detected[key] if item["d"] < item["D"] - TOL]
    counter = 1
    for category in CATEGORIES:
        for item in relations[category]:
            item["relation_id"] = f"G{counter:03d}"
            counter += 1
    return relations, detected, duplicates


def add_weight(weights, vertex, value):
    weights[vertex] = weights.get(vertex, 0.0) + value


def relation_witness(relation, coordinates):
    """Return baseline p, q and affine displacement weights for p-q."""
    kind = relation["type"]
    weights = {}
    detail = {}
    if kind in ("POLYGON", "NODE_NODE"):
        a, b = relation.get("pair", relation["affected"])
        p, q = coordinates[a], coordinates[b]
        add_weight(weights, a, 1.0); add_weight(weights, b, -1.0)
        detail = {"p_vertex": a, "q_vertex": b}
    elif kind == "NODE_EDGE":
        v = relation["affected"]["vertex"]
        a, b = relation["affected"]["edge"]
        p = coordinates[v]
        q = tuple(relation["closest_point"])
        t = segment_parameter(q, coordinates[a], coordinates[b])
        add_weight(weights, v, 1.0); add_weight(weights, a, -(1.0 - t)); add_weight(weights, b, -t)
        detail = {"p_vertex": v, "q_edge": [a, b], "q_segment_parameter": t}
    elif kind == "EDGE_EDGE":
        (a, b), (c, d) = relation["affected"]
        p, q = map(tuple, relation["closest_points"])
        s = segment_parameter(p, coordinates[a], coordinates[b])
        t = segment_parameter(q, coordinates[c], coordinates[d])
        add_weight(weights, a, 1.0 - s); add_weight(weights, b, s)
        add_weight(weights, c, -(1.0 - t)); add_weight(weights, d, -t)
        detail = {"p_edge": [a, b], "p_segment_parameter": s,
                  "q_edge": [c, d], "q_segment_parameter": t}
    elif kind == "SHARED_OVERLAY":
        v = relation["shared_vertex"]
        a, b = relation["nonshared_endpoints"]
        s1, s2 = relation["endpoint_to_other_segment_distances"]
        if s1 <= s2 + TOL:
            p = coordinates[a]
            q = cm.base.cp(p, coordinates[v], coordinates[b])
            t = segment_parameter(q, coordinates[v], coordinates[b])
            add_weight(weights, a, 1.0); add_weight(weights, v, -(1.0 - t)); add_weight(weights, b, -t)
            detail = {"responsible_witness": "A_TO_SEGMENT_V_B", "p_vertex": a,
                      "q_edge": [v, b], "q_segment_parameter": t}
        else:
            p = coordinates[b]
            q = cm.base.cp(p, coordinates[v], coordinates[a])
            t = segment_parameter(q, coordinates[v], coordinates[a])
            add_weight(weights, b, 1.0); add_weight(weights, v, -(1.0 - t)); add_weight(weights, a, -t)
            detail = {"responsible_witness": "B_TO_SEGMENT_V_A", "p_vertex": b,
                      "q_edge": [v, a], "q_segment_parameter": t}
    else:
        raise AssertionError(kind)
    return tuple(p), tuple(q), weights, detail


def assemble(relations, coordinates, vertex_order):
    rows, skipped = [], []
    index = {vertex: i for i, vertex in enumerate(vertex_order)}
    for category in CATEGORIES:
        for relation in relations[category]:
            p, q, weights, witness = relation_witness(relation, coordinates)
            dx, dy = p[0] - q[0], p[1] - q[1]
            d = math.hypot(dx, dy)
            if d <= TOL:
                skipped.append({"relation_id": relation["relation_id"], "type": relation["type"],
                                "affected": relation["affected"], "d": d, "D": relation["D"],
                                "reason": "undefined normal and no detector-established deterministic separation direction"})
                continue
            n = (dx / d, dy / d)
            row = np.zeros(2 * len(vertex_order), dtype=float)
            for vertex, weight in weights.items():
                offset = 2 * index[vertex]
                row[offset] += weight * n[0] / relation["D"]
                row[offset + 1] += weight * n[1] / relation["D"]
            b = 1.0 - relation["d"] / relation["D"]
            rows.append({"relation": relation, "p": list(p), "q": list(q), "normal": list(n),
                         "affine_displacement_weights": weights, "witness": witness,
                         "matrix_row": row, "rhs": b})
    A = np.vstack([item["matrix_row"] for item in rows]) if rows else np.empty((0, 2 * len(vertex_order)))
    b = np.array([item["rhs"] for item in rows], dtype=float)
    return A, b, rows, skipped


def cyclic_order(coordinates, edges):
    adjacency = {v: [] for v in coordinates}
    for a, b in edges:
        adjacency[a].append(b); adjacency[b].append(a)
    return {v: sorted(adjacency[v], key=lambda q: (math.atan2(coordinates[q][1] - coordinates[v][1], coordinates[q][0] - coordinates[v][0]), q))
            for v in sorted(adjacency)}


def cyclic_equal(a, b):
    if len(a) != len(b):
        return False
    if len(a) < 2:
        return True
    return any(a == b[i:] + b[:i] for i in range(len(b)))


def embedding_audit(before, after, edges):
    old, new = cyclic_order(before, edges), cyclic_order(after, edges)
    rows = []
    for vertex in sorted(old):
        relation = "same" if cyclic_equal(new[vertex], old[vertex]) else ("reversed" if cyclic_equal(new[vertex], list(reversed(old[vertex]))) else "changed")
        rows.append({"vertex": vertex, "before": old[vertex], "after": new[vertex], "relation": relation})
    changed = [x["vertex"] for x in rows if x["relation"] == "changed"]
    return {"valid": not changed, "changed_vertices": changed, "vertices": rows}


def ratio_summary(crowding_map):
    values = sorted(item["ratio"] for category in CATEGORIES for item in crowding_map[category])
    return {"count": len(values), "minimum": min(values) if values else None,
            "maximum": max(values) if values else None,
            "mean": sum(values) / len(values) if values else None,
            "median": (values[len(values)//2] if len(values) % 2 else (values[len(values)//2-1] + values[len(values)//2]) / 2) if values else None}


def relations_involving(crowding_map, names):
    return [item for category in CATEGORIES for item in crowding_map[category]
            if any(name in str(item.get("affected")) for name in names)]


def tpb_r3_measurement(coordinates):
    shared, first, second = "net:NODE_B", "component:R3", "component:TPB"
    d = guard.overlay_distance(coordinates, shared, first, second)
    D, detail = cm.local_d(coordinates, [shared, first, second])
    return {"relation_type": "SHARED_OVERLAY", "shared_vertex": shared,
            "nonshared_endpoints": [first, second], "d": d, "D": D,
            "ratio": d / D, "local_node_spacing": detail,
            "active": d < D - TOL}


def solve_once(vertices, edges, faces, start):
    vertex_order = sorted(vertices)
    relations, initial_map, duplicates = active_relations(faces, start, edges)
    A, b, row_records, skipped = assemble(relations, start, vertex_order)
    delta, residuals, rank, singular = np.linalg.lstsq(A, b, rcond=None)
    prediction = A @ delta
    residual = prediction - b
    for item, achieved, error in zip(row_records, prediction, residual):
        rel = item["relation"]
        item["audit"] = {"relation_id": rel["relation_id"], "relation_type": rel["type"],
                         "affected": rel["affected"], "baseline_d": rel["d"], "baseline_D": rel["D"],
                         "baseline_r": rel["ratio"], "requested_normalized_improvement": item["rhs"],
                         "predicted_achieved_row_value": float(achieved), "normalized_residual": float(error)}
    final = {vertex: (start[vertex][0] + float(delta[2*i]), start[vertex][1] + float(delta[2*i+1]))
             for i, vertex in enumerate(vertex_order)}
    final_map = cm.corrected_map(faces, final, edges)
    diag = geometry.graph_geometry_diagnostics(final, edges)
    rebuilt = geometry.rebuild_straight_edge_segments(final, edges)
    reconnect = geometry.validate_rebuilt_edge_segments(final, edges, rebuilt)
    displacements = {vertex: {"dx": float(delta[2*i]), "dy": float(delta[2*i+1]),
                              "magnitude": math.hypot(float(delta[2*i]), float(delta[2*i+1]))}
                     for i, vertex in enumerate(vertex_order)}
    requested = len(row_records)
    row_stats = {"fully_met_or_exceeded": sum(float(x) >= float(y) - TOL for x, y in zip(prediction, b)),
                 "partially_met": sum(TOL < float(x) < float(y) - TOL for x, y in zip(prediction, b)),
                 "predicted_to_worsen": sum(float(x) < -TOL for x in prediction)}
    cluster = ["component:U1", "net:-9V", "net:EB", "net:NODE_B"]
    return {"vertex_order": vertex_order, "relations": relations, "initial_map": initial_map,
            "deduplicated_relations": duplicates, "skipped_rows": skipped,
            "A_shape": list(A.shape), "matrix_rank": int(rank),
            "singular_values": [float(x) for x in singular],
            "condition_number_nonzero": float(singular[0] / singular[rank-1]) if rank else None,
            "least_squares_residual_array": [float(x) for x in residuals],
            "mean_dx": float(np.mean(delta[0::2])), "mean_dy": float(np.mean(delta[1::2])),
            "row_audit": [item["audit"] | {"witness": item["witness"], "normal": item["normal"],
                                                   "affine_displacement_weights": item["affine_displacement_weights"]}
                          for item in row_records],
            "residual_statistics": {"l2_norm": float(np.linalg.norm(residual)),
                                    "rms": float(math.sqrt(float(np.mean(residual * residual)))) if requested else 0.0,
                                    "maximum_absolute": float(np.max(np.abs(residual))) if requested else 0.0,
                                    **row_stats},
            "displacements": displacements, "final_coordinates": final,
            "maximum_vertex_displacement": max(x["magnitude"] for x in displacements.values()),
            "rms_vertex_displacement": math.sqrt(sum(x["magnitude"]**2 for x in displacements.values()) / len(displacements)),
            "final_map": final_map, "final_diagnostics": diag, "reconnection": reconnect,
            "embedding_audit": embedding_audit(start, final, edges),
            "central_cluster_pair_distances": {f"{a}--{b}": {"before": math.dist(start[a], start[b]), "after": math.dist(final[a], final[b])}
                                               for i, a in enumerate(cluster) for b in cluster[i+1:]},
            "central_cluster_events_before": relations_involving(initial_map, cluster),
            "central_cluster_events_after": relations_involving(final_map, cluster)}


def diagnostic_svg(vertices, edges, start, final, path):
    both = dict(start)
    both.update({f"final:{v}": p for v, p in final.items()})
    x0, y0, x1, y1 = bounds(both); pad = 75
    z = [f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="{x0-pad} {y0-pad} {x1-x0+2*pad} {y1-y0+2*pad}" width="1500" height="1100">',
         f'<rect x="{x0-pad}" y="{y0-pad}" width="{x1-x0+2*pad}" height="{y1-y0+2*pad}" fill="white"/>',
         f'<text x="{x0-pad+12}" y="{y0-pad+25}" font-size="18" font-weight="bold">IAMP GLOBAL LINEAR — SIMULTANEOUS DISPLACEMENT VECTORS</text>',
         '<defs><marker id="arrow" markerWidth="7" markerHeight="7" refX="6" refY="3.5" orient="auto"><path d="M0,0 L7,3.5 L0,7 z" fill="#d84315"/></marker></defs>']
    z += graph_elements(vertices, edges, start, 10)
    for vertex in sorted(vertices):
        a, b = start[vertex], final[vertex]
        z.append(f'<line x1="{a[0]:.9f}" y1="{a[1]:.9f}" x2="{b[0]:.9f}" y2="{b[1]:.9f}" stroke="#d84315" stroke-width="1.2" marker-end="url(#arrow)" opacity=".8"/>')
        z.append(f'<circle cx="{b[0]:.9f}" cy="{b[1]:.9f}" r="2.8" fill="#d84315"/>')
    z.append(f'<text x="{x0-pad+12}" y="{y1+pad-15}" font-size="12" fill="#d84315">Orange vectors: one minimum-norm global solution, applied simultaneously.</text>')
    z.append('</svg>')
    path.write_text('\n'.join(z) + '\n')


def serializable_result(result):
    copy = dict(result)
    copy["final_coordinates"] = {v: list(p) for v, p in sorted(result["final_coordinates"].items())}
    return copy


def run(graph_path, direct_path, outdir):
    graph = json.loads(graph_path.read_text()); direct = json.loads(direct_path.read_text())
    vertices = {item["id"]: item for item in graph["vertices"]}
    edges = [(item["source"], item["target"]) for item in graph["edges"]]
    endpoint_pairs = guard.edgeset(edges)
    start = {v: tuple(p) for v, p in direct["final_coordinates"].items()}
    valid, start_diag = cm.hard_valid(start, edges)
    if not valid or len(vertices) != 34 or len(edges) != 44:
        raise AssertionError("invalid authoritative post-uncrossing baseline")
    face_list, darts = inventory.enumerate_faces(vertices, edges, start)
    inventory.add_classification_and_adjacency(face_list, edges, darts)
    first = solve_once(vertices, edges, face_list, start)
    second = solve_once(vertices, edges, face_list, start)
    first_serial = serializable_result(first); second_serial = serializable_result(second)
    deterministic = json.dumps(first_serial, sort_keys=True, separators=(",", ":"), allow_nan=False) == json.dumps(second_serial, sort_keys=True, separators=(",", ":"), allow_nan=False)
    if not deterministic:
        raise AssertionError("global linear solve was not deterministic")
    final = first["final_coordinates"]
    before_penalty = derived.global_penalty(first["initial_map"])
    after_penalty = derived.global_penalty(first["final_map"])
    active_counts = {key: len(first["relations"][key]) for key in CATEGORIES}
    counts_before = {key: len(first["initial_map"][key]) for key in CATEGORIES}
    counts_after = {key: len(first["final_map"][key]) for key in CATEGORIES}
    min_before = {key: min((x["ratio"] for x in first["initial_map"][key]), default=None) for key in CATEGORIES}
    min_after = {key: min((x["ratio"] for x in first["final_map"][key]), default=None) for key in CATEGORIES}
    tpb_before = [x for x in first["initial_map"]["shared_endpoint_overlay"] if "component:TPB" in x["nonshared_endpoints"] and x["shared_vertex"] == "net:NODE_B"]
    tpb_after = [x for x in first["final_map"]["shared_endpoint_overlay"] if "component:TPB" in x["nonshared_endpoints"] and x["shared_vertex"] == "net:NODE_B"]
    outdir.mkdir(parents=True, exist_ok=True)
    start_svg = outdir / "iamp_global_linear_start.svg"
    final_svg = outdir / "iamp_global_linear_final.svg"
    diagnostic = outdir / "iamp_global_linear_diagnostic.svg"
    clean_svg(vertices, edges, start, "IAMP GLOBAL LINEAR — START", start_svg)
    clean_svg(vertices, edges, final, "IAMP GLOBAL LINEAR — FINAL", final_svg)
    diagnostic_svg(vertices, edges, start, final, diagnostic)
    report = {"schema": "graph-relax.global-linear-crowding.v1",
              "source_coordinates": f"{direct_path}#final_coordinates",
              "method": "one frozen detector snapshot; normalized linear crowding rows; numpy.linalg.lstsq minimum-norm solution; one simultaneous displacement",
              "starting_validation": {"vertices": 34, "edges": 44, "endpoint_pairs_identical": endpoint_pairs == guard.edgeset(edges),
                                      "proper_crossings": start_diag["proper_unrelated_edge_crossing_count"],
                                      "coincidences": start_diag["coincident_vertex_pair_count"],
                                      "vertex_on_unrelated_edge": start_diag["vertex_on_unrelated_edge_interior_count"]},
              "vertex_order": first["vertex_order"], "unknown_count": 2 * len(vertices),
              "active_relation_counts": active_counts, "deduplicated_relations": first["deduplicated_relations"],
              "skipped_rows": first["skipped_rows"], "matrix": {"shape": first["A_shape"], "rank": first["matrix_rank"],
              "singular_values": first["singular_values"], "condition_number_nonzero": first["condition_number_nonzero"]},
              "solution": {"mean_dx": first["mean_dx"], "mean_dy": first["mean_dy"],
                           "maximum_vertex_displacement": first["maximum_vertex_displacement"],
                           "rms_vertex_displacement": first["rms_vertex_displacement"],
                           "least_squares_residual_array": first["least_squares_residual_array"],
                           "residual_statistics": first["residual_statistics"], "per_vertex": first["displacements"]},
              "row_audit": first["row_audit"],
              "crowding": {"penalty_before": before_penalty, "penalty_after": after_penalty,
                           "violation_count_before": sum(counts_before.values()), "violation_count_after": sum(counts_after.values()),
                           "counts_before": counts_before, "counts_after": counts_after,
                           "minimum_ratio_before": min_before, "minimum_ratio_after": min_after,
                           "ratio_summary_before": ratio_summary(first["initial_map"]),
                           "ratio_summary_after": ratio_summary(first["final_map"])},
              "TPB_NODE_B_R3_overlay": {"direct_measurement_before": tpb_r3_measurement(start),
                                         "direct_measurement_after": tpb_r3_measurement(final),
                                         "all_active_TPB_overlays_before": tpb_before,
                                         "all_active_TPB_overlays_after": tpb_after},
              "J_EB_C2_distance": {"before": math.dist(start["component:J_EB"], start["component:C2"]),
                                    "after": math.dist(final["component:J_EB"], final["component:C2"])},
              "central_cluster_pair_distances": first["central_cluster_pair_distances"],
              "central_cluster_events_before": first["central_cluster_events_before"],
              "central_cluster_events_after": first["central_cluster_events_after"],
              "final_validation": {"vertices": len(vertices), "edges": len(edges),
                                   "endpoint_pairs_identical": endpoint_pairs == guard.edgeset(edges),
                                   "reconnected_edges_valid": first["reconnection"]["valid"],
                                   "proper_crossings": first["final_diagnostics"]["proper_unrelated_edge_crossing_count"],
                                   "crossing_set": first["final_diagnostics"]["proper_unrelated_edge_crossings"],
                                   "coincidences": first["final_diagnostics"]["coincident_vertex_pair_count"],
                                   "coincident_pairs": first["final_diagnostics"]["coincident_vertex_pairs"],
                                   "vertex_on_unrelated_edge": first["final_diagnostics"]["vertex_on_unrelated_edge_interior_count"],
                                   "vertex_on_edge_events": first["final_diagnostics"]["vertices_on_unrelated_edge_interiors"],
                                   "embedding_cyclic_order": first["embedding_audit"]},
              "determinism": {"runs": 2, "identical_serialized_calculation": deterministic},
              "final_coordinates": {v: list(p) for v, p in sorted(final.items())},
              "svg_files": [p.name for p in (start_svg, final_svg, diagnostic)]}
    report_path = outdir / "iamp_global_linear_report.json"
    report_path.write_text(json.dumps(report, indent=2, sort_keys=True, allow_nan=False) + "\n")
    return report


def main():
    parser = argparse.ArgumentParser(); base = Path("output/graph_first")
    parser.add_argument("--graph", type=Path, default=base / "iamp_graph.json")
    parser.add_argument("--direct", type=Path, default=base / "iamp_nearness_rotation_direct_report.json")
    parser.add_argument("--output", type=Path, default=base)
    args = parser.parse_args(); report = run(args.graph, args.direct, args.output)
    print(json.dumps({"active_relation_counts": report["active_relation_counts"], "matrix": report["matrix"],
                      "solution": {k: report["solution"][k] for k in ("mean_dx", "mean_dy", "maximum_vertex_displacement", "rms_vertex_displacement", "residual_statistics")},
                      "crowding": report["crowding"], "TPB": report["TPB_NODE_B_R3_overlay"],
                      "J_EB_C2": report["J_EB_C2_distance"], "final_validation": report["final_validation"],
                      "determinism": report["determinism"]}, indent=2))


if __name__ == "__main__":
    main()
