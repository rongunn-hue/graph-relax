#!/usr/bin/env python3
"""One normalized global crowding/preservation matrix operator for IAMP."""
from __future__ import annotations

import argparse
import html
import json
import math
from pathlib import Path

import numpy as np

import graph_first_experiment_common_metric as cm
import graph_first_experiment_common_metric_guard as guard
import graph_first_experiment_global_linear as linear
import graph_first_experiment_planar_region_inventory as inventory
import placement_geometry as geometry
from render_iamp_common_metric_direct_clean import clean_svg, bounds, graph_elements

TOL = 1e-7
CATEGORIES = linear.CATEGORIES


def detector_snapshot(faces, coordinates, edges):
    relations, raw, duplicates = linear.active_relations(faces, coordinates, edges)
    return relations, raw, duplicates


def construct_matrices(relations, coordinates, edges, vertex_order):
    index = {vertex: i for i, vertex in enumerate(vertex_order)}
    c_rows, normals, audits, skipped = [], [], [], []
    for category in CATEGORIES:
        for relation in relations[category]:
            p, q, weights, witness = linear.relation_witness(relation, coordinates)
            dx, dy = p[0] - q[0], p[1] - q[1]
            d = math.hypot(dx, dy)
            if d <= TOL:
                skipped.append({"relation_id": relation["relation_id"], "relation_type": relation["type"],
                                "affected": relation["affected"], "reason": "undefined witness normal"})
                continue
            normal = (dx / d, dy / d)
            row = np.zeros(len(vertex_order))
            for vertex, weight in weights.items():
                row[index[vertex]] += weight / relation["D"]
            c_rows.append(row); normals.append(normal)
            audits.append({"relation_id": relation["relation_id"], "relation_type": relation["type"],
                           "category": category, "affected": relation["affected"], "d": relation["d"],
                           "D": relation["D"], "r": relation["ratio"], "witness_p": list(p),
                           "witness_q": list(q), "normal": list(normal), "affine_difference_weights": weights,
                           "witness": witness})
    C = np.vstack(c_rows) if c_rows else np.empty((0, len(vertex_order)))
    N = np.asarray(normals, dtype=float).reshape((-1, 2))
    q_rows = []
    edge_audit = []
    for a, b in sorted((linear.edge(a, b) for a, b in edges)):
        length = math.dist(coordinates[a], coordinates[b])
        row = np.zeros(len(vertex_order)); row[index[a]] = 1.0 / length; row[index[b]] = -1.0 / length
        q_rows.append(row); edge_audit.append({"edge": [a, b], "baseline_length": length})
    Q = np.vstack(q_rows)
    return C, N, Q, audits, edge_audit, skipped


def solve_operator(vertices, edges, faces, start):
    order = sorted(vertices)
    X = np.asarray([start[v] for v in order], dtype=float)
    relations, raw_map, duplicates = detector_snapshot(faces, start, edges)
    C, N, Q, row_audit, edge_audit, skipped = construct_matrices(relations, start, edges, order)
    M = C.T @ C + Q.T @ Q
    B = C.T @ N + Q.T @ Q @ X
    ones = np.ones((len(order), 1))
    K = np.block([[M, ones], [ones.T, np.zeros((1, 1))]])
    rhs = np.vstack([B, np.sum(X, axis=0, keepdims=True)])
    solved = np.linalg.solve(K, rhs)
    Xp = solved[:-1]
    coordinates = {v: (float(Xp[i, 0]), float(Xp[i, 1])) for i, v in enumerate(order)}
    singular_M = np.linalg.svd(M, compute_uv=False)
    singular_K = np.linalg.svd(K, compute_uv=False)
    rank_M = int(np.linalg.matrix_rank(M)); rank_K = int(np.linalg.matrix_rank(K))
    before_fit = float(np.sum((C @ X - N) ** 2))
    fit = float(np.sum((C @ Xp - N) ** 2))
    preserve = float(np.sum((Q @ (Xp - X)) ** 2))
    for item, achieved in zip(row_audit, C @ Xp):
        item["solved_witness_vector_normalized"] = [float(x) for x in achieved]
        item["target_unit_vector"] = item["normal"]
        item["vector_residual"] = [float(achieved[i] - item["normal"][i]) for i in range(2)]
        item["squared_residual"] = sum(x*x for x in item["vector_residual"])
    delta = Xp - X
    displacement = {v: {"dx": float(delta[i, 0]), "dy": float(delta[i, 1]),
                        "magnitude": float(np.linalg.norm(delta[i]))} for i, v in enumerate(order)}
    top = sorted(({"vertex": v, **x} for v, x in displacement.items()), key=lambda x: (-x["magnitude"], x["vertex"]))[:10]
    final_relations, final_raw, final_duplicates = detector_snapshot(faces, coordinates, edges)
    diag = geometry.graph_geometry_diagnostics(coordinates, edges)
    rebuilt = geometry.rebuild_straight_edge_segments(coordinates, edges)
    reconnect = geometry.validate_rebuilt_edge_segments(coordinates, edges, rebuilt)
    original_centroid = np.mean(X, axis=0); solved_centroid = np.mean(Xp, axis=0)
    return {"vertex_order": order, "relations": relations, "raw_map": raw_map,
            "duplicates": duplicates, "final_relations": final_relations, "final_raw_map": final_raw,
            "final_duplicates": final_duplicates, "skipped": skipped, "C": C, "N": N, "Q": Q,
            "M": M, "K": K, "row_audit": row_audit, "edge_audit": edge_audit,
            "rank_M": rank_M, "rank_K": rank_K, "singular_M": singular_M,
            "singular_K": singular_K, "coordinates": coordinates, "displacement": displacement,
            "top_displacements": top, "maximum_displacement": top[0]["magnitude"],
            "maximum_displacement_vertex": top[0]["vertex"],
            "rms_displacement": float(math.sqrt(np.mean(np.sum(delta*delta, axis=1)))),
            "energy": {"crowding_fit_before": before_fit, "crowding_fit_after": fit,
                       "preservation": preserve, "total_after": fit + preserve},
            "original_centroid": [float(x) for x in original_centroid],
            "solved_centroid": [float(x) for x in solved_centroid],
            "centroid_max_error": float(np.max(np.abs(solved_centroid-original_centroid))),
            "diag": diag, "reconnect": reconnect,
            "embedding": linear.embedding_audit(start, coordinates, edges)}


def category_stats(relations):
    result = {}
    for category in CATEGORIES:
        ratios = sorted(x["ratio"] for x in relations[category])
        result[category] = {"count": len(ratios), "minimum_r": min(ratios) if ratios else None,
                            "mean_r": sum(ratios)/len(ratios) if ratios else None,
                            "median_r": (ratios[len(ratios)//2] if len(ratios)%2 else (ratios[len(ratios)//2-1]+ratios[len(ratios)//2])/2) if ratios else None,
                            "penalty": sum(1-r for r in ratios)}
    return result


def total_penalty(relations):
    return sum(max(0.0, 1.0-x["ratio"]) for category in CATEGORIES for x in relations[category])


def tpb_audit(coordinates, relations):
    specific = linear.tpb_r3_measurement(coordinates)
    active = [x for x in relations["shared_endpoint_overlay"] if "component:TPB" in x["nonshared_endpoints"]]
    return {"NODE_B_TPB_vs_NODE_B_R3": specific, "all_active_TPB_overlays": active}


def central_distances(coordinates):
    pairs = [("component:U1", "net:-9V"), ("component:U1", "net:EB"),
             ("component:U1", "net:NODE_B"), ("net:-9V", "net:EB"),
             ("net:-9V", "net:NODE_B"), ("net:EB", "net:NODE_B")]
    return {f"{a}--{b}": math.dist(coordinates[a], coordinates[b]) for a, b in pairs}


def diagnostic_svg(vertices, edges, coordinates, rows, path):
    x0, y0, x1, y1 = bounds(coordinates); pad = 70
    colors = {"polygon": "#ef6c00", "node_node": "#7b1fa2", "node_edge": "#00838f",
              "unrelated_edge_edge": "#d32f2f", "shared_endpoint_overlay": "#2e7d32"}
    z = [f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="{x0-pad} {y0-pad} {x1-x0+2*pad} {y1-y0+2*pad}" width="1500" height="1100">',
         f'<rect x="{x0-pad}" y="{y0-pad}" width="{x1-x0+2*pad}" height="{y1-y0+2*pad}" fill="white"/>',
         f'<text x="{x0-pad+12}" y="{y0-pad+25}" font-size="18" font-weight="bold">IAMP GLOBAL MATRIX — FROZEN ACTIVE RELATION WITNESSES</text>']
    z += graph_elements(vertices, edges, coordinates, 10)
    for item in rows:
        p, q = item["witness_p"], item["witness_q"]; color = colors[item["category"]]
        z.append(f'<line x1="{p[0]}" y1="{p[1]}" x2="{q[0]}" y2="{q[1]}" stroke="{color}" stroke-width="1" opacity=".38"/>')
        mx, my = (p[0]+q[0])/2, (p[1]+q[1])/2
        z.append(f'<text x="{mx+2}" y="{my-2}" font-size="5.5" fill="{color}" opacity=".8">{html.escape(item["relation_id"])}</text>')
    z.append(f'<text x="{x0-pad+12}" y="{y1+pad-12}" font-size="10">Orange polygon · purple node-node · teal node-edge · red unrelated edges · green shared overlay</text>')
    z.append('</svg>'); path.write_text('\n'.join(z)+'\n')


def serial_coordinates(result):
    return {v: [float(f"{p[0]:.17g}"), float(f"{p[1]:.17g}")] for v, p in sorted(result["coordinates"].items())}


def run(graph_path, direct_path, outdir):
    graph = json.loads(graph_path.read_text()); direct = json.loads(direct_path.read_text())
    vertices = {x["id"]: x for x in graph["vertices"]}; edges = [(x["source"], x["target"]) for x in graph["edges"]]
    endpoints = guard.edgeset(edges); start = {v: tuple(p) for v, p in direct["final_coordinates"].items()}
    valid, start_diag = cm.hard_valid(start, edges)
    if not valid or len(vertices) != 34 or len(edges) != 44:
        raise AssertionError("invalid authoritative zero-crossing baseline")
    faces, darts = inventory.enumerate_faces(vertices, edges, start)
    inventory.add_classification_and_adjacency(faces, edges, darts)
    first = solve_operator(vertices, edges, faces, start)
    second = solve_operator(vertices, edges, faces, start)
    deterministic = serial_coordinates(first) == serial_coordinates(second)
    if not deterministic:
        raise AssertionError("matrix construction/solve not deterministic at serialization precision")
    final = first["coordinates"]
    before_stats = category_stats(first["relations"]); after_stats = category_stats(first["final_relations"])
    outdir.mkdir(parents=True, exist_ok=True)
    paths = [outdir/"iamp_global_matrix_start.svg", outdir/"iamp_global_matrix_final.svg", outdir/"iamp_global_matrix_diagnostic.svg"]
    clean_svg(vertices, edges, start, "IAMP GLOBAL MATRIX — START", paths[0])
    clean_svg(vertices, edges, final, "IAMP GLOBAL MATRIX — FINAL", paths[1])
    diagnostic_svg(vertices, edges, start, first["row_audit"], paths[2])
    report = {"schema": "graph-relax.global-matrix-crowding.v1",
              "source_coordinates": f"{direct_path}#final_coordinates",
              "method": "one frozen normalized C/Q construction; one centroid-constrained KKT solve; one simultaneous plant",
              "starting_validation": {"vertices": 34, "edges": 44, "endpoint_pairs_identical": endpoints == guard.edgeset(edges),
                                      "proper_crossings": start_diag["proper_unrelated_edge_crossing_count"],
                                      "coincidences": start_diag["coincident_vertex_pair_count"],
                                      "vertex_on_unrelated_edge": start_diag["vertex_on_unrelated_edge_interior_count"]},
              "active_relation_counts": {k: len(first["relations"][k]) for k in CATEGORIES},
              "deduplicated_relations": first["duplicates"], "skipped_relations": first["skipped"],
              "matrices": {"C_dimensions": list(first["C"].shape), "Q_dimensions": list(first["Q"].shape),
                           "M_dimensions": list(first["M"].shape), "M_rank": first["rank_M"],
                           "M_nullity": len(vertices)-first["rank_M"],
                           "M_singular_values": [float(x) for x in first["singular_M"]],
                           "M_condition_number_nonzero": float(first["singular_M"][0]/first["singular_M"][first["rank_M"]-1]),
                           "KKT_dimensions": list(first["K"].shape), "KKT_rank": first["rank_K"],
                           "KKT_singular_values": [float(x) for x in first["singular_K"]],
                           "KKT_condition_number": float(first["singular_K"][0]/first["singular_K"][-1])},
              "centroid_constraint": {"original": first["original_centroid"], "solved": first["solved_centroid"],
                                      "maximum_absolute_error": first["centroid_max_error"]},
              "energy": first["energy"], "row_audit": first["row_audit"],
              "edge_preservation_rows": first["edge_audit"],
              "displacement": {"per_vertex": first["displacement"], "maximum": first["maximum_displacement"],
                               "maximum_vertex": first["maximum_displacement_vertex"],
                               "rms": first["rms_displacement"], "ten_largest": first["top_displacements"]},
              "crowding": {"before": before_stats, "after": after_stats,
                           "total_penalty_before": total_penalty(first["relations"]),
                           "total_penalty_after": total_penalty(first["final_relations"]),
                           "total_active_before": sum(x["count"] for x in before_stats.values()),
                           "total_active_after": sum(x["count"] for x in after_stats.values())},
              "TPB_audit": {"before": tpb_audit(start, first["relations"]), "after": tpb_audit(final, first["final_relations"])},
              "J_EB_C2_distance": {"before": math.dist(start["component:J_EB"], start["component:C2"]),
                                    "after": math.dist(final["component:J_EB"], final["component:C2"])},
              "central_distances": {"before": central_distances(start), "after": central_distances(final)},
              "final_validation": {"vertices": 34, "edges": 44, "endpoint_pairs_identical": endpoints == guard.edgeset(edges),
                                   "reconnected_edges_valid": first["reconnect"]["valid"],
                                   "proper_crossings": first["diag"]["proper_unrelated_edge_crossing_count"],
                                   "crossing_set": first["diag"]["proper_unrelated_edge_crossings"],
                                   "coincidences": first["diag"]["coincident_vertex_pair_count"],
                                   "coincident_pairs": first["diag"]["coincident_vertex_pairs"],
                                   "vertex_on_unrelated_edge": first["diag"]["vertex_on_unrelated_edge_interior_count"],
                                   "vertex_on_edge_events": first["diag"]["vertices_on_unrelated_edge_interiors"],
                                   "embedding_cyclic_order": first["embedding"]},
              "determinism": {"runs": 2, "coordinates_identical_at_17_digit_serialization": deterministic},
              "final_coordinates": {v: list(p) for v, p in sorted(final.items())},
              "svg_files": [p.name for p in paths]}
    (outdir/"iamp_global_matrix_report.json").write_text(json.dumps(report, indent=2, sort_keys=True, allow_nan=False)+'\n')
    return report


def main():
    parser = argparse.ArgumentParser(); base = Path("output/graph_first")
    parser.add_argument("--graph", type=Path, default=base/"iamp_graph.json")
    parser.add_argument("--direct", type=Path, default=base/"iamp_nearness_rotation_direct_report.json")
    parser.add_argument("--output", type=Path, default=base)
    args = parser.parse_args(); report = run(args.graph, args.direct, args.output)
    print(json.dumps({k: report[k] for k in ("active_relation_counts", "matrices", "centroid_constraint", "energy", "displacement", "crowding", "TPB_audit", "J_EB_C2_distance", "central_distances", "final_validation", "determinism")}, indent=2))


if __name__ == "__main__":
    main()
