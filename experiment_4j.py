#!/usr/bin/env python3
"""Experiment 4J: exhaustive integer uniform-scale sweep of stored 4I geometry."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import networkx as nx

import experiment_4i as e4i
import experiment_4p as e4p
import graph_relax as gr


SCALES = tuple(range(72, 0, -1))
SOURCE_TOLERANCE = 2e-5


def load_stored(circuit, output_dir):
    coordinates = json.loads((output_dir / "planar_node_initial_coordinates.json").read_text())
    embedding_document = json.loads((output_dir / "node_planar_embedding.json").read_text())
    validation = json.loads((output_dir / "planar_node_initial_validation.json").read_text())
    if coordinates["scale"] != 72 or validation["selected_scale"] != 72:
        raise ValueError("stored 4I baseline scale is not 72")
    graph = e4p.incidence_graph(circuit)
    if e4p.graph_document(graph) != embedding_document["graph"]:
        raise ValueError("stored 4I graph differs from source incidence graph")
    embedding = nx.PlanarEmbedding()
    embedding.set_data(embedding_document["cyclic_neighbor_order_clockwise"])
    embedding.check_structure()
    components = {ref: (coordinates["components"][ref]["x"] / 72.0,
                        coordinates["components"][ref]["y"] / 72.0) for ref in circuit.refs}
    nets = {net.name: (coordinates["net_junctions"][net.name]["x"] / 72.0,
                       coordinates["net_junctions"][net.name]["y"] / 72.0) for net in circuit.nets}
    return coordinates, embedding_document, validation, graph, embedding, components, nets


def geometry_at_scale(circuit, unscaled_components, unscaled_nets, scale):
    components = {ref: (point[0] * scale, point[1] * scale)
                  for ref, point in unscaled_components.items()}
    nets = {name: (point[0] * scale, point[1] * scale)
            for name, point in unscaled_nets.items()}
    attachments, edges = {}, []
    for net in circuit.nets:
        for terminal in net.terminals:
            ref = terminal.rsplit(".", 1)[0]
            attachment = e4i.ray_rectangle_attachment(components[ref], nets[net.name])
            attachments[terminal] = attachment
            edges.append({"terminal": terminal, "component": ref, "net": net.name,
                          "points": [attachment, nets[net.name]]})
    return components, nets, attachments, edges


def geometry_extents(components, nets):
    left = min([point[0] for point in nets.values()] +
               [point[0] - gr.WIDTH / 2 for point in components.values()])
    right = max([point[0] for point in nets.values()] +
                [point[0] + gr.WIDTH / 2 for point in components.values()])
    top = min([point[1] for point in nets.values()] +
              [point[1] - gr.HEIGHT / 2 for point in components.values()])
    bottom = max([point[1] for point in nets.values()] +
                 [point[1] + gr.HEIGHT / 2 for point in components.values()])
    return {"min_x": gr.rounded(left), "min_y": gr.rounded(top),
            "max_x": gr.rounded(right), "max_y": gr.rounded(bottom),
            "width": gr.rounded(right-left), "height": gr.rounded(bottom-top)}


def sweep(circuit, graph, embedding, unscaled_components, unscaled_nets):
    records, geometries = [], {}
    for scale in SCALES:
        geometry = geometry_at_scale(circuit, unscaled_components, unscaled_nets, scale)
        validation = e4i.validate_geometry(circuit, graph, embedding, *geometry)
        metrics = e4i.measure(circuit, geometry[0], geometry[3])
        fully_valid = (validation["component_overlaps"] == 0 and
                       validation["clearance_violations"] == 0 and
                       validation["different_net_crossings"] == 0 and
                       validation["unrelated_component_body_intersections"] == 0 and
                       validation["electrical_incidence_exact"] and
                       validation["expected_incidence_count"] == 44 and
                       validation["realized_incidence_count"] == 44 and
                       validation["cyclic_order_orientation"] in ("same", "global_reflection"))
        records.append({"scale": scale,
                        "total_connection_length": gr.rounded(metrics["total_connection_length"]),
                        "component_overlaps": validation["component_overlaps"],
                        "clearance_violations": validation["clearance_violations"],
                        "different_net_crossings": validation["different_net_crossings"],
                        "unrelated_component_body_intersections": validation["unrelated_component_body_intersections"],
                        "electrical_incidence_exact": validation["electrical_incidence_exact"],
                        "expected_incidences": validation["expected_incidence_count"],
                        "realized_incidences": validation["realized_incidence_count"],
                        "cyclic_order_orientation": validation["cyclic_order_orientation"],
                        "crowding": gr.rounded(metrics["energy"]["crowding"]),
                        "energy": {key: gr.rounded(value) for key, value in metrics["energy"].items()},
                        "objective": gr.rounded(metrics["objective"]),
                        "fully_valid": fully_valid})
        geometries[scale] = (geometry, validation, metrics)
    valid_scales = [record["scale"] for record in records if record["fully_valid"]]
    smallest = min(valid_scales)
    first_invalid = next((record for record in records if record["scale"] == smallest-1), None)
    return records, geometries, smallest, first_invalid


def verify_baseline(circuit, stored_validation, graph, embedding,
                    unscaled_components, unscaled_nets, output_dir):
    geometry = geometry_at_scale(circuit, unscaled_components, unscaled_nets, 72)
    validation = e4i.validate_geometry(circuit, graph, embedding, *geometry)
    for key in ("component_overlaps", "clearance_violations", "different_net_crossings",
                "unrelated_component_body_intersections"):
        if validation[key] != 0 or validation[key] != stored_validation[key]:
            raise ValueError(f"stored 4I baseline discrepancy: {key}")
    if not validation["electrical_incidence_exact"] or validation["realized_incidence_count"] != 44:
        raise ValueError("stored 4I incidence baseline discrepancy")
    expected_metrics = json.loads((output_dir / "planar_node_initial_metrics.json").read_text())["metrics"]
    actual_metrics = e4i.measure(circuit, geometry[0], geometry[3])
    if abs(actual_metrics["total_connection_length"] - expected_metrics["total_connection_length"]) > SOURCE_TOLERANCE:
        raise ValueError("stored 4I connection length discrepancy")
    return geometry, validation, actual_metrics


def coordinate_document(circuit, scale, geometry, validation):
    components, nets, attachments, edges = geometry
    return {"circuit": circuit.name, "source": "exact stored Experiment 4I point realization",
            "scale": scale, "component_dimensions": {"width": gr.WIDTH, "height": gr.HEIGHT},
            "components": {ref: e4i.point_doc(components[ref]) for ref in circuit.refs},
            "net_junctions": {name: e4i.point_doc(nets[name]) for name in sorted(nets)},
            "attachments": {terminal: e4i.point_doc(attachments[terminal]) for terminal in sorted(attachments)},
            "edges": [{"terminal": edge["terminal"], "component": edge["component"],
                       "net": edge["net"], "points": [e4i.point_doc(point) for point in edge["points"]]}
                      for edge in edges],
            "cyclic_order_orientation": validation["cyclic_order_orientation"]}


def failure_reasons(record):
    labels = (("component_overlaps", "component overlaps"),
              ("clearance_violations", "clearance violations"),
              ("different_net_crossings", "different-net crossings"),
              ("unrelated_component_body_intersections", "unrelated body intersections"))
    reasons = [f"{label}: {record[key]}" for key, label in labels if record[key] != 0]
    if not record["electrical_incidence_exact"]:
        reasons.append("electrical incidence mismatch")
    if record["cyclic_order_orientation"] not in ("same", "global_reflection"):
        reasons.append("cyclic-order inconsistency")
    return reasons


def run(input_path, output_dir):
    circuit = gr.parse_circuit(input_path)
    stored_coordinates, embedding_document, stored_validation, graph, embedding, uc, un = load_stored(circuit, output_dir)
    baseline_geometry, baseline_validation, baseline_metrics = verify_baseline(
        circuit, stored_validation, graph, embedding, uc, un, output_dir)
    records, geometries, smallest, first_invalid = sweep(circuit, graph, embedding, uc, un)
    compact_geometry, compact_validation, compact_metrics = geometries[smallest]
    sweep_path = output_dir / "planar_scale_sweep.json"
    gr.write_json(sweep_path, {"fixed_scale_sequence": list(SCALES), "tested_all_scales": True,
                               "records": records, "smallest_fully_valid_scale": smallest,
                               "first_invalid_lower_scale": first_invalid,
                               "first_invalid_reasons": failure_reasons(first_invalid) if first_invalid else []})
    coordinate_path = output_dir / "planar_scale_compact_coordinates.json"
    gr.write_json(coordinate_path, coordinate_document(circuit, smallest, compact_geometry, compact_validation))
    original_extents = geometry_extents(baseline_geometry[0], baseline_geometry[1])
    compact_extents = geometry_extents(compact_geometry[0], compact_geometry[1])
    metrics_document = {"original_4i": gr.rounded_metrics(baseline_metrics),
                        "experiment_4j": gr.rounded_metrics(compact_metrics),
                        "original_scale": 72, "smallest_fully_valid_scale": smallest,
                        "scale_reduction_percent": gr.rounded(100*(72-smallest)/72),
                        "connection_length_reduction_percent": gr.rounded(
                            100*(baseline_metrics["total_connection_length"]-compact_metrics["total_connection_length"])/baseline_metrics["total_connection_length"]),
                        "original_geometry_extents": original_extents,
                        "compacted_geometry_extents": compact_extents,
                        "first_invalid_lower_scale": first_invalid["scale"] if first_invalid else None,
                        "first_invalid_reasons": failure_reasons(first_invalid) if first_invalid else [],
                        "measurement_only": True}
    metrics_path = output_dir / "planar_scale_compact_metrics.json"; gr.write_json(metrics_path, metrics_document)
    validation_document = {"baseline_reproduced": True, "embedding_unchanged": True,
                           "incidence_graph_unchanged": True, "component_dimensions_unchanged": True,
                           "attachment_algorithm": "Experiment 4I ray/rectangle boundary intersection",
                           "selected_scale": smallest, "fully_valid": True,
                           "validation": compact_validation,
                           "first_invalid_lower_scale": first_invalid,
                           "first_invalid_reasons": failure_reasons(first_invalid) if first_invalid else []}
    validation_path = output_dir / "planar_scale_compact_validation.json"; gr.write_json(validation_path, validation_document)
    svg_path = output_dir / "planar_scale_compact.svg"; e4i.write_svg(svg_path, circuit, *compact_geometry)
    return {"smallest_fully_valid_scale": smallest, "first_invalid_lower_scale": first_invalid["scale"],
            "metrics": metrics_document,
            "hashes": {path.name: hashlib.sha256(path.read_bytes()).hexdigest()
                       for path in (sweep_path, coordinate_path, metrics_path, validation_path, svg_path)}}


def main():
    parser=argparse.ArgumentParser();parser.add_argument("input",type=Path);parser.add_argument("output_dir",type=Path);args=parser.parse_args()
    print(json.dumps(run(args.input,args.output_dir),sort_keys=True,indent=2))


if __name__=="__main__": main()
