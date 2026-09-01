#!/usr/bin/env python3
"""Experiment 4Q (tent): simultaneous fixed-cable radial contraction of 4K."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import time
from pathlib import Path

import experiment_4i as e4i
import experiment_4j as e4j
import experiment_4n as e4n
import experiment_4o as e4o
import graph_relax as gr


FRACTIONS = (0.00, 0.10, 0.20, 0.30, 0.40, 0.50,
             0.60, 0.70, 0.80, 0.90, 0.95, 0.98)
HEIGHT_TOLERANCE = 0.01
MAX_BISECTION_ITERATIONS = 32
CENTER_TOLERANCE = 1e-7


def tent_model(initial_positions):
    vertices = sorted(initial_positions)
    center = (sum(initial_positions[v][0] for v in vertices)/len(vertices),
              sum(initial_positions[v][1] for v in vertices)/len(vertices))
    cables, directions, center_vertices = {}, {}, []
    for vertex in vertices:
        dx = initial_positions[vertex][0]-center[0]
        dy = initial_positions[vertex][1]-center[1]
        length = math.hypot(dx, dy)
        cables[vertex] = length
        if length <= CENTER_TOLERANCE:
            center_vertices.append(vertex)
        else:
            directions[vertex] = (dx/length, dy/length)
    nonzero = [(length, vertex) for vertex, length in cables.items()
               if vertex not in center_vertices]
    minimum, owner = min(nonzero, key=lambda row: (row[0], row[1]))
    return center, cables, directions, center_vertices, minimum, owner


def positions_at_height(initial_positions, center, cables, directions, center_vertices, height):
    result = {}
    centered = set(center_vertices)
    for vertex in sorted(initial_positions):
        if vertex in centered:
            result[vertex] = (gr.rounded(center[0]), gr.rounded(center[1]))
            continue
        length = cables[vertex]
        if not 0.0 <= height < length:
            raise ValueError(f"height {height} is outside cable limit for {vertex}")
        planar_radius = math.sqrt(max(0.0, length*length-height*height))
        unit = directions[vertex]
        result[vertex] = (gr.rounded(center[0]+unit[0]*planar_radius),
                          gr.rounded(center[1]+unit[1]*planar_radius))
    return result


def evaluate(circuit, graph, embedding, initial_positions, model, height, fraction):
    center, cables, directions, center_vertices, _, _ = model
    positions = positions_at_height(initial_positions, center, cables, directions,
                                   center_vertices, height)
    validation, realization = e4n.complete_validate(circuit, graph, embedding, positions)
    if realization is None:
        metrics = angular = extents = None
        record = {"height": gr.rounded(height), "fraction_of_geometric_limit": gr.rounded(fraction),
                  "width": None, "drawing_height": None, "area": None,
                  "total_connection_length": None,
                  "global_minimum_incident_angular_gap": None,
                  "mean_per_vertex_minimum_angular_gap": None,
                  "mean_angular_error": None,
                  "different_net_crossings": None, "component_overlaps": None,
                  "clearance_violations": None,
                  "unrelated_component_body_intersections": None,
                  "electrical_incidence_exact": False,
                  "cyclic_order_orientation": None, "valid": False}
    else:
        metrics = e4i.measure(circuit, realization[0], realization[3])
        angular = e4o.angular_summary(graph, positions)
        extents = e4j.geometry_extents(realization[0], realization[1])
        record = {"height": gr.rounded(height), "fraction_of_geometric_limit": gr.rounded(fraction),
                  "width": extents["width"], "drawing_height": extents["height"],
                  "area": gr.rounded(extents["width"]*extents["height"]),
                  "total_connection_length": gr.rounded(metrics["total_connection_length"]),
                  "global_minimum_incident_angular_gap": angular["global_minimum_incident_angular_gap"],
                  "mean_per_vertex_minimum_angular_gap": angular["mean_per_vertex_minimum_angular_gap"],
                  "mean_angular_error": angular["mean_angular_error"],
                  "different_net_crossings": validation["different_net_crossings"],
                  "component_overlaps": validation["component_overlaps"],
                  "clearance_violations": validation["clearance_violations"],
                  "unrelated_component_body_intersections": validation["unrelated_component_body_intersections"],
                  "electrical_incidence_exact": validation["electrical_incidence_exact"],
                  "cyclic_order_orientation": validation["all_vertex_cyclic_order_orientation"],
                  "valid": validation["valid"]}
    return record, positions, validation, realization, metrics, angular, extents


def failure_reasons(record):
    fields = (("different_net_crossings", "different-net crossings"),
              ("component_overlaps", "component overlaps"),
              ("clearance_violations", "clearance violations"),
              ("unrelated_component_body_intersections", "unrelated body intersections"))
    reasons = [f"{label}: {record[key]}" for key, label in fields
               if record.get(key) not in (None, 0)]
    if not record.get("electrical_incidence_exact", False):
        reasons.append("electrical incidence mismatch or unrealizable geometry")
    if record.get("cyclic_order_orientation") not in ("same", "global_reflection"):
        reasons.append("cyclic-order inconsistency")
    return reasons


def diagnostic_svg(normal_path, diagnostic_path, center, positions):
    content = normal_path.read_text()
    guides = ['<g stroke="#4b82b8" stroke-width="0.7" opacity="0.25">']
    for vertex in sorted(positions):
        point = positions[vertex]
        guides.append(f'<line x1="{center[0]:.3f}" y1="{center[1]:.3f}" x2="{point[0]:.3f}" y2="{point[1]:.3f}"/>')
    guides.append('</g>')
    guides.append(f'<g><circle cx="{center[0]:.3f}" cy="{center[1]:.3f}" r="7" fill="#e31a1c" stroke="white" stroke-width="2"/><text x="{center[0]+10:.3f}" y="{center[1]-10:.3f}" font-family="monospace" font-size="12" fill="#e31a1c">C</text></g>')
    marker = '<g stroke="#777" stroke-width="1" fill="none">'
    content = content.replace(marker, "\n".join(guides)+"\n"+marker, 1)
    diagnostic_path.write_text(content, encoding="utf-8")


def write_csv(path, records):
    fields = ["height", "fraction_of_geometric_limit", "width", "drawing_height", "area",
              "total_connection_length", "global_minimum_incident_angular_gap",
              "mean_per_vertex_minimum_angular_gap", "mean_angular_error",
              "different_net_crossings", "component_overlaps", "clearance_violations",
              "unrelated_component_body_intersections", "electrical_incidence_exact",
              "cyclic_order_orientation", "valid"]
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        for record in records:
            writer.writerow({field: record[field] for field in fields})


def percent_change(new, old): return gr.rounded(100*(new-old)/old)


def run(input_path, output_dir):
    started = time.perf_counter()
    circuit = gr.parse_circuit(input_path)
    graph, embedding, initial_positions, order = e4n.load_start(circuit, output_dir)
    baseline, baseline_geometry = e4n.complete_validate(circuit, graph, embedding, initial_positions)
    if not baseline["valid"]: raise ValueError("stored Experiment 4K starting geometry failed reproduction")
    model = tent_model(initial_positions)
    center, cables, directions, center_vertices, minimum, owner = model
    coarse, payloads = [], {}
    for fraction in FRACTIONS:
        height = fraction*minimum
        payload = evaluate(circuit, graph, embedding, initial_positions, model, height, fraction)
        coarse.append(payload[0]); payloads[height] = payload
    positive_valid = [record for record in coarse if record["height"] > 0 and record["valid"]]
    if not positive_valid:
        selected = payloads[0.0]
        selected_record, positions, validation, realization, metrics, angular, extents = selected
        first_positive = next(record for record in coarse if record["height"] > 0)
        components,nets,attachments,edges=realization
        coordinate_path=output_dir/"planar_tent_compact_coordinates.json"
        gr.write_json(coordinate_path,{"circuit":circuit.name,"source":"exact stored Experiment 4K geometry",
          "status":"NO_POSITIVE_COARSE_HEIGHT_VALID",
          "center":{"x":gr.rounded(center[0]),"y":gr.rounded(center[1])},
          "minimum_nonzero_cable_length":gr.rounded(minimum),"minimum_cable_vertex":owner,
          "center_vertices":center_vertices,"selected_height":0.0,"selected_height_fraction":0.0,
          "component_dimensions":{"width":gr.WIDTH,"height":gr.HEIGHT},
          "graph_vertex_positions":{v:e4i.point_doc(positions[v]) for v in sorted(positions)},
          "components":{r:e4i.point_doc(components[r]) for r in circuit.refs},
          "net_junctions":{n:e4i.point_doc(nets[n]) for n in sorted(nets)},
          "attachments":{t:e4i.point_doc(attachments[t]) for t in sorted(attachments)},
          "edges":[{"terminal":e["terminal"],"component":e["component"],"net":e["net"],"points":[e4i.point_doc(p) for p in e["points"]]} for e in edges]})
        metrics_path=output_dir/"planar_tent_compact_metrics.json"
        metrics_document={"status":"NO_POSITIVE_COARSE_HEIGHT_VALID",
          "center":{"x":gr.rounded(center[0]),"y":gr.rounded(center[1])},
          "minimum_nonzero_cable_length":gr.rounded(minimum),"minimum_cable_vertex":owner,
          "coarse_fractions":list(FRACTIONS),"coarse_sweep":coarse,"selected_height":0.0,
          "selected_height_fraction":0.0,"bisection_performed":False,"bisection_iterations":0,
          "geometry_extents":extents,"drawing_area":gr.rounded(extents["width"]*extents["height"]),
          "metrics":gr.rounded_metrics(metrics),"angular":angular,
          "limiting_height":first_positive["height"],"limiting_constraint":failure_reasons(first_positive)}
        gr.write_json(metrics_path,metrics_document)
        validation_path=output_dir/"planar_tent_compact_validation.json"
        gr.write_json(validation_path,{"status":"NO_POSITIVE_COARSE_HEIGHT_VALID",
          "baseline_4k_reproduced":True,"embedding_unchanged":True,
          "direct_from_original_at_every_height":True,"coarse_validations":coarse,
          "selected_height":0.0,"selected_complete_validation":validation,
          "first_positive_coarse_failure":first_positive})
        svg_path=output_dir/"planar_tent_compact.svg";e4i.write_svg(svg_path,circuit,*realization)
        diagnostic_path=output_dir/"planar_tent_compact_diagnostic.svg";diagnostic_svg(svg_path,diagnostic_path,center,positions)
        csv_path=output_dir/"planar_tent_sweep.csv";write_csv(csv_path,coarse)
        paths=(coordinate_path,metrics_path,validation_path,svg_path,diagnostic_path,csv_path)
        return {"runtime_seconds":gr.rounded(time.perf_counter()-started),"metrics":metrics_document,
          "validation":validation,"hashes":{p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}}
    highest = max(positive_valid, key=lambda record: record["height"])
    low = highest["height"]
    higher_invalid = sorted((record for record in coarse
                             if record["height"] > low and not record["valid"]),
                            key=lambda record: record["height"])
    high = (higher_invalid[0]["height"] if higher_invalid else
            (1.0-1e-9)*minimum)
    low_payload = evaluate(circuit, graph, embedding, initial_positions, model,
                           low, low/minimum)
    high_payload = evaluate(circuit, graph, embedding, initial_positions, model,
                            high, high/minimum)
    if high_payload[0]["valid"]:
        low, low_payload = high, high_payload
        limiting_payload = None
        iterations = 0
    else:
        limiting_payload = high_payload
        iterations = 0
        while high-low > HEIGHT_TOLERANCE and iterations < MAX_BISECTION_ITERATIONS:
            middle = (low+high)/2.0
            payload = evaluate(circuit, graph, embedding, initial_positions, model,
                               middle, middle/minimum)
            iterations += 1
            if payload[0]["valid"]:
                low, low_payload = middle, payload
            else:
                high, limiting_payload = middle, payload
    selected_record, positions, validation, realization, metrics, angular, extents = low_payload
    if not validation["valid"]: raise RuntimeError("selected tent geometry is invalid")
    baseline_metrics=e4i.measure(circuit,baseline_geometry[0],baseline_geometry[3]);baseline_extents=e4j.geometry_extents(baseline_geometry[0],baseline_geometry[1]);baseline_angular=e4o.angular_summary(graph,initial_positions)
    comparison={"experiment_4k":{"width":baseline_extents["width"],"height":baseline_extents["height"],
      "area":gr.rounded(baseline_extents["width"]*baseline_extents["height"]),
      "total_connection_length":gr.rounded(baseline_metrics["total_connection_length"]),"angular":baseline_angular},
      "percent_change":{"width":percent_change(extents["width"],baseline_extents["width"]),
       "height":percent_change(extents["height"],baseline_extents["height"]),
       "area":percent_change(extents["width"]*extents["height"],baseline_extents["width"]*baseline_extents["height"]),
       "total_connection_length":percent_change(metrics["total_connection_length"],baseline_metrics["total_connection_length"]),
       "global_minimum_incident_angular_gap":percent_change(angular["global_minimum_incident_angular_gap"],baseline_angular["global_minimum_incident_angular_gap"]),
       "mean_per_vertex_minimum_angular_gap":percent_change(angular["mean_per_vertex_minimum_angular_gap"],baseline_angular["mean_per_vertex_minimum_angular_gap"]),
       "mean_angular_error":percent_change(angular["mean_angular_error"],baseline_angular["mean_angular_error"])}}
    components,nets,attachments,edges=realization
    coordinate_path=output_dir/"planar_tent_compact_coordinates.json"
    gr.write_json(coordinate_path,{"circuit":circuit.name,"source":"exact stored Experiment 4K geometry",
      "center":{"x":gr.rounded(center[0]),"y":gr.rounded(center[1])},
      "minimum_nonzero_cable_length":gr.rounded(minimum),"minimum_cable_vertex":owner,
      "center_vertices":center_vertices,"selected_height":gr.rounded(low),
      "selected_height_fraction":gr.rounded(low/minimum),"component_dimensions":{"width":gr.WIDTH,"height":gr.HEIGHT},
      "graph_vertex_positions":{v:e4i.point_doc(positions[v]) for v in sorted(positions)},
      "components":{r:e4i.point_doc(components[r]) for r in circuit.refs},
      "net_junctions":{n:e4i.point_doc(nets[n]) for n in sorted(nets)},
      "attachments":{t:e4i.point_doc(attachments[t]) for t in sorted(attachments)},
      "edges":[{"terminal":e["terminal"],"component":e["component"],"net":e["net"],"points":[e4i.point_doc(p) for p in e["points"]]} for e in edges]})
    metrics_path=output_dir/"planar_tent_compact_metrics.json"
    metrics_document={"center":{"x":gr.rounded(center[0]),"y":gr.rounded(center[1])},
      "minimum_nonzero_cable_length":gr.rounded(minimum),"minimum_cable_vertex":owner,
      "coarse_fractions":list(FRACTIONS),"coarse_sweep":coarse,
      "selected_height":gr.rounded(low),"selected_height_fraction":gr.rounded(low/minimum),
      "bisection_iterations":iterations,"height_tolerance":HEIGHT_TOLERANCE,
      "geometry_extents":extents,"drawing_area":gr.rounded(extents["width"]*extents["height"]),
      "metrics":gr.rounded_metrics(metrics),"angular":angular,"comparison":comparison,
      "limiting_height":None if limiting_payload is None else limiting_payload[0]["height"],
      "limiting_constraint":[] if limiting_payload is None else failure_reasons(limiting_payload[0])}
    gr.write_json(metrics_path,metrics_document)
    validation_path=output_dir/"planar_tent_compact_validation.json"
    gr.write_json(validation_path,{"baseline_4k_reproduced":True,"embedding_unchanged":True,
      "direct_from_original_at_every_height":True,"coarse_validations":[p[0] for p in [payloads[k] for k in payloads]],
      "selected_height":gr.rounded(low),"selected_complete_validation":validation,
      "limiting_evaluation":None if limiting_payload is None else limiting_payload[0]})
    svg_path=output_dir/"planar_tent_compact.svg";e4i.write_svg(svg_path,circuit,*realization)
    diagnostic_path=output_dir/"planar_tent_compact_diagnostic.svg";diagnostic_svg(svg_path,diagnostic_path,center,positions)
    csv_path=output_dir/"planar_tent_sweep.csv";write_csv(csv_path,coarse)
    paths=(coordinate_path,metrics_path,validation_path,svg_path,diagnostic_path,csv_path)
    return {"runtime_seconds":gr.rounded(time.perf_counter()-started),"metrics":metrics_document,
      "validation":validation,"hashes":{p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}}


def main():
    parser=argparse.ArgumentParser();parser.add_argument("input",type=Path);parser.add_argument("output_dir",type=Path);args=parser.parse_args()
    print(json.dumps(run(args.input,args.output_dir),sort_keys=True,indent=2))
if __name__=="__main__":main()
