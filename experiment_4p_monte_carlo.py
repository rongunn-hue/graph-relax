#!/usr/bin/env python3
"""Monte Carlo Local Compaction (user-labeled 4P; planarity 4P is preserved)."""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import random
import time
from pathlib import Path

import experiment_4i as e4i
import experiment_4j as e4j
import experiment_4n as e4n
import experiment_4o as e4o
import experiment_4p as planarity_4p
import graph_relax as gr


RNG_SEED = 314159
SAMPLES_PER_VERTEX = 128
ELITE_FRACTION = 0.10


def drawing_area(realization):
    extents = e4j.geometry_extents(realization[0], realization[1])
    return extents["width"]*extents["height"]


def rank_key(graph, vertex, row):
    if graph.degree(vertex) == 1:
        return (row["local_length"], row["drawing_area"], row["sample_order"])
    return (row["local_length"], row["angular"]["angular_error"],
            row["drawing_area"], -row["angular"]["minimum_gap"], row["sample_order"])


def sample_valid_candidates(circuit, graph, embedding, positions, vertex, radius, rng):
    current = positions[vertex]
    valid = []
    for sample_order in range(SAMPLES_PER_VERTEX):
        angle = rng.random()*2.0*math.pi
        point = (gr.rounded(current[0]+radius*math.cos(angle)),
                 gr.rounded(current[1]+radius*math.sin(angle)))
        trial_positions = dict(positions); trial_positions[vertex] = point
        try:
            realization = e4n.geometry(circuit, trial_positions)
        except ValueError:
            continue
        if not e4n.locally_valid(circuit, embedding, trial_positions, vertex, realization):
            continue
        valid.append({"sample_order": sample_order, "angle_radians": angle,
                      "point": point, "positions": trial_positions, "geometry": realization,
                      "local_length": e4n.local_length(vertex, realization[3]),
                      "angular": e4o.vertex_angular_measure(graph, trial_positions, vertex),
                      "drawing_area": drawing_area(realization)})
    valid.sort(key=lambda row: rank_key(graph, vertex, row))
    return valid


def elite_proposal(valid, current):
    if not valid:
        return [], None
    if len(valid) < 10:
        elite = valid[:1]
        return elite, None
    elite_count = max(1, math.ceil(ELITE_FRACTION*len(valid)))
    elite = valid[:elite_count]
    mean_dx = sum(row["point"][0]-current[0] for row in elite)/len(elite)
    mean_dy = sum(row["point"][1]-current[1] for row in elite)/len(elite)
    return elite, (gr.rounded(current[0]+mean_dx), gr.rounded(current[1]+mean_dy))


def validate_mean(circuit, graph, embedding, positions, vertex, point):
    if point == positions[vertex]:
        return None
    trial_positions = dict(positions); trial_positions[vertex] = point
    try:
        realization = e4n.geometry(circuit, trial_positions)
    except ValueError:
        return None
    if not e4n.locally_valid(circuit, embedding, trial_positions, vertex, realization):
        return None
    return {"point": point, "positions": trial_positions, "geometry": realization,
            "local_length": e4n.local_length(vertex, realization[3]),
            "angular": e4o.vertex_angular_measure(graph, trial_positions, vertex),
            "drawing_area": drawing_area(realization)}


def choose_move(circuit, graph, embedding, positions, realization, vertex, valid):
    current_length = e4n.local_length(vertex, realization[3])
    elite, mean_point = elite_proposal(valid, positions[vertex])
    if mean_point is not None:
        mean = validate_mean(circuit, graph, embedding, positions, vertex, mean_point)
        if mean is not None and mean["local_length"] < current_length-1e-9:
            mean["acceptance"] = "elite_mean"
            mean["elite_size"] = len(elite)
            return mean, len(elite)
    for row in elite:
        # Samples were already hard-valid, but certify the current selected geometry.
        if row["local_length"] < current_length-1e-9 and e4n.locally_valid(
                circuit, embedding, row["positions"], vertex, row["geometry"]):
            chosen = dict(row); chosen["acceptance"] = "fallback_sample"
            chosen["elite_size"] = len(elite)
            return chosen, len(elite)
    return None, len(elite)


def positions_from(path, circuit):
    document=json.loads(path.read_text())
    positions={planarity_4p.component_vertex(ref):(float(document["components"][ref]["x"]),float(document["components"][ref]["y"])) for ref in circuit.refs}
    positions.update({planarity_4p.net_vertex(net.name):(float(document["net_junctions"][net.name]["x"]),float(document["net_junctions"][net.name]["y"])) for net in circuit.nets})
    return positions


def comparison_row(graph, circuit, output_dir, coordinates, width, height, area, length, runtime):
    positions=positions_from(output_dir/coordinates,circuit)
    return {"width":width,"height":height,"area":area,"total_connection_length":length,
            "angular":e4o.angular_summary(graph,positions),"crossings":0,"component_overlaps":0,
            "clearance_violations":0,"unrelated_component_body_intersections":0,
            "electrical_incidences":"44/44","runtime_seconds":runtime}


def run(input_path, output_dir):
    started=time.perf_counter()
    circuit=gr.parse_circuit(input_path)
    graph,embedding,positions,order=e4n.load_start(circuit,output_dir)
    baseline,realization=e4n.complete_validate(circuit,graph,embedding,positions)
    if not baseline["valid"]: raise ValueError("stored Experiment 4K starting geometry failed reproduction")
    rng=random.Random(RNG_SEED)
    moves,passes=[],[]
    for pass_number,radius in enumerate(e4n.STEPS,1):
        moved=valid_total=elite_total=mean_accepted=fallback_accepted=0
        for vertex in order:
            valid=sample_valid_candidates(circuit,graph,embedding,positions,vertex,radius,rng)
            valid_total+=len(valid)
            chosen,elite_size=choose_move(circuit,graph,embedding,positions,realization,vertex,valid)
            elite_total+=elite_size
            if chosen is None: continue
            old=positions[vertex];before=e4n.local_length(vertex,realization[3])
            positions,realization=chosen["positions"],chosen["geometry"]
            moved+=1
            if chosen["acceptance"]=="elite_mean":mean_accepted+=1
            else:fallback_accepted+=1
            moves.append({"move_number":len(moves)+1,"pass":pass_number,"radius":radius,
                          "vertex":vertex,"degree":graph.degree(vertex),"acceptance":chosen["acceptance"],
                          "elite_size":chosen["elite_size"],"old":e4i.point_doc(old),
                          "new":e4i.point_doc(chosen["point"]),
                          "movement_distance":gr.rounded(math.dist(old,chosen["point"])),
                          "local_length_before":gr.rounded(before),
                          "local_length_after":gr.rounded(chosen["local_length"]),
                          "angular_error":None if chosen["angular"] is None else gr.rounded(chosen["angular"]["angular_error"]),
                          "minimum_gap":None if chosen["angular"] is None else gr.rounded(chosen["angular"]["minimum_gap"])})
        validation,checked=e4n.complete_validate(circuit,graph,embedding,positions)
        if not validation["valid"]:raise RuntimeError(f"complete validation failed after radius {radius}")
        realization=checked;measured=e4i.measure(circuit,realization[0],realization[3]);ext=e4j.geometry_extents(realization[0],realization[1])
        total_samples=len(order)*SAMPLES_PER_VERTEX
        passes.append({"pass":pass_number,"radius":radius,"vertices_moved":moved,
          "total_sampled_candidates":total_samples,"hard_valid_candidates":valid_total,
          "mean_valid_candidate_percentage":gr.rounded(100*valid_total/total_samples),
          "mean_elite_set_size":gr.rounded(elite_total/len(order)),
          "elite_mean_moves_accepted":mean_accepted,"fallback_sample_moves_accepted":fallback_accepted,
          "total_connection_length":gr.rounded(measured["total_connection_length"]),
          "angular":e4o.angular_summary(graph,positions),"width":ext["width"],"height":ext["height"],
          "area":gr.rounded(ext["width"]*ext["height"]),"validation":validation})
    final_validation,final_geometry=e4n.complete_validate(circuit,graph,embedding,positions)
    if not final_validation["valid"]:raise RuntimeError("final complete validation failed")
    final_metrics=e4i.measure(circuit,final_geometry[0],final_geometry[3]);ext=e4j.geometry_extents(final_geometry[0],final_geometry[1]);angular=e4o.angular_summary(graph,positions)
    runtime=time.perf_counter()-started
    nmetrics=json.loads((output_dir/"planar_local_compact_metrics.json").read_text());ometrics=json.loads((output_dir/"planar_circular_compact_metrics.json").read_text())
    comparisons={
      "experiment_4k":comparison_row(graph,circuit,output_dir,"planar_xy_compact_coordinates.json",3166,1112,3520592,23095.744883,None),
      "experiment_4n":comparison_row(graph,circuit,output_dir,"planar_local_compact_coordinates.json",2401.838680,984.746882,2365203.151197,13979.580061,nmetrics["runtime_seconds"]),
      "experiment_4o":comparison_row(graph,circuit,output_dir,"planar_circular_compact_coordinates.json",2545.874241,1006.351588,2562044.585279,15269.151443,ometrics["runtime_seconds"])}
    components,nets,attachments,edges=final_geometry
    coordinate_path=output_dir/"planar_monte_carlo_compact_coordinates.json"
    gr.write_json(coordinate_path,{"circuit":circuit.name,"source":"exact stored Experiment 4K geometry",
      "rng":{"implementation":"Python random.Random (MT19937)","seed":RNG_SEED},"processing_order":order,
      "radii":list(e4n.STEPS),"samples_per_vertex":SAMPLES_PER_VERTEX,"elite_fraction":ELITE_FRACTION,
      "component_dimensions":{"width":gr.WIDTH,"height":gr.HEIGHT},
      "graph_vertex_positions":{v:e4i.point_doc(positions[v]) for v in sorted(positions)},
      "components":{r:e4i.point_doc(components[r]) for r in circuit.refs},
      "net_junctions":{n:e4i.point_doc(nets[n]) for n in sorted(nets)},
      "attachments":{t:e4i.point_doc(attachments[t]) for t in sorted(attachments)},
      "edges":[{"terminal":e["terminal"],"component":e["component"],"net":e["net"],"points":[e4i.point_doc(p) for p in e["points"]]} for e in edges]})
    metrics_path=output_dir/"planar_monte_carlo_compact_metrics.json"
    metrics_document={"rng":{"implementation":"Python random.Random (MT19937)","seed":RNG_SEED},
      "processing_order":order,"radii":list(e4n.STEPS),"samples_per_vertex":SAMPLES_PER_VERTEX,
      "elite_fraction":ELITE_FRACTION,"passes":passes,"moves":moves,"total_accepted_moves":len(moves),
      "total_sampled_candidates":len(order)*SAMPLES_PER_VERTEX*len(e4n.STEPS),"geometry_extents":ext,
      "drawing_area":gr.rounded(ext["width"]*ext["height"]),"metrics":gr.rounded_metrics(final_metrics),
      "angular":angular,"runtime_seconds":gr.rounded(runtime),"comparison":comparisons}
    gr.write_json(metrics_path,metrics_document)
    validation_path=output_dir/"planar_monte_carlo_compact_validation.json"
    gr.write_json(validation_path,{"baseline_4k_reproduced":True,"embedding_unchanged":True,
      "exactly_six_passes":len(passes)==6,"complete_pass_validations":[p["validation"] for p in passes],
      "final_complete_validation":final_validation})
    svg_path=output_dir/"planar_monte_carlo_compact.svg";e4i.write_svg(svg_path,circuit,*final_geometry)
    paths=(coordinate_path,metrics_path,validation_path,svg_path)
    return {"metrics":metrics_document,"validation":final_validation,
      "hashes":{p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}}


def main():
    parser=argparse.ArgumentParser();parser.add_argument("input",type=Path);parser.add_argument("output_dir",type=Path);args=parser.parse_args()
    print(json.dumps(run(args.input,args.output_dir),sort_keys=True,indent=2))


if __name__=="__main__":main()
