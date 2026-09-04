#!/usr/bin/env python3
"""Local planarization using only fixed-length single-edge endpoint rotation."""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import networkx as nx

import placement_geometry as geometry


COARSE_DEGREES = 1
REFINE_TOLERANCE_DEGREES = 1e-7
REFINE_ITERATIONS = 40


def load_start(graph_path: Path, coordinate_path: Path, mds_report_path: Path):
    graph_doc = json.loads(graph_path.read_text(encoding="utf-8"))
    coordinate_doc = json.loads(coordinate_path.read_text(encoding="utf-8"))
    mds_report = json.loads(mds_report_path.read_text(encoding="utf-8"))
    graph = nx.Graph()
    graph.add_nodes_from((item["id"], {"type": item["type"], "name": item["name"]})
                         for item in graph_doc["vertices"])
    graph.add_edges_from((item["source"], item["target"]) for item in graph_doc["edges"])
    positions = {item["id"]: (float(item["x"]), float(item["y"]))
                 for item in coordinate_doc["vertices"]}
    if set(positions) != set(graph) or graph.number_of_nodes() != 34 or graph.number_of_edges() != 44:
        raise ValueError("MDS start does not match authoritative graph")
    spacing = mds_report["minimum_spacing_enforcement"]
    diagnostics = geometry.graph_geometry_diagnostics(positions, graph.edges)
    if spacing["status"] != "SATISFIED":
        raise ValueError("MDS start did not satisfy permanent spacing pipeline")
    if diagnostics["coincident_vertex_pair_count"] or \
            diagnostics["minimum_node_distance"] + geometry.GEOMETRY_TOLERANCE < geometry.MIN_NODE_DISTANCE:
        raise ValueError("MDS start violates permanent node spacing")
    if diagnostics["vertex_on_unrelated_edge_interior_count"]:
        raise ValueError("MDS start has a vertex on an unrelated edge")
    return graph, positions, diagnostics, spacing


def crossing_key(crossing: dict):
    return (tuple(crossing["edge_a"]), tuple(crossing["edge_b"]))


def rotated_position(moving_point, pivot_point, signed_degrees: float):
    radius = math.dist(moving_point, pivot_point)
    angle = math.atan2(moving_point[1] - pivot_point[1],
                       moving_point[0] - pivot_point[0]) + math.radians(signed_degrees)
    return (pivot_point[0] + radius * math.cos(angle),
            pivot_point[1] + radius * math.sin(angle))


def evaluate(graph, positions, moving, pivot, signed_degrees, target_key, stats=None):
    if stats is not None:
        stats["candidate_evaluations"] = stats.get("candidate_evaluations", 0) + 1
    candidate = dict(positions)
    candidate[moving] = rotated_position(positions[moving], positions[pivot], signed_degrees)
    diagnostics = geometry.graph_geometry_diagnostics(candidate, graph.edges)
    admissible = (
        diagnostics["coincident_vertex_pair_count"] == 0
        and diagnostics["vertex_on_unrelated_edge_interior_count"] == 0
        and diagnostics["minimum_node_distance"] + geometry.GEOMETRY_TOLERANCE
            >= geometry.MIN_NODE_DISTANCE
    )
    keys = {crossing_key(item) for item in diagnostics["proper_unrelated_edge_crossings"]}
    target_removed = target_key not in keys
    return candidate, diagnostics, admissible and target_removed


def incremental_diagnostics(graph, positions, candidate, moving, current_diagnostics):
    """Exactly update complete diagnostics for one changed vertex."""
    edges = sorted((min(a, b), max(a, b)) for a, b in graph.edges)
    affected = {edge for edge in edges if moving in edge}
    crossing_map = {crossing_key(item): item for item in
                    current_diagnostics["proper_unrelated_edge_crossings"]
                    if tuple(item["edge_a"]) not in affected and tuple(item["edge_b"]) not in affected}
    for first, second in sorted(affected):
        for third, fourth in edges:
            other = (third, fourth)
            if other == (first, second) or {first, second} & {third, fourth}:
                continue
            edge_a, edge_b = sorted(((first, second), other))
            key = (edge_a, edge_b)
            if geometry.proper_segment_crossing(candidate[edge_a[0]], candidate[edge_a[1]],
                                                candidate[edge_b[0]], candidate[edge_b[1]]):
                crossing_map[key] = {"edge_a": list(edge_a), "edge_b": list(edge_b)}
            else:
                crossing_map.pop(key, None)
    distances = [(node, math.dist(candidate[moving], candidate[node]))
                 for node in sorted(graph) if node != moving]
    unchanged_minimum = min(
        (math.dist(positions[first], positions[second])
         for index, first in enumerate(sorted(graph)) if first != moving
         for second in sorted(graph)[index + 1:] if second != moving),
        default=math.inf)
    minimum = min(unchanged_minimum, min(distance for _, distance in distances))
    coincidences = [[min(moving, node), max(moving, node)] for node, distance in distances
                    if distance <= geometry.GEOMETRY_TOLERANCE]
    on_edges = []
    for first, second in edges:
        if moving not in {first, second} and geometry.point_on_segment_interior(
                candidate[moving], candidate[first], candidate[second]):
            on_edges.append({"vertex": moving, "edge": [first, second]})
    for first, second in sorted(affected):
        for node in sorted(graph):
            if node not in {first, second} and geometry.point_on_segment_interior(
                    candidate[node], candidate[first], candidate[second]):
                on_edges.append({"vertex": node, "edge": [first, second]})
    on_edges.sort(key=lambda item: (item["vertex"], item["edge"]))
    crossings = [crossing_map[key] for key in sorted(crossing_map)]
    return {
        "minimum_node_distance": minimum,
        "coincident_vertex_pair_count": len(coincidences),
        "coincident_vertex_pairs": coincidences,
        "proper_unrelated_edge_crossing_count": len(crossings),
        "proper_unrelated_edge_crossings": crossings,
        "vertex_on_unrelated_edge_interior_count": len(on_edges),
        "vertices_on_unrelated_edge_interiors": on_edges,
    }


def refine_candidate(graph, positions, moving, pivot, signed_degrees, target_key,
                     desired_crossings, stats=None):
    sign = 1.0 if signed_degrees > 0 else -1.0
    upper = abs(signed_degrees)
    lower = max(0.0, upper - COARSE_DEGREES)
    best = None
    for _ in range(REFINE_ITERATIONS):
        if upper - lower <= REFINE_TOLERANCE_DEGREES:
            break
        middle = (lower + upper) / 2.0
        candidate, diagnostics, useful = evaluate(
            graph, positions, moving, pivot, sign * middle, target_key, stats)
        if useful and diagnostics["proper_unrelated_edge_crossing_count"] <= desired_crossings:
            upper = middle
            best = (candidate, diagnostics, sign * middle)
        else:
            lower = middle
    if best is None:
        candidate, diagnostics, useful = evaluate(
            graph, positions, moving, pivot, signed_degrees, target_key, stats)
        if not useful:
            raise AssertionError("coarse candidate became unusable during refinement")
        return candidate, diagnostics, signed_degrees
    candidate, diagnostics, useful = evaluate(
        graph, positions, moving, pivot, sign * upper, target_key, stats)
    return candidate, diagnostics, sign * upper


def best_move(graph, positions, current_diagnostics, stats=None):
    before_count = current_diagnostics["proper_unrelated_edge_crossing_count"]
    options = []
    evaluation_cache = {}
    signed_angles = [float(value) for magnitude in range(COARSE_DEGREES, 180, COARSE_DEGREES)
                     for value in (-magnitude, magnitude)] + [180.0]
    for crossing in current_diagnostics["proper_unrelated_edge_crossings"]:
        target = crossing_key(crossing)
        endpoints = [(crossing["edge_a"][0], crossing["edge_a"][1]),
                     (crossing["edge_a"][1], crossing["edge_a"][0]),
                     (crossing["edge_b"][0], crossing["edge_b"][1]),
                     (crossing["edge_b"][1], crossing["edge_b"][0])]
        for moving, pivot in endpoints:
            best_coarse = None
            for signed in signed_angles:
                    cache_key = (moving, pivot, signed)
                    if cache_key not in evaluation_cache:
                        if stats is not None:
                            stats["candidate_evaluations"] = stats.get("candidate_evaluations", 0) + 1
                        candidate = dict(positions)
                        candidate[moving] = rotated_position(
                            positions[moving], positions[pivot], signed)
                        candidate_diagnostics = incremental_diagnostics(
                            graph, positions, candidate, moving, current_diagnostics)
                        admissible = (
                            candidate_diagnostics["coincident_vertex_pair_count"] == 0
                            and candidate_diagnostics["vertex_on_unrelated_edge_interior_count"] == 0
                            and candidate_diagnostics["minimum_node_distance"]
                                + geometry.GEOMETRY_TOLERANCE >= geometry.MIN_NODE_DISTANCE
                        )
                        keys = {crossing_key(item) for item in
                                candidate_diagnostics["proper_unrelated_edge_crossings"]}
                        evaluation_cache[cache_key] = (
                            candidate, candidate_diagnostics, admissible, keys)
                    candidate, diagnostics, admissible, keys = evaluation_cache[cache_key]
                    useful = admissible and target not in keys
                    count = diagnostics["proper_unrelated_edge_crossing_count"]
                    if not useful or count >= before_count:
                        continue
                    displacement = math.dist(positions[moving], candidate[moving])
                    key = (count, abs(signed), displacement, signed)
                    if best_coarse is None or key < best_coarse[0]:
                        best_coarse = (key, signed, candidate, diagnostics)
            if best_coarse is None:
                continue
            _, signed, _, coarse_diagnostics = best_coarse
            candidate, diagnostics, refined_signed = refine_candidate(
                graph, positions, moving, pivot, signed, target,
                coarse_diagnostics["proper_unrelated_edge_crossing_count"], stats)
            displacement = math.dist(positions[moving], candidate[moving])
            count = diagnostics["proper_unrelated_edge_crossing_count"]
            options.append({
                "target": target,
                "moving": moving,
                "pivot": pivot,
                "signed_degrees": refined_signed,
                "displacement": displacement,
                "positions": candidate,
                "diagnostics": diagnostics,
                "selection_key": (count, -(before_count - count), abs(refined_signed),
                                  displacement, moving, pivot, refined_signed),
            })
    return min(options, key=lambda item: item["selection_key"]) if options else None


def run_search(graph, initial_positions, stats=None):
    positions = dict(initial_positions)
    diagnostics = geometry.graph_geometry_diagnostics(positions, graph.edges)
    states = [dict(positions)]
    moves = []
    crossing_history = [diagnostics["proper_unrelated_edge_crossings"]]
    while diagnostics["proper_unrelated_edge_crossing_count"]:
        move = best_move(graph, positions, diagnostics, stats)
        if move is None:
            break
        complete = geometry.graph_geometry_diagnostics(move["positions"], graph.edges)
        if complete != move["diagnostics"]:
            raise AssertionError("incremental candidate validation disagrees with complete validator")
        move["diagnostics"] = complete
        before_keys = {crossing_key(item) for item in diagnostics["proper_unrelated_edge_crossings"]}
        after_keys = {crossing_key(item) for item in
                      move["diagnostics"]["proper_unrelated_edge_crossings"]}
        old_length = math.dist(positions[move["moving"]], positions[move["pivot"]])
        new_length = math.dist(move["positions"][move["moving"]],
                               move["positions"][move["pivot"]])
        if not math.isclose(old_length, new_length, rel_tol=1e-12, abs_tol=1e-9):
            raise AssertionError("rotation failed to preserve pivot-edge length")
        record = {
            "iteration": len(moves) + 1,
            "crossing_count_before": diagnostics["proper_unrelated_edge_crossing_count"],
            "targeted_crossing": {"edge_a": list(move["target"][0]),
                                  "edge_b": list(move["target"][1])},
            "moving_vertex": move["moving"],
            "pivot_vertex": move["pivot"],
            "preserved_edge_length_before": old_length,
            "preserved_edge_length_after": new_length,
            "signed_rotation_degrees": move["signed_degrees"],
            "absolute_rotation_degrees": abs(move["signed_degrees"]),
            "moving_vertex_displacement": move["displacement"],
            "crossing_count_after": move["diagnostics"]["proper_unrelated_edge_crossing_count"],
            "crossings_removed": [
                {"edge_a": list(key[0]), "edge_b": list(key[1])}
                for key in sorted(before_keys - after_keys)
            ],
            "crossings_newly_created": [
                {"edge_a": list(key[0]), "edge_b": list(key[1])}
                for key in sorted(after_keys - before_keys)
            ],
            "minimum_node_distance_after": move["diagnostics"]["minimum_node_distance"],
            "start_position": list(positions[move["moving"]]),
            "end_position": list(move["positions"][move["moving"]]),
            "pivot_position": list(positions[move["pivot"]]),
        }
        moves.append(record)
        positions = move["positions"]
        diagnostics = complete
        states.append(dict(positions))
        crossing_history.append(diagnostics["proper_unrelated_edge_crossings"])
    return positions, diagnostics, moves, states, crossing_history


def svg_frame(graph, positions, title, viewbox, highlighted=None):
    min_x, min_y, width, height = viewbox
    parts = [f'<rect x="{min_x}" y="{min_y}" width="{width}" height="{height}" fill="white"/>',
             f'<text x="{min_x+16}" y="{min_y+28}" font-family="sans-serif" font-size="18">{title}</text>']
    for first, second in sorted((min(a, b), max(a, b)) for a, b in graph.edges):
        parts.append(f'<line x1="{positions[first][0]:.6f}" y1="{positions[first][1]:.6f}" '
                     f'x2="{positions[second][0]:.6f}" y2="{positions[second][1]:.6f}" '
                     'stroke="#777" stroke-width="1.5"/>')
    for node in sorted(graph):
        x, y = positions[node]
        color = "#1565c0" if graph.nodes[node]["type"] == "COMPONENT" else "#d84315"
        radius = 6 if node != highlighted else 10
        parts.append(f'<circle cx="{x:.6f}" cy="{y:.6f}" r="{radius}" fill="{color}" '
                     'stroke="white" stroke-width="1"/>')
        parts.append(f'<text x="{x+8:.6f}" y="{y-8:.6f}" font-family="sans-serif" '
                     f'font-size="11">{graph.nodes[node]["name"]}</text>')
    return parts


def viewbox_for(states):
    xs = [point[0] for state in states for point in state.values()]
    ys = [point[1] for state in states for point in state.values()]
    margin = 100.0
    return (min(xs)-margin, min(ys)-margin,
            max(xs)-min(xs)+2*margin, max(ys)-min(ys)+2*margin)


def static_svg(graph, positions, title, viewbox):
    min_x, min_y, width, height = viewbox
    return (f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="{min_x} {min_y} {width} {height}" '
            'width="1400" height="1000">\n' +
            "\n".join(svg_frame(graph, positions, title, viewbox)) + '\n</svg>\n')


def motion_svg(graph, states, moves, viewbox):
    min_x, min_y, width, height = viewbox
    duration = max(1, len(states)) * 2
    parts = [f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="{min_x} {min_y} {width} {height}" '
             'width="1400" height="1000">']
    for index, state in enumerate(states):
        move = moves[index-1] if index else None
        parts.append(f'<g visibility="hidden">')
        parts.append(f'<set attributeName="visibility" to="visible" begin="{index*2}s" dur="2s" '
                     f'repeatCount="indefinite" repeatDur="{duration}s"/>')
        parts.extend(svg_frame(graph, state, f'Endpoint rotation stage {index}', viewbox,
                               move["moving_vertex"] if move else None))
        if move:
            px, py = move["pivot_position"]
            radius = move["preserved_edge_length_before"]
            parts.append(f'<circle cx="{px:.6f}" cy="{py:.6f}" r="{radius:.6f}" fill="none" '
                         'stroke="#8e24aa" stroke-width="2" stroke-dasharray="7 5" opacity="0.65"/>')
            parts.append(f'<circle cx="{px:.6f}" cy="{py:.6f}" r="9" fill="none" '
                         'stroke="#8e24aa" stroke-width="3"/>')
            sx, sy = move["start_position"]
            parts.append(f'<circle cx="{sx:.6f}" cy="{sy:.6f}" r="7" fill="none" '
                         'stroke="#00897b" stroke-width="3"/>')
            parts.append(f'<text x="{min_x+16}" y="{min_y+52}" font-family="sans-serif" '
                         f'font-size="15">move {move["moving_vertex"]} around {move["pivot_vertex"]}; '
                         f'{move["signed_rotation_degrees"]:.6f} degrees</text>')
        parts.append('</g>')
    parts.append('</svg>')
    return "\n".join(parts) + "\n"


def write_json(path, document):
    path.write_text(json.dumps(document, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def run(graph_path: Path, coordinate_path: Path, mds_report_path: Path, output_dir: Path):
    graph, start, initial, spacing = load_start(graph_path, coordinate_path, mds_report_path)
    final, final_diagnostics, moves, states, history = run_search(graph, start)
    report = {
        "schema": "graph-relax.endpoint-rotation-planarization.v1",
        "authoritative_start": coordinate_path.as_posix(),
        "search": {"coarse_step_degrees": COARSE_DEGREES,
                   "refinement_tolerance_degrees": REFINE_TOLERANCE_DEGREES,
                   "refinement_iteration_limit": REFINE_ITERATIONS},
        "permanent_spacing_start": spacing,
        "initial_crossing_count": initial["proper_unrelated_edge_crossing_count"],
        "initial_crossings": initial["proper_unrelated_edge_crossings"],
        "accepted_move_count": len(moves),
        "moves": moves,
        "crossing_history": history,
        "final_crossing_count": final_diagnostics["proper_unrelated_edge_crossing_count"],
        "final_crossings": final_diagnostics["proper_unrelated_edge_crossings"],
        "final_geometry_diagnostics": final_diagnostics,
        "stopping_reason": ("ZERO_CROSSINGS" if not final_diagnostics["proper_unrelated_edge_crossing_count"]
                            else "NO_ADMISSIBLE_IMPROVING_SINGLE_ENDPOINT_ROTATION"),
        "final_coordinates": {node: list(final[node]) for node in sorted(final)},
    }
    output_dir.mkdir(parents=True, exist_ok=True)
    viewbox = viewbox_for(states)
    (output_dir / "iamp_nearness_rotation_start.svg").write_text(
        static_svg(graph, start, "Nearness MDS — endpoint rotation start", viewbox), encoding="utf-8")
    (output_dir / "iamp_nearness_rotation_end.svg").write_text(
        static_svg(graph, final, "Nearness MDS — endpoint rotation end", viewbox), encoding="utf-8")
    (output_dir / "iamp_nearness_rotation_motion.svg").write_text(
        motion_svg(graph, states, moves, viewbox), encoding="utf-8")
    write_json(output_dir / "iamp_nearness_rotation_report.json", report)
    return report


def main():
    parser = argparse.ArgumentParser()
    base = Path("output/graph_first")
    parser.add_argument("--graph", type=Path, default=base / "iamp_graph.json")
    parser.add_argument("--coordinates", type=Path, default=base / "iamp_nearness_mds.json")
    parser.add_argument("--mds-report", type=Path, default=base / "iamp_nearness_mds_report.json")
    parser.add_argument("--output", type=Path, default=base)
    args = parser.parse_args()
    report = run(args.graph, args.coordinates, args.mds_report, args.output)
    print(json.dumps(report, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
