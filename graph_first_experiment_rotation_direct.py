#!/usr/bin/env python3
"""Direct geometric single-endpoint rotation without angular search."""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import placement_geometry as geometry
import graph_first_experiment_endpoint_rotation as sweep
import direct_crossing_repair as generic_repair
import direct_clearance_geometry as clearance_geometry


TAU = 2.0 * math.pi
ANGLE_TOLERANCE = 1e-12
NUMERICAL_SAFETY_FACTOR = 16.0


def svg_number(value):
    """Serialize a coordinate with enough precision to round-trip its float."""
    return format(float(value), ".17g")


def normalize(angle):
    return angle % TAU


def ccw_span(start, end):
    return normalize(end - start)


def point_on_ccw_arc(angle, start, end, tolerance=ANGLE_TOLERANCE):
    return ccw_span(start, angle) <= ccw_span(start, end) + tolerance


def segment_disk_clip(first, second, center, radius):
    """Return the segment portion inside a closed disk, analytically."""
    dx, dy = second[0] - first[0], second[1] - first[1]
    fx, fy = first[0] - center[0], first[1] - center[1]
    a = dx*dx + dy*dy
    tolerance = geometry.GEOMETRY_TOLERANCE
    if a <= tolerance*tolerance:
        return ([first, first] if math.dist(first, center) <= radius + tolerance else [])
    b = 2.0 * (fx*dx + fy*dy)
    c = fx*fx + fy*fy - radius*radius
    discriminant = b*b - 4*a*c
    inside_first = c <= tolerance
    ex, ey = second[0] - center[0], second[1] - center[1]
    inside_second = ex*ex + ey*ey - radius*radius <= tolerance
    parameters = [0.0] if inside_first else []
    if discriminant >= -tolerance:
        root = math.sqrt(max(0.0, discriminant))
        parameters.extend(value for value in ((-b-root)/(2*a), (-b+root)/(2*a))
                          if -tolerance <= value <= 1.0+tolerance)
    if inside_second:
        parameters.append(1.0)
    parameters = sorted(max(0.0, min(1.0, value)) for value in parameters)
    unique = []
    for value in parameters:
        if not unique or abs(value-unique[-1]) > tolerance:
            unique.append(value)
    if not unique:
        return []
    low, high = unique[0], unique[-1]
    return [(first[0]+low*dx, first[1]+low*dy),
            (first[0]+high*dx, first[1]+high*dy)]


def segments_intersect_inclusive(a, b, c, d):
    tolerance = geometry.GEOMETRY_TOLERANCE
    values = (geometry.cross(a, b, c), geometry.cross(a, b, d),
              geometry.cross(c, d, a), geometry.cross(c, d, b))
    if ((values[0] > tolerance and values[1] < -tolerance) or
        (values[0] < -tolerance and values[1] > tolerance)) and \
       ((values[2] > tolerance and values[3] < -tolerance) or
        (values[2] < -tolerance and values[3] > tolerance)):
        return True
    def on(point, first, second):
        if abs(geometry.cross(first, second, point)) > tolerance:
            return False
        return (min(first[0], second[0])-tolerance <= point[0] <= max(first[0], second[0])+tolerance
                and min(first[1], second[1])-tolerance <= point[1] <= max(first[1], second[1])+tolerance)
    return on(c, a, b) or on(d, a, b) or on(a, c, d) or on(b, c, d)


def target_forbidden_arc(pivot, radius, clipped, theta_old):
    if len(clipped) != 2 or math.dist(clipped[0], clipped[1]) <= geometry.GEOMETRY_TOLERANCE:
        raise ValueError("target clip is not a nonzero segment")
    angles = [math.atan2(point[1]-pivot[1], point[0]-pivot[0]) for point in clipped]
    first, second = angles
    midpoint_first = normalize(first + ccw_span(first, second)/2.0)
    midpoint_second = normalize(second + ccw_span(second, first)/2.0)
    ray_first = (pivot[0]+radius*math.cos(midpoint_first),
                 pivot[1]+radius*math.sin(midpoint_first))
    ray_second = (pivot[0]+radius*math.cos(midpoint_second),
                  pivot[1]+radius*math.sin(midpoint_second))
    first_hits = segments_intersect_inclusive(pivot, ray_first, clipped[0], clipped[1])
    second_hits = segments_intersect_inclusive(pivot, ray_second, clipped[0], clipped[1])
    if first_hits == second_hits:
        raise ValueError("could not uniquely determine target forbidden angular interval")
    start, end = (first, second) if first_hits else (second, first)
    if not point_on_ccw_arc(theta_old, start, end):
        raise ValueError("current moving-edge angle is outside calculated target interval")
    return normalize(start), normalize(end)


def direct_solution(graph, positions, moving, pivot, crossed_edge, target_key, counters):
    positions_evaluated = 0
    pivot_point, moving_point = positions[pivot], positions[moving]
    radius = math.dist(pivot_point, moving_point)
    theta_old = math.atan2(moving_point[1]-pivot_point[1], moving_point[0]-pivot_point[0])
    clipped = segment_disk_clip(positions[crossed_edge[0]], positions[crossed_edge[1]],
                                pivot_point, radius)
    if len(clipped) != 2 or math.dist(clipped[0], clipped[1]) <= geometry.GEOMETRY_TOLERANCE:
        return {"valid": False, "reason": "EMPTY_OR_POINT_TARGET_CLIP", "moving": moving,
                "pivot": pivot, "clipped": clipped}
    try:
        mathematical_target_arc = target_forbidden_arc(pivot_point, radius, clipped, theta_old)
    except ValueError as error:
        return {"valid": False, "reason": str(error), "moving": moving,
                "pivot": pivot, "clipped": clipped}
    mathematical_deltas = [clearance_geometry.signed_delta(angle, theta_old)
                           for angle in mathematical_target_arc]
    mathematical_boundary_delta = min(mathematical_deltas, key=abs)
    try:
        clearance_answer = clearance_geometry.nearest_clearance_escape(
            pivot_point, radius, moving_point,
            positions[crossed_edge[0]], positions[crossed_edge[1]],
            geometry.REQUIRED_CENTERLINE_CLEARANCE)
    except ValueError as error:
        return {"valid": False, "reason": str(error), "moving": moving,
                "pivot": pivot, "clipped": clipped,
                "mathematical_target_interval": mathematical_target_arc}
    target_arc = clearance_answer["forbidden_interval"]
    boundary_delta = clearance_geometry.signed_delta(
        clearance_answer["boundary_angle"], theta_old)
    direction = -1.0 if boundary_delta <= 0 else 1.0
    drawing_scale = max(math.dist(pivot_point, point) for point in positions.values())
    epsilon = NUMERICAL_SAFETY_FACTOR * geometry.GEOMETRY_TOLERANCE * \
        max(1.0, drawing_scale) / max(radius, geometry.GEOMETRY_TOLERANCE)
    signed_delta = boundary_delta + direction*epsilon
    theta_new = theta_old + signed_delta
    candidate = dict(positions)
    candidate[moving] = (pivot_point[0]+radius*math.cos(theta_new),
                         pivot_point[1]+radius*math.sin(theta_new))
    # Plant first. The coordinate map is authoritative; every graph segment is
    # then rebuilt from current endpoints before any geometry is validated.
    rebuilt_segments = geometry.rebuild_straight_edge_segments(candidate, graph.edges)
    reconnection = geometry.validate_rebuilt_edge_segments(
        candidate, graph.edges, rebuilt_segments)
    if not reconnection["valid"]:
        raise AssertionError(f"plant-and-reconnect invariant failed: {reconnection['failures']}")
    positions_evaluated += 1
    if positions_evaluated > 1:
        raise AssertionError("more than one position evaluated for endpoint/pivot choice")
    counters["positions_per_choice"].append({
        "crossing": target_key,
        "moving": moving,
        "pivot": pivot,
        "positions_evaluated": positions_evaluated,
    })
    counters["candidate_evaluations"] += 1
    diagnostics = geometry.graph_geometry_diagnostics(candidate, graph.edges)
    target_centerline_distance = geometry.segment_segment_distance(
        pivot_point, candidate[moving], positions[crossed_edge[0]], positions[crossed_edge[1]])
    admissible = (diagnostics["coincident_vertex_pair_count"] == 0
                  and diagnostics["vertex_on_unrelated_edge_interior_count"] == 0
                  and target_centerline_distance + geometry.GEOMETRY_TOLERANCE
                      >= geometry.REQUIRED_CENTERLINE_CLEARANCE)
    return {
        "valid": admissible, "reason": None if admissible else "DIRECT_POSITION_FAILED_VALIDATION",
        "moving": moving, "pivot": pivot, "radius": radius,
        "clipped": [list(point) for point in clipped], "target_interval": list(target_arc),
        "mathematical_target_interval": list(mathematical_target_arc),
        "mathematical_crossing_boundary_angle": theta_old + mathematical_boundary_delta,
        "visible_clearance_boundary_angle": clearance_answer["boundary_angle"],
        "visible_clearance_boundary_feature": clearance_answer["boundary_feature"],
        "additional_boundary_rotation_for_visible_clearance_radians": (
            abs(clearance_geometry.signed_delta(clearance_answer["boundary_angle"], theta_old))
            - abs(mathematical_boundary_delta)),
        "required_centerline_clearance": geometry.REQUIRED_CENTERLINE_CLEARANCE,
        "target_centerline_distance": target_centerline_distance,
        "resulting_minimum_node_distance": diagnostics["minimum_node_distance"],
        "numerical_legal_side_offset_radians": direction*epsilon,
        "theta_old": theta_old, "theta_new": normalize(theta_new),
        "signed_delta_radians": signed_delta, "signed_delta_degrees": math.degrees(signed_delta),
        "displacement": math.dist(moving_point, candidate[moving]),
        "candidate": candidate, "diagnostics": diagnostics,
        "rebuilt_segments": rebuilt_segments, "reconnection": reconnection,
    }


def best_direct_move(graph, positions, current, counters):
    before = current["proper_unrelated_edge_crossing_count"]
    def solver_factory(target):
        def solve(moving, pivot, crossed_edge):
            solution = direct_solution(
                graph, positions, moving, pivot, crossed_edge, target, counters)
            solution["positions_evaluated"] = 1 if "candidate" in solution else 0
            solution["resulting_total_crossings"] = (
                solution["diagnostics"]["proper_unrelated_edge_crossing_count"]
                if solution.get("valid") else None)
            return solution
        return solve

    selected, audits = generic_repair.best_move_across_crossings(
        current["proper_unrelated_edge_crossings"], before, solver_factory)
    for audit in audits:
        counters["solutions_per_crossing"].append(
            {"crossing": audit["crossing"], "calculated": len(audit["solutions"])})
        serializable_solutions = []
        for solution in audit["solutions"]:
            outcome = {key: value for key, value in solution.items()
                       if key not in {"candidate", "diagnostics", "rebuilt_segments"}}
            counters["all_solution_outcomes"].append(outcome)
            serializable_solutions.append(outcome)
        counters["orientation_audits"].append(
            {"crossing": audit["crossing"], "solutions": serializable_solutions})
    return selected


def direct_search(graph, start):
    positions = dict(start)
    current = geometry.graph_geometry_diagnostics(positions, graph.edges)
    counters = {"candidate_evaluations": 0, "positions_per_choice": [],
                "solutions_per_crossing": [], "all_solution_outcomes": [],
                "orientation_audits": []}
    moves, states, history = [], [dict(positions)], [current["proper_unrelated_edge_crossings"]]
    while current["proper_unrelated_edge_crossing_count"]:
        solution = best_direct_move(graph, positions, current, counters)
        if solution is None:
            break
        complete = geometry.graph_geometry_diagnostics(solution["candidate"], graph.edges)
        if complete != solution["diagnostics"]:
            raise AssertionError("direct solution validation mismatch")
        if (complete["proper_unrelated_edge_crossing_count"] >=
                current["proper_unrelated_edge_crossing_count"]):
            raise AssertionError("refusing to commit a non-reducing direct move")
        rebuilt = geometry.rebuild_straight_edge_segments(solution["candidate"], graph.edges)
        reconnection = geometry.validate_rebuilt_edge_segments(
            solution["candidate"], graph.edges, rebuilt)
        if not reconnection["valid"] or rebuilt != solution["rebuilt_segments"]:
            raise AssertionError("accepted move was not validated with its reconnected graph")
        neighbors = sorted(graph.neighbors(solution["moving"]))
        incident = [item for item in rebuilt if solution["moving"] in item["edge"]]
        before_keys = {sweep.crossing_key(item) for item in current["proper_unrelated_edge_crossings"]}
        after_keys = {sweep.crossing_key(item) for item in complete["proper_unrelated_edge_crossings"]}
        moves.append({
            "iteration": len(moves)+1,
            "crossing_being_repaired": {"edge_a": list(solution["target"][0]),
                                        "edge_b": list(solution["target"][1])},
            "moving_vertex": solution["moving"], "pivot_vertex": solution["pivot"],
            "radius": solution["radius"], "clipped_crossed_segment": solution["clipped"],
            "target_forbidden_angular_interval": solution["target_interval"],
            "mathematical_crossing_boundary_angle": solution["mathematical_crossing_boundary_angle"],
            "visible_clearance_boundary_angle": solution["visible_clearance_boundary_angle"],
            "visible_clearance_boundary_feature": solution["visible_clearance_boundary_feature"],
            "additional_boundary_rotation_for_visible_clearance_radians": solution[
                "additional_boundary_rotation_for_visible_clearance_radians"],
            "required_centerline_clearance": solution["required_centerline_clearance"],
            "resulting_target_centerline_distance": solution["target_centerline_distance"],
            "theta_old": solution["theta_old"], "theta_new": solution["theta_new"],
            "signed_rotation_radians": solution["signed_delta_radians"],
            "signed_rotation_degrees": solution["signed_delta_degrees"],
            "absolute_rotation_degrees": abs(solution["signed_delta_degrees"]),
            "moving_vertex_displacement": solution["displacement"],
            "crossing_count_before": current["proper_unrelated_edge_crossing_count"],
            "crossing_count_after": complete["proper_unrelated_edge_crossing_count"],
            "crossings_removed": [{"edge_a": list(key[0]), "edge_b": list(key[1])}
                                  for key in sorted(before_keys-after_keys)],
            "crossings_newly_created": [{"edge_a": list(key[0]), "edge_b": list(key[1])}
                                        for key in sorted(after_keys-before_keys)],
            "minimum_node_distance_after": complete["minimum_node_distance"],
            "start_position": list(positions[solution["moving"]]),
            "end_position": list(solution["candidate"][solution["moving"]]),
            "pivot_position": list(positions[solution["pivot"]]),
            "preserved_edge_length_before": solution["radius"],
            "moved_vertex_degree": graph.degree[solution["moving"]],
            "moved_vertex_neighbors": neighbors,
            "incident_edges_regenerated": incident,
            "all_incident_edges_regenerated": len(incident) == graph.degree[solution["moving"]],
            "reconnection_validation": reconnection,
        })
        positions, current = solution["candidate"], complete
        states.append(dict(positions)); history.append(current["proper_unrelated_edge_crossings"])
    maximum_positions = max((item["positions_evaluated"]
                             for item in counters["positions_per_choice"]), default=0)
    maximum_solutions = max((item["calculated"] for item in counters["solutions_per_crossing"]), default=0)
    if maximum_positions > 1 or maximum_solutions > 4:
        raise AssertionError("direct-solution cardinality invariant failed")
    return positions, current, moves, states, history, counters


def rendered_edge_records(graph, positions):
    """Return the authoritative active edge geometry used by direct SVGs."""
    records = geometry.rebuild_straight_edge_segments(positions, graph.edges)
    validation = geometry.validate_rebuilt_edge_segments(positions, graph.edges, records)
    if not validation["valid"]:
        raise AssertionError(f"refusing to render stale graph edges: {validation['failures']}")
    return records


def direct_static_svg(graph, positions, title, viewbox):
    """Render direct results without discarding topology-significant digits."""
    min_x, min_y, width, height = viewbox
    records = rendered_edge_records(graph, positions)
    parts = [f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="{svg_number(min_x)} '
             f'{svg_number(min_y)} {svg_number(width)} {svg_number(height)}" '
             'width="1400" height="1000">',
             f'<rect x="{svg_number(min_x)}" y="{svg_number(min_y)}" '
             f'width="{svg_number(width)}" height="{svg_number(height)}" fill="white"/>',
             f'<text x="{svg_number(min_x+16)}" y="{svg_number(min_y+28)}" '
             f'font-family="sans-serif" font-size="18">{title}</text>']
    for record in records:
        first, second = record["edge"]
        start, end = record["segment_start"], record["segment_end"]
        parts.append(f'<line data-edge-u="{first}" data-edge-v="{second}" '
                     f'x1="{svg_number(start[0])}" y1="{svg_number(start[1])}" '
                     f'x2="{svg_number(end[0])}" y2="{svg_number(end[1])}" '
                     'stroke="#777" stroke-width="1.5"/>')
    for node in sorted(graph):
        x, y = positions[node]
        color = "#1565c0" if graph.nodes[node]["type"] == "COMPONENT" else "#d84315"
        parts.append(f'<circle data-active-vertex="{node}" cx="{svg_number(x)}" '
                     f'cy="{svg_number(y)}" r="6" fill="{color}" stroke="white" stroke-width="1"/>')
        parts.append(f'<text x="{svg_number(x+8)}" y="{svg_number(y-8)}" '
                     f'font-family="sans-serif" font-size="11">{graph.nodes[node]["name"]}</text>')
    parts.append('</svg>')
    return "\n".join(parts) + "\n"


def direct_motion_svg(graph, states, moves, viewbox):
    """Render complete cumulative post-move states with reconnected stars."""
    if len(states) != len(moves) + 1:
        raise AssertionError("motion states do not correspond one-to-one with accepted moves")
    min_x, min_y, width, height = viewbox
    stage_duration = 3
    final_begin = len(states) * stage_duration
    parts = [f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="{min_x} {min_y} {width} {height}" '
             'width="1400" height="1000">',
             '<defs><marker id="move-arrow" markerWidth="8" markerHeight="8" refX="7" refY="4" '
             'orient="auto"><path d="M0,0 L8,4 L0,8 z" fill="#00897b"/></marker></defs>']
    final_active_records = None
    previous_serialized_crossings = None
    for index, state in enumerate(states):
        move = moves[index-1] if index else None
        active_records = rendered_edge_records(graph, state)
        serialized_state = {node: (float(svg_number(point[0])), float(svg_number(point[1])))
                            for node, point in state.items()}
        serialized_diagnostics = geometry.graph_geometry_diagnostics(serialized_state, graph.edges)
        if index and (serialized_diagnostics["proper_unrelated_edge_crossing_count"] >=
                      previous_serialized_crossings):
            raise AssertionError("serialized SVG stage does not strictly reduce total crossings")
        previous_serialized_crossings = serialized_diagnostics["proper_unrelated_edge_crossing_count"]
        if move:
            if tuple(state[move["moving_vertex"]]) != tuple(move["end_position"]):
                raise AssertionError("post-move stage does not contain planted final coordinate")
            incident = [record for record in active_records
                        if move["moving_vertex"] in record["edge"]]
            if len(incident) != move["moved_vertex_degree"]:
                raise AssertionError("rendered reconnected star has wrong degree")
            for record in incident:
                first, second = record["edge"]
                endpoint = (record["segment_start"] if first == move["moving_vertex"]
                            else record["segment_end"])
                if math.dist(tuple(endpoint), tuple(move["end_position"])) > geometry.GEOMETRY_TOLERANCE:
                    raise AssertionError("rendered incident edge terminates at stale node coordinate")
        parts.append('<g visibility="hidden">')
        parts.append(f'<set attributeName="visibility" to="visible" begin="{index*stage_duration}s" '
                     f'dur="{stage_duration}s"/>')
        parts.append(f'<rect x="{min_x}" y="{min_y}" width="{width}" height="{height}" fill="white"/>')
        label = "Stage 0 — original graph" if not move else \
            f'Stage {index} — POST-MOVE graph: {move["moving_vertex"]} planted and reconnected'
        parts.append(f'<text x="{min_x+16}" y="{min_y+28}" font-family="sans-serif" font-size="18">{label}</text>')
        for record in active_records:
            first, second = record["edge"]
            highlighted = bool(move and move["moving_vertex"] in {first, second})
            color, stroke_width = ("#e53935", 4.5) if highlighted else ("#777", 1.5)
            start, end = record["segment_start"], record["segment_end"]
            parts.append(f'<line data-edge-u="{first}" data-edge-v="{second}" '
                         f'x1="{svg_number(start[0])}" y1="{svg_number(start[1])}" '
                         f'x2="{svg_number(end[0])}" y2="{svg_number(end[1])}" '
                         f'stroke="{color}" stroke-width="{stroke_width}"/>')
        for node in sorted(graph):
            x, y = state[node]
            highlighted = bool(move and node == move["moving_vertex"])
            color = "#1565c0" if graph.nodes[node]["type"] == "COMPONENT" else "#d84315"
            radius = 11 if highlighted else 6
            stroke = "#ffeb3b" if highlighted else "white"
            stroke_width = 4 if highlighted else 1
            parts.append(f'<circle data-active-vertex="{node}" cx="{svg_number(x)}" cy="{svg_number(y)}" '
                         f'r="{radius}" fill="{color}" stroke="{stroke}" stroke-width="{stroke_width}"/>')
            parts.append(f'<text x="{svg_number(x+8)}" y="{svg_number(y-8)}" font-family="sans-serif" '
                         f'font-size="11">{graph.nodes[node]["name"]}</text>')
        if move:
            old_x, old_y = move["start_position"]
            new_x, new_y = move["end_position"]
            pivot_x, pivot_y = move["pivot_position"]
            radius = move["preserved_edge_length_before"]
            parts.append(f'<circle cx="{svg_number(pivot_x)}" cy="{svg_number(pivot_y)}" r="{svg_number(radius)}" fill="none" '
                         'stroke="#8e24aa" stroke-width="1.5" stroke-dasharray="7 5" opacity="0.35"/>')
            parts.append(f'<circle cx="{svg_number(old_x)}" cy="{svg_number(old_y)}" r="7" fill="none" '
                         'stroke="#00897b" stroke-width="2" stroke-dasharray="3 3" opacity="0.55"/>')
            parts.append(f'<line x1="{svg_number(old_x)}" y1="{svg_number(old_y)}" x2="{svg_number(new_x)}" y2="{svg_number(new_y)}" '
                         'stroke="#00897b" stroke-width="2" opacity="0.65" marker-end="url(#move-arrow)"/>')
            parts.append(f'<text x="{min_x+16}" y="{min_y+52}" font-family="sans-serif" font-size="15">'
                         f'Highlighted active star: degree {move["moved_vertex_degree"]}; old position is annotation only</text>')
        parts.append('</g>')
        if index == len(states)-1:
            final_active_records = active_records
    expected_final = rendered_edge_records(graph, states[-1])
    if final_active_records != expected_final:
        raise AssertionError("motion final active graph differs from direct end graph")
    # After the last move demonstration, replace its annotations with one clean
    # final graph and freeze that graph on screen permanently.
    parts.append('<g id="final-persistent-graph" visibility="hidden">')
    parts.append(f'<set attributeName="visibility" to="visible" begin="{final_begin}s" '
                 'dur="indefinite" fill="freeze"/>')
    parts.append(f'<rect x="{min_x}" y="{min_y}" width="{width}" height="{height}" fill="white"/>')
    parts.append(f'<text x="{min_x+16}" y="{min_y+28}" font-family="sans-serif" font-size="18">'
                 'FINAL — direct geometric crossing repair</text>')
    for record in expected_final:
        first, second = record["edge"]
        start, end = record["segment_start"], record["segment_end"]
        parts.append(f'<line data-final-edge-u="{first}" data-final-edge-v="{second}" '
                     f'x1="{svg_number(start[0])}" y1="{svg_number(start[1])}" '
                     f'x2="{svg_number(end[0])}" y2="{svg_number(end[1])}" '
                     'stroke="#777" stroke-width="1.5"/>')
    for node in sorted(graph):
        x, y = states[-1][node]
        color = "#1565c0" if graph.nodes[node]["type"] == "COMPONENT" else "#d84315"
        parts.append(f'<circle data-final-vertex="{node}" cx="{svg_number(x)}" cy="{svg_number(y)}" '
                     f'r="6" fill="{color}" stroke="white" stroke-width="1"/>')
        parts.append(f'<text x="{svg_number(x+8)}" y="{svg_number(y-8)}" font-family="sans-serif" '
                     f'font-size="11">{graph.nodes[node]["name"]}</text>')
    parts.append('</g>')
    parts.append('</svg>')
    return "\n".join(parts) + "\n"


def regenerate_motion_from_saved_report(graph_path, coordinate_path, mds_report_path,
                                        direct_report_path, motion_path):
    """Regenerate animation only; never recalculate accepted direct moves."""
    graph, start, _, _ = sweep.load_start(graph_path, coordinate_path, mds_report_path)
    report = json.loads(direct_report_path.read_text(encoding="utf-8"))
    current = dict(start)
    states = [dict(current)]
    for move in report["moves"]:
        current[move["moving_vertex"]] = tuple(move["end_position"])
        states.append(dict(current))
    expected = {node: tuple(point) for node, point in report["final_coordinates"].items()}
    if current != expected:
        raise AssertionError("saved cumulative moves do not reconstruct saved final coordinates")
    motion_path.write_text(
        direct_motion_svg(graph, states, report["moves"], sweep.viewbox_for(states)),
        encoding="utf-8")


def run(graph_path, coordinate_path, mds_report_path, output_dir):
    graph, start, initial, spacing = sweep.load_start(graph_path, coordinate_path, mds_report_path)
    final, diagnostics, moves, states, history, counters = direct_search(graph, start)
    report = {
        "schema": "graph-relax.direct-geometric-endpoint-rotation.v1",
        "initial_crossing_count": initial["proper_unrelated_edge_crossing_count"],
        "initial_crossings": initial["proper_unrelated_edge_crossings"],
        "final_crossing_count": diagnostics["proper_unrelated_edge_crossing_count"],
        "final_crossings": diagnostics["proper_unrelated_edge_crossings"],
        "accepted_move_count": len(moves), "moves": moves, "crossing_history": history,
        "final_geometry_diagnostics": diagnostics,
        "stopping_reason": ("ZERO_CROSSINGS" if not diagnostics["proper_unrelated_edge_crossing_count"]
                            else "NO_REDUCING_DIRECT_GEOMETRIC_SOLUTION"),
        "direct_solution_assertions": {
            "maximum_positions_evaluated_per_endpoint_pivot": max(
                (item["positions_evaluated"] for item in counters["positions_per_choice"]),
                default=0),
            "maximum_direct_solutions_calculated_per_crossing": max(
                (item["calculated"] for item in counters["solutions_per_crossing"]), default=0),
            "candidate_evaluations": counters["candidate_evaluations"],
            "passed": True,
        },
        "all_direct_solution_outcomes": counters["all_solution_outcomes"],
        "generic_orientation_audits": counters["orientation_audits"],
        "generic_mechanism": "direct_crossing_repair.best_move_across_crossings",
        "final_coordinates": {node: list(final[node]) for node in sorted(final)},
    }
    output_dir.mkdir(parents=True, exist_ok=True)
    viewbox = sweep.viewbox_for(states)
    (output_dir/"iamp_nearness_rotation_direct_start.svg").write_text(
        direct_static_svg(graph, start, "Direct geometric endpoint rotation start", viewbox), encoding="utf-8")
    (output_dir/"iamp_nearness_rotation_direct_end.svg").write_text(
        direct_static_svg(graph, final, "Direct geometric endpoint rotation end", viewbox), encoding="utf-8")
    (output_dir/"iamp_nearness_rotation_direct_motion.svg").write_text(
        direct_motion_svg(graph, states, moves, viewbox), encoding="utf-8")
    sweep.write_json(output_dir/"iamp_nearness_rotation_direct_report.json", report)
    return report


def main():
    parser = argparse.ArgumentParser()
    base = Path("output/graph_first")
    parser.add_argument("--graph", type=Path, default=base/"iamp_graph.json")
    parser.add_argument("--coordinates", type=Path, default=base/"iamp_nearness_mds.json")
    parser.add_argument("--mds-report", type=Path, default=base/"iamp_nearness_mds_report.json")
    parser.add_argument("--output", type=Path, default=base)
    args = parser.parse_args()
    print(json.dumps(run(args.graph, args.coordinates, args.mds_report, args.output),
                     indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
