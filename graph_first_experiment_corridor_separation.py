#!/usr/bin/env python3
"""One-shot, topology-derived separation of the central IAMP corridor."""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import graph_first_experiment_crowding_analysis as crowding
import graph_first_experiment_crowding_points as crowding_points
import graph_first_experiment_planar_region_inventory as inventory
import placement_geometry as geometry


DELTA = 0.5
SOURCE_STAGE = "AFTER face:007"
RMIN1_FACE = "face:007"
J_OUT_FACE = "face:008"


def graph_data(path):
    doc = json.loads(path.read_text())
    vertices = {item["id"]: item for item in doc["vertices"]}
    edges = [(item["source"], item["target"]) for item in doc["edges"]]
    return vertices, edges


def split_cycle(cycle, first, second):
    """Return the two directed paths between two vertices in a simple cycle."""
    i, j = cycle.index(first), cycle.index(second)
    if i > j:
        i, j = j, i
    path1 = cycle[i:j + 1]
    path2 = cycle[j:] + cycle[:i + 1]
    return path1, path2


def edge_key(a, b):
    return tuple(sorted((a, b)))


def path_edges(path):
    return [edge_key(a, b) for a, b in zip(path, path[1:])]


def point_polyline_distance(point, path, coordinates):
    return min(geometry.point_segment_distance(point, coordinates[a], coordinates[b])
               for a, b in zip(path, path[1:]))


def boundary_metrics(path_a, path_b, coordinates):
    distances_a = [point_polyline_distance(coordinates[v], path_b, coordinates)
                   for v in path_a]
    distances_b = [point_polyline_distance(coordinates[v], path_a, coordinates)
                   for v in path_b]
    minimum_nonincident = min(
        geometry.segment_segment_distance(coordinates[a], coordinates[b],
                                          coordinates[c], coordinates[d])
        for a, b in zip(path_a, path_a[1:])
        for c, d in zip(path_b, path_b[1:])
        if not ({a, b} & {c, d})
    )
    # The shared endpoints are excluded from the representative width: their
    # zero separation describes the corridor termini, not the crowded band.
    interior = distances_a[1:-1] + distances_b[1:-1]
    return {"minimum_polyline_separation_including_shared_termini": 0.0,
            "minimum_nonincident_segment_separation": minimum_nonincident,
            "mean_boundary_vertex_to_opposite_polyline_separation": sum(interior) / len(interior),
            "samples": {"boundary_a": distances_a, "boundary_b": distances_b}}


def metrics(vertices, edges, coordinates):
    diagnostic = geometry.graph_geometry_diagnostics(coordinates, edges)
    analysis = crowding.analyze(vertices, edges, coordinates)
    events = crowding_points.localize_events(analysis)
    clusters, _, _ = crowding_points.cluster_events(events)
    unique = []
    for event in events:
        p = tuple(event["crowding_point"])
        if not any(math.dist(p, q) <= geometry.GEOMETRY_TOLERANCE for q in unique):
            unique.append(p)
    return diagnostic, analysis, events, clusters, {
        "proper_crossings": diagnostic["proper_unrelated_edge_crossing_count"],
        "coincident_vertices": diagnostic["coincident_vertex_pair_count"],
        "vertex_on_unrelated_edge_events": diagnostic["vertex_on_unrelated_edge_interior_count"],
        "minimum_node_distance": diagnostic["minimum_node_distance"],
        "clearance_violations": len(events),
        "unique_crowding_event_locations": len(unique),
        "crowding_clusters": len(clusters),
        "worst_clearance": min((e["clearance"] for e in events), default=None),
    }


def relation_key(item):
    return item.get("geometry_type", item.get("violation_type")), item["item_a"], item["item_b"]


def event_changes(before_analysis, after_analysis, before_events):
    before = {relation_key(x): x for x in before_analysis["relations"]}
    after = {relation_key(x): x for x in after_analysis["relations"]}
    event_id = {relation_key(x): x["violation_id"] for x in before_events}
    rows = []
    for key, identifier in sorted(event_id.items(), key=lambda x: x[1]):
        old, new = before[key]["distance"], after[key]["distance"]
        change = new - old
        if new >= geometry.READABILITY_CLEARANCE:
            status = "disappeared"
        elif change > geometry.GEOMETRY_TOLERANCE:
            status = "improved"
        elif change < -geometry.GEOMETRY_TOLERANCE:
            status = "worsened"
        else:
            status = "unchanged"
        rows.append({"event_id": identifier, "item_a": before[key]["item_a"],
                     "item_b": before[key]["item_b"], "clearance_before": old,
                     "clearance_after": new, "change": change, "classification": status})
    old_bad = set(event_id)
    new = []
    for item in after_analysis["violations"]:
        key = relation_key(item)
        if key not in old_bad:
            new.append({"item_a": item["item_a"], "item_b": item["item_b"],
                        "geometry_type": item["geometry_type"], "clearance": item["distance"],
                        "deficit": item["deficit"]})
    return rows, sorted(new, key=lambda x: (x["clearance"], x["item_a"], x["item_b"]))


def render(vertices, edges, before, after, path_a, path_b, centerline, normal,
           moved_a, moved_b, before_events, after_events, result, output):
    allx = [p[0] for p in before.values()]; ally = [p[1] for p in before.values()]
    minx, maxx, miny, maxy = min(allx), max(allx), min(ally), max(ally)
    margin = 90.; panel_w = maxx-minx+2*margin; height=maxy-miny+2*margin
    total_w=2*panel_w; view_minx=minx-margin; view_miny=miny-margin
    parts=[f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="{view_minx:.17g} {view_miny:.17g} {total_w:.17g} {height:.17g}" width="1800" height="900">',
           f'<rect x="{view_minx:.17g}" y="{view_miny:.17g}" width="{total_w:.17g}" height="{height:.17g}" fill="white"/>']
    for panel,(title,co,events) in enumerate((("BEFORE",before,before_events),(f"AFTER — {result}",after,after_events))):
        shift=panel*panel_w
        parts.append(f'<g id="{title.split()[0].lower()}" transform="translate({shift:.17g} 0)">')
        for a,b in edges:
            parts.append(f'<line data-edge="{a}|{b}" x1="{co[a][0]:.17g}" y1="{co[a][1]:.17g}" x2="{co[b][0]:.17g}" y2="{co[b][1]:.17g}" stroke="#777" stroke-width="1.5"/>')
        for path,color in ((path_a,"#00897b"),(path_b,"#8e24aa")):
            parts.append(f'<polyline points="'+" ".join(f'{co[v][0]:.17g},{co[v][1]:.17g}' for v in path)+f'" fill="none" stroke="{color}" stroke-width="5" opacity=".72"/>')
        a,b=centerline
        parts.append(f'<line x1="{a[0]:.17g}" y1="{a[1]:.17g}" x2="{b[0]:.17g}" y2="{b[1]:.17g}" stroke="#1565c0" stroke-width="2" stroke-dasharray="8 5"/>')
        mid=((a[0]+b[0])/2,(a[1]+b[1])/2); scale=42
        for sign,label,color in ((1,"+n RMIN1 side","#00897b"),(-1,"-n J_OUT side","#8e24aa")):
            end=(mid[0]+sign*scale*normal[0],mid[1]+sign*scale*normal[1])
            parts.append(f'<line x1="{mid[0]:.17g}" y1="{mid[1]:.17g}" x2="{end[0]:.17g}" y2="{end[1]:.17g}" stroke="{color}" stroke-width="3"/><text x="{end[0]+5:.17g}" y="{end[1]-5:.17g}" font-size="11" fill="{color}">{label}</text>')
        if panel:
            for v,color in [(v,"#00897b") for v in moved_a]+[(v,"#8e24aa") for v in moved_b]:
                p,q=before[v],after[v]
                parts.append(f'<line x1="{p[0]:.17g}" y1="{p[1]:.17g}" x2="{q[0]:.17g}" y2="{q[1]:.17g}" stroke="{color}" stroke-width="3"/>')
        for v in sorted(vertices):
            x,y=co[v]; color="#00897b" if v in moved_a else "#8e24aa" if v in moved_b else "#1565c0" if vertices[v]["type"]=="COMPONENT" else "#d84315"
            parts.append(f'<circle data-vertex="{v}" cx="{x:.17g}" cy="{y:.17g}" r="{8 if v in moved_a+moved_b else 6}" fill="{color}" stroke="white"/><text x="{x+8:.17g}" y="{y-8:.17g}" font-family="sans-serif" font-size="10">{vertices[v]["name"]}</text>')
        for event in events:
            x,y=event["crowding_point"]
            parts.append(f'<circle cx="{x:.17g}" cy="{y:.17g}" r="4" fill="#d50000"/><text x="{x+5:.17g}" y="{y+13:.17g}" font-size="9" fill="#b71c1c">{event["violation_id"]}</text>')
        parts.append(f'<text x="{minx-margin+15:.17g}" y="{miny-margin+26:.17g}" font-family="sans-serif" font-size="18" font-weight="bold">{title}</text></g>')
    parts.append('</svg>');output.write_text('\n'.join(parts)+'\n')


def run(graph_path, sequential_path, output_dir):
    vertices, edges = graph_data(graph_path)
    history = json.loads(sequential_path.read_text())
    source = next(stage for stage in history["stages"] if stage["stage"] == SOURCE_STAGE)
    original = {v: tuple(p) for v,p in source["coordinates"].items()}
    before_diag,before_analysis,before_events,before_clusters,before_metrics=metrics(vertices,edges,original)
    if (len(vertices),len(edges),before_metrics["proper_crossings"],before_metrics["clearance_violations"]) != (34,44,0,6):
        raise AssertionError(f"historical baseline discrepancy: {before_metrics}")

    faces,dart_face=inventory.enumerate_faces(vertices,edges,original)
    inventory.add_classification_and_adjacency(faces,edges,dart_face)
    by_id={f["face_id"]:f for f in faces}
    face_a,face_b=by_id[RMIN1_FACE],by_id[J_OUT_FACE]
    shared=sorted(set(face_a["distinct_boundary_vertices"]) & set(face_b["distinct_boundary_vertices"]))
    if shared != ["component:U1","net:GND"]:
        raise AssertionError(f"corridor side faces do not share expected termini: {shared}")
    paths_a=split_cycle(face_a["ordered_boundary_walk"],*shared)
    path_a=next(p for p in paths_a if "component:RMIN1" not in p)
    paths_b=split_cycle(face_b["ordered_boundary_walk"],*shared)
    path_b=next(p for p in paths_b if "component:J_OUT" in p)
    # Orient both boundaries U1 -> GND.
    for path in (path_a,path_b):
        if path[0] != "component:U1": path.reverse()

    start,end=original["component:U1"],original["net:GND"]
    length=math.dist(start,end); tangent=((end[0]-start[0])/length,(end[1]-start[1])/length)
    base_normal=(-tangent[1],tangent[0])
    signed=lambda p:(p[0]-start[0])*base_normal[0]+(p[1]-start[1])*base_normal[1]
    mean_a=sum(signed(original[v]) for v in path_a[1:-1])/len(path_a[1:-1])
    mean_b=sum(signed(original[v]) for v in path_b[1:-1])/len(path_b[1:-1])
    if mean_a*mean_b >= 0: raise AssertionError("measured boundaries are not on opposite centerline sides")
    normal=base_normal if mean_a>0 else (-base_normal[0],-base_normal[1])

    # Corridor-side boundary vertices may move only if the complete graph has
    # exactly the two incident boundary edges. This prevents dragging termini,
    # external attachments, or pendant structures.
    incidence={v:[] for v in vertices}
    for edge in edges:
        incidence[edge[0]].append(edge_key(*edge));incidence[edge[1]].append(edge_key(*edge))
    moved_a=[];moved_b=[];classification={}
    boundary_sets={"RMIN1":set(path_edges(path_a)),"J_OUT":set(path_edges(path_b))}
    corridor_vertices=sorted(set(path_a)|set(path_b))
    for v in sorted(vertices):
        side="RMIN1" if v in path_a else "J_OUT" if v in path_b else None
        h=(original[v][0]-start[0])*normal[0]+(original[v][1]-start[1])*normal[1]
        if side and abs(h)<=geometry.GEOMETRY_TOLERANCE:
            reason="CENTER / FIXED: corridor terminus on centerline"
        elif side and len(incidence[v])==2 and set(incidence[v]).issubset(boundary_sets[side]):
            reason=f"MOVABLE: topology-free internal {side}-side boundary vertex"
            (moved_a if side=="RMIN1" else moved_b).append(v)
        elif side:
            reason="FIXED: external anchor or non-boundary incidence"
        else:
            reason="FIXED: not on either measured corridor boundary"
        classification[v]={"signed_distance_before":h,"corridor_side":side,"movement":reason,
                           "degree":len(incidence[v]),"incident_edges":[list(e) for e in sorted(incidence[v])]}
    if set(moved_a)&set(moved_b): raise AssertionError("duplicate movement membership")

    final=dict(original)
    for v in moved_a: final[v]=(original[v][0]+DELTA*normal[0],original[v][1]+DELTA*normal[1])
    for v in moved_b: final[v]=(original[v][0]-DELTA*normal[0],original[v][1]-DELTA*normal[1])
    movement=[]
    for v in moved_a+moved_b:
        h0=classification[v]["signed_distance_before"]
        h1=(final[v][0]-start[0])*normal[0]+(final[v][1]-start[1])*normal[1]
        if abs(h1)<=abs(h0)+geometry.GEOMETRY_TOLERANCE: raise AssertionError(f"{v} moved toward corridor")
        classification[v]["signed_distance_after"]=h1
        movement.append({"vertex":v,"side":classification[v]["corridor_side"],"before":list(original[v]),"after":list(final[v]),"displacement":[final[v][0]-original[v][0],final[v][1]-original[v][1]],"absolute_signed_distance_increased":True})

    segments=geometry.rebuild_straight_edge_segments(final,edges)
    reconnect=geometry.validate_rebuilt_edge_segments(final,edges,segments)
    if not reconnect["valid"]:raise AssertionError(reconnect)
    after_diag,after_analysis,after_events,after_clusters,after_metrics=metrics(vertices,edges,final)
    event_rows,new_events=event_changes(before_analysis,after_analysis,before_events)
    boundary_before=boundary_metrics(path_a,path_b,original);boundary_after=boundary_metrics(path_a,path_b,final)
    opened=(boundary_after["mean_boundary_vertex_to_opposite_polyline_separation"] > boundary_before["mean_boundary_vertex_to_opposite_polyline_separation"] + geometry.GEOMETRY_TOLERANCE)
    hard_valid=(after_metrics["proper_crossings"]==0 and after_metrics["coincident_vertices"]==0 and after_metrics["vertex_on_unrelated_edge_events"]==0 and after_metrics["minimum_node_distance"]+geometry.GEOMETRY_TOLERANCE>=geometry.MIN_NODE_DISTANCE)
    if not hard_valid: interpretation="C. GEOMETRIC FAILURE"
    elif not opened: interpretation="D. NO EFFECT"
    elif after_metrics["clearance_violations"]<=before_metrics["clearance_violations"]:interpretation="A. SUCCESS"
    else:interpretation="B. PARTIAL SUCCESS"

    band_events=[e["violation_id"] for e in before_events if e["violation_id"]!="V01"]
    report={"schema":"graph-relax.corridor-separation.v1","source_coordinates":f"{sequential_path}#{SOURCE_STAGE}",
      "diagnostic_delta":DELTA,"readability_clearance":geometry.READABILITY_CLEARANCE,"vertex_count":len(vertices),"electrical_incidence_count":len(edges),
      "corridor":{"events_included":band_events,"excluded_events":["V01"],"rmin1_side_face":RMIN1_FACE,"j_out_side_face":J_OUT_FACE,
        "faces_touching_rmin1_boundary":sorted({dart_face.get((a,b)) for a,b in zip(path_a,path_a[1:])}|{dart_face.get((b,a)) for a,b in zip(path_a,path_a[1:])}),
        "faces_touching_j_out_boundary":sorted({dart_face.get((a,b)) for a,b in zip(path_b,path_b[1:])}|{dart_face.get((b,a)) for a,b in zip(path_b,path_b[1:])}),
        "rmin1_side_boundary":{"ordered_vertices":path_a,"ordered_edges":[list(e) for e in path_edges(path_a)]},
        "j_out_side_boundary":{"ordered_vertices":path_b,"ordered_edges":[list(e) for e in path_edges(path_b)]}},
      "centerline":{"definition":"straight chord joining the two topology-shared corridor termini component:U1 and net:GND","endpoints":[list(start),list(end)],"unit_tangent":list(tangent),"unit_normal_toward_rmin1_side":list(normal),"mean_signed_distance_rmin1_boundary_internal":abs(mean_a),"mean_signed_distance_j_out_boundary_internal":-abs(mean_b)},
      "vertex_classification":classification,"rmin1_side_movable_vertices":moved_a,"j_out_side_movable_vertices":moved_b,"movement":movement,
      "before_metrics":before_metrics,"after_metrics":after_metrics,"before_crowding_events":before_events,"after_crowding_events":after_events,
      "before_crowding_clusters":before_clusters,"after_crowding_clusters":after_clusters,"original_event_outcomes":event_rows,"new_violations":new_events,
      "boundary_separation_before":boundary_before,"boundary_separation_after":boundary_after,"corridor_mean_separation_increased":opened,
      "edge_reconnection":reconnect,"hard_geometry_valid":hard_valid,"interpretation":interpretation,
      "coordinates_before":{v:list(original[v]) for v in sorted(original)},"coordinates_after":{v:list(final[v]) for v in sorted(final)}}
    output_dir.mkdir(parents=True,exist_ok=True)
    render(vertices,edges,original,final,path_a,path_b,(start,end),normal,moved_a,moved_b,before_events,after_events,interpretation,output_dir/'iamp_corridor_separation.svg')
    (output_dir/'iamp_corridor_separation_report.json').write_text(json.dumps(report,indent=2,sort_keys=True)+'\n')
    return report


def main():
    p=argparse.ArgumentParser();b=Path('output/graph_first')
    p.add_argument('--graph',type=Path,default=b/'iamp_graph.json');p.add_argument('--sequential-report',type=Path,default=b/'iamp_sequential_local_face_pulls_report.json');p.add_argument('--output',type=Path,default=b)
    a=p.parse_args();r=run(a.graph,a.sequential_report,a.output)
    print(r['source_coordinates']);print('moved',r['rmin1_side_movable_vertices'],r['j_out_side_movable_vertices'])
    print('before',r['before_metrics']);print('after',r['after_metrics']);print(r['interpretation'])

if __name__=='__main__':main()
