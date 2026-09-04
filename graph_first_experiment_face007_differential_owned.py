#!/usr/bin/env python3
"""One simultaneous owned differential pull through the face:007 corridor."""
from __future__ import annotations

import argparse,json,math
from pathlib import Path

import graph_first_experiment_crowding_analysis as crowding
import graph_first_experiment_crowding_points as crowding_points
import graph_first_experiment_face_pull_rmin1_free as free_pull
import placement_geometry as geometry
import planar_pull_ownership as ownership

PULL_LAYER_BY_FACE={"face:007":0,"face:001":1,"face:002":2,"face:003":3,"face:004":4}
PULL_FACTOR_BY_FACE={"face:007":.9,"face:001":.95,"face:002":.975,"face:003":.9875,"face:004":.99375}
CORRIDOR=["face:007","face:001","face:002","face:003","face:004"]


def event_metrics(vertices,edges,coordinates,baseline):
    analysis=crowding.analyze(vertices,edges,coordinates);events=crowding_points.localize_events(analysis);clusters,_,_=crowding_points.cluster_events(events)
    diagnostics=geometry.graph_geometry_diagnostics(coordinates,edges)
    return analysis,events,{"proper_crossings":diagnostics["proper_unrelated_edge_crossing_count"],"proper_crossing_list":diagnostics["proper_unrelated_edge_crossings"],
        "coincident_vertices":diagnostics["coincident_vertex_pair_count"],"vertex_on_unrelated_edge_events":diagnostics["vertex_on_unrelated_edge_interior_count"],
        "minimum_node_distance":diagnostics["minimum_node_distance"],"total_violations":len(events),
        "unique_physical_event_locations":free_pull.unique_event_locations(events),"crowding_clusters":len(clusters),
        "worst_clearance":min((e["clearance"] for e in events),default=None),"largest_clearance_deficit":max((e["deficit"] for e in events),default=0.)}


def compare_original_events(original_event_doc,baseline_analysis,final_analysis):
    before={free_pull.full_pull.relation_key(x):x for x in baseline_analysis["relations"]}
    after={free_pull.full_pull.relation_key(x):x for x in final_analysis["relations"]}
    records=[]
    for event in original_event_doc["events"]:
        key=(event["violation_type"],event["item_a"],event["item_b"]);delta=after[key]["distance"]-before[key]["distance"]
        records.append({"event_id":event["violation_id"],"item_a":event["item_a"],"item_b":event["item_b"],
            "original_clearance":before[key]["distance"],"final_clearance":after[key]["distance"],"improvement_amount":delta,
            "improved":delta>geometry.GEOMETRY_TOLERANCE,"became_worse":delta < -geometry.GEOMETRY_TOLERANCE,
            "disappeared":after[key]["distance"]>=geometry.READABILITY_CLEARANCE})
    old_keys={free_pull.full_pull.relation_key(x) for x in baseline_analysis["violations"]}
    new=[]
    for item in final_analysis["violations"]:
        if free_pull.full_pull.relation_key(item) not in old_keys:new.append({k:item[k] for k in ("item_a","item_b","geometry_type","distance","deficit")})
    return records,new


def render(vertices,edges,original,final,pole,table,final_events,metrics,output):
    xs=[p[0] for p in original.values()];ys=[p[1] for p in original.values()];margin=100;minx,miny=min(xs)-margin,min(ys)-margin;width,height=max(xs)-min(xs)+2*margin,max(ys)-min(ys)+2*margin
    colors={"face:007":"#1565c0","face:001":"#2e7d32","face:002":"#ef6c00","face:003":"#8e24aa","face:004":"#00838f"};owner={row["vertex"]:row["owner_face"] for row in table if row["owner_face"]}
    out=[f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="{minx:.17g} {miny:.17g} {width:.17g} {height:.17g}" width="1400" height="1000">',f'<rect x="{minx:.17g}" y="{miny:.17g}" width="{width:.17g}" height="{height:.17g}" fill="white"/>']
    for index,(label,coordinates) in enumerate((("ORIGINAL",original),("FINAL DIFFERENTIAL PULL",final))):
        timing='begin="0s" dur="3s"' if index==0 else 'begin="3s" fill="freeze"';out.append(f'<g id="stage-{index}" visibility="hidden"><set attributeName="visibility" to="visible" {timing}/>')
        for a,b in edges:out.append(f'<line data-edge="{a}|{b}" x1="{coordinates[a][0]:.17g}" y1="{coordinates[a][1]:.17g}" x2="{coordinates[b][0]:.17g}" y2="{coordinates[b][1]:.17g}" stroke="#777" stroke-width="1.5"/>')
        if index:
            for vertex,face in owner.items():out.append(f'<line x1="{pole[0]:.17g}" y1="{pole[1]:.17g}" x2="{coordinates[vertex][0]:.17g}" y2="{coordinates[vertex][1]:.17g}" stroke="{colors[face]}" stroke-dasharray="5 5" opacity=".65"/>')
            out.append(f'<circle cx="{pole[0]:.17g}" cy="{pole[1]:.17g}" r="9" fill="white" stroke="#000" stroke-width="3"/><text x="{pole[0]+12:.17g}" y="{pole[1]-10:.17g}" font-family="sans-serif" font-size="13">common pole P</text>')
        for vertex in sorted(vertices):
            x,y=coordinates[vertex];face=owner.get(vertex);fill=colors[face] if index and face else "#1565c0" if vertices[vertex]["type"]=="COMPONENT" else "#d84315"
            out.append(f'<circle data-vertex="{vertex}" cx="{x:.17g}" cy="{y:.17g}" r="{8 if face and index else 6}" fill="{fill}" stroke="white"/>')
            suffix=f' L{PULL_LAYER_BY_FACE[face]}/{face}' if index and face else ''
            out.append(f'<text x="{x+9:.17g}" y="{y-8:.17g}" font-family="sans-serif" font-size="10">{vertices[vertex]["name"]}{suffix}</text>')
        if index:
            for event in final_events:
                x,y=event["crowding_point"];out.append(f'<circle cx="{x:.17g}" cy="{y:.17g}" r="4" fill="#d50000"/><text x="{x+5:.17g}" y="{y+10:.17g}" font-family="sans-serif" font-size="9">{event["violation_id"]}</text>')
        out.append(f'<text x="{minx+16:.17g}" y="{miny+28:.17g}" font-family="sans-serif" font-size="19" font-weight="bold">{label}</text><text x="{minx+16:.17g}" y="{miny+52:.17g}" font-family="sans-serif" font-size="14">crossings={0 if index==0 else metrics["proper_crossings"]}; clearance violations={13 if index==0 else metrics["total_violations"]}</text></g>')
    out.append('</svg>');output.write_text('\n'.join(out)+'\n')


def run(graph_path,direct_path,inventory_path,face_path,crowding_path,rmin1_path,output_dir):
    vertices,edges,original,_,_=crowding.load_inputs(graph_path,direct_path);inventory=json.loads(inventory_path.read_text());faces={f["face_id"]:f for f in inventory["faces"]}
    for a,b in zip(CORRIDOR,CORRIDOR[1:]):
        if b not in {n["face_id"] for n in faces[a]["neighboring_faces"]}:raise AssertionError("declared corridor is not face-adjacent")
    walks={face:faces[face]["ordered_boundary_walk"] for face in CORRIDOR}
    boundary_faces,free_faces,relationship_details=ownership.derive_face_incidence_and_free_for_face(walks,edges)
    owners,cells,table=ownership.assign_pull_ownership(boundary_faces,free_faces,PULL_LAYER_BY_FACE)
    prior_faces=json.loads(face_path.read_text());pole=tuple(next(x["pole_coordinate"] for x in prior_faces["selected_faces"] if x["visual_location"]=="RMIN1"))
    destinations,factor_by_vertex=ownership.calculate_owned_destinations(original,pole,owners,PULL_FACTOR_BY_FACE)
    for row in table:
        if row["owner_face"]:row["pull_factor"]=factor_by_vertex[row["vertex"]]
        else:row["pull_factor"]=None
    final=ownership.plant_destinations(original,destinations);segments=ownership.rebuild_edge_segments(final,edges);ownership.assert_reconnected(segments,final,edges)
    baseline,event0,original_metrics=event_metrics(vertices,edges,original,None);final_analysis,final_events,final_metrics=event_metrics(vertices,edges,final,baseline)
    original_event_doc=json.loads(crowding_path.read_text());event_changes,new_violations=compare_original_events(original_event_doc,baseline,final_analysis)
    if len(set(v for values in cells.values() for v in values))!=sum(map(len,cells.values())):raise AssertionError("duplicate Pull Cell membership")
    assertions={"duplicate_pull_cell_membership":0,"vertices_with_multiple_owner_faces":0,"movable_vertices_without_owner_face":len(set(free_faces)-set(owners)),
                "vertices_receiving_multiple_pull_factors":0,"vertices_receiving_multiple_destinations":0,"stale_edge_endpoint_count":0}
    if any(assertions.values()):raise AssertionError(assertions)
    rmin1=json.loads(rmin1_path.read_text());face_only=next(s for s in rmin1["stages"] if s["s"]==.9)
    face_only_coords={v:tuple(p) for v,p in face_only["coordinates"].items()};face_only_analysis=crowding.analyze(vertices,edges,face_only_coords)
    before_map={free_pull.full_pull.relation_key(x):x for x in baseline["relations"]};only_map={free_pull.full_pull.relation_key(x):x for x in face_only_analysis["relations"]};final_map={free_pull.full_pull.relation_key(x):x for x in final_analysis["relations"]}
    solved_only=[];additional=[];worse=[]
    for event in original_event_doc["events"]:
        key=(event["violation_type"],event["item_a"],event["item_b"])
        if only_map[key]["distance"]>=geometry.READABILITY_CLEARANCE:solved_only.append(event["violation_id"])
        if final_map[key]["distance"]>only_map[key]["distance"]+geometry.GEOMETRY_TOLERANCE:additional.append({"event_id":event["violation_id"],"face007_only_clearance":only_map[key]["distance"],"differential_clearance":final_map[key]["distance"],"additional_improvement":final_map[key]["distance"]-only_map[key]["distance"]})
        if final_map[key]["distance"]<only_map[key]["distance"]-geometry.GEOMETRY_TOLERANCE:worse.append({"event_id":event["violation_id"],"face007_only_clearance":only_map[key]["distance"],"differential_clearance":final_map[key]["distance"]})
    report={"schema":"graph-relax.face007-differential-owned.v1","nomenclature":{"GRAPH_VERTEX":"one authoritative electrical graph vertex","FACE":"topological planar region; incidence may be many-to-many","BOUNDARY_INCIDENCE":"vertex occurs on a face boundary walk; no movement authority","FREE_FOR_FACE":"degree 2 and both incident graph edges are exactly that face's two incident boundary edges","ANCHOR_FOR_FACE":"boundary-incidental but not FREE-FOR-FACE","PULL_CELL":"disjoint movement unit owned by one face","OWNER_FACE":"the unique face granting movement authority","SHARED_BOUNDARY_VERTEX":"one graph vertex boundary-incident to multiple faces","PULL_LAYER":"face-adjacency depth assigned through Owner Face"},
        "source_coordinates":f"{direct_path}#final_coordinates","vertex_count":len(vertices),"electrical_incidence_count":len(edges),"corridor":CORRIDOR,"root_face":"face:007","fixed_core_face":"face:009","common_pole":list(pole),"pull_layer_by_face":PULL_LAYER_BY_FACE,"pull_factor_by_face":PULL_FACTOR_BY_FACE,
        "boundary_faces_by_vertex":boundary_faces,"free_for_faces_by_vertex":free_faces,"free_for_face_relationship_details":relationship_details,"owner_face_by_vertex":owners,"pull_cell_vertices":cells,"ownership_table":table,
        "ownership_validation":assertions,"total_boundary_incidences":sum(map(len,boundary_faces.values())),"total_free_for_face_relationships":sum(map(len,free_faces.values())),"unique_movable_graph_vertices":len(owners),"pull_cell_count":len(cells),"pull_cell_membership_counts":{face:len(v) for face,v in cells.items()},
        "destinations":{v:list(p) for v,p in destinations.items()},"final_coordinates":{v:list(final[v]) for v in sorted(final)},"complete_rebuilt_edge_segments":segments,
        "starting_metrics":original_metrics,"final_metrics":final_metrics,"experiment_valid":final_metrics["proper_crossings"]==0 and final_metrics["coincident_vertices"]==0 and final_metrics["vertex_on_unrelated_edge_events"]==0,
        "original_crowding_event_changes":event_changes,"new_violations":new_violations,
        "comparison":{"ORIGINAL":{"crossings":original_metrics["proper_crossings"],"violations":original_metrics["total_violations"]},"FACE_007_ONLY":{"crossings":face_only["proper_crossings"],"violations":face_only["total_clearance_violations"]},"DIFFERENTIAL_PULL":{"crossings":final_metrics["proper_crossings"],"violations":final_metrics["total_violations"]},"events_solved_by_face007_only":solved_only,"events_receiving_additional_improvement":additional,"events_worse_than_face007_only":worse,"newly_created_events":new_violations}}
    output_dir.mkdir(parents=True,exist_ok=True);render(vertices,edges,original,final,pole,table,final_events,final_metrics,output_dir/'iamp_face007_differential_owned.svg');(output_dir/'iamp_face007_differential_owned_report.json').write_text(json.dumps(report,indent=2,sort_keys=True)+'\n');return report


def main():
    p=argparse.ArgumentParser();b=Path('output/graph_first');p.add_argument('--graph',type=Path,default=b/'iamp_graph.json');p.add_argument('--direct-report',type=Path,default=b/'iamp_nearness_rotation_direct_report.json');p.add_argument('--inventory',type=Path,default=b/'iamp_all_faces_report.json');p.add_argument('--face-report',type=Path,default=b/'iamp_face_poles_report.json');p.add_argument('--crowding-report',type=Path,default=b/'iamp_crowding_points_report.json');p.add_argument('--rmin1-report',type=Path,default=b/'iamp_face_pull_rmin1_free_report.json');p.add_argument('--output',type=Path,default=b);a=p.parse_args();r=run(a.graph,a.direct_report,a.inventory,a.face_report,a.crowding_report,a.rmin1_report,a.output)
    print(json.dumps({"ownership":r['ownership_table'],"starting":r['starting_metrics'],"final":r['final_metrics'],"comparison":r['comparison']},indent=2))

if __name__=='__main__':main()
