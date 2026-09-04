#!/usr/bin/env python3
"""Sequential local-face pulls with persistent exclusive Pull-Cell ownership."""
from __future__ import annotations

import argparse,json,math
from pathlib import Path

import graph_first_experiment_crowding_analysis as crowding
import graph_first_experiment_crowding_points as crowding_points
import graph_first_experiment_face_poles as face_poles
import graph_first_experiment_planar_region_inventory as inventory
import placement_geometry as geometry
import planar_pull_ownership as ownership

FACE_ORDER=("face:007","face:001","face:002","face:003","face:004")
PULL_FACTOR=.90


def key(item):return item["geometry_type"],item["item_a"],item["item_b"]


def geometry_metrics(vertices,edges,coordinates):
    diagnostics=geometry.graph_geometry_diagnostics(coordinates,edges);analysis=crowding.analyze(vertices,edges,coordinates)
    events=crowding_points.localize_events(analysis);clusters,_,_=crowding_points.cluster_events(events)
    points=[]
    for event in events:
        point=tuple(event["crowding_point"])
        if not any(math.dist(point,p)<=geometry.GEOMETRY_TOLERANCE for p in points):points.append(point)
    return diagnostics,analysis,events,{"proper_crossings":diagnostics["proper_unrelated_edge_crossing_count"],
        "proper_crossing_list":diagnostics["proper_unrelated_edge_crossings"],"coincident_vertices":diagnostics["coincident_vertex_pair_count"],
        "vertex_on_unrelated_edge_events":diagnostics["vertex_on_unrelated_edge_interior_count"],"minimum_node_distance":diagnostics["minimum_node_distance"],
        "clearance_violations":len(events),"unique_crowding_event_locations":len(points),"crowding_clusters":len(clusters),
        "worst_clearance":min((e["clearance"] for e in events),default=None),"largest_clearance_deficit":max((e["deficit"] for e in events),default=0.)}


def current_face(face_id,vertices,edges,coordinates,reference_walk):
    faces,_=inventory.enumerate_faces(vertices,edges,coordinates)
    matches=[face for face in faces if inventory.undirected_walk_key(face["ordered_boundary_walk"])==inventory.undirected_walk_key(reference_walk)]
    if len(matches)!=1:raise AssertionError(f"could not preserve identity of {face_id}")
    face=matches[0]
    if face["face_id"]!=face_id:raise AssertionError(f"stable face ID changed: {face_id} -> {face['face_id']}")
    if not face["simple_boundary"]:raise AssertionError(f"selected face {face_id} is not a simple polygon")
    return face


def ownership_snapshot(table,owner,movement_stage):
    result=[]
    for row in table:
        vertex=row["vertex"];face=owner.get(vertex)
        result.append({"vertex":vertex,"boundary_incidental_faces":row["all_boundary_incidental_selected_faces"],
            "free_for_face_relationships":row["all_free_for_face_qualifying_faces"],"owner_face":face,
            "pull_cell":f"cell:{face.split(':')[1]}" if face else None,"moved":vertex in movement_stage,
            "movement_stage":movement_stage.get(vertex)})
    return result


def relation_changes(before,after):
    old={key(x):x for x in before["relations"]};new={key(x):x for x in after["relations"]}
    old_bad={k for k,x in old.items() if x["distance"]<geometry.READABILITY_CLEARANCE};new_bad={k for k,x in new.items() if x["distance"]<geometry.READABILITY_CLEARANCE}
    def record(k):return {"item_a":new[k]["item_a"],"item_b":new[k]["item_b"],"geometry_type":new[k]["geometry_type"],"clearance_before":old[k]["distance"],"clearance_after":new[k]["distance"],"change":new[k]["distance"]-old[k]["distance"]}
    return {"removed":[record(k) for k in sorted(old_bad-new_bad)],"new":[record(k) for k in sorted(new_bad-old_bad)],
        "improved":[record(k) for k in sorted(old_bad) if new[k]["distance"]>old[k]["distance"]+geometry.GEOMETRY_TOLERANCE],
        "worsened":[record(k) for k in sorted(old_bad) if new[k]["distance"]<old[k]["distance"]-geometry.GEOMETRY_TOLERANCE]}


def tracked_events(original_event_doc,baseline,current,new_ids):
    base={key(x):x for x in baseline["relations"]};now={key(x):x for x in current["relations"]};original_keys={}
    records=[]
    for event in original_event_doc["events"]:
        k=(event["violation_type"],event["item_a"],event["item_b"]);original_keys[k]=event["violation_id"];delta=now[k]["distance"]-base[k]["distance"]
        records.append({"event_id":event["violation_id"],"item_a":event["item_a"],"item_b":event["item_b"],"original_clearance":base[k]["distance"],"current_clearance":now[k]["distance"],
            "remains":now[k]["distance"]<geometry.READABILITY_CLEARANCE,"improved":delta>geometry.GEOMETRY_TOLERANCE,"disappeared":now[k]["distance"]>=geometry.READABILITY_CLEARANCE})
    current_bad={key(x):x for x in current["violations"]}
    for k in sorted(set(current_bad)-set(original_keys)):
        if k not in new_ids:new_ids[k]=f"N{len(new_ids)+1:02d}"
    new_records=[]
    for k,event_id in sorted(new_ids.items(),key=lambda pair:pair[1]):
        item=now[k];new_records.append({"event_id":event_id,"item_a":item["item_a"],"item_b":item["item_b"],"current_clearance":item["distance"],"active_violation":k in current_bad})
    return records,new_records


def render(vertices,edges,original,stages,faces,output):
    xs=[p[0] for p in original.values()];ys=[p[1] for p in original.values()];margin=100;minx,miny=min(xs)-margin,min(ys)-margin;width,height=max(xs)-min(xs)+2*margin,max(ys)-min(ys)+2*margin
    out=[f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="{minx:.17g} {miny:.17g} {width:.17g} {height:.17g}" width="1400" height="1000">',f'<rect x="{minx:.17g}" y="{miny:.17g}" width="{width:.17g}" height="{height:.17g}" fill="white"/>']
    for index,stage in enumerate(stages):
        begin=index*2.5;timing=f'begin="{begin}s" dur="2.5s"' if index<len(stages)-1 else f'begin="{begin}s" fill="freeze"';out.append(f'<g id="stage-{index}" visibility="hidden"><set attributeName="visibility" to="visible" {timing}/>')
        co={v:tuple(p) for v,p in stage["coordinates"].items()};face=faces.get(stage["face"]);walk=face["ordered_boundary_walk"] if face else []
        if face:out.append('<polygon points="'+' '.join(f'{co[v][0]:.17g},{co[v][1]:.17g}' for v in walk)+'" fill="#90caf9" fill-opacity=".2"/>')
        for a,b in edges:out.append(f'<line data-edge="{a}|{b}" x1="{co[a][0]:.17g}" y1="{co[a][1]:.17g}" x2="{co[b][0]:.17g}" y2="{co[b][1]:.17g}" stroke="#777" stroke-width="1.5"/>')
        if stage["pole_coordinate"]:
            px,py=stage["pole_coordinate"]
            for v in stage["new_vertices_moved"]:out.append(f'<line x1="{px:.17g}" y1="{py:.17g}" x2="{co[v][0]:.17g}" y2="{co[v][1]:.17g}" stroke="#2e7d32" stroke-dasharray="5 5"/>')
            out.append(f'<circle cx="{px:.17g}" cy="{py:.17g}" r="9" fill="white" stroke="#6a1b9a" stroke-width="4"/>')
        for v in sorted(vertices):
            x,y=co[v];fill="#2e7d32" if v in stage["new_vertices_moved"] else "#795548" if v in stage["owned_vertices"] else "#1565c0" if vertices[v]["type"]=="COMPONENT" else "#d84315"
            out.append(f'<circle data-vertex="{v}" cx="{x:.17g}" cy="{y:.17g}" r="{8 if v in stage["new_vertices_moved"] else 6}" fill="{fill}" stroke="white"/><text x="{x+9:.17g}" y="{y-8:.17g}" font-family="sans-serif" font-size="10">{vertices[v]["name"]}</text>')
        for event in stage["crowding_event_points"]:
            x,y=event["crowding_point"];out.append(f'<circle cx="{x:.17g}" cy="{y:.17g}" r="4" fill="#d50000"/>')
        out.append(f'<text x="{minx+16:.17g}" y="{miny+28:.17g}" font-family="sans-serif" font-size="19" font-weight="bold">{stage["stage"]} — {stage["result"]}</text><text x="{minx+16:.17g}" y="{miny+52:.17g}" font-family="sans-serif" font-size="14">crossings={stage["proper_crossings"]}; violations={stage["clearance_violations"]}; green=new; brown=previously owned</text></g>')
    out.append('</svg>');output.write_text('\n'.join(out)+'\n')


def run(graph_path,direct_path,inventory_path,crowding_path,face007_result_path,output_dir):
    vertices,edges,original,_,_=crowding.load_inputs(graph_path,direct_path);inventory_doc=json.loads(inventory_path.read_text());reference={f["face_id"]:f for f in inventory_doc["faces"]};walks={f:reference[f]["ordered_boundary_walk"] for f in FACE_ORDER}
    boundary_faces,free_faces,details=ownership.derive_face_incidence_and_free_for_face(walks,edges);_,_,base_table=ownership.assign_pull_ownership(boundary_faces,free_faces,{face:i for i,face in enumerate(FACE_ORDER)})
    original_events=json.loads(crowding_path.read_text());_,baseline,events0,metrics0=geometry_metrics(vertices,edges,original);current=dict(original);owners={};cells={face:[] for face in FACE_ORDER};movement_stage={};new_ids={}
    original_tracking,new_tracking=tracked_events(original_events,baseline,baseline,new_ids)
    stages=[{"stage":"ORIGINAL","face":None,"result":"BASELINE","accepted":True,"new_vertices_moved":[],"owned_vertices":[],"pole_coordinate":None,"coordinates":{v:list(current[v]) for v in sorted(current)},"crowding_event_points":events0,"ownership_table":ownership_snapshot(base_table,owners,movement_stage),"original_events":original_tracking,"new_events":new_tracking,**metrics0}]
    operations=[]
    established=json.loads(face007_result_path.read_text());expected007=next(s for s in established["stages"] if s["s"]==.9)
    for face_id in FACE_ORDER:
        before=dict(current);before_diag,before_analysis,before_events,before_metrics=geometry_metrics(vertices,edges,before)
        face=current_face(face_id,vertices,edges,current,walks[face_id]);points=[current[v] for v in face["ordered_boundary_walk"]];area=abs(face_poles.polygon_area(points));pole,clearance,_=face_poles.safe_pole(points)
        qualifying=sorted(free_faces.get(v,[]) and v for v in free_faces if face_id in free_faces[v]);qualifying=[v for v in qualifying if v]
        already=sorted(v for v in qualifying if v in owners);newly=sorted(v for v in qualifying if v not in owners)
        tentative_owner=dict(owners);tentative_cells={f:list(vs) for f,vs in cells.items()};dest={}
        for v in newly:
            tentative_owner[v]=face_id;tentative_cells[face_id].append(v);x,y=current[v];dest[v]=(pole[0]+PULL_FACTOR*(x-pole[0]),pole[1]+PULL_FACTOR*(y-pole[1]))
        candidate=ownership.plant_destinations(current,dest);segments=ownership.rebuild_edge_segments(candidate,edges);ownership.assert_reconnected(segments,candidate,edges)
        cand_diag,cand_analysis,cand_events,cand_metrics=geometry_metrics(vertices,edges,candidate)
        accepted=not newly or (cand_metrics["proper_crossings"]==0 and cand_metrics["coincident_vertices"]==0 and cand_metrics["vertex_on_unrelated_edge_events"]==0)
        if accepted:
            current=candidate;owners=tentative_owner;cells=tentative_cells
            for v in newly:movement_stage[v]=face_id
            result="SKIPPED — NO UNOWNED FREE VERTICES" if not newly else ("IMPROVED CROWDING" if cand_metrics["clearance_violations"]<before_metrics["clearance_violations"] else "WORSENED CROWDING" if cand_metrics["clearance_violations"]>before_metrics["clearance_violations"] else "NO CHANGE")
        else:result="REJECTED — INVALID GEOMETRY"
        after_diag,after_analysis,after_events,after_metrics=geometry_metrics(vertices,edges,current);changes=relation_changes(before_analysis,after_analysis);orig_track,new_track=tracked_events(original_events,baseline,after_analysis,new_ids)
        if face_id=="face:007":
            max_delta=max(math.dist(current[v],tuple(expected007["coordinates"][v])) for v in current)
            if max_delta>geometry.GEOMETRY_TOLERANCE or after_metrics["proper_crossings"]!=0 or after_metrics["clearance_violations"]!=6:raise AssertionError(f"face:007 reproduction discrepancy {max_delta} {after_metrics}")
        memberships=[v for values in cells.values() for v in values]
        assertions={"duplicate_pull_cell_membership":len(memberships)-len(set(memberships)),"multiple_owner_faces":0,"vertices_moved_more_than_once":len(movement_stage)-len(set(movement_stage)),"multiple_destination_coordinates":0,"stale_edge_endpoints":0}
        if any(assertions.values()):raise AssertionError(assertions)
        stage={"stage":f"AFTER {face_id}","face":face_id,"result":result,"accepted":accepted,"new_vertices_moved":newly if accepted else [],"owned_vertices":sorted(owners),"pole_coordinate":list(pole),"coordinates":{v:list(current[v]) for v in sorted(current)},"crowding_event_points":after_events,
            "current_face_area_before_pull":area,"free_for_face_count":len(qualifying),"already_owned_free_count":len(already),"newly_movable_free_count":len(newly),"current_safe_pole_clearance":clearance,"total_displacement":sum(math.dist(before[v],current[v]) for v in newly) if accepted else 0.,
            "violations_before":before_metrics["clearance_violations"],"violations_after":after_metrics["clearance_violations"],"new_violations_this_stage":changes["new"],"removed_violations_this_stage":changes["removed"],"improved_violations_this_stage":changes["improved"],"worsened_violations_this_stage":changes["worsened"],
            "candidate_metrics":cand_metrics,"ownership_assertions":assertions,"ownership_table":ownership_snapshot(base_table,owners,movement_stage),"original_events":orig_track,"new_events":new_track,**after_metrics}
        stages.append(stage);operations.append(stage)
    report={"schema":"graph-relax.sequential-local-face-pulls.v1","source_coordinates":f"{direct_path}#final_coordinates","face_order":list(FACE_ORDER),"pull_factor":PULL_FACTOR,"vertex_count":len(vertices),"electrical_incidence_count":len(edges),
        "boundary_faces_by_vertex":boundary_faces,"free_for_faces_by_vertex":free_faces,"free_for_face_relationship_details":details,"persistent_owner_face_by_vertex":owners,"persistent_pull_cell_vertices":cells,"movement_stage_by_vertex":movement_stage,"stages":stages,"operations":operations,
        "final_metrics":{k:stages[-1][k] for k in metrics0},"final_coordinates":stages[-1]["coordinates"]}
    output_dir.mkdir(parents=True,exist_ok=True);render(vertices,edges,original,stages,reference,output_dir/'iamp_sequential_local_face_pulls.svg');(output_dir/'iamp_sequential_local_face_pulls_report.json').write_text(json.dumps(report,indent=2,sort_keys=True)+'\n');return report


def main():
    p=argparse.ArgumentParser();b=Path('output/graph_first');p.add_argument('--graph',type=Path,default=b/'iamp_graph.json');p.add_argument('--direct-report',type=Path,default=b/'iamp_nearness_rotation_direct_report.json');p.add_argument('--inventory',type=Path,default=b/'iamp_all_faces_report.json');p.add_argument('--crowding-report',type=Path,default=b/'iamp_crowding_points_report.json');p.add_argument('--face007-result',type=Path,default=b/'iamp_face_pull_rmin1_free_report.json');p.add_argument('--output',type=Path,default=b);a=p.parse_args();r=run(a.graph,a.direct_report,a.inventory,a.crowding_report,a.face007_result,a.output)
    for s in r['stages']:print(s['stage'],s.get('face'),s.get('new_vertices_moved'),len(s.get('owned_vertices',[])),s['proper_crossings'],s.get('violations_before',s['clearance_violations']),s.get('violations_after',s['clearance_violations']),len(s.get('new_violations_this_stage',[])),len(s.get('removed_violations_this_stage',[])),s['worst_clearance'],s['result'])

if __name__=='__main__':main()
