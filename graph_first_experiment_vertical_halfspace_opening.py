#!/usr/bin/env python3
"""IAMP-only vertical opening of the two face regions around the central band."""
from __future__ import annotations

import argparse,json,math
from pathlib import Path

import graph_first_experiment_corridor_separation as corridor
import graph_first_experiment_crowding_analysis as crowding
import graph_first_experiment_crowding_points as crowding_points
import graph_first_experiment_planar_region_inventory as inventory
import placement_geometry as geometry

UPPER_FACE="face:007";LOWER_FACE="face:008"

def internal_mean_y(path,coordinates):
    values=[coordinates[v][1] for v in path[1:-1]]
    if not values:raise AssertionError("boundary path has no internal vertices")
    return sum(values)/len(values)

def split_selected(face,first,second,required_vertex):
    paths=corridor.split_cycle(face["ordered_boundary_walk"],first,second)
    selected=next(path for path in paths if required_vertex in path)
    if selected[0]!=first:selected.reverse()
    other=next(path for path in paths if path is not selected)
    if other[0]!=first:other.reverse()
    return selected,other

def band_width(upper,lower,coordinates):
    return corridor.boundary_metrics(upper,lower,coordinates)

def event_outcomes(before_analysis,after_analysis,before_events):
    return corridor.event_changes(before_analysis,after_analysis,before_events)

def render(vertices,edges,before,after,report,output):
    all_points=list(before.values())+list(after.values());xs=[p[0] for p in all_points];ys=[p[1] for p in all_points]
    margin=90.;minx,maxx,miny,maxy=min(xs),max(xs),min(ys),max(ys);pw=maxx-minx+2*margin;h=maxy-miny+2*margin
    parts=[f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="{minx-margin:.17g} {miny-margin:.17g} {2*pw:.17g} {h:.17g}" width="1800" height="900">',f'<rect x="{minx-margin:.17g}" y="{miny-margin:.17g}" width="{2*pw:.17g}" height="{h:.17g}" fill="white"/>']
    upper=set(report["upper_movable_vertices"]);lower=set(report["lower_movable_vertices"])
    paths=(report["boundaries"]["upper_corridor_facing"],report["boundaries"]["upper_opposite"],report["boundaries"]["lower_corridor_facing"],report["boundaries"]["lower_opposite"])
    for panel,(title,co,events) in enumerate((("BEFORE",before,report["crowding_events_before"]),("AFTER",after,report["crowding_events_after"]))):
        shift=panel*pw;parts.append(f'<g id="{title.lower()}" transform="translate({shift:.17g} 0)">')
        y=report["Y_CENTER"]
        parts.append(f'<line x1="{minx:.17g}" y1="{y:.17g}" x2="{maxx:.17g}" y2="{y:.17g}" stroke="#1565c0" stroke-width="2" stroke-dasharray="8 5"/><text x="{minx+5:.17g}" y="{y-7:.17g}" font-size="11" fill="#1565c0">Y_CENTER</text>')
        for a,b in edges:parts.append(f'<line data-edge="{a}|{b}" x1="{co[a][0]:.17g}" y1="{co[a][1]:.17g}" x2="{co[b][0]:.17g}" y2="{co[b][1]:.17g}" stroke="#777" stroke-width="1.5"/>')
        for path,color in zip(paths,("#00897b","#80cbc4","#8e24aa","#ce93d8")):
            parts.append('<polyline points="'+' '.join(f'{co[v][0]:.17g},{co[v][1]:.17g}' for v in path)+f'" fill="none" stroke="{color}" stroke-width="4" opacity=".75"/>')
        if panel:
            for v,color in [(v,"#00897b") for v in upper]+[(v,"#8e24aa") for v in lower]:
                p,q=before[v],after[v];parts.append(f'<line x1="{p[0]:.17g}" y1="{p[1]:.17g}" x2="{q[0]:.17g}" y2="{q[1]:.17g}" stroke="{color}" stroke-width="3"/>')
        for v in sorted(vertices):
            x,yv=co[v];color="#00897b" if v in upper else "#8e24aa" if v in lower else "#1565c0" if vertices[v]["type"]=="COMPONENT" else "#d84315"
            parts.append(f'<circle data-vertex="{v}" cx="{x:.17g}" cy="{yv:.17g}" r="{8 if v in upper|lower else 6}" fill="{color}" stroke="white"/><text x="{x+8:.17g}" y="{yv-8:.17g}" font-size="10">{vertices[v]["name"]}</text>')
        for e in events:
            x,yv=e["crowding_point"];parts.append(f'<circle cx="{x:.17g}" cy="{yv:.17g}" r="4" fill="#d50000"/>')
        parts.append(f'<text x="{minx-margin+15:.17g}" y="{miny-margin+25:.17g}" font-size="18" font-weight="bold">{title}</text><text x="{minx-margin+15:.17g}" y="{miny-margin+47:.17g}" font-size="13">crossings={report[title.lower()+"_metrics"]["proper_crossings"]}; violations={report[title.lower()+"_metrics"]["clearance_violations"]}</text></g>')
    parts.append('</svg>');output.write_text('\n'.join(parts)+'\n')

def run(graph_path,corridor_report_path,output_dir):
    base=json.loads(corridor_report_path.read_text());g=json.loads(graph_path.read_text());vertices={x["id"]:x for x in g["vertices"]};edges=[(x["source"],x["target"]) for x in g["edges"]];before={v:tuple(p) for v,p in base["coordinates_before"].items()}
    bd,ba,be,bc,bm=corridor.metrics(vertices,edges,before)
    if (len(vertices),len(edges),bm["proper_crossings"],bm["clearance_violations"])!=(34,44,0,6):raise AssertionError(f"baseline discrepancy {bm}")
    band_ids=set(base["corridor"]["events_included"]);band_events=[e for e in be if e["violation_id"] in band_ids]
    y_center=sum(e["crowding_point"][1] for e in band_events)/len(band_events)
    faces,darts=inventory.enumerate_faces(vertices,edges,before);inventory.add_classification_and_adjacency(faces,edges,darts);by={f["face_id"]:f for f in faces}
    upper_corridor,upper_outer=split_selected(by[UPPER_FACE],"component:U1","net:GND","net:EA")
    lower_outer,lower_corridor=split_selected(by[LOWER_FACE],"component:U1","net:GND","component:J_OUT")
    upper_corr_y=internal_mean_y(upper_corridor,before);upper_outer_y=internal_mean_y(upper_outer,before)
    lower_corr_y=internal_mean_y(lower_corridor,before);lower_outer_y=internal_mean_y(lower_outer,before)
    h_upper=upper_corr_y-upper_outer_y;h_lower=lower_outer_y-lower_corr_y
    if h_upper<=0 or h_lower<=0:raise AssertionError("face-region depth does not point outward")
    move_up=h_upper/2;move_down=h_lower/2
    upper_free=sorted(by[UPPER_FACE]["free_vertices"]);lower_free=sorted(by[LOWER_FACE]["free_vertices"])
    overlap=set(upper_free)&set(lower_free);upper_move=[v for v in upper_free if v not in overlap];lower_move=[v for v in lower_free if v not in overlap]
    final=dict(before)
    for v in upper_move:final[v]=(before[v][0],before[v][1]-move_up)
    for v in lower_move:final[v]=(before[v][0],before[v][1]+move_down)
    for v in final:
        if final[v][0]!=before[v][0]:raise AssertionError("X coordinate changed")
        if v in upper_move and not final[v][1]<before[v][1]:raise AssertionError("upper vertex did not move up")
        if v in lower_move and not final[v][1]>before[v][1]:raise AssertionError("lower vertex did not move down")
    segments=geometry.rebuild_straight_edge_segments(final,edges);reconnect=geometry.validate_rebuilt_edge_segments(final,edges,segments)
    if not reconnect["valid"]:raise AssertionError(reconnect)
    ad,aa,ae,ac,am=corridor.metrics(vertices,edges,final);outcomes,new=event_outcomes(ba,aa,be)
    width_before=band_width(upper_corridor,lower_corridor,before);width_after=band_width(upper_corridor,lower_corridor,final)
    hard_valid=am["proper_crossings"]==0 and am["coincident_vertices"]==0 and am["vertex_on_unrelated_edge_events"]==0 and am["minimum_node_distance"]+geometry.GEOMETRY_TOLERANCE>=geometry.MIN_NODE_DISTANCE
    classification={}
    for v in sorted(vertices):
        memberships=[]
        if v in by[UPPER_FACE]["distinct_boundary_vertices"]:memberships.append(UPPER_FACE)
        if v in by[LOWER_FACE]["distinct_boundary_vertices"]:memberships.append(LOWER_FACE)
        if v in overlap:reason="FIXED: shared by upper and lower face descriptions"
        elif v in upper_move:reason="MOVABLE UPPER: FREE-FOR-FACE face:007"
        elif v in lower_move:reason="MOVABLE LOWER: FREE-FOR-FACE face:008"
        elif memberships:reason="FIXED: face anchor or external/pendant attachment anchor"
        else:reason="FIXED: unrelated to participating face regions"
        classification[v]={"participating_faces":memberships,"reason":reason,"coordinate_before":list(before[v]),"coordinate_after":list(final[v])}
    if not hard_valid:interpretation="C. GEOMETRIC FAILURE"
    elif width_after["mean_boundary_vertex_to_opposite_polyline_separation"]<=width_before["mean_boundary_vertex_to_opposite_polyline_separation"]+geometry.GEOMETRY_TOLERANCE:interpretation="D. NO EFFECT"
    elif am["clearance_violations"]<=bm["clearance_violations"]:interpretation="A. SUCCESS"
    else:interpretation="B. PARTIAL SUCCESS"
    report={"schema":"graph-relax.iamp-vertical-halfspace-opening.v1","source_coordinates":base["source_coordinates"],"vertex_count":len(vertices),"electrical_incidence_count":len(edges),"readability_clearance":geometry.READABILITY_CLEARANCE,
      "Y_CENTER":y_center,"central_band_event_ids":sorted(band_ids),"upper_face":UPPER_FACE,"lower_face":LOWER_FACE,
      "boundaries":{"upper_corridor_facing":upper_corridor,"upper_opposite":upper_outer,"lower_corridor_facing":lower_corridor,"lower_opposite":lower_outer},
      "face_region_depth":{"upper_corridor_internal_mean_y":upper_corr_y,"upper_opposite_internal_mean_y":upper_outer_y,"H_UPPER":h_upper,"MOVE_UP":move_up,"lower_corridor_internal_mean_y":lower_corr_y,"lower_opposite_internal_mean_y":lower_outer_y,"H_LOWER":h_lower,"MOVE_DOWN":move_down},
      "upper_movable_vertices":upper_move,"lower_movable_vertices":lower_move,"shared_vertices_fixed":sorted(overlap),"vertex_classification":classification,
      "total_up_movement":move_up*len(upper_move),"total_down_movement":move_down*len(lower_move),"before_metrics":bm,"after_metrics":am,
      "before_geometry_diagnostics":bd,"after_geometry_diagnostics":ad,
      "central_gap_before":width_before,"central_gap_after":width_after,"central_gap_mean_increase":width_after["mean_boundary_vertex_to_opposite_polyline_separation"]-width_before["mean_boundary_vertex_to_opposite_polyline_separation"],"central_gap_minimum_nonincident_increase":width_after["minimum_nonincident_segment_separation"]-width_before["minimum_nonincident_segment_separation"],
      "original_event_outcomes":outcomes,"new_violations":new,"crowding_events_before":be,"crowding_events_after":ae,"crowding_clusters_before":bc,"crowding_clusters_after":ac,
      "edge_reconnection":reconnect,"hard_geometry_valid":hard_valid,"interpretation":interpretation,"coordinates_before":{v:list(p) for v,p in sorted(before.items())},"coordinates_after":{v:list(p) for v,p in sorted(final.items())}}
    output_dir.mkdir(parents=True,exist_ok=True);render(vertices,edges,before,final,report,output_dir/'iamp_vertical_halfspace_opening.svg');(output_dir/'iamp_vertical_halfspace_opening_report.json').write_text(json.dumps(report,indent=2,sort_keys=True)+'\n');return report

def main():
    p=argparse.ArgumentParser();b=Path('output/graph_first');p.add_argument('--graph',type=Path,default=b/'iamp_graph.json');p.add_argument('--corridor-report',type=Path,default=b/'iamp_corridor_separation_report.json');p.add_argument('--output',type=Path,default=b);a=p.parse_args();r=run(a.graph,a.corridor_report,a.output)
    print('Y_CENTER',r['Y_CENTER']);print('H/MOVE',r['face_region_depth']);print('moved',r['upper_movable_vertices'],r['lower_movable_vertices']);print('before',r['before_metrics']);print('after',r['after_metrics']);print(r['interpretation'])
if __name__=='__main__':main()
