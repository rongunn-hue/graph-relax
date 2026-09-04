#!/usr/bin/env python3
"""Deterministic outward wave over precomputed horizontal face boundaries."""
from __future__ import annotations
import argparse,json,math
from pathlib import Path

import graph_first_experiment_corridor_separation as corridor
import graph_first_experiment_planar_region_inventory as inventory
import placement_geometry as geometry

TOP_FACE="face:007";BOTTOM_FACE="face:008"
VALIDITY_SCAN_STEPS=4096
VALIDITY_BISECTION_STEPS=64

def y_values(path,x,coordinates):
    values=[]
    for a,b in zip(path,path[1:]):
        p,q=coordinates[a],coordinates[b];lo,hi=sorted((p[0],q[0]))
        if x < lo-geometry.GEOMETRY_TOLERANCE or x > hi+geometry.GEOMETRY_TOLERANCE:continue
        if abs(q[0]-p[0])<=geometry.GEOMETRY_TOLERANCE:
            values.extend((p[1],q[1]))
        else:
            t=(x-p[0])/(q[0]-p[0]);values.append(p[1]+t*(q[1]-p[1]))
    return values

def vertical_gap(current,next_path,coordinates):
    """Minimum vertical gap on the internal-vertex X overlap."""
    current_x=[coordinates[v][0] for v in current[1:-1]];next_x=[coordinates[v][0] for v in next_path[1:-1]]
    if not current_x or not next_x:return None
    lo=max(min(current_x),min(next_x));hi=min(max(current_x),max(next_x))
    if hi-lo<=geometry.GEOMETRY_TOLERANCE:return None
    candidates={lo,hi}
    for path in (current,next_path):
        for v in path[1:-1]:
            x=coordinates[v][0]
            if lo-geometry.GEOMETRY_TOLERANCE<=x<=hi+geometry.GEOMETRY_TOLERANCE:candidates.add(x)
    samples=[]
    for x in sorted(candidates):
        first=y_values(current,x,coordinates);second=y_values(next_path,x,coordinates)
        if first and second:samples.append({"x":x,"gap":min(abs(a-b) for a in first for b in second),"current_y":first,"next_y":second})
    if not samples:return None
    limiting=min(samples,key=lambda z:(z["gap"],z["x"]))
    return {"overlapping_x_range":[lo,hi],"minimum_vertical_separation":limiting["gap"],
            "mean_vertical_separation":sum(s["gap"] for s in samples)/len(samples),
            "limiting_sample":limiting,"samples":samples}

def prepare_sequences(vertices,edges,coordinates):
    faces,darts=inventory.enumerate_faces(vertices,edges,coordinates);inventory.add_classification_and_adjacency(faces,edges,darts);by={f["face_id"]:f for f in faces}
    def paths(face_id,inner_marker,inner_is_marker=True):
        ps=corridor.split_cycle(by[face_id]["ordered_boundary_walk"],"component:U1","net:GND")
        inner=next(p for p in ps if (inner_marker in p)==inner_is_marker);outer=next(p for p in ps if p is not inner)
        for p in (inner,outer):
            if p[0]!="component:U1":p.reverse()
        return inner,outer
    top0,top1=paths(TOP_FACE,"net:EA");bottom0,bottom1=paths(BOTTOM_FACE,"component:J_OUT",False)
    return {"top":[top0,top1],"bottom":[bottom0,bottom1],"top_face_sequence":["face:001","face:007"],"bottom_face_sequence":["face:004","face:008"],"constructed_from_original_embedding":True,"frozen_before_movement":True,
      "termination":{"top":"TOP_1 reaches exterior/branched adjacency; no next single bounded outward boundary","bottom":"BOTTOM_1 reaches exterior/branched adjacency; no next single bounded outward boundary"}}

def valid(vertices,edges,coordinates):
    d,a,e,c,m=corridor.metrics(vertices,edges,coordinates)
    ok=d["proper_unrelated_edge_crossing_count"]==0 and d["coincident_vertex_pair_count"]==0 and d["vertex_on_unrelated_edge_interior_count"]==0 and d["minimum_node_distance"]+geometry.GEOMETRY_TOLERANCE>=geometry.MIN_NODE_DISTANCE
    segments=geometry.rebuild_straight_edge_segments(coordinates,edges);re=geometry.validate_rebuilt_edge_segments(coordinates,edges,segments)
    return ok and re["valid"],d,a,e,c,m,re

def maximum_valid_move(vertices,edges,coordinates,movers,direction,half_gap):
    """Greatest valid translation in [0, half_gap], deterministically."""
    def candidate(distance):
        result=dict(coordinates)
        for v in movers:result[v]=(coordinates[v][0],coordinates[v][1]+direction*distance)
        return result
    def hard_ok(points):
        diagnostic=geometry.graph_geometry_diagnostics(points,edges)
        if (diagnostic["proper_unrelated_edge_crossing_count"] or diagnostic["coincident_vertex_pair_count"] or
                diagnostic["vertex_on_unrelated_edge_interior_count"] or
                diagnostic["minimum_node_distance"]+geometry.GEOMETRY_TOLERANCE<geometry.MIN_NODE_DISTANCE):return False
        segments=geometry.rebuild_straight_edge_segments(points,edges)
        return geometry.validate_rebuilt_edge_segments(points,edges,segments)["valid"]
    if not movers:return 0.0,candidate(0.0),{"cap_is_valid":True,"candidate_evaluations":1,"scan_steps":VALIDITY_SCAN_STEPS,"bisection_steps":0}
    evaluations=1;at_cap=candidate(half_gap)
    if hard_ok(at_cap):return half_gap,at_cap,{"cap_is_valid":True,"candidate_evaluations":evaluations,"scan_steps":VALIDITY_SCAN_STEPS,"bisection_steps":0}
    # Descend from the absolute cap to find the highest valid bracket. This
    # seeks the greatest valid final displacement, not an optimized crowding score.
    invalid_hi=half_gap;valid_lo=None
    for index in range(VALIDITY_SCAN_STEPS-1,-1,-1):
        distance=half_gap*index/VALIDITY_SCAN_STEPS;evaluations+=1
        if hard_ok(candidate(distance)):valid_lo=distance;break
        invalid_hi=distance
    if valid_lo is None:raise AssertionError("zero displacement unexpectedly invalid")
    for _ in range(VALIDITY_BISECTION_STEPS):
        mid=(valid_lo+invalid_hi)/2;evaluations+=1
        if hard_ok(candidate(mid)):valid_lo=mid
        else:invalid_hi=mid
    result=candidate(valid_lo)
    return valid_lo,result,{"cap_is_valid":False,"candidate_evaluations":evaluations,"scan_steps":VALIDITY_SCAN_STEPS,"bisection_steps":VALIDITY_BISECTION_STEPS,"first_invalid_above":invalid_hi}

def render(vertices,edges,stages,sequences,output):
    points=[tuple(p) for s in stages for p in s["coordinates"].values()];xs=[p[0] for p in points];ys=[p[1] for p in points];margin=90;minx,maxx,miny,maxy=min(xs),max(xs),min(ys),max(ys);pw=maxx-minx+2*margin;h=maxy-miny+2*margin
    parts=[f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="{minx-margin:.17g} {miny-margin:.17g} {len(stages)*pw:.17g} {h:.17g}" width="{900*len(stages)}" height="900">',f'<rect x="{minx-margin:.17g}" y="{miny-margin:.17g}" width="{len(stages)*pw:.17g}" height="{h:.17g}" fill="white"/>']
    for i,s in enumerate(stages):
        co={v:tuple(p) for v,p in s["coordinates"].items()};shift=i*pw;parts.append(f'<g id="stage-{i}" transform="translate({shift:.17g} 0)">')
        for a,b in edges:parts.append(f'<line data-edge="{a}|{b}" x1="{co[a][0]:.17g}" y1="{co[a][1]:.17g}" x2="{co[b][0]:.17g}" y2="{co[b][1]:.17g}" stroke="#777" stroke-width="1.5"/>')
        if i:
            side=stages[i]["operation"]["side"];color="#00897b" if side=="top" else "#8e24aa";path=stages[i]["operation"]["boundary"]
            parts.append('<polyline points="'+' '.join(f'{co[v][0]:.17g},{co[v][1]:.17g}' for v in path)+f'" fill="none" stroke="{color}" stroke-width="5"/>')
        for v in sorted(vertices):
            x,y=co[v];color="#1565c0" if vertices[v]["type"]=="COMPONENT" else "#d84315";parts.append(f'<circle data-vertex="{v}" cx="{x:.17g}" cy="{y:.17g}" r="6" fill="{color}" stroke="white"/><text x="{x+8:.17g}" y="{y-8:.17g}" font-size="10">{vertices[v]["name"]}</text>')
        op=s.get("operation",{});label="ORIGINAL" if not i else f'{op["side"].upper()}_{op["index"]}: move={op["actual_move"]:.6g}; max={op["max_valid_move"]:.6g}; half-gap={op["half_gap"]:.6g}'
        parts.append(f'<text x="{minx-margin+12:.17g}" y="{miny-margin+24:.17g}" font-size="16" font-weight="bold">{label}</text><text x="{minx-margin+12:.17g}" y="{miny-margin+45:.17g}" font-size="12">crossings={s["metrics"]["proper_crossings"]}; violations={s["metrics"]["clearance_violations"]}</text></g>')
    parts.append('</svg>');output.write_text('\n'.join(parts)+'\n')

def run(graph_path,corridor_report_path,output_dir):
    base=json.loads(corridor_report_path.read_text());g=json.loads(graph_path.read_text());vertices={x["id"]:x for x in g["vertices"]};edges=[(x["source"],x["target"]) for x in g["edges"]];original={v:tuple(p) for v,p in base["coordinates_before"].items()}
    ok,d,a,e,c,m,re=valid(vertices,edges,original)
    if not ok or (len(vertices),len(edges),m["proper_crossings"])!=(34,44,0):raise AssertionError("invalid starting state")
    sequences=prepare_sequences(vertices,edges,original);frozen=json.dumps(sequences,sort_keys=True);current=dict(original);moved=set();processed={"top":set(),"bottom":set()}
    top0=sequences["top"][0];bottom0=sequences["bottom"][0]
    initial_gap=vertical_gap(top0,bottom0,current)
    stages=[{"stage":"ORIGINAL","coordinates":{v:list(p) for v,p in sorted(current.items())},"metrics":m,"central_gap":{"minimum":initial_gap["minimum_vertical_separation"],"mean":initial_gap["mean_vertical_separation"]}}];operations=[]
    max_transitions=max(len(sequences["top"]),len(sequences["bottom"]))-1
    for k in range(max_transitions):
        for side,sign in (("top",-1),("bottom",1)):
            if k>=len(sequences[side])-1:continue
            before=dict(current);_,_,before_analysis,before_events,_,before_metrics,_=valid(vertices,edges,before)
            current_path,next_path=sequences[side][k],sequences[side][k+1];processed[side].add(k);gap=vertical_gap(current_path,next_path,before)
            if gap is None:
                operations.append({"side":side,"index":k,"face_id":TOP_FACE if side=="top" else BOTTOM_FACE,"status":"STOPPED — NO MEANINGFUL X OVERLAP","boundary":current_path,"next_boundary":next_path});continue
            half_gap=gap["minimum_vertical_separation"]/2;movers=[v for v in current_path[1:-1] if v not in moved]
            max_move,candidate,search=maximum_valid_move(vertices,edges,before,movers,sign,half_gap);actual=min(max_move,half_gap)
            if abs(actual-max_move)>geometry.GEOMETRY_TOLERANCE:raise AssertionError("actual move does not equal min(max valid, half gap)")
            current=candidate;moved.update(movers)
            ok2,d2,a2,e2,c2,m2,re2=valid(vertices,edges,current)
            if not ok2:raise AssertionError("limited move invalid")
            event_rows,new_events=corridor.event_changes(before_analysis,a2,before_events);gap_now=vertical_gap(top0,bottom0,current)
            face_sequence=sequences[f"{side}_face_sequence"]
            record={"side":side,"index":k,"face_id":face_sequence[k],"next_outward_face_id":face_sequence[k+1],"status":"ACCEPTED","boundary":current_path,"next_boundary":next_path,"anchors":[current_path[0],current_path[-1]],"non_anchor_vertices":current_path[1:-1],"eligible_unmoved_vertices":movers,"full_gap":gap["minimum_vertical_separation"],"half_gap":half_gap,"max_valid_move":max_move,"actual_move":actual,"limiting_validity_search":search,"crossings_after":m2["proper_crossings"],"clearance_violations_after":m2["clearance_violations"]}
            stage={"stage":f"AFTER {side.upper()}_{k}","coordinates":{v:list(p) for v,p in sorted(current.items())},"operation":record,"metrics":m2,"edge_reconnection":re2,"central_gap":{"minimum":gap_now["minimum_vertical_separation"],"mean":gap_now["mean_vertical_separation"]},"original_event_outcomes_this_move":event_rows,"new_violations_this_move":new_events}
            stages.append(stage);operations.append(record)
    assertions={"sequences_frozen":json.dumps(sequences,sort_keys=True)==frozen,"exterior_processed":False,"duplicate_boundary_processing":any(len(x)!=len(set(x)) for x in processed.values()),"vertex_moved_more_than_once":False,"x_coordinates_changed":any(current[v][0]!=original[v][0] for v in current),"final_edge_reconnection_valid":stages[-1].get("edge_reconnection",re)["valid"]}
    if not assertions["sequences_frozen"] or assertions["duplicate_boundary_processing"] or assertions["x_coordinates_changed"] or not assertions["final_edge_reconnection_valid"]:raise AssertionError(assertions)
    side_progress={side:{"available_transitions":len(sequences[side])-1,
                         "processed_transitions":len(processed[side]),
                         "accepted_transitions":sum(stage.get("side")==side and stage.get("status")=="ACCEPTED" for stage in operations),
                         "rejected_transitions":sum(stage.get("side")==side and stage.get("status","").startswith("REJECTED") for stage in operations),
                         "exhausted":len(processed[side])==len(sequences[side])-1}
                   for side in ("top","bottom")}
    report={"schema":"graph-relax.outward-boundary-wave.v1","source_coordinates":base["source_coordinates"],"vertex_count":len(vertices),"electrical_incidence_count":len(edges),"readability_clearance":geometry.READABILITY_CLEARANCE,"boundary_sequences":sequences,"sequence_assertions":assertions,"side_progress":side_progress,"stages":stages,"operations":operations,"moved_vertices":sorted(moved),"final_coordinates":stages[-1]["coordinates"]}
    output_dir.mkdir(parents=True,exist_ok=True);render(vertices,edges,stages,sequences,output_dir/'iamp_outward_boundary_wave.svg');(output_dir/'iamp_outward_boundary_wave_report.json').write_text(json.dumps(report,indent=2,sort_keys=True)+'\n');return report

def main():
    p=argparse.ArgumentParser();b=Path('output/graph_first');p.add_argument('--graph',type=Path,default=b/'iamp_graph.json');p.add_argument('--corridor-report',type=Path,default=b/'iamp_corridor_separation_report.json');p.add_argument('--output',type=Path,default=b);a=p.parse_args();r=run(a.graph,a.corridor_report,a.output);print(json.dumps(r['boundary_sequences'],indent=2));
    for s in r['stages']:print(s['stage'],s['metrics']['proper_crossings'],s['metrics']['clearance_violations'],s['central_gap'])
if __name__=='__main__':main()
