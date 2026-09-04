#!/usr/bin/env python3
"""One-pass IAMP crowding experiment with one node-spacing reference metric."""
import argparse, html, json, math
from itertools import combinations
from pathlib import Path

import graph_first_experiment_corrected_crowding as base
import graph_first_experiment_planar_region_inventory as inventory
import graph_first_experiment_polygon_area_expansion as polygon
import placement_geometry as geometry

TOL = 1e-7

def local_d(co, participants):
    """Mean leave-event-out nearest-node spacing for an event's vertices."""
    group = set(participants)
    vals = {v: min(math.dist(co[v], co[q]) for q in co if q not in group) for v in sorted(group)}
    return sum(vals.values()) / len(vals), vals

def event(kind, affected, d, D, where, **extra):
    return {"type": kind, "affected": affected, "d": d, "D": D,
            "ratio": d / D, "location": list(where), **extra}

def corrected_map(faces, co, edges):
    out = {k: [] for k in ("polygon", "node_node", "node_edge", "unrelated_edge_edge", "shared_endpoint_overlay")}
    ranks = polygon.ranking(faces, co)
    byface = {f["face_id"]: f for f in faces}
    for r in ranks:
        f = byface[r["face_id"]]; w = f["ordered_boundary_walk"]
        out["polygon"].append(event("POLYGON", r["face_id"], r["minimum_nearest_boundary_node_distance"],
            r["average_nearest_boundary_node_distance"], polygon_centroid(w, co),
            face_id=r["face_id"], boundary=w, area=r["area"], crowding_rank=r["crowding_rank"],
            minimum_NN=r["minimum_nearest_boundary_node_distance"], average_NN=r["average_nearest_boundary_node_distance"]))
    # One nearest node partner per node; pairs deduplicated.
    seen = set()
    for v in sorted(co):
        q = min((q for q in co if q != v), key=lambda q: (math.dist(co[v], co[q]), q)); pair = tuple(sorted((v, q)))
        if pair in seen: continue
        seen.add(pair); d = math.dist(co[v], co[q]); D, detail = local_d(co, pair)
        if d < D - TOL:
            out["node_node"].append(event("NODE_NODE", list(pair), d, D, midpoint(co[pair[0]], co[pair[1]]), local_node_spacing=detail))
    # Nearest eligible unrelated edge per node. Incident edges are absent by construction.
    all_ne = base.node_edge(co, edges)
    for v in sorted(co):
        x = min((x for x in all_ne if x["vertex"] == v), key=lambda x: (x["clearance"], x["edge"])); D, detail = local_d(co, [v, *x["edge"]])
        if x["clearance"] < D - TOL:
            out["node_edge"].append(event("NODE_EDGE", {"vertex": v, "edge": x["edge"]}, x["clearance"], D, x["closest_point"], closest_point=x["closest_point"], local_node_spacing=detail))
    # Nearest unrelated edge partner per edge; pairs deduplicated.
    all_ee = base.edge_edge(co, edges); seen = set()
    for e in sorted(tuple(sorted(x)) for x in edges):
        xs = [x for x in all_ee if e in (tuple(x["edges"][0]), tuple(x["edges"][1]))]
        x = min(xs, key=lambda x: (x["clearance"], x["edges"])); key = tuple(tuple(z) for z in x["edges"])
        if key in seen: continue
        seen.add(key); D, detail = local_d(co, sum(x["edges"], []))
        if x["clearance"] < D - TOL:
            p, q = x["closest_points"]; out["unrelated_edge_edge"].append(event("EDGE_EDGE", x["edges"], x["clearance"], D, midpoint(p, q), closest_points=x["closest_points"], local_node_spacing=detail))
    incidence={v:[] for v in co}
    for a,b in edges:incidence[a].append(b);incidence[b].append(a)
    raw=[]
    for v in sorted(co):
      for a,b in combinations(sorted(incidence[v]),2):
        s1=geometry.point_segment_distance(co[a],co[v],co[b]);s2=geometry.point_segment_distance(co[b],co[v],co[a]);d=min(s1,s2)
        ua=(co[a][0]-co[v][0],co[a][1]-co[v][1]);ub=(co[b][0]-co[v][0],co[b][1]-co[v][1]);la,lb=math.hypot(*ua),math.hypot(*ub);ang=math.degrees(math.acos(max(-1,min(1,(ua[0]*ub[0]+ua[1]*ub[1])/(la*lb)))))
        raw.append((d,v,a,b,ang,[la,lb]))
    for d,v,a,b,ang,lengths in sorted(raw):
        members = [v,a,b]; D, detail = local_d(co, members)
        if d < D - TOL:
            out["shared_endpoint_overlay"].append(event("SHARED_OVERLAY", [[v,a],[v,b]], d, D, co[v],
                shared_vertex=v, nonshared_endpoints=[a,b], angle_degrees=ang,
                available_segment_lengths=lengths, endpoint_to_other_segment_distances=[geometry.point_segment_distance(co[a],co[v],co[b]),geometry.point_segment_distance(co[b],co[v],co[a])],local_node_spacing=detail))
    for key in out:
        out[key].sort(key=lambda x: (x["ratio"], str(x["affected"])))
    n = 1
    for key in out:
        for x in out[key]: x["event_id"] = f"E{n:03d}"; n += 1
    return out

def midpoint(a,b): return [(a[0]+b[0])/2, (a[1]+b[1])/2]
def polygon_centroid(w,co): return [sum(co[v][0] for v in w)/len(w), sum(co[v][1] for v in w)/len(w)]
def hard_valid(co,edges):
    d=geometry.graph_geometry_diagnostics(co,edges); r=geometry.validate_rebuilt_edge_segments(co,edges,geometry.rebuild_straight_edge_segments(co,edges))
    return r["valid"] and not d["proper_unrelated_edge_crossing_count"] and not d["coincident_vertex_pair_count"] and not d["vertex_on_unrelated_edge_interior_count"],d

def direct_candidate(co,edges,v,target,score):
    c=dict(co);c[v]=target;ok,diag=hard_valid(c,edges);after=score(c)
    return c,ok and after>score(co)+TOL,after,{"hard_valid":ok,"diagnostics":diag}

def direct_polygon_expand(edges,face,co):
    w=face["ordered_boundary_walk"];floating=face["free_vertices"];fixed=[v for v in sorted(set(w)) if v not in floating];bm=polygon.metrics(w,co);ba=abs(polygon.area(w,co));dirs=polygon.directions(w,co,floating) if floating else {}
    # One scale-free candidate: outward displacement equals this polygon's current minimum boundary NN.
    t=bm["minimum_nearest_boundary_node_distance"];cand=polygon.proposed(co,dirs,t);ok,diag=hard_valid(cand,edges);am=abs(polygon.area(w,cand));mm=polygon.metrics(w,cand);accepted=bool(floating) and ok and polygon.shell_simple(w,cand) and am>ba+TOL and mm["minimum_nearest_boundary_node_distance"]>bm["minimum_nearest_boundary_node_distance"]+TOL and mm["average_nearest_boundary_node_distance"]>bm["average_nearest_boundary_node_distance"]+TOL
    return (cand if accepted else co),{"status":"ACCEPTED" if accepted else "REJECTED","single_shot":True,"candidate_count":1 if floating else 0,"target_displacement":t,"movement":t if accepted else 0.0,"floating_vertices":floating,"fixed_vertices":fixed,"directions":dirs,"area_before":ba,"area_after":am if accepted else ba,"metrics_before":bm,"metrics_after":mm if accepted else bm,"candidate_validation":diag,"rejection_reason":None if accepted else "HARD_INVALID_OR_AREA/CROWDING_NOT_IMPROVED"}

def node_pass(co,edges,m):
    ops=[]
    candidates=sorted(m["node_node"]+m["node_edge"],key=lambda x:(x["ratio"],x["event_id"]))
    for x in candidates:
        if x["type"]=="NODE_NODE":
            a,b=x["affected"]; choices=[]
            # Existing deterministic local rule: move lower-degree endpoint, then stable ID.
            degree={z:sum(z in e for e in edges) for z in (a,b)};v=sorted((a,b),key=lambda z:(degree[z],z))[0];q=b if v==a else a
            u=base.unit(co[v],co[q],v); score=lambda c,v=v,q=q:math.dist(c[v],c[q]); choices=[(v,u,score)]
        else:
            v=x["affected"]["vertex"]
            a,b=x["affected"]["edge"]; q=base.cp(co[v],co[a],co[b]); u=base.unit(co[v],q,v)
            score=lambda c,v=v,a=a,b=b:geometry.point_segment_distance(c[v],c[a],c[b]);choices=[(v,u,score)]
        participants = x["affected"] if x["type"]=="NODE_NODE" else [x["affected"]["vertex"], *x["affected"]["edge"]]
        D,_ = local_d(co, participants); before = choices[0][2](co) if choices else x["d"]
        v,u,score=choices[0];target=(co[v][0]+(D-before)*u[0],co[v][1]+(D-before)*u[1]);new,accepted,after,validation=direct_candidate(co,edges,v,target,score);disp=math.dist(co[v],target)
        if accepted:co=new;decision="ACCEPTED"
        else:after=before;decision="REJECTED"
        ops.append({"event_id":x["event_id"],"event_type":x["type"],"offending_geometry":x["affected"],"local_D":D,"d_before":before,"ratio_before":before/D,"selected_mover":v,"calculated_target":list(target),"displacement":disp,"d_after":after,"ratio_after":after/D,"decision":decision,"candidate_validation":validation,"candidate_count":1})
    return co,ops

def edge_pass(co,edges,m):
    ops=[]; candidates=sorted(m["unrelated_edge_edge"]+m["shared_endpoint_overlay"],key=lambda x:(x["ratio"],x["event_id"]))
    for x in candidates:
        choices=[]
        if x["type"]=="EDGE_EDGE":
            e,f=x["affected"]
            def score(c,e=e,f=f):p,q=base.segcp(c[e[0]],c[e[1]],c[f[0]],c[f[1]]);return math.dist(p,q)
            for v in e+f:
                p,q=base.segcp(co[e[0]],co[e[1]],co[f[0]],co[f[1]]); ref=q if v in e else p;choices.append((v,ref,score))
        else:
            s=x["shared_vertex"];a,b=x["nonshared_endpoints"]
            def score(c,s=s,a=a,b=b):return min(geometry.point_segment_distance(c[a],c[s],c[b]),geometry.point_segment_distance(c[b],c[s],c[a]))
            for v,o in ((a,b),(b,a)):
                q=base.cp(co[v],co[s],co[o]);choices.append((v,q,score))
        participants = sum(x["affected"],[]) if x["type"]=="EDGE_EDGE" else [x["shared_vertex"],*x["nonshared_endpoints"]]
        D,_ = local_d(co, participants); before = score(co)
        results=[]
        for v,q,score in choices:
            u=base.unit(co[v],q,v);target=(q[0]+D*u[0],q[1]+D*u[1]);new,accepted,after,validation=direct_candidate(co,edges,v,target,score)
            results.append({"mover":v,"target":list(target),"displacement":math.dist(co[v],target),"valid_and_improving":accepted,"d_after":after,"validation":validation,"coordinates":new})
        legal=[r for r in results if r["valid_and_improving"]]
        if legal:
            chosen=max(legal,key=lambda r:(r["d_after"],r["mover"]));co=chosen.pop("coordinates");v=chosen["mover"];after=chosen["d_after"];decision="ACCEPTED"
        else:v=None;after=before;decision="REJECTED"
        for r in results:r.pop("coordinates",None)
        ops.append({"event_id":x["event_id"],"event_type":x["type"],"offending_geometry":x["affected"],"local_D":D,"d_before":before,"ratio_before":before/D,"selected_mover":v,"d_after":after,"ratio_after":after/D,"decision":decision,"candidate_results":results,"candidate_count":len(results)})
    return co,ops

def audit(name,co,edges,m):
    inc=[list(e) for e in edges if name in e]
    relevant=[x for k in ("node_node","node_edge","shared_endpoint_overlay") for x in m[k] if name in str(x["affected"])]
    return {"vertex":name,"degree":len(inc),"incident_edges":inc,"flagged":bool(relevant),"events":relevant,
            "own_incident_edges_excluded_from_node_edge":all(name not in x["affected"]["edge"] for x in m["node_edge"] if x["affected"]["vertex"]==name)}

def render(vertices,edges,co,m,title,path):
    xs=[p[0] for p in co.values()];ys=[p[1] for p in co.values()];pad=50;panel=620;X=min(xs)-pad;Y=min(ys)-pad;W=max(xs)-min(xs)+2*pad;H=max(max(ys)-min(ys)+2*pad,420)
    colors={"POLYGON":"#ef6c00","NODE_NODE":"#7b1fa2","NODE_EDGE":"#00838f","EDGE_EDGE":"#d32f2f","SHARED_OVERLAY":"#2e7d32"};events=[x for k in m for x in m[k]]
    z=[f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="{X} {Y} {W+panel} {H}" width="1800" height="1100"><rect x="{X}" y="{Y}" width="{W+panel}" height="{H}" fill="white"/><text x="{X+10}" y="{Y+20}" font-size="16">{html.escape(title)}</text>']
    for a,b in edges:z.append(f'<line x1="{co[a][0]}" y1="{co[a][1]}" x2="{co[b][0]}" y2="{co[b][1]}" stroke="#999"/>')
    for x in events:
        px,py=x["location"];c=colors[x["type"]];z.append(f'<circle cx="{px}" cy="{py}" r="7" fill="none" stroke="{c}" stroke-width="2"/><text x="{px+8}" y="{py+4}" font-size="8" fill="{c}">{x["event_id"]}</text>')
    for v in sorted(vertices):x,y=co[v];z.append(f'<circle cx="{x}" cy="{y}" r="4" fill="#1565c0"/><text x="{x+5}" y="{y-5}" font-size="8">{html.escape(vertices[v]["name"])}</text>')
    tx=X+W+10;ty=Y+18;z.append(f'<text x="{tx}" y="{ty}" font-size="11">Legend: POLYGON orange; NODE_NODE purple; NODE_EDGE teal;</text>');ty+=14;z.append(f'<text x="{tx}" y="{ty}" font-size="11">EDGE_EDGE red; SHARED_OVERLAY green. All use d/D.</text>');ty+=18
    for x in events:
        affected=html.escape(str(x["affected"]));line=f'{x["event_id"]} {x["type"]} {affected} d={x["d"]:.3f} D={x["D"]:.3f} d/D={x["ratio"]:.4f}'
        z.append(f'<text x="{tx}" y="{ty}" font-size="8" fill="{colors[x["type"]]}">{line}</text>');ty+=10
    z.append(f'<!-- EVENT_ANNOTATION_COUNT={len(events)} --></svg>');path.write_text('\n'.join(z)+'\n');return len(events)

def run(graph_path,direct_path,outdir):
    gd=json.loads(graph_path.read_text());dr=json.loads(direct_path.read_text());V={x["id"]:x for x in gd["vertices"]};E=[(x["source"],x["target"]) for x in gd["edges"]];start={v:tuple(p) for v,p in dr["final_coordinates"].items()};ok,sd=hard_valid(start,E)
    if not ok or (len(V),len(E))!=(34,44):raise AssertionError("authoritative start invalid")
    faces,darts=inventory.enumerate_faces(V,E,start);inventory.add_classification_and_adjacency(faces,E,darts);initial=corrected_map(faces,start,E)
    tpo=audit("component:TPO",start,E,initial);tpb=audit("component:TPB",start,E,initial)
    # Reservoir-first expansion of the most crowded polygon only (ranking, no unit threshold).
    target=initial["polygon"][0]["face_id"];_,_,paths,order=polygon.dependencies(faces,polygon.ranking(faces,start),[target]);by={f["face_id"]:f for f in faces};co=dict(start);pops=[]
    for fid in order:
        co,r=direct_polygon_expand(E,by[fid],co);r.update({"face_id":fid,"boundary_walk":by[fid]["ordered_boundary_walk"]});pops.append(r)
    polygon_co=dict(co);after_polygon=corrected_map(faces,co,E)
    co,nops=node_pass(co,E,after_polygon);node_co=dict(co);after_nodes=corrected_map(faces,co,E)
    co,eops=edge_pass(co,E,after_nodes);final=corrected_map(faces,co,E);ok,fd=hard_valid(co,E)
    if not ok:raise AssertionError("final graph invalid")
    outdir.mkdir(parents=True,exist_ok=True);files=[]
    for name,c,m,title in (("iamp_common_metric_direct_start.svg",start,initial,"COMMON-METRIC DIRECT CROWDING — START"),("iamp_common_metric_direct_polygon.svg",polygon_co,after_polygon,"AFTER DIRECT POLYGON PASS"),("iamp_common_metric_direct_nodes.svg",node_co,after_nodes,"AFTER DIRECT NODE PASS"),("iamp_common_metric_direct_final.svg",co,final,"FINAL AFTER DIRECT EDGE PASS")):
        p=outdir/name;count=render(V,E,c,m,title,p);assert count==sum(len(x) for x in m.values());files.append(name)
    report={"schema":"graph-relax.common-node-spacing-crowding.v1","source_coordinates":f"{direct_path}#final_coordinates","reference_D_definition":"mean leave-event-out nearest-neighbor graph-node distance across the event's participating vertices; polygon D is its average boundary nearest-neighbor distance","starting_validation":{"vertices":34,"edges":44,"connectivity_unchanged":True,"proper_crossings":sd["proper_unrelated_edge_crossing_count"],"coincidences":sd["coincident_vertex_pair_count"],"vertex_on_unrelated_edge":sd["vertex_on_unrelated_edge_interior_count"]},"initial_map":initial,"audits":{"TPO":tpo,"TPB":tpb},"false_positive_rules":{"incident_node_edges_excluded":True,"shared_endpoint_zero_ignored":True},"false_positives_eliminated":[{"vertex":"component:TPO","excluded_geometry":["component:TPO","net:EOUT"],"reason":"TPO-EOUT is incident to TPO and therefore cannot be a NODE_EDGE event"},{"category":"SHARED_ENDPOINT_EDGE_OVERLAY","reason":"the common endpoint's zero distance is ignored; d is min(distance(A,V-B), distance(B,V-A))"}],"polygon_stage":{"target":target,"dependency_paths":paths,"processing_order":order,"operations":pops,"map_after":after_polygon},"node_stage":{"operations":nops,"map_after":after_nodes},"edge_stage":{"operations":eops},"final_map":final,"final_validation":{"vertices":34,"edges":44,"connectivity_unchanged":True,"proper_crossings":fd["proper_unrelated_edge_crossing_count"],"coincidences":fd["coincident_vertex_pair_count"],"vertex_on_unrelated_edge":fd["vertex_on_unrelated_edge_interior_count"]},"svg_files":files,"final_coordinates":{v:list(co[v]) for v in sorted(co)}}
    (outdir/"iamp_common_metric_direct_report.json").write_text(json.dumps(report,indent=2,sort_keys=True)+'\n');return report

def main():
    p=argparse.ArgumentParser();b=Path("output/graph_first");p.add_argument("--graph",type=Path,default=b/"iamp_graph.json");p.add_argument("--direct",type=Path,default=b/"iamp_nearness_rotation_direct_report.json");p.add_argument("--output",type=Path,default=b);a=p.parse_args();r=run(a.graph,a.direct,a.output)
    print(json.dumps({"audits":r["audits"],"counts_start":{k:len(v) for k,v in r["initial_map"].items()},"polygon_operations":len(r["polygon_stage"]["operations"]),"node_operations":len(r["node_stage"]["operations"]),"edge_operations":len(r["edge_stage"]["operations"]),"final_validation":r["final_validation"]},indent=2))
if __name__=="__main__":main()
