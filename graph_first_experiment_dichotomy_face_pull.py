#!/usr/bin/env python3
"""Two-sided large-open-face pulls toward existing outward graph vertices."""
import argparse,json,math
from pathlib import Path
import graph_first_experiment_corridor_separation as corridor
import graph_first_experiment_face_poles as polygons
import graph_first_experiment_outward_boundary_wave as hard
import graph_first_experiment_planar_region_inventory as inventory
import placement_geometry as geometry

SIDES=(("SIDE_A","face:007","face:001"),("SIDE_B","face:008","face:004"))

def faces_now(vertices,edges,co):
 fs,d=inventory.enumerate_faces(vertices,edges,co);inventory.add_classification_and_adjacency(fs,edges,d);return {f["face_id"]:f for f in fs}

def shared_chain(face,inside_face):
 shared={tuple(x) for x in face["distinct_boundary_edges"]}&{tuple(x) for x in inside_face["distinct_boundary_edges"]}
 adjacency={}
 for a,b in shared:adjacency.setdefault(a,[]).append(b);adjacency.setdefault(b,[]).append(a)
 ends=sorted(v for v,n in adjacency.items() if len(n)==1)
 if len(ends)!=2:return None
 path=[ends[0]];previous=None
 while path[-1]!=ends[1]:
  options=sorted(x for x in adjacency[path[-1]] if x!=previous)
  if not options:return None
  previous,path_end=path[-1],options[0];path.append(path_end)
 return path

def point_path_distance(point,path,co):
 return min(geometry.point_segment_distance(point,co[a],co[b]) for a,b in zip(path,path[1:]))

def chord_inside_face(v,pole,face,co):
 a,b=co[v],co[pole];boundary=[tuple(x) for x in face["distinct_boundary_edges"]]
 for e in boundary:
  if v in e or pole in e:continue
  if geometry.proper_segment_crossing(a,b,co[e[0]],co[e[1]]):return False,"pull chord exits face across "+str(e)
 mid=((a[0]+b[0])/2,(a[1]+b[1])/2)
 pts=[co[x] for x in face["ordered_boundary_walk"]]
 # Deterministic odd-even containment for the chord midpoint.
 inside=False
 for x,y in zip(pts,pts[1:]+pts[:1]):
  if ((x[1]>mid[1])!=(y[1]>mid[1])) and mid[0] < (y[0]-x[0])*(mid[1]-x[1])/(y[1]-x[1])+x[0]:inside=not inside
 return (inside,"safe chord lies in selected face" if inside else "pull chord midpoint outside face")

def scaled(co,movers,pole,scale):
 out=dict(co);p=co[pole]
 for v in movers:out[v]=(p[0]+scale*(co[v][0]-p[0]),p[1]+scale*(co[v][1]-p[1]))
 return out

def maximum_pull(vertices,edges,co,movers,pole):
 target=scaled(co,movers,pole,.5);ok,diag,_,_,_,_,rec=hard.valid(vertices,edges,target)
 failures=[]
 if not ok:
  if diag["proper_unrelated_edge_crossing_count"]:failures.append({"type":"PROPER_EDGE_CROSSING","details":diag["proper_unrelated_edge_crossings"]})
  if diag["coincident_vertex_pair_count"]:failures.append({"type":"COINCIDENT_VERTICES","details":diag["coincident_vertex_pairs"]})
  if diag["vertex_on_unrelated_edge_interior_count"]:failures.append({"type":"VERTEX_ON_UNRELATED_EDGE","details":diag["vertices_on_unrelated_edge_interiors"]})
  if not rec["valid"]:failures.append({"type":"STALE_OR_CONNECTIVITY","details":rec["failures"]})
 if ok:return .5,target,True,failures
 lo=.5;hi=1.
 # Find the smallest valid scale (largest movement) on the same pull trajectory.
 for _ in range(80):
  mid=(lo+hi)/2
  if hard.valid(vertices,edges,scaled(co,movers,pole,mid))[0]:hi=mid
  else:lo=mid
 return hi,scaled(co,movers,pole,hi),False,failures

def render(vertices,edges,stages,out,clean):
 pts=[p for s in stages for p in s["coordinates"].values()];xs=[p[0] for p in pts];ys=[p[1] for p in pts];m=90;mnx,mxx,mny,mxy=min(xs),max(xs),min(ys),max(ys);count=1 if clean else len(stages);pw=mxx-mnx+2*m;h=mxy-mny+2*m
 z=[f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="{mnx-m:.17g} {mny-m:.17g} {count*pw:.17g} {h:.17g}" width="{1000*count}" height="900">',f'<rect x="{mnx-m:.17g}" y="{mny-m:.17g}" width="{count*pw:.17g}" height="{h:.17g}" fill="white"/>']
 shown=[stages[-1]] if clean else stages
 for i,s in enumerate(shown):
  co={v:tuple(p) for v,p in s["coordinates"].items()};z.append(f'<g transform="translate({i*pw:.17g} 0)">')
  for a,b in edges:z.append(f'<line data-edge="{a}|{b}" x1="{co[a][0]:.17g}" y1="{co[a][1]:.17g}" x2="{co[b][0]:.17g}" y2="{co[b][1]:.17g}" stroke="#777" stroke-width="1.5"/>')
  if not clean and s.get("pole"):
   p=co[s["pole"]];z.append(f'<circle cx="{p[0]:.17g}" cy="{p[1]:.17g}" r="10" fill="none" stroke="#6a1b9a" stroke-width="3"/>')
   for v in s["moved_vertices"]:z.append(f'<line x1="{p[0]:.17g}" y1="{p[1]:.17g}" x2="{co[v][0]:.17g}" y2="{co[v][1]:.17g}" stroke="#6a1b9a" stroke-dasharray="6 4"/>')
  moved=set(s.get("moved_vertices",[]))
  for v in sorted(vertices):
   x,y=co[v];fill="#2e7d32" if v in moved else "#1565c0" if vertices[v]["type"]=="COMPONENT" else "#d84315";z.append(f'<circle data-vertex="{v}" cx="{x:.17g}" cy="{y:.17g}" r="{8 if v in moved else 6}" fill="{fill}" stroke="white"/><text x="{x+8:.17g}" y="{y-8:.17g}" font-size="10">{vertices[v]["name"]}</text>')
  title="FINAL CLEAN" if clean else s["stage"];z.append(f'<text x="{mnx-m+12:.17g}" y="{mny-m+25:.17g}" font-size="17" font-weight="bold">{title}</text></g>')
 z.append("</svg>");out.write_text("\n".join(z)+"\n")

def run(graph_path,base_path,outdir):
 base=json.loads(base_path.read_text());g=json.loads(graph_path.read_text());vertices={x["id"]:x for x in g["vertices"]};edges=[(x["source"],x["target"]) for x in g["edges"]];co={v:tuple(p) for v,p in base["coordinates_before"].items()}
 ok,d,a,e,c,m,re=hard.valid(vertices,edges,co)
 if not ok or (len(vertices),len(edges),m["proper_crossings"])!=(34,44,0):raise AssertionError("invalid start")
 initial=faces_now(vertices,edges,co);areas=sorted(f["geometric_area"] for f in initial.values() if f["is_bounded"] and f["simple_boundary"]);large_floor=areas[len(areas)//2]
 consumed=set();moved=set();stages=[{"stage":"ORIGINAL","coordinates":{v:list(p) for v,p in sorted(co.items())},"metrics":m,"moved_vertices":[]}];trace=[];considered=[]
 for iteration,(side,fid,inward) in enumerate(SIDES,1):
  fs=faces_now(vertices,edges,co);face=fs[fid];chain=shared_chain(face,fs[inward]);area=face["geometric_area"];useful=area is not None and area>=large_floor and fid not in consumed
  record={"iteration":iteration,"side":side,"current_crowded_boundary_face":inward,"selected_open_polygon_face_id":fid,"polygon_area":area,"large_area_floor_median":large_floor,"sufficiently_open_useful":useful,"boundary_vertex_ids":face["distinct_boundary_vertices"],"facing_boundary_vertices":chain,"planned_pull_factor":.5}
  if not useful or not chain:
   record.update({"status":"BLOCKED","rejection_reason":"not a large unused simple polygon or no continuous shared facing boundary"});trace.append(record);considered.append(record);continue
  candidates=chain[1:-1];poles=[v for v in face["distinct_boundary_vertices"] if v not in chain]
  pole=max(poles,key=lambda v:(point_path_distance(co[v],chain,co),v));eligible=[];blocked={}
  for v in candidates:
   if v in moved:blocked[v]="vertex already moved";continue
   chord_ok,reason=chord_inside_face(v,pole,face,co)
   one=scaled(co,[v],pole,.5);hard_ok=hard.valid(vertices,edges,one)[0]
   if chord_ok and hard_ok:eligible.append(v)
   else:blocked[v]=reason if not chord_ok else "individual half-pull creates hard-invalid graph"
  record.update({"pull_pole_vertex_id":pole,"pull_pole_coordinates":list(co[pole]),"facing_boundary_candidate_vertex_ids":candidates,"safe_movable_vertex_ids":eligible,"blocked_vertex_ids_and_reasons":blocked})
  if not eligible:
   record.update({"status":"BLOCKED","rejection_reason":"no safely pullable facing-boundary vertices"});trace.append(record);considered.append(record);continue
  before=dict(co);scale,co,full,failures=maximum_pull(vertices,edges,co,eligible,pole);ok2,d2,a2,e2,c2,m2,re2=hard.valid(vertices,edges,co)
  if not ok2:raise AssertionError("accepted pull invalid")
  movements=[]
  for v in eligible:
   movements.append({"vertex":v,"original_coordinates":list(before[v]),"new_coordinates":list(co[v]),"actual_movement_distance":math.dist(before[v],co[v]),"full_0_5_pull_accepted":full,"accepted_scale":scale,"shortening_hard_validity_reason":failures or None});moved.add(v)
  consumed.add(fid);record.update({"status":"ACCEPTED","accepted_scale":scale,"full_0_5_pull_accepted":full,"shortening_hard_validity_reason":failures or None,"movements":movements,"crossings_after_reconnect":m2["proper_crossings"],"crowding_violations_after_reconnect":m2["clearance_violations"],"worst_clearance_after":m2["worst_clearance"],"metrics_after":m2,"edge_reconnection":re2});trace.append(record);considered.append(record)
  stages.append({"stage":f"AFTER {side} {fid}","coordinates":{v:list(p) for v,p in sorted(co.items())},"metrics":m2,"pole":pole,"moved_vertices":eligible})
 # Reevaluate outward neighbors; none qualifies as a further large, continuously facing unused polygon.
 final_faces=faces_now(vertices,edges,co)
 for side,fid,_ in SIDES:
  adjacent=[x["face_id"] for x in final_faces[fid]["neighboring_faces"] if final_faces[x["face_id"]]["is_bounded"] and x["face_id"] not in consumed]
  considered.append({"side":side,"from_consumed_face":fid,"reevaluated_adjacent_bounded_faces":adjacent,"status":"EXHAUSTED — NO NEXT UNUSED LARGE POLYGON WITH A CONTINUOUS CROWD-FACING BOUNDARY"})
 report={"schema":"graph-relax.dichotomy-face-pull.v1","source_coordinates":base["source_coordinates"],"vertex_count":34,"electrical_incidence_count":44,"principal_crowded_region_event_ids":base["corridor"]["events_included"],"side_definition":{"SIDE_A":"RMIN1-side topology","SIDE_B":"J_OUT-side topology"},"large_open_area_rule":"simple bounded face area >= median current simple bounded-face area","large_area_floor":large_floor,"trace":trace,"considered_polygons":considered,"consumed_faces":sorted(consumed),"stages":stages,"final_coordinates":stages[-1]["coordinates"],"final_metrics":stages[-1]["metrics"]}
 outdir.mkdir(parents=True,exist_ok=True);render(vertices,edges,stages,outdir/"iamp_dichotomy_face_pull.svg",False);render(vertices,edges,stages,outdir/"iamp_dichotomy_face_pull_clean.svg",True);(outdir/"iamp_dichotomy_face_pull_report.json").write_text(json.dumps(report,indent=2,sort_keys=True)+"\n");return report

def main():
 p=argparse.ArgumentParser();b=Path("output/graph_first");p.add_argument("--graph",type=Path,default=b/"iamp_graph.json");p.add_argument("--base",type=Path,default=b/"iamp_corridor_separation_report.json");p.add_argument("--output",type=Path,default=b);x=p.parse_args();r=run(x.graph,x.base,x.output)
 for t in r["trace"]:print(json.dumps(t,sort_keys=True))
if __name__=="__main__":main()
