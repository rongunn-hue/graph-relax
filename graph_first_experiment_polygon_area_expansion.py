#!/usr/bin/env python3
"""Reservoir-first deterministic polygon area expansion on the zero-crossing IAMP graph."""
import argparse,json,math
from collections import deque
from pathlib import Path
import graph_first_experiment_planar_region_inventory as inventory
import placement_geometry as geometry

TOL=1e-7; BISECT=60; MAX_DOUBLINGS=40
def area(w,co):return .5*sum(co[a][0]*co[b][1]-co[b][0]*co[a][1] for a,b in zip(w,w[1:]+w[:1]))
def metrics(w,co):
 vs=sorted(set(w));nearest={v:min(math.dist(co[v],co[u]) for u in vs if u!=v) for v in vs}
 return {"minimum_nearest_boundary_node_distance":min(nearest.values()),"average_nearest_boundary_node_distance":sum(nearest.values())/len(nearest),"nearest_by_vertex":nearest}
def directions(w,co,floating):
 sign=1 if area(w,co)>0 else -1;out={};n=len(w)
 for i,v in enumerate(w):
  if v not in floating:continue
  vectors=[]
  for a,b in ((co[w[i-1]],co[v]),(co[v],co[w[(i+1)%n]])):
   dx,dy=b[0]-a[0],b[1]-a[1];q=math.hypot(dx,dy)
   # Bounded walks are CCW: right normal is outward.
   vectors.append((sign*dy/q,-sign*dx/q))
  x,y=vectors[0][0]+vectors[1][0],vectors[0][1]+vectors[1][1];q=math.hypot(x,y)
  out[v]=(x/q,y/q) if q>TOL else vectors[0]
 return out
def shell_simple(w,co):
 es=list(zip(w,w[1:]+w[:1]))
 return not any(geometry.proper_segment_crossing(co[a],co[b],co[c],co[d]) for i,(a,b) in enumerate(es) for c,d in es[i+1:] if not {a,b}&{c,d})
def proposed(co,dirs,t):
 x=dict(co)
 for v,(dx,dy) in dirs.items():x[v]=(co[v][0]+t*dx,co[v][1]+t*dy)
 return x
def valid(graph_edges,w,base_area,base_metrics,co):
 diag=geometry.graph_geometry_diagnostics(co,graph_edges);re=geometry.validate_rebuilt_edge_segments(co,graph_edges,geometry.rebuild_straight_edge_segments(co,graph_edges));m=metrics(w,co);a=abs(area(w,co))
 ok=(re["valid"] and not diag["proper_unrelated_edge_crossing_count"] and not diag["coincident_vertex_pair_count"] and not diag["vertex_on_unrelated_edge_interior_count"] and shell_simple(w,co) and a>base_area+TOL and m["minimum_nearest_boundary_node_distance"]>base_metrics["minimum_nearest_boundary_node_distance"]+TOL and m["average_nearest_boundary_node_distance"]>base_metrics["average_nearest_boundary_node_distance"]+TOL)
 return ok,{"diagnostics":diag,"reconnection":re,"area":a,"metrics":m}
def expand(edges,face,co):
 w=face["ordered_boundary_walk"];floating=face["free_vertices"];fixed=[v for v in sorted(set(w)) if v not in floating];bm=metrics(w,co);ba=abs(area(w,co))
 if not floating:return co,{"status":"NO_FLOATING_VERTICES","movement":0.0,"floating_vertices":[],"fixed_vertices":fixed,"area_before":ba,"area_after":ba,"metrics_before":bm,"metrics_after":bm}
 dirs=directions(w,co,floating);step=max(bm["minimum_nearest_boundary_node_distance"],1e-6);low=0.;high=step;last=None
 for _ in range(MAX_DOUBLINGS):
  ok,info=valid(edges,w,ba,bm,proposed(co,dirs,high))
  if not ok:last=info;break
  low=high;high*=2
 else:raise AssertionError("expansion did not encounter first geometric/crowding limit")
 for _ in range(BISECT):
  mid=(low+high)/2;ok,info=valid(edges,w,ba,bm,proposed(co,dirs,mid))
  if ok:low=mid
  else:high=mid;last=info
 final=proposed(co,dirs,low);ok,fi=valid(edges,w,ba,bm,final) if low>TOL else (False,{"area":ba,"metrics":bm,"diagnostics":geometry.graph_geometry_diagnostics(co,edges)})
 stop=[]
 if last:
  d=last["diagnostics"]
  if d["proper_unrelated_edge_crossing_count"]:stop.append("PROPER_CROSSING")
  if d["coincident_vertex_pair_count"]:stop.append("COINCIDENCE")
  if d["vertex_on_unrelated_edge_interior_count"]:stop.append("VERTEX_ON_UNRELATED_EDGE")
  if last["area"]<=ba+TOL:stop.append("AREA_NOT_INCREASING")
  if last["metrics"]["minimum_nearest_boundary_node_distance"]<=bm["minimum_nearest_boundary_node_distance"]+TOL:stop.append("MINIMUM_BOUNDARY_SPACING_NOT_IMPROVING")
  if last["metrics"]["average_nearest_boundary_node_distance"]<=bm["average_nearest_boundary_node_distance"]+TOL:stop.append("AVERAGE_BOUNDARY_SPACING_NOT_IMPROVING")
 return final,{"status":"EXPANDED" if low>TOL else "BLOCKED","movement":low,"floating_vertices":floating,"fixed_vertices":fixed,"directions":dirs,"area_before":ba,"area_after":fi["area"],"metrics_before":bm,"metrics_after":fi["metrics"],"first_limit":sorted(set(stop)),"crossing_limited":"PROPER_CROSSING" in stop,"fixed_coordinates_unchanged":all(final[v]==co[v] for v in fixed)}
def dependencies(faces,ranked,targets):
 by={f["face_id"]:f for f in faces};ext=next(f["face_id"] for f in faces if not f["is_bounded"]);open_ids={x["face_id"] for x in ranked[len(ranked)*3//4:]};paths=[]
 for target in targets:
  q=deque([(target,[target])]);seen={target};found=None
  while q:
   cur,path=q.popleft()
   if cur!=target and (cur in open_ids or cur==ext):found=path;break
   for n in sorted(x["face_id"] for x in by[cur]["neighboring_faces"]):
    if n not in seen:seen.add(n);q.append((n,path+[n]))
  paths.append({"target":target,"crowded_to_reservoir_path":found,"processing_order":[x for x in reversed(found or []) if x!=ext]})
 order=[]
 for depth in range(max(map(lambda p:len(p["processing_order"]),paths),default=0)):
  for p in paths:
   if depth<len(p["processing_order"]) and p["processing_order"][depth] not in order:order.append(p["processing_order"][depth])
 return ext,sorted(open_ids),paths,order
def ranking(faces,co):
 rows=[]
 for f in faces:
  if f["is_bounded"] and f["simple_boundary"]:
   m=metrics(f["ordered_boundary_walk"],co);rows.append({"face_id":f["face_id"],"area":abs(area(f["ordered_boundary_walk"],co)),**{k:v for k,v in m.items() if k!="nearest_by_vertex"}})
 rows.sort(key=lambda x:(x["minimum_nearest_boundary_node_distance"],x["average_nearest_boundary_node_distance"],x["face_id"]))
 for i,x in enumerate(rows,1):x["crowding_rank"]=i
 return rows
def render(vertices,edges,co,faces,ranks,title,path):
 xs=[p[0] for p in co.values()];ys=[p[1] for p in co.values()];m=90;minx,miny=min(xs)-m,min(ys)-m;w=max(xs)-min(xs)+2*m;h=max(ys)-min(ys)+2*m;by={x["face_id"]:x for x in ranks}
 z=[f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="{minx} {miny} {w} {h}" width="1400" height="1000"><rect x="{minx}" y="{miny}" width="{w}" height="{h}" fill="white"/><text x="{minx+15}" y="{miny+25}" font-size="18">{title}</text>']
 for f in faces:
  if f["face_id"] in by:z.append('<polygon points="'+' '.join(f'{co[v][0]},{co[v][1]}' for v in f["ordered_boundary_walk"])+f'" fill="#ef5350" fill-opacity="{.22 if by[f["face_id"]]["crowding_rank"]<=3 else .05}" stroke="#aaa"/>')
 for a,b in edges:z.append(f'<line x1="{co[a][0]}" y1="{co[a][1]}" x2="{co[b][0]}" y2="{co[b][1]}" stroke="#666" stroke-width="1.4"/>')
 for v in sorted(vertices):
  x,y=co[v];z.append(f'<circle cx="{x}" cy="{y}" r="5" fill="#1565c0"/><text x="{x+7}" y="{y-7}" font-size="10">{vertices[v]["name"]}</text>')
 z.append('</svg>');path.write_text('\n'.join(z)+'\n')
def run(graph_path,direct_path,outdir):
 gd=json.loads(graph_path.read_text());dr=json.loads(direct_path.read_text());vertices={x["id"]:x for x in gd["vertices"]};edges=[(x["source"],x["target"]) for x in gd["edges"]];start={v:tuple(p) for v,p in dr["final_coordinates"].items()};sd=geometry.graph_geometry_diagnostics(start,edges)
 if (len(vertices),len(edges),sd["proper_unrelated_edge_crossing_count"])!=(34,44,0):raise AssertionError("wrong authoritative start")
 faces,darts=inventory.enumerate_faces(vertices,edges,start);inventory.add_classification_and_adjacency(faces,edges,darts);before=ranking(faces,start);targets=[x["face_id"] for x in before[:3]];ext,openids,paths,order=dependencies(faces,before,targets);by={f["face_id"]:f for f in faces};co=dict(start);ops=[]
 for fid in order:
  old=dict(co);co,result=expand(edges,by[fid],co);result.update({"face_id":fid,"boundary_walk":by[fid]["ordered_boundary_walk"],"crowding_rank_before":next(x["crowding_rank"] for x in before if x["face_id"]==fid)})
  result["adjacent_face_effects"]={n["face_id"]:{"before":metrics(by[n["face_id"]]["ordered_boundary_walk"],old),"after":metrics(by[n["face_id"]]["ordered_boundary_walk"],co)} for n in by[fid]["neighboring_faces"] if by[n["face_id"]]["is_bounded"] and by[n["face_id"]]["simple_boundary"]};ops.append(result)
 final=geometry.graph_geometry_diagnostics(co,edges);after=ranking(faces,co);re=geometry.validate_rebuilt_edge_segments(co,edges,geometry.rebuild_straight_edge_segments(co,edges));outdir.mkdir(parents=True,exist_ok=True)
 render(vertices,edges,start,faces,before,"IAMP polygon crowding — START",outdir/'iamp_polygon_area_expansion_start.svg');render(vertices,edges,co,faces,after,"IAMP polygon area expansion — FINAL",outdir/'iamp_polygon_area_expansion_final.svg')
 report={"schema":"graph-relax.polygon-area-expansion.v1","source_coordinates":f"{direct_path}#final_coordinates","vertex_count":34,"edge_count":44,"starting_crossings":0,"before_ranking":before,"crowded_target_faces":targets,"open_reservoir_faces":openids,"exterior_face":ext,"dependency_paths":paths,"processing_order":order,"operations":ops,"after_ranking":after,"final_validation":{"proper_crossings":final["proper_unrelated_edge_crossing_count"],"coincidences":final["coincident_vertex_pair_count"],"vertex_on_unrelated_edge_events":final["vertex_on_unrelated_edge_interior_count"],"vertices":34,"edges":44,"connectivity_unchanged":re["valid"]},"final_coordinates":{v:list(co[v]) for v in sorted(co)}};(outdir/'iamp_polygon_area_expansion_report.json').write_text(json.dumps(report,indent=2,sort_keys=True)+'\n');return report
def main():
 p=argparse.ArgumentParser();b=Path('output/graph_first');p.add_argument('--graph',type=Path,default=b/'iamp_graph.json');p.add_argument('--direct',type=Path,default=b/'iamp_nearness_rotation_direct_report.json');p.add_argument('--output',type=Path,default=b);a=p.parse_args();r=run(a.graph,a.direct,a.output);print(json.dumps({k:r[k] for k in ('crowded_target_faces','dependency_paths','processing_order','final_validation')},indent=2));[print(x['face_id'],x['status'],x['movement'],x['area_before'],x['area_after'],x['first_limit']) for x in r['operations']]
if __name__=='__main__':main()
