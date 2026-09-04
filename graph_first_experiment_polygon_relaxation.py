#!/usr/bin/env python3
"""Independent fixed-area, fixed-anchor perimeter relaxations of three IAMP shells."""
import argparse,json,math
from pathlib import Path
import numpy as np
from scipy.optimize import minimize
import graph_first_experiment_corridor_separation as corridor
import graph_first_experiment_planar_region_inventory as inventory
import placement_geometry as geometry

REGIONS={"R1_EXT":("face:007","net:R1_EXT"),"RA_SW":("face:006","net:RA_SW"),"RMIN4":("face:010","component:RMIN4")}

def area(points):return .5*sum(a[0]*b[1]-b[0]*a[1] for a,b in zip(points,points[1:]+points[:1]))
def perimeter(points):return sum(math.dist(a,b) for a,b in zip(points,points[1:]+points[:1]))

def simple_shell(face):
 walk=face["ordered_boundary_walk"];excluded=[]
 # Remove exact out-and-back bridge excursions A -> leaf -> A.
 result=[];i=0
 while i<len(walk):
  if i+2<len(walk) and walk[i]==walk[i+2]:
   excluded.append({"walk":[walk[i],walk[i+1],walk[i+2]],"bridge_vertex":walk[i+1]})
   if not result or result[-1]!=walk[i]:result.append(walk[i])
   i+=3
  else:
   result.append(walk[i]);i+=1
 # The cyclic rotation can put the repeated base across the list boundary.
 changed=True
 while changed:
  changed=False
  n=len(result)
  for j in range(n):
   if n>=3 and result[j]==result[(j+2)%n]:
    excluded.append({"walk":[result[j],result[(j+1)%n],result[(j+2)%n]],"bridge_vertex":result[(j+1)%n]})
    result=[result[k] for k in range(n) if k!=(j+1)%n];changed=True;break
 return result,excluded

def classify(shell,edges):
 inc={v:[] for v in shell}
 for a,b in edges:
  if a in inc:inc[a].append(tuple(sorted((a,b))))
  if b in inc:inc[b].append(tuple(sorted((a,b))))
 boundary={tuple(sorted((a,b))) for a,b in zip(shell,shell[1:]+shell[:1])}
 floating=[v for v in shell if len(inc[v])==2 and set(inc[v]).issubset(boundary)]
 return floating,[v for v in shell if v not in floating],inc

def solve(shell,floating,coordinates):
 original={v:coordinates[v] for v in shell};a0=area([original[v] for v in shell]);sign=1 if a0>=0 else -1
 x0=np.array([q for v in floating for q in original[v]],dtype=float)
 def points(x):
  co=dict(original)
  for i,v in enumerate(floating):co[v]=(float(x[2*i]),float(x[2*i+1]))
  return [co[v] for v in shell],co
 def objective(x):return perimeter(points(x)[0])
 def equality(x):return area(points(x)[0])-a0
 result=minimize(objective,x0,method="SLSQP",constraints=[{"type":"eq","fun":equality}],options={"ftol":1e-12,"maxiter":10000,"disp":False})
 pts,co=points(result.x)
 return co,{"method":"SciPy SLSQP; analytic objective evaluated deterministically from one authoritative start; signed-area equality constraint","success":bool(result.success),"message":result.message,"iterations":int(result.nit),"objective_evaluations":int(result.nfev),"signed_area_target":a0,"signed_area_result":area(pts),"perimeter_result":perimeter(pts)}

def shell_simple(shell,co):
 es=list(zip(shell,shell[1:]+shell[:1]))
 bad=[]
 for i,(a,b) in enumerate(es):
  for j,(c,d) in enumerate(es):
   if j<=i or {a,b}&{c,d}:continue
   if geometry.proper_segment_crossing(co[a],co[b],co[c],co[d]):bad.append([[a,b],[c,d]])
 return not bad,bad

def render(vertices,edges,original,results,out):
 allco=[original]+[x["coordinates"] for x in results];pts=[p for co in allco for p in co.values()];xs=[p[0] for p in pts];ys=[p[1] for p in pts];m=80;mnx,mxx,mny,mxy=min(xs),max(xs),min(ys),max(ys);pw=mxx-mnx+2*m;h=mxy-mny+2*m
 panels=[]
 for r in results:panels.extend([(r["region"]+" ORIGINAL",original,r),(r["region"]+" RELAXED",r["coordinates"],r)])
 z=[f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="{mnx-m:.17g} {mny-m:.17g} {len(panels)*pw:.17g} {h:.17g}" width="{650*len(panels)}" height="750">',f'<rect x="{mnx-m:.17g}" y="{mny-m:.17g}" width="{len(panels)*pw:.17g}" height="{h:.17g}" fill="white"/>']
 for i,(title,co,r) in enumerate(panels):
  z.append(f'<g transform="translate({i*pw:.17g} 0)">')
  for a,b in edges:z.append(f'<line x1="{co[a][0]:.17g}" y1="{co[a][1]:.17g}" x2="{co[b][0]:.17g}" y2="{co[b][1]:.17g}" stroke="#aaa" stroke-width="1"/>')
  z.append('<polygon points="'+' '.join(f'{co[v][0]:.17g},{co[v][1]:.17g}' for v in r["shell"])+f'" fill="none" stroke="#6a1b9a" stroke-width="4"/>')
  for v in sorted(vertices):
   x,y=co[v];fill="#2e7d32" if v in r["floating_vertices"] else "#c62828" if v in r["fixed_vertices"] else "#78909c";z.append(f'<circle cx="{x:.17g}" cy="{y:.17g}" r="{7 if v in r["shell"] else 3}" fill="{fill}"/><text x="{x+7:.17g}" y="{y-7:.17g}" font-size="9">{vertices[v]["name"]}</text>')
  z.append(f'<text x="{mnx-m+10:.17g}" y="{mny-m+22:.17g}" font-size="15" font-weight="bold">{title}</text></g>')
 z.append("</svg>");out.write_text("\n".join(z)+"\n")

def clean_svg(vertices,edges,co,r,out):
 # Reuse the comparison renderer with one result and retain only a compact standalone view.
 xs=[p[0] for p in co.values()];ys=[p[1] for p in co.values()];m=70;mnx,mxx,mny,mxy=min(xs),max(xs),min(ys),max(ys)
 z=[f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="{mnx-m:.17g} {mny-m:.17g} {mxx-mnx+2*m:.17g} {mxy-mny+2*m:.17g}" width="1000" height="850"><rect x="{mnx-m:.17g}" y="{mny-m:.17g}" width="{mxx-mnx+2*m:.17g}" height="{mxy-mny+2*m:.17g}" fill="white"/>']
 for a,b in edges:z.append(f'<line x1="{co[a][0]:.17g}" y1="{co[a][1]:.17g}" x2="{co[b][0]:.17g}" y2="{co[b][1]:.17g}" stroke="#777" stroke-width="1.5"/>')
 for v in sorted(vertices):
  x,y=co[v];z.append(f'<circle cx="{x:.17g}" cy="{y:.17g}" r="5" fill="#1565c0"/><text x="{x+7:.17g}" y="{y-7:.17g}" font-size="10">{vertices[v]["name"]}</text>')
 z.append(f'<text x="{mnx-m+10:.17g}" y="{mny-m+22:.17g}" font-size="17">{r["region"]} relaxed</text></svg>');out.write_text("\n".join(z)+"\n")

def run(graph_path,base_path,outdir):
 base=json.loads(base_path.read_text());g=json.loads(graph_path.read_text());vertices={x["id"]:x for x in g["vertices"]};edges=[(x["source"],x["target"]) for x in g["edges"]];original={v:tuple(p) for v,p in base["coordinates_before"].items()}
 faces,_=inventory.enumerate_faces(vertices,edges,original);by={f["face_id"]:f for f in faces};results=[]
 for region,(fid,marker) in REGIONS.items():
  face=by[fid]
  if marker not in face["distinct_boundary_vertices"]:raise AssertionError("face identity mismatch")
  shell,excursions=simple_shell(face);floating,fixed,inc=classify(shell,edges);relaxed_shell,opt=solve(shell,floating,original);co=dict(original);co.update(relaxed_shell)
  p0=[original[v] for v in shell];p1=[co[v] for v in shell];a0=abs(area(p0));a1=abs(area(p1));l0=perimeter(p0);l1=perimeter(p1);simple,bad=shell_simple(shell,co)
  diag=geometry.graph_geometry_diagnostics(co,edges);segments=geometry.rebuild_straight_edge_segments(co,edges);re=geometry.validate_rebuilt_edge_segments(co,edges,segments)
  fixed_checks={v:{"original":list(original[v]),"relaxed":list(co[v]),"unchanged":co[v]==original[v]} for v in fixed}
  movements={v:{"original":list(original[v]),"relaxed":list(co[v]),"displacement":math.dist(original[v],co[v])} for v in floating}
  r={"region":region,"face_id":fid,"formal_face_walk":face["ordered_boundary_walk"],"shell":shell,"excluded_bridge_excursions":excursions,"original_area":a0,"relaxed_area":a1,"area_error":a1-a0,"original_perimeter":l0,"relaxed_perimeter":l1,"absolute_perimeter_reduction":l0-l1,"percentage_perimeter_reduction":100*(l0-l1)/l0,"floating_vertices":floating,"fixed_vertices":fixed,"floating_movements":movements,"fixed_coordinate_checks":fixed_checks,"optimizer":opt,"relaxed_shell_simple":simple,"shell_crossings":bad,"whole_graph_diagnostics":diag,"edge_reconnection":re,"embedding_remained_valid":simple and diag["proper_unrelated_edge_crossing_count"]==0 and diag["coincident_vertex_pair_count"]==0 and diag["vertex_on_unrelated_edge_interior_count"]==0,"coordinates":co}
  results.append(r)
 outdir.mkdir(parents=True,exist_ok=True);render(vertices,edges,original,results,outdir/"iamp_polygon_relaxation_test.svg")
 for r in results:clean_svg(vertices,edges,r["coordinates"],r,outdir/f'iamp_relaxed_{r["region"]}.svg')
 report={"schema":"graph-relax.polygon-perimeter-relaxation.v1","source_coordinates":base["source_coordinates"],"vertex_count":34,"electrical_incidence_count":44,"independent_start_for_every_region":True,"method":"deterministic SLSQP constrained minimization of straight-edge perimeter with signed area equality and fixed anchor coordinates","results":[{k:v for k,v in r.items() if k!="coordinates"} for r in results]}
 (outdir/"iamp_polygon_relaxation_test_report.json").write_text(json.dumps(report,indent=2,sort_keys=True)+"\n");return report

def main():
 p=argparse.ArgumentParser();b=Path("output/graph_first");p.add_argument("--graph",type=Path,default=b/"iamp_graph.json");p.add_argument("--base",type=Path,default=b/"iamp_corridor_separation_report.json");p.add_argument("--output",type=Path,default=b);a=p.parse_args();r=run(a.graph,a.base,a.output)
 for x in r["results"]:print(x["region"],x["face_id"],x["shell"],"floating",x["floating_vertices"],"fixed",x["fixed_vertices"],"area",x["original_area"],x["relaxed_area"],"perimeter",x["original_perimeter"],x["relaxed_perimeter"],"crossings",x["whole_graph_diagnostics"]["proper_unrelated_edge_crossing_count"])
if __name__=="__main__":main()
