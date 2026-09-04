#!/usr/bin/env python3
"""Measurement-only clipped Voronoi node-density diagnostic."""
import argparse,json,math,statistics
from pathlib import Path
import graph_first_experiment_crowding_analysis as crowding
import placement_geometry as geometry

def clip_halfplane(poly,a,b,c,tol=1e-10):
 """Keep ax+by<=c."""
 out=[]
 for p,q in zip(poly,poly[1:]+poly[:1]):
  fp=a*p[0]+b*p[1]-c;fq=a*q[0]+b*q[1]-c;pin=fp<=tol;qin=fq<=tol
  if pin:out.append(p)
  if pin!=qin:
   t=fp/(fp-fq);out.append((p[0]+t*(q[0]-p[0]),p[1]+t*(q[1]-p[1])))
 return out

def polygon_area(poly):
 return abs(sum(p[0]*q[1]-q[0]*p[1] for p,q in zip(poly,poly[1:]+poly[:1])))/2 if len(poly)>=3 else 0.

def cells(co):
 xs=[p[0] for p in co.values()];ys=[p[1] for p in co.values()];box=(min(xs),max(xs),min(ys),max(ys));xmin,xmax,ymin,ymax=box
 rectangle=[(xmin,ymin),(xmax,ymin),(xmax,ymax),(xmin,ymax)];result={}
 for v in sorted(co):
  x,y=co[v];poly=list(rectangle)
  for w in sorted(co):
   if w==v:continue
   u,z=co[w]
   # |p-v|^2 <= |p-w|^2 => 2(w-v).p <= |w|^2-|v|^2
   poly=clip_halfplane(poly,2*(u-x),2*(z-y),u*u+z*z-x*x-y*y)
   if not poly:raise AssertionError(f"empty Voronoi cell {v}")
  result[v]=poly
 return box,result

def render(vertices,edges,co,records,cellmap,box,out):
 xmin,xmax,ymin,ymax=box;width=xmax-xmin;height=ymax-ymin
 ratio={x["vertex_id"]:x["crowding_ratio"] for x in records}
 z=[f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="{xmin:.17g} {ymin:.17g} {width:.17g} {height:.17g}" width="1500" height="1000">',f'<rect x="{xmin:.17g}" y="{ymin:.17g}" width="{width:.17g}" height="{height:.17g}" fill="white" stroke="#111" stroke-width="2"/>']
 for i,v in enumerate(sorted(cellmap)):
  poly=cellmap[v];color=f'hsl({(i*137)%360} 55% 70%)'
  z.append('<polygon points="'+' '.join(f'{x:.17g},{y:.17g}' for x,y in poly)+f'" fill="{color}" fill-opacity=".28" stroke="#607d8b" stroke-width=".8"/>')
 for a,b in edges:z.append(f'<line x1="{co[a][0]:.17g}" y1="{co[a][1]:.17g}" x2="{co[b][0]:.17g}" y2="{co[b][1]:.17g}" stroke="#555" stroke-width="1.5"/>')
 for v in sorted(vertices):
  x,y=co[v];fill="#1565c0" if vertices[v]["type"]=="COMPONENT" else "#d84315"
  z.append(f'<circle cx="{x:.17g}" cy="{y:.17g}" r="6" fill="{fill}" stroke="white"/><text x="{x+8:.17g}" y="{y-8:.17g}" font-family="sans-serif" font-size="10">{vertices[v]["name"]}  CR={ratio[v]:.3f}</text>')
 z.append('</svg>');out.write_text("\n".join(z)+"\n")

def run(graph_path,base_path,outdir):
 g=json.loads(graph_path.read_text());base=json.loads(base_path.read_text());vertices={x["id"]:x for x in g["vertices"]};edges=[(x["source"],x["target"]) for x in g["edges"]];co={v:tuple(p) for v,p in base["coordinates_before"].items()};original=dict(co)
 diag=geometry.graph_geometry_diagnostics(co,edges)
 if (len(vertices),len(edges),diag["proper_unrelated_edge_crossing_count"])!=(34,44,0):raise AssertionError("bad source geometry")
 box,cellmap=cells(co);xmin,xmax,ymin,ymax=box;width=xmax-xmin;height=ymax-ymin;domain=width*height;target=domain/len(co);density=len(co)/domain
 areas={v:polygon_area(p) for v,p in cellmap.items()};total=sum(areas.values());error=total-domain
 if abs(error)>1e-8*max(1.,domain):raise AssertionError(f"Voronoi area partition error {error}")
 ordered=sorted(areas,key=lambda v:(areas[v],v));records=[]
 for rank,v in enumerate(ordered,1):
  a=areas[v];records.append({"rank":rank,"vertex_id":v,"voronoi_area":a,"area_over_target":a/target,"crowding_ratio":target/a,"percent_difference_from_target":100*(a/target-1),"interpretation":"denser / more crowded than average" if a<target else "less dense / more spatially open than average" if a>target else "average"})
 values=list(areas.values());stats={"minimum_area":min(values),"maximum_area":max(values),"median_area":statistics.median(values),"mean_area":statistics.mean(values),"population_standard_deviation":statistics.pstdev(values),"max_min_area_ratio":max(values)/min(values)}
 names=["component:TPB","component:U1","net:+9V","net:-9V","net:EA","net:EB","component:J_PWR","component:J_EA","component:J_EB","net:GND"]
 by={x["vertex_id"]:x for x in records};report={"schema":"graph-relax.global-node-density.v1","source_coordinates":base["source_coordinates"],"vertex_count":34,"electrical_incidence_count":44,"coordinates_unchanged":co==original,"proper_crossings":0,"domain":{"xmin":xmin,"xmax":xmax,"ymin":ymin,"ymax":ymax,"width":width,"height":height,"area":domain},"target_area":target,"target_density":density,"sum_clipped_cell_areas":total,"partition_area_error":error,"every_vertex_has_one_cell":len(cellmap)==34,"negative_cell_count":sum(a<0 for a in values),"ranked_vertices":records,"requested_sanity_vertices":{v:by[v] for v in names},"five_smallest_cells":[by[v] for v in ordered[:5]],"five_largest_cells":[by[v] for v in reversed(ordered[-5:])],"area_statistics":stats,"cells":{v:{"area":areas[v],"polygon":[list(p) for p in cellmap[v]]} for v in sorted(cellmap)}}
 outdir.mkdir(parents=True,exist_ok=True);render(vertices,edges,co,records,cellmap,box,outdir/"iamp_global_node_density.svg");(outdir/"iamp_global_node_density_report.json").write_text(json.dumps(report,indent=2,sort_keys=True)+"\n");return report

def main():
 p=argparse.ArgumentParser();b=Path("output/graph_first");p.add_argument("--graph",type=Path,default=b/"iamp_graph.json");p.add_argument("--base",type=Path,default=b/"iamp_corridor_separation_report.json");p.add_argument("--output",type=Path,default=b);a=p.parse_args();r=run(a.graph,a.base,a.output)
 print("domain",r["domain"],"target area",r["target_area"],"target density",r["target_density"],"sum error",r["partition_area_error"])
 for x in r["ranked_vertices"]:print(x["rank"],x["vertex_id"],x["voronoi_area"],x["area_over_target"],x["crowding_ratio"],x["percent_difference_from_target"])
if __name__=="__main__":main()
