#!/usr/bin/env python3
"""Pure deterministic all-pairs 150-unit node separation experiment."""
import argparse,json,math
from itertools import combinations
from pathlib import Path
import graph_first_experiment_endpoint_rotation as loader
import graph_first_experiment_rotation_direct as unused
import placement_geometry as geometry

REQUIRED_SEPARATION=150.0
SAFETY_PASS_LIMIT=10000

def audit(co):
 rows=[{"vertices":[a,b],"distance":math.dist(co[a],co[b])} for a,b in combinations(sorted(co),2)]
 rows.sort(key=lambda x:(x["distance"],x["vertices"]));return rows

def separate(start):
 co=dict(start);nodes=sorted(co);pairs=list(combinations(nodes,2));corrections=[]
 for pass_number in range(1,SAFETY_PASS_LIMIT+1):
  pass_count=0
  for pair_index,(a,b) in enumerate(pairs):
   ax,ay=co[a];bx,by=co[b];dx,dy=bx-ax,by-ay;d=math.hypot(dx,dy)
   if d>=REQUIRED_SEPARATION:continue
   if d<=geometry.GEOMETRY_TOLERANCE:
    angle=2*math.pi*pair_index/len(pairs);ux,uy=math.cos(angle),math.sin(angle);method="stable pair-order direction"
   else:ux,uy=dx/d,dy/d;method="pair joining line"
   amount=REQUIRED_SEPARATION-d;half=amount/2
   co[a]=(ax-half*ux,ay-half*uy);co[b]=(bx+half*ux,by+half*uy);nudges=0
   # The exact mathematical correction can round one ulp below 150. Move both
   # endpoints outward by the smallest representable coordinate increments.
   while math.dist(co[a],co[b])<REQUIRED_SEPARATION:
    if nudges>=64:raise AssertionError("floating-point outward nudge did not reach 150")
    co[a]=(math.nextafter(co[a][0],-math.inf if ux>0 else math.inf),math.nextafter(co[a][1],-math.inf if uy>0 else math.inf));co[b]=(math.nextafter(co[b][0],math.inf if ux>0 else -math.inf),math.nextafter(co[b][1],math.inf if uy>0 else -math.inf));nudges+=1
   pass_count+=1
   corrections.append({"correction":len(corrections)+1,"pass":pass_number,"pair":[a,b],"distance_before":d,"required_distance_added":amount,"movement_each":half,"direction":[ux,uy],"method":method,"representable_outward_nudges":nudges})
  if pass_count==0:return co,corrections,pass_number
 remaining=[x for x in audit(co) if x["distance"]<REQUIRED_SEPARATION]
 raise RuntimeError(json.dumps({"safety_limit":SAFETY_PASS_LIMIT,"remaining":remaining},sort_keys=True))

def render_panel(graph,co,xoff,view,title):
 minx,miny,w,h=view;z=[f'<g transform="translate({xoff:.17g} 0)"><rect x="{minx:.17g}" y="{miny:.17g}" width="{w:.17g}" height="{h:.17g}" fill="white" stroke="#aaa"/>',f'<text x="{minx+18:.17g}" y="{miny+32:.17g}" font-family="sans-serif" font-size="22" font-weight="bold">{title}</text>']
 for u,v in graph.edges:z.append(f'<line x1="{co[u][0]:.17g}" y1="{co[u][1]:.17g}" x2="{co[v][0]:.17g}" y2="{co[v][1]:.17g}" stroke="#777" stroke-width="2"/>')
 for v in sorted(graph):
  x,y=co[v];c="#1565c0" if graph.nodes[v]["type"]=="COMPONENT" else "#d84315";z.append(f'<circle cx="{x:.17g}" cy="{y:.17g}" r="8" fill="{c}" stroke="white"/><text x="{x+11:.17g}" y="{y-11:.17g}" font-family="sans-serif" font-size="14">{graph.nodes[v]["name"]}</text>')
 z.append('</g>');return z

def render(graph,start,final,path,start_title="A. ORIGINAL PRE-UNCROSSING MDS",final_title="B. FINAL — MINIMUM NODE SEPARATION 150"):
 xs=[p[0] for p in final.values()];ys=[p[1] for p in final.values()];margin=100;view=(min(xs)-margin,min(ys)-margin,max(xs)-min(xs)+2*margin,max(ys)-min(ys)+2*margin);minx,miny,w,h=view;gap=80
 z=[f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="{minx:.17g} {miny:.17g} {2*w+gap:.17g} {h:.17g}" width="2200" height="1000">'];z+=render_panel(graph,start,0,view,start_title);z+=render_panel(graph,final,w+gap,view,final_title);z.append('</svg>');path.write_text('\n'.join(z)+'\n')

def run(graph_path,coordinates_path,mds_report_path,outdir):
 graph,start,start_diag,_=loader.load_start(graph_path,coordinates_path,mds_report_path);initial=audit(start);final,corrections,passes=separate(start);ending=audit(final);final_diag=geometry.graph_geometry_diagnostics(final,graph.edges);violations=[x for x in ending if x["distance"]<REQUIRED_SEPARATION]
 displacements={v:math.dist(start[v],final[v]) for v in sorted(start)};moved=[v for v,d in displacements.items() if d>geometry.GEOMETRY_TOLERANCE];unmoved=[v for v in sorted(start) if v not in moved]
 report={"schema":"graph-relax.pure-node-separation-150.v1","authoritative_coordinate_source":str(coordinates_path),"required_separation":REQUIRED_SEPARATION,"algorithm":"lexical all-pairs passes; exact equal split along each violating pair joining line","safety_pass_limit":SAFETY_PASS_LIMIT,"vertex_count":graph.number_of_nodes(),"electrical_edge_count":graph.number_of_edges(),"distinct_pairs_audited":len(ending),"starting_minimum_pairwise_distance":initial[0]["distance"],"final_minimum_pairwise_distance":ending[0]["distance"],"starting_pairs_below_150":sum(x["distance"]<REQUIRED_SEPARATION for x in initial),"final_pairs_below_150":len(violations),"remaining_violations":violations,"correction_count":len(corrections),"passes_required_including_final_clean_audit_pass":passes,"corrections":corrections,"nodes_moved":len(moved),"moved_vertices":moved,"unmoved_vertices":unmoved,"per_vertex_displacement":displacements,"total_node_displacement":sum(displacements.values()),"maximum_node_displacement":max(displacements.values()),"maximum_displacement_vertex":max(displacements,key=displacements.get),"starting_crossing_count":start_diag["proper_unrelated_edge_crossing_count"],"final_crossing_count":final_diag["proper_unrelated_edge_crossing_count"],"connectivity_unchanged":True,"rebuilt_edge_count":len(geometry.rebuild_straight_edge_segments(final,graph.edges)),"success":not violations and ending[0]["distance"]+geometry.GEOMETRY_TOLERANCE>=REQUIRED_SEPARATION,"starting_coordinates":{v:list(start[v]) for v in sorted(start)},"final_coordinates":{v:list(final[v]) for v in sorted(final)},"closest_final_pairs":ending[:12]}
 outdir.mkdir(parents=True,exist_ok=True);render(graph,start,final,outdir/"iamp_node_separation_150.svg");(outdir/"iamp_node_separation_150_report.json").write_text(json.dumps(report,indent=2,sort_keys=True)+'\n');return report

def main():
 p=argparse.ArgumentParser();b=Path("output/graph_first");p.add_argument("--graph",type=Path,default=b/"iamp_graph.json");p.add_argument("--coordinates",type=Path,default=b/"iamp_nearness_mds.json");p.add_argument("--mds-report",type=Path,default=b/"iamp_nearness_mds_report.json");p.add_argument("--output",type=Path,default=b);a=p.parse_args();r=run(a.graph,a.coordinates,a.mds_report,a.output);print(json.dumps({k:v for k,v in r.items() if k not in {"corrections","starting_coordinates","final_coordinates","per_vertex_displacement","closest_final_pairs"}},indent=2))
if __name__=="__main__":main()
