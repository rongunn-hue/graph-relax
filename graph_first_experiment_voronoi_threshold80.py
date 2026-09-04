#!/usr/bin/env python3
"""One 80%-thresholded frozen-cell Voronoi-centroid pass."""
import argparse,json,math
from pathlib import Path
import graph_first_experiment_voronoi_centroid_move as core

IMPORTANT=("net:EB","net:+9V","net:-9V","component:U1","net:NODE_B","component:J_PWR","component:J_EB","component:TPB","component:C1","net:EA","component:J_EA","net:GND","net:R1_EXT","component:J_OUT")

def render(vertices,edges,initial,final,final_cells,rows,attempts,active,box,path):
 xmin,xmax,ymin,ymax=box;w=xmax-xmin;h=ymax-ymin;ratio={r["vertex_id"]:r["crowding_ratio"] for r in rows};attempt={a["vertex_id"]:a for a in attempts}
 z=[f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="{xmin:.17g} {ymin:.17g} {w:.17g} {h:.17g}" width="1500" height="1000">',f'<rect x="{xmin:.17g}" y="{ymin:.17g}" width="{w:.17g}" height="{h:.17g}" fill="white" stroke="#111" stroke-width="2"/>']
 for i,v in enumerate(sorted(final_cells)):z.append('<polygon points="'+' '.join(f'{x:.17g},{y:.17g}' for x,y in final_cells[v])+f'" fill="hsl({(i*137)%360} 55% 70%)" fill-opacity=".25" stroke="#78909c" stroke-width=".8"/>')
 for v,a in attempt.items():
  p=a["starting_coordinate_at_attempt"];q=a["frozen_cell_centroid"];ok=a["decision"]=="ACCEPTED";color="#2e7d32" if ok else "#c62828"
  z.append(f'<line x1="{p[0]:.17g}" y1="{p[1]:.17g}" x2="{q[0]:.17g}" y2="{q[1]:.17g}" stroke="{color}" stroke-width="2" stroke-dasharray="{5 if ok else 2} 4"/>')
  if not ok:z.append(f'<path d="M {q[0]-5:.17g} {q[1]-5:.17g} L {q[0]+5:.17g} {q[1]+5:.17g} M {q[0]-5:.17g} {q[1]+5:.17g} L {q[0]+5:.17g} {q[1]-5:.17g}" stroke="#c62828" stroke-width="2"/>')
 for a,b in edges:z.append(f'<line x1="{final[a][0]:.17g}" y1="{final[a][1]:.17g}" x2="{final[b][0]:.17g}" y2="{final[b][1]:.17g}" stroke="#444" stroke-width="1.5"/>')
 for v in sorted(vertices):
  x,y=final[v];isactive=v in active;fill="#ff8f00" if isactive else "#546e7a";stroke="#111" if isactive else "white"
  z.append(f'<circle cx="{x:.17g}" cy="{y:.17g}" r="{7 if isactive else 5:.17g}" fill="{fill}" stroke="{stroke}"/><text x="{x+8:.17g}" y="{y-8:.17g}" font-family="sans-serif" font-size="10">{vertices[v]["name"]} {"ACTIVE" if isactive else "FROZEN"} CR={ratio[v]:.3f}</text>')
 z.append(f'<text x="{xmin+12:.17g}" y="{ymin+24:.17g}" font-family="sans-serif" font-size="18" font-weight="bold">Threshold 80% — final cells (orange=ACTIVE, gray=FROZEN)</text></svg>');path.write_text('\n'.join(z)+'\n')

def run(graph_path,pass1_path,outdir):
 g=json.loads(graph_path.read_text());p1=json.loads(pass1_path.read_text());vertices={v["id"]:v for v in g["vertices"]};edges=[(e["source"],e["target"]) for e in g["edges"]]
 initial={v:tuple(p) for v,p in p1["final_coordinates"].items()};box=tuple(p1["domain"][k] for k in ("xmin","xmax","ymin","ymax"));target=p1["target_area"];threshold=.8*target
 d0=core.diagnostics(initial,edges);frozen_cells=core.cells_in_fixed_box(initial,box);starting=core.records(frozen_cells,target)
 if starting!=p1["after_density"] or any(d0[k] for k in ("proper_crossings","coincident_vertices","vertex_on_unrelated_edge_events")):raise AssertionError("Pass-1 baseline mismatch")
 active=[r["vertex_id"] for r in starting if r["voronoi_area"]<threshold];frozen=[r["vertex_id"] for r in starting if r["voronoi_area"]>=threshold];co=dict(initial);attempts=[]
 for row in starting:
  v=row["vertex_id"]
  if v not in active:continue
  start=co[v];c=core.centroid(frozen_cells[v]);co[v]=c;d=core.diagnostics(co,edges);reasons=[]
  if d["proper_crossings"]:reasons.append({"type":"proper_crossings","details":d["proper_crossing_list"]})
  if d["coincident_vertices"]:reasons.append({"type":"coincident_vertices","details":d["coincident_vertex_pairs"]})
  if d["vertex_on_unrelated_edge_events"]:reasons.append({"type":"vertex_on_unrelated_edge","details":d["vertex_on_unrelated_edge_list"]})
  decision="REJECTED" if reasons else "ACCEPTED"
  if reasons:co[v]=start
  attempts.append({"attempt_order":len(attempts)+1,"vertex_id":v,"starting_area":row["voronoi_area"],"starting_crowding_ratio":row["crowding_ratio"],"starting_coordinate_at_attempt":list(start),"frozen_cell_centroid":list(c),"proposed_displacement_distance":math.dist(start,c),"decision":decision,"rejection_reasons":reasons,"complete_edge_regeneration_count":44})
 final_diag=core.diagnostics(co,edges);final_cells=core.cells_in_fixed_box(co,box);final=core.records(final_cells,target);area_error=sum(x["voronoi_area"] for x in final)-p1["domain"]["area"]
 if any(final_diag[k] for k in ("proper_crossings","coincident_vertices","vertex_on_unrelated_edge_events")) or abs(area_error)>1e-7:raise AssertionError("invalid threshold final")
 before={r["vertex_id"]:r for r in starting};after={r["vertex_id"]:r for r in final};attempt_by={a["vertex_id"]:a for a in attempts};per=[]
 for v in sorted(vertices):
  cls="ACTIVE" if v in active else "FROZEN";a=attempt_by.get(v);per.append({"vertex_id":v,"pass1_starting_area":before[v]["voronoi_area"],"classification":cls,"centroid_attempted":a is not None,"decision":a["decision"] if a else "NOT_ATTEMPTED","final_area":after[v]["voronoi_area"],"final_crowding_ratio":after[v]["crowding_ratio"],"final_density_rank":after[v]["rank"]})
 accepted=[a for a in attempts if a["decision"]=="ACCEPTED"]
 report={"schema":"graph-relax.voronoi-threshold80.v1","source_pass1_report":str(pass1_path),"domain":p1["domain"],"target_area":target,"minimum_acceptable_area":threshold,"threshold_fraction":.8,"vertex_count":34,"electrical_incidence_count":44,"active_vertices":active,"frozen_vertices":frozen,"active_count":len(active),"frozen_count":len(frozen),"attempts":attempts,"accepted_active_moves":len(accepted),"rejected_active_moves":len(attempts)-len(accepted),"total_accepted_displacement":sum(a["proposed_displacement_distance"] for a in accepted),"maximum_accepted_displacement":max(a["proposed_displacement_distance"] for a in accepted),"final_coordinates":{v:list(co[v]) for v in sorted(co)},"final_validation":final_diag,"before_statistics":p1["after_statistics"],"after_statistics":core.stats(final),"cells_below_threshold_before":len(active),"cells_below_threshold_after":sum(r["voronoi_area"]<threshold for r in final),"cells_at_or_above_threshold_after":sum(r["voronoi_area"]>=threshold for r in final),"per_vertex":per,"important_vertices":{x["vertex_id"]:x for x in per if x["vertex_id"] in IMPORTANT},"final_density":final,"partition_area_error":area_error}
 outdir.mkdir(parents=True,exist_ok=True);render(vertices,edges,initial,co,final_cells,final,attempts,set(active),box,outdir/"iamp_voronoi_threshold80.svg");(outdir/"iamp_voronoi_threshold80_report.json").write_text(json.dumps(report,indent=2,sort_keys=True)+'\n');return report

def main():
 p=argparse.ArgumentParser();b=Path("output/graph_first");p.add_argument("--graph",type=Path,default=b/"iamp_graph.json");p.add_argument("--pass1",type=Path,default=b/"iamp_voronoi_centroid_move_report.json");p.add_argument("--output",type=Path,default=b);a=p.parse_args();r=run(a.graph,a.pass1,a.output);print("ACTIVE",r["active_vertices"]);print("accepted",r["accepted_active_moves"],"rejected",r["rejected_active_moves"]);[print(x["attempt_order"],x["vertex_id"],x["decision"],x["rejection_reasons"]) for x in r["attempts"]]
if __name__=="__main__":main()
