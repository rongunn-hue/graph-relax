#!/usr/bin/env python3
"""Threshold-aware frozen-cell centroid pass using total deficiency acceptance."""
import argparse,json,math
from pathlib import Path
import graph_first_experiment_voronoi_centroid_move as core

FRACTION=.8
IMPORTANT=("net:EB","net:+9V","net:-9V","component:U1","net:NODE_B","component:J_PWR","component:J_EB","component:TPB","component:C1","net:EA","component:J_EA","net:GND","net:R1_EXT","component:J_OUT")
def deficiency(rows,threshold):return sum(max(0.,threshold-r["voronoi_area"]) for r in rows)

def render(vertices,edges,co,cells,rows,attempts,active,box,path):
 xmin,xmax,ymin,ymax=box;w=xmax-xmin;h=ymax-ymin;rat={r["vertex_id"]:r["crowding_ratio"] for r in rows};a={x["vertex_id"]:x for x in attempts};z=[f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="{xmin:.17g} {ymin:.17g} {w:.17g} {h:.17g}" width="1500" height="1000">',f'<rect x="{xmin:.17g}" y="{ymin:.17g}" width="{w:.17g}" height="{h:.17g}" fill="white" stroke="#111" stroke-width="2"/>']
 for i,v in enumerate(sorted(cells)):z.append('<polygon points="'+' '.join(f'{x:.17g},{y:.17g}' for x,y in cells[v])+f'" fill="hsl({(i*137)%360} 55% 70%)" fill-opacity=".25" stroke="#78909c" stroke-width=".8"/>')
 for v,x in a.items():
  p=x["starting_coordinate_at_attempt"];q=x["frozen_cell_centroid"];color="#2e7d32" if x["decision"]=="ACCEPTED" else "#c62828" if not x["graph_valid"] else "#7b1fa2"
  z.append(f'<line x1="{p[0]:.17g}" y1="{p[1]:.17g}" x2="{q[0]:.17g}" y2="{q[1]:.17g}" stroke="{color}" stroke-width="2" stroke-dasharray="4 4"/>')
  if x["decision"]!="ACCEPTED":z.append(f'<path d="M {q[0]-5:.17g} {q[1]-5:.17g} L {q[0]+5:.17g} {q[1]+5:.17g} M {q[0]-5:.17g} {q[1]+5:.17g} L {q[0]+5:.17g} {q[1]-5:.17g}" stroke="{color}" stroke-width="2"/>')
 for u,v in edges:z.append(f'<line x1="{co[u][0]:.17g}" y1="{co[u][1]:.17g}" x2="{co[v][0]:.17g}" y2="{co[v][1]:.17g}" stroke="#444" stroke-width="1.5"/>')
 for v in sorted(vertices):
  x,y=co[v];act=v in active;z.append(f'<circle cx="{x:.17g}" cy="{y:.17g}" r="{7 if act else 5}" fill="{"#ff8f00" if act else "#546e7a"}" stroke="white"/><text x="{x+8:.17g}" y="{y-8:.17g}" font-family="sans-serif" font-size="10">{vertices[v]["name"]} {"ACTIVE" if act else "FROZEN"} CR={rat[v]:.3f}</text>')
 z.append(f'<text x="{xmin+12:.17g}" y="{ymin+24:.17g}" font-family="sans-serif" font-size="17" font-weight="bold">80% deficiency-aware: green accepted, red graph reject, purple deficiency reject</text></svg>');path.write_text('\n'.join(z)+'\n')

def run(graph_path,pass1_path,previous_path,outdir):
 g=json.loads(graph_path.read_text());p1=json.loads(pass1_path.read_text());prev=json.loads(previous_path.read_text());vertices={v["id"]:v for v in g["vertices"]};edges=[(e["source"],e["target"]) for e in g["edges"]]
 initial={v:tuple(p) for v,p in p1["final_coordinates"].items()};box=tuple(p1["domain"][k] for k in ("xmin","xmax","ymin","ymax"));target=p1["target_area"];threshold=FRACTION*target
 frozen_cells=core.cells_in_fixed_box(initial,box);starting=core.records(frozen_cells,target);d0=core.diagnostics(initial,edges)
 if starting!=p1["after_density"] or any(d0[k] for k in ("proper_crossings","coincident_vertices","vertex_on_unrelated_edge_events")):raise AssertionError("baseline mismatch")
 active=[r["vertex_id"] for r in starting if r["voronoi_area"]<threshold];frozen=[r["vertex_id"] for r in starting if r["voronoi_area"]>=threshold];co=dict(initial);current_cells=frozen_cells;current_rows=starting;current_def=deficiency(starting,threshold);attempts=[]
 for row in starting:
  v=row["vertex_id"]
  if v not in active:continue
  start=co[v];cent=core.centroid(frozen_cells[v]);before=current_def;co[v]=cent;diag=core.diagnostics(co,edges);reasons=[]
  if diag["proper_crossings"]:reasons.append({"type":"proper_crossings","details":diag["proper_crossing_list"]})
  if diag["coincident_vertices"]:reasons.append({"type":"coincident_vertices","details":diag["coincident_vertex_pairs"]})
  if diag["vertex_on_unrelated_edge_events"]:reasons.append({"type":"vertex_on_unrelated_edge","details":diag["vertex_on_unrelated_edge_list"]})
  graph_valid=not reasons;proposed_def=None;change=None;improves=False
  if graph_valid:
   proposed_cells=core.cells_in_fixed_box(co,box);proposed_rows=core.records(proposed_cells,target);proposed_def=deficiency(proposed_rows,threshold);change=proposed_def-before;tol=1e-10*max(1.,before);improves=proposed_def<before-tol
  accepted=graph_valid and improves
  if accepted:current_cells=proposed_cells;current_rows=proposed_rows;current_def=proposed_def
  else:co[v]=start
  attempts.append({"attempt_order":len(attempts)+1,"vertex_id":v,"frozen_starting_area":row["voronoi_area"],"frozen_cell_centroid":list(cent),"starting_coordinate_at_attempt":list(start),"proposed_displacement_distance":math.dist(start,cent),"total_deficiency_before":before,"graph_valid":graph_valid,"graph_invalid_reasons":reasons,"proposed_total_deficiency":proposed_def,"deficiency_change":change,"deficiency_improving":improves,"decision":"ACCEPTED" if accepted else "REJECTED","rejection_category":None if accepted else "GRAPH_INVALID" if not graph_valid else "NON_IMPROVING_DEFICIENCY","complete_edge_regeneration_count":44})
 final_cells=core.cells_in_fixed_box(co,box);final=core.records(final_cells,target);final_def=deficiency(final,threshold);diag=core.diagnostics(co,edges);area_error=sum(r["voronoi_area"] for r in final)-p1["domain"]["area"]
 if abs(final_def-current_def)>1e-7 or abs(area_error)>1e-7 or any(diag[k] for k in ("proper_crossings","coincident_vertices","vertex_on_unrelated_edge_events")):raise AssertionError("invalid final")
 before={r["vertex_id"]:r for r in starting};old={x["vertex_id"]:x for x in prev["per_vertex"]};new={r["vertex_id"]:r for r in final};att={x["vertex_id"]:x for x in attempts};per=[]
 for v in sorted(vertices):
  a=att.get(v);per.append({"vertex_id":v,"pass1_area":before[v]["voronoi_area"],"previous_threshold80_area":old[v]["final_area"],"threshold_aware_area":new[v]["voronoi_area"],"pass1_crowding_ratio":before[v]["crowding_ratio"],"previous_threshold80_crowding_ratio":old[v]["final_crowding_ratio"],"threshold_aware_crowding_ratio":new[v]["crowding_ratio"],"pass1_rank":before[v]["rank"],"previous_threshold80_rank":old[v]["final_density_rank"],"threshold_aware_rank":new[v]["rank"],"classification":"ACTIVE" if v in active else "FROZEN","proposed":a is not None,"graph_valid":a["graph_valid"] if a else None,"deficiency_improving":a["deficiency_improving"] if a else None,"decision":a["decision"] if a else "NOT_ATTEMPTED"})
 def state(stat,rows,total):return {**stat,"cells_below_threshold":sum(r["voronoi_area"]<threshold for r in rows),"total_deficiency":total}
 prev_rows=[{"voronoi_area":x["final_area"]} for x in prev["per_vertex"]]
 accepted=[x for x in attempts if x["decision"]=="ACCEPTED"]
 report={"schema":"graph-relax.voronoi-threshold80-deficiency.v1","source_pass1_report":str(pass1_path),"fixed_domain":p1["domain"],"target_area":target,"threshold":threshold,"active_vertices":active,"frozen_vertices":frozen,"active_count":len(active),"frozen_count":len(frozen),"attempts":attempts,"accepted_moves":len(accepted),"rejected_graph_invalid":sum(x["rejection_category"]=="GRAPH_INVALID" for x in attempts),"rejected_non_improving_deficiency":sum(x["rejection_category"]=="NON_IMPROVING_DEFICIENCY" for x in attempts),"total_accepted_displacement":sum(x["proposed_displacement_distance"] for x in accepted),"maximum_accepted_displacement":max(x["proposed_displacement_distance"] for x in accepted),"starting_total_deficiency":deficiency(starting,threshold),"final_total_deficiency":final_def,"final_coordinates":{v:list(co[v]) for v in sorted(co)},"final_validation":diag,"comparison":{"PASS1_BASELINE":state(p1["after_statistics"],starting,deficiency(starting,threshold)),"PREVIOUS_THRESHOLD80":state(prev["after_statistics"],prev_rows,sum(max(0.,threshold-x["final_area"]) for x in prev["per_vertex"])),"THRESHOLD_AWARE":state(core.stats(final),final,final_def)},"per_vertex":per,"important_vertices":{x["vertex_id"]:x for x in per if x["vertex_id"] in IMPORTANT},"final_density":final,"partition_area_error":area_error}
 outdir.mkdir(parents=True,exist_ok=True);render(vertices,edges,co,final_cells,final,attempts,set(active),box,outdir/"iamp_voronoi_threshold80_deficiency.svg");(outdir/"iamp_voronoi_threshold80_deficiency_report.json").write_text(json.dumps(report,indent=2,sort_keys=True)+'\n');return report

def main():
 p=argparse.ArgumentParser();b=Path("output/graph_first");p.add_argument("--graph",type=Path,default=b/"iamp_graph.json");p.add_argument("--pass1",type=Path,default=b/"iamp_voronoi_centroid_move_report.json");p.add_argument("--previous",type=Path,default=b/"iamp_voronoi_threshold80_report.json");p.add_argument("--output",type=Path,default=b);a=p.parse_args();r=run(a.graph,a.pass1,a.previous,a.output);print("ACTIVE",r["active_vertices"]);print("deficiency",r["starting_total_deficiency"],"->",r["final_total_deficiency"],"accepted",r["accepted_moves"]);[print(x["attempt_order"],x["vertex_id"],x["decision"],x["rejection_category"],x["total_deficiency_before"],x["proposed_total_deficiency"]) for x in r["attempts"]]
if __name__=="__main__":main()
