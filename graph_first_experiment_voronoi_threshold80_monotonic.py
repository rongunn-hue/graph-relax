#!/usr/bin/env python3
"""Crowding-monotonic 80%-threshold frozen-centroid experiment."""
import argparse,json,math
from pathlib import Path
import graph_first_experiment_voronoi_centroid_move as core
import graph_first_experiment_voronoi_threshold80_deficiency as aware

IMPORTANT=aware.IMPORTANT
def deficiencies(rows,t):return {r["vertex_id"]:max(0.,t-r["voronoi_area"]) for r in rows}

def render(vertices,edges,co,cells,rows,attempts,active,box,path):
 xmin,xmax,ymin,ymax=box;w=xmax-xmin;h=ymax-ymin;rat={r["vertex_id"]:r["crowding_ratio"] for r in rows};z=[f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="{xmin:.17g} {ymin:.17g} {w:.17g} {h:.17g}" width="1500" height="1000">',f'<rect x="{xmin:.17g}" y="{ymin:.17g}" width="{w:.17g}" height="{h:.17g}" fill="white" stroke="#111" stroke-width="2"/>']
 for i,v in enumerate(sorted(cells)):z.append('<polygon points="'+' '.join(f'{x:.17g},{y:.17g}' for x,y in cells[v])+f'" fill="hsl({(i*137)%360} 55% 70%)" fill-opacity=".25" stroke="#78909c" stroke-width=".8"/>')
 colors={"ACCEPTED":"#2e7d32","REJECTED_GRAPH_INVALID":"#c62828","REJECTED_TOTAL_DEFICIENCY":"#7b1fa2","REJECTED_DEFICIENT_VERTEX_WORSENED":"#ef6c00"}
 for a in attempts:
  p=a["starting_coordinate_at_attempt"];q=a["frozen_cell_centroid"];c=colors[a["decision"]];z.append(f'<line x1="{p[0]:.17g}" y1="{p[1]:.17g}" x2="{q[0]:.17g}" y2="{q[1]:.17g}" stroke="{c}" stroke-width="2" stroke-dasharray="4 4"/>')
  if a["decision"]!="ACCEPTED":z.append(f'<path d="M {q[0]-5:.17g} {q[1]-5:.17g} L {q[0]+5:.17g} {q[1]+5:.17g} M {q[0]-5:.17g} {q[1]+5:.17g} L {q[0]+5:.17g} {q[1]-5:.17g}" stroke="{c}" stroke-width="2"/>')
 for u,v in edges:z.append(f'<line x1="{co[u][0]:.17g}" y1="{co[u][1]:.17g}" x2="{co[v][0]:.17g}" y2="{co[v][1]:.17g}" stroke="#444" stroke-width="1.5"/>')
 for v in sorted(vertices):
  x,y=co[v];act=v in active;z.append(f'<circle cx="{x:.17g}" cy="{y:.17g}" r="{7 if act else 5}" fill="{"#ff8f00" if act else "#546e7a"}" stroke="white"/><text x="{x+8:.17g}" y="{y-8:.17g}" font-family="sans-serif" font-size="10">{vertices[v]["name"]} {"ACTIVE" if act else "FROZEN"} CR={rat[v]:.3f}</text>')
 z.append(f'<text x="{xmin+10:.17g}" y="{ymin+22:.17g}" font-family="sans-serif" font-size="15" font-weight="bold">Monotonic 80%: green accept; red graph; purple total; orange protected-cell reject</text></svg>');path.write_text('\n'.join(z)+'\n')

def run(graph_path,pass1_path,simple_path,aware_path,outdir):
 g=json.loads(graph_path.read_text());p1=json.loads(pass1_path.read_text());simple=json.loads(simple_path.read_text());prior=json.loads(aware_path.read_text());vertices={v["id"]:v for v in g["vertices"]};edges=[(e["source"],e["target"]) for e in g["edges"]]
 initial={v:tuple(p) for v,p in p1["final_coordinates"].items()};box=tuple(p1["domain"][k] for k in ("xmin","xmax","ymin","ymax"));target=p1["target_area"];threshold=.8*target;start_cells=core.cells_in_fixed_box(initial,box);start_rows=core.records(start_cells,target)
 if start_rows!=p1["after_density"]:raise AssertionError("baseline mismatch")
 active=[r["vertex_id"] for r in start_rows if r["voronoi_area"]<threshold];frozen=[r["vertex_id"] for r in start_rows if r["voronoi_area"]>=threshold];co=dict(initial);current_rows=start_rows;current_cells=start_cells;current_def=deficiencies(current_rows,threshold);attempts=[]
 for row in start_rows:
  v=row["vertex_id"]
  if v not in active:continue
  start=co[v];cent=core.centroid(start_cells[v]);before_total=sum(current_def.values());protected=[x for x,d in current_def.items() if d>0];co[v]=cent;diag=core.diagnostics(co,edges);reasons=[]
  if diag["proper_crossings"]:reasons.append({"type":"proper_crossings","details":diag["proper_crossing_list"]})
  if diag["coincident_vertices"]:reasons.append({"type":"coincident_vertices","details":diag["coincident_vertex_pairs"]})
  if diag["vertex_on_unrelated_edge_events"]:reasons.append({"type":"vertex_on_unrelated_edge","details":diag["vertex_on_unrelated_edge_list"]})
  graph_valid=not reasons;proposed_total=None;change=None;total_ok=False;monotonic_ok=False;worsened=[]
  if graph_valid:
   proposed_cells=core.cells_in_fixed_box(co,box);proposed_rows=core.records(proposed_cells,target);proposed_def=deficiencies(proposed_rows,threshold);proposed_total=sum(proposed_def.values());change=proposed_total-before_total;total_tol=1e-10*max(1.,before_total);total_ok=proposed_total<before_total-total_tol;area_before={r["vertex_id"]:r["voronoi_area"] for r in current_rows};area_after={r["vertex_id"]:r["voronoi_area"] for r in proposed_rows};area_tol=1e-10*max(1.,threshold)
   for x in protected:
    if proposed_def[x]>current_def[x]+area_tol:worsened.append({"vertex_id":x,"area_before":area_before[x],"area_after_proposed":area_after[x],"deficiency_before":current_def[x],"deficiency_after_proposed":proposed_def[x],"worsening_amount":proposed_def[x]-current_def[x]})
   monotonic_ok=not worsened
  accepted=graph_valid and total_ok and monotonic_ok
  if accepted:current_cells=proposed_cells;current_rows=proposed_rows;current_def=proposed_def;decision="ACCEPTED"
  else:
   co[v]=start
   decision="REJECTED_GRAPH_INVALID" if not graph_valid else "REJECTED_TOTAL_DEFICIENCY" if not total_ok else "REJECTED_DEFICIENT_VERTEX_WORSENED"
  attempts.append({"attempt_order":len(attempts)+1,"vertex_id":v,"frozen_starting_area":row["voronoi_area"],"frozen_cell_centroid":list(cent),"starting_coordinate_at_attempt":list(start),"proposed_displacement_distance":math.dist(start,cent),"current_total_deficiency_before":before_total,"currently_deficient_count":len(protected),"currently_deficient_vertices":protected,"graph_valid":graph_valid,"graph_invalid_reasons":reasons,"proposed_total_deficiency":proposed_total,"total_deficiency_change":change,"total_deficiency_decreased":total_ok,"all_currently_deficient_non_worsening":monotonic_ok,"worsened_deficient_vertices":worsened,"decision":decision,"complete_edge_regeneration_count":44})
 final_cells=core.cells_in_fixed_box(co,box);final_rows=core.records(final_cells,target);final_def=sum(deficiencies(final_rows,threshold).values());diag=core.diagnostics(co,edges);err=sum(r["voronoi_area"] for r in final_rows)-p1["domain"]["area"]
 if any(diag[k] for k in ("proper_crossings","coincident_vertices","vertex_on_unrelated_edge_events")) or abs(err)>1e-7:raise AssertionError("invalid final")
 by0={r["vertex_id"]:r for r in start_rows};bys={x["vertex_id"]:x for x in simple["per_vertex"]};bya={x["vertex_id"]:x for x in prior["per_vertex"]};byf={r["vertex_id"]:r for r in final_rows};att={a["vertex_id"]:a for a in attempts};per=[]
 for v in sorted(vertices):
  x=att.get(v);per.append({"vertex_id":v,"pass1_area":by0[v]["voronoi_area"],"simple_threshold_area":bys[v]["final_area"],"deficiency_aware_area":bya[v]["threshold_aware_area"],"monotonic_area":byf[v]["voronoi_area"],"final_crowding_ratio":byf[v]["crowding_ratio"],"final_density_rank":byf[v]["rank"],"classification":"ACTIVE" if v in active else "FROZEN","proposed":x is not None,"decision":x["decision"] if x else "NOT_ATTEMPTED"})
 def state(stats,areas):return {**stats,"cells_below_threshold":sum(a<threshold for a in areas),"total_deficiency":sum(max(0.,threshold-a) for a in areas)}
 states={"PASS1_BASELINE":state(p1["after_statistics"],[r["voronoi_area"] for r in start_rows]),"SIMPLE_THRESHOLD80":state(simple["after_statistics"],[x["final_area"] for x in simple["per_vertex"]]),"DEFICIENCY_AWARE":prior["comparison"]["THRESHOLD_AWARE"],"CROWDING_MONOTONIC":state(core.stats(final_rows),[r["voronoi_area"] for r in final_rows])};accepted=[a for a in attempts if a["decision"]=="ACCEPTED"]
 report={"schema":"graph-relax.voronoi-threshold80-monotonic.v1","source_pass1_report":str(pass1_path),"domain":p1["domain"],"target_area":target,"threshold":threshold,"active_vertices":active,"frozen_vertices":frozen,"active_count":len(active),"frozen_count":len(frozen),"attempts":attempts,"accepted_moves":len(accepted),"rejected_graph_invalid":sum(a["decision"]=="REJECTED_GRAPH_INVALID" for a in attempts),"rejected_total_deficiency":sum(a["decision"]=="REJECTED_TOTAL_DEFICIENCY" for a in attempts),"rejected_deficient_vertex_worsened":sum(a["decision"]=="REJECTED_DEFICIENT_VERTEX_WORSENED" for a in attempts),"starting_total_deficiency":sum(deficiencies(start_rows,threshold).values()),"final_total_deficiency":final_def,"total_accepted_displacement":sum(a["proposed_displacement_distance"] for a in accepted),"maximum_accepted_displacement":max(a["proposed_displacement_distance"] for a in accepted) if accepted else 0.,"final_coordinates":{v:list(co[v]) for v in sorted(co)},"final_validation":diag,"four_state_comparison":states,"per_vertex":per,"important_vertices":{x["vertex_id"]:x for x in per if x["vertex_id"] in IMPORTANT},"final_density":final_rows,"partition_area_error":err}
 outdir.mkdir(parents=True,exist_ok=True);render(vertices,edges,co,final_cells,final_rows,attempts,set(active),box,outdir/"iamp_voronoi_threshold80_monotonic.svg");(outdir/"iamp_voronoi_threshold80_monotonic_report.json").write_text(json.dumps(report,indent=2,sort_keys=True)+'\n');return report

def main():
 p=argparse.ArgumentParser();b=Path("output/graph_first");p.add_argument("--graph",type=Path,default=b/"iamp_graph.json");p.add_argument("--pass1",type=Path,default=b/"iamp_voronoi_centroid_move_report.json");p.add_argument("--simple",type=Path,default=b/"iamp_voronoi_threshold80_report.json");p.add_argument("--aware",type=Path,default=b/"iamp_voronoi_threshold80_deficiency_report.json");p.add_argument("--output",type=Path,default=b);a=p.parse_args();r=run(a.graph,a.pass1,a.simple,a.aware,a.output);print("ACTIVE",r["active_vertices"]);print("accepted",r["accepted_moves"],"deficiency",r["starting_total_deficiency"],r["final_total_deficiency"]);[print(x["attempt_order"],x["vertex_id"],x["decision"],x["worsened_deficient_vertices"]) for x in r["attempts"]]
if __name__=="__main__":main()
