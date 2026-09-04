#!/usr/bin/env python3
"""Third frozen-cell Voronoi-centroid movement pass on the IAMP graph."""
import argparse,json,math
from pathlib import Path
import graph_first_experiment_voronoi_centroid_move as core

TRACK=("net:EB","component:U1","net:-9V","component:J_PWR","component:C1","component:TPB","net:NODE_B","net:+9V")

def run(graph_path,pass2_path,outdir):
 g=json.loads(graph_path.read_text()); p2=json.loads(pass2_path.read_text()); vertices={v["id"]:v for v in g["vertices"]}; edges=[(e["source"],e["target"]) for e in g["edges"]]
 initial={v:tuple(p) for v,p in p2["pass2_final_coordinates"].items()}; box=tuple(p2["domain"][k] for k in ("xmin","xmax","ymin","ymax")); target=p2["target_area"]
 d0=core.diagnostics(initial,edges)
 if (len(vertices),len(edges),d0["proper_crossings"],d0["coincident_vertices"],d0["vertex_on_unrelated_edge_events"])!=(34,44,0,0,0):raise AssertionError("invalid Pass-2 source")
 frozen=core.cells_in_fixed_box(initial,box); starting=core.records(frozen,target)
 if starting!=p2["pass2_final_density"]:raise AssertionError("Pass-3 starting partition differs from Pass-2 final partition")
 co=dict(initial);attempts=[]
 for row in starting:
  v=row["vertex_id"];start=co[v];c=core.centroid(frozen[v]);co[v]=c;d=core.diagnostics(co,edges);reasons=[]
  if d["proper_crossings"]:reasons.append({"type":"proper_crossings","details":d["proper_crossing_list"]})
  if d["coincident_vertices"]:reasons.append({"type":"coincident_vertices","details":d["coincident_vertex_pairs"]})
  if d["vertex_on_unrelated_edge_events"]:reasons.append({"type":"vertex_on_unrelated_edge","details":d["vertex_on_unrelated_edge_list"]})
  decision="REJECTED" if reasons else "ACCEPTED"
  if reasons:co[v]=start
  attempts.append({"attempt_order":row["rank"],"vertex_id":v,"pass3_frozen_voronoi_area":row["voronoi_area"],"pass3_starting_crowding_ratio":row["crowding_ratio"],"starting_coordinate_at_attempt":list(start),"frozen_cell_centroid":list(c),"proposed_displacement_distance":math.dist(start,c),"decision":decision,"rejection_reasons":reasons,"complete_edge_regeneration_count":44})
 final_diag=core.diagnostics(co,edges)
 if any(final_diag[k] for k in ("proper_crossings","coincident_vertices","vertex_on_unrelated_edge_events")):raise AssertionError("invalid Pass-3 final graph")
 final_cells=core.cells_in_fixed_box(co,box);final=core.records(final_cells,target);area_error=sum(x["voronoi_area"] for x in final)-p2["domain"]["area"]
 if abs(area_error)>1e-7:raise AssertionError("Pass-3 partition area mismatch")
 old={x["vertex_id"]:x for x in p2["three_state_by_vertex"]};three={x["vertex_id"]:x for x in p2["pass2_attempts"]};now={x["vertex_id"]:x for x in final};a3={x["vertex_id"]:x for x in attempts};four=[]
 for v in sorted(vertices):
  x=old[v];four.append({"vertex_id":v,"original":x["original"],"after_pass1":x["after_pass1"],"after_pass2":x["after_pass2"],"after_pass3":{"rank":now[v]["rank"],"area":now[v]["voronoi_area"],"crowding_ratio":now[v]["crowding_ratio"]},"pass1_decision":"REJECTED" if x["rejected_pass1"] else "ACCEPTED","pass2_decision":three[v]["decision"],"pass3_decision":a3[v]["decision"]})
 accepted=[x for x in attempts if x["decision"]=="ACCEPTED"];stats={**p2["three_state_statistics"],"AFTER_PASS_3":core.stats(final)};sd=[stats[k]["population_standard_deviation"] for k in ("ORIGINAL","AFTER_PASS_1","AFTER_PASS_2","AFTER_PASS_3")]
 report={"schema":"graph-relax.voronoi-centroid-pass3.v1","source_pass2_report":str(pass2_path),"domain":p2["domain"],"target_area":target,"vertex_count":34,"electrical_incidence_count":44,"pass3_initial_validation":d0,"frozen_pass3_partition":True,"pass3_attempts":attempts,"pass3_accepted_moves":len(accepted),"pass3_rejected_moves":34-len(accepted),"pass3_total_accepted_displacement":sum(x["proposed_displacement_distance"] for x in accepted),"pass3_maximum_accepted_displacement":max(x["proposed_displacement_distance"] for x in accepted),"pass3_final_coordinates":{v:list(co[v]) for v in sorted(co)},"pass3_final_validation":final_diag,"pass3_starting_density":starting,"pass3_final_density":final,"pass3_partition_area_error":area_error,"four_state_statistics":stats,"standard_deviation_percent_changes":{"ORIGINAL_TO_PASS_1":100*(sd[1]/sd[0]-1),"PASS_1_TO_PASS_2":100*(sd[2]/sd[1]-1),"PASS_2_TO_PASS_3":100*(sd[3]/sd[2]-1)},"accepted_displacement_by_pass":{"PASS_1":p2["three_state_statistics"] and json.loads(Path(p2["source_pass1_report"]).read_text())["total_vertex_displacement"],"PASS_2":p2["pass2_total_accepted_displacement"],"PASS_3":sum(x["proposed_displacement_distance"] for x in accepted)},"four_state_by_vertex":four,"tracked_vertices":{x["vertex_id"]:x for x in four if x["vertex_id"] in TRACK}}
 outdir.mkdir(parents=True,exist_ok=True);core.render(vertices,edges,initial,co,frozen,final_cells,starting,final,box,attempts,outdir/"iamp_voronoi_centroid_pass3.svg");(outdir/"iamp_voronoi_centroid_pass3_report.json").write_text(json.dumps(report,indent=2,sort_keys=True)+'\n');return report

def main():
 p=argparse.ArgumentParser();b=Path("output/graph_first");p.add_argument("--graph",type=Path,default=b/"iamp_graph.json");p.add_argument("--pass2",type=Path,default=b/"iamp_voronoi_centroid_pass2_report.json");p.add_argument("--output",type=Path,default=b);a=p.parse_args();r=run(a.graph,a.pass2,a.output);print("accepted",r["pass3_accepted_moves"],"rejected",r["pass3_rejected_moves"],"displacement",r["pass3_total_accepted_displacement"]);[print(x["attempt_order"],x["vertex_id"],x["decision"],x["proposed_displacement_distance"],x["rejection_reasons"]) for x in r["pass3_attempts"]]
if __name__=="__main__":main()
