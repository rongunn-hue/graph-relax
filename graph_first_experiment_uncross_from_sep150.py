#!/usr/bin/env python3
"""Direct analytic uncrossing from the authoritative 150-separated geometry."""
import argparse,json
from itertools import combinations
from pathlib import Path
import networkx as nx
import placement_geometry as geometry
import graph_first_experiment_rotation_direct as direct
import graph_first_experiment_node_separation_150 as sep

S_MIN=150.0
def audit(co):
 rows=[{"vertices":[a,b],"distance":__import__('math').dist(co[a],co[b])} for a,b in combinations(sorted(co),2)];rows.sort(key=lambda x:(x["distance"],x["vertices"]));return rows

def load(graph_path,source_path):
 gd=json.loads(graph_path.read_text());sd=json.loads(source_path.read_text());g=nx.Graph();g.add_nodes_from((v["id"],{"type":v["type"],"name":v["name"]}) for v in gd["vertices"]);g.add_edges_from((e["source"],e["target"]) for e in gd["edges"]);co={v:tuple(p) for v,p in sd["final_coordinates"].items()}
 if (g.number_of_nodes(),g.number_of_edges(),set(g),set(co))!=(34,44,set(co),set(co)):raise AssertionError("source connectivity mismatch")
 return g,co,sd

def run(graph_path,source_path,outdir):
 graph,start,source=load(graph_path,source_path);pairs=audit(start);viol=[x for x in pairs if x["distance"]<S_MIN]
 if len(pairs)!=561 or viol:raise AssertionError(f"authoritative sep150 input failed literal audit: {viol}")
 start_diag=geometry.graph_geometry_diagnostics(start,graph.edges)
 final,final_diag,moves,states,history,counters=direct.direct_search(graph,start,S_MIN,0.0)
 accepted=[]
 for i,m in enumerate(moves):
  before=audit(states[i]);after=audit(states[i+1]);bad=[x for x in after if x["distance"]<S_MIN]
  if bad:raise AssertionError(f"accepted move violates literal 150 audit: {bad}")
  accepted.append({"repair_number":i+1,"crossing_being_addressed":m["crossing_being_repaired"],"moving_endpoint":m["moving_vertex"],"pivot":m["pivot_vertex"],"angular_displacement_radians":m["signed_rotation_radians"],"angular_displacement_degrees":m["signed_rotation_degrees"],"crossings_before":m["crossing_count_before"],"crossings_after":m["crossing_count_after"],"minimum_pairwise_distance_before":before[0]["distance"],"minimum_pairwise_distance_after":after[0]["distance"],"pair_audit":{"pairs_audited":len(after),"pairs_below_150":len(bad),"passed":not bad}})
 final_pairs=audit(final);final_bad=[x for x in final_pairs if x["distance"]<S_MIN]
 if final_bad:raise AssertionError("final literal 150 audit failed")
 # One limiting analytic relaxation, not a parameter sweep: remove spacing
 # exclusions entirely and retain the unchanged direct candidate mechanism.
 relaxed_counters={"candidate_evaluations":0,"positions_per_choice":[],"solutions_per_crossing":[],"all_solution_outcomes":[],"orientation_audits":[]}
 relaxed_choice=None
 if final_diag["proper_unrelated_edge_crossing_count"]:
  relaxed_choice=direct.best_direct_move(graph,final,final_diag,relaxed_counters,0.0,0.0)
 relaxed_reducing=[x for x in relaxed_counters["all_solution_outcomes"] if x.get("valid") and x.get("resulting_total_crossings",10**9)<final_diag["proper_unrelated_edge_crossing_count"]]
 if relaxed_choice is not None or relaxed_reducing:
  raise AssertionError("spacing-relaxed continuation unexpectedly exists")
 blocked=[]
 for crossing in final_diag["proper_unrelated_edge_crossings"]:
  matching=[a for a in counters["orientation_audits"] if a["crossing"]==crossing]
  audit_record=matching[-1] if matching else None;solutions=[]
  for s in (audit_record or {"solutions":[]})["solutions"]:
   spacing=s.get("spacing_intervals",[]);solutions.append({"moving_endpoint":s.get("moving"),"pivot":s.get("pivot"),"current_radius":s.get("radius"),"crossing_forbidden_interval":s.get("target_interval"),"spacing_exclusion_interval_count":sum("interval" in x for x in spacing),"spacing_interval_vertices":[x["vertex"] for x in spacing if "interval" in x],"full_circle_infeasible":s.get("reason") in {"ENTIRE_MOTION_CIRCLE_FORBIDDEN_BY_SPACING","FORBIDDEN_INTERVAL_UNION_COVERS_CIRCLE"},"full_circle_blocking_vertices":[x["vertex"] for x in spacing if x["status"]=="ENTIRE"],"candidate_crossing_count":s.get("resulting_total_crossings"),"candidate_minimum_node_distance":s.get("resulting_minimum_node_distance"),"rejection_reason":s.get("reason") or s.get("selection_status")})
  blocked.append({"crossing":crossing,"orientations":solutions})
 attempted=len(counters["all_solution_outcomes"]);infeasible=sum(not x.get("valid",False) for x in counters["all_solution_outcomes"])
 report={"schema":"graph-relax.uncross-from-sep150.v1","authoritative_input":str(source_path),"node_separation_rerun":False,"S_MIN":S_MIN,"starting_audit":{"pairs_audited":561,"minimum_pairwise_distance":pairs[0]["distance"],"pairs_below_150":0,"violating_pairs":[],"proper_crossing_count":start_diag["proper_unrelated_edge_crossing_count"]},"accepted_moves":accepted,"crossing_history":history,"blocked_crossings":blocked,"orientation_audits":counters["orientation_audits"],"attempted_orientation_count":attempted,"infeasible_orientation_count":infeasible,"maximum_positions_per_orientation":max((x["positions_evaluated"] for x in counters["positions_per_choice"]),default=0),"maximum_orientations_per_crossing":max((x["calculated"] for x in counters["solutions_per_crossing"]),default=0),"final_audit":{"pairs_audited":561,"minimum_pairwise_distance":final_pairs[0]["distance"],"pairs_below_150":0,"violating_pairs":[],"proper_crossing_count":final_diag["proper_unrelated_edge_crossing_count"],"coincident_vertex_count":final_diag["coincident_vertex_pair_count"],"vertex_on_unrelated_edge_count":final_diag["vertex_on_unrelated_edge_interior_count"],"vertex_count":34,"edge_count":44,"connectivity_unchanged":True},"accepted_move_count":len(moves),"zero_crossings_achieved":final_diag["proper_unrelated_edge_crossing_count"]==0,"full_150_spacing_preserved":not final_bad and all(x["pair_audit"]["passed"] for x in accepted),"final_coordinates":{v:list(final[v]) for v in sorted(final)}}
 report["schema"]="graph-relax.uncross-from-sep150.v2"
 report["spacing_policy"]={"target":S_MIN,"zero_crossings_required":True,"target_attempted_first":True,"fallback":{"kind":"analytic limiting relaxation","minimum_spacing_floor":0.0,"parameter_sweep":False}}
 report["relaxed_spacing_fallback_audits"]=relaxed_counters["orientation_audits"]
 report["relaxed_reducing_candidate_count"]=len(relaxed_reducing)
 report["spacing_relaxation_outcome"]="NO_CROSSING_REDUCING_DIRECT_CANDIDATE_EVEN_WITH_ALL_SPACING_EXCLUSIONS_REMOVED" if final_diag["proper_unrelated_edge_crossing_count"] else "NOT_NEEDED"
 report["final_audit"]["amount_below_150_target"]=max(0.0,S_MIN-final_pairs[0]["distance"])
 report["final_audit"]["limiting_node_pair"]=final_pairs[0]["vertices"]
 outdir.mkdir(parents=True,exist_ok=True);sep.render(graph,start,final,outdir/"iamp_uncross_from_sep150.svg","A. AUTHORITATIVE 150-SEPARATED START","B. FINAL — TARGET-SPACING ANALYTIC UNCROSSING");(outdir/"iamp_uncross_from_sep150_report.json").write_text(json.dumps(report,indent=2,sort_keys=True)+'\n');return report

def main():
 p=argparse.ArgumentParser();b=Path("output/graph_first");p.add_argument("--graph",type=Path,default=b/"iamp_graph.json");p.add_argument("--source",type=Path,default=b/"iamp_node_separation_150_report.json");p.add_argument("--output",type=Path,default=b);a=p.parse_args();r=run(a.graph,a.source,a.output);print(json.dumps({k:v for k,v in r.items() if k not in {"orientation_audits","final_coordinates"}},indent=2))
if __name__=="__main__":main()
