#!/usr/bin/env python3
"""Retry the direct analytic uncrosser from its current saved final geometry."""
import argparse, json, math
from itertools import combinations
from pathlib import Path
import networkx as nx
import placement_geometry as geometry
import graph_first_experiment_rotation_direct as direct
import graph_first_experiment_node_separation_150 as sep

TARGET=150.0

def pair_audit(coordinates):
 rows=[{"vertices":[a,b],"distance":math.dist(coordinates[a],coordinates[b])}
       for a,b in combinations(sorted(coordinates),2)]
 return sorted(rows,key=lambda x:(x["distance"],x["vertices"]))

def counters():
 return {"candidate_evaluations":0,"positions_per_choice":[],"solutions_per_crossing":[],
         "all_solution_outcomes":[],"orientation_audits":[]}

def run(graph_path,source_path,outdir):
 gd=json.loads(graph_path.read_text()); source=json.loads(source_path.read_text())
 graph=nx.Graph(); graph.add_nodes_from((v["id"],{"type":v["type"],"name":v["name"]}) for v in gd["vertices"])
 graph.add_edges_from((e["source"],e["target"]) for e in gd["edges"])
 start={v:tuple(p) for v,p in source["final_coordinates"].items()}
 if graph.number_of_nodes()!=34 or graph.number_of_edges()!=44 or set(start)!=set(graph):
  raise AssertionError("authoritative retry input mismatch")
 positions=dict(start); accepted=[]; iterations=[]
 while True:
  diagnostic=geometry.graph_geometry_diagnostics(positions,graph.edges)
  if not diagnostic["proper_unrelated_edge_crossing_count"]: break
  target_counters=counters()
  target_choice=direct.best_direct_move(graph,positions,diagnostic,target_counters,TARGET,0.0)
  relaxed_counters=None; choice=target_choice; mode="TARGET_150"
  if choice is None:
   relaxed_counters=counters()
   choice=direct.best_direct_move(graph,positions,diagnostic,relaxed_counters,0.0,0.0)
   mode="LIMITING_RELAXATION_TO_ZERO_FLOOR"
  record={"iteration":len(iterations)+1,"crossings_before":diagnostic["proper_unrelated_edge_crossing_count"],
          "minimum_node_distance_before":diagnostic["minimum_node_distance"],
          "target_150_orientation_audits":target_counters["orientation_audits"],
          "relaxed_orientation_audits":None if relaxed_counters is None else relaxed_counters["orientation_audits"],
          "selection_mode":mode,"accepted":choice is not None}
  iterations.append(record)
  if choice is None: break
  candidate=choice["candidate"]
  after=geometry.graph_geometry_diagnostics(candidate,graph.edges)
  if after["proper_unrelated_edge_crossing_count"]>=diagnostic["proper_unrelated_edge_crossing_count"]:
   raise AssertionError("retry accepted non-reducing candidate")
  rebuilt=geometry.rebuild_straight_edge_segments(candidate,graph.edges)
  reconnect=geometry.validate_rebuilt_edge_segments(candidate,graph.edges,rebuilt)
  if not reconnect["valid"] or after["coincident_vertex_pair_count"] or after["vertex_on_unrelated_edge_interior_count"]:
   raise AssertionError("retry accepted hard-invalid candidate")
  pairs=pair_audit(candidate)
  accepted.append({"repair_number":len(accepted)+1,"moving_endpoint":choice["moving"],"pivot":choice["pivot"],
                   "crossing_being_addressed":{"edge_a":list(choice["target"][0]),"edge_b":list(choice["target"][1])},
                   "angular_displacement_radians":choice["signed_delta_radians"],
                   "crossings_before":diagnostic["proper_unrelated_edge_crossing_count"],
                   "crossings_after":after["proper_unrelated_edge_crossing_count"],
                   "minimum_node_distance_before":diagnostic["minimum_node_distance"],
                   "minimum_node_distance_after":pairs[0]["distance"],
                   "preserved_150":pairs[0]["distance"]>=TARGET,"limiting_pair":pairs[0]["vertices"]})
  positions=candidate
 final=geometry.graph_geometry_diagnostics(positions,graph.edges); pairs=pair_audit(positions)
 last=iterations[-1] if iterations else None
 blocked=[]
 if last and not last["accepted"]:
  relaxed=last["relaxed_orientation_audits"] or []
  for crossing in relaxed:
   orientations=[]
   for s in crossing["solutions"]:
    after=s.get("resulting_total_crossings")
    if not s.get("valid"): classification="GEOMETRICALLY_INVALID"
    elif after is not None and after<last["crossings_before"]: classification="CROSSING_REDUCING"
    elif after==last["crossings_before"]: classification="NO_NET_REDUCTION"
    else: classification="CREATES_ADDITIONAL_CROSSINGS"
    orientations.append({"moving_endpoint":s.get("moving"),"pivot":s.get("pivot"),"pivot_radius":s.get("radius"),
      "current_angle":s.get("theta_old"),"crossing_reducing_position_exists":classification=="CROSSING_REDUCING",
      "maximum_minimum_node_spacing_achievable_while_removing_crossing":s.get("resulting_minimum_node_distance") if classification=="CROSSING_REDUCING" else None,
      "resulting_crossing_count":after,"candidate_minimum_node_distance":s.get("resulting_minimum_node_distance"),
      "classification":classification,"reason":s.get("reason") or s.get("selection_status")})
   blocked.append({"crossing":crossing["crossing"],"orientations":orientations})
 report={"schema":"graph-relax.uncross-retry.v1","authoritative_input":str(source_path),
  "starting_coordinates_are_source_final_coordinates":start=={v:tuple(p) for v,p in source["final_coordinates"].items()},
  "spacing_target":TARGET,"spacing_is_hard_minimum":False,"spacing_sweep_used":False,
  "starting":{"crossings":geometry.graph_geometry_diagnostics(start,graph.edges)["proper_unrelated_edge_crossing_count"],"minimum_node_distance":pair_audit(start)[0]["distance"]},
  "iterations":iterations,"accepted_moves":accepted,"blocked_crossings":blocked,
  "final":{"crossings":final["proper_unrelated_edge_crossing_count"],"minimum_node_distance":pairs[0]["distance"],
   "amount_below_150_target":max(0.0,TARGET-pairs[0]["distance"]),"limiting_node_pairs":[x for x in pairs if abs(x["distance"]-pairs[0]["distance"])<=1e-9],
   "coincident_vertex_count":final["coincident_vertex_pair_count"],"vertex_on_unrelated_edge_count":final["vertex_on_unrelated_edge_interior_count"],
   "vertex_count":34,"edge_count":44,"connectivity_unchanged":True,"zero_crossings_achieved":final["proper_unrelated_edge_crossing_count"]==0},
  "final_coordinates":{v:list(positions[v]) for v in sorted(positions)}}
 outdir.mkdir(parents=True,exist_ok=True)
 sep.render(graph,start,positions,outdir/"iamp_uncross_retry_from_current.svg","A. CURRENT UNCROSSER FINAL — RETRY START","B. AFTER DIRECT ANALYTIC RETRY")
 (outdir/"iamp_uncross_retry_from_current_report.json").write_text(json.dumps(report,indent=2,sort_keys=True)+'\n')
 return report

def main():
 p=argparse.ArgumentParser(); b=Path("output/graph_first")
 p.add_argument("--graph",type=Path,default=b/"iamp_graph.json")
 p.add_argument("--source",type=Path,default=b/"iamp_uncross_from_sep150_report.json")
 p.add_argument("--output",type=Path,default=b); a=p.parse_args()
 r=run(a.graph,a.source,a.output)
 print(json.dumps({"starting":r["starting"],"accepted_moves":r["accepted_moves"],"final":r["final"],"blocked_crossings":r["blocked_crossings"]},indent=2))
if __name__=="__main__": main()
