#!/usr/bin/env python3
"""Second frozen-cell Voronoi-centroid movement pass on the IAMP graph."""
import argparse,json,math
from pathlib import Path
import placement_geometry as geometry
import graph_first_experiment_voronoi_centroid_move as pass1mod

IMPORTANT=("component:TPB","net:EB","net:+9V","net:-9V","component:U1","net:NODE_B","component:J_PWR","component:J_EB","net:EA","component:J_EA","net:GND","net:R1_EXT","component:J_OUT")

def run(graph_path,pass1_path,outdir):
    g=json.loads(graph_path.read_text()); p1=json.loads(pass1_path.read_text())
    vertices={v["id"]:v for v in g["vertices"]}; edges=[(e["source"],e["target"]) for e in g["edges"]]
    initial={v:tuple(p) for v,p in p1["final_coordinates"].items()}; box=tuple(p1["domain"][k] for k in ("xmin","xmax","ymin","ymax")); target=p1["target_area"]
    d0=pass1mod.diagnostics(initial,edges)
    if (len(vertices),len(edges),d0["proper_crossings"],d0["coincident_vertices"],d0["vertex_on_unrelated_edge_events"])!=(34,44,0,0,0): raise AssertionError("invalid Pass-1 source")
    frozen=pass1mod.cells_in_fixed_box(initial,box); starting=pass1mod.records(frozen,target)
    if starting!=p1["after_density"]: raise AssertionError("Pass-2 starting partition differs from Pass-1 final partition")
    co=dict(initial); attempts=[]
    for row in starting:
        v=row["vertex_id"]; start=co[v]; c=pass1mod.centroid(frozen[v]); co[v]=c; d=pass1mod.diagnostics(co,edges); reasons=[]
        if d["proper_crossings"]: reasons.append({"type":"proper_crossings","details":d["proper_crossing_list"]})
        if d["coincident_vertices"]: reasons.append({"type":"coincident_vertices","details":d["coincident_vertex_pairs"]})
        if d["vertex_on_unrelated_edge_events"]: reasons.append({"type":"vertex_on_unrelated_edge","details":d["vertex_on_unrelated_edge_list"]})
        decision="REJECTED" if reasons else "ACCEPTED"
        if reasons: co[v]=start
        attempts.append({"attempt_order":row["rank"],"vertex_id":v,"pass2_frozen_voronoi_area":row["voronoi_area"],"pass2_starting_crowding_ratio":row["crowding_ratio"],
          "starting_coordinate_at_attempt":list(start),"frozen_cell_centroid":list(c),"proposed_displacement_distance":math.dist(start,c),"decision":decision,
          "rejection_reasons":reasons,"complete_edge_regeneration_count":len(edges)})
    final_diag=pass1mod.diagnostics(co,edges)
    if any(final_diag[k] for k in ("proper_crossings","coincident_vertices","vertex_on_unrelated_edge_events")): raise AssertionError("invalid Pass-2 final graph")
    final_cells=pass1mod.cells_in_fixed_box(co,box); final=pass1mod.records(final_cells,target)
    area_error=sum(x["voronoi_area"] for x in final)-p1["domain"]["area"]
    if abs(area_error)>1e-7: raise AssertionError("Pass-2 partition area mismatch")
    original={x["vertex_id"]:x for x in p1["before_density"]}; one={x["vertex_id"]:x for x in p1["after_density"]}; two={x["vertex_id"]:x for x in final}
    p1attempt={x["vertex_id"]:x for x in p1["attempts"]}; p2attempt={x["vertex_id"]:x for x in attempts}
    three=[]
    for v in sorted(vertices):
        three.append({"vertex_id":v,
          "original":{"rank":original[v]["rank"],"area":original[v]["voronoi_area"],"crowding_ratio":original[v]["crowding_ratio"]},
          "after_pass1":{"rank":one[v]["rank"],"area":one[v]["voronoi_area"],"crowding_ratio":one[v]["crowding_ratio"]},
          "after_pass2":{"rank":two[v]["rank"],"area":two[v]["voronoi_area"],"crowding_ratio":two[v]["crowding_ratio"]},
          "moved_pass1":p1attempt[v]["decision"]=="ACCEPTED","moved_pass2":p2attempt[v]["decision"]=="ACCEPTED",
          "rejected_pass1":p1attempt[v]["decision"]=="REJECTED","rejected_pass2":p2attempt[v]["decision"]=="REJECTED"})
    accepted=[x for x in attempts if x["decision"]=="ACCEPTED"]
    report={"schema":"graph-relax.voronoi-centroid-pass2.v1","source_pass1_report":str(pass1_path),"domain":p1["domain"],"target_area":target,
      "vertex_count":34,"electrical_incidence_count":44,"pass2_initial_validation":d0,"frozen_pass2_partition":True,"pass2_attempts":attempts,
      "pass2_accepted_moves":len(accepted),"pass2_rejected_moves":34-len(accepted),"pass2_total_accepted_displacement":sum(x["proposed_displacement_distance"] for x in accepted),
      "pass2_maximum_individual_displacement":max(x["proposed_displacement_distance"] for x in accepted),"pass2_final_coordinates":{v:list(co[v]) for v in sorted(co)},
      "pass2_final_validation":final_diag,"pass2_starting_density":starting,"pass2_final_density":final,"pass2_final_partition_area_error":area_error,
      "three_state_statistics":{"ORIGINAL":p1["before_statistics"],"AFTER_PASS_1":p1["after_statistics"],"AFTER_PASS_2":pass1mod.stats(final)},
      "three_state_by_vertex":three,"important_vertices":{x["vertex_id"]:x for x in three if x["vertex_id"] in IMPORTANT}}
    outdir.mkdir(parents=True,exist_ok=True)
    pass1mod.render(vertices,edges,initial,co,frozen,final_cells,starting,final,box,attempts,outdir/"iamp_voronoi_centroid_pass2.svg")
    (outdir/"iamp_voronoi_centroid_pass2_report.json").write_text(json.dumps(report,indent=2,sort_keys=True)+'\n'); return report

def main():
    p=argparse.ArgumentParser(); b=Path("output/graph_first"); p.add_argument("--graph",type=Path,default=b/"iamp_graph.json"); p.add_argument("--pass1",type=Path,default=b/"iamp_voronoi_centroid_move_report.json"); p.add_argument("--output",type=Path,default=b); a=p.parse_args(); r=run(a.graph,a.pass1,a.output)
    print("accepted",r["pass2_accepted_moves"],"rejected",r["pass2_rejected_moves"],"total displacement",r["pass2_total_accepted_displacement"])
    for x in r["pass2_attempts"]: print(x["attempt_order"],x["vertex_id"],x["decision"],x["proposed_displacement_distance"],x["rejection_reasons"])
if __name__=="__main__": main()
