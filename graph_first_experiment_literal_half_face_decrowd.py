#!/usr/bin/env python3
"""Run the literal D/2 alternating IAMP face de-crowding experiment."""
import argparse,json
from pathlib import Path
import graph_first_experiment_outward_face_wave as core

def run(graph_path,base_path,outdir):
    report=core.run(graph_path,base_path,outdir)
    report["schema"]="graph-relax.literal-half-face-decrowd.v1"
    report["movement_definition"]={"D":"minimum actual vertical separation between the non-shared boundaries of the current and dynamically selected next adjacent face","TARGET":"D / 2","readability_clearance_used_for_target":False,"readability_clearance_role":"post-move reporting only"}
    moves=[x for x in report["trace"] if x.get("status") in ("MOVED_UP","MOVED_DOWN")]
    report["literal_half_assertions"]={
      "every_target_equals_D_over_2":all(abs(x["half_gap"]-x["full_gap"]/2)<=1e-9 for x in moves),
      "every_valid_half_was_accepted_exactly":all((not x["exact_half_hard_valid"]) or abs(x["actual_move"]-x["half_gap"])<=1e-9 for x in moves),
      "every_reduction_has_hard_failure":all(x["actual_move"]+1e-9>=x["half_gap"] or bool(x["exact_half_hard_failures"]) for x in moves),
      "d6_never_used_for_target":all(not x["target_uses_readability_clearance"] for x in moves),
      "alternating_execution":[x["side"] for x in moves],
      "second_top_contains_C1":report["second_top_reaches_c1_face"],
      "second_bottom_contains_C2":report["second_bottom_reaches_c2_face"]}
    if not all(v for k,v in report["literal_half_assertions"].items() if k!="alternating_execution"):raise AssertionError(report["literal_half_assertions"])
    outdir.mkdir(parents=True,exist_ok=True)
    source=outdir/"iamp_outward_face_wave.svg";target=outdir/"iamp_literal_half_face_decrowd.svg"
    target.write_text(source.read_text())
    (outdir/"iamp_literal_half_face_decrowd_report.json").write_text(json.dumps(report,indent=2,sort_keys=True)+"\n")
    return report

def main():
    p=argparse.ArgumentParser();b=Path("output/graph_first");p.add_argument("--graph",type=Path,default=b/"iamp_graph.json");p.add_argument("--base",type=Path,default=b/"iamp_corridor_separation_report.json");p.add_argument("--output",type=Path,default=b);a=p.parse_args();r=run(a.graph,a.base,a.output)
    for x in r["trace"]:
      print(x["side"],"current =",x["face_id"],"next =",x.get("next_face_id"),"movable =",x.get("moved_vertices"),"D =",x.get("full_gap"),"HALF =",x.get("half_gap"),"HALF_VALID =",x.get("exact_half_hard_valid"),"ACTUAL_MOVE =",x.get("actual_move"),"reason =",x.get("reduction_reason"),"crossings =",x.get("metrics_after",{}).get("proper_crossings"),"violations =",x.get("metrics_after",{}).get("clearance_violations"))
if __name__=="__main__":main()
