#!/usr/bin/env python3
"""Metric-derived aggregate pass with locks on explicitly satisfied targets."""
import argparse,json,math
from pathlib import Path
import graph_first_experiment_common_metric as cm
import graph_first_experiment_common_metric_guard as guard
import graph_first_experiment_metric_derived as derived
import graph_first_experiment_planar_region_inventory as inventory
import graph_first_experiment_polygon_area_expansion as polygon
from render_iamp_common_metric_direct_clean import clean_svg

TOL=1e-7
def signature(co):return tuple((v,co[v][0],co[v][1]) for v in sorted(co))
def canonical(k):return listify(k)
def listify(x):return [listify(v) for v in x] if isinstance(x,tuple) else x

class ProtectionGate:
 def __init__(self):self.locks={};self.creation=[];self.pending={};self.stage='polygon';self.rejections=0
 def sync(self,co,faces):
  sig=signature(co);p=self.pending.get(sig)
  if p:
   k=p['key'];old=self.locks.get(k)
   if old is None:
    item={**p,'creation_index':len(self.creation)+1,'minimum_d_after_protection':p['d_at_lock']};self.locks[k]=item;self.creation.append(item)
   elif p['D_lock']>old['D_lock']+TOL:old['D_lock']=p['D_lock'];old['d_at_lock']=p['d_at_lock'];old['updated_by_later_target']=True
  for k,x in self.locks.items():x['minimum_d_after_protection']=min(x['minimum_d_after_protection'],guard.rel_value(k,co,faces))
 def __call__(self,co,candidate,edges,faces,movers,target_key,target_before,target_after,exclude=()):
  self.sync(co,faces);accepted,decision,result=derived.aggregate_gate(co,candidate,edges,faces,movers,target_key,target_before,target_after,exclude)
  result.update({'protected_count_before':len(self.locks),'protected_relation_violations':[],'candidate_created_protection':False,'candidate_protection_identity':canonical(target_key) if target_key else None})
  if not accepted:return accepted,decision,result
  for k,x in self.locks.items():
   d=guard.rel_value(k,candidate,faces)
   if d<x['D_lock']-TOL:result['protected_relation_violations'].append({'canonical_relation_identity':canonical(k),'relation_type':k[0],'d_candidate':d,'D_lock':x['D_lock'],'shortfall':x['D_lock']-d,'created_by_event':x.get('event_created')})
  if result['protected_relation_violations']:
   self.rejections+=1;return False,'REJECT_PROTECTED_RELATION_REGRESSION',result
  if target_key is not None and target_before<1-TOL and target_after>=1-TOL:
   D=guard.affected_snapshot(co,edges,faces,movers).get(target_key,{}).get('D')
   if D is None:
    # The target is always in the affected set; retained only as a defensive assertion.
    raise AssertionError(('target absent from affected set',target_key))
   d=guard.rel_value(target_key,candidate,faces);old=self.locks.get(target_key);lock=max(D,old['D_lock'] if old else D)
   p={'key':target_key,'canonical_relation_identity':canonical(target_key),'relation_type':target_key[0],'D_lock':lock,'d_at_lock':d,'r_before':target_before,'r_after':target_after,'stage_created':self.stage,'event_created':None,'remained_satisfied_through_final':None}
   self.pending[signature(candidate)]=p;result['candidate_created_protection']=True;result['candidate_D_lock']=lock
  return True,'ACCEPTED',result

def summary(ops):
 c=[]
 for o in ops:c+=o.get('candidate_results',[o])
 return {'attempted':len(c),'accepted':sum(x.get('decision')=='ACCEPTED' for x in c),'rejected_target_not_improved':sum(x.get('decision')=='REJECTED_TARGET_NOT_IMPROVED' for x in c),'rejected_aggregate_crowding':sum(x.get('decision')=='REJECTED_AGGREGATE_CROWDING' for x in c),'rejected_graph_validity':sum(x.get('decision')=='REJECTED_GRAPH_INVALID' for x in c),'rejected_protected_relation_regression':sum(x.get('decision')=='REJECT_PROTECTED_RELATION_REGRESSION' for x in c),'skipped_no_longer_crowded':sum(x.get('decision')=='SKIPPED_NO_LONGER_CROWDED' for x in c)}

def attach_creation_events(manager,nops,eops):
 accepted=[]
 for o in nops:
  if o.get('decision')=='ACCEPTED' and o.get('gate',{}).get('candidate_created_protection'):accepted.append(('node',o.get('event_id'),o))
 for o in eops:
  if o.get('decision')=='ACCEPTED':
   c=next((x for x in o.get('candidate_results',[]) if x.get('mover')==o.get('selected_mover') and x.get('decision')=='ACCEPTED'),None)
   if c and c.get('gate',{}).get('candidate_created_protection'):accepted.append(('edge',o.get('event_id'),c))
 for item in manager.creation:
  match=next(((stage,eid,c) for stage,eid,c in accepted if c.get('gate',{}).get('candidate_protection_identity')==item['canonical_relation_identity']),None)
  if match:
   stage,eid,c=match;item['event_created']=eid;item['stage_created']=stage;c['created_protection']=True

def audits(nops,eops):
 tpb=[];jeb=[]
 for stage,ops in [('node',nops),('edge',eops)]:
  for o in ops:
   text=str(o);relations=o.get('candidate_results',[o])
   for c in relations:
    g=c.get('gate',o.get('gate',{}));row={'event_id':o.get('event_id'),'stage':stage,'relation_type':o.get('event_type'),'candidate_mover':c.get('mover',o.get('selected_mover')),'d_before':o.get('d_before'),'D':o.get('D'),'r_before':o.get('r_before'),'d_after':c.get('d_after',o.get('d_after_candidate')),'r_after':c.get('r_after',o.get('r_after_candidate')),'P_before':g.get('P_before'),'P_after':g.get('P_after'),'decision':c.get('decision',o.get('decision')),'protected_before_candidate':g.get('protected_count_before',0)>0,'created_protection':c.get('created_protection',False),'D_lock_after_event':g.get('candidate_D_lock'),'protected_regression_violations':g.get('protected_relation_violations',[])}
    if 'component:TPB' in text and 'net:NODE_B' in text:tpb.append(row)
    if row['candidate_mover'] in ('component:J_EB','component:C2') or 'component:J_EB' in text or 'component:C2' in text:jeb.append(row)
 return tpb,jeb

def run(graph_path,direct_path,outdir):
 gd=json.loads(graph_path.read_text());dr=json.loads(direct_path.read_text());V={x['id']:x for x in gd['vertices']};E=[(x['source'],x['target']) for x in gd['edges']];pairs=guard.edgeset(E);start={v:tuple(p) for v,p in dr['final_coordinates'].items()};ok,sd=cm.hard_valid(start,E)
 if not ok or (len(V),len(E))!=(34,44):raise AssertionError('invalid authoritative start')
 fl,darts=inventory.enumerate_faces(V,E,start);inventory.add_classification_and_adjacency(fl,E,darts);faces={f['face_id']:f for f in fl};initial=cm.corrected_map(fl,start,E);manager=ProtectionGate();old=guard.gate;guard.gate=manager
 try:
  co,pstage=derived.polygon_pass(dict(start),E,faces,initial['polygon'][0]);manager.sync(co,faces);polygon_co=dict(co);after_polygon=cm.corrected_map(fl,co,E);manager.stage='node';co,nops=guard.node_pass(co,E,faces,after_polygon);manager.sync(co,faces);node_co=dict(co);after_nodes=cm.corrected_map(fl,co,E);manager.stage='edge';co,eops=guard.edge_pass(co,E,faces,after_nodes);manager.sync(co,faces)
 finally:guard.gate=old
 attach_creation_events(manager,nops,eops);final=cm.corrected_map(fl,co,E);manager.sync(co,faces);ok,fd=cm.hard_valid(co,E);assert ok and pairs==guard.edgeset(E)
 for k,x in manager.locks.items():x['final_d']=guard.rel_value(k,co,faces);x['remained_satisfied_through_final']=x['final_d']>=x['D_lock']-TOL
 stages={'start':start,'polygon':polygon_co,'nodes':node_co,'final':co};outdir.mkdir(parents=True,exist_ok=True);svgs=[]
 for name,c in stages.items():p=outdir/f'iamp_metric_protected_{name}.svg';clean_svg(V,E,c,f'IAMP METRIC-PROTECTED — {name.upper()}',p);svgs.append(p.name)
 tpb,jeb=audits(nops,eops);report={'schema':'graph-relax.metric-protected.v1','source_coordinates':f'{direct_path}#final_coordinates','starting_validation':{'vertices':34,'edges':44,'endpoint_pairs_identical':pairs==guard.edgeset(E),'proper_crossings':sd['proper_unrelated_edge_crossing_count'],'coincidences':sd['coincident_vertex_pair_count'],'vertex_on_unrelated_edge':sd['vertex_on_unrelated_edge_interior_count']},'initial_map':initial,'start_detector_penalty':derived.global_penalty(initial),'polygon_stage':{**pstage,'summary':summary(pstage['operations']),'map_after':after_polygon},'node_stage':{'operations':nops,'summary':summary(nops),'map_after':after_nodes},'edge_stage':{'operations':eops,'summary':summary(eops)},'TPB_audit':tpb,'J_EB_C2_distances':{k:math.dist(c['component:J_EB'],c['component:C2']) for k,c in stages.items()},'J_EB_C2_candidate_audit':jeb,'protection_report':{'creation_order':manager.creation,'total_created':len(manager.creation),'protection_regression_rejections':manager.rejections,'all_satisfied_at_final':all(x['remained_satisfied_through_final'] for x in manager.creation)},'final_map':final,'final_detector_penalty':derived.global_penalty(final),'final_validation':{'vertices':34,'edges':44,'endpoint_pairs_identical':pairs==guard.edgeset(E),'proper_crossings':fd['proper_unrelated_edge_crossing_count'],'coincidences':fd['coincident_vertex_pair_count'],'vertex_on_unrelated_edge':fd['vertex_on_unrelated_edge_interior_count'],'all_protected_relations_satisfied':all(x['remained_satisfied_through_final'] for x in manager.creation)},'svg_files':svgs,'stage_coordinates':{k:{v:list(p) for v,p in sorted(c.items())} for k,c in stages.items()}}
 (outdir/'iamp_metric_protected_report.json').write_text(json.dumps(report,indent=2,sort_keys=True)+'\n');return report
def main():
 p=argparse.ArgumentParser();b=Path('output/graph_first');p.add_argument('--graph',type=Path,default=b/'iamp_graph.json');p.add_argument('--direct',type=Path,default=b/'iamp_nearness_rotation_direct_report.json');p.add_argument('--output',type=Path,default=b);a=p.parse_args();r=run(a.graph,a.direct,a.output);print(json.dumps({'start_penalty':r['start_detector_penalty'],'final_penalty':r['final_detector_penalty'],'node':r['node_stage']['summary'],'edge':r['edge_stage']['summary'],'protections':r['protection_report'],'J_EB_C2':r['J_EB_C2_distances'],'final_validation':r['final_validation']},indent=2))
if __name__=='__main__':main()
