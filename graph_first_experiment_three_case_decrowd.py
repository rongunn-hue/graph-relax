#!/usr/bin/env python3
"""One polygon, non-polygon-node, and edge-edge de-crowding pass."""
import argparse,json,math
from pathlib import Path
import graph_first_experiment_planar_region_inventory as inventory
import graph_first_experiment_polygon_area_expansion as poly
import placement_geometry as geo

TOL=1e-7
def closest_point(p,a,b):
 dx,dy=b[0]-a[0],b[1]-a[1];q=dx*dx+dy*dy
 t=0 if q==0 else max(0,min(1,((p[0]-a[0])*dx+(p[1]-a[1])*dy)/q))
 return (a[0]+t*dx,a[1]+t*dy)
def segment_closest(a,b,c,d):
 candidates=[(a,closest_point(a,c,d)),(b,closest_point(b,c,d)),(closest_point(c,a,b),c),(closest_point(d,a,b),d)]
 return min(candidates,key=lambda z:(math.dist(*z),z))
def hard_valid(vertices,edges,co):
 d=geo.graph_geometry_diagnostics(co,edges);r=geo.validate_rebuilt_edge_segments(co,edges,geo.rebuild_straight_edge_segments(co,edges))
 return len(co)==34 and len(edges)==44 and r['valid'] and not d['proper_unrelated_edge_crossing_count'] and not d['coincident_vertex_pair_count'] and not d['vertex_on_unrelated_edge_interior_count'],d,r
def node_event(v,co,edges):
 neighbors={b for a,b in edges if a==v}|{a for a,b in edges if b==v};items=[]
 for q in sorted(co):
  if q!=v and q not in neighbors:items.append((math.dist(co[v],co[q]),'VERTEX',q,co[q]))
 for a,b in sorted(tuple(sorted(e)) for e in edges):
  if v not in (a,b):
   p=closest_point(co[v],co[a],co[b]);items.append((math.dist(co[v],p),'EDGE',[a,b],p))
 distance,kind,item,point=min(items,key=lambda x:(x[0],x[1],str(x[2])))
 return {'vertex':v,'clearance':distance,'offender_type':kind,'offender':item,'closest_point':list(point)}
def nonpolygon_map(vertices,edges,co,faces):
 represented={v for f in faces if f['is_bounded'] and f['simple_boundary'] for v in f['distinct_boundary_vertices']}
 rows=[node_event(v,co,edges) for v in sorted(vertices) if v not in represented]
 rows.sort(key=lambda x:(x['clearance'],x['vertex']));
 for i,x in enumerate(rows,1):x['rank']=i
 return rows
def edge_map(edges,co):
 es=sorted(tuple(sorted(e)) for e in edges);rows=[]
 for i,e in enumerate(es):
  for f in es[i+1:]:
   if set(e)&set(f):continue
   p,q=segment_closest(co[e[0]],co[e[1]],co[f[0]],co[f[1]])
   rows.append({'edge_a':list(e),'edge_b':list(f),'clearance':math.dist(p,q),'closest_point_a':list(p),'closest_point_b':list(q)})
 rows.sort(key=lambda x:(x['clearance'],x['edge_a'],x['edge_b']))
 for i,x in enumerate(rows,1):x['rank']=i
 return rows
def unit_away(point,other,stable):
 dx,dy=point[0]-other[0],point[1]-other[1];q=math.hypot(dx,dy)
 if q>TOL:return dx/q,dy/q
 angle=(sum(map(ord,stable))%360)*math.pi/180;return math.cos(angle),math.sin(angle)
def node_pass(vertices,edges,co,faces):
 initial=nonpolygon_map(vertices,edges,co,faces);attempts=[]
 for base in initial:
  v=base['vertex'];current=node_event(v,co,edges);dx,dy=unit_away(co[v],tuple(current['closest_point']),v);distance=current['clearance'];candidate=dict(co);candidate[v]=(co[v][0]+distance*dx,co[v][1]+distance*dy);ok,d,r=hard_valid(vertices,edges,candidate);after=node_event(v,candidate,edges)
  accepted=ok and after['clearance']>current['clearance']+TOL
  attempts.append({'vertex':v,'before':current,'proposed_displacement':distance,'after_if_proposed':after,'accepted':accepted,'rejection_reason':None if accepted else ('HARD_GEOMETRY_INVALID' if not ok else 'TARGET_CLEARANCE_NOT_IMPROVED')})
  if accepted:co=candidate
 return co,initial,attempts,nonpolygon_map(vertices,edges,co,faces)
def edge_pass(vertices,edges,co):
 initial=edge_map(edges,co);minimum=initial[0]['clearance'];targets=[x for x in initial if abs(x['clearance']-minimum)<=TOL];attempts=[]
 for target in targets:
  current=next((x for x in edge_map(edges,co) if x['edge_a']==target['edge_a'] and x['edge_b']==target['edge_b']),None)
  if current is None:continue
  chosen=None;trials=[]
  for side,edge,otherpoint in [('A',current['edge_a'],current['closest_point_b']),('B',current['edge_b'],current['closest_point_a'])]:
   for v in edge:
    dx,dy=unit_away(co[v],tuple(otherpoint),v);candidate=dict(co);candidate[v]=(co[v][0]+current['clearance']*dx,co[v][1]+current['clearance']*dy);ok,d,r=hard_valid(vertices,edges,candidate);after=next(x for x in edge_map(edges,candidate) if x['edge_a']==current['edge_a'] and x['edge_b']==current['edge_b']);accepted=ok and after['clearance']>current['clearance']+TOL;trials.append({'endpoint':v,'side':side,'clearance_after':after['clearance'],'hard_valid':ok,'accepted':accepted})
    if accepted:chosen=(v,candidate,after);break
   if chosen:break
  attempts.append({'target':current,'endpoint_trials':trials,'accepted':chosen is not None,'endpoint_moved':None if chosen is None else chosen[0],'clearance_after':current['clearance'] if chosen is None else chosen[2]['clearance']})
  if chosen:co=chosen[1]
 return co,initial,attempts,edge_map(edges,co)
def render(vertices,edges,co,faces,polygon_rows,node_rows,edge_rows,title,path):
 xs=[p[0] for p in co.values()];ys=[p[1] for p in co.values()];m=90;minx,miny=min(xs)-m,min(ys)-m;w=max(xs)-min(xs)+2*m;h=max(ys)-min(ys)+2*m;crowded={x['face_id'] for x in polygon_rows[:3]}
 z=[f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="{minx} {miny} {w} {h}" width="1400" height="1000"><rect x="{minx}" y="{miny}" width="{w}" height="{h}" fill="white"/><text x="{minx+15}" y="{miny+25}" font-size="18">{title}</text>']
 for f in faces:
  if f['face_id'] in crowded:z.append('<polygon points="'+' '.join(f'{co[v][0]},{co[v][1]}' for v in f['ordered_boundary_walk'])+'" fill="#ffca28" fill-opacity=".22" stroke="#f57f17"/>')
 for a,b in edges:z.append(f'<line x1="{co[a][0]}" y1="{co[a][1]}" x2="{co[b][0]}" y2="{co[b][1]}" stroke="#666" stroke-width="1.4"/>')
 for event in node_rows[:3]:
  x,y=co[event['vertex']];z.append(f'<circle cx="{x}" cy="{y}" r="11" fill="none" stroke="#8e24aa" stroke-width="2"/>')
 for event in edge_rows[:1]:
  for p in (event['closest_point_a'],event['closest_point_b']):z.append(f'<circle cx="{p[0]}" cy="{p[1]}" r="5" fill="#d50000"/>')
 for v in sorted(vertices):
  x,y=co[v];z.append(f'<circle cx="{x}" cy="{y}" r="5" fill="#1565c0"/><text x="{x+7}" y="{y-7}" font-size="10">{vertices[v]["name"]}</text>')
 z.append('</svg>');path.write_text('\n'.join(z)+'\n')
def run(graph_path,direct_path,outdir):
 gd=json.loads(graph_path.read_text());dr=json.loads(direct_path.read_text());vertices={x['id']:x for x in gd['vertices']};edges=[(x['source'],x['target']) for x in gd['edges']];start={v:tuple(p) for v,p in dr['final_coordinates'].items()};ok,sd,re=hard_valid(vertices,edges,start)
 if not ok:raise AssertionError('authoritative zero-crossing start invalid')
 faces,darts=inventory.enumerate_faces(vertices,edges,start);inventory.add_classification_and_adjacency(faces,edges,darts);p_before=poly.ranking(faces,start);targets=[x['face_id'] for x in p_before[:3]];ext,opens,paths,order=poly.dependencies(faces,p_before,targets);by={f['face_id']:f for f in faces};co=dict(start);pops=[]
 node_before=nonpolygon_map(vertices,edges,start,faces);edge_before=edge_map(edges,start)
 for fid in order:
  co,result=poly.expand(edges,by[fid],co);result.update({'face_id':fid,'boundary_walk':by[fid]['ordered_boundary_walk']});pops.append(result)
 polygon_co=dict(co);p_after_polygon=poly.ranking(faces,co);node_after_polygon=nonpolygon_map(vertices,edges,co,faces);edge_after_polygon=edge_map(edges,co)
 co,node_stage_before,node_attempts,node_after=node_pass(vertices,edges,co,faces);node_co=dict(co);edge_before_stage=edge_map(edges,co)
 co,edge_stage_initial,edge_attempts,edge_after=edge_pass(vertices,edges,co);final=co;p_final=poly.ranking(faces,final);ok,fd,re=hard_valid(vertices,edges,final)
 if not ok:raise AssertionError('final graph invalid')
 outdir.mkdir(parents=True,exist_ok=True);render(vertices,edges,start,faces,p_before,node_before,edge_before,'START — THREE CROWDING MAPS',outdir/'iamp_three_case_decrowd_start.svg');render(vertices,edges,polygon_co,faces,p_after_polygon,node_after_polygon,edge_after_polygon,'AFTER POLYGON AREA PASS',outdir/'iamp_three_case_decrowd_polygon.svg');render(vertices,edges,node_co,faces,poly.ranking(faces,node_co),node_after,edge_before_stage,'AFTER NON-POLYGON NODE PASS',outdir/'iamp_three_case_decrowd_nodes.svg');render(vertices,edges,final,faces,p_final,nonpolygon_map(vertices,edges,final,faces),edge_after,'FINAL — AFTER EDGE PASS',outdir/'iamp_three_case_decrowd_final.svg')
 report={'schema':'graph-relax.three-case-decrowd.v1','source_coordinates':f'{direct_path}#final_coordinates','starting_validation':{'vertices':34,'edges':44,'proper_crossings':0,'coincidences':0,'vertex_on_unrelated_edge':0,'connectivity_unchanged':True},'initial_crowding_map':{'polygon':p_before,'non_polygon_nodes':node_before,'edge_edge':edge_before},'polygon_stage':{'targets':targets,'dependency_paths':paths,'processing_order':order,'operations':pops,'ranking_after':p_after_polygon},'non_polygon_node_stage':{'map_before':node_stage_before,'attempts':node_attempts,'map_after':node_after},'edge_stage':{'map_before':edge_stage_initial,'worst_events_attempted':edge_attempts,'map_after':edge_after},'final_polygon_ranking':p_final,'final_validation':{'vertices':34,'edges':44,'proper_crossings':fd['proper_unrelated_edge_crossing_count'],'coincidences':fd['coincident_vertex_pair_count'],'vertex_on_unrelated_edge':fd['vertex_on_unrelated_edge_interior_count'],'connectivity_unchanged':re['valid']},'worst_remaining':{'polygon':p_final[0],'non_polygon_node':nonpolygon_map(vertices,edges,final,faces)[0],'edge_edge':edge_after[0]},'final_coordinates':{v:list(final[v]) for v in sorted(final)}};(outdir/'iamp_three_case_decrowd_report.json').write_text(json.dumps(report,indent=2,sort_keys=True)+'\n');return report
def main():
 p=argparse.ArgumentParser();b=Path('output/graph_first');p.add_argument('--graph',type=Path,default=b/'iamp_graph.json');p.add_argument('--direct',type=Path,default=b/'iamp_nearness_rotation_direct_report.json');p.add_argument('--output',type=Path,default=b);a=p.parse_args();r=run(a.graph,a.direct,a.output);print(json.dumps({'polygon_order':r['polygon_stage']['processing_order'],'node_attempts':[(x['vertex'],x['accepted']) for x in r['non_polygon_node_stage']['attempts']],'edge_attempts':r['edge_stage']['worst_events_attempted'],'worst_remaining':r['worst_remaining'],'final_validation':r['final_validation']},indent=2))
if __name__=='__main__':main()
