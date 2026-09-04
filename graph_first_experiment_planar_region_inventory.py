#!/usr/bin/env python3
"""Measurement-only dart-based inventory of the current IAMP planar regions."""
from __future__ import annotations

import argparse
import json
import math
from collections import defaultdict
from pathlib import Path

import graph_first_experiment_crowding_analysis as crowding
import graph_first_experiment_face_poles as old_faces
import graph_first_experiment_face_pull_three as face_pull
import placement_geometry as geometry


def rotate_canonical(walk):
    return list(min(tuple(walk[i:]+walk[:i]) for i in range(len(walk))))


def undirected_walk_key(walk):
    return tuple(old_faces.canonical_cycle(list(walk)))


def enumerate_faces(vertices,edges,coordinates):
    adjacency={vertex:[] for vertex in vertices}
    for first,second in edges:
        adjacency[first].append(second);adjacency[second].append(first)
    for vertex in adjacency:
        adjacency[vertex].sort(key=lambda other:(math.atan2(coordinates[other][1]-coordinates[vertex][1],coordinates[other][0]-coordinates[vertex][0]),other))
    visited=set();raw=[];dart_raw={}
    for first in sorted(adjacency):
        for second in adjacency[first]:
            if (first,second) in visited:continue
            walk=[];darts=[];previous,current=first,second
            while (previous,current) not in visited:
                visited.add((previous,current));walk.append(previous);darts.append((previous,current))
                neighbors=adjacency[current];following=neighbors[(neighbors.index(previous)-1)%len(neighbors)]
                previous,current=current,following
            index=min(range(len(walk)),key=lambda i:tuple(walk[i:]+walk[:i]))
            walk=walk[index:]+walk[:index];darts=darts[index:]+darts[:index]
            signed=old_faces.polygon_area([coordinates[v] for v in walk])
            raw.append({"walk":walk,"darts":darts,"signed_area":signed})
    # IDs use lexical canonical walks, not a geometric or pull-utility ranking.
    raw.sort(key=lambda item:tuple(item["walk"]))
    faces=[];dart_face={}
    for index,item in enumerate(raw,1):
        face_id=f"face:{index:03d}"
        for dart in item["darts"]:dart_face[dart]=face_id
        walk=item["walk"];distinct_vertices=sorted(set(walk));distinct_edges=sorted({tuple(sorted(dart)) for dart in item["darts"]})
        repeated_vertices=len(distinct_vertices)!=len(walk);repeated_edges=len(distinct_edges)!=len(item["darts"])
        simple=not repeated_vertices and not repeated_edges
        exterior=item["signed_area"]<0
        record={"face_id":face_id,"face_type":"EXTERIOR" if exterior else "BOUNDED","is_bounded":not exterior,
                "ordered_boundary_walk":walk,"ordered_boundary_darts":[list(dart) for dart in item["darts"]],
                "distinct_boundary_vertices":distinct_vertices,"distinct_boundary_edges":[list(edge) for edge in distinct_edges],
                "distinct_vertex_count":len(distinct_vertices),"distinct_edge_count":len(distinct_edges),
                "boundary_walk_length":len(walk),"simple_boundary":simple,"repeated_vertices_present":repeated_vertices,
                "repeated_bridge_edges_present":repeated_edges,
                "non_simple_reason":([reason for reason,condition in (("repeated boundary vertices",repeated_vertices),("repeated undirected boundary edges / bridge sides",repeated_edges)) if condition] or None),
                "geometric_area":None,"perimeter":None,"bounding_box":None}
        if not exterior and simple:
            points=[coordinates[v] for v in walk];xs=[p[0] for p in points];ys=[p[1] for p in points]
            record.update({"geometric_area":abs(item["signed_area"]),"perimeter":old_faces.polygon_perimeter(points),
                           "bounding_box":{"min_x":min(xs),"min_y":min(ys),"max_x":max(xs),"max_y":max(ys)}})
        faces.append(record)
    exterior=[face for face in faces if face["face_type"]=="EXTERIOR"]
    if len(exterior)!=1:raise AssertionError(f"expected one exterior face, found {len(exterior)}")
    return faces,dart_face


def add_classification_and_adjacency(faces,edges,dart_face):
    by_id={face["face_id"]:face for face in faces}
    edge_sides=[];neighbor_edges=defaultdict(lambda:defaultdict(list))
    for edge in edges:
        first,second=edge;left=dart_face[(first,second)];right=dart_face[(second,first)]
        record={"edge":[first,second],"face_on_first_to_second_side":left,"face_on_second_to_first_side":right,
                "same_face_on_both_sides":left==right}
        edge_sides.append(record)
        if left!=right:
            neighbor_edges[left][right].append(list(edge));neighbor_edges[right][left].append(list(edge))
    for face in faces:
        face["neighboring_faces"]=[{"face_id":neighbor,"shared_graph_edges":sorted(shared)} for neighbor,shared in sorted(neighbor_edges[face["face_id"]].items())]
        if face["is_bounded"]:
            classification,free,anchors=face_pull.classify_face(face["ordered_boundary_walk"],edges)
            face.update({"boundary_vertex_classification":classification,"free_vertices":free,"anchor_vertices":anchors,
                         "free_count":len(free),"anchor_count":len(anchors)})
        else:
            face.update({"boundary_vertex_classification":[],"free_vertices":[],"anchor_vertices":[],"free_count":None,"anchor_count":None})
    return edge_sides


def match_known(faces,diagnostic):
    by_key={undirected_walk_key(face["ordered_boundary_walk"]):face for face in faces}
    matches=[]
    for prior in diagnostic["selected_faces"]:
        current=by_key.get(undirected_walk_key(prior["ordered_boundary_cycle"]))
        if current is None:raise AssertionError(f"could not match prior {prior['visual_location']} face")
        matches.append({"visual_region":prior["visual_location"],"inventory_face_id":current["face_id"],
                        "prior_diagnostic_face_id":prior["face_id"],"exact_boundary_walk":current["ordered_boundary_walk"],
                        "is_one_actual_bounded_face":current["is_bounded"],"is_composite_cycle":False,
                        "constituent_faces":[current["face_id"]],"simple_boundary":current["simple_boundary"],
                        "area":current["geometric_area"] if current["geometric_area"] is not None else prior["area"],
                        "area_source":"ordinary simple-polygon boundary" if current["simple_boundary"] else "prior diagnostic walk shoelace; current inventory leaves ordinary area undefined",
                        "free_vertices":current["free_vertices"],"anchor_vertices":current["anchor_vertices"],
                        "adjacent_faces":current["neighboring_faces"],"pole_coordinate":prior["pole_coordinate"]})
    return matches


def crowding_face_references(events,faces,edge_sides):
    vertex_faces=defaultdict(list)
    for face in faces:
        for vertex in face["distinct_boundary_vertices"]:vertex_faces[vertex].append(face["face_id"])
    sides={tuple(item["edge"]):(item["face_on_first_to_second_side"],item["face_on_second_to_first_side"]) for item in edge_sides}
    references=[]
    for event in events:
        item_refs=[]
        for item in (event["item_a"],event["item_b"]):
            if item.startswith("vertex:"):
                vertex=item[len("vertex:"):];face_ids=sorted(vertex_faces[vertex])
            else:
                edge=tuple(sorted(item[len("edge:"):].split("|")));face_ids=sorted(set(sides[edge]))
            item_refs.append({"item":item,"boundary_face_ids":face_ids})
        references.append({"violation_id":event["violation_id"],"coordinate":event["crowding_point"],"item_a":item_refs[0],"item_b":item_refs[1]})
    return references


def render(vertices,edges,coordinates,faces,known,crowding_events):
    xs=[p[0] for p in coordinates.values()];ys=[p[1] for p in coordinates.values()];margin=100
    minx,miny=min(xs)-margin,min(ys)-margin;width,height=max(xs)-min(xs)+2*margin,max(ys)-min(ys)+2*margin
    palette=["#bbdefb","#c8e6c9","#ffe0b2","#e1bee7","#b2dfdb","#f8bbd0","#dcedc8","#fff9c4","#d1c4e9","#b3e5fc","#ffccbc"]
    out=[f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="{minx:.17g} {miny:.17g} {width:.17g} {height:.17g}" width="1400" height="1000">',f'<rect x="{minx:.17g}" y="{miny:.17g}" width="{width:.17g}" height="{height:.17g}" fill="white"/>',f'<text x="{minx+15:.17g}" y="{miny+28:.17g}" font-family="sans-serif" font-size="18">IAMP global planar-region inventory — geometry unchanged</text>']
    known_by_id={item["inventory_face_id"]:item["visual_region"] for item in known}
    for index,face in enumerate(face for face in faces if face["is_bounded"] and face["simple_boundary"]):
        points=' '.join(f'{coordinates[v][0]:.17g},{coordinates[v][1]:.17g}' for v in face["ordered_boundary_walk"])
        out.append(f'<polygon points="{points}" fill="{palette[index%len(palette)]}" fill-opacity=".25" stroke="#aaa" stroke-width=".6"/>')
    for a,b in edges:out.append(f'<line x1="{coordinates[a][0]:.17g}" y1="{coordinates[a][1]:.17g}" x2="{coordinates[b][0]:.17g}" y2="{coordinates[b][1]:.17g}" stroke="#666" stroke-width="1.5"/>')
    # Label simple faces at a deterministic safe interior point; use the prior
    # safe pole for the known non-simple RA_SW face.
    prior_poles={item["inventory_face_id"]:item["pole_coordinate"] for item in known}
    for face in faces:
        if face["face_type"]=="EXTERIOR":continue
        if face["face_id"] in prior_poles:point=prior_poles[face["face_id"]]
        elif face["simple_boundary"]:point,_,_=old_faces.safe_pole([coordinates[v] for v in face["ordered_boundary_walk"]])
        else:
            points=[coordinates[v] for v in face["distinct_boundary_vertices"]];point=[sum(x for x,_ in points)/len(points),sum(y for _,y in points)/len(points)]
        suffix=f' ({known_by_id[face["face_id"]]})' if face["face_id"] in known_by_id else ''
        out.append(f'<text x="{point[0]:.17g}" y="{point[1]:.17g}" text-anchor="middle" font-family="sans-serif" font-size="13" font-weight="bold" fill="#222">{face["face_id"]}{suffix}</text>')
    # Topological classification is consistent at a vertex across incident
    # bounded faces here; green ring denotes FREE somewhere, orange ANCHOR.
    free_any={v for face in faces if face["is_bounded"] for v in face["free_vertices"]};anchor_any={v for face in faces if face["is_bounded"] for v in face["anchor_vertices"]}-free_any
    for vertex in sorted(vertices):
        x,y=coordinates[vertex];fill="#1565c0" if vertices[vertex]["type"]=="COMPONENT" else "#d84315"
        out.append(f'<circle cx="{x:.17g}" cy="{y:.17g}" r="6" fill="{fill}" stroke="white"/>')
        if vertex in free_any:out.append(f'<circle cx="{x:.17g}" cy="{y:.17g}" r="10" fill="none" stroke="#2e7d32" stroke-width="2.5"/>')
        elif vertex in anchor_any:out.append(f'<rect x="{x-9:.17g}" y="{y-9:.17g}" width="18" height="18" fill="none" stroke="#ef6c00" stroke-width="2"/>')
        out.append(f'<text x="{x+10:.17g}" y="{y-8:.17g}" font-family="sans-serif" font-size="11">{vertices[vertex]["name"]}</text>')
    for event in crowding_events:
        x,y=event["coordinate"];out.append(f'<circle cx="{x:.17g}" cy="{y:.17g}" r="4" fill="#d50000" stroke="white"/><text x="{x+5:.17g}" y="{y+12:.17g}" font-family="sans-serif" font-size="9" fill="#d50000">{event["violation_id"]}</text>')
    exterior=next(face for face in faces if not face["is_bounded"])
    out.append(f'<text x="{minx+18:.17g}" y="{miny+55:.17g}" font-family="sans-serif" font-size="14" font-weight="bold">{exterior["face_id"]} — EXTERIOR</text>')
    out.append('</svg>');return '\n'.join(out)+'\n'


def run(graph_path,direct_path,face_path,crowding_path,output_dir):
    vertices,edges,coordinates,_,_=crowding.load_inputs(graph_path,direct_path);original=dict(coordinates)
    diagnostics=geometry.graph_geometry_diagnostics(coordinates,edges)
    if diagnostics["proper_unrelated_edge_crossing_count"]:raise AssertionError("input is not zero-crossing")
    faces,dart_face=enumerate_faces(vertices,edges,coordinates);edge_sides=add_classification_and_adjacency(faces,edges,dart_face)
    known=match_known(faces,json.loads(face_path.read_text()));event_doc=json.loads(crowding_path.read_text());event_refs=crowding_face_references(event_doc["events"],faces,edge_sides)
    report={"schema":"graph-relax.planar-region-inventory.v1","source_coordinates":f"{direct_path}#final_coordinates","vertex_count":len(vertices),"electrical_incidence_count":len(edges),"proper_crossing_count":diagnostics["proper_unrelated_edge_crossing_count"],"coordinates_unchanged":coordinates==original,
            "face_id_assignment":"lexicographic canonical directed boundary walk; no area/utility ranking","total_faces":len(faces),"bounded_faces":sum(face["is_bounded"] for face in faces),"exterior_face_id":next(face["face_id"] for face in faces if not face["is_bounded"]),
            "bounded_faces_with_free_vertices":sum(face["is_bounded"] and face["free_count"]>0 for face in faces),"bounded_faces_without_free_vertices":sum(face["is_bounded"] and face["free_count"]==0 for face in faces),
            "faces":faces,"graph_edge_face_sides":edge_sides,"known_region_matches":known,"existing_crowding_event_face_references":event_refs}
    output_dir.mkdir(parents=True,exist_ok=True);(output_dir/'iamp_all_faces.svg').write_text(render(vertices,edges,coordinates,faces,known,event_refs));(output_dir/'iamp_all_faces_report.json').write_text(json.dumps(report,indent=2,sort_keys=True)+'\n');return report


def main():
    parser=argparse.ArgumentParser();base=Path('output/graph_first');parser.add_argument('--graph',type=Path,default=base/'iamp_graph.json');parser.add_argument('--direct-report',type=Path,default=base/'iamp_nearness_rotation_direct_report.json');parser.add_argument('--face-report',type=Path,default=base/'iamp_face_poles_report.json');parser.add_argument('--crowding-report',type=Path,default=base/'iamp_crowding_points_report.json');parser.add_argument('--output',type=Path,default=base);a=parser.parse_args();r=run(a.graph,a.direct_report,a.face_report,a.crowding_report,a.output)
    for f in r['faces']:print(f['face_id'],f['face_type'],f['distinct_vertex_count'],f['distinct_edge_count'],f['simple_boundary'],f['geometric_area'],f['free_count'],f['anchor_count'],[n['face_id'] for n in f['neighboring_faces']])
    print('total',r['total_faces'],'bounded',r['bounded_faces'],'with FREE',r['bounded_faces_with_free_vertices'],'without FREE',r['bounded_faces_without_free_vertices']);print([(x['visual_region'],x['inventory_face_id']) for x in r['known_region_matches']])

if __name__=='__main__':main()
