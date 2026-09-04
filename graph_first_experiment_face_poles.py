#!/usr/bin/env python3
"""Measurement-only face-pole diagnostic for the current planar IAMP drawing."""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import placement_geometry as geometry


def canonical_cycle(cycle):
    variants=[]
    for seq in (cycle, list(reversed(cycle))):
        for i in range(len(seq)):
            variants.append(tuple(seq[i:]+seq[:i]))
    return list(min(variants))


def polygon_area(points):
    return .5*sum(points[i][0]*points[(i+1)%len(points)][1]-points[(i+1)%len(points)][0]*points[i][1]
                    for i in range(len(points)))


def polygon_perimeter(points):
    return sum(math.dist(points[i],points[(i+1)%len(points)]) for i in range(len(points)))


def enumerate_geometric_faces(vertices, edges, coordinates):
    adjacency={node:[] for node in vertices}
    for first,second in edges:
        adjacency[first].append(second); adjacency[second].append(first)
    for node in adjacency:
        adjacency[node].sort(key=lambda other:(math.atan2(coordinates[other][1]-coordinates[node][1],
                                                         coordinates[other][0]-coordinates[node][0]),other))
    visited=set(); raw=[]
    for first in sorted(adjacency):
        for second in adjacency[first]:
            if (first,second) in visited: continue
            walk=[]; previous,current=first,second
            while (previous,current) not in visited:
                visited.add((previous,current)); walk.append(previous)
                neighbors=adjacency[current]
                following=neighbors[(neighbors.index(previous)-1)%len(neighbors)]
                previous,current=current,following
            raw.append(walk)
    faces=[]
    for walk in raw:
        points=[coordinates[node] for node in walk]
        signed=polygon_area(points)
        faces.append({"boundary_cycle":canonical_cycle(walk),"signed_area":signed,
                      "area":abs(signed),"is_bounded":signed>0})
    faces.sort(key=lambda face:(not face["is_bounded"],-face["area"],face["boundary_cycle"]))
    for index,face in enumerate(faces,1):
        points=[coordinates[node] for node in face["boundary_cycle"]]
        xs=[point[0] for point in points]; ys=[point[1] for point in points]
        face.update({"face_id":f"face:{index:03d}","perimeter":polygon_perimeter(points),
                     "bounding_box":{"min_x":min(xs),"min_y":min(ys),"max_x":max(xs),"max_y":max(ys)}})
    return faces


def cross(a,b,c):
    return (b[0]-a[0])*(c[1]-a[1])-(b[1]-a[1])*(c[0]-a[0])


def point_strictly_in_triangle(p,a,b,c,tol=1e-12):
    values=(cross(a,b,p),cross(b,c,p),cross(c,a,p))
    return all(value>tol for value in values) or all(value < -tol for value in values)


def triangulate(points):
    """Deterministic ear clipping; returns triangles wholly inside the polygon."""
    indices=list(range(len(points))); triangles=[]
    orientation=1 if polygon_area(points)>0 else -1
    while len(indices)>3:
        ear=None
        for position,index in enumerate(indices):
            before=indices[position-1]; after=indices[(position+1)%len(indices)]
            if orientation*cross(points[before],points[index],points[after]) <= geometry.GEOMETRY_TOLERANCE:
                continue
            if any(point_strictly_in_triangle(points[other],points[before],points[index],points[after])
                   for other in indices if other not in (before,index,after)):
                continue
            ear=(position,(before,index,after)); break
        if ear is None: raise ValueError("face boundary could not be triangulated")
        position,triangle=ear; triangles.append(triangle); indices.pop(position)
    triangles.append(tuple(indices)); return triangles


def triangle_incenter(a,b,c):
    wa=math.dist(b,c); wb=math.dist(a,c); wc=math.dist(a,b); total=wa+wb+wc
    return ((wa*a[0]+wb*b[0]+wc*c[0])/total,(wa*a[1]+wb*b[1]+wc*c[1])/total)


def boundary_clearance(point,points):
    return min(geometry.point_segment_distance(point,points[i],points[(i+1)%len(points)])
               for i in range(len(points)))


def safe_pole(points):
    candidates=[]
    for triangle in triangulate(points):
        point=triangle_incenter(*(points[index] for index in triangle))
        candidates.append((boundary_clearance(point,points),triangle,point))
    clearance,triangle,point=max(candidates,key=lambda item:(item[0],tuple(-i for i in item[1])))
    if clearance<=geometry.GEOMETRY_TOLERANCE: raise ValueError("failed to produce strict interior pole")
    return point,clearance,list(triangle)


def analyze(graph_path,direct_report_path):
    graph=json.loads(graph_path.read_text()); direct=json.loads(direct_report_path.read_text())
    vertices={item["id"]:item for item in graph["vertices"]}
    edges=[(item["source"],item["target"]) for item in graph["edges"]]
    coordinates={node:tuple(point) for node,point in direct["final_coordinates"].items()}
    diagnostics=geometry.graph_geometry_diagnostics(coordinates,edges)
    if len(vertices)!=34 or len(edges)!=44 or diagnostics["proper_unrelated_edge_crossing_count"]!=0:
        raise ValueError("source is not the required 34-vertex/44-incidence zero-crossing drawing")
    faces=enumerate_geometric_faces(vertices,edges,coordinates)
    if len(faces)!=len(edges)-len(vertices)+2: raise AssertionError("face count violates Euler formula")
    targets={"RMIN1":"component:RMIN1","RA_SW":"net:RA_SW","RMIN4":"component:RMIN4"}
    selected=[]
    for label,node in targets.items():
        matches=[face for face in faces if face["is_bounded"] and node in face["boundary_cycle"]]
        if len(matches)!=1: raise ValueError(f"{label} belongs to {len(matches)} bounded face boundaries")
        face=matches[0]; points=[coordinates[item] for item in face["boundary_cycle"]]
        pole,clearance,triangle=safe_pole(points)
        cables=[]
        for vertex in face["boundary_cycle"]:
            point=coordinates[vertex]; dx=pole[0]-point[0]; dy=pole[1]-point[1]; distance=math.hypot(dx,dy)
            cables.append({"vertex":vertex,"coordinate":list(point),"distance_to_pole":distance,
                           "direction_toward_pole":[dx/distance,dy/distance]})
        boundary_edges={frozenset((face["boundary_cycle"][i],face["boundary_cycle"][(i+1)%len(points)]))
                        for i in range(len(points))}
        inside_vertices=[]
        # A face cannot contain graph geometry. Verify using strict ray casting
        # and edge midpoints, excluding its own boundary.
        def inside(point):
            result=False
            for i,a in enumerate(points):
                b=points[(i+1)%len(points)]
                if (a[1]>point[1]) != (b[1]>point[1]):
                    x=a[0]+(point[1]-a[1])*(b[0]-a[0])/(b[1]-a[1])
                    if x>point[0]: result=not result
            return result and boundary_clearance(point,points)>geometry.GEOMETRY_TOLERANCE
        boundary_set=set(face["boundary_cycle"])
        for vertex,point in coordinates.items():
            if vertex not in boundary_set and inside(point): inside_vertices.append(vertex)
        inside_edges=[]
        for edge in edges:
            if frozenset(edge) in boundary_edges: continue
            midpoint=((coordinates[edge[0]][0]+coordinates[edge[1]][0])/2,
                      (coordinates[edge[0]][1]+coordinates[edge[1]][1])/2)
            if inside(midpoint): inside_edges.append(list(edge))
        selected.append({"visual_location":label,"face_id":face["face_id"],
                         "ordered_boundary_cycle":face["boundary_cycle"],
                         "boundary_vertex_count":len(face["boundary_cycle"]),
                         "area":face["area"],"perimeter":face["perimeter"],
                         "bounding_box":face["bounding_box"],"pole_coordinate":list(pole),
                         "pole_to_boundary_minimum_clearance":clearance,
                         "pole_source_triangle_indices":triangle,"boundary_vertex_cables":cables,
                         "strictly_interior_graph_vertices":inside_vertices,
                         "strictly_interior_graph_edges":inside_edges})
    return vertices,edges,coordinates,faces,selected,diagnostics


def render(vertices,edges,coordinates,selected):
    xs=[p[0] for p in coordinates.values()]; ys=[p[1] for p in coordinates.values()]
    margin=100; minx=min(xs)-margin; miny=min(ys)-margin; width=max(xs)-min(xs)+2*margin; height=max(ys)-min(ys)+2*margin
    colors=[("#1565c0","#90caf9"),("#2e7d32","#a5d6a7"),("#8e24aa","#ce93d8")]
    out=[f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="{minx:.17g} {miny:.17g} {width:.17g} {height:.17g}" width="1400" height="1000">',
         f'<rect x="{minx:.17g}" y="{miny:.17g}" width="{width:.17g}" height="{height:.17g}" fill="white"/>',
         f'<text x="{minx+15:.17g}" y="{miny+27:.17g}" font-family="sans-serif" font-size="18">IAMP face-pole diagnostic — graph unchanged</text>']
    for face,(_,fill) in zip(selected,colors):
        points=' '.join(f'{coordinates[v][0]:.17g},{coordinates[v][1]:.17g}' for v in face["ordered_boundary_cycle"])
        out.append(f'<polygon points="{points}" fill="{fill}" fill-opacity="0.20" stroke="none"/>')
    for a,b in edges:
        out.append(f'<line x1="{coordinates[a][0]:.17g}" y1="{coordinates[a][1]:.17g}" x2="{coordinates[b][0]:.17g}" y2="{coordinates[b][1]:.17g}" stroke="#777" stroke-width="1.5"/>')
    for face,(stroke,_) in zip(selected,colors):
        cycle=face["ordered_boundary_cycle"]
        for i,a in enumerate(cycle):
            b=cycle[(i+1)%len(cycle)]
            out.append(f'<line x1="{coordinates[a][0]:.17g}" y1="{coordinates[a][1]:.17g}" x2="{coordinates[b][0]:.17g}" y2="{coordinates[b][1]:.17g}" stroke="{stroke}" stroke-width="4" stroke-opacity="0.75"/>')
        px,py=face["pole_coordinate"]
        for vertex in dict.fromkeys(cycle):
            x,y=coordinates[vertex]
            out.append(f'<line x1="{px:.17g}" y1="{py:.17g}" x2="{x:.17g}" y2="{y:.17g}" stroke="{stroke}" stroke-width="1" stroke-dasharray="6 5" opacity="0.6"/>')
        out.extend([f'<circle cx="{px:.17g}" cy="{py:.17g}" r="9" fill="white" stroke="{stroke}" stroke-width="4"/>',
                    f'<path d="M {px-5:.17g} {py:.17g} H {px+5:.17g} M {px:.17g} {py-5:.17g} V {py+5:.17g}" stroke="{stroke}" stroke-width="2"/>',
                    f'<text x="{px+12:.17g}" y="{py-10:.17g}" font-family="sans-serif" font-size="14" font-weight="bold" fill="{stroke}">{face["visual_location"]} / {face["face_id"]} pole</text>'])
    for node in sorted(vertices):
        x,y=coordinates[node]; fill="#1565c0" if vertices[node]["type"]=="COMPONENT" else "#d84315"
        out.append(f'<circle cx="{x:.17g}" cy="{y:.17g}" r="5" fill="{fill}" stroke="white" stroke-width="1"/>')
        out.append(f'<text x="{x+8:.17g}" y="{y-7:.17g}" font-family="sans-serif" font-size="11">{node.split(":",1)[1]}</text>')
    out.append('</svg>'); return '\n'.join(out)+'\n'


def run(graph_path,direct_report_path,output_dir):
    vertices,edges,coordinates,faces,selected,diagnostics=analyze(graph_path,direct_report_path)
    report={"schema":"graph-relax.face-poles.v1","source_coordinates":f"{direct_report_path}#final_coordinates",
            "vertex_count":len(vertices),"electrical_incidence_count":len(edges),
            "proper_crossing_count":diagnostics["proper_unrelated_edge_crossing_count"],
            "coordinates_unchanged":True,"face_count":len(faces),"bounded_face_count":sum(f["is_bounded"] for f in faces),
            "faces":faces,"selected_faces":selected}
    output_dir.mkdir(parents=True,exist_ok=True)
    (output_dir/'iamp_face_poles.svg').write_text(render(vertices,edges,coordinates,selected))
    (output_dir/'iamp_face_poles_report.json').write_text(json.dumps(report,indent=2,sort_keys=True)+'\n')
    return report


def main():
    parser=argparse.ArgumentParser(); base=Path('output/graph_first')
    parser.add_argument('--graph',type=Path,default=base/'iamp_graph.json')
    parser.add_argument('--direct-report',type=Path,default=base/'iamp_nearness_rotation_direct_report.json')
    parser.add_argument('--output',type=Path,default=base); args=parser.parse_args()
    report=run(args.graph,args.direct_report,args.output)
    for face in report['selected_faces']:
        print(face['visual_location'],face['face_id'],face['ordered_boundary_cycle'],face['area'],face['perimeter'],face['pole_coordinate'],face['pole_to_boundary_minimum_clearance'])

if __name__=='__main__': main()
