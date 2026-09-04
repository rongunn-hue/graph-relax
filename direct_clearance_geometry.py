"""Analytic finite-segment clearance boundary for direct endpoint rotation."""

from __future__ import annotations

import math

import placement_geometry as geometry


TAU = 2.0 * math.pi


def normalize(angle):
    return angle % TAU


def signed_delta(angle, origin):
    return (angle-origin+math.pi) % TAU - math.pi


def circle_circle_points(first, first_radius, second, second_radius, tolerance):
    dx, dy = second[0]-first[0], second[1]-first[1]
    distance = math.hypot(dx, dy)
    if distance <= tolerance:
        return []
    if distance > first_radius+second_radius+tolerance:
        return []
    if distance < abs(first_radius-second_radius)-tolerance:
        return []
    along = (first_radius**2-second_radius**2+distance**2)/(2*distance)
    height_squared = first_radius**2-along**2
    if height_squared < -tolerance:
        return []
    height = math.sqrt(max(0.0, height_squared))
    base = (first[0]+along*dx/distance, first[1]+along*dy/distance)
    perpendicular = (-dy/distance, dx/distance)
    points = [(base[0]+height*perpendicular[0], base[1]+height*perpendicular[1])]
    if height > tolerance:
        points.append((base[0]-height*perpendicular[0], base[1]-height*perpendicular[1]))
    return points


def circle_segment_points(center, radius, first, second, tolerance):
    dx, dy = second[0]-first[0], second[1]-first[1]
    fx, fy = first[0]-center[0], first[1]-center[1]
    a = dx*dx+dy*dy
    if a <= tolerance*tolerance:
        return [first] if abs(math.dist(first, center)-radius) <= tolerance else []
    b = 2*(fx*dx+fy*dy)
    c = fx*fx+fy*fy-radius*radius
    discriminant = b*b-4*a*c
    if discriminant < -tolerance:
        return []
    root = math.sqrt(max(0.0, discriminant))
    parameters = [(-b-root)/(2*a), (-b+root)/(2*a)]
    points = []
    for value in parameters:
        if -tolerance <= value <= 1+tolerance:
            value = max(0.0, min(1.0, value))
            point = (first[0]+value*dx, first[1]+value*dy)
            if not any(math.dist(point, old) <= tolerance for old in points):
                points.append(point)
    return points


def _cap_contains(point, endpoint, other, tolerance):
    # Endpoint cap is the semicircle facing away from the segment interior.
    return ((point[0]-endpoint[0])*(other[0]-endpoint[0]) +
            (point[1]-endpoint[1])*(other[1]-endpoint[1])) <= tolerance


def capsule_boundary_angles(pivot, motion_radius, first, second, clearance,
                            tolerance=geometry.GEOMETRY_TOLERANCE):
    """Return exact silhouette/finite-radius boundary angles and their source."""
    dx, dy = second[0]-first[0], second[1]-first[1]
    length = math.hypot(dx, dy)
    if length <= tolerance:
        raise ValueError("crossed segment is degenerate")
    normal = (-dy/length, dx/length)
    side_segments = [
        ((first[0]+normal[0]*clearance, first[1]+normal[1]*clearance),
         (second[0]+normal[0]*clearance, second[1]+normal[1]*clearance)),
        ((first[0]-normal[0]*clearance, first[1]-normal[1]*clearance),
         (second[0]-normal[0]*clearance, second[1]-normal[1]*clearance)),
    ]
    records = []
    for label, endpoint, other in (("CAP_A", first, second), ("CAP_B", second, first)):
        q = math.dist(pivot, endpoint)
        if q > clearance+tolerance:
            phi = math.atan2(endpoint[1]-pivot[1], endpoint[0]-pivot[0])
            offset = math.asin(clearance/q)
            tangent_distance = math.sqrt(max(0.0, q*q-clearance*clearance))
            if tangent_distance <= motion_radius+tolerance:
                for angle in (phi-offset, phi+offset):
                    point = (pivot[0]+tangent_distance*math.cos(angle),
                             pivot[1]+tangent_distance*math.sin(angle))
                    if _cap_contains(point, endpoint, other, tolerance):
                        records.append({"angle": normalize(angle), "feature": label,
                                        "point": point})
        for point in circle_circle_points(pivot, motion_radius, endpoint, clearance, tolerance):
            if _cap_contains(point, endpoint, other, tolerance):
                records.append({"angle": normalize(math.atan2(point[1]-pivot[1],
                                                               point[0]-pivot[0])),
                                "feature": label, "point": point})
    for side_index, (side_first, side_second) in enumerate(side_segments):
        for point in (side_first, side_second):
            if math.dist(pivot, point) <= motion_radius+tolerance:
                records.append({"angle": normalize(math.atan2(point[1]-pivot[1],
                                                               point[0]-pivot[0])),
                                "feature": f"STRIP_{side_index}", "point": point})
        for point in circle_segment_points(pivot, motion_radius, side_first, side_second, tolerance):
            records.append({"angle": normalize(math.atan2(point[1]-pivot[1],
                                                           point[0]-pivot[0])),
                            "feature": f"STRIP_{side_index}", "point": point})
    unique = []
    for record in sorted(records, key=lambda item: item["angle"]):
        if not unique or abs(signed_delta(record["angle"], unique[-1]["angle"])) > tolerance:
            unique.append(record)
    return unique


def nearest_clearance_escape(pivot, motion_radius, moving, first, second, clearance,
                             tolerance=geometry.GEOMETRY_TOLERANCE):
    """Find the nearest exact capsule boundary leading outside, without search."""
    theta_old = math.atan2(moving[1]-pivot[1], moving[0]-pivot[0])
    current_distance = geometry.segment_segment_distance(pivot, moving, first, second)
    if current_distance > clearance+tolerance:
        raise ValueError("current radial segment is outside target clearance capsule")
    records = capsule_boundary_angles(
        pivot, motion_radius, first, second, clearance, tolerance)
    drawing_scale = max(1.0, motion_radius, math.dist(pivot, first), math.dist(pivot, second))
    angular_epsilon = 16.0*tolerance*drawing_scale/max(motion_radius, tolerance)
    escapes = []
    for record in records:
        delta = signed_delta(record["angle"], theta_old)
        if abs(delta) <= tolerance:
            continue
        direction = 1.0 if delta > 0 else -1.0
        legal_angle = theta_old + delta + direction*angular_epsilon
        endpoint = (pivot[0]+motion_radius*math.cos(legal_angle),
                    pivot[1]+motion_radius*math.sin(legal_angle))
        distance = geometry.segment_segment_distance(pivot, endpoint, first, second)
        if distance+tolerance >= clearance:
            escapes.append((abs(delta), delta, legal_angle, endpoint, distance, record,
                            direction*angular_epsilon))
    if not escapes:
        raise ValueError("no analytic finite-segment clearance escape exists")
    negative_escapes = [item for item in escapes if item[1] < 0]
    positive_escapes = [item for item in escapes if item[1] > 0]
    if not negative_escapes or not positive_escapes:
        raise ValueError("clearance capsule does not yield two-sided angular escape")
    negative_escape = min(negative_escapes, key=lambda item: abs(item[1]))
    positive_escape = min(positive_escapes, key=lambda item: abs(item[1]))
    _, boundary_delta, legal_angle, endpoint, distance, record, epsilon = min(
        escapes, key=lambda item: (item[0], item[1]))
    return {"theta_old": theta_old, "boundary_angle": theta_old+boundary_delta,
            "theta_new": legal_angle, "new_point": endpoint,
            "boundary_feature": record["feature"], "boundary_point": record["point"],
            "legal_side_offset": epsilon, "resulting_distance": distance,
            "candidate_boundary_count": len(records),
            "forbidden_interval": (normalize(theta_old+negative_escape[1]),
                                   normalize(theta_old+positive_escape[1]))}
