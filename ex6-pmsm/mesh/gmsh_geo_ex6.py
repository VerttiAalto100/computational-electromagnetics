"""
2D PMSM (surface-mounted magnet) cross-section geometry, built with the
Gmsh Python API (OpenCASCADE kernel).

Regions produced (as Gmsh Physical Groups on 2D surfaces):
    Shaft
    Magnet_N1, Magnet_S1, Magnet_N2, Magnet_S2   (4 magnets, 2 pole pairs)
    AirGap
    Stator                                       (core ring, 6 bolt holes cut out)
    BoltHole_1 ... BoltHole_6
    AirBox                                       (surrounding air, outside the stator)
Physical curve:
    OuterBoundary                                (outer edge of the air box, for BCs)

Edit the PARAMETERS block below to match your machine.
"""

import gmsh
import math

# ----------------------------------------------------------------------
# PARAMETERS (mm)
# ----------------------------------------------------------------------
r_shaft      = 12.0    # shaft radius
r_rotor      = 32.0    # rotor / magnet outer radius
r_stator_in  = 34.0    # stator bore radius  -> airgap = r_stator_in - r_rotor
r_stator_out = 65.0    # stator outer radius
r_bolt_pcd   = 50.0    # bolt hole pitch circle radius
r_bolt       = 4.0     # bolt hole radius
n_bolts      = 6
n_magnets    = 4       # 2 pole pairs

box_half     = 100.0   # surrounding square air box, half side length

gmsh.initialize()
gmsh.model.add("pmsm_2d")
occ = gmsh.model.occ

# ----------------------------------------------------------------------
# 1) Basic shapes: nested full disks + outer square + bolt disks + spokes
# ----------------------------------------------------------------------
shaft_disk     = occ.addDisk(0, 0, 0, r_shaft,      r_shaft)
rotor_disk     = occ.addDisk(0, 0, 0, r_rotor,      r_rotor)
statorin_disk  = occ.addDisk(0, 0, 0, r_stator_in,  r_stator_in)
statorout_disk = occ.addDisk(0, 0, 0, r_stator_out, r_stator_out)
air_box        = occ.addRectangle(-box_half, -box_half, 0, 2 * box_half, 2 * box_half)

# 6 bolt holes (evenly spaced, first one on the +x axis, like the reference picture)
bolt_disks, bolt_centers = [], []
for i in range(n_bolts):
    ang = 2 * math.pi * i / n_bolts
    bx, by = r_bolt_pcd * math.cos(ang), r_bolt_pcd * math.sin(ang)
    bolt_centers.append((bx, by))
    bolt_disks.append(occ.addDisk(bx, by, 0, r_bolt, r_bolt))

# 4 short radial spokes (0/90/180/270 deg) that split the magnet ring into quadrants
spoke_lines = []
for k in range(n_magnets):
    ang = 2 * math.pi * k / n_magnets
    p1 = occ.addPoint(r_shaft * math.cos(ang), r_shaft * math.sin(ang), 0)
    p2 = occ.addPoint(r_rotor * math.cos(ang), r_rotor * math.sin(ang), 0)
    spoke_lines.append(occ.addLine(p1, p2))

# ----------------------------------------------------------------------
# 2) Fragment everything together -> conformal, non-overlapping regions
# ----------------------------------------------------------------------
surfaces = [(2, air_box), (2, statorout_disk), (2, statorin_disk),
            (2, rotor_disk), (2, shaft_disk)] + [(2, d) for d in bolt_disks]
curves = [(1, l) for l in spoke_lines]

occ.fragment(surfaces, curves)
occ.synchronize()

# ----------------------------------------------------------------------
# 3) Classify every resulting surface
#
# NOTE: a full, rotationally-symmetric ring (Shaft disk, AirGap, Stator)
# has its centre of mass exactly at the origin no matter its radius, so
# centre-of-mass distance can't tell those apart. Bounding-box extent
# (how far the surface reaches from the origin) is used instead - it
# scales with the ring's outer radius the way we need. Bolt holes (which
# are NOT symmetric about the origin) are still identified by proximity
# of their centroid to a known bolt centre.
# ----------------------------------------------------------------------
tol = 1e-6
groups = {"RotorYoke": [], "AirGap": [], "StatorYoke": [], "AirBox": [],
          "Magnet_N1": [], "Magnet_S1": [], "Magnet_N2": [], "Magnet_S2": []}
for i in range(n_bolts):
    groups[f"Slot_{i + 1}"] = []

# midpoint thresholds between successive radii
t_shaft_rotor  = (r_shaft + r_rotor) / 2
t_rotor_sin    = (r_rotor + r_stator_in) / 2
t_sin_sout     = (r_stator_in + r_stator_out) / 2
t_sout_box     = (r_stator_out + box_half) / 2

for dim, tag in gmsh.model.getEntities(2):
    x, y, _ = occ.getCenterOfMass(dim, tag)

    # is this one of the small bolt-hole disks? (checked first: bolt holes
    # sit at the same radial band as the stator, so radius alone won't do)
    bolt_hit = False
    for i, (bx, by) in enumerate(bolt_centers):
        if math.hypot(x - bx, y - by) < r_bolt * 0.9:
            groups[f"Slot_{i + 1}"].append(tag)
            bolt_hit = True
            break
    if bolt_hit:
        continue

    xmin, ymin, _, xmax, ymax, _ = gmsh.model.getBoundingBox(dim, tag)
    extent = max(abs(xmin), abs(xmax), abs(ymin), abs(ymax))

    if extent < t_shaft_rotor:
        groups["RotorYoke"].append(tag)
    elif extent < t_rotor_sin:
        ang = math.atan2(y, x) % (2 * math.pi)
        quadrant = int(ang // (math.pi / 2)) % 4          # 0..3
        name = ["Magnet_N1", "Magnet_S1", "Magnet_N2", "Magnet_S2"][quadrant]
        groups[name].append(tag)
    elif extent < t_sin_sout:
        groups["AirGap"].append(tag)
    elif extent < t_sout_box:
        groups["StatorYoke"].append(tag)
    else:
        groups["AirBox"].append(tag)

for name, tags in groups.items():
    if tags:
        gmsh.model.addPhysicalGroup(2, tags, name=name)

# outer boundary of the air box (the 4 straight edges of the square), for BCs
edges = gmsh.model.getBoundary([(2, t) for t in groups["AirBox"]], combined=False, oriented=False)
outer_edges = []
for dim, tag in edges:
    xmin, ymin, _, xmax, ymax, _ = gmsh.model.getBoundingBox(dim, tag)
    if (abs(xmin + box_half) < tol or abs(xmax - box_half) < tol or
            abs(ymin + box_half) < tol or abs(ymax - box_half) < tol):
        outer_edges.append(tag)
gmsh.model.addPhysicalGroup(1, outer_edges, name="OuterBoundary")

# ----------------------------------------------------------------------
# 4) Mesh size field: refine near the airgap and the bolt holes
#
# A "Distance" field measures, for every point in the domain, how close
# it is to a chosen set of curves. A "Threshold" field turns that
# distance into an element size: SizeMin right on the curves, ramping up
# to SizeMax once you're DistMax away. Feeding this in as the background
# mesh (instead of a single global size) lets the airgap/bolt holes be
# fine while the rest of the model stays coarse and fast to mesh.
# ----------------------------------------------------------------------
lc_fine   = 0.6   # element size in the airgap / around bolt holes (mm)
lc_coarse = 4.0   # element size far from those features (mm)
dist_min  = 1.0   # stay at lc_fine within this distance of the curves (mm)
dist_max  = 9.0   # grow linearly out to lc_coarse by this distance (mm)

# disable the size heuristics so only our field controls element size
gmsh.option.setNumber("Mesh.MeshSizeExtendFromBoundary", 0)
gmsh.option.setNumber("Mesh.MeshSizeFromPoints", 0)
gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", 0)

refine_curves = gmsh.model.getBoundary(
    [(2, t) for t in groups["AirGap"]] +
    [(2, t) for i in range(n_bolts) for t in groups[f"Slot_{i + 1}"]],
    oriented=False, recursive=False)
refine_curve_tags = [tag for dim, tag in refine_curves if dim == 1]

dist_field = gmsh.model.mesh.field.add("Distance")
gmsh.model.mesh.field.setNumbers(dist_field, "CurvesList", refine_curve_tags)
gmsh.model.mesh.field.setNumber(dist_field, "Sampling", 100)

thresh_field = gmsh.model.mesh.field.add("Threshold")
gmsh.model.mesh.field.setNumber(thresh_field, "InField", dist_field)
gmsh.model.mesh.field.setNumber(thresh_field, "SizeMin", lc_fine)
gmsh.model.mesh.field.setNumber(thresh_field, "SizeMax", lc_coarse)
gmsh.model.mesh.field.setNumber(thresh_field, "DistMin", dist_min)
gmsh.model.mesh.field.setNumber(thresh_field, "DistMax", dist_max)

gmsh.model.mesh.field.setAsBackgroundMesh(thresh_field)

# ----------------------------------------------------------------------
# 5) Mesh & save
# ----------------------------------------------------------------------
gmsh.model.mesh.generate(2)

from pathlib import Path
# wirte the .msh to same folder
mesh_file = Path(__file__).with_name("mesh_ex6.msh")

gmsh.write(str(mesh_file))

gmsh.fltk.run()   # uncomment to inspect interactively (requires a display)
gmsh.finalize()
