"""
2D PMSM (surface-mounted magnet) cross-section geometry, built with the
Gmsh Python API (OpenCASCADE kernel).

Thermal-study version: no surrounding air box and no Maxwell stress
tensor integration circle (those were only needed for the EM/torque
study). The stator's outer radius is now the true outer edge of the
model.

Regions produced (as Gmsh Physical Groups on 2D surfaces):
    RotorYoke
    Magnet_N, Magnet_S                           (2 magnets, 1 pole pair)
    AirGap                                       (single ring, no longer split)
    StatorYoke                                   (core ring, 6 bolt holes cut out)
    Slot_1 ... Slot_6
Physical curves:
    OuterBoundary                                (outer edge of the stator, for
                                                   convection/ambient BCs)

Edit the PARAMETERS block below to match your machine.
"""

import gmsh
import math

# ----------------------------------------------------------------------
# PARAMETERS (m)
# ----------------------------------------------------------------------
mm = 1e-3  # nastran needs to be x1000 to match comsol. This is done to get correct mesh


r_shaft      = 12.0*mm    # shaft radius
r_rotor      = 32.0*mm   # rotor / magnet outer radius
r_stator_in  = 34.0*mm    # stator bore radius  -> airgap = r_stator_in - r_rotor
r_stator_out = 65.0*mm    # stator outer radius
r_bolt_pcd   = 50.0*mm    # bolt hole pitch circle radius
r_bolt       = 4.0*mm     # bolt hole radius
n_bolts      = 6
n_magnets    = 1       # 1 pole pair

gmsh.initialize()
gmsh.model.add("pmsm_2d")
occ = gmsh.model.occ

# ----------------------------------------------------------------------
# 1) Basic shapes: nested full disks + bolt disks + spokes
# ----------------------------------------------------------------------
shaft_disk     = occ.addDisk(0, 0, 0, r_shaft,      r_shaft)
rotor_disk     = occ.addDisk(0, 0, 0, r_rotor,      r_rotor)
statorin_disk  = occ.addDisk(0, 0, 0, r_stator_in,  r_stator_in)
statorout_disk = occ.addDisk(0, 0, 0, r_stator_out, r_stator_out)

# 6 bolt holes (evenly spaced, first one on the +x axis, like the reference picture)
bolt_disks, bolt_centers = [], []
for i in range(n_bolts):
    ang = 2 * math.pi * i / n_bolts
    bx, by = r_bolt_pcd * math.cos(ang), r_bolt_pcd * math.sin(ang)
    bolt_centers.append((bx, by))
    bolt_disks.append(occ.addDisk(bx, by, 0, r_bolt, r_bolt))

# 2 short radial spokes (0/180 deg) that split the magnet ring into quadrants
spoke_lines = []
#
#for k in range(n_magnets):
#    ang = 2 * math.pi * k / n_magnets
#    p1 = occ.addPoint(r_shaft * math.cos(ang), r_shaft * math.sin(ang), 0)
#    p2 = occ.addPoint(r_rotor * math.cos(ang), r_rotor * math.sin(ang), 0)
#    spoke_lines.append(occ.addLine(p1, p2))

# ----------------------------------------------------------------------
# 2) Fragment everything together -> conformal, non-overlapping regions
# ----------------------------------------------------------------------
surfaces = [(2, statorout_disk), (2, statorin_disk),
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
tol = max(r_stator_out * 1e-4, 1e-6)
groups = {"RotorYoke": [], "AirGap": [], "StatorYoke": [], "Magnet": []}
for i in range(n_bolts):
    groups[f"Slot_{i + 1}"] = []

# midpoint thresholds between successive radii
t_shaft_rotor  = (r_shaft + r_rotor) / 2
t_rotor_sin    = (r_rotor + r_stator_in) / 2
t_sin_sout     = (r_stator_in + r_stator_out) / 2

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
        groups["Magnet"].append(tag)
    elif extent < t_sin_sout:
        groups["AirGap"].append(tag)
    else:
        groups["StatorYoke"].append(tag)

for name, tags in groups.items():
    if tags:
        gmsh.model.addPhysicalGroup(2, tags, name=name)

# ----------------------------------------------------------------------
# Outer boundary of the whole model = outer edge of the stator (no air
# box in the thermal study, so this is the true domain boundary; use it
# for a convection / ambient-temperature BC).
# ----------------------------------------------------------------------
outer_edges = []
for dim, tag in gmsh.model.getEntities(1):
    xmin, ymin, _, xmax, ymax, _ = gmsh.model.getBoundingBox(dim, tag)
    extent = max(abs(xmin), abs(xmax), abs(ymin), abs(ymax))
    if abs(extent - r_stator_out) < tol:
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
lc_fine   = 0.6*mm   # element size in the airgap / around bolt holes
lc_coarse = 4.0*mm   # element size far from those features
dist_min  = 1.0*mm   # stay at lc_fine within this distance of the curves
dist_max  = 9.0*mm   # grow linearly out to lc_coarse by this distance

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
# write the .msh to same folder
mesh_file = Path(__file__).with_name("mesh_pmsm_2d.bdf")


gmsh.write(str(mesh_file))

gmsh.fltk.run()   # uncomment to inspect interactively (requires a display)
gmsh.finalize()