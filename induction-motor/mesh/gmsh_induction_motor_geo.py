import gmsh
import math

gmsh.initialize()
gmsh.model.add("induction_motor")

Rotor_Yoke_tag      = 1
Rotor_Aluminium_tag = 2
Airgap_tags         = [3, 17]
Air_Holes_tags      = [4, 5, 6, 7, 8, 9]
Slots_tags          = [10, 11, 12, 13, 14, 15]
Stator_Yoke_tag     = 16

gmsh.model.occ.importShapes("induction-motor/mesh/compumag-induction-motor.step")
gmsh.model.occ.synchronize()

xmin, ymin, zmin, xmax, ymax, zmax = gmsh.model.getBoundingBox(2, Stator_Yoke_tag)
cx = 0.5 * (xmin + xmax)
cy = 0.5 * (ymin + ymax)
stator_outer_r = 0.5 * max(xmax - xmin, ymax - ymin)

AIRBOX_FACTOR = 5
airbox_r = AIRBOX_FACTOR * stator_outer_r

airbox_disk_tag = gmsh.model.occ.addDisk(cx, cy, 0, airbox_r, airbox_r)
gmsh.model.occ.synchronize()

name_by_old_tag = {Rotor_Yoke_tag: "Rotor_Yoke", Rotor_Aluminium_tag: "Rotor_Aluminium",
                   Stator_Yoke_tag: "Stator_Yoke"}
for t in Airgap_tags:
    name_by_old_tag[t] = "AirGap"
for i, t in enumerate(Air_Holes_tags):
    name_by_old_tag[t] = f"Air_Hole"
for i, t in enumerate(Slots_tags):
    name_by_old_tag[t] = f"Slot_{i+1}"

motor_surface_tags = list(name_by_old_tag.keys())

pre_frag_info = {}
for t in motor_surface_tags:
    x, y, _ = gmsh.model.occ.getCenterOfMass(2, t)
    area = gmsh.model.occ.getMass(2, t)
    pre_frag_info[t] = (x, y, area)

all_shapes = [(2, airbox_disk_tag)] + [(2, t) for t in motor_surface_tags]
gmsh.model.occ.fragment(all_shapes, [])
gmsh.model.occ.synchronize()

AREA_REL_TOL = 1e-6
LEN_TOL = stator_outer_r * 1e-6

def match_original(x, y, area):
    for old_tag, (ox, oy, oarea) in pre_frag_info.items():
        if oarea > 0 and abs(area - oarea) / oarea < AREA_REL_TOL:
            if math.hypot(x - ox, y - oy) < LEN_TOL:
                return old_tag
    return None

groups = {name: [] for name in name_by_old_tag.values()}
groups["Air_Box"] = []

for dim, tag in gmsh.model.getEntities(2):
    x, y, _ = gmsh.model.occ.getCenterOfMass(2, tag)
    area = gmsh.model.occ.getMass(2, tag)
    old_tag = match_original(x, y, area)
    if old_tag is not None:
        groups[name_by_old_tag[old_tag]].append(tag)
    else:
        groups["Air_Box"].append(tag)

Air_Box_tags = groups["Air_Box"]

boundary_curves = gmsh.model.getBoundary(
    [(2, t) for t in Air_Box_tags], oriented=False, combined=False
)

outer_curves = []
for dim, tag in boundary_curves:
    cxmin, cymin, _, cxmax, cymax, _ = gmsh.model.getBoundingBox(dim, tag)
    r_est = 0.5 * max(cxmax - cxmin, cymax - cymin)
    if abs(r_est - airbox_r) < 1e-3 * airbox_r:
        outer_curves.append(tag)
outer_curves = list(set(outer_curves))

gmsh.model.occ.synchronize()

airgap_curves = []
for s in groups["AirGap"]:
    b = gmsh.model.getBoundary([(2, s)], oriented=False, combined=False)
    airgap_curves.extend(tag for dim, tag in b)
airgap_curves = list(set(airgap_curves))

# ------------------------------------------------------------------
# Scale geometry to meters
# ------------------------------------------------------------------
gmsh.model.occ.dilate(
    gmsh.model.occ.getEntities(),
    0, 0, 0,
    0.001, 0.001, 0.001
)
gmsh.model.occ.synchronize()

# ------------------------------------------------------------------
# NOW add physical groups -- ONLY ONCE, after everything is finalized
# ------------------------------------------------------------------
for name, tags in groups.items():
    if tags:
        gmsh.model.addPhysicalGroup(2, tags, name=name)

gmsh.model.addPhysicalGroup(1, outer_curves, name="Outer_Boundary")
gmsh.model.addPhysicalGroup(1, [34], name="MST_Curve")


# ------------------------------------------------------------------
# 6. Mesh sizing: fine in the air gap, coarse at the air box's outer edge
# ------------------------------------------------------------------
mm = 1e-3   # geometry is already in meters after the dilate() call above

# --- Measure the real airgap width instead of assuming it ---
# Airgap surfaces are thin annular rings; their bounding box's
# radial extent (max - min, along either axis, whichever is smaller
# since it's a thin ring) approximates the physical gap width.
airgap_xmin, airgap_ymin = float("inf"), float("inf")
airgap_xmax, airgap_ymax = float("-inf"), float("-inf")
for s in groups["AirGap"]:
    x0, y0, _, x1, y1, _ = gmsh.model.getBoundingBox(2, s)
    airgap_xmin, airgap_ymin = min(airgap_xmin, x0), min(airgap_ymin, y0)
    airgap_xmax, airgap_ymax = max(airgap_xmax, x1), max(airgap_ymax, y1)

# For an annulus split into pieces, use rotor OD / stator ID directly if
# you have them as variables -- more reliable than bounding boxes on
# fragmented slivers. Otherwise this is a reasonable estimate:
airgap_outer_r = 0.5 * max(airgap_xmax - airgap_xmin, airgap_ymax - airgap_ymin)
# If you already have rotor_outer_r / stator_inner_r variables from
# earlier in your script, prefer: airgap_width = stator_inner_r - rotor_outer_r
airgap_width = 2 * mm   # <-- replace with a real measured value if possible
print(f"Using airgap_width = {airgap_width/mm:.3f} mm")


lc_fine   = airgap_width / 5     # ~5 elements across the gap radially
lc_coarse = 8 * mm                 # bulk stator/airbox size -- TUNE vs slot/tooth size
dist_min  = 1 * airgap_width     # covers the full radial span of the gap
dist_max  = 2.0 * airgap_width     # grade out into rotor/stator over a short distance

print(f"lc_fine={lc_fine/mm:.4f}mm  lc_coarse={lc_coarse/mm:.2f}mm  "
      f"dist_min={dist_min/mm:.3f}mm  dist_max={dist_max/mm:.3f}mm")

# disable the size heuristics so only our field controls element size
gmsh.option.setNumber("Mesh.MeshSizeExtendFromBoundary", 0)
gmsh.option.setNumber("Mesh.MeshSizeFromPoints", 0)
gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", 0)
gmsh.option.setNumber("Mesh.Algorithm", 6)  # Frontal-Delaunay, plays well with size fields

refine_curves = gmsh.model.getBoundary(
    [(2, t) for t in groups["AirGap"]],
    oriented=False, recursive=False)
refine_curve_tags = [tag for dim, tag in refine_curves if dim == 1]

dist_field = gmsh.model.mesh.field.add("Distance")
gmsh.model.mesh.field.setNumbers(dist_field, "CurvesList", refine_curve_tags)
gmsh.model.mesh.field.setNumber(dist_field, "Sampling", 50)

thresh_field = gmsh.model.mesh.field.add("Threshold")
gmsh.model.mesh.field.setNumber(thresh_field, "InField", dist_field)
gmsh.model.mesh.field.setNumber(thresh_field, "SizeMin", lc_fine)
gmsh.model.mesh.field.setNumber(thresh_field, "SizeMax", lc_coarse)
gmsh.model.mesh.field.setNumber(thresh_field, "DistMin", dist_min)
gmsh.model.mesh.field.setNumber(thresh_field, "DistMax", dist_max)

gmsh.model.mesh.field.setAsBackgroundMesh(thresh_field)



for dim, tag in gmsh.model.getPhysicalGroups():
    name = gmsh.model.getPhysicalName(dim, tag)
    ents = gmsh.model.getEntitiesForPhysicalGroup(dim, tag)
    print(f"Physical group '{name}' (dim={dim}, tag={tag}): entities={ents}")
# ------------------------------------------------------------------
# 7. Generate & save the mesh
# ------------------------------------------------------------------

# Curved mesh elements (2nd order) are needed for accurate FEM results in the airgap.
gmsh.option.setNumber("Mesh.ElementOrder", 2)
gmsh.option.setNumber("Mesh.SecondOrderLinear", 0)
gmsh.option.setNumber("Mesh.HighOrderOptimize", 2)

gmsh.model.mesh.generate(2)

from pathlib import Path
# write the .msh to same folder
mesh_file = Path(__file__).with_name("induction_motor.bdf")


gmsh.write(str(mesh_file))

gmsh.fltk.run()
gmsh.finalize()