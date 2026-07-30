import gmsh
import math

gmsh.initialize()
gmsh.model.add("induction_motor")

# ------------------------------------------------------------------
# 1. Import geometry
# ------------------------------------------------------------------
gmsh.model.occ.importShapes("induction-motor/mesh/compumag-induction-motor.step")
gmsh.model.occ.synchronize()

# Original surface tags coming straight out of the STEP file
Rotor_Yoke_tag      = 1
Rotor_Aluminium_tag = 2
Airgap_tags         = [3, 17]
Air_Holes_tags      = [4, 5, 6, 7, 8, 9]
Slots_tags          = [10, 11, 12, 13, 14, 15]
Stator_Yoke_tag     = 16   # the surface that borders the outside air

# Map every original tag -> the physical-group name it belongs to
name_by_old_tag = {Rotor_Yoke_tag: "Rotor_Yoke", Rotor_Aluminium_tag: "Rotor_Aluminium",
                   Stator_Yoke_tag: "Stator_Yoke"}
for t in Airgap_tags:
    name_by_old_tag[t] = "AirGap"
for i, t in enumerate(Air_Holes_tags):
    name_by_old_tag[t] = f"Air_Hole_{i+1}"
for i, t in enumerate(Slots_tags):
    name_by_old_tag[t] = f"Slot_{i+1}"

motor_surface_tags = list(name_by_old_tag.keys())

print("\nSurfaces before airbox:")
for dim, tag in gmsh.model.getEntities(dim=2):
    print(f"  Surface tag = {tag}")

# Record area + centroid of every original surface BEFORE the boolean op.
# We'll use these to recognise each surface again afterwards, since a
# boolean fragment can renumber tags even for pieces whose shape didn't
# actually change (see step 3 below).
pre_frag_info = {}
for t in motor_surface_tags:
    x, y, _ = gmsh.model.occ.getCenterOfMass(2, t)
    area = gmsh.model.occ.getMass(2, t)
    pre_frag_info[t] = (x, y, area)

# ------------------------------------------------------------------
# 2. Build a circular air box around the motor
# ------------------------------------------------------------------
xmin, ymin, zmin, xmax, ymax, zmax = gmsh.model.getBoundingBox(2, Stator_Yoke_tag)
cx = 0.5 * (xmin + xmax)
cy = 0.5 * (ymin + ymax)
stator_outer_r = 0.5 * max(xmax - xmin, ymax - ymin)

AIRBOX_FACTOR = 1.5          # air box outer radius = 1.5 x stator outer radius -- TUNE
airbox_r = AIRBOX_FACTOR * stator_outer_r

airbox_disk_tag = gmsh.model.occ.addDisk(cx, cy, 0, airbox_r, airbox_r)
gmsh.model.occ.synchronize()

# ------------------------------------------------------------------
# 3. Fragment everything together, then re-classify by geometry
# ------------------------------------------------------------------
# The disk is solid from the centre out to airbox_r, so it overlaps
# every motor surface, not just the stator. Fragment all of them
# together so the whole assembly becomes conformal (shared boundaries
# everywhere), same technique as the working PMSM script.
all_shapes = [(2, airbox_disk_tag)] + [(2, t) for t in motor_surface_tags]
gmsh.model.occ.fragment(all_shapes, [])
gmsh.model.occ.synchronize()

# Do NOT trust out_map index bookkeeping here -- with real (imported)
# B-rep geometry a surface can legitimately come back as more than one
# piece, and mesh.generate() meshes every 2D entity that exists in the
# model regardless of physical groups. Any surface silently dropped by
# index-based tracking still gets meshed and shows up as an overlapping
# "phantom" -- which is what caused the duplicate mesh last time.
#
# Instead, classify every surface that exists AFTER the fragment by
# matching it back to the pre-fragment (area, centroid) we recorded.
# Area alone would be ambiguous for congruent slots (same area,
# different position); centroid alone is ambiguous for concentric
# rings (same centre, different radius). Together they're unique.
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
        groups["Air_Box"].append(tag)   # doesn't match any original region -> new air

print("\nClassification after fragment:")
for name, tags in groups.items():
    print(f"  {name:16s}: {tags}")

total_in = len(gmsh.model.getEntities(2))
total_classified = sum(len(v) for v in groups.values())
print(f"\nTotal surfaces: {total_in}, classified: {total_classified}")
if total_in != total_classified:
    print("  WARNING: mismatch -- some surface was double counted or missed.")
if not groups["Air_Box"]:
    print("  WARNING: no Air_Box surface found -- check AIRBOX_FACTOR / geometry.")

# ------------------------------------------------------------------
# 4. Physical groups (surfaces)
# ------------------------------------------------------------------
for name, tags in groups.items():
    if tags:
        gmsh.model.addPhysicalGroup(2, tags, name=name)

Air_Box_tags = groups["Air_Box"]

# ------------------------------------------------------------------
# 5. Physical group for the exterior boundary of the air box (curves)
# ------------------------------------------------------------------
# The Air_Box surface has (at least) two boundary curves: the brand-new
# outer circle (radius = airbox_r) and the curve(s) shared with the
# stator (radius = stator_outer_r). We only want the outer one.
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
gmsh.model.addPhysicalGroup(1, outer_curves, name="Outer_Boundary")

# ------------------------------------------------------------------
# 6. Mesh sizing: fine in the air gap, coarse at the air box's outer edge
# ------------------------------------------------------------------
airgap_curves = []
for s in groups["AirGap"]:
    b = gmsh.model.getBoundary([(2, s)], oriented=False, combined=False)
    airgap_curves.extend(tag for dim, tag in b)
airgap_curves = list(set(airgap_curves))

AIRGAP_SIZE  = stator_outer_r / 100.0     # fine size in/around the air gap -- TUNE
COARSE_SIZE  = stator_outer_r / 10.0       # coarse size at the air box outer edge -- TUNE
GRADING_DIST = airbox_r - stator_outer_r  # distance over which size grows -- TUNE

gmsh.model.mesh.field.add("Distance", 1)
gmsh.model.mesh.field.setNumbers(1, "CurvesList", airgap_curves)
gmsh.model.mesh.field.setNumber(1, "Sampling", 200)

gmsh.model.mesh.field.add("Threshold", 2)
gmsh.model.mesh.field.setNumber(2, "InField", 1)
gmsh.model.mesh.field.setNumber(2, "SizeMin", AIRGAP_SIZE)
gmsh.model.mesh.field.setNumber(2, "SizeMax", COARSE_SIZE)
gmsh.model.mesh.field.setNumber(2, "DistMin", 0)
gmsh.model.mesh.field.setNumber(2, "DistMax", GRADING_DIST)

gmsh.model.mesh.field.setAsBackgroundMesh(2)

# Let the field be the only thing controlling element size
gmsh.option.setNumber("Mesh.MeshSizeExtendFromBoundary", 0)
gmsh.option.setNumber("Mesh.MeshSizeFromPoints", 0)
gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", 0)
gmsh.option.setNumber("Mesh.Algorithm", 6)  # Frontal-Delaunay, plays well with size fields

# ------------------------------------------------------------------
# 7. Generate & save the mesh
# ------------------------------------------------------------------
gmsh.model.mesh.generate(2)

from pathlib import Path
# write the .msh to same folder
mesh_file = Path(__file__).with_name("induction_motor.msh")

gmsh.write(str(mesh_file))

gmsh.fltk.run()
gmsh.finalize()