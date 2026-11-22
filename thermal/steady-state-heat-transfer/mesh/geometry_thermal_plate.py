import gmsh

gmsh.initialize()
gmsh.model.add("thermal_plate")

# ============================================================
# Geometry
# ============================================================

L = 1.0       # Length [m]
H = 0.5       # Height [m]

# Corner points
p1 = gmsh.model.geo.addPoint(0, 0, 0, 0.05)
p2 = gmsh.model.geo.addPoint(L, 0, 0, 0.05)
p3 = gmsh.model.geo.addPoint(L, H, 0, 0.05)
p4 = gmsh.model.geo.addPoint(0, H, 0, 0.05)

# Lines
bottom = gmsh.model.geo.addLine(p1, p2)
right  = gmsh.model.geo.addLine(p2, p3)
top    = gmsh.model.geo.addLine(p3, p4)
left   = gmsh.model.geo.addLine(p4, p1)

# Closed boundary
loop = gmsh.model.geo.addCurveLoop([
    bottom,
    right,
    top,
    left
])

# Surface
surface = gmsh.model.geo.addPlaneSurface([loop])

# Synchronize geometry
gmsh.model.geo.synchronize()

# ============================================================
# Physical groups
# ============================================================

# Thermal domain
domain = gmsh.model.addPhysicalGroup(
    2,
    [surface]
)
gmsh.model.setPhysicalName(
    2,
    domain,
    "ThermalDomain"
)

# Left boundary: T = 100 C
left_boundary = gmsh.model.addPhysicalGroup(
    1,
    [left]
)
gmsh.model.setPhysicalName(
    1,
    left_boundary,
    "Temperature_Left"
)

# Right boundary: T = 20 C
right_boundary = gmsh.model.addPhysicalGroup(
    1,
    [right]
)
gmsh.model.setPhysicalName(
    1,
    right_boundary,
    "Temperature_Right"
)

# Top and bottom: insulated
insulated = gmsh.model.addPhysicalGroup(
    1,
    [top, bottom]
)
gmsh.model.setPhysicalName(
    1,
    insulated,
    "Insulated"
)

# ============================================================
# Mesh
# ============================================================

gmsh.option.setNumber(
    "Mesh.CharacteristicLengthMin",
    0.05
)

gmsh.option.setNumber(
    "Mesh.CharacteristicLengthMax",
    0.05
)

gmsh.model.mesh.generate(2)

# Save mesh
from pathlib import Path
# write the .msh to same folder
mesh_file = Path(__file__).with_name("thermal_plate.msh")
gmsh.write(str(mesh_file))
# Open Gmsh GUI
gmsh.fltk.run()

gmsh.finalize()