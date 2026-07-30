import gmsh

gmsh.initialize()

gmsh.model.add("induction_motor")

# Import STEP file
gmsh.model.occ.importShapes("induction-motor/mesh/compumag-induction-motor.step")

# Synchronize OCC kernel with Gmsh model
gmsh.model.occ.synchronize()

# ------------------------------------
# Print all surfaces
# ------------------------------------

# Print all surfaces
surfaces = gmsh.model.getEntities(dim=2)

print("\nSurfaces:")
for dim, tag in surfaces:
    print(f"Surface tag = {tag}")

# Open gmsh and see the surfaces in the GUI, assign correct physical groups to the surfaces.

# Air gap
air = gmsh.model.addPhysicalGroup(2, [1])
gmsh.model.setPhysicalName(2, air, "AirGap")

# Rotor
rotor = gmsh.model.addPhysicalGroup(2, [2])
gmsh.model.setPhysicalName(2, rotor, "Rotor")

# Stator
stator = gmsh.model.addPhysicalGroup(2, [3])
gmsh.model.setPhysicalName(2, stator, "Stator")

# Coil A+
Aplus = gmsh.model.addPhysicalGroup(2, [4, 5, 6, 7, 8, 9])
gmsh.model.setPhysicalName(2, Aplus, "A_plus")

# Coil A-
Aminus = gmsh.model.addPhysicalGroup(2, [10, 11, 12, 13, 14, 15])
gmsh.model.setPhysicalName(2, Aminus, "A_minus")

gmsh.fltk.run()