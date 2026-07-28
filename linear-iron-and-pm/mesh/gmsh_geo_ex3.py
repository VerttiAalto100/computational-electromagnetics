import gmsh

gmsh.initialize()
gmsh.model.add("c_core_pm_air")

# --------------------------------------------------
# Units
# --------------------------------------------------

cm = 1.0 / 100.0
lc = 0.2 * cm


# --------------------------------------------------
# Dimensions
# --------------------------------------------------

w = 6 * cm
h = w

hy = 2 * cm
wmag = 2 * cm

# Air padding
air_pad = 1 * w


# --------------------------------------------------
# Points for C-core
# --------------------------------------------------

pts = []

coords = [
    (0, 0),
    (w, 0),
    (w, hy),
    (w - wmag, hy),
    (hy, hy),
    (hy, h - hy),
    (w - wmag, h - hy),
    (w, h - hy),
    (w, h),
    (0, h),
]

for x, y in coords:
    pts.append(
        gmsh.model.occ.addPoint(x, y, 0, lc)
    )


# --------------------------------------------------
# C-core boundary
# --------------------------------------------------

lines = []

for i in range(10):
    lines.append(
        gmsh.model.occ.addLine(
            pts[i],
            pts[(i+1)%10]
        )
    )

l0,l1,l2,l3,l4,l5,l6,l7,l8,l9 = lines
lm = gmsh.model.occ.addLine(pts[2], pts[7])


# --------------------------------------------------
# Surfaces
# --------------------------------------------------

core_loop = gmsh.model.occ.addCurveLoop(
    [l0, l1, l2, l3, l4, l5, l6, l7, l8, l9]
)

core_surface = gmsh.model.occ.addPlaneSurface([core_loop])


# Single rectangular magnet
mag_loop = gmsh.model.occ.addCurveLoop([
    l2,     # p2 -> p3
    l3,     # p3 -> p4
    l4,     # p4 -> p5
    l5,     # p5 -> p6
    l6,     # p6 -> p7
    -lm     # p7 -> p2
])

magnet = gmsh.model.occ.addPlaneSurface([mag_loop])


# --------------------------------------------------
# Air box
# --------------------------------------------------

air_box = gmsh.model.occ.addRectangle(
    -air_pad,
    -air_pad,
    0,
    w + 2*air_pad,
    h + 2*air_pad
)


gmsh.model.occ.synchronize()


# --------------------------------------------------
# Fragment air with core/magnets
# --------------------------------------------------

objects = [
    (2, air_box),
    (2, core_surface),
    (2, magnet)
]


gmsh.model.occ.fragment(
    [(2,air_box)],
    [(2,core_surface),
     (2,magnet)]
)


gmsh.model.occ.synchronize()



# --------------------------------------------------
# Identify domains
# --------------------------------------------------

# --------------------------------------------------
# Fragment air with core/magnets
# --------------------------------------------------

result, frag_map = gmsh.model.occ.fragment(
    [(2,air_box)],
    [
        (2,core_surface),
        (2,magnet)
    ]
)

gmsh.model.occ.synchronize()# --------------------------------------------------
# Recover physical regions
# --------------------------------------------------

core_surfaces = []
mag_surfaces = []



# frag_map order follows the input list:
#
# 0 -> air_box
# 1 -> core_surface
# 2 -> magnet


# Core
core_surfaces = [
    e[1] for e in frag_map[1]
]


# North magnet
mag_surfaces = [
    e[1] for e in frag_map[2]
]


# --------------------------------------------------
# Air = all remaining surfaces
# --------------------------------------------------

all_surfaces = [
    s[1] for s in gmsh.model.getEntities(2)
]


solid_surfaces = (
    core_surfaces +
    mag_surfaces
)


air_surfaces = list(
    set(all_surfaces) -
    set(solid_surfaces)
)

# --------------------------------------------------
# Physical groups
# --------------------------------------------------

gmsh.model.addPhysicalGroup(
    2,
    core_surfaces,
    name="core"
)

gmsh.model.addPhysicalGroup(
    2,
    mag_surfaces,
    name="magnet"
)

gmsh.model.addPhysicalGroup(
    2,
    air_surfaces,
    name="air"
)

# --------------------------------------------------
# Find only external air boundary
# --------------------------------------------------
# --------------------------------------------------
# Find outer air rectangle edges geometrically
# --------------------------------------------------

outer_boundary = []

xmin = -air_pad
xmax = w + air_pad
ymin = -air_pad
ymax = h + air_pad

tol = 1e-8   # increase to 1e-6 if needed


for dim, tag in gmsh.model.getEntities(1):

    x, y, z = gmsh.model.occ.getCenterOfMass(
        dim,
        tag
    )

    # vertical left edge
    if abs(x - xmin) < tol:
        outer_boundary.append(tag)

    # vertical right edge
    elif abs(x - xmax) < tol:
        outer_boundary.append(tag)

    # horizontal bottom edge
    elif abs(y - ymin) < tol:
        outer_boundary.append(tag)

    # horizontal top edge
    elif abs(y - ymax) < tol:
        outer_boundary.append(tag)


outer_boundary = list(set(outer_boundary))



gmsh.model.addPhysicalGroup(
    1,
    outer_boundary,
    name="Outer_Air_Boundary"
)

print("Outer boundary edges:", outer_boundary)# --------------------------------------------------
# Mesh
# --------------------------------------------------
# Disable automatic mesh size adaptation
gmsh.option.setNumber("Mesh.MeshSizeFromPoints", 0)
gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", 0)
gmsh.option.setNumber("Mesh.MeshSizeExtendFromBoundary", 0)

# Force a uniform element size
gmsh.option.setNumber("Mesh.MeshSizeMin", lc)
gmsh.option.setNumber("Mesh.MeshSizeMax", lc)
gmsh.model.mesh.generate(2)

from pathlib import Path

# wirte the .msh to same folder
mesh_file = Path(__file__).with_name("mesh_ex3.msh")
gmsh.write(str(mesh_file))

gmsh.fltk.run()

gmsh.finalize()
