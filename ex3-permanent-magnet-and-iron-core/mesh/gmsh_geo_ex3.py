import gmsh

gmsh.initialize()
gmsh.model.add("c_core_pm_air")

# --------------------------------------------------
# Units
# --------------------------------------------------

cm = 1.0 / 100.0
lc = 0.1 * cm


# --------------------------------------------------
# Dimensions
# --------------------------------------------------

w = 6 * cm
h = w

hy = 2 * cm
wmag = 2 * cm

# Air padding
air_pad = 3 * w


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


# Internal magnet separation
lm1 = gmsh.model.occ.addLine(pts[2], pts[7])
lm2 = gmsh.model.occ.addLine(pts[3], pts[6])


# --------------------------------------------------
# Surfaces
# --------------------------------------------------

core_loop = gmsh.model.occ.addCurveLoop(
    [l0,l1,l2,l3,l4,l5,l6,l7,l8,l9]
)

core_surface = gmsh.model.occ.addPlaneSurface(
    [core_loop]
)


mag_loop_south = gmsh.model.occ.addCurveLoop(
    [-l2,lm1,-l6,-lm2]
)

mag_loop_north = gmsh.model.occ.addCurveLoop(
    [l3,l4,l5,-lm2]
)

magnet_S = gmsh.model.occ.addPlaneSurface(
    [mag_loop_south]
)

magnet_N = gmsh.model.occ.addPlaneSurface(
    [mag_loop_north]
)



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
    (2, magnet_N),
    (2, magnet_S)
]


gmsh.model.occ.fragment(
    [(2,air_box)],
    [(2,core_surface),
     (2,magnet_N),
     (2,magnet_S)]
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
        (2,magnet_N),
        (2,magnet_S)
    ]
)

gmsh.model.occ.synchronize()# --------------------------------------------------
# Recover physical regions
# --------------------------------------------------

core_surfaces = []
magN_surfaces = []
magS_surfaces = []


# frag_map order follows the input list:
#
# 0 -> air_box
# 1 -> core_surface
# 2 -> magnet_N
# 3 -> magnet_S


# Core
core_surfaces = [
    e[1] for e in frag_map[1]
]


# North magnet
magN_surfaces = [
    e[1] for e in frag_map[2]
]


# South magnet
magS_surfaces = [
    e[1] for e in frag_map[3]
]


# --------------------------------------------------
# Air = all remaining surfaces
# --------------------------------------------------

all_surfaces = [
    s[1] for s in gmsh.model.getEntities(2)
]


solid_surfaces = (
    core_surfaces +
    magN_surfaces +
    magS_surfaces
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
    magN_surfaces,
    name="magnet_N"
)

gmsh.model.addPhysicalGroup(
    2,
    magS_surfaces,
    name="magnet_S"
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

gmsh.model.mesh.generate(2)

gmsh.write("c_core_pm.msh")


gmsh.fltk.run()

gmsh.finalize()