import numpy as np

from skfem import Basis, ElementTriP1
from skfem.utils import project


def eval_Az_at_points(A_z, basis, points):
    """Evaluate the primary DOF field A_z at arbitrary (x,y) points."""
    P = basis.probes(points.T)
    return P @ A_z


def elem_to_nodal(mesh, elem_values):
    """
    Recover a continuous nodal (P1) field from elementwise-constant values
    by simple averaging over all elements touching each vertex.
    """
    t = mesh.t          # shape (3, n_elems) - vertex indices per triangle
    n_vertices = mesh.p.shape[1]

    nodal_sum = np.zeros(n_vertices)
    nodal_count = np.zeros(n_vertices)

    for local_node in range(3):
        vertex_ids = t[local_node]
        nodal_sum[vertex_ids] += elem_values
        nodal_count[vertex_ids] += 1

    return nodal_sum / nodal_count


def eval_B_at_points(A_z, basis, points, smooth=True):
    """
    Evaluate B = curl(A_z ẑ) = (dAz/dy, -dAz/dx) at arbitrary (x,y) points.
    """
    uh = basis.interpolate(A_z)
    Bx_elem = uh.grad[1].mean(axis=1)     # (n_elems,)
    By_elem = -uh.grad[0].mean(axis=1)    # (n_elems,)

    if smooth:
        mesh = basis.mesh
        Bx_nodal = elem_to_nodal(mesh, Bx_elem)
        By_nodal = elem_to_nodal(mesh, By_elem)

        basis_p1 = Basis(mesh, ElementTriP1())
        P = basis_p1.probes(points.T)
        Bx_pts = P @ Bx_nodal
        By_pts = P @ By_nodal
    else:
        finder = basis.mesh.element_finder()
        tris = finder(points[:, 0], points[:, 1])
        Bx_pts = Bx_elem[tris]
        By_pts = By_elem[tris]

    Bnorm_pts = np.sqrt(Bx_pts**2 + By_pts**2)
    return Bx_pts, By_pts, Bnorm_pts

def eval_Br_Bt_at_points(A_z, basis, points):
    
    #Evaluate Br, Bt (polar components of B = curl(A_z ẑ)) at arbitrary
    #(x,y) points, using elementwise-constant gradient + nearest-element lookup.
    
    uh = basis.interpolate(A_z)
    Bx_elem = uh.grad[1].mean(axis=1)     # dAz/dy per element
    By_elem = -uh.grad[0].mean(axis=1)    # -dAz/dx per element

    finder = basis.mesh.element_finder()
    tris = finder(points[:, 0], points[:, 1])

    Bx_pts = Bx_elem[tris]
    By_pts = By_elem[tris]

    x, y = points[:, 0], points[:, 1]
    r = np.sqrt(x**2 + y**2) + 1e-15

    Br_pts = (Bx_pts*x + By_pts*y) / r
    Bt_pts = (-Bx_pts*y + By_pts*x) / r

    return Br_pts, Bt_pts

