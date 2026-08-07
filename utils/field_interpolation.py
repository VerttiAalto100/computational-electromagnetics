import numpy as np

from skfem import Basis, ElementTriP1
from skfem.utils import project


def eval_scalar_field_at_points(scalar_field, basis, points):
    """Evaluate the primary DOF field at arbitrary (x,y) points."""
    P = basis.probes(points.T)
    return P @ scalar_field


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


def eval_B_at_points(A_z, basis, points, average_qp=True):
    """
    Evaluate B = curl(Az zhat) at arbitrary points.

    Parameters
    ----------
    A_z : ndarray
        FE solution.
    basis : Basis
        Basis used to solve A_z.
    points : (N,2) ndarray
        Query points.
    average_qp : bool
        Average quadrature values inside each element.
        Recommended for visualization.

    Returns
    -------
    Bx, By, Bnorm
        Arrays of length N.
    """

    uh = basis.interpolate(A_z)

    gradx = uh.grad[0]
    grady = uh.grad[1]

    # P2 -> (nelems,nqp)
    # P1 -> (nelems,)
    if gradx.ndim == 2:
        if average_qp:
            dAdx = gradx.mean(axis=1)
            dAdy = grady.mean(axis=1)
        else:
            dAdx = gradx[:, 0]
            dAdy = grady[:, 0]
    else:
        dAdx = gradx
        dAdy = grady

    Bx_elem = dAdy
    By_elem = -dAdx

    finder = basis.mesh.element_finder()
    elem = finder(points[:, 0], points[:, 1])

    if np.any(elem < 0):
        raise ValueError("Some points lie outside the mesh.")

    Bx = Bx_elem[elem]
    By = By_elem[elem]

    return Bx, By, np.sqrt(Bx**2 + By**2)

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

