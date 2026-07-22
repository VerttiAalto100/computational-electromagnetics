import numpy as np


def compare_point_values(
    basis,
    solution,
    points,
    reference,
    quantity_name="Value",
):
    """
    Compare FEM solution against reference values at specified points.

    Parameters
    ----------
    basis : skfem Basis
        FEM basis.
    solution : ndarray
        FEM solution vector.
    points : (N, 2) ndarray
        Evaluation points [[x1, y1], ...].
    reference : ndarray
        Reference values (e.g. COMSOL or analytical solution).
    quantity_name : str
        Name of the compared quantity.

    Returns
    -------
    values : ndarray
        FEM values.
    abs_error : ndarray
        Absolute error.
    rel_error : ndarray
        Relative error in percent.
    """

    values = np.empty(len(points))

    for i, (x, y) in enumerate(points):
        P = basis.probes(np.array([[x], [y]]))
        values[i] = (P @ solution)[0]

    abs_error = np.abs(values - reference)
    rel_error = 100 * abs_error / np.abs(reference)

    print(
        f"{'x':>7} {'y':>7}"
        f"{'FEM':>15}"
        f"{'Reference':>15}"
        f"{'Abs. Error':>15}"
        f"{'Rel. Error %':>15}"
    )
    print("-" * 80)

    for i, (x, y) in enumerate(points):
        print(
            f"{x:7.3f}"
            f"{y:7.3f}"
            f"{values[i]:15.6e}"
            f"{reference[i]:15.6e}"
            f"{abs_error[i]:15.6e}"
            f"{rel_error[i]:15.3f}"
        )

    print("\nSummary")
    print("-------")
    print(f"{quantity_name} mean relative error : {rel_error.mean():.3f} %")
    print(f"{quantity_name} max relative error  : {rel_error.max():.3f} %")
    print(
        f"{quantity_name} RMS relative error  : "
        f"{np.sqrt(np.mean(rel_error**2)):.3f} %"
    )

    return values, abs_error, rel_error