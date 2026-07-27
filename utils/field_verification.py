import numpy as np


def compare_point_values(
    own_values,
    ref_values,
    points,
    quantity_name="field",
    n_show=20,
    sort_by="abs_error",
    plot=True,
):
    """
    Compare two arrays of already-evaluated point values (any scalar field:
    A_z, Bnorm, Bx, etc.) against each other at the same coordinates.

    Does NOT know about basis, DOFs, or FE machinery — pure post-processing.
    """
    own_values = np.asarray(own_values)
    ref_values = np.asarray(ref_values)

    abs_error = own_values - ref_values
    ref_scale = np.maximum(np.abs(ref_values), 1e-12 * np.max(np.abs(ref_values)))
    rel_error = np.abs(abs_error) / ref_scale

    L2_abs = np.sqrt(np.mean(abs_error**2))
    L2_ref = np.sqrt(np.mean(ref_values**2))
    L2_rel = L2_abs / L2_ref

    print(f"[{quantity_name}] N points           : {len(ref_values)}")
    print(f"[{quantity_name}] Absolute L2 (RMS)  : {L2_abs:.6e}")
    print(f"[{quantity_name}] Reference RMS norm : {L2_ref:.6e}")
    print(f"[{quantity_name}] Relative L2 error  : {L2_rel:.4%}\n")

    order = np.arange(len(ref_values))
    if sort_by == "abs_error":
        order = np.argsort(-np.abs(abs_error))
    elif sort_by == "rel_error":
        order = np.argsort(-rel_error)

    n_show = len(order) if n_show is None else min(n_show, len(order))
    print(f"{'x [m]':>10} {'y [m]':>10} {quantity_name+'_own':>14} "
          f"{quantity_name+'_ref':>14} {'abs err':>12} {'rel err %':>10}")
    for i in order[:n_show]:
        x, y = points[i]
        print(f"{x:10.4e} {y:10.4e} {own_values[i]:14.6e} "
              f"{ref_values[i]:14.6e} {abs_error[i]:12.4e} {100*rel_error[i]:10.4f}")

    if plot:
        import matplotlib.pyplot as plt
        fig, axes = plt.subplots(1, 3, figsize=(18, 5))
        for ax, data, title, cmap in zip(
            axes,
            [own_values, ref_values, np.abs(abs_error)],
            [f"{quantity_name} (own)", f"{quantity_name} (reference)", "Absolute error"],
            ["viridis", "viridis", "Reds"],
        ):
            sc = ax.scatter(points[:, 0], points[:, 1], c=data, cmap=cmap)
            ax.set_title(title)
            ax.set_aspect("equal")
            ax.set_xlabel("x [m]")
            ax.set_ylabel("y [m]")
            plt.colorbar(sc, ax=ax)
        plt.tight_layout()
        plt.show()

    return {"own": own_values, "ref": ref_values,
            "abs_error": abs_error, "rel_error": rel_error,
            "L2_abs": L2_abs, "L2_rel": L2_rel}