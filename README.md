Finite element examples for computational physics implemented in Python using **scikit-fem** and the **Gmsh Python API**.

- **scikit-fem** – https://github.com/kinnala/scikit-fem
- **Gmsh** – https://gmsh.info/

Implementations of FEM problems for different physics domains. The examples are not optimized in any way, but they are light enough to run on a laptop.

**Author:** Vertti Aalto  
**Email:** vertti.aalto@gmail.com

---

## Repository Structure

Examples are grouped by physics domain:

```
.
├── elmag/                      # Electromagnetics
│   ├── back-emf/
│   ├── induction-motor/
│   ├── induction-motor-transient/
│   ├── linear-iron-and-pm/
│   ├── nonlinear-iron-and-pm/
│   ├── pmsm-mechanical-dynamics/
│   ├── pmsm-torque-evaluation/
│   └── rotating-pmsm/
│
├── thermal/                    # Heat transfer
│   ├── steady-state-heat-transfer/
│   ├── steady-state-convection/
│   └── transient-convection/
│
├── cfd/                        # fluid dynamics
│   └── poiseuille_flow/
│
├── utils/                      # Helper functions shared by all examples
│
└── README.md
```

### Example layout

Each example follows the same structure:

```
example_name/
├── mesh/
│   ├── geometry.py             # Gmsh geometry file
│   └── mesh.msh or .bdf        # Generated mesh
├── verification/
│   ├── model.mph               # COMSOL model (if available)
│   └── data                    # Reference data (if available)
└── example_notebook.ipynb      # FEM implementation
```

Each example typically includes:

- Geometry generation using the Gmsh Python API
- Generated mesh (`.msh`) and/or (`.bdf`)
- FEM implementation in a Jupyter notebook
- Solution visualization
- Verification against analytical solutions and/or COMSOL (when available)

---

## Branches

- **main** – Stable & verified examples
- **development** – Work in progress
