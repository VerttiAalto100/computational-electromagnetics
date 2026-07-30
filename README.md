# Computational Electromagnetics with FEM

Finite element method (FEM) examples for computational electromagnetics implemented in Python using **scikit-fem** and the **Gmsh Python API**.

- **scikit-fem** – https://github.com/kinnala/scikit-fem
- **Gmsh** – https://gmsh.info/

These problems emphasize clarity and understanding of the FEM workflow, rather than optimization and computational efficiency.

**Author:** Vertti Aalto  
**Email:** vertti.aalto@gmail.com

---

## Repository Structure

```
.
├── example_name/
│   ├── mesh/
│   │   ├── geometry.py         # Gmsh geometry file
│   │   └── mesh.msh or .bdf    # Generated mesh
│   ├── verification/
│   │   ├── model.mph           # COMSOL model (if available)
│   │   └── data                # Reference data (if available)
│   └── example_notebook.ipynb  # FEM implementation
│
├── another_example/
│   └── ...
│
├── utils/
│   # Helper functions
│
└── README.md
```
Each example typically includes:

- Geometry generation using the Gmsh Python API
- Generated mesh (`.msh`) and/or (`.bdf`)
- FEM implementation in a Jupyter notebook
- Solution visualization
- Verification against analytical solutions and/or COMSOL (when available)

---

## Branches

- **master** – Stable & verified examples
- **development** – Work in progress
