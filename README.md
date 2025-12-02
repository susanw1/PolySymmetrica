---

# **PolySymmetrica**

### *A parametric polyhedral symmetry engine for OpenSCAD*

PolySymmetrica is a geometry engine that generates and manipulates polyhedra using symmetry, duals, and parametric structure. It provides a clean framework for constructing 3D-printable polyhedral frames, exploring duals and truncations, and attaching arbitrary geometry to faces, edges, or vertices of any polyhedron.

---

## 🔶 Features

* **Primitive symmetry families**

  * tetrahedral
  * octahedral
  * icosahedral

* **Derived solids through duality**

  * cube (dual of octahedron)
  * dodecahedron (dual of icosahedron)

* **Generic placement operators**
  Attach arbitrary geometry to:

  * faces
  * edges
  * vertices

* **Automatic scaling using inter-radius**
  Ensures consistent geometry across derived solids.

* **Fully combinatorial dual operator**
  Works for any convex polyhedron defined in the system.

* **Contextual `$ph_*` variables**
  Child modules automatically receive local coordinate frames and geometry metadata.

---

## 🔶 Quick Start

```scad
use <src/polysymmetrica/polysymmetrica.scad>;

place_on_faces_ir(octahedron(), 30)
    circle(r = $ph_face_poly_radius, $fn = 3);
```

Attach geometry to edges:

```scad
place_on_edges_ir(icosahedron(), 40)
    edge_mount($ph_edge_len);
```

Create dual polyhedra:

```scad
dual_dodeca = poly_dual(icosahedron());
```

---

## 🔶 Directory Structure

```
src/polysymmetrica/
    core/        # math, placement, duals
    models/      # tetra, octa, icosa…
    examples/    # 3D printable frames
    polysymmetrica.scad
docs/
    developer_guide.md
```

---

## 🔶 Documentation

See:
📄 **docs/developer_guide.md**
for a full explanation of:

* polyhedral descriptors
* placement operators
* dual construction
* scaling system
* extending with new solids

---

## 🔶 Roadmap

* Truncation engine
* Archimedean solids
* Catalan solids (via duals)
* Stellation framework
* Compound construction
* More examples and printable models

---

## 🔶 License

MIT — permissive and open for community use.

---


