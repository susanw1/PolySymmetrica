Absolutely — here is the **updated, polished, fully renamed Developer Guide** for **PolySymmetrica**.
---

# **PolySymmetrica – Developer Guide**

*A parametric geometry & symmetry engine for OpenSCAD*

---

## **1. Introduction**

**PolySymmetrica** is a parametric, symmetry-driven polyhedral geometry engine for **OpenSCAD**.

It provides:

* **Primitive symmetry groups**
  Tetrahedral, octahedral, and icosahedral symmetries as foundational structures.

* **Canonical polyhedral descriptors**
  Tetrahedron, octahedron, icosahedron — and via duals: cube and dodecahedron.

* **Generic placement operators**
  Attach arbitrary geometry to *faces*, *edges*, or *vertices* of any polyhedron.

* **A uniform scaling system** via **inter-radius**
  Ensures consistent sizing for primitives, duals, truncations, and compounds.

* A robust, fully combinatorial **dual operator**
  Which correctly generates convex duals of any polyhedron described in the system.

* A clean mechanism for contextual special variables (`$ph_*`)
  allowing child modules to respond to their placement environment.

PolySymmetrica’s design philosophy:

> Polyhedra should be generated, transformed, and explored through **symmetry**, not by hardcoded rotations.

This enables a wide range of applications:

* 3D-printed **polyhedral frames**, edge mounts, and decorative solids
* Dual solids (e.g., cube from octahedron, dodecahedron from icosahedron)
* Future: **Archimedean and Catalan solids**, truncations, rectifications
* Custom geometric elements attached to any symmetry site
* Explorations of symmetry, stellation, and polyhedral combinatorics

PolySymmetrica is as much a *research tool* as a *maker tool*.

---

## **2. Repository Structure**

A suggested layout for clarity:

```
PolySymmetrica/
│
├─ src/
│   ├─ core/
│   │   ├─ funcs.scad           # vector ops, dot, cross, sum, etc.
│   │   ├─ placement.scad       # face/edge/vertex placement operators
│   │   ├─ dual.scad            # poly_dual implementation
│   │   ├─ descriptors.scad     # tetrahedron(), octahedron(), icosahedron(), etc.
│   │   └─ util.scad            # centroids, normals, edges, orientation helpers
│   │
│   ├─ printers/
│   │   └─ frame_mounts.scad    # physical connection parts
│   │
│   └─ examples/
│       ├─ show_octahedron.scad
│       ├─ show_dual_cube.scad
│       ├─ show_dodecahedron.scad
│       └─ …
│
└─ docs/
    ├─ developer_guide.md       # **THIS FILE**
    └─ diagrams/                # optional figures and debug images
```

Each module is deliberately small so that users can take only what they need.

---

## **3. Core Concepts**

### **3.1 Polyhedral Descriptor**

Every polyhedron in PolySymmetrica is defined as:

```
[ verts, faces, unit_edge, e_over_ir ]
```

Where:

| Field       | Meaning                                                   |
| ----------- | --------------------------------------------------------- |
| `verts`     | List of 3D vertex coordinates                             |
| `faces`     | List of faces, each face = ordered list of vertex indices |
| `unit_edge` | Edge length in the canonical embedding                    |
| `e_over_ir` | Ratio of edge length to **inter-radius**                  |

This compact representation allows:

* perfect scaling,
* exact dual construction,
* clean symmetry placement,
* compatibility with truncations and derived families.

---

### **3.2 Inter-Radius Scaling**

PolySymmetrica uses **inter-radius** = radius to the midpoint of an edge as its main input scale.

```
edge_len = inter_radius * e_over_ir
```

This has multiple advantages:

* symmetry-friendly metric (same for duals),
* compatible across truncations,
* makes solids *fit* if they share the same inter-radius,
* intuitive for printing: “how long should an edge be?”

---

### **3.3 Generic Placement Operators**

These operators attach arbitrary geometry to each face/edge/vertex of a polyhedron:

* `place_on_faces(poly, edge_len)`
* `place_on_edges(poly, edge_len)`
* `place_on_vertices(poly, edge_len)`

…or their inter-radius siblings:

* `place_on_faces_ir(poly, inter_radius)`
* `place_on_edges_ir(poly, inter_radius)`
* `place_on_vertices_ir(poly, inter_radius)`

Each operator:

* constructs a local orthonormal frame,
* aligns child geometry accordingly,
* exposes special contextual variables:

| Variable               | Meaning                                               |
| ---------------------- | ----------------------------------------------------- |
| `$ph_facet_idx`        | Index of the face being placed                        |
| `$ph_edge_idx`         | Index of the edge                                     |
| `$ph_vertex_idx`       | Index of the vertex                                   |
| `$ph_edge_len`         | Edge length after scaling                             |
| `$ph_face_poly_radius` | Radius of the face polygon                            |
| `$ph_face_midradius`   | Distance from center to face center                   |
| `$ph_center_local`     | Vector from face center to poly center (local coords) |

These make the system extremely expressive.

---

### **3.4 Dual Operator**

PolySymmetrica implements a **fully correct convex dual construction**:

```
dual_poly = poly_dual(poly);
```

Algorithm:

1. **Dual vertices** = centroids of original faces
2. **Dual faces** = cyclic lists of faces around original vertices
   (computed combinatorially)
3. **Face orientation** corrected via normals
4. **unit_edge** and **e_over_ir** computed automatically

This gives:

* tetra ↔ tetra
* octa ↔ hexa
* icosa ↔ dodeca

…and any future solid will have a correct dual.

---

## **4. Included Polyhedra**

### **4.1 Primitive Symmetry Sets**

```
tetrahedron()
octahedron()
icosahedron()
```

These are the fundamental rotational symmetry families.

### **4.2 Dual Derived Solids**

```
hexahedron() = poly_dual(octahedron());
dodecahedron() = poly_dual(icosahedron());
```

Tetrahedron is self-dual.

---

## **5. Usage Examples**

### **5.1 Octahedron face mounts**

```scad
place_on_faces_ir(octahedron(), 30)
    circle(r = $ph_face_poly_radius, $fn = 3);
```

### **5.2 Edge frames on an icosahedron**

```scad
place_on_edges_ir(icosahedron(), 40)
    edge_mount($ph_edge_len);
```

### **5.3 Dodecahedron from dual of icosa**

```scad
d = poly_dual(icosahedron());

place_on_faces_ir(d, 40)
    circle(r = $ph_face_poly_radius, $fn = 5);
```

### **5.4 Visual debugging: original + dual overlay**

```scad
color("green")
place_on_faces_ir(octahedron(), 30)
    circle(r = $ph_face_poly_radius, $fn = 3);

color("gold")
place_on_faces_ir(hexahedron(), 60)    // note scaling for vertex/face alignment
    square([$ph_edge_len, $ph_edge_len], center = true);
```

---

## **6. Extending PolySymmetrica**

### **6.1 Adding New Polyhedra**

Any convex polyhedron may be added by supplying:

* list of vertices,
* list of faces (cyclic order),
* or even just vertices + faces and letting PolySymmetrica compute:

  * unit_edge,
  * e_over_ir.

The descriptor is then compatible with all placement operators.

### **6.2 Derived Families (future modules)**

Planned:

* **Truncation**: `poly_truncate(poly, t)`
* **Rectification**
* **Cantellation / expansion**
* **Stellation frameworks**
* **Catalan solids** via `poly_dual(poly_truncate(...))`

### **6.3 Custom placement modes**

Just follow the structure of `place_on_faces`:

* compute a frame,
* set `$ph_*` variables,
* descend into `children()`.

---

## **7. Debugging Guide**

If geometry misbehaves:

### **Face centres wrong?**

Ensure centroid uses *all* vertices of the face.

### **Normals pointing inward?**

Use `orient_all_faces_outward`.

### **Dual faces scrambled?**

Check:

* `edges_from_faces`
* `edge_faces_table`
* `faces_around_vertex_rec`

### **Scaling strange?**

Confirm:

* `unit_edge`,
* `e_over_ir`,
* inter_radius argument.

### **Child geometry rotated oddly?**

Use `$ph_*` variables (especially `$ph_center_local`) to confirm the local frame orientation.

---

## **8. Roadmap**

* ✔ Convex dual operator
* ✔ Face/edge/vertex placement operators
* ✔ Primitive Platonic descriptors
* ✔ Dual-generated cube and dodecahedron

**Planned enhancements:**

* ☐ General truncation operator
* ☐ Archimedean generation via truncation
* ☐ Catalan generation via duals
* ☐ Stellation tools
* ☐ Better debugging visualizers
* ☐ PolySymmetrica web documentation / diagrams
* ☐ Interactive examples

---

## **9. Conclusion**

**PolySymmetrica** is a powerful, symmetry-driven geometric engine built for OpenSCAD.
It aims to provide:

* clean mathematical foundations,
* elegant parametric control,
* deep symmetry manipulation,
* and physical print-ready robustness.

It is ideal for:

* makers building 3D polyhedral frameworks,
* educators demonstrating symmetry,
* researchers exploring polyhedral combinatorics,
* hobbyists designing duals, compounds, or stellations.

PolySymmetrica turns the beautiful world of polyhedral geometry into a programmable, extensible toolkit.

---

If you want next steps:

* polish the README,
* generate a logo,
* create example diagrams,
* implement the truncation operator,
* or add Archimedean family support,

just let me know!

