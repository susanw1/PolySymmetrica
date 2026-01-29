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

* A clean mechanism for contextual special variables (`$ps_*`)
  allowing child modules to respond to their placement environment.

PolySymmetrica’s design philosophy:

> Polyhedra should be generated, transformed, and explored through **symmetry**, not by hardcoded rotations.

This enables a wide range of applications:

* 3D-printed **polyhedral frames**, edge mounts, and decorative solids
* Dual solids (e.g., cube from octahedron, dodecahedron from icosahedron)
* **Archimedean and Catalan solids** via truncation/duals are now supported; stellations are next
* Custom geometric elements attached to any symmetry site
* Explorations of symmetry, stellation, and polyhedral combinatorics


---

## **2. Repository Structure**

```
PolySymmetrica/
│
├─ README.md
├─ src/
│   ├─ polysymmetrica/
│   │   ├─ core/
│   │   │   ├─ funcs.scad           # math, vector, centroid, helpers
│   │   │   ├─ placement.scad       # face/edge/vertex placement
│   │   │   ├─ truncation.scad      # poly_truncate and helpers
│   │   │   └─ duals.scad           # poly_dual and helpers
│   │   │
│   │   ├─ models/
│   │   │   ├─ tetrahedron.scad
│   │   │   ├─ octahedron.scad
│   │   │   └─ icosahedron.scad
│   │   │
│   │   └─ examples/
│   │       ├─ basics/                 # core placement demos
│   │       ├─ truncation/             # trunc/rectify/cantellate demos
│   │       ├─ poly-frame/             # frame + mount experiments
│   │       └─ printing/               # printable face/frame experiments
│   └─ tests/...
│
└─ docs/
    ├─ developer_guide.md
    └─ diagrams/
```

Each module is deliberately small so that users can take only what they need.

---

## **3. Core Concepts**

### **3.1 Polyhedral Descriptor**

Every polyhedron in PolySymmetrica is defined as:

```
[ verts, faces, e_over_ir ]
```

Where:

| Field       | Meaning                                                                  |
| ----------- | ---------------------------------------------------------                |
| `verts`     | List of 3D vertex coordinates (unit-edge by convention; not required)     |
| `faces`     | List of faces, each face = ordered list of vertex indices                |
| `e_over_ir` | Ratio of edge length to **inter-radius**                                 |

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

* `place_on_faces(poly, inter_radius)`
* `place_on_edges(poly, inter_radius)`
* `place_on_vertices(poly, inter_radius)`

Or calculate it based on edge length:

* `place_on_faces(poly, edge_len = 2.5)`
* `place_on_edges(poly, edge_len = 2.5)`
* `place_on_vertices(poly, edge_len = 2.5)`


Each operator:

* constructs a local orthonormal frame,
* aligns child geometry accordingly,
* exposes special contextual variables:

PolySymmetrica exposes per-placement metadata via `$ps_*` variables.

(Conventions: ✅ = available, ☐ = not set)

| Variable                        | Faces | Verts  | Edges | Meaning                                                                                                                                                        |
| ------------------------------- | :---: | :---:  | :---: | -------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `$ps_face_idx`                 |   ✅  |   ☐   |   ☐   | Index of the face being placed (0..N-1)                                                                                                                        |
| `$ps_vertex_count`              |   ✅  |   ☐   |   ☐   | Number of vertices of this face (length of `$ps_face_pts2d`)                                                                                                   |
| `$ps_face_midradius`            |   ✅  |   ☐   |   ☐   | Distance from poly centre to face centre (world units; scale-derived)                                                                                          |
| `$ps_face_radius`              |   ✅  |   ☐   |   ☐   | Mean distance from face centre to its vertices (useful even for irregular faces; scale-derived)                                                                |
| `$ps_face_pts2d`                |   ✅  |   ☐   |   ☐   | Face polygon vertices in **face-local 2D coords** `[[x,y]...]` suitable for `polygon(points=...)`                                                              |
| `$ps_face_neighbors_idx`       |   ✅  |   ☐   |   ☐   | Adjacent face indices per face edge (aligned with `$ps_face_pts2d` order; `undef` for boundary edges)                                                         |
| `$ps_face_dihedrals`            |   ✅  |   ☐   |   ☐   | Dihedral angles per face edge, degrees (aligned with `$ps_face_pts2d` order; internal dihedral for closed manifolds)                                           |
| `$ps_vertex_valence`            |   ☐   |   ✅   |   ☐   | Number of incident edges at this vertex                                                                                                                       |
| `$ps_vertex_neighbors_idx`      |   ☐   |   ✅   |   ☐   | Neighbor vertex indices incident to this vertex (order currently unspecified)                                                                                 |
| `$ps_vertex_neighbor_pts_local` |   ☐   |   ✅   |   ☐   | Vectors from vertex to each neighbor in **vertex-local coords** `[[x,y,z]...]`                                                                                |
| `$ps_vert_radius`               |   ☐   |   ✅   |   ☐   | Distance from poly centre to this vertex (world units; scale-derived)                                                                                         |
| `$ps_vertex_idx`                |   ☐   |   ✅   |   ☐   | Index of the vertex being placed (0..K-1)                                                                                                                     |
| `$ps_edge_idx`                  |   ☐   |   ☐   |   ✅  | Index of the edge being placed (0..M-1)                                                                                                                        |
| `$ps_edge_midradius`            |   ☐   |   ☐   |   ✅  | Distance from poly centre to edge midpoint (world units; scale-derived)                                                                                        |
| `$ps_edge_pts_local`            |   ☐   |   ☐   |   ✅  | Edge endpoints in **edge-local coords**, typically `[[ -L/2,0,0 ], [ +L/2,0,0 ]]`                                                                              |
| `$ps_edge_verts_idx`            |   ☐   |   ☐   |   ✅  | Vertex indices of this edge `[v0, v1]`                                                                                                                         |
| `$ps_edge_adj_faces_idx`        |   ☐   |   ☐   |   ✅  | Face indices adjacent to this edge (usually 2 for closed manifold polys)                                                                                       |
| `$ps_edge_len`                  |   ✅  |   ✅   |   ✅  | **Faces:** mean edge length for the face (scale-derived); **Edges:** *actual* length of this edge; **Verts:** target edge length parameter passed to placement |
| `$ps_poly_center_local`         |   ✅  |   ✅   |   ✅  | Polyhedral centre vector expressed in the current local coord frame (Faces/Edges/Verts). For vertices: `[0,0,-$ps_vert_radius]` by construction                |

These make the system extremely expressive.

#### Naming conventions for `$ps_*` variables

The following conventions are used consistently to make their meaning predictable:

##### Index vs geometry

* `*_idx`
  An **index** into the polyhedral descriptor (e.g. face index, edge index, vertex index).

  * Examples: `$ps_face_idx`, `$ps_edge_idx`, `$ps_vertex_idx`
* `*_idx` (plural values)
  A list of indices.

  * Examples: `$ps_edge_verts_idx`, `$ps_edge_adj_faces_idx`, `$ps_vertex_neighbors_idx`

##### Coordinate systems

* `*_local`
  Values expressed in the **local coordinate frame** of the current placement (face / edge / vertex).

  * Examples: `$ps_poly_center_local`, `$ps_edge_pts_local`, `$ps_vertex_neighbor_pts_local`
* `*_pts2d`
  2D points in the local XY plane, suitable for `polygon(points=...)` and **ordered LHR (CW from outside)**.

  * Example: `$ps_face_pts2d`

##### Distances and radii

* `*_radius`
  A distance measured from the polyhedral centre to a geometric feature.

  * Examples: `$ps_face_midradius`, `$ps_edge_midradius`, `$ps_vert_radius`
* `*_len`
  A linear length (edge length or scale-derived length).
  Note: `$ps_edge_len` may represent either an actual edge length (edge placement) or a scale-derived / target value (face or vertex placement), depending on context.

##### Counts

* `*_count`
  Cardinality of a local feature.

  * Example: `$ps_vertex_count` (number of vertices in a face)

---

##### Design intent

These conventions aim to:

* keep `$ps_*` variables **self-describing**
* avoid ambiguity between indices and coordinates
* make it obvious which coordinate system a value belongs to
* allow child modules to remain geometry-agnostic and reusable



---

### **3.4 Dual Operator**

PolySymmetrica implements a **fully correct convex dual construction**:

```
dual_poly = poly_dual(poly);
```

Algorithm:

1. **Dual vertices** = polar face normals (plane normals divided by offset)
2. **Dual faces** = cyclic lists of faces around original vertices
   (computed combinatorially)
3. **Face orientation** corrected via normals (LHR; clockwise when viewed from outside)
4. **e_over_ir** computed automatically
5. Scaling helper function to match dual size to original for overlay

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
place_on_faces(octahedron(), 30)
    circle(r = $ps_face_radius, $fn = $ps_vertex_count);
```
<img src="images/5.1 octa.png" alt="Simple octahedron" width="200px"/>

### **5.2 Edge frames on an icosahedron**

```scad
place_on_edges(icosahedron(), 40)
    cube([$ps_edge_len, 5, 1], center = true);
```
<img src="images/5.2 icosa.png" alt="Edge frames on an icosahedron" width="200px"/>

### **5.3 Exploded dodecahedron from dual of icosa**

```scad
d = poly_dual(icosahedron());

place_on_faces(d, 40)
    translate([0,0, -$ps_face_midradius + 20])
        cylinder(r1 = 0, r2 = $ps_face_radius, h = $ps_face_midradius, $fn = $ps_vertex_count);
```
<img src="images/5.3 dodeca-exploded.png" alt="Exploded dodecahedron" width="200px"/>

### **5.4 Truncated icosa**

```scad
!place_on_faces(poly_truncate(icosahedron()), 40)
    color($ps_vertex_count == 5? "blue" : "orange")
        cylinder(r = $ps_face_radius, $fn = $ps_vertex_count, h = 0.2);
```
<img src="images/5.4 trunc-icosa.png" alt="Truncated icosahedron" width="200px"/>

### **5.5 Catalan from dual of truncated icosahedron**

```scad
d = poly_dual(poly_truncate(icosahedron()));

place_on_faces(d, 40)
    linear_extrude(height = 0.2) polygon(points = $ps_face_pts2d);
```
<img src="images/5.5 dual-trunc-icosa.png" alt="Dual of truncated icosahedron" width="200px"/>

### **5.6 Visual debugging: original + dual overlay**

```scad
color("red", alpha = 0.5)
place_on_faces(octahedron(), 30)
    cylinder(r = $ps_face_radius, $fn = $ps_vertex_count, h = 0.2);

color("gold", alpha = 0.3)
place_on_faces(poly_dual(octahedron()), 30)    // note scaling for vertex/face alignment
    cylinder(r = $ps_face_radius, $fn = $ps_vertex_count, h = 0.2);
```
<img src="images/5.6 cube-octa.png" alt="Overlayed cube on octahedron" width="200px"/>

---

## **6. Extending PolySymmetrica**

### **6.1 Adding New Polyhedra**

Any convex polyhedron may be added by supplying:

* list of vertices,
* list of faces (cyclic order),
* or even just vertices + faces and letting PolySymmetrica compute *e_over_ir*.

The descriptor is then compatible with all placement operators.

### **6.2 Derived Families**

Implemented:

* **Rectification**
* **Cantellation / expansion**

Planned:

* **Stellation frameworks**

### **6.4 Transform parameter sign conventions**

Operator parameters follow their construction, so the sign may differ by operator.

* **Chamfer**: positive `t` offsets face planes **inward** (toward the poly center); negative `t` produces an anti‑chamfer.
* **Cantellate/expand**: positive offsets typically push new edge/vertex features **outward** from the original poly.

### **6.3 Custom placement modes**

Just follow the structure of `place_on_faces`:

* compute a frame,
* set `$ps_*` variables,
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

* `e_over_ir`,
* inter_radius argument,
* use scale_dual(poly, dual) to calculate dual's inter_radius.

---

## **8. Roadmap**

* ✔ Face/edge/vertex placement operators
* ✔ Primitive Platonic descriptors
* ✔ Dual-generated cube and dodecahedron
* ✔ Convex dual operator
* ✔ General truncation/rectification/cantellation operators (`core/truncation.scad`)
* ✔ Archimedean generation via truncation/cantellation
* ✔ Catalan generation via duals
* ✔ Catalan scaling to ensure dual overlay

**Planned enhancements:**

* ☐ Stellation tools
* ☐ Better debugging visualizers
* ☐ PolySymmetrica web documentation / diagrams

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
