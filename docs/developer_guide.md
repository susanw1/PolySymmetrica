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
│   │   │   ├─ truncation.scad      # truncate/rectify/chamfer/cantellate/cantitruncate
│   │   │   ├─ params.scad          # shared params_overrides helpers
│   │   │   ├─ solvers.scad         # parameter solvers (cantitruncate, etc.)
│   │   │   ├─ transform.scad       # site-based mesh assembly
│   │   │   ├─ transform_util.scad  # shared transform helpers
│   │   │   ├─ validate.scad        # poly validity checks
│   │   │   ├─ cleanup.scad         # structural cleanup/normalization
│   │   │   ├─ construction.scad    # delete/cap/slice construction helpers
│   │   │   ├─ construction.scad    # slice/cap/attach construction helpers
│   │   │   ├─ prisms.scad          # prism/antiprism generators
│   │   │   ├─ render.scad          # poly_describe, render helpers
│   │   │   └─ duals.scad           # poly_dual and helpers
│   │   │
│   │   ├─ models/
│   │   │   ├─ platonics_all.scad
│   │   │   ├─ archimedians_all.scad
│   │   │   ├─ catalans_all.scad
│   │   │   ├─ johnsons_all.scad
│   │   │   ├─ tetrahedron.scad
│   │   │   ├─ octahedron.scad
│   │   │   └─ icosahedron.scad
│   │   │
│   │   └─ examples/
│   │       ├─ basics/                 # core placement + archimedean/catalan demos
│   │       ├─ truncation/             # trunc/rectify/cantellate/cantitruncate demos
│   │       ├─ poly-frame/             # frame + mount experiments
│   │       └─ printing/               # printable face/frame experiments
│   └─ tests/...
│
└─ docs/
    ├─ developer_guide.md
    ├─ params_overrides.md
    ├─ cantitruncation.md
    ├─ cantellation.md
    ├─ snubs.md
    ├─ attach.md
    └─ images/
```

Each module is deliberately small so that users can take only what they need.

Related deep-dive notes:

- [Cantellation notes](cantellation.md)
- [Cantitruncation notes](cantitruncation.md)
- [Snub notes](snubs.md)
- [Prisms and antiprisms](prisms.md)
- [Params overrides](params_overrides.md)
- [Face attachment](attach.md)
- [Construction helpers](construction.md)

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

**Face winding:** PolySymmetrica follows OpenSCAD’s LHR convention.
Faces are ordered **clockwise when viewed from outside**.

---

### **3.2 Classification (Families)**

`poly_classify(poly, detail=1, eps=1e-6, radius=1, include_geom=false)` groups faces, edges, and vertices into families.

- Returns `[face_fams, edge_fams, vert_fams]`, each as `[key, idxs]`.
- `detail`:
  - `0`: topology only (face size, edge-adjacent face sizes, vertex valence + cyclic face sizes).
  - `1`: topology + one-pass neighbour refinement.
  - `2`: iterated neighbour refinement (propagates by `radius`).
- `include_geom`: when true, appends average edge length to keys (rounded by `eps`).
- `radius`: controls propagation depth for `detail=2` refinement.

Helper for printing: `show_poly(poly, detail, eps, radius, include_geom)`.

---

### **3.3 Inter-Radius Scaling**

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

### **3.4 Generic Placement Operators**

These operators attach arbitrary geometry to each face/edge/vertex of a polyhedron:

* `place_on_faces(poly, inter_radius, edge_len=undef, classify=undef, classify_opts=undef)`
* `place_on_edges(poly, inter_radius, edge_len=undef, classify=undef, classify_opts=undef)`
* `place_on_vertices(poly, inter_radius, edge_len=undef, classify=undef, classify_opts=undef)`

Or calculate it based on edge length:

* `place_on_faces(poly, edge_len = 2.5, classify=cls)`
* `place_on_edges(poly, edge_len = 2.5, classify=cls)`
* `place_on_vertices(poly, edge_len = 2.5, classify=cls)`

Classification controls:

* `classify`: optional precomputed value from `poly_classify(...)`; preferred for consistency and speed.
* `classify_opts`: optional `[detail, eps, radius, include_geom]`; used only when `classify` is not passed.
* if both are omitted, placement remains geometry-only (no classification is performed and family vars are `undef`).


Each operator:

* constructs a local orthonormal frame,
* aligns child geometry accordingly,
* exposes special contextual variables:

PolySymmetrica exposes per-placement metadata via `$ps_*` variables.

### **3.5 Stable public placement context metadata**

The variables documented in this section are the **stable public placement
context metadata** for the core placement operators:

- `place_on_faces(...)`
- `place_on_edges(...)`
- `place_on_vertices(...)`

If a `$ps_*` variable is **not** documented here, it should currently be
treated as internal or experimental rather than part of the public contract.

Likewise, nested iterator modules in other files, such as `segments.scad` or
`face_regions.scad`, may expose additional `$ps_*` variables of their own; they
are governed by their own docs and are not automatically part of this core
placement context metadata surface.

The contract is about **semantics**, not implementation details. Future helper
functions such as site-record builders should preserve these meanings.

#### Common guarantees

All three core placement operators guarantee that children are evaluated in a
local orthonormal frame centered on the current placement site.

The exact axis-construction details are part of each placement operator's
geometry and may evolve. Consumers should rely on the documented `$ps_*`
metadata contract rather than inferring hidden frame semantics from the current
implementation.

Every placement site also exposes the following common metadata:

| Variable | Faces | Verts | Edges | Meaning |
| --- | :---: | :---: | :---: | --- |
| `$ps_edge_len` | ✅ | ✅ | ✅ | Scale-driving edge length for this placement context. On faces this is the scale-derived face edge length; on edges it is the actual length of the current edge; on vertices it is the target edge length used to scale the poly. |
| `$ps_poly_center_local` | ✅ | ✅ | ✅ | Vector from the placement-site origin to the polyhedral center, expressed in the current local frame. |
| `$ps_face_family_count` | ✅ | ✅ | ✅ | Total number of face families in the active classification, or `undef` when no classification context is active. |
| `$ps_edge_family_count` | ✅ | ✅ | ✅ | Total number of edge families in the active classification, or `undef` when no classification context is active. |
| `$ps_vertex_family_count` | ✅ | ✅ | ✅ | Total number of vertex families in the active classification, or `undef` when no classification context is active. |

#### Face placement contract

When inside `place_on_faces(...)`, the following face-specific metadata is
public:

| Variable | Meaning |
| --- | --- |
| `$ps_face_idx` | Index of the current face in the poly descriptor. |
| `$ps_vertex_count` | Number of vertices in the current face. |
| `$ps_face_midradius` | Distance from poly center to face center, in world units after scaling. |
| `$ps_face_radius` | Mean distance from the face center to its vertices, in world units after scaling. |
| `$ps_face_pts2d` | Face vertices in face-local 2D coordinates, suitable for `polygon(points=...)`. |
| `$ps_face_pts3d_local` | Face vertices in face-local 3D coordinates, with local `z` mean-centered. |
| `$ps_face_planarity_err` | Maximum local-z deviation across the current face vertices. |
| `$ps_face_is_planar` | Boolean planarity flag derived from `$ps_face_planarity_err`. |
| `$ps_face_family_id` | Family id of the current face from `poly_classify(...)`, or `undef` if no classification is active. |
| `$ps_face_neighbors_idx` | Adjacent face indices per face edge, aligned with `$ps_face_pts2d` order. Boundary edges are `undef`. |
| `$ps_face_dihedrals` | Dihedral angles per face edge, aligned with `$ps_face_pts2d` order. Boundary edges use a fallback self-normal value and should not be treated as meaningful two-face dihedrals. |

#### Vertex placement contract

When inside `place_on_vertices(...)`, the following vertex-specific metadata is
public:

| Variable | Meaning |
| --- | --- |
| `$ps_vertex_idx` | Index of the current vertex in the poly descriptor. |
| `$ps_vertex_valence` | Number of incident edges at the current vertex. |
| `$ps_vertex_neighbors_idx` | Neighbor vertex indices incident to the current vertex. Their order is available but not currently guaranteed as a public ordering contract. |
| `$ps_vertex_neighbor_pts_local` | Vectors from the current vertex to each neighbor, expressed in vertex-local coordinates. |
| `$ps_vert_radius` | Distance from the poly center to the current vertex, in world units after scaling. |
| `$ps_vertex_family_id` | Family id of the current vertex from `poly_classify(...)`, or `undef` if no classification is active. |

#### Edge placement contract

When inside `place_on_edges(...)`, the following edge-specific metadata is
public:

The edge frame itself is currently dihedral-centered when a usable adjacent
face pair exists: local `+X` runs along the edge and local `+Z` follows the
projected bisector of the two adjacent outward face normals. Boundary or
degenerate edges fall back to the older radial frame.

| Variable | Meaning |
| --- | --- |
| `$ps_edge_idx` | Index of the current edge in the poly descriptor. |
| `$ps_edge_midradius` | Distance from the poly center to the current edge midpoint, in world units after scaling. |
| `$ps_edge_pts_local` | Edge endpoints in edge-local coordinates, typically `[[ -L/2, 0, 0 ], [ +L/2, 0, 0 ]]`. |
| `$ps_edge_verts_idx` | Vertex indices of the current edge `[v0, v1]`. |
| `$ps_edge_adj_faces_idx` | Face indices adjacent to the current edge. Closed manifold edges usually have two. |
| `$ps_edge_family_id` | Family id of the current edge from `poly_classify(...)`, or `undef` if no classification is active. |

#### Classification note

Family ids are only meaningful relative to one classification context. If you
need stable matching between placement and transform/override logic, pass the
same precomputed `classify` object to all consumers.

These make the system extremely expressive.

### **3.6 Face Segment Notes**

See [segments.md](segments.md) for the face-local analysis layer:

- `ps_face_arrangement(...)`
- `ps_face_segments(...)`
- `ps_face_geom_cut_entries(...)`
- `ps_face_geom_cut_segments(...)`
- `ps_face_foreign_intrusion_records(...)`
- `ps_face_visible_segments(...)`

Use this layer when you want analyzable face-local cells, cuts, and visibility,
not arbitrary 3D clipping.

See [face_regions.md](face_regions.md) for positive face-local volumes built
from those boundary spans:

- `ps_face_anti_interference_shells(...)`
- `ps_face_anti_interference_volume(...)`

See [face_arrangement.md](face_arrangement.md) for the planned next layer:

- `ps_face_arrangement(...)`
- `ps_face_boundary_model(...)`
- `ps_face_cells(...)`
- `ps_face_atoms(...)`

### **3.7 Cantitruncation Notes**

See [cantitruncation.md](cantitruncation.md) for current parameterization, trig solver, and dominant‑family notes.

### **3.8 Cantellation Notes**

See [cantellation.md](cantellation.md) for current parameterization, helpers, and planarity notes.

### **3.9 Snub Notes**

See [snubs.md](snubs.md) for snub usage, default solving, helper notes, and current caveats.

### **3.10 Shared Params Overrides**

PolySymmetrica supports operator-agnostic structured overrides using shared row formats:

```scad
["face"|"vert"|"edge", "all", ["key", value], ...]
["face"|"vert"|"edge", "family", family_id, ["key", value], ...]
["face"|"vert"|"edge", "id", id_or_ids, ["key", value], ...]
```

Implemented in `core/params.scad` and used by operators such as `poly_snub(...)`.
Rows can set multiple keys, and precedence is `id > family > all`.
For hot paths, compile rows once with:

- `ps_params_compile_key(...)`
- `ps_params_compile_specs(...)`

and read in O(1) with:

- `ps_compiled_param_get(...)`

For schema details, compile-spec format, and examples, see:

- [params_overrides.md](params_overrides.md)
- [`src/polysymmetrica/examples/basics/main_params.scad`](../src/polysymmetrica/examples/basics/main_params.scad)

When using family-targeted overrides (`"family"` rows), keep classification context consistent (`detail`, `radius`, `include_geom`) across all stages; best practice is to classify once and reuse that object.

### **3.11 Face Attachment**

`poly_attach(...)` composes two closed polys by welding selected faces:

- selected seam faces are aligned and opposed in normal direction,
- seam faces are removed from each input,
- seam vertices are merged using cleanup tolerance,
- output is asserted closed-valid.

API:

```scad
poly_attach(
    p1, p2,
    f1=0, f2=0,        // f1 accepts index or [idx...]
    rotate_step=0,
    scale_mode="fit_edge",
    mirror=false,
    eps=1e-8,
    cleanup=true,
    cleanup_eps=1e-8
)
```

See:

- [attach.md](attach.md)
- [`src/polysymmetrica/examples/basics/main_attach.scad`](../src/polysymmetrica/examples/basics/main_attach.scad)

### **3.12 Validation Modes**

`poly_valid(poly, mode, eps)` supports progressively stricter checks:

| Mode       | Structural checks | Manifold closedness | Self-intersections | Outward + convex |
| ---------- | ----------------- | ------------------- | ------------------ | ---------------- |
| `"struct"` | ✅                | ☐                   | ☐                  | ☐                |
| `"closed"` | ✅                | ✅                  | disallowed         | ☐                |
| `"star_ok"`| ✅                | ✅                  | allowed            | ☐                |
| `"convex"` | ✅                | ✅                  | disallowed         | ✅               |

Use `"struct"` for early pipeline/debug stages, `"closed"` for printable manifold outputs, `"star_ok"` for stellated/self-intersecting faces, and `"convex"` when strict outward convex geometry is required.

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
* `*_pts3d_local`
  3D points in the local frame, useful when a face may be non-planar.

  * Example: `$ps_face_pts3d_local`

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

### **3.13 Dual Operator**

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

### **4.3 Aggregated Sets**

```scad
platonics_all();      // [["tetrahedron", fn], ...]
archimedians_all();   // 13 Archimedeans as [name, fn]
catalans_all();       // Catalan duals as [name, fn]
johnsons_all();       // current Johnson previews as [name, fn]
```

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

* **Truncation**
* **Rectification**
* **Chamfer**
* **Cantellation / expansion**
* **Cantitruncation**
* **Catalans** via `poly_dual(...)` of Archimedeans

Planned:

* **Stellation frameworks**

### **6.3 Transform parameter sign conventions**

Operator parameters follow their construction, so the sign may differ by operator.

* **Chamfer**: positive `t` offsets face planes **inward** (toward the poly center); negative `t` produces an anti‑chamfer.
* **Cantellate/expand**: positive offsets typically push new edge/vertex features **outward** from the original poly.
* **Cantitruncate**: `t` is a truncation-style edge fraction; `c` offsets face/edge planes by `±c * ir` (see [cantitruncation.md](cantitruncation.md)).
* **Snub**: `angle` controls in-face twist, with `df/de/c` controlling offsets; defaults are solver-derived (see [snubs.md](snubs.md)).

### **6.4 Custom placement modes**

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
* ✔ Truncation/rectification/chamfer/cantellate/cantitruncation/snub (`core/truncation.scad`)
* ✔ Solver helpers (`core/solvers.scad`)
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
