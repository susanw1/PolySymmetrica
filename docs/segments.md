# Face Segments

`src/polysymmetrica/core/segments.scad` is the face-local analysis layer for:

- splitting concave or self-intersecting face loops into simple arrangement cells
- triangulating those arrangement cells consistently
- deriving geometry cut segments where other faces cross the current face plane
- determining which split arrangement cells are actually visible from the face-local `+Z` side

This file is about data and iteration, not arbitrary 3D clipping.

See also [face_arrangement.md](face_arrangement.md) for the next planned layer:
explicit arrangement, boundary-model, face-cell, and atom APIs for self-crossing
faces.

## Mental Model

Use `segments.scad` when you want answers to questions like:

- "What is the raw arrangement induced by this self-crossing face?"
- "What is the true filled boundary of this self-crossing face?"
- "How does this self-crossing face split under `nonzero` fill?"
- "What cut segments cross this face?"
- "Which face cells of this face are still visible?"
- "Iterate the retained pieces and label them / color them / build custom logic."

In short:

- `segments.scad` = analyze

## Definitions

- `face loop`
  The original ordered walk of the face in face-local coordinates. This may be concave or self-crossing.

- `simple cell`
  A simple polygonal region traced from that face loop after splitting at crossings.
  Here "simple" means non-self-intersecting with one closed boundary loop.
  It does not imply convex, triangular, minimal, or canonical.
  In this document, a simple cell is a 2D face-local arrangement region.

- `parent edge`
  An edge of a derived cell that lies on one of the original face-loop edges.

- `cut edge`
  An edge of a derived cell created by splitting the face with crossings or with geometry-derived cut segments.

- `filled cell`
  A simple 2D face cell that is kept by the chosen fill rule (`nonzero`, `evenodd`, or `all`).

- `atom`
  A convex 2D sub-piece of a filled face cell, introduced when a filled cell is too complex or non-convex for robust clipping or cutter generation.
  An atom is derived from a filled cell.
  It is not the same thing as a 3D body/volume.

- `face boundary`
  The boundary of the original face loop as drawn, before any fill-rule interpretation.

- `filled boundary`
  The boundary of the region that remains after applying the fill rule. For a self-crossing face, this can differ from the raw face loop.

- `geometry cut`
  A segment formed where another face crosses the current face plane, expressed in the current face-local 2D coordinates.

- `visible cell`
  A filled 2D face cell, or sub-cell after further splitting by geometry cuts, that is still visible from the face-local `+Z` side.

## Fill Modes

Most user-facing solid-face work should use `mode="nonzero"`.

- `nonzero`
  Treats self-intersecting loops as filled solids under winding-number rules.
  This is the normal choice for star faces intended to behave as solid polyhedron faces.

- `evenodd`
  Treats alternating regions as in/out under parity rules.
  Keep this for debugging, comparison, or when you intentionally want the central overlap of a star to be hollow.

- `all`
  Keeps every extracted cycle for inspection/debugging.

## Main APIs

### `ps_face_segments(face_pts3d_local, mode="nonzero", eps=1e-8)`

Splits one face loop into simple 2D face cells.

Returns:

```scad
[
    [pts2d, pts3d_local, edge_ids, edge_kinds],
    ...
]
```

Where:

- `pts2d`
  simple 2D loop in face-local coordinates
- `pts3d_local`
  matching 3D points in face-local coordinates
- `edge_ids`
  source face-loop edge indices for the loop edges
- `edge_kinds`
  current edge lineage tags for the loop edges

Currently, arrangement-backed `ps_face_segments(...)` emits `"source"` for inherited
spans. Later segment layers may also emit tags such as `"cut"` or `"synthetic"`.

For ordinary unsplit faces, this usually returns one face cell.

### `ps_face_arrangement(face_pts3d_local, eps=1e-8)`

Builds the raw planar arrangement induced by the face loop.

Returns:

```scad
[
    face_pts2d,
    crossings,
    nodes,
    spans,
    cells
]
```

This is the unfiltered arrangement-level product.
It does not apply a fill rule and does not imply convex decomposition.

Sub-records are:

```scad
crossing = [source_edge_a, source_t_a, source_edge_b, source_t_b, pt2d, node_idx]
node     = [pt2d, kind]
span     = [seg2d, node_a, node_b, source_edge_idx, source_t0, source_t1, kind]
cell     = [pts2d, pts3d_local, node_ids, span_ids, signed_area]
```

Where:

- `crossings`
  one record per proper self-crossing of the original face loop.
- `source_edge_a`, `source_edge_b`
  indices of the original face-loop edges that cross.
- `source_t_a`, `source_t_b`
  parametric positions of the crossing on those source edges.
- `pt2d`
  crossing point in face-local 2D.
- `node_idx`
  index into `nodes` for that crossing point.

- `nodes`
  unique arrangement nodes, including original face vertices and inserted crossings.
- `kind`
  currently `source_vertex` or `crossing`.

- `spans`
  split edge spans between arrangement nodes.
- `seg2d`
  span endpoints in face-local 2D.
- `node_a`, `node_b`
  indices into `nodes` for the span endpoints.
- `source_edge_idx`
  original face-loop edge that this span came from.
- `source_t0`, `source_t1`
  parametric interval on that original source edge.
- span `kind`
  currently `"source"` for arrangement spans derived directly from the face loop.

- `cells`
  all traced simple arrangement cells before any fill-rule filtering.
- `pts2d`, `pts3d_local`
  ordered loop points for that cell in 2D and matching face-local 3D.
- `node_ids`
  indices into `nodes` around the cell boundary.
- `span_ids`
  indices into `spans` around the cell boundary.
- `signed_area`
  signed 2D loop area, useful for orientation and fill-rule logic.

### `ps_face_boundary_model(face_pts3d_local, mode="nonzero", eps=1e-8)`

Derives the true filled boundary from the arrangement for the chosen fill rule.

Returns:

```scad
[
    mode,
    filled_cell_ids,
    boundary_loops,
    boundary_spans
]
```

Where:

- `filled_cell_ids`
  indices into the arrangement `cells` array for the cells kept by the fill rule
- `boundary_loops`
  the reconstructed filled-boundary loops
- `boundary_spans`
  oriented boundary spans carrying lineage back to source edges where applicable

This is the face-like boundary product to use when you need the resulting
perimeter rather than just the arrangement cells.

Sub-records are:

```scad
boundary_loop = [pts2d, span_ids]
boundary_span = [seg2d, loop_idx, source_edge_idx, source_t0, source_t1, kind, filled_cell_idx, other_cell_idx]
```

Where:

- `boundary_loops`
  one record per traced filled-boundary loop.
- loop `pts2d`
  ordered 2D vertices around that filled-boundary loop.
- loop `span_ids`
  indices into `boundary_spans` in loop order.

- `boundary_spans`
  one record per oriented span on the filled boundary.
- `seg2d`
  oriented 2D segment in loop order.
- `loop_idx`
  which `boundary_loop` this span belongs to.
- `source_edge_idx`
  original face-loop edge that this span descends from, when applicable.
- `source_t0`, `source_t1`
  parametric interval on that source edge, oriented to match `seg2d`.
- `kind`
  currently `"source"` for arrangement-derived spans; later layers may add `"cut"` or `"synthetic"`.
- `filled_cell_idx`, `other_cell_idx`
  arrangement cell indices adjacent to the span, with `other_cell_idx` possibly
  `undef` on the outside.

### `ps_face_filled_boundary_source_edges(face_pts3d_local, mode="nonzero", eps=1e-8)`

Groups true filled-boundary spans by the original source edge they came from.

Returns:

```scad
[
    [source_edge_idx, source_seg2d, source_boundary_spans],
    ...
]
```

Where:

- `source_edge_idx`
  original face-loop edge index.
- `source_seg2d`
  the original source edge segment in face-local 2D.
- `source_boundary_spans`
  all surviving filled-boundary spans descended from that source edge.

Only source edges with at least one surviving filled-boundary span are returned.
The returned source-edge records are ordered by `source_edge_idx`.

Sub-records are:

```scad
source_boundary_span = [
    boundary_span_idx,
    seg2d,
    loop_idx,
    source_t0,
    source_t1,
    kind,
    filled_cell_idx,
    other_cell_idx,
    filled_side
]
```

Where:

- `boundary_span_idx`
  index into `ps_face_boundary_model(...)[3]`.
- `seg2d`
  oriented span segment in filled-boundary traversal order.
- `loop_idx`
  filled-boundary loop containing this span.
- `source_t0`, `source_t1`
  source-edge parameter range, oriented to match `seg2d`.
- `kind`
  currently `"source"` for arrangement-derived spans.
- `filled_cell_idx`, `other_cell_idx`
  arrangement cells on either side of the span.
- `filled_side`
  `+1` when the filled region lies on the left of `seg2d`, `-1` when it lies
  on the right, and `0` only for degenerate/ambiguous spans.

Use this when source-edge lineage is the primary structure. Use
`place_on_face_boundary_spans(...)` instead when each surviving span needs its
own local frame or adjacent-face/dihedral context.

### `place_on_face_segments(mode="nonzero", eps=1e-8)`

Iterator wrapper over `ps_face_segments(...)` for use inside `place_on_faces(...)`.

Provides:

- `$ps_seg_idx`
- `$ps_seg_count`
- `$ps_seg_vertex_count`
- `$ps_seg_pts2d`
- `$ps_seg_pts3d_local`
- `$ps_seg_parent_face_edge_idx`
- `$ps_seg_edge_kind`

### `place_on_face_filled_boundary_source_edges(mode="nonzero", eps=1e-8, coords="element")`

Iterator wrapper over `ps_face_filled_boundary_source_edges(...)` for use inside
`place_on_faces(...)`.

`coords` controls the child coordinate system:

- `"element"` (default): place children in a normalized source-edge coordinate
  system.
- `"parent"`: leave children in the enclosing face-local coordinate system.

In element coordinates:

- local `+X` runs along the source edge, reversed when needed
- local `+Y` is the in-face left normal of local `+X`
- local `+Z` is face-local `+Z`
- the filled region is on local `-Y` when the source edge has an unambiguous
  filled side

Provides:

- `$ps_boundary_source_edge_idx`
- `$ps_boundary_source_edge_count`
- `$ps_boundary_source_edge_len`
- `$ps_boundary_source_edge_segment2d_local`
- `$ps_boundary_source_edge_span_count`
- `$ps_boundary_source_edge_spans`
- `$ps_boundary_source_edge_boundary_span_idxs`
- `$ps_boundary_source_edge_span_segments2d_local`
- `$ps_boundary_source_edge_span_t_ranges`
- `$ps_boundary_source_edge_sides`
- `$ps_boundary_source_edge_frame_reversed`
- `$ps_boundary_source_edge_frame_source_side`
- `$ps_boundary_source_edge_frame_side`
- `$ps_boundary_source_edge_span_t_ranges_local`
- `$ps_boundary_source_edge_span_sides_local`

The raw record metadata remains source-oriented. The `*_frame_*` variables
describe the child placement frame, and `*_span_t_ranges_local` is the
same span interval expressed in that element coordinate system. If one source edge contributes
spans with mixed filled sides, use
`$ps_boundary_source_edge_span_sides_local` per span rather than assuming every
span has the representative frame side.

When `coords="parent"`, the exposed segment records such as
`$ps_boundary_source_edge_segment2d_local` and
`$ps_boundary_source_edge_span_segments2d_local` can be drawn directly in the
current face-local coordinates.

### `place_on_face_boundary_spans(mode="nonzero", eps=1e-8, coords="element")`

Iterator wrapper over internal dihedral-aware boundary-span site records for
use inside `place_on_faces(...)`.

Use this when later geometry needs both:

- the face-local filled boundary span itself
- the adjacent-face / dihedral context inherited from the original source edge
- the per-span filled side in the face plane

This is the key distinction from source-edge-level metadata:
the same original source edge can contribute multiple surviving boundary spans,
and those spans can legitimately differ in direction/order.

`coords` controls the child coordinate system:

- `"element"` (default): place children in the current boundary-span coordinate
  system.
- `"parent"`: leave children in the enclosing face-local coordinate system.

Provides:

- `$ps_boundary_span_idx`
- `$ps_boundary_span_count`
- `$ps_boundary_span_len`
- `$ps_boundary_span_segment2d_local`
- `$ps_boundary_span_loop_idx`
- `$ps_boundary_span_source_edge_idx`
- `$ps_boundary_span_source_t0`
- `$ps_boundary_span_source_t1`
- `$ps_boundary_span_kind`
- `$ps_boundary_span_filled_cell_idx`
- `$ps_boundary_span_other_cell_idx`
- `$ps_boundary_span_adj_face_idx`
- `$ps_boundary_span_dihedral`
- `$ps_boundary_span_adj_face_normal_local`
- `$ps_boundary_span_filled_side`
- `$ps_boundary_span_adj_face_dir_span_local`

Where:

- local frame:
  - `+X` runs along the oriented boundary span
  - `+Y` is the in-face left normal of that span
  - `+Z` is face-local `+Z`
- `$ps_boundary_span_source_t0/$ps_boundary_span_source_t1`
  are oriented to match `$ps_boundary_span_segment2d_local`
- `$ps_boundary_span_filled_side`
  is `+1` when the filled region lies on the left side of the oriented span,
  `-1` when it lies on the right, and `0` only for degenerate/ambiguous cases
- `$ps_boundary_span_adj_face_dir_span_local`
  is the direction, in the current span-local frame, of the adjacent face plane
  perpendicular to the shared span and oriented toward the current face's local
  `+Z` side. This is often more useful for visualization and anti-interference
  work than the adjacent-face normal alone, especially when a single source edge
  contributes multiple surviving boundary spans. The underlying source-edge
  direction is projected into the current face plane before this direction is
  derived, so warped/non-planar faces still preserve the span-frame `Y/Z`
  contract.

### `ps_polygon(points, mode="nonzero", eps=1e-8)`

Safe 2D polygon builder for concave or self-intersecting face loops.

This is the user-facing replacement for assuming raw `polygon(points=...)`
will behave well on star/self-intersecting loops.

### `ps_face_geom_cut_entries(face_pts2d, face_idx, poly_faces_idx, poly_verts_local, eps=1e-8, mode="nonzero", filter_parent=true)`

Derives the geometry cuts made by other faces crossing the current face plane.

Returns:

```scad
[
    [seg2d, cutter_face_idx, cut_dihed],
    ...
]
```

Where:

- `seg2d`
  the cut segment in face-local 2D
- `cutter_face_idx`
  the source face that generated this cut
- `cut_dihed`
  the corresponding cut dihedral

### `ps_face_geom_cut_segments(...)`

Convenience wrapper returning only the `seg2d` part of
`ps_face_geom_cut_entries(...)`.

### `ps_face_foreign_intrusion_records(face_pts2d, face_idx, poly_faces_idx, poly_verts_local, eps=1e-8, mode="nonzero", filter_parent=true)`

Derives exact face-local foreign intrusion records for other faces that cross
the current face plane.

Returns:

```scad
[
    [record_kind, target_face_idx, foreign_kind, foreign_idx, seg2d_local, cut_dihed, confidence],
    ...
]
```

Where:

- `record_kind`
  currently `"face_plane_cut"`.
- `target_face_idx`
  the face being tested for intrusion.
- `foreign_kind`
  currently `"face"`.
- `foreign_idx`
  the source face that generated this intrusion.
- `seg2d_local`
  the intrusion segment in target face-local 2D.
- `cut_dihed`
  the face-plane cut angle inherited from `ps_face_geom_cut_entries(...)`.
- `confidence`
  currently `"exact"`; later clearance/envelope candidates should use a
  different classification such as `"conservative"`.

Accessor helpers:

- `ps_intrusion_kind(record)`
- `ps_intrusion_target_face_idx(record)`
- `ps_intrusion_foreign_kind(record)`
- `ps_intrusion_foreign_idx(record)`
- `ps_intrusion_segment2d_local(record)`
- `ps_intrusion_dihedral(record)`
- `ps_intrusion_confidence(record)`

This API intentionally does not answer clearance questions such as “does a
4mm-wide edge site come within 2mm of this face?”. Those are conservative
placement-envelope candidates and belong in a later geometry/proxy layer.

### `place_on_face_geom_cut_segments(...)`

Iterator wrapper over geometry-derived cut segments.

Provides:

- `$ps_face_cut_idx`
- `$ps_face_cut_count`
- `$ps_face_cut_segment2d_local`
- `$ps_face_cut_segments2d_local`

### `place_on_face_foreign_intrusions(...)`

Iterator wrapper over `ps_face_foreign_intrusion_records(...)`.

Provides:

- `$ps_intrusion_idx`
- `$ps_intrusion_count`
- `$ps_intrusion_record`
- `$ps_intrusion_kind`
- `$ps_intrusion_target_face_idx`
- `$ps_intrusion_foreign_kind`
- `$ps_intrusion_foreign_idx`
- `$ps_intrusion_segment2d_local`
- `$ps_intrusion_dihedral`
- `$ps_intrusion_confidence`

### `ps_face_visible_segments(face_pts2d, face_idx, poly_faces_idx, poly_verts_local, eps=1e-8, mode="nonzero", filter_parent=true)`

Splits the current face by geometry cuts and keeps only 2D face cells visible from the
face-local `+Z` side.

Returns:

```scad
[
    [cell_pts2d, cell_pts3d_local, cell_edge_ids, cell_edge_kinds],
    ...
]
```

Where:

- `cell_pts2d`
  visible cell boundary in face-local 2D.
- `cell_pts3d_local`
  matching face-local 3D points.
- `cell_edge_ids`
  source-edge ids for inherited edges where applicable.
- `cell_edge_kinds`
  currently `"parent"` or `"cut"` depending on whether that visible-cell edge comes from the pre-cut cell boundary or from a geometry cut.

### `place_on_face_visible_segments(...)`

Iterator wrapper over `ps_face_visible_segments(...)`.

Provides:

- `$ps_vis_seg_idx`
- `$ps_vis_seg_count`
- `$ps_vis_seg_vertex_count`
- `$ps_vis_seg_pts2d`
- `$ps_vis_seg_pts3d_local`
- `$ps_vis_seg_edge_ids`
- `$ps_vis_seg_edge_kinds`

## Notes

- `segments.scad` is intentionally face-local.
- Most printing/punch-through workflows should use `mode="nonzero"`.
- These APIs are intended to stay analyzable as data. Arbitrary 3D clipping belongs in later geometry layers.
