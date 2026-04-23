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
boundary_span = [seg2d, loop_idx, source_edge_idx, source_t0, source_t1, kind, left_cell_idx, right_cell_idx]
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
  parametric interval on that source edge.
- `kind`
  currently `"source"` for arrangement-derived spans; later layers may add `"cut"` or `"synthetic"`.
- `left_cell_idx`, `right_cell_idx`
  arrangement cell indices adjacent to the span, with one side possibly `undef` on the outside.

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

### `place_on_face_geom_cut_segments(...)`

Iterator wrapper over geometry-derived cut segments.

Provides:

- `$ps_face_cut_idx`
- `$ps_face_cut_count`
- `$ps_face_cut_segment2d_local`
- `$ps_face_cut_segments2d_local`

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
