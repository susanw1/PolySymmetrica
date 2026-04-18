# Face Segments

`src/polysymmetrica/core/segments.scad` is the face-local analysis layer for:

- splitting concave or self-intersecting face loops into simple cells
- triangulating those cells consistently
- deriving geometry cut segments where other faces cross the current face plane
- determining which split cells are actually visible from the face-local `+Z` side

This file is about **data and iteration**, not clipping arbitrary 3D geometry.

If you want to clip child geometry to a face or to visible face pieces, use
[`face_regions.md`](face_regions.md) and `src/polysymmetrica/core/face_regions.scad`.

## Mental Model

Use `segments.scad` when you want answers to questions like:

- "How does this star face split under `nonzero` fill?"
- "What cut segments cross this face?"
- "Which cells of this face are still visible?"
- "Iterate the retained pieces and label them / color them / build custom logic."

Use `face_regions.scad` when you want:

- "Take my arbitrary geometry and keep only what lies inside the face region."
- "Clip my geometry to the visible segmented pieces of this face."

In short:

- `segments.scad` = analyze
- `face_regions.scad` = build/clip

## Fill Modes

Most user-facing solid-face work should use `mode="nonzero"`.

- `nonzero`
  Treats self-intersecting loops as filled solids under winding-number rules.
  This is the normal choice for star faces that are intended to behave as solid polyhedron faces.

- `evenodd`
  Treats alternating regions as in/out under parity rules.
  Keep this for debugging, comparison, or when you intentionally want the central overlap of a star to be hollow.

The default in `segments.scad` is now `mode="nonzero"`.

## Main APIs

### `ps_face_segments(face_pts3d_local, mode="nonzero", eps=1e-8)`

Splits a single face loop into simple face cells.

Returns:

```scad
[
    [pts2d, pts3d_local, edge_ids, edge_kinds, cut_entry_ids],
    ...
]
```

Where:

- `pts2d`
  Simple 2D loop in face-local coordinates
- `pts3d_local`
  Matching 3D points in face-local coordinates
- `edge_ids`
  Parent edge indices for the loop edges
- `edge_kinds`
  `"parent"` or `"cut"`
- `cut_entry_ids`
  `undef` for parent edges, otherwise the originating cut-entry id

For ordinary unsplit faces, this usually returns one cell.

### `ps_face_filled_cells(face_pts3d_local, eps=1e-8)`

Semantic wrapper for the normal filled face arrangement under the `nonzero`
fill rule.

This is the usual starting point for true filled-boundary and atom analysis.

### `ps_face_filled_boundary_segments(face_pts3d_local, eps=1e-8)`

Extracts the true filled-boundary subsegments of a possibly self-intersecting
face.

Returns:

```scad
[
    [seg2d, source_edge_idx, source_t0, source_t1, filled_cell_idx],
    ...
]
```

Where:

- `seg2d`
  one true boundary subsegment in face-local 2D
- `source_edge_idx`
  the original face-walk edge this subsegment came from
- `source_t0`, `source_t1`
  parameters along that original source edge
- `filled_cell_idx`
  the filled cell that contributed this boundary edge

Shared edges between two filled cells are internal and are filtered out here.

### `ps_face_filled_boundary_source_edges(face_pts3d_local, eps=1e-8)`

Groups the true filled-boundary subsegments back onto their original face-walk
source edges.

Returns:

```scad
[
    [seg2d, source_edge_idx, spans, inward_sign],
    ...
]
```

Where:

- `seg2d`
  the full original source edge in face-local 2D
- `source_edge_idx`
  the original face-walk edge index
- `spans`
  `[[source_t0, source_t1, filled_cell_idx], ...]` for the true filled spans
  that survive on that source edge
- `inward_sign`
  the agreed filled-side sign for that source edge in its local 2D frame

This is the canonical data API behind
`place_on_face_filled_boundary_source_edges(...)`.

### `place_on_face_segments(mode="nonzero", eps=1e-8)`

Iterator wrapper over `ps_face_segments(...)` for use inside `place_on_faces(...)`.

Provides:

- `$ps_seg_idx`
- `$ps_seg_count`
- `$ps_seg_vertex_count`
- `$ps_seg_pts2d`
- `$ps_seg_pts3d_local`
- `$ps_seg_edge_ids`
- `$ps_seg_edge_kinds`
- `$ps_seg_cut_entry_ids`

### `place_on_face_filled_boundary_segments(eps=1e-8)`

Iterator wrapper over `ps_face_filled_boundary_segments(...)`.

This is the true filled-boundary traversal surface for star/self-crossing faces.
Use it when you want one record per actual surviving boundary subsegment.

### `place_on_face_filled_boundary_source_edges(source_edge_indices=undef, eps=1e-8)`

Iterator wrapper over `ps_face_filled_boundary_source_edges(...)`.

This is the preferred iterator when you want:

- one record per original source edge
- access to the surviving filled boundary spans on that source edge
- a dihedral-centered edge-style frame without artificial subsegment end caps

Provides the usual boundary-edge context plus:

- `$ps_face_boundary_source_spans`
- `$ps_face_boundary_source_span_count`

### `ps_polygon(points, mode="nonzero", eps=1e-8)`

Safe 2D polygon builder for concave or self-intersecting face loops.

This is the user-facing replacement for assuming that raw `polygon(points=...)`
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

This is the provenance consumed by `face_regions.scad`.

### `ps_face_geom_cut_pair_ids(cut_entries, face_idx, face_center_world, face_ex_world, face_ey_world, pair_eps=1e-6)`

Builds a world-stable join id for each cut entry.

This is primarily for:

- cross-face debugging
- verifying that two neighboring faces are talking about the same geometric join
- future exact-region work that needs stable join pairing

Each id is derived from:

- the two face ids involved in the cut
- the world-space intersection line, represented by a quantized anchor point and direction

So corresponding cut edges on two different faces get the same `cut_pair_id`
even though their local edge indices, local cut-entry ids, and local cropped
segment spans can differ.

### `ps_face_geom_cut_segments(...)`

Convenience wrapper returning only the `seg2d` part of
`ps_face_geom_cut_entries(...)`.

Useful for labeling/debug drawing when you do not need the cutter metadata.

### `place_on_face_geom_cut_segments(...)`

Iterator wrapper over geometry-derived cut segments.

Provides:

- `$ps_face_cut_idx`
- `$ps_face_cut_count`
- `$ps_face_cut_segment2d_local`
- `$ps_face_cut_segments2d_local`

### `ps_face_visible_segments(face_pts2d, face_idx, poly_faces_idx, poly_verts_local, eps=1e-8, mode="nonzero", filter_parent=true)`

Splits the current face by geometry cuts and keeps only cells visible from the
face-local `+Z` side.

Returns the same cell structure as `ps_face_segments(...)`:

```scad
[
    [cell_pts2d, cell_pts3d_local, cell_edge_ids, cell_edge_kinds, cell_cut_entry_ids, cell_cut_run_ids],
    ...
]
```

This is the main handoff from analysis to clipping.

Important:

- a visible cell returned here may still be **nonconvex**
- `face_regions.scad` is responsible for decomposing such cells into convex
  atoms when it needs convex local region construction

So `ps_face_visible_segments(...)` should keep returning the true visible cell
topology, not a speculative convexified approximation.

`cell_cut_run_ids` is parallel to `cell_edge_kinds` and `cell_cut_entry_ids`:

- `undef` for parent edges
- a local run id for cut edges

This is how split cut spans are distinguished. If the same cutter contributes
multiple disjoint spans to the same visible cell, those spans keep distinct
`cut_run_id`s even when they share the same `cut_entry_id`.

### `place_on_face_visible_segments(...)`

Iterator wrapper over visible cells.

Provides:

- `$ps_vis_seg_idx`
- `$ps_vis_seg_count`
- `$ps_vis_seg_vertex_count`
- `$ps_vis_seg_pts2d`
- `$ps_vis_seg_pts3d_local`
- `$ps_vis_seg_edge_ids`
- `$ps_vis_seg_edge_kinds`
- `$ps_vis_seg_cut_entry_ids`
- `$ps_vis_seg_cut_pair_ids`
- `$ps_vis_seg_cut_run_ids`

### `face_cut_stencil(face_thk, kerf=0.2, extend=0.5, z_pad=0.2, mode="nonzero", eps=1e-8, filter_parent=true)`

Builds a subtraction body from the geometry-derived cut segments of the current face.

This is still useful as a simple/debug helper and for the compact example in
`main_star_face_render.scad`, but it is **not** the preferred generic path for
clipping arbitrary geometry. For that, use `face_regions.scad`.

## Usage Patterns

### 1. Fill a star face safely

```scad
place_on_faces(poly_antiprism(5, p=2), 30) {
    if ($ps_face_idx == 0)
        linear_extrude(height = 0.4, center = true)
            ps_polygon(points = $ps_face_pts2d, mode = "nonzero");
}
```

### 2. Iterate visible face pieces

```scad
place_on_faces(poly_prism(n=7, p=3), 30) {
    if ($ps_face_idx == 2)
        place_on_face_visible_segments(mode = "nonzero") {
            linear_extrude(height = 0.4, center = true)
                polygon(points = $ps_vis_seg_pts2d);
        }
}
```

### 3. Show raw geometry cut segments

```scad
place_on_faces(poly_antiprism(5, p=2, angle=15), 30) {
    if ($ps_face_idx == 2)
        place_on_face_geom_cut_segments(mode = "nonzero", filter_parent = true) {
            seg = $ps_face_cut_segment2d_local;
            translate([(seg[0][0] + seg[1][0]) / 2, (seg[0][1] + seg[1][1]) / 2, 0.2])
                linear_extrude(height = 0.1)
                    text(str($ps_face_cut_idx), size = 1.5, halign = "center", valign = "center");
        }
}
```

## When To Use What

Use `segments.scad` when:

- you want split cells, cut lines, or visible cells as data
- you want to label or inspect those cells
- you want to write your own custom geometry per cell

Use `face_regions.scad` when:

- you already have arbitrary geometry
- you want to clip that geometry to a face region or visible segmented cells
- you want cut-band subtraction handled as a volume operation

One practical boundary:

- `segments.scad` should not try to pre-convexify or otherwise distort visible
  cells just to make downstream region construction easier
- if a consumer needs convex pieces, that decomposition belongs in
  `face_regions.scad` or the consumer, not here

## Join Metadata

There are now two useful levels of cut metadata:

- `cut_entry_id`
  Local to one face. Good for indexing back into
  `ps_face_geom_cut_entries(...)`.

- `cut_pair_id`
  Stable across the two faces that share the same geometric join.
  Use this when you need to compare or debug the same segmentation edge from
  both sides.

See also:

- [Face region volumes](face_regions.md)
- `src/polysymmetrica/examples/basics/main_star_face_render.scad`
