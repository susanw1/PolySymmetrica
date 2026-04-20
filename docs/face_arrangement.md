# Face Arrangement and Boundary Models

This note defines the intended next layer for self-crossing face work.

It sits one level below the current `segments.scad` user-facing APIs and is
meant to separate four ideas that are currently too entangled:

- arrangement extraction
- filled-boundary reconstruction
- filled simple cells
- convex atoms used for clipping/CSG

The goal is to make self-crossing faces behave like analyzable data first,
before any later printing or punch-through logic consumes them.

## Definitions

- `arrangement`
  The planar graph obtained by splitting a face loop at all of its crossings.

- `simple cell`
  A non-self-intersecting polygonal region traced from that arrangement.
  A simple cell is not necessarily convex.

- `filled cell`
  A simple cell kept by the chosen fill rule, such as `nonzero` or `evenodd`.

- `boundary model`
  The true filled boundary derived from the arrangement and the chosen fill
  rule. This is the face-like product we would use when we want the resulting
  shape to behave like a replacement face with explicit lineage.

- `boundary span`
  One boundary segment of the boundary model. A boundary span may come from an
  original source edge, a split of a source edge, or a synthetic/cut edge.

- `atom`
  A convex sub-piece of a filled cell, introduced when a filled cell is too
  complex or non-convex for robust clipping or cutter generation.
  An atom is always convex and is derived from a filled cell.
  It is not the same thing as a simple cell.

- `visible cell`
  A filled cell, or later an atom, that remains visible after face-local
  occlusion/cut filtering.

## Intended Layering

The intended progression is:

1. `ps_face_arrangement(...)`
2. `ps_face_boundary_model(...)`
3. `ps_face_cells(...)`
4. `ps_face_atoms(...)`
5. `ps_face_visible_cells(...)`

That means:

- arrangement first
- boundary reconstruction second
- filled simple cells third
- convex decomposition only when needed
- visibility as a separate later filter

## Relationship To Current `segments.scad`

The current `segments.scad` layer already provides useful face-local products,
but it is still organized mainly around fill modes and visible segments.

Roughly:

- `ps_face_segments(...)`
  currently behaves closest to `ps_face_cells(...)`
- `mode="all"`
  means “return all simple cells from the arrangement”, not “return convex atoms”
- `ps_face_visible_segments(...)`
  is a later visibility-stage result, not the canonical base layer

So this new arrangement/boundary-model layer is not meant to replace the
current segment APIs immediately. It is meant to give them a clearer foundation.

## Proposed Root Data Products

### `ps_face_arrangement(face_pts3d_local, eps=1e-8)`

Build the planar arrangement induced by the input face loop.

Proposed root record:

```scad
[
    face_pts2d,
    crossings,
    nodes,
    spans,
    cells
]
```

Where:

- `face_pts2d`
  The original face loop in face-local 2D.

- `crossings`
  One record per proper self-crossing.

- `nodes`
  Unique arrangement nodes, including original vertices and inserted crossings.

- `spans`
  Split edge spans between arrangement nodes.

- `cells`
  All traced simple cells in the arrangement, before any fill-rule filtering.

Suggested sub-records:

```scad
crossing = [source_edge_a, source_t_a, source_edge_b, source_t_b, pt2d, node_idx]
node     = [pt2d, kind]
span     = [seg2d, node_a, node_b, source_edge_idx, source_t0, source_t1, kind]
cell     = [pts2d, node_ids, span_ids, signed_area]
```

Where `kind` is intentionally explicit and small, for example:

- node kind: `source_vertex`, `crossing`
- span kind: `source`, `synthetic`, `cut`

## Proposed Boundary Product

### `ps_face_boundary_model(face_pts3d_local, mode="nonzero", eps=1e-8)`

Derive the true filled boundary from the arrangement.

This is the product to use when we want the result to behave like a proper
replacement face perimeter rather than just a collection of cells.

Proposed root record:

```scad
[
    mode,
    filled_cell_ids,
    boundary_loops,
    boundary_spans
]
```

Suggested sub-records:

```scad
boundary_loop = [pts2d, span_ids]
boundary_span = [seg2d, loop_idx, source_edge_idx, source_t0, source_t1, kind, left_cell_idx, right_cell_idx]
```

Boundary spans should preserve lineage where possible:

- `source_edge_idx`
  original face-loop edge index, or `undef`
- `source_t0`, `source_t1`
  parametric interval on the source edge when applicable
- `kind`
  for example `source`, `cut`, `synthetic`
- `left_cell_idx`, `right_cell_idx`
  or equivalent ownership/adjacency information

## Proposed Filled-Cell Product

### `ps_face_cells(face_pts3d_local, mode="nonzero", eps=1e-8)`

Return the filled simple cells selected by the fill rule.

This is the topological cell layer.
Cells are simple, but not necessarily convex.

Proposed record:

```scad
[
    [pts2d, span_ids, source_edge_ids],
    ...
]
```

The exact record shape can be refined later, but it should be clearly derived
from the arrangement/boundary model rather than recomputed ad hoc.

## Proposed Atom Product

### `ps_face_atoms(face_pts3d_local, mode="nonzero", eps=1e-8)`

Return a convex decomposition of the filled cells.

This is the geometry-work-unit layer.
Atoms should be the preferred substrate for:

- clipping
- cutter generation
- robust CSG/interference logic

Atoms should carry lineage back to their parent filled cell:

```scad
atom = [pts2d, parent_cell_idx, edge_ids, edge_kinds]
```

The exact decomposition method can vary, but the API meaning should stay:

- atoms are convex
- atoms are derived from filled cells
- atoms are not the same thing as arrangement cells

## Likely Iterator Surfaces

Once the data products exist, the thin iterator wrappers become much clearer:

- `place_on_face_boundary_loops(...)`
- `place_on_face_boundary_spans(...)`
- `place_on_face_cells(...)`
- `place_on_face_atoms(...)`
- possibly `place_on_face_cell_edges(...)`
- possibly `place_on_face_atom_edges(...)`

## Non-Goals

- Do not treat current `mode="all"` as if it already means convex
  decomposition.
- Do not fold visible-cell logic into the first arrangement/boundary-model
  pass.
- Do not drag `proxy_interaction` or printable carve logic into this layer.

## Practical Direction

The first implementation step should be:

1. define `ps_face_arrangement(...)`
2. define `ps_face_boundary_model(...)`
3. only then add `ps_face_cells(...)` and `ps_face_atoms(...)`

That keeps terminology and data shapes stable before the next wave of
geometry/code changes.
