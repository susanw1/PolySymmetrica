# Construction Notes

This document records the current polyhedral construction helpers in
`src/polysymmetrica/core/construction.scad`.

These are explicit topology tools intended to support later Johnson-style
construction, slicing, opening, and capping workflows.

## Current API

```scad
poly_delete_faces(poly, fids, cap=false, cleanup=true, cleanup_eps=1e-8)
poly_boundary_loops(poly)
poly_cap_loops(poly, loops=undef, cleanup=true, cleanup_eps=1e-8)
poly_slice(poly, plane_pt, plane_n, keep="above", cap=true, cleanup=true, cleanup_eps=1e-8)
poly_pyramid(n=4, p=1, edge=1, height=undef, height_scale=1)
poly_cupola(n=3, edge=1, height=undef, height_scale=1)
poly_attach(p1, p2, f1=0, f2=0, rotate_step=0, scale_mode="fit_edge", mirror=false, eps=1e-8, cleanup=true, cleanup_eps=1e-8)
```

## Design Intent

This layer is deliberately explicit.

Rather than starting with high-level "remove vertex and repair" style helpers,
the construction API currently exposes the lower-level operations that those
wrappers would need anyway:

1. remove faces,
2. recover open boundary loops,
3. cap selected loops,
4. slice by a plane and optionally cap the result.

The first direct Johnson-oriented constructor now also lives here:

5. build a regular/star pyramid from a `{n,p}` base.
6. build an exact n-gonal cupola from concentric top/base polygons.

That keeps the topology changes visible and makes later Johnson/cupola/pyramid
operations easier to reason about.

## Function Notes

### `poly_delete_faces(...)`

Deletes one face or a list of faces.

- `fids` may be a scalar or list.
- if `cap=false`, the result is generally open and should not be expected to
  satisfy `poly_valid(..., "closed")`
- if `cap=true`, the deleted openings are recovered via `poly_boundary_loops(...)`
  and capped automatically

Example:

```scad
q = poly_delete_faces(hexahedron(), 0, cap=false);
r = poly_delete_faces(hexahedron(), [0,1], cap=true);
```

### `poly_boundary_loops(...)`

Returns the current open boundary loops of a poly as vertex-index cycles.

- loops do not repeat the first vertex at the end
- this is intended for open meshes
- on a closed-valid poly it should return `[]`

Boundary detection currently uses a simple and reliable rule:

- build directed face edges from the current face list
- an edge is a boundary edge if its undirected form appears in exactly one face
- keep the directed occurrence from that surviving face

This is not heavily optimized; correctness is preferred here.

### `poly_cap_loops(...)`

Adds caps over boundary loops.

- if `loops=undef`, all current boundary loops are capped
- if `loops` is supplied, it should be a list of vertex-index cycles

The result is passed through `poly_cleanup(..., fix_winding=true, ...)`, so the
caller does not need to determine cap winding manually.

Example:

```scad
open_poly = poly_delete_faces(hexahedron(), 0, cap=false, cleanup=false);
closed_poly = poly_cap_loops(open_poly);
```

### `poly_slice(...)`

Slices a closed polyhedron by a plane and keeps one side.

Parameters:

- `plane_pt`: a point on the plane
- `plane_n`: plane normal
- `keep`: `"above"` or `"below"` relative to `plane_n`
- `cap`: whether to cap the cut opening

Implementation outline:

1. clip each face polygon against the chosen half-space
2. rebuild an open poly from the clipped face polygons
3. optionally cap the resulting boundary loops

Example:

```scad
top_half_open  = poly_slice(hexahedron(), [0,0,0], [0,0,1], keep="above", cap=false);
top_half_closed = poly_slice(hexahedron(), [0,0,0], [0,0,1], keep="above", cap=true);
```

### `poly_pyramid(...)`

Builds a regular or star pyramid over a `{n,p}` base.

Parameters:

- `n`: number of base vertices
- `p`: polygon step (`p=1` for an ordinary pyramid)
- `edge`: base edge length
- `height`: explicit apex-to-base-plane height
- `height_scale`: multiplier on the chosen height

If `height=undef`, the constructor solves the apex height so that the side edges
match `edge`. That gives the exact Johnson J1/J2 forms for `n=4` and `n=5`.

Examples:

```scad
j1 = poly_pyramid(4);
j2 = poly_pyramid(5);
star = poly_pyramid(5, 2);
```

Runnable demo:

- `src/polysymmetrica/examples/basics/main_johnsons.scad`

### `poly_cupola(...)`

Builds an exact n-gonal cupola with:

- a top n-gon,
- a base 2n-gon,
- alternating square and triangular side faces.

Parameters:

- `n`: top polygon arity
- `edge`: common target edge length
- `height`: explicit cupola height
- `height_scale`: multiplier on the chosen height

If `height=undef`, the constructor solves the exact cupola height from the
apothem difference between the top n-gon and the base 2n-gon. This gives the
exact Johnson cupolae:

- `poly_cupola(3)` => J3 triangular cupola
- `poly_cupola(4)` => J4 square cupola
- `poly_cupola(5)` => J5 pentagonal cupola

Examples:

```scad
j3 = poly_cupola(3);
j4 = poly_cupola(4);
j5 = poly_cupola(5);
```

### `poly_attach(...)`

`poly_attach(...)` is also part of the construction layer and is documented in
[attach.md](attach.md). It is expected to work together with slicing, capping,
and later cap/belt constructors for Johnson-style assemblies.

## Current Assumptions

- `poly_slice(...)` expects a closed-valid source poly:
  `poly_valid(poly, "closed")`
- cleanup is structural only; this is not a full self-intersection or manifold
  surgery framework
- capping currently creates one cap face per recovered loop

## Current Tests

Coverage lives in:

- `src/tests/core/TestConstruction.scad`

Current tests cover:

- deleting one or more cube faces
- recovering resulting boundary loops
- recapping those loops
- slicing a cube through its midplane, both open and capped

Each `Test*.scad` file is directly runnable as well as via `run_all.scad`.

## Likely Next Steps

This construction layer is intended to support later helpers such as:

- `poly_delete_vertices(..., repair=true)`
- `poly_delete_edges(..., repair=true)`
- higher-level Johnson constructors
- slice-driven cap extraction for pyramids/cupolae/rotundae

The likely next practical use is building pyramids and cupola-style solids on
top of `poly_slice(...)`, `poly_attach(...)`, and prism/antiprism belts.
