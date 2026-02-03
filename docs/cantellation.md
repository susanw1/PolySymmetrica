# Cantellation Notes

This document captures current cantellation behavior, parameterization, and solver helpers. It is a working note, not final user documentation.

## Parameterization

`poly_cantellate(poly, df=undef, c=undef, df_max=undef, steps=16, family_edge_idx=0, eps, len_eps)`

- **df** controls face offsets (how far original faces move along their normals).
- **c** provides a normalized knob; `c=0.5` targets square edge faces (via `cantellate_square_df`).
- If `df` is omitted, `c` (or a default `c=0.5`) is used to derive a `df`.
- `df_max` bounds the normalized mapping; `steps` controls the square‑target search.

Topology (current):
- Face cycles: original faces (expanded) become larger polygons.
- Edge cycles: quads derived from each original edge.
- Vertex cycles: polygons derived from each original vertex (valence‑gons).

## Components & Current Limitations

Cantellation is constructed from three components:

1) **Face cycles**
   - Derived by shifting original face planes by `df` and intersecting with edge‑bisector planes.
   - Intended to remain planar by construction.

2) **Edge cycles** (quads)
   - Built from the two edge points per original edge (one per adjacent face).
   - These are the “rectangular/square” faces of cantellation.

3) **Vertex cycles** (valence‑gons)
   - Built from edge‑adjacent points around each original vertex.

Current limitations:
- Mixed‑family bases can show slight warping depending on how offsets are balanced.
- Extreme offsets (large `df`) can yield self‑intersections or degenerate faces; use with care.

## Solvers / Helpers

- `cantellate_square_df(poly, df_min, df_max, steps, family_edge_idx, eps)`
  - Searches for a `df` that makes a chosen edge‑family as square as possible.

- `poly_cantellate_norm(poly, c, df_max=undef, steps=16, family_edge_idx=0, eps, len_eps)`
  - Normalized cantellation (maps `c in [0,1]` to a `df` range).
  - Intended as the user‑friendly entry point for many examples.

## Inspecting Planarity

Use `poly_describe(poly, detail=3)` to echo `max_plane_err` per face. This is the max distance of a face’s vertices from its best‑fit plane.

## Usage Examples

### Basic cantellation (explicit df)
```
p = poly_cantellate(hexahedron(), df = 0.2);
```

### Normalized cantellation
```
p = poly_cantellate_norm(hexahedron(), 0.5);
```

### Square‑targeted face family
```
base = poly_rectify(octahedron()); // cuboctahedron
// pick a representative edge from the family you want squared
sol_df = cantellate_square_df(base, 0.0, 1.0, 40, family_edge_idx=0);
p = poly_cantellate(base, sol_df);
```

### Inspect planarity
```
poly_describe(p, detail=3);
```
