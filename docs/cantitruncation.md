# Cantitruncation Notes

This document captures current cantitruncation behavior, parameterization, and the non‑grid (trig) solvers. It is a working note, not final user documentation.

## Parameterization

`poly_cantitruncate(poly, t, c, eps, len_eps)`

- **t** controls edge‑point placement along each original edge.
  - `t=0` at original vertices, `t=0.5` at mid‑edge.
  - `t<0` / `t>1` are allowed (anti/hyper).
- **c** controls face/edge expansion.
  - `d_f = -c * ir` (face plane shift along normal)
  - `d_e =  c * ir` (edge‑bisector plane offset)

Topology (current):
- Face cycles: 2n‑gons built from *face‑edge points*.
- Edge cycles: quads built from the same face‑edge points.
- Vertex cycles: 2·valence‑gons built from face‑edge points (hex for valence 3).

## Trig Solver (regular bases)

`solve_cantitruncate_trig(poly, face_idx=0, edge_idx=undef)` returns `[t, c]` directly:

- Let **φ** be the interior angle of the selected face.
- Let **α** be the angle between outward normals of two adjacent faces.
- Let **a** be the original edge length.

Closed‑form:

```
t = 1 / (2 * (1 + sin(φ/2)))
d_f = (1 - 2t) * a / (2 * sin(α/2))
c = |d_f| / ir
```

This avoids any grid search for regular bases (cube, dodecahedron).

## Dominant‑Family Solver (mixed face sizes)

For bases with multiple face families (e.g. cuboctahedron), one global `c` cannot keep all vertex faces planar. Use the dominant‑family path:

- `solve_cantitruncate_dominant(poly, dominant_size)`
  - returns `[t, c_by_size]`
  - `c_by_size` is a list of `[face_size, c]` pairs.

- `solve_cantitruncate_dominant_edges(poly, dominant_size)`
  - returns `[t, c_by_size, c_edge_by_pair]`
  - `c_edge_by_pair` is a list of `[a, b, c]` for edge pairs between face sizes `a` and `b`
  - improves planarity when a single `c_by_size` still warps vertex faces.

- `poly_cantitruncate_families(poly, t, c_by_size, default_c)`
  - applies per‑face‑family `c` values.
  - optional `c_edge_by_pair` overrides edge‑bisector offsets per family pair.

This keeps the dominant family planar and lets secondary families follow.

## Inspecting Planarity

Use `poly_describe(poly, detail=3)` to echo `max_plane_err` per face. This is the max distance of a face’s vertices from its best‑fit plane (0 means perfectly planar).

## Example (cuboctahedron)

```
base = poly_rectify(octahedron()); // cuboctahedron
sol = solve_cantitruncate_dominant_edges(base, 4); // squares dominate
p = poly_cantitruncate_families(base, sol[0], sol[1], c_edge_by_pair=sol[2]);
```
