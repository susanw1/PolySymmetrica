# Notes

## Poly Validity Checks (Draft Spec)

### Core structural checks (always-on)
- Index bounds: every face index is within `[0, len(verts)-1]`.
- No duplicate indices per face (prevents zero-area degeneracy).
- No zero-length edges: each consecutive edge length > EPS.
- Face planarity: all vertices of a face lie on a plane within EPS.
- Edge manifoldness: each undirected edge appears in exactly 2 faces for closed solids.

### Geometric checks (mode-specific)
- Convex mode: for each face plane `n·x = d`, all vertices satisfy `n·x <= d + eps`.
- No-self-intersection mode: faces must not self-intersect after projection onto their plane.
- Star/stellated mode: allow self-intersecting faces but keep planarity + edge manifoldness.
- Outward orientation: require consistent outward normals in convex mode.

### Topological checks (optional)
- Euler characteristic: `V - E + F == 2` for closed genus-0 solids.
- No duplicate faces (same vertex cycle, possibly rotated).

### Why this matters
Some visually incorrect geometries can still satisfy generic validity checks, so
operation-specific tests (e.g., cantellation face-family counts or adjacency) are
needed to catch issues early.

## Cantellation-Specific Test Ideas
- Face-family counts: for cantellated tetra, expect 4 hex + 6 quad + 4 tri.
- Edge-face adjacency: edge faces should connect the two adjacent vertex-face points
  at each endpoint (prevents crossing across the vertex hex).
- Uniform edge-length checks (optional, for uniform parameter sets).

## Cantellation Normalization (Why We Search)
`poly_cantellate_norm` maps a 0–1 slider to a geometric offset `df` and chooses the
midpoint where edge-faces are squares. There is no simple closed-form formula that
works for arbitrary convex polys because the square condition depends on local
dihedral geometry and edge families. A small numeric search is the most robust
general solution. Closed-form solutions are possible for specific shapes (e.g., cube)
but do not generalize; a search keeps the API uniform across convex polyhedra and
Johnson solids.
