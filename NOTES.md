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
