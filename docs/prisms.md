# Prisms And Antiprisms

`core/prisms.scad` provides two generators:

- `poly_prism(n=3, p=1, edge=1, height=undef, height_scale=1)`
- `poly_antiprism(n=3, p=1, edge=1, angle=0, height=undef, height_scale=1)`

Both return a standard PolySymmetrica poly descriptor: `[verts, faces, e_over_ir]`.

## Parameters

### Shared

- `n`: polygon side count (`n >= 3`)
- `p`: polygon step (`1 <= p < n/2`, `gcd(n,p)=1`)
  - `p=1` gives regular n-gon caps
  - `p>1` gives star/compound-style cap ordering `{n,p}` (single-cycle only in v1)
- `edge`: target edge length
- `height`: explicit cap-to-cap separation (if `undef`, regular default is solved)
- `height_scale`: multiplier applied to the chosen base height (`height` if supplied, else regular default)

### Antiprism only

- `angle`: additive top-ring twist offset in degrees relative to the regular antiprism twist
  - actual twist = `180*p/n + angle`
  - `angle=0` gives the exact regular/star antiprism for `{n,p}`

## Regular defaults

### Prism

When `height=undef`, default is `height=edge`, so base edges and vertical edges match.

### Antiprism

When `height=undef`, height is solved from:

- base ring radius: `R = edge / (2*sin(180/n))`
- base ring radius (regular/star): `R = edge / (2*sin(180*p/n))`
- twist: `theta = 180/n + angle`
- twist (regular/star): `theta = 180*p/n + angle`
- lateral edge condition: `edge^2 = 2*R^2*(1-cos(theta)) + h^2`

So:

`h = sqrt(edge^2 - 2*R^2*(1-cos(theta)))`

If that radicand is negative, the requested `(n, edge, angle)` is infeasible and the function asserts.

## Usage examples

```scad
use <src/polysymmetrica/core/prisms.scad>
use <src/polysymmetrica/examples/truncation/util_demo.scad>

demo(poly_prism(6));                         // regular hexagonal prism
demo(poly_prism(5, height=0.8));             // explicit height
demo(poly_prism(5, p=2));                    // pentagrammic prism {5/2}
demo(poly_antiprism(5));                     // regular pentagonal antiprism
demo(poly_antiprism(5, angle=8));            // extra twist
demo(poly_antiprism(6, height=1.0));         // explicit antiprism height
demo(poly_antiprism(5, p=2));                // pentagrammic antiprism {5/2}
```

Example scene:

- `src/polysymmetrica/examples/basics/main_prisms.scad`

Negative-test runners (expected to fail with assertions):

- `src/tests/negative/prism_bad_p0.scad`
- `src/tests/negative/prism_bad_half.scad`
- `src/tests/negative/prism_bad_coprime.scad`
- or run all negatives: `src/tests/run_negative_all.sh`

## Known rendering limitation (star faces)

For star caps (`p > 1`), faces are self-intersecting loops in 3D.  
`poly_render(...)` currently delegates to `polyhedron(points, faces)`, and OpenSCAD F6 backends may triangulate/fill those loops differently than 2D `polygon()` even-odd behavior. This can appear as unexpected filled regions on star faces.

Current workaround for visual/debug output:

- render face-local plates from placement (`$ps_face_pts2d`) via `polygon()`/`linear_extrude`, instead of a single `polyhedron(...)` cap for those faces.

Repro/inspection scene:

- `src/polysymmetrica/examples/basics/main_star_face_render.scad`
