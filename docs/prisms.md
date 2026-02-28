# Prisms And Antiprisms

`core/prisms.scad` provides two generators:

- `poly_prism(n=3, edge=1, height=undef, height_scale=1)`
- `poly_antiprism(n=3, edge=1, angle=0, height=undef, height_scale=1)`

Both return a standard PolySymmetrica poly descriptor: `[verts, faces, e_over_ir]`.

## Parameters

### Shared

- `n`: polygon side count (`n >= 3`)
- `edge`: target edge length
- `height`: explicit cap-to-cap separation (if `undef`, regular default is solved)
- `height_scale`: multiplier applied to the chosen base height (`height` if supplied, else regular default)

### Antiprism only

- `angle`: additive top-ring twist offset in degrees relative to the regular antiprism twist
  - actual twist = `180/n + angle`
  - `angle=0` gives the exact regular antiprism

## Regular defaults

### Prism

When `height=undef`, default is `height=edge`, so base edges and vertical edges match.

### Antiprism

When `height=undef`, height is solved from:

- base ring radius: `R = edge / (2*sin(180/n))`
- twist: `theta = 180/n + angle`
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
demo(poly_antiprism(5));                     // regular pentagonal antiprism
demo(poly_antiprism(5, angle=8));            // extra twist
demo(poly_antiprism(6, height=1.0));         // explicit antiprism height
```

Example scene:

- `src/polysymmetrica/examples/basics/main_prisms.scad`

