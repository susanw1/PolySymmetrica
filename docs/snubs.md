# Snub Notes

This document captures current `poly_snub` behavior, solver strategy, parameters,
and practical usage/debug notes. It is a working note (implementation-focused),
not final polished user docs.

## API

```scad
poly_snub(
    poly,
    angle = undef,
    c = undef,
    df = undef,
    de = undef,
    handedness = 1,
    eps = 1e-8,
    len_eps = 1e-6,
    params_overrides = undef
)
```

Related helper:

```scad
ps_snub_default_params_overrides(poly, handedness=1, eps=1e-9)
```

## Parameter Semantics

- `angle`: in-face twist angle (degrees). `undef` means "solve default angle".
- `df`: face-plane offset distance (caller convention: positive is outward).
- `de`: edge/vertex inset distance (same outward convention as `df`).
- `c`: normalized control for `de` via cantellation mapping.
  - Internally mapped by `_ps_cantellate_df_from_c_linear(...)`.
  - If `c` is supplied and `de` is not, then `de` comes from `c`.
- `handedness`: diagonal split/chirality selector (`>=0` vs `<0`).
- `params_overrides`: structured per-element overrides (see `docs/params_overrides.md`).

Important fallback rules:

- If all of `angle/c/df/de` are `undef`, full auto defaults are solved.
- If `angle` is `undef` but other controls are provided:
  - with explicit `de`: angle is solved against `(df,de)`,
  - else with explicit `c`: angle solve uses `c` mapping,
  - else with explicit `df`: angle solve uses `(df,df)`.

## Current Build Model

Snub is constructed in the standard two-phase style:

1. **Structural cycles** from per-face transformed points:
   - original face cycles,
   - edge cycles (quads),
   - vertex cycles.
2. **Final face set**:
   - keep face cycles,
   - split each edge quad into two triangles according to `handedness`,
   - keep vertex cycles.

Face orientation is normalized with `ps_orient_all_faces_outward(...)`.

## Solver/Objective Notes

Snub default solving currently uses:

- `_ps_snub_eval_errors_base(...)` returning:
  - `ev[0]`: global edge-length spread,
  - `ev[1]`: edge-triangle equilateral error,
  - `ev[2]`: vertex-face planarity error,
  - `ev[3]`: edge-triangle isosceles error.

- Regular path (`_ps_is_regular_base(poly)`):
  - uses `_ps_snub_default_params_full(...)`,
  - searches a small `(c, angle)` region with local refinement,
  - currently keeps `df` coupled to `de` (`r=1`) in this path.

- Non-regular path:
  - uses heuristic search with stronger planarity penalty and capped `c` range.

### Handedness Consistency Fix

Snub now canonicalizes face-pair order per edge before scoring and building edge cycles:

- `_ps_snub_oriented_edge_faces(edge, fpair, face_n, verts0)`

This prevents dependence on incidental `edge_faces` ordering and keeps handedness application consistent.

## Structured Overrides (Recommended for Advanced Use)

`poly_snub` compiles and applies overrides for:

- face keys: `df`, `angle`
- vertex keys: `c`, `de`

Example:

```scad
params = [
    ["face", "all", ["angle", 14], ["df", 0.08]],
    ["vert", "family", 0, ["c", 0.05]]
];

q = poly_snub(p, params_overrides=params);
```

See full schema in `docs/params_overrides.md`.

## Usage Examples

### 1) Plain auto snub (regular base)

```scad
q = poly_snub(hexahedron());
```

### 2) Explicit controls

```scad
q = poly_snub(hexahedron(), c=0.07, angle=16.5);
```

### 3) Solve only angle for fixed offsets

```scad
q = poly_snub(hexahedron(), df=0.10, de=0.10, angle=undef);
```

### 4) Use structured default rows directly

```scad
rows = ps_snub_default_params_overrides(hexahedron());
q = poly_snub(hexahedron(), params_overrides=rows);
```

### 5) Family-specific tweak (requires family IDs from classify)

```scad
cls = poly_classify(poly_rectify(octahedron()), detail=1);
// family id chosen from cls inspection
q = poly_snub(
    poly_rectify(octahedron()),
    params_overrides=[
        ["face", "family", 0, ["df", 0.06], ["angle", 12]],
        ["vert", "all", ["c", 0.04]]
    ]
);
```

## Tests Covering Snub

Snub behavior is currently pinned in `src/tests/core/TestTruncation.scad`, including:

- combinatorics/counts (`cube`, `dodeca`),
- twist actually moving geometry,
- default edge-triangle quality checks,
- auto-default ranges,
- fixed-`c` angle solve sanity,
- structured override precedence,
- canonical edge-face orientation regression.

## Practical Caveats

- Auto defaults are best on regular bases.
- Multi-family bases may need `params_overrides` (especially `df/angle` by face family) for good planarity/shape quality.
- The solver remains search-based and can be expensive on complex bases.
- Current objective balances can trade off global edge uniformity vs triangle quality; tests reflect the current compromise.

## Performance Notes

- Default solves are fast enough for basic regular cases but still search-based.
- If you are generating many variants, prefer:
  1. solve once (`ps_snub_default_params_overrides`),
  2. reuse rows in repeated `poly_snub` calls.

## Re-entry Notes

- If defaults regress toward `angle=0`, first verify objective inputs:
  - are the same diagonal splits scored that are emitted?
  - is edge face-pair ordering canonicalized?
- Validate against tests before changing search granularity.
- Prefer fixing objective correctness before tuning step counts.

