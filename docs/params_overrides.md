# Params Overrides Notes

This document describes the shared `params_overrides` mechanism and helper functions in `core/params.scad`.

It is intentionally operator-agnostic, so the same mechanism can be reused by truncation/cantellation/chamfer/snub/cantitruncation and future operators.

## Row Schema

Rows can target all elements, families, or explicit element indices.

### Family rows (first-class mode)

```scad
["face"|"vert"|"edge", "family", family_id, ["key", value], ...]
```

### Additional selector scopes

```scad
["face"|"vert"|"edge", "all", ["key", value], ...]
["face"|"vert"|"edge", "id", id_or_ids, ["key", value], ...]
```

- `kind` is usually `"face"`, `"vert"`, `"edge"` (operators may define others)
- `id_or_ids` can be a single index or a list of indices
- keys are operator-defined (e.g. `df`, `angle`, `c`, `de`, `t`)
- each row can set **multiple keys** for the same target

Example with multiple params per row:

```scad
params_overrides = [
    ["face", "all",    ["angle", 15], ["df", 0.03]],
    ["face", "family", 1, ["angle", 20], ["df", 0.04]],
    ["vert", "family", 0, ["c", 0.06], ["de", 0.05]],
    ["face", "id", [2, 7], ["angle", 19]]
];
```

## Override Precedence

For a given `(kind, element_id, family_id, key)`, precedence is fixed:

1. `id`
2. `family`
3. `all`

Within a scope, later rows win; within a row, later duplicate keys win.

## Helper API (`core/params.scad`)

Generic helpers:

- `ps_params_get(params_overrides, kind, key, element_id=undef, family_id=undef)`
- `ps_params_count_kind(params_overrides, kind)`
- `ps_params_compile_key(params_overrides, kind, key, count, family_ids=undef)`
- `ps_params_compile_specs(params_overrides, specs)`
- `ps_compiled_param_get(arr, idx)`

## Compile Specs

When performance matters, compile sparse rows into dense arrays once.

Spec row formats:

```scad
[kind, key, count]
[kind, key, count, family_ids]  // optional parallel array of family ids by element idx
```

Important:
- Family-scoped rows require a known family id at lookup time.
- In compile mode, that means you must provide `family_ids` for that spec if you expect family rows to match.
- If `family_ids` is omitted, family rows do not match compiled lookups (only `all` and `id` rows apply).

Example (per-element compile with family-aware matching):

```scad
params_compiled = ps_params_compile_specs(params_overrides, [
    ["face", "df", len(faces0), face_fid],
    ["face", "angle", len(faces0), face_fid],
    ["vert", "c", len(verts0), vert_fid],
    ["vert", "de", len(verts0), vert_fid]
]);

face_df_by_idx    = params_compiled[0];
face_angle_by_idx = params_compiled[1];
vert_c_by_idx     = params_compiled[2];
vert_de_by_idx    = params_compiled[3];
```

Then lookup by element idx in O(1):

```scad
df_ov = ps_compiled_param_get(face_df_by_idx, face_idx);
```

If you are not compiling and call `ps_params_get(...)` directly, family matching still works if you pass `family_id` explicitly.

## Usage Pattern

1. Classify once (`poly_classify`) if family selectors are used.
2. Compile requested keys once via `ps_params_compile_specs`.
3. During geometry loops, use `ps_compiled_param_get` and operator fallback precedence:
   - selector overrides (`id/family/all`)
   - scalar/operator explicit args
   - derived auto/default values

This keeps operator code fast and keeps parameterization semantics consistent across operators.

## Classification Context Consistency

Family selectors depend on how classification is computed. If two stages use
different classify settings, family ids may differ.

Recommended pattern:

1. Compute `cls = poly_classify(...)` once.
2. Derive any `*_fid` arrays and compiled params from that same `cls`.
3. Pass the same `cls` into placement modules (`classify=cls`) and any
   transform logic that consumes family ids.

This avoids family-id drift between overrides, transforms, and placement.

## Example File

Runnable example:

- [`src/polysymmetrica/examples/basics/main_params.scad`](../src/polysymmetrica/examples/basics/main_params.scad)

It demonstrates:
- `face_fid` as family-id-by-face-index
- compiled dense arrays for `df` and `angle`
- `ps_params_print(...)` output

## Snub Defaults As Structured Params

Snub now has a structured default solver output:

- `ps_snub_default_params_overrides(poly, handedness=1, eps=1e-9)`

This returns override rows directly (face/vert scopes), so it can be passed
straight into `poly_snub`:

```scad
rows = ps_snub_default_params_overrides(hexahedron());
q = poly_snub(hexahedron(), params_overrides=rows);
```
