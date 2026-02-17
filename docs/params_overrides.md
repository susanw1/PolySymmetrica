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

## Usage Pattern

1. Classify once (`poly_classify`) if family selectors are used.
2. Compile requested keys once via `ps_params_compile_specs`.
3. During geometry loops, use `ps_compiled_param_get` and operator fallback precedence:
   - selector overrides (`id/family/all`)
   - scalar/operator explicit args
   - derived auto/default values

This keeps operator code fast and keeps parameterization semantics consistent across operators.
