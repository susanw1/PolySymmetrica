// Shared parameter override helpers for transform operators.
//
// INPUT FORMAT (`params_overrides`)
// ---------------------------------
// `params_overrides` is a list of rows. Each row targets an element kind
// ("face" | "vert" | "edge"), a selector scope, and one or more key/value
// pairs.
//
// Row schemas:
//   ["face"|"vert"|"edge", "all", ["key", value], ...]
//   ["face"|"vert"|"edge", "family", family_id, ["key", value], ...]
//   ["face"|"vert"|"edge", "id", id_or_ids, ["key", value], ...]
//
// Where:
// - `id_or_ids` may be a single index or a list of indices.
// - keys are operator-defined (examples: "df", "angle", "c", "de", "t").
//
// Precedence is fixed and deterministic:
//   id > family > all
// For the same scope/key, later rows win.
// For duplicate keys in one row, the last key/value wins.
//
// EXAMPLE INPUT
// -------------
// rows = [
//   ["face", "all", ["angle", 15]],
//   ["face", "family", 1, ["df", 0.04]],
//   ["face", "id", [2, 7], ["angle", 19]],
//   ["vert", "family", 0, ["c", 0.06]]
// ];
//
// OUTPUT FORMAT (compiled arrays)
// --------------------------------
// `ps_params_compile_key(...)` returns a dense per-index array of length
// `count`, containing either a resolved value or `undef`.
//
// `ps_params_compile_specs(...)` returns a list of dense arrays, one per spec,
// in the same order as `specs`.
//
// Example:
// face_fid = [0, 1, 1, 0];
// face_df = ps_params_compile_key(rows, "face", "df", 4, face_fid);
// // => [undef, 0.04, 0.04, undef]
//
// face_angle = ps_params_compile_key(rows, "face", "angle", 4, face_fid);
// // => [15, 15, 19, 15]
//
// specs = [
//   ["face", "df", 4, face_fid],
//   ["face", "angle", 4, face_fid]
// ];
// compiled = ps_params_compile_specs(rows, specs);
// // => [face_df, face_angle]
//
// ACCESSOR USAGE
// --------------
// v = ps_compiled_param_get(compiled[0], face_idx); // safe bounds + undef
//
// DEBUG PRINT
// -----------
// Use `ps_params_print(rows);` to dump normalized row interpretation.

function _ps_params_row_kind(row) =
    (!is_list(row) || len(row) < 1) ? undef : row[0];

function _ps_params_row_scope(row) =
    (!is_list(row) || len(row) < 2) ? undef
        : (is_string(row[1]) && (row[1] == "all" || row[1] == "family" || row[1] == "id")) ? row[1]
        : undef;

function _ps_params_row_selector(row) =
    let(scope = _ps_params_row_scope(row))
    (scope == "all") ? undef
        : (scope == "family") ? row[2]
        : (scope == "id") ? row[2]
        : undef;

function _ps_params_row_kv_start(row) =
    let(scope = _ps_params_row_scope(row))
    (scope == "family" || scope == "id") ? 3
        : (scope == "all") ? 2
        : 999999;

function _ps_params_selector_has_id(selector, id) =
    is_undef(id) ? false
        : is_list(selector)
            ? (len([for (x = selector) if (x == id) 1]) > 0)
            : (selector == id);

function _ps_params_row_applies(row, kind, scope, element_id=undef, family_id=undef) =
    let(
        row_kind = _ps_params_row_kind(row),
        row_scope = _ps_params_row_scope(row),
        selector = _ps_params_row_selector(row)
    )
    (row_kind == kind) && (row_scope == scope) && (
        (scope == "all")
        || (scope == "family" && !is_undef(family_id) && selector == family_id)
        || (scope == "id" && _ps_params_selector_has_id(selector, element_id))
    );

function _ps_params_row_get(row, key) =
    let(
        start = _ps_params_row_kv_start(row),
        vals = [
            for (i = [start:1:len(row)-1])
                let(p = row[i])
                if (is_list(p) && len(p) >= 2 && p[0] == key) p[1]
        ],
        vals2 = [for (v = vals) if (!is_undef(v)) v]
    )
    (len(vals2) == 0) ? undef : vals2[len(vals2)-1];

function ps_params_get(params_overrides, kind, key, element_id=undef, family_id=undef) =
    let(
        rows = is_undef(params_overrides) ? [] : params_overrides,
        vals_all = [for (row = rows) if (_ps_params_row_applies(row, kind, "all", element_id, family_id)) _ps_params_row_get(row, key)],
        vals_family = [for (row = rows) if (_ps_params_row_applies(row, kind, "family", element_id, family_id)) _ps_params_row_get(row, key)],
        vals_id = [for (row = rows) if (_ps_params_row_applies(row, kind, "id", element_id, family_id)) _ps_params_row_get(row, key)],
        all2 = [for (v = vals_all) if (!is_undef(v)) v],
        fam2 = [for (v = vals_family) if (!is_undef(v)) v],
        id2 = [for (v = vals_id) if (!is_undef(v)) v],
        v_all = (len(all2) == 0) ? undef : all2[len(all2)-1],
        v_family = (len(fam2) == 0) ? undef : fam2[len(fam2)-1],
        v_id = (len(id2) == 0) ? undef : id2[len(id2)-1]
    )
    !is_undef(v_id) ? v_id
        : !is_undef(v_family) ? v_family
        : v_all;

function ps_params_count_kind(params_overrides, kind) =
    is_undef(params_overrides) ? 0
        : len([for (row = params_overrides) if (_ps_params_row_kind(row) == kind) 1]);

// Compile one kind/key into a dense lookup array, optionally using a parallel
// family-id array for family-scoped matches.
function ps_params_compile_key(params_overrides, kind, key, count, family_ids=undef) =
    [for (idx = [0:1:max(0, count-1)])
        ps_params_get(
            params_overrides,
            kind,
            key,
            idx,
            (is_undef(family_ids) || idx >= len(family_ids)) ? undef : family_ids[idx]
        )
    ];

// Compile multiple specs into dense lookup arrays.
// Spec row format:
//   [kind, key, count]
//   [kind, key, count, family_ids]
function ps_params_compile_specs(params_overrides, specs) =
    [for (spec = specs)
        ps_params_compile_key(
            params_overrides,
            spec[0],
            spec[1],
            spec[2],
            (len(spec) >= 4) ? spec[3] : undef
        )
    ];

function ps_compiled_param_get(arr, idx) =
    (is_undef(arr) || idx < 0 || idx >= len(arr)) ? undef : arr[idx];

module ps_params_print(params_overrides) {
    rows = is_undef(params_overrides) ? [] : params_overrides;
    echo(str("params_overrides: rows=", len(rows)));
    for (ri = [0:1:len(rows)-1]) {
        row = rows[ri];
        kind = _ps_params_row_kind(row);
        scope = _ps_params_row_scope(row);
        sel = _ps_params_row_selector(row);
        start = _ps_params_row_kv_start(row);
        kv = (start > len(row)-1) ? [] : [for (i = [start:1:len(row)-1]) row[i]];
        echo(str("  row#", ri, " kind=", kind, " scope=", scope, " selector=", sel, " kv=", kv));
    }
}
