# Iterator APIs

This note catalogs the **iterator-style modules** in PolySymmetrica: modules
that iterate over a set of derived entities and expose `$ps_*` context for each
item.

The intended design rule is:

1. a **function** builds the entity list
2. a thin **iterator module** loops over that list
3. the module sets placement/context variables and calls `children()`

That makes the entities first-class API objects instead of burying the logic
inside nested modules.

This note is about **iterator/traversal surfaces** only. It does **not** try to
catalog plain CSG emitters such as `ps_clip_*`, `ps_partition_*`, or local
debug-only modules.

## API Pattern

Prefer this shape:

```scad
function ps_some_entities(...) = [record0, record1, ...];

module place_on_some_entities(...) {
    ents = ps_some_entities(...);
    for (ei = [0:1:len(ents)-1]) {
        ent = ents[ei];
        // set $ps_* vars
        children();
    }
}
```

Guidelines:

- the function is the canonical source of truth
- the module should stay thin and mostly do record unpacking / frame setup
- if a record shape is nontrivial, add accessors instead of relying on raw
  positional indices everywhere
- do not force this pattern onto modules that are really just geometry emitters

## Current Families

### 1. Global Site Iterators

These are core API surfaces already, but they are still **module-first** today.
They do not yet have matching public “site record” functions.

| Iterator module | Entity type | Backing function today | Status | Likely future data API |
| --- | --- | --- | --- | --- |
| `place_on_faces(...)` | poly face placement sites | none | module-first | `ps_face_sites(poly, edge_len=..., classify=...)` |
| `place_on_edges(...)` | poly edge placement sites | none | module-first | `ps_edge_sites(poly, edge_len=..., classify=...)` |
| `place_on_vertices(...)` | poly vertex placement sites | none | module-first | `ps_vertex_sites(poly, edge_len=..., classify=...)` |

These are the biggest remaining gaps if iterator modules are to become a more
explicit part of the public API.

### 2. Construction / Topology Iterators

These are already on solid ground: the **function** is the public API, and no
iterator module is strictly required unless a future consumer wants one.

| Data function | Entity type | Notes |
| --- | --- | --- |
| `poly_boundary_loops(poly)` | boundary vertex loops of an open mesh | already first-class API |
| `poly_cap_loops(poly, ...)` | cap workflow over boundary loops | consumes loops rather than iterating with a nested module |

This is the cleanest existing example of the “entities first, module second”
pattern.

### 3. Face-Local Segment / Cell Iterators

This family is already mostly in the desired shape.

| Iterator module | Backing function | Entity type | Status |
| --- | --- | --- | --- |
| `place_on_face_segments(...)` | `ps_face_segments(...)` | simple face cells under chosen fill rule | good |
| `place_on_face_geom_cut_segments(...)` | `ps_face_geom_cut_segments(...)` | geometry-derived cut segments | good |
| `place_on_face_visible_segments(...)` | `ps_face_visible_segments(...)` | retained visible face cells | good |
| `place_on_face_filled_boundary_segments(...)` | `ps_face_filled_boundary_segments(...)` | true filled boundary subsegments | good |
| `place_on_face_filled_boundary_edges(...)` | `ps_face_filled_boundary_segments(...)` | same subsegments, but in edge-style bisector frames | acceptable, but record type is implicit |
| `place_on_face_filled_boundary_source_edges(...)` | `ps_face_filled_boundary_source_edges(...)` | unique source edges that still contribute true boundary spans | good |

Related data helpers in the same family:

| Data function | Entity type |
| --- | --- |
| `ps_face_filled_cells(...)` | non-zero filled arrangement cells |
| `ps_face_filled_atoms(...)` | convex atoms derived from filled cells |

The main remaining cleanup target here is:

- `place_on_face_filled_boundary_edges(...)` still derives its edge-style
  records implicitly inside the module rather than iterating an explicit
  edge-record function

## What Should *Not* Be Folded Into This Family

These modules may use iterator-style helpers internally, but they are **not**
themselves iterator APIs:

- `core/face_regions.scad` volume/cutter emitters such as
  `ps_face_interference_volume_ctx(...)`
- `core/proxy_interaction.scad` CSG/context wrappers such as
  `ps_clip_*`, `ps_partition_*`, and `_ps_proxy_*influence(...)`

They should prefer function helpers internally where sensible, but they do not
need the same “record builder + traversal module” treatment as true iterator
surfaces.

## Recommended Next Steps

In order of payoff:

1. add explicit site-record builders for global placement:
   - `ps_face_sites(...)`
   - `ps_edge_sites(...)`
   - `ps_vertex_sites(...)`
2. add lightweight accessors for any record shape that proves durable enough to
   be reused outside its original file

That would make the iterator family much easier to document, test, and reuse.
