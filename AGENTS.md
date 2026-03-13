# Repository Guidelines

## Project Structure & Module Organization
- `src/polysymmetrica/core/`: core math, placement, duals, truncation utilities.
- `src/polysymmetrica/models/`: base polyhedra (tetrahedron, octahedron, icosahedron) and derived solids.
- `src/polysymmetrica/examples/`: runnable OpenSCAD examples (basics, truncation, poly-frame).
- `src/tests/`, `src/tests/core/`, and `src/tests/negative/`: OpenSCAD unit tests and runner files.
- `docs/`: developer guide and images used in documentation.

## Build, Test, and Development Commands
This repo is OpenSCAD-first; there is no separate build system.
- Prefer `openscad-nightly` for command-line renders/tests when available (significantly faster F6/CGAL path on this machine).
- Render an example locally:
  `openscad -o /tmp/ps-preview.stl src/polysymmetrica/examples/basics/main-basics.scad`
- Run the full test suite (check console for PASS):
  `openscad -o /tmp/ps-tests.stl src/tests/run_all.scad`
- Negative test (expects a failure):
  `openscad -o /tmp/ps-neg.stl src/tests/run_negative.scad`
- Run all negative tests (expects each file in `src/tests/negative/` to fail):
  `src/tests/run_negative_all.sh`
- Scratch/probe `.scad` files should be created in `/tmp` (for example `/tmp/tmp_probe.scad`), not in the repo root.
- Generated outputs (`.stl`, logs, screenshots) should also go to `/tmp` unless they are intentional docs/examples assets.
- Process safety: never kill `openscad-nightly` broadly (`pkill openscad*` etc.). The user keeps an interactive `openscad-nightly` session running.
  Only target exact non-interactive command lines started for this task.

## Coding Style & Naming Conventions
- Indentation: 4 spaces; keep blank lines between logical blocks.
- Functions/modules: `snake_case` (e.g., `poly_dual`, `place_on_faces`).
- Files: lower snake case (e.g., `placement.scad`, `test_duals.scad`).
- Poly descriptor shape is `[verts, faces, e_over_ir]`; use the accessor helpers in `core/funcs.scad`.
- Preserve `$ps_*` variable naming conventions for placement metadata (see `docs/developer_guide.md`).

## Testing Guidelines
- Tests are plain OpenSCAD modules that call `assert(...)` in `src/tests/core/`.
- Add new test modules with a `test_*` prefix and register them in `src/tests/run_all.scad`.
- Keep numeric tolerances explicit (see `EPS` in `TestFuncs.scad`).
- Do not rely on `poly_valid(...)` alone for transform correctness; add operation-specific tests (counts, adjacency, family behavior) for new geometry operators.

## Scaling & Dual Alignment Notes
- Dual overlays are sensitive to which geometric feature you choose to align (edge midpoints vs face families).
- `scale_dual()` aligns edge-mid spheres; for edge-crossing alignment of a specific edge family, use `scale_dual_edge_cross(poly, dual, face_idx, edge_pos)`.
- For non-uniform truncations, no single scale aligns every edge family; pick the family you want to prioritize.

## Validation Modes
- `poly_valid(poly, "struct")`: structural + planarity only.
- `poly_valid(poly, "closed")`: adds manifoldness + no self-intersections.
- `poly_valid(poly, "star_ok")`: allows self-intersections but keeps manifoldness.
- `poly_valid(poly, "convex")`: adds outward orientation + convexity.

## Commit & Pull Request Guidelines
- Commit messages are short, imperative, and capitalized (e.g., "Consolidate utility funcs").
- PRs should include a brief summary, affected `.scad` paths, and screenshots or renders for geometry changes.
- Note any new tests or expected output changes in the PR description.

## Session Notes (Cantellation Debugging)
- `poly_cantellate` is sensitive to how face, edge, and vertex faces are ordered. Visual crossings can occur even when `assert_poly_valid` passes.
- The test failure we saw was duplicate indices in vertex faces; it was fixed by constructing vertex faces as intersections of the vertex plane and the two incident edge planes per face (see `poly_cantellate` in `src/polysymmetrica/core/truncation.scad`).
- Debugging helpers live in `src/polysymmetrica/examples/truncation/debug_cant_vertex.scad` and can isolate face families or print edge/face mappings.
- A generic validity suite is not enough for cantellation; operation-specific tests (face-family counts, edge-face adjacency) should be added once design stabilizes.
- Cantellation operators now live in `src/polysymmetrica/core/truncation.scad`; solver helpers live in `src/polysymmetrica/core/solvers.scad`.
- **Orientation now uses OpenSCAD LHR** (CW from outside). Models are LHR; `ps_face_normal` and `poly_face_ez` follow LHR. `poly_render` no longer flips faces.
- `face_plate` expects LHR/CW 2D input order for `pts` and aligned `diheds`; bevel sign must match LHR.
- `poly_chamfer` now builds hex edge faces (true chamfer): edge faces include original vertices; vertex faces are omitted.
- Shared mesh build helper: `_ps_poly_from_face_points(...)` dedups points, orients faces, and rescales to unit edge.
- Placement modules now accept optional classification context:
  - `place_on_faces/edges/vertices(..., classify=cls, classify_opts=[detail,eps,radius,include_geom])`
  - New `$ps_*` family vars are exposed per element (`*_family_id`) plus global family counts.
  - To avoid family-id drift between placement and overrides/solvers, classify once and reuse the same `cls`.
- Dependency hygiene: `faces_around_vertex` helpers now live in `core/funcs.scad` (shared primitive). Avoid pulling them from `duals.scad` to prevent `placement -> classify -> duals -> placement` style use-cycles.

## Session Notes (Recent Cleanup Insights)
- Prefer shared scalar helpers in `funcs.scad` when used across core files (for example `ps_clamp(...)`), rather than duplicating private variants per file.
- Remove thin pass-through wrappers when they add no semantic value; call the canonical helper directly.
- For inert cleanup passes, include comment-only/doc-only normalization together with dead-local/dead-helper removal, then always run:
  `openscad -o /tmp/ps-tests.stl src/tests/run_all.scad`
- `classify.scad` now uses `_ps_*_keys_from(...)` forms directly; legacy wrapper variants were removed as dead code.
- Keep debug/probe and generated artifacts in `/tmp`; do not leave temporary `.scad` probes in repo root.
- `poly_attach(...)` lives in `core/construction.scad` and requires `poly_valid(..., "closed")` for both inputs.
  It aligns two selected planar faces (same arity), drops seam faces, then always seam-merges via `poly_cleanup(...)`.
  Use `rotate_step` for cyclic vertex correspondence and `scale_mode="fit_edge"` when input face sizes differ.
  `f1` now accepts either a scalar face index or a list (`f1=[...]`) to attach one copy of `p2` per listed face in a single pass.
  Attach mapping now defaults to chirality-preserving orientation (`mirror=false`); legacy reflected behavior is opt-in via `mirror=true`.

## Session Notes (Construction Layer)
- Construction primitives now live in `src/polysymmetrica/core/construction.scad`, including:
  - `poly_delete_faces(...)`
  - `poly_boundary_loops(...)`
  - `poly_cap_loops(...)`
  - `poly_slice(...)`
  - `poly_attach(...)`
  - `poly_pyramid(...)`
- `src/polysymmetrica/core/attach.scad` was deleted; construction is now the single home for attach/slice/cap workflows.
- Generic `{n,p}` polygon/polygram helpers were hoisted into `core/funcs.scad` so prisms/antiprisms and future Johnson constructors share one geometry basis.
- `poly_pyramid(n, p, edge, height=undef, height_scale=1)` defaults to the regular-equal-edge height when `height` is omitted; this now backs exact `j1_square_pyramid()` and `j2_pentagonal_pyramid()`.
- `poly_cupola(n, edge, height=undef, height_scale=1)` now provides exact J3/J4/J5 cupolae from a direct concentric-polygon construction; `j3_triangular_cupola()`, `j4_square_cupola()`, and `j5_pentagonal_cupola()` are wrappers onto it.
- `poly_rotunda(edge=1)` now provides exact J6 by slicing `icosidodecahedron()` through the origin using a pentagon-face normal and capping the decagonal cut.
- `poly_elongate(...)` and `poly_gyroelongate(...)` are thin construction wrappers over `poly_attach(...)` plus `poly_prism(...)` / `poly_antiprism(...)`; current exact exemplars are elongated cupolae and a gyroelongated triangular cupola.
- `models/johnsons_all.scad` now exports `johnsons_all()` as `[name, fn]` like the other aggregate model files, though the set is still mixed exact/approximate/WIP.
- `src/polysymmetrica/examples/basics/main_johnsons.scad` is the current runnable Johnson/construction demo surface; keep new direct constructors visible there as they are added.

## Session Notes (Non-Planar Face Frames)
- `place_on_faces(...)` now uses a frame normal intended for placement (`ps_face_frame_normal(...)`) rather than relying only on the first triangle normal.
- `ps_face_frame_normal(...)` uses a Newell-style best-fit normal for non-planar faces, then aligns sign to `ps_face_normal(...)` so winding/orientation semantics remain consistent.
- `ps_face_normal(...)` is still the topological/orientation normal and should remain unchanged for validation/duals logic.
- For non-planar faces, `$ps_poly_center_local` may legitimately have significant local X/Y components; this is expected and indicates the face center is not radially aligned to poly center.
- `poly_face_ex(...)` must be projected onto the local face plane (perpendicular to face EZ) before normalization; otherwise local-dot projections drift and `$ps_poly_center_local` placement can be wrong.
- Coverage added in `src/tests/core/TestFuncs.scad`:
  - `test_face_frame_normal__planar_matches_face_normal`
  - `test_face_frame_normal__nonplanar_is_unit_and_oriented`
  - `test_poly_face_ez__uses_frame_normal_for_nonplanar_face`

## Session Notes (Printing Face Segmentation)
- `face_plate_visible(...)` should only use the segmented visible-cell path when `ps_face_geom_cut_entries(...)` returns actual geometry cuts; otherwise it must fall back to plain `face_plate(...)` so regular/star faces keep their known-good bevel and pillow behavior.
- Cutter/cut-entry logic must thread the same fill mode (`"nonzero"` vs `"evenodd"`) that the face geometry uses; otherwise star/self-intersecting cutters silently segment against the wrong filled region.
- `ps_face_visible_segments(...)` must reorient kept cells to match the parent face winding before handing them to bevel code; otherwise parent edges bevel outward.
- For segmented printable pieces, keep original parent edges beveled and let cut edges use cut-derived metadata; do not clip finished 3D plates with `intersection()`, because the normalized CSG tree explodes badly for printing demos.

## Session Notes (Construction Toolkit)
- Johnson-style construction should build on explicit topology tools first: delete faces, recover boundary loops, then cap them. Keep “repair” semantics visible instead of burying them in magical element-deletion helpers.
- `core/construction.scad` now starts with:
  - `poly_delete_faces(...)`
  - `poly_boundary_loops(...)`
  - `poly_cap_loops(...)`
  - `poly_slice(...)`
  - `poly_attach(...)`
- Boundary detection there uses undirected edge multiplicity (`ps_face_has_edge(...) == 1`) and then keeps the directed occurrence from the surviving face. It is intentionally simple and reliable rather than optimized.
- `poly_attach(...)` now lives in `core/construction.scad`; do not reintroduce a separate `core/attach.scad` split unless there is a very strong reason.

## Session Notes (Printing / Visible Face Pieces)
- `face_plate.scad` now handles self-intersecting faces by first segmenting the 2D loop into simple loops and then building body/roof/clearance/pillow from those loops. Do not regress to lofting or hulling the raw self-intersecting loop.
- For star/self-intersecting pillows, `hull()` convexifies the shape and is wrong. Use topology-preserving segmented loops plus loft/extrude instead.
- `segments.scad` now contains a direct “keep visible face cells” path:
  - `ps_face_visible_segments(...)`
  - `place_on_face_visible_segments(...)`
  - `face_visible_mask(...)`
- This path works by splitting a face by geometry-derived cut segments, then discarding cells occluded by other local triangles. It is preferable to `hull() face_cut_stencil(...)` for printing because it avoids magic-number cutters and over-subtraction.
- When splitting a simple face by cut segments, do not require child cell winding to match the outer boundary. The traversal can yield mixed orientation for valid cells.
- Endpoint-on-boundary intersections matter for cut segmentation. A strict interior/interior-only segment intersection test will miss essential split nodes.

## Pre-PR Inert Cleanup Checklist
- Confirm the work is functionally inert (no geometry/output-intent changes).
- Remove clearly unused locals/helpers/wrappers only when references are zero.
- Prefer shared helpers in `funcs.scad` over duplicate per-file utilities.
- Keep comments/docs aligned with actual behavior and parameter semantics.
- Run full tests and verify PASS:
  `openscad -o /tmp/ps-tests.stl src/tests/run_all.scad`
- Ensure scratch/probe files are in `/tmp`, not in repo tree.
- Refresh this `AGENTS.md` with any useful insights from the work unit.
