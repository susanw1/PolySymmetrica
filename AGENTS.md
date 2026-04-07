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
- `openscad-nightly` is installed via snap on this machine, so it is sandboxed more tightly than the apt `openscad`. When feeding `openscad-nightly` a scratch `.scad`, prefer a repo-local `.tmp/` path over `/tmp`, otherwise snap confinement may block reading the input file.
- This environment has `python3` available but not `python`. For quick probe-file generation or text rewriting, use `python3`.
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
  - `place_on_faces/edges/vertices(..., indices=[...])` filters to an exact site list while reusing the same placement frames and `$ps_*` metadata.
  - `core/proxy_interaction.scad` should build on those exact index lists rather than recreating feature transforms by hand.
  - New `$ps_*` family vars are exposed per element (`*_family_id`) plus global family counts.
  - To avoid family-id drift between placement and overrides/solvers, classify once and reuse the same `cls`.
- Dependency hygiene: `faces_around_vertex` helpers now live in `core/funcs.scad` (shared primitive). Avoid pulling them from `duals.scad` to prevent `placement -> classify -> duals -> placement` style use-cycles.

## Session Notes (Recent Cleanup Insights)
- Prefer shared scalar helpers in `funcs.scad` when used across core files (for example `ps_clamp(...)`), rather than duplicating private variants per file.
- Promote only genuinely generic geometry primitives into `funcs.scad` (for example 2D orientation, convexity, or line intersection); keep face-region-specific half-plane clipping and offset machinery local until it has a second real consumer.
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
- Public Johnson model wrappers in `models/johnsons_all.scad` should use `jNN_*` names consistently; keep descriptive names for generic construction helpers, not the model surface.
- `src/polysymmetrica/examples/basics/main_johnsons.scad` is the current runnable Johnson/construction demo surface; keep new direct constructors visible there as they are added.

## Session Notes (Non-Planar Face Frames)
- `place_on_faces(...)` now uses a frame normal intended for placement (`ps_face_frame_normal(...)`) rather than relying only on the first triangle normal.
- `ps_face_frame_normal(...)` uses a Newell-style best-fit normal for non-planar faces, then aligns sign to `ps_face_normal(...)` so winding/orientation semantics remain consistent.
- `ps_face_normal(...)` is still the topological/orientation normal and should remain unchanged for validation/duals logic.
- `place_on_edges(...)` now uses the adjacent-face normal bisector as edge-local `+Z` when a usable face pair exists; this is the preferred dihedral-centered frame for edge structure/proxy work. Boundary/degenerate edges fall back to the older radial frame.
- For non-planar faces, `$ps_poly_center_local` may legitimately have significant local X/Y components; this is expected and indicates the face center is not radially aligned to poly center.
- `poly_face_ex(...)` must be projected onto the local face plane (perpendicular to face EZ) before normalization; otherwise local-dot projections drift and `$ps_poly_center_local` placement can be wrong.
- Coverage added in `src/tests/core/TestFuncs.scad`:
  - `test_face_frame_normal__planar_matches_face_normal`
  - `test_face_frame_normal__nonplanar_is_unit_and_oriented`
  - `test_poly_face_ez__uses_frame_normal_for_nonplanar_face`

## Session Notes (Face Region Volumes)
- `src/polysymmetrica/core/face_regions.scad` is the new home for local face-region and segmentation volume helpers.
- These helpers are intentionally face-local: they depend on face geometry, cut geometry, and dihedral/2 rules, not on any global polyhedral centre convention.
- Keep `examples/printing/face_plate.scad` thin. If segmented-face or printable join geometry needs more machinery, move that machinery into `core/face_regions.scad`, not back into the example.
- Helper tests now live in `src/tests/core/TestFaceRegions.scad`, not under `src/tests/examples/`.
- Generic child-clipping wrappers now live there too:
  - `ps_clip_to_face_region(_ctx)`
  - `ps_clip_to_visible_face_cell_ctx(...)`
  - `ps_clip_to_visible_face_segments_ctx(...)`
- `core/face_regions.scad` now owns the real segmented-face geometry story. `segments.scad` finds visible cells and cut provenance; `face_regions.scad` turns that into admissible 3D regions.
- The supported contract is: derive keep regions from the underlying polyhedron alone, then clip arbitrary child geometry to them. Do not try to infer a perfect result recursively from already-decorated neighboring faces.
- `examples/printing/face_plate.scad` should stay a thin consumer of those core regions. Avoid reintroducing segmentation or join solving into the example layer.
- For segmented visible cells, the important model is direct region construction from **all** side constraints together, not "keep prism minus a union of independent cut bands". This matters when multiple cut edges meet at one corner, for example where a true edge punches through another face.
- Parent edges use the ordinary face dihedral/2 rule. Cut edges use cutter-derived join geometry from the matched cut entries and cutter face normals.
- `face_plate_visible(..., seg_apply_cut_bands=true)` now selects the direct visible-cell region path. `seg_along_pad` and the old subtraction-band tuning were removed as dead baggage.
- Consumer rule: pass the full intended segmentation join clearance into the core region path. Do not halve `edge_inset` before passing it in.
- The exact 2D clipped loop helper (`ps_face_visible_cell_loop_at_z_clipped(...)`) is now trusted as cross-section truth for bad punch-through cases such as `poly_antiprism(7,3, angle=15)` `f2/c0`.
- Important lesson: a visible cell may still be nonconvex. In that case the clipped half-plane loop is the convex kernel of the cell, not the whole visible region. Do not replace a nonconvex live cell loop with the clipped loop directly.
- The current live fix is: decompose nonconvex visible cells into convex atoms in `face_regions.scad`, build one atom region per convex atom, then union them.
- Important refinement: a convex atom created by triangulating a nonconvex cell is not automatically a valid half-plane-bounded region for clipped-loop lofting. If an atom contains any triangulation `"inner"` edge, the clipped-loop path can blow up far beyond the true atom footprint. For now, only atoms bounded entirely by real cell edges may use clipped-loop lofting; inner-edge atoms must stay on the adjacent-line path.
- Do not switch live geometry over just because the clipped top loop looks right. The acceptance rule for future exact-region work is stricter: an exact convex-atom region model must reproduce the accepted local cross-sections across the actual occupied `z` span before the live region builder is replaced.
- The raw exact 3D CSG probe/debug path was removed as a dead end. Keep the useful diagnostics only:
  - face/cell/cut labels
  - orange-vs-cyan loop comparison
  - plane-derived 2D cross-section checks
- For the canonical troublesome `poly_antiprism(7,3, angle=15)` join, the matched cut pair `f2/c0/e1/c1/f5` and `f5/c1/e1/c0/f2` already agree numerically:
  - `cut_dihed ~= 94.2341`
  - `join_dihed ~= 265.7659`
  - identical profile over the face slab span `[[0.55, -0.4], [2.03592, 1.2]]`
  So if that join still looks wrong, the bug is in whole-cell region construction, not in the per-edge dihedral math.

## Session Notes (Printing Segmentation Metadata)
- `face_plate_visible(...)` uses `ps_clip_to_visible_face_segments_ctx(...)` as its segmented path and should remain a thin consumer of the core region model.
- `ps_face_visible_segments(...)` now returns `cell_cut_entry_ids` as a fifth parallel vector alongside `cell_edge_kinds`. Values are `undef` for parent edges and stable indices into `ps_face_geom_cut_entries(...)` for cut edges.
- `place_on_face_visible_segments(...)` exposes the same data as `$ps_vis_seg_cut_entry_ids`; `place_on_face_segments(...)` exposes `$ps_seg_cut_entry_ids`.
- `cut_entry_id` is face-local; use it to index back into `ps_face_geom_cut_entries(...)`.
- `cut_pair_id` is the world-stable join id for the same geometric cut seen from two faces. It is derived from the two face ids plus a quantized world-space line signature (anchor point + direction), and is exposed as `$ps_vis_seg_cut_pair_ids` when `place_on_faces(...)` world-frame context is available.
- `cut_run_id` is local to one visible-cell boundary and groups one continuous cutter contribution. Use it to distinguish split cut spans when the same cutter re-enters the same visible cell in multiple disjoint runs.
- Core region semantics: `z0`/`z1` are only admissible-region bounds. Parent/cut side faces should be full-depth planes across that span. If a cut span is finite, the principled fix is full-depth run-end planes, not top-only notches or cap-specific heuristics.
- Related rule: the target model is that the live cut-edge path should use the full active `z0..z1` span for side-plane construction. Treat any internal `band_z0` / `band_z1` as suspect legacy/debug plumbing rather than a real top/bottom distinction in core geometry. A direct swap to full-span live construction was attempted and backed out because it caused broader region-builder regressions.
- `ps_face_visible_cell_cut_run_end_entries(...)` is the first principled helper for finite split cut spans. Keep it as analysis/debug metadata until a full-depth run-end-plane live path is ready; do not reintroduce heuristic top-only endpoint notches.
- Keep `_ps_fr_atom_is_enabled(...)` minimal and justified: all-inner atoms stay off, explicit debug skip lists are honored, but avoid heuristic atom-role suppressions in the live path unless they have a general geometric proof.
- The key rule for future cut-edge relief work is: use these propagated cut-entry ids, not fuzzy segment-equality matching, when tying visible-cell edges back to cutter geometry.
- A plausible future fabrication path is **proxy interaction**: keep the analytic layer for visibility/provenance, but let user-supplied raw face/edge/vertex proxy solids drive interference subtraction with ordinary CSG. If pursued, proxies must be raw proposals only, bounded by explicit local influence extents, and must not recurse on already-clipped neighbors.
- For face-face interference in that proxy path, a useful early mode is `face_proxy_mode="sweep_to_bounds"`: project the neighboring face proxy to its local XY footprint and extrude it through `face_bounds` to form a bounded subtractor volume.
- Signature sketches for that proxy path live in `src/polysymmetrica/core/proxy_interaction.scad`; keep that file at the API-contract level until there is a clear first implementation.
- Proxy interaction needs a hard semantic split:
  - **occupancy proxies** describe material that really exists and should cut foreign intersecting polys
  - **clearance proxies** describe intentional fit/seat/tolerance gaps and should stay separate by default
- Do not reuse local face/edge seating clearance as the foreign inter-poly cutter. Inter-poly subtraction should be driven by foreign occupancy, with any desired extra fit gap represented explicitly as a separate clearance/dilation step.
- Proxy fabrication should now be framed as:
  - `face_shielded(i) = F_raw(i) ∩ Z(i) ∩ ⋂ B(i,e)` using one adjacent-edge bisector half-space per edge
  - `V(i) = ⋃ Occ(x)` over non-adjacent intersecting features only
  - `face_final(i) = face_shielded(i) - V(i) - C_local(i)`
- Adjacent faces should not be treated as foreign cutters once those bisector shields are active; their anti-interference role is already handled by the shield planes.
- Local edge-clearance strips still belong on the simple indexed `place_on_edges(...)` path in dihedral-centered frames. Do not reshape them with extra corridor/span clipping unless a concrete need survives review.
- For self-intersecting/star faces, do not apply those shields to the original self-crossing walk directly. First split at self-intersections, create pseudo-vertices, classify the filled arrangement, keep only the true filled-boundary subsegments, then apply shielding/clearance from those boundary subsegments.
- In other words: the next principled proxy step is face-local arrangement/boundary extraction for nonconvex/self-crossing faces, not more cutter/clearance tuning.
- Preferred proxy execution order is now:
  1. face-local arrangement extraction
  2. convex atomization of the filled face region
  3. shielded-face construction from true boundary subsegments
  4. non-adjacent `ps_intersections(...)` hit-list generation
  5. foreign occupied-volume builder
  6. final `face_shielded - foreign_occ - local_clearance` carve
  7. only then frame integration and later analytic run-end relief
- One concrete trap already found: do not pass a short `edge_length` influence clip (for example `IR`) into edge proxy subtraction unless you really intend to truncate the strip along its own `x` axis. That was the cause of the “middle part of the edge works, ends missing” regression.
- First landed cell-builder on that proxy baseline: `ps_partition_face_by_feature_proxies(...)` / `_ctx(...)` in `core/proxy_interaction.scad` recursively partition the target face proxy by the actual proxy cutters and emit the resulting cells as separate solids. Keep cutter counts explicit and small (`face_indices`, `edge_indices`, `vertex_indices`, `max_cutters`) because the split is exponential in the number of cutters.
- First landed composed fabrication helper on that same baseline: `ps_carve_face_by_feature_proxies(...)` clips face occupancy and local clearance separately by the same foreign occupancy cutters, then subtracts the clipped local clearance from the clipped face occupancy. Keep lower-level `clip`/`cells` helpers as internal debugging aids only; the example path should use the composed carve.
- Pure cut-edge cross-section helpers for printing live in `examples/printing/face_plate.scad` for now (`ps_face_cut_join_dihed`, `ps_face_cut_relief_u_at_z`, `ps_face_cut_relief_profile2d`); their tests live under `src/tests/examples/`, not `src/tests/core/`.
- Failed geometry experiments should be deleted rather than left around dead. Keep the abstraction boundary clean: segmentation metadata in `segments.scad`, admissible regions in `face_regions.scad`, example styling in `face_plate.scad`.

## Session Notes (Printing Segmentation)
- The generic segmented path now builds visible-cell regions directly from all edge side lines together. This is the right abstraction for punch-through cases and future stellation-style intersections.

## Session Notes (Printing Face Segmentation)
- `face_plate_visible(...)` should only use the segmented visible-cell path when `ps_face_geom_cut_entries(...)` returns actual geometry cuts; otherwise it must fall back to plain `face_plate(...)` so regular/star faces keep their known-good bevel and pillow behavior.
- `examples/printing/face_plate.scad`: keep segmentation join clearance separate from ordinary outer-edge inset. `face_plate_visible(..., seg_cut_clearance=...)` controls the cut-edge gap; it defaults to `edge_inset` only for backward-compatible behavior.
- `core/face_regions.scad`: avoid reintroducing alternate cut-band/subtraction APIs unless they have a real consumer. The main segmented path should stay the direct visible-cell region build.
- `core/segments.scad`: `nonzero` is now the intended default for `ps_face_segments(...)` / `place_on_face_segments(...)`; keep `evenodd` as an explicit/debug option rather than the normal user path.
- Cutter/cut-entry logic must thread the same fill mode (`"nonzero"` vs `"evenodd"`) that the face geometry uses; otherwise star/self-intersecting cutters silently segment against the wrong filled region.
- Same-face cut fragments must merge by connected collinear intervals only. Do not collapse all fragments from one cutter face to the farthest pair, or disjoint spans get bridged across empty gaps.
- Docs split: `docs/segments.md` explains split-cell / cut-line analysis and iteration; `docs/face_regions.md` explains the 3D clipping/volume layer that consumes that data. Keep examples and docs pointing users at the right layer instead of re-explaining segmentation inside `face_plate`-style consumers.
- `ps_face_visible_segments(...)` must reorient kept cells to match the parent face winding before handing them to bevel code; otherwise parent edges bevel outward.
- `ps_face_visible_segments(...)` should only collapse back to the unsplit parent polygon when **all** split cells are visible. If a cut leaves multiple visible survivors, keep them all as separate cells.
- For segmented printable pieces, keep original parent edges beveled and let cut edges use cut-derived metadata from the region model. Avoid solving segmented joins in the example layer.

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
- `make_poly(...)` and `_ps_poly_ir(...)` now recenter geometry by the mean edge-midpoint center before computing/storing `e_over_ir`. This is required for asymmetric construction/transform outputs (e.g. elongated Johnsons, non-uniform truncations) to keep descriptor scaling meaningful.
- Open construction shells are now valid poly descriptors for internal workflows. `make_poly(...)` no longer insists on 4 faces / 6 edges, which lets `poly_slice(..., cap=true)` keep a 3-face intermediate shell alive until capping closes it.
- Construction editing APIs should reject non-integer face IDs up front. Do not silently `round()` computed/parameterized face indices before validating them, or topology can be corrupted without any error.

## Session Notes (Printing / Visible Face Pieces)
- `face_plate.scad` now handles self-intersecting faces by first segmenting the 2D loop into simple loops and then building body/roof/clearance/pillow from those loops. Do not regress to lofting or hulling the raw self-intersecting loop.
- For star/self-intersecting pillows, `hull()` convexifies the shape and is wrong. Use topology-preserving segmented loops plus loft/extrude instead.
- `segments.scad` now contains a direct “keep visible face cells” path:
  - `ps_face_visible_segments(...)`
  - `place_on_face_visible_segments(...)`
  - `face_visible_mask(...)`
- This path works by splitting a face by geometry-derived cut segments, then discarding cells occluded by other local triangles. It is preferable to `hull() face_cut_stencil(...)` for printing because it avoids magic-number cutters and over-subtraction.
- If cut lines split a face but every child cell remains visible, keep the original unsplit face cell.
- Star-prism side faces should resolve to the two visible side panels after spike occlusion; do not collapse that case to a single octagon or center strip.
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
