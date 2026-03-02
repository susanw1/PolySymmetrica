# Repository Guidelines

## Project Structure & Module Organization
- `src/polysymmetrica/core/`: core math, placement, duals, truncation utilities.
- `src/polysymmetrica/models/`: base polyhedra (tetrahedron, octahedron, icosahedron) and derived solids.
- `src/polysymmetrica/examples/`: runnable OpenSCAD examples (basics, truncation, poly-frame).
- `src/tests/` and `src/tests/core/`: OpenSCAD unit tests and test runner files.
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

## Pre-PR Inert Cleanup Checklist
- Confirm the work is functionally inert (no geometry/output-intent changes).
- Remove clearly unused locals/helpers/wrappers only when references are zero.
- Prefer shared helpers in `funcs.scad` over duplicate per-file utilities.
- Keep comments/docs aligned with actual behavior and parameter semantics.
- Run full tests and verify PASS:
  `openscad -o /tmp/ps-tests.stl src/tests/run_all.scad`
- Ensure scratch/probe files are in `/tmp`, not in repo tree.
- Refresh this `AGENTS.md` with any useful insights from the work unit.
