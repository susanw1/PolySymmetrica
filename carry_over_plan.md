# Carry-Over Plan

This branch is being frozen as an archive/reference point. The goal is to
restart from `main` and bring forward only the parts that are clear, durable,
and worth keeping.

## Carry Over

- Iterator/data-first API direction.
- `place_on_edges(...)` as a dihedral-centered placement surface.
- Useful new `$ps_*` metadata that is simple and durable.
- Face-local segment APIs from `segments.scad`:
  - `ps_face_segments(...)`
  - `ps_face_geom_cut_segments(...)`
  - `ps_face_visible_segments(...)`
  - `ps_face_filled_cells(...)`
  - `ps_face_filled_boundary_segments(...)`
  - `ps_face_filled_boundary_source_edges(...)`
  - iterator wrappers over those
- Documentation improvements:
  - [docs/iterators.md](docs/iterators.md)
  - tighter API comments / naming discipline
  - experiments under `examples/experiments/`

## Drop / Do Not Salvage Whole

- Current live proxy/interference carve stack in `core/proxy_interaction.scad`.
- Current phase-2 lobe machinery as implemented.
- Complex probe/debug scaffolding that only served this branch’s dead ends.
- Any geometry path that cannot be explained simply in terms of:
  - data
  - iterator
  - volume

## First Milestones On New Branch

1. Reintroduce iterator-backed placement/data surfaces cleanly on `main`.
2. Reintroduce face-local segment/boundary APIs with tests on counts, spans,
   and visibility.
3. Build a narrow, trustworthy punch-through probe for
   `poly_antiprism(7,3, angle=15)` before any live printable/carve
   integration.
