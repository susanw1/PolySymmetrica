# Segmentation Region Plan

This plan is the working design for segmentation / face-region handling in
Polysymmetrica, especially for stellations and other heavy face-intersection
cases.

## Core Decision

The exact core object is still the **admissible region**, not the clipped mesh.

That gives us a stable boundary:

- `segments.scad` answers:
  - what visible cells exist
  - what edges bound each cell
  - which edges are true parent edges vs cut edges
  - which cutter face produced each cut edge
- `face_regions.scad` answers:
  - what 3D region is legal for that cell between `z_below` and `z_above`

Then any consumer:

- `face_plate`
- future edge frames
- arbitrary `children()`

just gets clipped to that region.

That is the only generic promise we can actually keep.

## Region Model

The region is defined by:

- top slab plane
- bottom slab plane
- one side half-space per parent edge
- one side half-space per cut edge

That is the underlying truth.

The important implementation lesson is that visible cells may still be
nonconvex. So:

- the exact half-plane-clipped loop at one `z` can be the convex kernel of a
  nonconvex visible cell
- directly replacing a live nonconvex cell loop with that clipped loop is wrong

The current live path therefore works like this:

- keep visible-cell topology from `segments.scad`
- if a visible cell is nonconvex, decompose it into convex atoms
- build one admissible region per convex atom
- union the atom regions

The harder future exact-region work is still about replacing the sampled/lofted
convex-atom builder with a more exact one, not about collapsing nonconvex cells
to their kernels.

The sampled path is hard because it still approximates the convex-atom truth by:

- computing 2D loops at a few `z` levels
- lofting them

That fails when:

- multiple cuts meet at one corner
- a true edge punches through another face
- loop arity changes with `z`

which is exactly what stellations will do all the time.

## Plan That Should Work

1. Keep visible cells as true topology; decompose later if needed.

- `segments.scad` should keep returning the real visible cells.
- `face_regions.scad` should decompose nonconvex visible cells into convex atoms
  only when building local regions.

Why:

- the analysis layer should preserve topology
- convex decomposition is a geometry-construction concern, not a visibility one

2. Add stable world-space cut pairing.

- local labels like `e1/c1/f5` are not enough
- add a `join_id` / `cut_pair_id` for each visible cut edge
- derive it from the world-space cut line and the two face ids

Why:

- lets us test "these two neighboring pieces should use the same join plane"
- makes debugging much less ambiguous

3. Represent every cell boundary edge as an exact 3D plane.

For each edge of one visible cell:

- parent edge:
  - plane from face dihedral/2
- cut edge:
  - plane from cutter face / join bisector
  - plus `seg_cut_clearance`

This should be a real 3D plane equation, not a sampled `u(z)` profile.

4. Keep the current convex-atom sampled builder as the stable live path.

That is the working result now.

Future exact-region work can still target exact half-space region construction,
but it should do so per convex atom, not per nonconvex visible cell.

5. Keep the exact 2D helper as the cross-section verifier.

The clipped-loop helper in `face_regions.scad` is still valuable. Use it to:

- verify the cross-section at chosen `z`
- debug the region builder
- write tests

But it should not be the main construction path by itself.

6. Only after an exact convex-atom region path works, consider replacing the
   sampled live path.

Once the exact region path is visually correct for the canonical cases:

- `5/2`
- `7/2`
- `7/3`
- `7/3 angle=15`

then we can decide whether to:

- keep it as the generic path
- or add a faster mesh builder for examples like `face_plate`

The fast path must match the exact path, not replace it speculatively.

## What We Learned From Recent Attempts

1. The topology side is mostly fine.

`segments.scad` is not the blocker now. Visible cells, cut provenance, and
`cut_pair_id` are doing their job.

2. The exact 2D clipped loop is trustworthy as a local cross-section tool.

`ps_face_visible_cell_loop_at_z_clipped(...)` correctly identifies the bad
`f2/c0` corner in the `poly_antiprism(n=7, p=3, angle=15)` case.

3. The active 3D builder was only part of the problem.

The sampled-loop loft in `face_regions.scad` was too weak for punch-through /
multi-cut corners, but the deeper issue was that `f2/c0` itself was nonconvex.

4. "Use clipped loops directly" is not enough.

For a nonconvex visible cell, the clipped loop is the convex kernel. So
replacing the live loop with it collapses real visible area and makes the piece
worse.

5. Consumer `z` range matters a lot.

The live result changed drastically depending on the `z` span fed into the
region builder. The exact region path cannot be validated only at one
cross-section; it has to be checked over the actual occupied `z` span.

6. The raw exact 3D CSG probe was not useful.

The useful debug artifacts were:

- the 2D orange-vs-cyan loop overlay
- plane-derived cross-sections compared against the clipped-loop truth

The abandoned exact 3D CSG probe was removed rather than kept around as
misleading machinery.

## Revised Staging Rule

Do not replace the live convex-atom region builder until an exact convex-atom
region model reproduces the accepted local cross-sections across the occupied
`z` span for the canonical bad cases.

This remains the acceptance gate.

That means we must keep the work separated into:

1. exact 2D cross-section truth
2. convexity diagnosis of the visible cell
3. exact convex-atom region construction
4. consumer integration

The earlier failed attempts jumped too early from 1 to 3.

## Practical Execution Order

1. `segments.scad`

- keep the visible-cell / cut-pair metadata stable
- avoid changing topology unless a real visibility bug is found

2. `face_regions.scad`

- keep the current convex-atom sampled live path until a replacement is proven
- keep exact 2D / plane-derived helpers as acceptance tools
- validate future exact convex-atom builders by slicing at chosen `z` values
  and comparing with the accepted local cross-sections

3. `main_demo`

- keep debug labels
- keep loop-comparison overlays
- use pair-id / cell-id overlays when discussing a problematic join
- do not reintroduce the abandoned exact-volume blob debug

4. Tests

- sample-point tests on exact region membership
- pair consistency tests across matched seg edges
- canonical visual/probe cases:
  - `poly_prism(5,2)`
  - `poly_prism(7,2)`
  - `poly_prism(7,3)`
  - `poly_antiprism(7,3, angle=15)`

## Recommendation

Do not touch `face_plate` next, except as a thin consumer or debug surface.

The next real work unit should be:

- keep `segments.scad` stable
- keep the current live geometry path stable
- keep the convex-atom decomposition
- only pursue future exact-region work per convex atom
- accept any future replacement only when its horizontal slices match the
  accepted local cross-sections across the actual occupied `z` span

That is the route that still looks robust enough for stellations instead of
just patching the current antiprism bug.
