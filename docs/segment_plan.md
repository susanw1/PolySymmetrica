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

## Proxy Fabrication Model

The current proxy design should be described in three distinct layers:

- **anti-interference shield**
  - one bisector half-space per adjacent edge
- **foreign occupied volume**
  - non-adjacent intersecting faces / edges / vertices only
- **local clearance**
  - the target face's own seating / inset subtraction

The important correction is that adjacent faces should not be treated as
foreign cutters once the bisector shields are active. Their role is already
accounted for by the target face's own anti-interference planes.

Using compact notation:

- `F_raw(i)` = supplied face geometry for face `i`
- `Z(i)` = face slab `z0..z1`
- `B(i,e)` = "stay on my side" bisector half-space for edge `e`
- `C_local(i)` = local edge / vertex clearance for face `i`
- `I(i)` = non-adjacent intersecting features for face `i`
- `Occ(x)` = occupied solid of feature `x`

Then the intended target model is:

- `F_shield(i) = F_raw(i) ∩ Z(i) ∩ ⋂[e in boundary(i)] B(i,e)`
- `V(i) = ⋃[x in I(i)] Occ(x)`
- `F_final(i) = F_shield(i) - V(i) - C_local(i)`

OpenSCAD-ish pseudocode:

```scad
module face_shielded(i) {
    intersection() {
        face_raw(i);
        z_slab(i);

        for (e = boundary_edges(i))
            bisector_halfspace(i, e, eps = EPS);
    }
}

module intersecting_vols(i) {
    union() {
        for (x = ps_intersections(i))   // non-adjacent only
            occupied_feature(x);
    }
}

module local_clearance(i) {
    union() {
        for (e = boundary_edges(i))
            edge_clearance(i, e);

        for (v = boundary_vertices(i))
            vertex_clearance(i, v);
    }
}

module face_final(i) {
    difference() {
        face_shielded(i);
        intersecting_vols(i);
        local_clearance(i);
    }
}
```

This is the basis for future proxy work. It deliberately avoids the older,
failed corridor/ownership experiments where extra face-local clipping reshaped
the clearance strips before subtraction.

Current implementation direction:

- target face occupancy should be intersected with a shield volume built from
  arrangement atoms
- adjacent faces should be excluded from the foreign face-cutter set
- local clearance should stay a separate final subtraction, placed on the true
  filled-boundary edge subsegments

### Edge Interference Principle

The edge-side anti-interference cutter should follow the same logic as the
planned face-cell removal path:

- start from a simple default edge cutter in the dihedral-centered boundary-edge
  frame
- then remove the parts of that cutter occupied by crossing face cells

So for one true filled-boundary segment `b` on face `i`:

- `E0(i,b)` = default edge interference cutter
- `owner_cell(i,b)` = the filled cell that owns boundary segment `b`
- `X(i,b)` = union of the other filled face cells of `i` that cross `E0(i,b)`
- `E(i,b) = E0(i,b) - X(i,b)`

For execution, the simpler first implementation is:

- build `owner_cell_prism(i,b)` by extruding the owning filled cell through the
  active face slab
- then define:
  - `E(i,b) = E0(i,b) ∩ owner_cell_prism(i,b)`

Within the face slab, that is equivalent to subtracting the other crossing
cells, but it avoids explicit crossing-edge detection in the first pass.

This is the preferred next execution model for the remaining "armpit hair"
corner artifacts. The correct shape should come from clipping the default edge
cutter by crossing face-cell volumes, not by inventing ad hoc tapered end
geometry.

Consequences:

- for ordinary convex faces, `X(i,b)` is empty, so `E(i,b) = E0(i,b)`
- for star/self-intersecting faces, crossing cells automatically trim the edge
  cutter into the correct corner shape
- the edge cutter should be treated as another cell-removal problem, not as a
  special-case bevel/taper problem
- there is no need to reason about crossing edges explicitly in the first
  implementation; the owning filled cell already encodes the right keep region

### What Proved Hard

The hard part was not "draw a sloping cutter". It was assigning the right
constraints to a nonconvex/self-crossing face without collapsing it or leaving
local seams.

What failed:

- a single global intersection of all boundary-derived bisector half-spaces
  collapsed star/heptagram faces to their convex kernels
- broadening each convex atom to use all parent-cell boundary cutters also
  over-constrained the face and collapsed the star again
- keeping cutters atom-local preserved the arms and fixed the main bevel
  direction, but left the remaining interior-corner "armpit hair" because the
  default rectangular edge cutter was still blind to crossing face cells
- treating the problem as an `eps` or strip-height issue was a dead end; those
  only changed the symptoms

What turned out to matter:

- the boundary-edge frame orientation had to be fixed at the source, not by
  blind sign-flipping of individual cutters
- the face-side interference problem is now good enough in isolation
  (`test_interference.scad`) and in the proxy face branch
- the remaining corner cleanup belongs to the edge cutter, not to more global
  shield heuristics

If the older ideas are revisited later, the real questions are:

- how to distribute boundary constraints across nonconvex filled cells without
  reintroducing kernel collapse
- how to make atom unions meet cleanly without leaving seam artifacts
- and how to do that without turning the model into a pile of mode flags or
  hand-designed taper cases

## Nonconvex / Self-Intersecting Faces

The shielded-face model above needs one important refinement for star and other
self-intersecting faces.

For a self-intersecting face, the naive global form:

- `F_shield(i) = F_raw(i) ∩ Z(i) ∩ ⋂ B(i,e)`

is wrong, because intersecting all edge half-spaces collapses the face to its
convex kernel and destroys the protruding filled arms.

So the next real starting point is not more cutter tuning. It is a face-local
arrangement step:

1. split the original face walk at self-intersections
2. create pseudo-vertices at those crossings
3. classify the filled planar regions using the chosen fill rule
   (`"nonzero"` is the intended normal case)
4. distinguish:
   - **true filled-boundary segments**
   - **internal crossing segments**
5. use only the true filled-boundary segments for:
   - anti-interference shields
   - local edge-clearance placement

That means the useful future helper is not just:

- `ps_face_self_intersectors(...)`

but something closer to:

- `ps_face_filled_arrangement(...)`
- or `ps_face_filled_boundary_segments(...)`

returning subsegments with inherited metadata from their source face edges.

Because those filled regions may still be nonconvex, the final shielded-face
construction should then be:

- decompose the filled face region into convex atoms/cells
- apply the relevant bisector shields to each convex atom
- union the results back together

So for star/self-intersecting faces, the target model becomes:

- `F_shield(i) = ⋃ atom_shield(atom_k)`

rather than one monolithic half-space intersection.

This is the right place to start next. It is more fundamental than any further
clearance/cutter ordering tweaks, because it determines what the real face
boundary even is for shielding and local seating.

Practical notes already established:

- local edge-clearance strips should come from real indexed
  `place_on_edges(...)` placement
- do not pass a short `edge_length` influence bound such as `IR` into that
  path unless you intentionally want to truncate the strip along its own `x`
  axis
- anti-interference is conceptually just the dihedral/2 half-space; it is not
  the same thing as local seat clearance

The first landed constructive/debug step on top of that baseline is still:

- `ps_partition_face_by_feature_proxies(...)` /
  `ps_partition_face_by_feature_proxies_ctx(...)` in
  `core/proxy_interaction.scad`

That path recursively splits one target face proxy by the selected proxy
cutters and emits the resulting cells as separate solids. It is diagnostic
rather than the final fabrication path.

There is now a second, narrower refinement problem too:

- even when the side planes are correct, the **endpoints of a cut span** can
  still misbehave because the region currently lacks finite run-end bounds
- this happens when a cutter face contributes multiple disjoint spans to the
  same visible cell, or when a cut span terminates into another cut / true-edge
  constraint
- the visible symptom can show up near one face boundary first, but there is no
  special top/bottom rule in the core model: the missing geometry is a full-
  depth run-end plane, not a cap-specific tweak

So the region problem is now split into two levels:

- **whole-cell region correctness**
  - solved enough for the current `f2/c0` / `f5/c1` join by convex-atom
    decomposition
- **cut-span endpoint / cut-vertex relief**
  - still unresolved for more complex cases like the `f2` / `f11` repeated-cut
    spans

The sampled path is hard because it still approximates the convex-atom truth by:

- computing 2D loops at a few `z` levels
- lofting them

That fails when:

- multiple cuts meet at one corner
- a true edge punches through another face
- loop arity changes with `z`

which is exactly what stellations will do all the time.

## Proxy Execution Plan

This is the implementation order that should now be followed. The point is to
land one geometrically meaningful layer at a time and avoid mixing adjacency
shielding, foreign cutters, and local clearance in one patch.

### Phase 0: Freeze The Known-Good Baseline

Goal:

- keep the current simple proxy carve path working while new pieces are built
  alongside it

Deliverables:

- preserve the raw indexed `place_on_edges(...)` local-clearance strip path
- keep `test_proxy.scad` able to show a known-good simple face carve
- do not add new corridor/span ownership shaping to the live strip path

Exit criteria:

- one ordinary convex face still carves correctly with local clearance only
- no return of the old "middle works, ends missing" regression

### Phase 1: Face Arrangement Extraction

Goal:

- turn one self-intersecting face walk into a usable filled planar arrangement

Deliverables:

- helper to split original face edges at self-intersections
- pseudo-vertices at crossing points
- filled-region classification using the chosen fill rule (`"nonzero"`)
- boundary-vs-internal segment classification

Suggested helper surface:

- `ps_face_filled_arrangement(...)`
- or `ps_face_filled_boundary_segments(...)`

First landed step:

- `ps_face_filled_cells(...)`
- `ps_face_filled_boundary_segments(...)`
- `ps_face_filled_atoms(...)`
- `place_on_face_filled_boundary_segments(...)`
- `place_on_face_filled_boundary_edges(...)`

These now expose the non-zero filled arrangement as simple cells plus the true
filled-boundary subsegments with inherited source-edge ids/parameters, plus a
conservative convex atomization where atom edges are marked as `"boundary"` or
`"inner"`. The segment-placement helper iterates those true boundary
subsegments in stable local frames, and the boundary-edge helper lifts them
into dihedral-centered edge frames using the inherited source-edge dihedral.
That gives later clearance/shield logic a way to follow the actual filled
perimeter rather than the original self-crossing walk.

Per-segment metadata should include:

- endpoints
- source original edge index
- inherited source dihedral / edge parameters
- boundary/internal classification

Exit criteria:

- a pentagram/star face yields true boundary subsegments for the filled arms
- internal crossing segments are no longer treated as face boundaries
- local edge-clearance placement can be driven from those true boundary
  subsegments instead of the original face walk

### Phase 2: Convex Atomization Of Filled Face Regions

Goal:

- convert the filled planar face result into convex pieces suitable for
  shielding

Deliverables:

- helper to decompose the filled 2D arrangement into convex atoms/cells
- stable mapping from each atom edge back to its source boundary metadata where
  applicable

Why this phase exists:

- even after arrangement extraction, the filled face may still be nonconvex
- one global half-space intersection is still wrong for shield construction

Exit criteria:

- one nonconvex filled face can be represented as `union(convex_atoms)`
- atom boundaries retain enough metadata to know which edges are real boundary
  segments

### Phase 3: Shielded Face Construction

Goal:

- build `F_shield(i)` correctly from the arrangement-derived boundary

Deliverables:

- `face_shielded(i)` built as:
  - face slab clip `Z(i)`
  - per-boundary-segment bisector half-spaces `B(i,e)`
  - applied per convex atom, then unioned

Important rule:

- adjacent-edge anti-interference belongs here
- adjacent faces are not yet in the foreign cutter set

Exit criteria:

- convex faces shield correctly
- star/self-intersecting faces keep their filled arms instead of collapsing to
  the convex kernel

Current status:

- face-side anti-interference now works in a dedicated probe and is wired into
  the proxy face branch
- the remaining unresolved part is the edge-side corner cleanup described in
  the edge-interference principle above

### Phase 4: Intersecting Feature Query

Goal:

- produce the exact non-adjacent hit-list for one target face

Deliverables:

- `ps_intersections(face_i, ...)` or equivalent helper
- returns candidate non-adjacent:
  - faces
  - edges
  - vertices
- bounded by:
  - target face slab thickness
  - maximum edge radius / vertex radius assumptions where needed

Important rule:

- exclude adjacent faces already handled by shield planes
- exclude the target face's own boundary edges/vertices from the foreign set

Exit criteria:

- for a target face, the helper returns the small finite set of truly relevant
  foreign features

### Phase 5: Foreign Occupied Volume Builder

Goal:

- build `V(i)` from actual occupied geometry, not raw adjacency

Deliverables:

- helper that instantiates and unions foreign:
  - face occupancy
  - edge occupancy
  - vertex occupancy
- initially, this can still be face-led if that is what is available, but the
  abstraction should be "occupied feature", not "raw face slab"

Important rule:

- this is foreign occupied volume only
- local clearance is not reused as the foreign cutter

Exit criteria:

- `V(i)` can be rendered independently for one target face and inspected

### Phase 6: Final Face Carve

Goal:

- build the finished printable face shape from the three agreed layers

Deliverables:

- `face_final(i) = face_shielded(i) - V(i) - C_local(i)`
- local edge/vertex clearance remains on the target face's own boundary only

Important rule:

- foreign cutters subtract from the shielded face
- local clearance is the final target-local subtraction
- do not reintroduce "clip local clearance by adjacent faces" behavior

Exit criteria:

- adjacent anti-interference is handled by shields alone
- truly penetrating non-adjacent features remove material
- local seating still works afterward

### Phase 7: Frame / Seat Consumer Integration

Goal:

- make the face result usable by the frame-building path

Deliverables:

- explicit distinction between:
  - `face_shielded(i)`
  - `face_final(i)` / seated face
- frame subtraction should use the smaller seated face where appropriate

Exit criteria:

- frame seats use the intended smaller face shape
- anti-interference and seating are no longer conflated

### Phase 8: Only Then Revisit Endpoint / Run-End Geometry

Goal:

- return to the harder analytic cut-span endpoint problem once the proxy face
  model is settled

Important rule:

- do not mix this back into Phases 1-7
- run-end planes and repeated-cut-span relief are a separate analytic problem

Why deferred:

- the current proxy effort is about the correct occupied/seated face shape
- run-end plane work belongs to admissible-region exactness in
  `face_regions.scad`

## Guardrails For Execution

At each phase:

- land one helper family at a time
- keep a visible example/probe for the current phase
- add focused tests for new helpers, not just render smoke where a geometric
  assertion is possible
- do not mix:
  - arrangement extraction
  - foreign hit-list generation
  - local clearance tuning
  - endpoint relief
  in one patch

The immediate start point is Phase 1, not more clearance/cutter tweaking.

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

2a. Keep split cut spans distinct on each visible cell boundary.

- `cut_entry_id` and `cut_pair_id` are not enough for repeated same-cutter spans
- add `cut_run_id` local to one visible cell boundary
- use it to distinguish split cut spans when the same cutter re-enters the
  same visible cell in multiple disjoint runs

Why:

- endpoint relief needs to know which local span endpoint it is acting on
- repeated spans from the same cutter must not be conflated

3. Represent every cell boundary edge as an exact 3D plane.

For each edge of one visible cell:

- parent edge:
  - plane from face dihedral/2
- cut edge:
  - plane from cutter face / join bisector
  - plus `seg_cut_clearance`

This should be a real 3D plane equation, not a sampled `u(z)` profile.

3a. Add a distinct cut-span endpoint / cut-vertex model.

This is now a separate requirement from the edge planes themselves.

For a visible cell, classify each vertex as:

- parent-parent
- parent-cut
- cut-cut
- repeated-cutter re-entry / punch-through neighbourhood

Then apply local endpoint relief rules:

- parent-parent:
  - no extra treatment
- parent-cut:
  - one-sided retreat aligned to the cut span endpoint
- cut-cut:
  - vertex bisector / chamfer rule derived from the two incident cut
    constraints
- repeated-cutter re-entry / punch-through:
  - treat the two cut-span endpoints independently; do not assume the whole
    same-face contribution is one continuous join

The important point is that edge planes and endpoint relief are different
layers:

- the edge planes define the correct side faces of the admissible region
- the run-end bounds stop a finite cut span from behaving like an infinite
  cutter join

Crucially:

- `z0` and `z1` only bound the admissible region
- they do **not** imply different geometric behavior at the "top" and
  "bottom"
- the principled target is therefore:
  - full-depth side planes for parent and cut edges
  - full-depth run-end planes for finite cut spans

So any future endpoint handling should be expressed as extra bounding planes,
not as top-only notches, tapers, or cap-specific heuristics.

The same principle should eventually apply to the live cut-edge path:

- active cut-edge side planes ought to be derived across the full active
  `z0..z1` span
- a privileged interior `band_z0` / `band_z1` is the wrong abstraction for
  core geometry

But this is still design intent, not a landed live-path guarantee: a direct
swap to full-span construction caused broader region-builder regressions and
was backed out.

Implementation guardrails:

- do **not** derive endpoint relief independently per convex atom
- derive it on the original visible cell first, then propagate it into the
  convex atoms used by the live region builder
- do **not** assume a triangulated convex atom is automatically a valid
  half-plane-bounded region for clipped-loop lofting; atoms containing
  triangulation `"inner"` edges are under-constrained and can blow up far
  beyond the true atom footprint

Why:

- one original visible-cell vertex can appear in multiple atoms after
  triangulation
- if each atom invents its own local endpoint chamfer, the same real endpoint
  can be cut twice or inconsistently

Also:

- the first endpoint-relief implementation must preserve the sampled loop arity
  used by the current live convex-atom loft builder
- otherwise we will just reintroduce the earlier "exact loop is right, loft is
  broken" failure mode under a different name

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

## Alternative Direction: Proxy Interaction Contract

The current analytic path is trying to infer all admissible face segmentation
from planes, loops, and lofted region approximations. That remains useful for:

- topology / visibility analysis
- cut provenance
- pairing / run metadata
- debug labeling

But there is a second plausible path for fabrication-oriented consumers:

- build **proxy occupancy solids**
- subtract those directly with ordinary OpenSCAD CSG

This should be treated as a separate interaction engine, not a replacement for
the topology layer.

### Occupancy vs Clearance

The proxy path must distinguish between two different kinds of solids:

- **occupancy proxies**
  - material that really exists in the finished part
  - face occupancy
  - edge/frame occupancy
  - vertex occupancy
- **clearance proxies**
  - empty space intentionally reserved for fit/tolerance
  - face inset seats
  - edge seating strips / skeleton gap
  - vertex relief

These solve different problems:

- **inter-poly subtraction**
  - subtract foreign occupancy
  - this is what should make one intersecting poly cut through another
- **intra-poly fitting**
  - subtract local clearance
  - this is what creates seats/gaps so parts fit around a local frame

So the intended model is:

- foreign cutters should be built from the other poly's occupied structure,
  not from local seat-cutting strips
- local inset/skeleton gaps should be created by dedicated clearance proxies,
  not baked into face occupancy by default

Do not mix occupancy and clearance by default. If an inter-poly fit gap is
needed, it should be represented explicitly as a separate clearance/dilation
step, not by reusing the local seating proxy.

### Core Idea

For one feature family, define canonical local occupancy proxy geometry:

- child `0`: face occupancy proxy
- child `1`: edge occupancy proxy
- child `2`: vertex occupancy proxy

These are **proposal solids**, not already-clipped final geometry.

For a target face:

1. instantiate the target face proxy in the target face frame
2. instantiate neighboring face / edge / vertex proxies that can reach that
   face
3. clip those neighbors to a local influence bound
4. subtract the union of those clipped neighbors from the target face proxy

That gives a fabrication interaction result directly from SCAD solids, rather
than inferring every detail analytically.

### Important Constraint

The proxy path must use **raw proposals only**.

Do **not** recursively build already-clipped neighbors and subtract those.
That creates a circular problem:

- face A depends on face B's final shape
- face B depends on face A's final shape

Instead:

- each feature contributes one local raw occupancy proxy
- interaction is resolved only at the target being built

### Why This Could Help

This path naturally handles:

- split cut spans
- punch-throughs
- future edge/vertex structure
- finite extents of real user geometry

because the proxy already encodes finite occupancy.

It also matches the actual fabrication question more directly:

- "what other geometry would physically intrude into this feature volume?"

instead of:

- "what plane/loop approximation should stand in for that intrusion?"

### Why This Does Not Replace the Analytic Path

The analytic path is still needed for:

- visible-cell enumeration
- world-stable cut pairing
- split cut span metadata (`cut_run_id`)
- feature classification
- debugging / labels / reasoning about where interference comes from

So the likely end state is:

- **analytic layer**: determine candidate interacting features and local bounds
- **proxy CSG layer**: perform the actual fabrication subtraction using proxy
  solids

### Proposed API Contract

Introduce a new family of occupancy-proxy interaction modules that consume up
to three children:

```scad
ps_proxy_interaction_face(... ) {
    // child 0: face occupancy proxy
    // child 1: edge occupancy proxy
    // child 2: vertex occupancy proxy
}
```

The current signature sketch for this work lives in:

- `src/polysymmetrica/core/proxy_interaction.scad`

Or, more explicitly:

```scad
ps_clip_face_by_feature_proxies(
    poly,
    face_idx,
    face_bounds = [z0, z1],
    face_proxy_mode = "raw" | "sweep_to_bounds",
    edge_radius = ...,
    edge_length = ...,
    vertex_radius = ...,
    include_faces = true,
    include_edges = true,
    include_vertices = true,
    face_indices = ...,
    edge_indices = ...,
    vertex_indices = ...,
    filter = ...
) {
    children();
}
```

Conceptually:

- child `0` is required if face-face occupancy interaction is desired
- child `1` is optional edge occupancy proxy
- child `2` is optional vertex occupancy proxy

The implementation would:

- instantiate the target face proxy in the target face frame
- instantiate interfering proxies on neighboring faces / edges / vertices
- clip each interfering proxy to its influence bound
- subtract the union of those clipped interfering proxies from the target face
  proxy

For face-face interaction, a useful first realization mode is:

- `face_proxy_mode = "sweep_to_bounds"`
  - project the neighboring face proxy to its local XY footprint
  - extrude that footprint through `face_bounds`
  - use the resulting bounded occupancy solid as the subtractor

Helpful supporting primitive:

- `place_on_faces/edges/vertices(..., indices=[...])`
  - lets the proxy path reuse the existing exact placement frames while
    instantiating only the specific neighboring features that matter

Exact site lists are first-class inputs to the proxy path:

- `face_indices`
- `edge_indices`
- `vertex_indices`

`undef` means “all candidate features of that family”, while a list means
“instantiate exactly these known influencing sites”.

A future clearance path may mirror that structure, but it should stay a
separate API or an explicitly separate mode:

- face clearance proxy
- edge clearance proxy
- vertex clearance proxy

### Ownership Zones

Ownership zones remain a plausible future direction, but they are currently
parked.

Reason:

- the first zone/corridor experiments repeatedly over-constrained the raw edge
  strip before subtraction
- that reshaping was exactly what broke the previously good edge-separation
  behavior

So for now:

- keep ownership/corridor work out of the active proxy baseline
- keep the idea documented as a future refinement
- do not let corridor logic reshape the raw edge proxy until there is a clear
  consumer and a proof that it improves rather than distorts the result

### Bounds Contract

Because proxy geometry is user-defined, the core API must ask for explicit
influence bounds:

- `face_bounds = [z0, z1]`
  - local face-normal extent of the face proxy
- `edge_radius`
  - how far an edge proxy may influence surrounding face/edge space
- `edge_length`
  - how far along the edge tangent the proxy should be considered
- `vertex_radius`
  - spherical / local bound for vertex proxies

These are not stylistic parameters; they are interaction bounds so that the
CSG problem stays local and finite.

### Scope Rules

The first proxy-interaction implementation should be conservative:

- build one target at a time
- use only raw neighboring proposals
- do not recurse
- do not attempt fixed-point convergence
- do not use final decorative geometry as the interaction proxy unless the user
  explicitly wants that

The recommended default is:

- use a simplified interaction proxy
- then separately render / clip the final decorative geometry through that
  result

### Practical Implication

If this path works, a large part of the current segmentation-specific math
becomes:

- analytic help for deciding **what interacts**

rather than:

- the sole mechanism for constructing the final removed volume

That is promising for future face/edge/vertex fabrication work, especially
where explicit SCAD solids are easier to reason about than ever-more-complicated
plane intersection logic.

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

7. There is a real distinction between edge-angle correctness and endpoint
   correctness.

The `f2` / `f11` case established that:

- the cut side faces can be parallel and correctly spaced at the lower `z`
  portion of the span
- the local cut dihedral can be right
- but the top surface can still overrun at the ends of that span

So once the side planes are trusted, the next bug class is not "wrong dihedral"
but "missing cut-span endpoint relief".

8. Repeated spans from the same cutter face are a real case.

The `f2/c0` chain

- `e2/c4/f11`
- `e3/c3/f7`
- `e4/c5/f13`
- `e5/c4/f11`

shows that one cutter face can re-enter the same visible cell in multiple
disjoint spans, with other cuts intervening.

That means:

- same cutter face does not imply one simple continuous join on the live piece
- `cut_entry_id` / `cut_pair_id` are not sufficient by themselves to identify a
  specific cut-span endpoint inside one visible cell
- endpoint handling must therefore be keyed per **cut run / span**, not just per
  cutter face or pair id globally

So the next metadata addition is likely to be:

- `cut_run_id` (or equivalent), local to one visible cell boundary, grouping
  one continuous cut span between its two endpoints

This should be derived after visible-cell construction, not from the raw global
cut-entry list.
- endpoint handling must work per span, not per cutter face globally

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
- add endpoint-focused debug views only if they show the local top-cap geometry
  clearly; avoid broad translucent CSG overlays
- when debugging endpoint relief, label cut spans/runs explicitly, not just
  cutter-face ids, because the same cutter can re-enter the same cell multiple
  times

4. Tests

- sample-point tests on exact region membership
- pair consistency tests across matched seg edges
- endpoint-specific tests for repeated same-cutter spans and cut-cut vertices
- canonical visual/probe cases:
  - `poly_prism(5,2)`
  - `poly_prism(7,2)`
  - `poly_prism(7,3)`
  - `poly_antiprism(7,3, angle=15)`

Additional canonical endpoint-relief cases:

- `poly_antiprism(7,3, angle=15)`:
  - `f2/c0` with the repeated `f11` spans
  - verify that the lower cut faces remain parallel and spaced
  - verify that the top cap does not overrun at either end of the `f11` spans
  - verify that the two `f11` spans are treated as distinct runs with distinct
    endpoints, not merged because they share a cutter face
- any future edge-mounted example where a true edge punches through a face:
  - use it to validate that endpoint relief composes with `place_on_edges(...)`
  - endpoint relief must not consume the true-edge neighbourhood so aggressively
    that edge-mounted geometry loses a predictable interface

## Recommendation

Do not touch `face_plate` next, except as a thin consumer or debug surface.

The next real work unit should be:

- keep `segments.scad` stable
- keep the current live geometry path stable
- keep the convex-atom decomposition
- add cut-run / span metadata on visible-cell boundaries
- derive endpoint relief on the original visible cell, then map it into convex
  atoms
- add cut-span endpoint / cut-vertex relief as a new layer in
  `face_regions.scad`, without revisiting the already-correct side-plane math
- preserve per-atom loop arity in the first implementation so the current loft
  builder remains valid
- only pursue future exact-region work per convex atom
- accept any future replacement only when its horizontal slices match the
  accepted local cross-sections across the actual occupied `z` span

That is the route that still looks robust enough for stellations instead of
just patching the current antiprism bug.
