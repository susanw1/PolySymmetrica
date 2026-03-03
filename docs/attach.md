# Face Attachment (`poly_attach`)

`poly_attach` attaches one closed polyhedron to another by aligning a selected face from each, removing those seam faces, and welding the seam.

Implemented in: `src/polysymmetrica/core/attach.scad`

## API

```scad
poly_attach(
    p1, p2,
    f1=0, f2=0,        // f1 accepts index or [idx...]
    rotate_step=0,
    scale_mode="fit_edge",   // "fit_edge" | "none"
    mirror=false,            // false=preserve chirality, true=mirror across seam
    eps=1e-8,
    cleanup=true,
    cleanup_eps=1e-8
)
```

## Preconditions

- `p1` and `p2` must pass `poly_valid(..., "closed")`.
- Selected faces must be planar (within `eps`).
- Selected faces must have equal arity (`len(face1) == len(face2)`).

## Parameter Meaning

- `f1`: seam face index (or list of indices) in `p1`.
  - scalar `f1=x` is shorthand for `f1=[x]`.
  - list mode attaches one transformed copy of `p2` per listed face.
- `f2`: seam face index in `p2` (same seam face used for each `f1` entry in list mode).
- `rotate_step`: cyclic correspondence shift for `p2` seam face vertices before alignment.
  - useful when same-arity faces have multiple valid edge matchings.
- `scale_mode`:
  - `"fit_edge"`: scales `p2` so mean edge length of face `f2` matches face `f1`.
  - `"none"`: no scaling; use only if seam faces already match size (otherwise seam welding may fail closed-valid checks).
- `mirror`:
  - `false` (default): orientation-preserving seam mapping (keeps `p2` chirality).
  - `true`: reflected seam mapping (legacy mirror behavior).
- `cleanup`: enables degenerate-face dropping in seam cleanup.
- `cleanup_eps`: tolerance for seam vertex merge and cleanup operations.

## How It Works

1. Outward-orients both input face sets.
2. Builds local frames for seam faces.
3. Maps all points of `p2` into `p1` seam frame, with opposite seam normal.
4. Drops all seam faces in `f1` from `p1` and seam face `f2` from each attached copy.
5. Concatenates meshes and seam-merges with `poly_cleanup(...)`.
6. Asserts final result is `poly_valid(..., "closed")`.

## Examples

Basic usage:

```scad
use <../src/polysymmetrica/core/attach.scad>
use <../src/polysymmetrica/models/platonics_all.scad>

q = poly_attach(hexahedron(), hexahedron(), f1=0, f2=0);
```

Edge correspondence shift:

```scad
q = poly_attach(hexahedron(), hexahedron(), f1=0, f2=0, rotate_step=1);
```

Attach one poly on multiple target faces:

```scad
q = poly_attach(hexahedron(), hexahedron(), f1=[0,1], f2=0);
```

Attach scaled poly with automatic seam scaling:

```scad
p2 = [
    [for (v = poly_verts(hexahedron())) v * 2],
    poly_faces(hexahedron()),
    poly_e_over_ir(hexahedron())
];
q = poly_attach(hexahedron(), p2, f1=0, f2=0, scale_mode="fit_edge");
```

See runnable demo: `src/polysymmetrica/examples/basics/main_attach.scad`

## Notes

- This is currently face-to-face only (no edge/vertex attach mode).
- It is strict by design: invalid seam topology fails early.
- For reproducible orientation/correspondence control, prefer explicit `f1/f2` and `rotate_step`.
- In list mode, `f1` entries must be unique and in-range.
