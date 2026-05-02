# Face-Interference Experiments

Developer-facing stress cases for self-crossing face analysis. These are not
printing examples; they keep arrangement, boundary, source-edge, dihedral, and
anti-interference data visible without reintroducing the proxy/carve stack.

- `test_7_3_15_punch_through_probe.scad`
  Focuses on the `poly_antiprism(7, 3, angle=15)` punch-through case, comparing
  a star face with a triangular face crossed by other arms. It also includes a
  simple proxy-replay panel using `place_on_face_foreign_proxy_sites(...)`.
- `test_antitruncated_tetrahedron_hex_probe.scad`
  Focuses on the self-crossing hex face produced by
  `poly_truncate(tetrahedron(), t=-0.5)`, where a single source edge can
  contribute boundary spans with opposite adjacent-face directions.
- `test_minimal_printable_punch_through_probe.scad`
  Minimal printable-style punch-through integration probe for
  `poly_antiprism(7,3, angle=15)`. It intersects raw face material with the
  visible-cell mask and positive anti-interference volume, then compares that
  keep-body with exact-intrusion strip clearance and an opt-in
  `face_plate_minus_foreign_proxies(...)` panel that subtracts replayed foreign
  face proxy bodies.
