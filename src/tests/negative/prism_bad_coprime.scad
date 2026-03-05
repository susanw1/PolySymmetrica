use <../../polysymmetrica/core/prisms.scad>

// EXPECT FAIL: n and p must be coprime for single-cycle star polygons
_ = poly_prism(6, p=2);

