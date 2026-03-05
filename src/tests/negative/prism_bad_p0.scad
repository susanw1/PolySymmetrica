use <../../polysymmetrica/core/prisms.scad>

// EXPECT FAIL: p must be >= 1
_ = poly_prism(5, p=0);

