use <../../polysymmetrica/core/prisms.scad>

// EXPECT FAIL: p must satisfy p < n/2
_ = poly_antiprism(6, p=3);

