// Handy functions (and aliases)

function v_add(a, b)   = a + b;
function v_sub(a, b)   = a - b;
function v_scale(a, k) = a * k;             // scalar multiplication
function v_dot(a, b)   = a * b;             // dot product
function v_cross(a, b) = cross(a, b);       // OpenSCAD built-in
function v_len(a)      = norm(a);           // built-in length
function v_norm(a)     = let(L = norm(a)) (L == 0 ? [0,0,0] : a / L);


/** Calculate polygon edge, given N and radius */ 
function calc_edge(n_vertex, rad) = 2 * rad * sin(180 / n_vertex);

/** Calculate polygon radius, given N and edge length */ 
function calc_radius(n_vertex, edge_len) = edge_len / (2 * sin(180 / n_vertex));

