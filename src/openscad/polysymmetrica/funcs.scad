// Handy functions (and aliases)

function v_add(a, b)   = a + b;
function v_sub(a, b)   = a - b;
function v_scale(a, k) = a * k;             // scalar multiplication
function v_dot(a, b)   = a * b;             // dot product
function v_cross(a, b) = cross(a, b);       // OpenSCAD built-in
function v_len(a)      = norm(a);           // built-in length
function v_norm(a)     = let(L = norm(a)) (L == 0 ? [0,0,0] : a / L);

// Edge equality
function edge_equal(e1, e2) = (e1[0] == e2[0] && e1[1] == e2[1]);

/** Calculate polygon edge, given N and radius */ 
function calc_edge(n_vertex, rad) = 2 * rad * sin(180 / n_vertex);

/** Calculate polygon radius, given N and edge length */ 
function calc_radius(n_vertex, edge_len) = edge_len / (2 * sin(180 / n_vertex));

/** Sum of vector */
function sum(a, i = 0) = 
    i >= len(a) ? 0 : a[i] + sum(a, i + 1);

/** Sum of vector list (3D) */
function v_sum(list) = [
    sum([for (v = list) v[0]]),
    sum([for (v = list) v[1] ]),
    sum([for (v = list) v[2] ])
];
    
/** Face centroid from verts + index list */
function face_centroid(verts, f) =
    len(f) == 0
        ? [0,0,0]
        : v_scale(v_sum([for (vid = f) verts[vid]]), 1 / len(f));

// Face normal (not scaled, just direction)
function face_normal(verts, f) =
    v_norm(v_cross(
        verts[f[1]] - verts[f[0]],
        verts[f[2]] - verts[f[0]]
    ));
    