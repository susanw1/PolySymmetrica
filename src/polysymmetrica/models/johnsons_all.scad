// ---------------------------------------------------------------------------
// PolySymmetrica - Johnson solids (more to come!)
// Version: 0.1.0
// Copyright 2026 Susan Witts
// SPDX-License-Identifier: MIT

use <../core/funcs.scad>

// *** WARNING! *** 
// These definitions are currently for experimentation only, to verify correct behaviour on non-uniform
// shapes! One day we will have a complete list of all the johnsons, created with construction.
// Even the names here are just placeholders.

// J1: square pyramid (unit edge length, LHR orientation).
function j1_square_pyramid() =
    let(h = sqrt(2) / 2)
    make_poly(
        [
            [ 0.5,  0.5, 0],  // 0
            [-0.5,  0.5, 0],  // 1
            [-0.5, -0.5, 0],  // 2
            [ 0.5, -0.5, 0],  // 3
            [ 0,    0,   h]   // 4 apex
        ],
        [
            [0, 1, 2, 3],  // base (outside faces -Z)
            [0, 4, 1],
            [1, 4, 2],
            [2, 4, 3],
            [3, 4, 0]
        ]
    );

// J3: triangular cupola (approximate coords; preview-only).
// Not unit-edge; intended for placement/orientation experiments.
function j3_triangular_cupola_approx() =
    let(
        r = 1,
        s = sqrt(3) / 2,
        h = 0.6,
        verts = [
            [ r, 0, 0],           // 0
            [ r/2, s, 0],         // 1
            [-r/2, s, 0],         // 2
            [-r, 0, 0],           // 3
            [-r/2, -s, 0],        // 4
            [ r/2, -s, 0],        // 5
            [ r, 0, h],           // 6 (top)
            [-r/2, s, h],         // 7 (top)
            [-r/2, -s, h]         // 8 (top)
        ],
        faces_raw = [
            [6, 8, 7],            // top triangle (LHR)
            [0, 1, 2, 3, 4, 5],   // bottom hex (LHR from below)
            [6, 1, 0, 5],         // square
            [7, 3, 2, 1],         // square
            [8, 5, 4, 3],         // square
            [6, 7, 1],            // triangle
            [7, 8, 3],            // triangle
            [8, 6, 5]             // triangle
        ],
        faces = ps_orient_all_faces_outward(verts, faces_raw)
    )
    make_poly(verts, faces);

// J40: elongated pentagonal orthocupolarotunda (imported coords; preview-only).
// Faces are reoriented outward for LHR compatibility.
function j40_elongated_pentagonal_orthocupolarotunda_approx() =
    let(
        verts = [
            [-1.05518, -0.061289, -0.047893],
            [-0.934164, 0.280612, 0.409939],
            [-0.859454, -0.241561, -0.56784],
            [-0.777777, -0.505581, 0.210572],
            [-0.776073, 0.311702, -0.400212],
            [-0.656761, -0.163679, 0.668404],
            [-0.655057, 0.653604, 0.05762],
            [-0.582051, -0.685853, -0.309375],
            [-0.542629, 0.653549, 0.63078],
            [-0.421745, -0.191346, -0.9513],
            [-0.400139, -0.685942, 0.618017],
            [-0.338365, 0.361918, -0.783672],
            [-0.265226, 0.209257, 0.889245],
            [-0.144342, -0.635637, -0.692835],
            [-0.142556, 0.915126, -0.042884],
            [-0.083446, -0.977628, -0.223275],
            [-0.030129, 0.915071, 0.530275],
            [0.028982, -0.977682, 0.349884],
            [0.05317, 0.734854, -0.562831],
            [0.090755, 0.070177, -1.051805],
            [0.149998, -0.635781, 0.807716],
            [0.233378, -0.082518, 0.975345],
            [0.247274, 0.47078, 0.788741],
            [0.368158, -0.374115, -0.793339],
            [0.40758, 0.965287, 0.146815],
            [0.46669, -0.927467, -0.033576],
            [0.48229, 0.443114, -0.830964],
            [0.603306, 0.785015, -0.373132],
            [0.662498, -0.374258, 0.707212],
            [0.684983, 0.520995, 0.405281],
            [0.745797, -0.554475, -0.385895],
            [0.759693, -0.001178, -0.572498],
            [0.858225, -0.55453, 0.187265],
            [0.880709, 0.340724, -0.114666],
            [0.941605, -0.001267, 0.354893]
        ],
        faces_raw = [
            [6, 1, 8],
            [14, 16, 24],
            [18, 27, 26],
            [11, 19, 9],
            [4, 2, 0],
            [20, 10, 17],
            [10, 5, 3],
            [17, 15, 25],
            [15, 7, 13],
            [25, 30, 32],
            [30, 23, 31],
            [32, 34, 28],
            [34, 33, 29],
            [28, 21, 20],
            [21, 22, 12],
            [6, 4, 0, 1],
            [14, 6, 8, 16],
            [18, 14, 24, 27],
            [11, 18, 26, 19],
            [4, 11, 9, 2],
            [26, 31, 23, 19],
            [19, 23, 13, 9],
            [9, 13, 7, 2],
            [2, 7, 3, 0],
            [0, 3, 5, 1],
            [1, 5, 12, 8],
            [8, 12, 22, 16],
            [16, 22, 29, 24],
            [24, 29, 33, 27],
            [27, 33, 31, 26],
            [4, 6, 14, 18, 11],
            [20, 17, 25, 32, 28],
            [20, 21, 12, 5, 10],
            [17, 10, 3, 7, 15],
            [25, 15, 13, 23, 30],
            [32, 30, 31, 33, 34],
            [28, 34, 29, 22, 21]
        ],
        faces = ps_orient_all_faces_outward(verts, faces_raw)
    )
    make_poly(verts, faces);
