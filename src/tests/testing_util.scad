use <../polysymmetrica/core/funcs.scad>
use <../polysymmetrica/core/placement.scad>

RAD = 40;

module placement_tester(poly, rad=RAD) {
    color("yellow") 
    place_on_faces(poly, rad)
        cylinder(r = $ps_facet_radius, h = 2, $fn = $ps_vertex_count);
        
    color("red") 
    place_on_vertices(poly, rad)
        cylinder(h = 5, r = 5, $fn = $ps_vertex_valence);
        
    color("blue") 
    place_on_edges(poly, rad)
        rotate([0,90,0]) cylinder(r1 = 5, r2 = 1, h = $ps_edge_len/2);
}


// --- facet equivalence helpers ---

// rotate list left by k (k can be >= len)
function rotl(v, k) =
    let(n = len(v), s = (n==0 ? 0 : k % n))
    [ for (i = [0:n-1]) v[(i+s) % n] ];

// true if a and b are equal up to rotation (same direction), false otherwise
function facet_eq_rot(a, b) =
    (len(a) == len(b)) &&
    (len(a) == 0 ? true
                 : let(n = len(a))
                   // try all rotations of b
                   (max([ for (k = [0:n-1]) (a == rotl(b, k)) ? 1 : 0 ]) == 1));

// true if every facet in A can be matched with a distinct facet in B (multiset match)
function facets_eq_rot(A, B) =
    let(n = len(A))
    (n == len(B)) ?
        let(M = [
                for (i = [0:n-1])
                    [ for (j = [0:n-1])
                        facet_eq_rot(A[i], B[j]) ? 1 : 0
                    ]
            ]
        )
        _perfect_match(M)
    : false;


// --- matching implementation (Kuhn algorithm) ---
// Returns true if there is a perfect matching in bipartite graph M (n x n of 0/1)
function _perfect_match(M) =
    let(n = len(M))
    // matchR[j] = i matched to right node j, or -1
    let(matchR = [ for (j=[0:n-1]) -1 ])
    // try to match each left node in turn
    _pm_iter(M, 0, matchR);

// iterative over left side
function _pm_iter(M, i, matchR) =
    let(n = len(M))
    (i == n) ? true
             : let(res = _try_augment(M, i, matchR, [for (k=[0:n-1]) 0]))
               // res[0]=success?, res[1]=updated matchR
               (res[0] ? _pm_iter(M, i+1, res[1]) : false);

// DFS augmenting path from left node i
function _try_augment(M, i, matchR, seen) =
    let(n = len(M))
    _ta_j(M, i, 0, matchR, seen);

// walk j across right nodes
function _ta_j(M, i, j, matchR, seen) =
    let(n = len(M))
    (j == n) ? [false, matchR] :
    (M[i][j] == 0 || seen[j] == 1) ?
        _ta_j(M, i, j+1, matchR, seen) :
        // mark seen[j]=1
        let(seen2 = [ for (k=[0:n-1]) (k==j ? 1 : seen[k]) ])
        (
          (matchR[j] == -1) ?
              // free right node: match it
              [true, [ for (k=[0:n-1]) (k==j ? i : matchR[k]) ]] :
              // occupied: try to re-route the existing match
              let(next = _try_augment(M, matchR[j], matchR, seen2))
              (next[0] ?
                  [true, [ for (k=[0:n-1]) (k==j ? i : next[1][k]) ]] :
                  _ta_j(M, i, j+1, matchR, seen2)
              )
        );

// --- example ---
 A = [[1,2,3,4], [5,6,7]];
 B = [[3,4,1,2], [6,7,5]];
 B1 = [[6,7,5], [3,4,1,2]];
 C = [[3,4,1,2], [6,5,7]];
 echo(facets_eq_rot(A,B));  // true
 echo(facets_eq_rot(A,B1)); // true
 echo(facets_eq_rot(A,C));  // false


function assert_facet_matches(p1, p2) =
    assert(facets_eq_rot(poly_faces(p1), poly_faces(p2)));

function assert_verts_matches(A, B, eps=1e-9) =
    let(
        same_len = (len(A) == len(B)),
        diffs_ok =
            same_len &&
            min([
                for (i = [0 : len(A)-1]) 
                    norm(v_sub(A[i], B[i])) <= eps ? 1 : 0
            ]) == 1
    )
    assert(same_len, "Vertex arrays have different lengths")
    assert(diffs_ok, "Vertex coordinates differ beyond epsilon")
    true;
