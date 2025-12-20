use <../polysymmetrica/core/funcs.scad>
use <../polysymmetrica/core/placement.scad>

use <../polysymmetrica/models/regular_all.scad>


// Demo: draw neighbor-direction spikes around each vertex
module demo_vertex_neighbors(spike_len = 8, spike_d = 1, tip_d = 2.0) {
    // In vertex-local coords, vertex is at origin.
    // For each neighbor vector p=[x,y,z], draw a spike in that direction.
    for (p = $ps_vertex_neighbor_pts_local) {
        L = norm(p);
        if (L > 1e-9) {
            dir = p / L;

            // Shaft (a cylinder from origin)
            rotate(a = acos(dir[2]), v = v_cross([0,0,1], dir))
                cylinder(h = spike_len, d = spike_d, center = false, $fn=30);

            // Tip
            translate(dir * spike_len)
                sphere(d = tip_d, $fn=30);
        }
    }

    // Show valence as a little marker size (optional)
    sphere(d = spike_d + 0.3 * $ps_vertex_valence, $fn=30);
    text(str($ps_vertex_valence), size=5);
}

place_on_vertices_ir(hexahedron(), 30)
    demo_vertex_neighbors(spike_d = 2);
