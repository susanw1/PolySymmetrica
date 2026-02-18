use <../../core/params.scad>

/*
Params-overrides basics example.

Shows that `face_fid` is a family-id array indexed by face index.
Here we use a synthetic 4-face case:

    face index:  0  1  2  3
    face_fid:   [0, 1, 1, 0]

Meaning:
- faces 0 and 3 belong to face family 0
- faces 1 and 2 belong to face family 1
*/

rows = [
    ["face", "all", ["angle", 15]],
    ["face", "family", 1, ["df", 0.04]],
    ["face", "id", [2, 3], ["angle", 19]],
    ["vert", "family", 0, ["c", 0.06]]
];

face_fid = [0, 1, 1, 0];

face_df = ps_params_compile_key(rows, "face", "df", 4, face_fid);
face_angle = ps_params_compile_key(rows, "face", "angle", 4, face_fid);

echo("params example: face_fid=", face_fid);
echo("params example: face_df=", face_df);         // [undef, 0.04, 0.04, undef]
echo("params example: face_angle=", face_angle);   // [15, 15, 19, 19]
ps_params_print(rows);

// Minimal geometry so this file renders to a valid STL target.
cube(1, center=true);

