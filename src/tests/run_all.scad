use <core/TestFuncs.scad>
use <core/TestDuals.scad>
use <core/TestCantellation.scad>
use <core/TestTruncation.scad>
use <core/TestValidity.scad>

echo("=== PolySymmetrica tests: START ===");

run_TestFuncs();
run_TestDuals();
run_TestCantellation();
run_TestTruncation();
run_TestValidity();

echo("=== PolySymmetrica tests: PASS ===");

color("green") 
linear_extrude(height = 2)
    text(str("Tests Passed!"));
