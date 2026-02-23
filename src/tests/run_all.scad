use <core/TestFuncs.scad>
use <core/TestDuals.scad>
use <core/TestCantellation.scad>
use <core/TestTruncation.scad>
use <core/TestValidity.scad>
use <core/TestClassify.scad>
use <core/TestPlacement.scad>

echo("=== PolySymmetrica tests: START ===");

run_TestFuncs();
run_TestDuals();
run_TestCantellation();
run_TestTruncation();
run_TestValidity();
run_TestClassify();
run_TestPlacement();

echo("=== PolySymmetrica tests: PASS ===");

color("green") 
linear_extrude(height = 2)
    text(str("Tests Passed!"));
