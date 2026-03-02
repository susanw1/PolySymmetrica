use <core/TestFuncs.scad>
use <core/TestDuals.scad>
use <core/TestCantellation.scad>
use <core/TestTruncation.scad>
use <core/TestCleanup.scad>
use <core/TestValidity.scad>
use <core/TestClassify.scad>
use <core/TestPlacement.scad>
use <core/TestPrisms.scad>

echo("=== PolySymmetrica tests: START ===");

run_TestFuncs();
run_TestDuals();
run_TestCantellation();
run_TestTruncation();
run_TestCleanup();
run_TestValidity();
run_TestClassify();
run_TestPlacement();
run_TestPrisms();

echo("=== PolySymmetrica tests: PASS ===");

color("green") 
linear_extrude(height = 2)
    text(str("Tests Passed!"));
