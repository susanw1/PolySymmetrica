use <core/TestFuncs.scad>
use <core/TestDuals.scad>
use <core/TestTruncation.scad>

echo("=== PolySymmetrica tests: START ===");

run_TestFuncs();
run_TestDuals();
run_TestTruncation();

echo("=== PolySymmetrica tests: PASS ===");
