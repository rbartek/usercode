// Unnamed ROOT Macro which compiles our code and runs the optimizer
// To use this, run root like this:
// root -b -q runBinSyst.C &> logfile
//
{
	// Compile and keep classes
	gSystem->CompileMacro("TMVAClassification_BDT.class.C","k");
	gSystem->CompileMacro("PropShapeUncert.C","k");

	// Now we run the optimizer
	PropShapeUncert();
}



