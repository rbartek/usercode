// Unnamed ROOT Macro which compiles our code and runs the optimizer
// To use this, run root like this:
// root -b -q runOptimizer.C &> logfile
//
{
	// Compile and keep classes
	gSystem->CompileMacro("TMVAClassification_BDT.class.C","k");
	gSystem->CompileMacro("SignificanceFinder.C","k");

	// Now we run the optimizer
	SignificanceFinder();
}



