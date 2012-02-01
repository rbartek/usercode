// @(#)root/tmva $Id: ClassBOOST_ZH.C,v 1.2 2011/12/14 14:50:10 madfish Exp $
/**********************************************************************************
 * Project   : TMVA - a ROOT-integrated toolkit for multivariate data analysis    *
 * Package   : TMVA                                                               *
 * Root Macro: TMVAClassification                                                 *
 *                                                                                *
 * This macro provides examples for the training and testing of the               *
 * TMVA classifiers.                                                              *
 *                                                                                *
 * As input data is used a toy-MC sample consisting of four Gaussian-distributed  *
 * and linearly correlated input variables.                                       *
 *                                                                                *
 * The methods to be used can be switched on and off by means of booleans, or     *
 * via the prompt command, for example:                                           *
 *                                                                                *
 *    root -l ./TMVAClassification.C\(\"Fisher,Likelihood\"\)                     *
 *                                                                                *
 * (note that the backslashes are mandatory)                                      *
 * If no method given, a default set of classifiers is used.                      *
 *                                                                                *
 * The output file "TMVA.root" can be analysed with the use of dedicated          *
 * macros (simply say: root -l <macro.C>), which can be conveniently              *
 * invoked through a GUI that will appear at the end of the run of this macro.    *
 * Launch the GUI via the command:                                                *
 *                                                                                *
 *    root -l ./TMVAGui.C                                                         *
 *                                                                                *
 **********************************************************************************/

#include <cstdlib>
#include <iostream>
#include <map>
#include <string>
#include "xsecnov10f.h"
#include "TChain.h"
#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TObjString.h"
#include "TSystem.h"
#include "TROOT.h"

#include "TMVAGui.C"

#if not defined(__CINT__) || defined(__MAKECINT__)
// needs to be included when makecint runs (ACLIC)
//gSystem->Load( "TMVA/lib/libTMVA.1" )
#include "TMVA/Factory.h"
#include "TMVA/Tools.h"
#endif

double max(double csv1, double  csv2)
{
if(csv1 > csv2) return csv1;
else return csv2;
}

double min(double csv1, double  csv2)
{
if(csv1 > csv2) return csv2;
else return csv1;
}

double deltaPhi(double phi1, double phi2) 
{ 
    double PI = 3.14159265;
    double result = phi1 - phi2;
    while (result > PI) result -= 2*PI;
    while (result <= -PI) result += 2*PI;
    return result;
}


void ClassBOOST_ZH( TString myMethodList = "" )
{

double lumi = 100.0;
int lnum = 115; 
std::cout << "Training " << lnum <<std::endl;


   //ONLY OF THESE MAY BE TRUE!

   bool trainALL   = 1;

   // The explicit loading of the shared libTMVA is done in TMVAlogon.C, defined in .rootrc
   // if you use your private .rootrc, or run from a different directory, please copy the
   // corresponding lines from .rootrc

   // methods to be processed can be given as an argument; use format:
   //
   // mylinux~> root -l TMVAClassification.C\(\"myMethod1,myMethod2,myMethod3\"\)
   //
   // if you like to use a method via the plugin mechanism, we recommend using
   //
   // mylinux~> root -l TMVAClassification.C\(\"P_myMethod\"\)
   // (an example is given for using the BDT as plugin (see below),
   // but of course the real application is when you write your own
   // method based)

   //---------------------------------------------------------------
   // This loads the library
   TMVA::Tools::Instance();

   // Default MVA methods to be trained + tested
   std::map<std::string,int> Use;

   // --- Cut optimisation
   Use["Cuts"]            = 0;
   Use["CutsD"]           = 0;
   Use["CutsPCA"]         = 0;
   Use["CutsGA"]          = 0;
   Use["CutsSA"]          = 0;
   // 
   // --- 1-dimensional likelihood ("naive Bayes estimator")
   Use["Likelihood"]      = 0;
   Use["LikelihoodD"]     = 0; // the "D" extension indicates decorrelated input variables (see option strings)
   Use["LikelihoodPCA"]   = 0; // the "PCA" extension indicates PCA-transformed input variables (see option strings)
   Use["LikelihoodKDE"]   = 0;
   Use["LikelihoodMIX"]   = 0;
   //
   // --- Mutidimensional likelihood and Nearest-Neighbour methods
   Use["PDERS"]           = 0;
   Use["PDERSD"]          = 0;
   Use["PDERSPCA"]        = 0;
   Use["PDEFoam"]         = 0;
   Use["PDEFoamBoost"]    = 0; // uses generalised MVA method boosting
   Use["KNN"]             = 0; // k-nearest neighbour method
   //
   // --- Linear Discriminant Analysis
   Use["LD"]              = 0; // Linear Discriminant identical to Fisher
   Use["Fisher"]          = 0;
   Use["FisherG"]         = 0;
   Use["BoostedFisher"]   = 0; // uses generalised MVA method boosting
   Use["HMatrix"]         = 0;
   //
   // --- Function Discriminant analysis
   Use["FDA_GA"]          = 0; // minimisation of user-defined function using Genetics Algorithm
   Use["FDA_SA"]          = 0;
   Use["FDA_MC"]          = 0;
   Use["FDA_MT"]          = 0;
   Use["FDA_GAMT"]        = 0;
   Use["FDA_MCMT"]        = 0;
   //
   // --- Neural Networks (all are feed-forward Multilayer Perceptrons)
   Use["MLP"]             = 0; // Recommended ANN
   Use["MLPBFGS"]         = 0; // Recommended ANN with optional training method
   Use["MLPBNN"]          = 0; // Recommended ANN with BFGS training method and bayesian regulator
   Use["CFMlpANN"]        = 0; // Depreciated ANN from ALEPH
   Use["TMlpANN"]         = 0; // ROOT's own ANN
   //
   // --- Support Vector Machine 
   Use["SVM"]             = 0;
   // 
   // --- Boosted Decision Trees
   Use["BDT"]             = 1; // uses Adaptive Boost
   Use["BDTG"]            = 0; // uses Gradient Boost
   Use["BDTB"]            = 0; // uses Bagging
   Use["BDTD"]            = 0; // decorrelation + Adaptive Boost
   // 
   // --- Friedman's RuleFit method, ie, an optimised series of cuts ("rules")
   Use["RuleFit"]         = 0;
   // ---------------------------------------------------------------

   std::cout << std::endl;
   std::cout << "==> Start TMVAClassification" << std::endl;

   // Select methods (don't look at this code - not of interest)
   if (myMethodList != "") {
      for (std::map<std::string,int>::iterator it = Use.begin(); it != Use.end(); it++) it->second = 0;

      std::vector<TString> mlist = TMVA::gTools().SplitString( myMethodList, ',' );
      for (UInt_t i=0; i<mlist.size(); i++) {
         std::string regMethod(mlist[i]);

         if (Use.find(regMethod) == Use.end()) {
            std::cout << "Method \"" << regMethod << "\" not known in TMVA under this name. Choose among the following:" << std::endl;
            for (std::map<std::string,int>::iterator it = Use.begin(); it != Use.end(); it++) std::cout << it->first << " ";
            std::cout << std::endl;
            return;
         }
         Use[regMethod] = 1;
      }
   }

   // --------------------------------------------------------------------------------------------------

   // --- Here the preparation phase begins

   // Create a ROOT output file where TMVA will store ntuples, histograms, etc.
   TString outfileName( "TMVA.root" );
   TFile* outputFile = TFile::Open( outfileName, "RECREATE" );

   // Create the factory object. Later you can choose the methods
   // whose performance you'd like to investigate. The factory is 
   // the only TMVA object you have to interact with
   //
   // The first argument is the base of the name of all the
   // weightfiles in the directory weight/
   //
   // The second argument is the output file for the training results
   // All TMVA output can be suppressed by removing the "!" (not) in
   // front of the "Silent" argument in the option string


   TMVA::Factory *factory = new TMVA::Factory( "TMVAClassification", outputFile,
                                               "!V:!Silent:Color:DrawProgressBar:Transformations=I;D;P;G,D:AnalysisType=Classification" );


   // If you wish to modify default settings
   // (please check "src/Config.h" to see all available global options)
   //    (TMVA::gConfig().GetVariablePlotting()).fTimesRMS = 8.0;
   //    (TMVA::gConfig().GetIONames()).fWeightFileDir = "myWeightDirectory";

   // Define the input variables that shall be used for the MVA training
   // note that you may also use variable expressions, such as: "3*var1/var2*abs(var3)"
   // [all types of expressions that can also be parsed by TTree::Draw( "expression" )]
//   factory->AddVariable( "myvar1 := var1+var2", 'F' );
//   factory->AddVariable( "myvar2 := var1-var2", "Expression 2", "", 'F' );
//   factory->AddVariable( "var3",                "Variable 3", "units", 'F' );


//bool lng = 0;

//if(lng)   factory->AddVariable( "cosRest := abs(Jet1.cosTheta) + abs(Jet2.cosTheta)", 'F' );
//   factory->AddVariable( "H.mass"    ,  "H Mass"    , "GeV"	, 'F' );
//   factory->AddVariable( "H.pt"      ,  "H PT"      , "GeV" 	, 'F' );
//   factory->AddVariable( "mincsv := min(hJet_csv[0],hJet_csv[1])"      ,  "CSV 1"     , ""		, 'F' );
//   factory->AddVariable( "maxcsv := max(hJet_csv[0],hJet_csv[1])"      ,  "CSV 2"     , "" 	, 'F' );

   factory->AddVariable( "Hmass"    ,  "H Mass"    , "GeV"	, 'F' );
   factory->AddVariable( "Hpt"      ,  "H PT"      , "GeV" 	, 'F' );
//   factory->AddVariable( "mincsv := min(hJet_csvCor[0],hJet_csvCor[1])"      ,  "CSV 1"     , ""		, 'F' );
//   factory->AddVariable( "maxcsv := max(hJet_csvCor[0],hJet_csvCor[1])"      ,  "CSV 2"     , "" 	, 'F' );
   factory->AddVariable( "CSV0"      ,  "CSV 1"     , ""		, 'F' );
   factory->AddVariable( "CSV1"      ,  "CSV 2"     , "" 	, 'F' );


   factory->AddVariable( "Zpt"        ,  "Z #rightarrow #mu#mu p_{T}"	     , "GeV"      , 'F' );
   factory->AddVariable( "hjjzmmdPhi := abs(deltaPhi(Hphi,Zphi))"   ,  "H-Z #Delta #varphi", "rad"      , 'F' );
   factory->AddVariable( "deltaEta := abs(jetEta1-jetEta0)"      ,  "J-J Delta Eta"     , ""		, 'F' );
   factory->AddVariable( "Naj"      ,  "Num Add Jets"     , ""		, 'F' );

//   factory->AddVariable( "Zmass"        ,  "Z #rightarrow #mu#mu Mass"	     , "GeV"      , 'F' );
 
     // You can add so-called "Spectator variables", which are not used in the MVA training,
   // but will appear in the final "TestTree" produced by TMVA. This TestTree will contain the
   // input variables, the response values of all trained MVAs, and the spectator variables
   //  factory->AddSpectator( "zeePt",  "Z pt", "GeV" , 'F' );
//   factory->AddSpectator( "spec1 := var1*2",  "Spectator 1", "units", 'F' );
 //  factory->AddSpectator( "spec2 := var1*3",  "Spectator 2", "units", 'F' );

   // Read training and test data
   // (it is also possible to use ASCII format as input -> see TMVA Users Guide)
   
//   if (gSystem->AccessPathName( fname ))  // file does not exist in local directory
//      gSystem->Exec("wget http://root.cern.ch/files/tmva_class_example.root");


 TFile *inputBtt  =   TFile::Open("/home/hep/wilken/Weights/CMSSW_4_2_8_patch3/src/UserCode/wilken/CSVShapeCorr/TTJets_60M.root");
 TFile *inputBzz  =   TFile::Open("/home/hep/wilken/Weights/CMSSW_4_2_8_patch3/src/UserCode/wilken/CSVShapeCorr/ZZ.root"            );
 TFile *inputts   =   TFile::Open("/home/hep/wilken/Weights/CMSSW_4_2_8_patch3/src/UserCode/wilken/CSVShapeCorr/T_schannel.root"         );
 TFile *inputTt   =   TFile::Open("/home/hep/wilken/Weights/CMSSW_4_2_8_patch3/src/UserCode/wilken/CSVShapeCorr/T_tchannel.root"         );
 TFile *inputtw   =   TFile::Open("/home/hep/wilken/Weights/CMSSW_4_2_8_patch3/src/UserCode/wilken/CSVShapeCorr/T_tWchannel_DR.root"     );
 TFile *inputtbs   =  TFile::Open("/home/hep/wilken/Weights/CMSSW_4_2_8_patch3/src/UserCode/wilken/CSVShapeCorr/Tbar_schannel.root"      );
 TFile *inputTtb   =  TFile::Open("/home/hep/wilken/Weights/CMSSW_4_2_8_patch3/src/UserCode/wilken/CSVShapeCorr/Tbar_tchannel.root"      );
 TFile *inputtbw   =  TFile::Open("/home/hep/wilken/Weights/CMSSW_4_2_8_patch3/src/UserCode/wilken/CSVShapeCorr/Tbar_tWchannel_DR.root"  );
 TFile *inputww   =   TFile::Open("/home/hep/wilken/Weights/CMSSW_4_2_8_patch3/src/UserCode/wilken/CSVShapeCorr/WW.root"		      );
 TFile *inputwz   =   TFile::Open("/home/hep/wilken/Weights/CMSSW_4_2_8_patch3/src/UserCode/wilken/CSVShapeCorr/WZ.root"                 );
 TFile *inputwj   =   TFile::Open("/home/hep/wilken/Weights/CMSSW_4_2_8_patch3/src/UserCode/wilken/CSVShapeCorr/WJetsToLNu_TuneZ2.root"        );
 TFile *inputBzj   =  TFile::Open("/home/hep/wilken/Weights/CMSSW_4_2_8_patch3/src/UserCode/wilken/CSVShapeCorr/DYJetsToLL_M50.root"   );
 TFile *inputBzj2   = TFile::Open("/home/hep/wilken/Weights/CMSSW_4_2_8_patch3/src/UserCode/wilken/CSVShapeCorr/DYJetsToLL_PtZ.root");



/*
 TFile *inputBtt  = TFile::Open("A-TTJets_TuneZ2_7TeV-madgraph-tauola.root");                                   
 TFile *inputBzz  = TFile::Open("A-ZZ_TuneZ2_7TeV_pythia6_tauola.root");
 TFile *inputtbs   = TFile::Open("A-Tbar_TuneZ2_s-channel_7TeV-powheg-tauola.root");
 TFile *inputts   = TFile::Open("A-T_TuneZ2_s-channel_7TeV-powheg-tauola.root");
 TFile *inputTt   = TFile::Open("A-T_TuneZ2_t-channel_7TeV-powheg-tauola.root");
 TFile *inputTtb   = TFile::Open("A-Tbar_TuneZ2_t-channel_7TeV-powheg-tauola.root");
 TFile *inputtw   = TFile::Open("A-T_TuneZ2_tW-channel-DS_7TeV-powheg-tauola.root");
 TFile *inputtbw   = TFile::Open("A-Tbar_TuneZ2_tW-channel-DR_7TeV-powheg-tauola.root");
 TFile *inputww   = TFile::Open("A-WW_TuneZ2_7TeV_pythia6_tauola.root");
 TFile *inputwz   = TFile::Open("A-WZ_TuneZ2_7TeV_pythia6_tauola.root");
 TFile *inputwj   = TFile::Open("A-WJetsToLNu_TuneZ2_7TeV-madgraph-tauola.root");
 TFile *inputBzj   = TFile::Open("A-DYJetsToLL_TuneZ2_M-50_7TeV-madgraph-tauola.root");
 TFile *inputBzj2   = TFile::Open("A-DYJetsToLL_PtZ-100_TuneZ2_7TeV-madgraph-tauola.root");
*/
 TFile *inputS = TFile::Open("/home/hep/wilken/Weights/CMSSW_4_2_8_patch3/src/UserCode/wilken/CSVShapeCorr/ZH_M115.root");

//   std::cout << "--- TMVAClassification       : Using input file: " << input->GetName() << std::endl;
   
   // --- Register the training and test trees


 TTree *signal        = (TTree*) inputS->Get("TMVA_tree"); 
 TTree *backgroundtt  = (TTree*) inputBtt ->Get("TMVA_tree");   
 TTree *backgroundzz  = (TTree*) inputBzz ->Get("TMVA_tree"); 
 TTree *backgroundzj  = (TTree*) inputBzj ->Get("TMVA_tree"); 
 TTree *backgroundzj2  = (TTree*) inputBzj2 ->Get("TMVA_tree"); 
 TTree *backgroundTt = (TTree*) inputTt->Get("TMVA_tree");
 TTree *backgroundts  = (TTree*) inputts ->Get("TMVA_tree");
 TTree *backgroundtw  = (TTree*) inputtw ->Get("TMVA_tree");

 TTree *backgroundTtb = (TTree*) inputTtb->Get("TMVA_tree");
 TTree *backgroundtbs  = (TTree*) inputtbs ->Get("TMVA_tree");
 TTree *backgroundtbw  = (TTree*) inputtbw ->Get("TMVA_tree");
 

 TTree *backgroundww  = (TTree*) inputww ->Get("TMVA_tree");
 TTree *backgroundwz  = (TTree*) inputwz ->Get("TMVA_tree");
 TTree *backgroundwj  = (TTree*) inputwj ->Get("TMVA_tree");

 
// global event weights per tree (see below for setting event-wise weights)

  
   //Double_t backgroundWeight3 = lumi/lumiZJets;
   
   /*
   Double_t signalWeight1     = 1;
   Double_t backgroundWeight1 =1;
   Double_t backgroundWeight2 =1;
   Double_t backgroundWeight3 =1;
   */
   
   // You can add an arbitrary number of signal or background trees
  factory->AddSignalTree    ( signal	    , lumi/(lumiZH115/2.0)  );



if(trainALL)
{   
if(backgroundTtb->GetEntries() > 0)    factory->AddBackgroundTree(backgroundTtb,lumi/(lumiTtb/2.0));
if(backgroundtbs->GetEntries() > 0)    factory->AddBackgroundTree(backgroundtbs, lumi/(lumiTsb/2.0));
if(backgroundtbw->GetEntries() > 0)    factory->AddBackgroundTree(backgroundtbw, lumi/(lumiTtWb/2.0));
if(backgroundts->GetEntries() > 0) factory->AddBackgroundTree(backgroundts,lumi/(lumiTs/2.0));
if(backgroundTt->GetEntries() > 0) factory->AddBackgroundTree(backgroundTt,lumi/(lumiTt/2.0));
if(backgroundtw->GetEntries() > 0) factory->AddBackgroundTree(backgroundtw,lumi/(lumiTtW/2.0));
if(backgroundww->GetEntries() > 0) factory->AddBackgroundTree(backgroundww,lumi/(lumiWW/2.0));
if(backgroundwz->GetEntries() > 0) factory->AddBackgroundTree(backgroundwz,lumi/(lumiWZ/2.0));
if(backgroundwj->GetEntries() > 0) factory->AddBackgroundTree(backgroundwj,lumi/(lumiWJ/2.0));
if(backgroundzz->GetEntries() > 0) factory->AddBackgroundTree(backgroundzz   , lumi/(lumiZZ/2.0));
if(backgroundtt->GetEntries() > 0)   factory->AddBackgroundTree(backgroundtt   , lumi/(lumiTT/2.0));
if(backgroundzj->GetEntries() > 0) factory->AddBackgroundTree(backgroundzj,lumi/(lumiZJL/2.0));
if(backgroundzj2->GetEntries() > 0) factory->AddBackgroundTree(backgroundzj2,lumi/(lumiZJH/2.0));

}

//  factory->SetSignalWeightExpression("weight*");
//  factory->SetBackgroundWeightExpression( "weight*" );



   // To give different trees for training and testing, do as follows:
   //    factory->AddSignalTree( signalTrainingTree, signalTrainWeight, "Training" );
   //    factory->AddSignalTree( signalTestTree,     signalTestWeight,  "Test" );
   
   // Use the following code instead of the above two or four lines to add signal and background
   // training and test events "by hand"
   // NOTE that in this case one should not give expressions (such as "var1+var2") in the input
   //      variable definition, but simply compute the expression before adding the event
   //
   //     // --- begin ----------------------------------------------------------
   //     std::vector<Double_t> vars( 4 ); // vector has size of number of input variables
   //     Float_t  treevars[4], weight;
   //     
   //     // Signal
   //     for (UInt_t ivar=0; ivar<4; ivar++) signal->SetBranchAddress( Form( "var%i", ivar+1 ), &(treevars[ivar]) );
   //     for (UInt_t i=0; i<signal->GetEntries(); i++) {
   //        signal->GetEntry(i);
   //        for (UInt_t ivar=0; ivar<4; ivar++) vars[ivar] = treevars[ivar];
   //        // add training and test events; here: first half is training, second is testing
   //        // note that the weight can also be event-wise
   //        if (i < signal->GetEntries()/2.0) factory->AddSignalTrainingEvent( vars, signalWeight );
   //        else                              factory->AddSignalTestEvent    ( vars, signalWeight );
   //     }
   //   
   //     // Background (has event weights)
   //     background->SetBranchAddress( "weight", &weight );
   //     for (UInt_t ivar=0; ivar<4; ivar++) background->SetBranchAddress( Form( "var%i", ivar+1 ), &(treevars[ivar]) );
   //     for (UInt_t i=0; i<background->GetEntries(); i++) {
   //        background->GetEntry(i);
   //        for (UInt_t ivar=0; ivar<4; ivar++) vars[ivar] = treevars[ivar];
   //        // add training and test events; here: first half is training, second is testing
   //        // note that the weight can also be event-wise
   //        if (i < background->GetEntries()/2) factory->AddBackgroundTrainingEvent( vars, backgroundWeight*weight );
   //        else                                factory->AddBackgroundTestEvent    ( vars, backgroundWeight*weight );
   //     }
         // --- end ------------------------------------------------------------
   //
   // --- end of tree registration 

   // Set individual event weights (the variables must exist in the original TTree)
   //    for signal    : factory->SetSignalWeightExpression    ("weight1*weight2");
   //    for background: factory->SetBackgroundWeightExpression("weight1*weight2");
//   factory->SetBackgroundWeightExpression( "weight" );

/*
   factory->AddVariable( "H.mass"    ,  "H Mass"  , "GeV"	, 'F' );
   factory->AddVariable( "nHjj"       ,  "Num H"   , ""		, 'F' );
   //   factory->AddVariable( "H.pt"      ,  "H PT"    , "GeV" 	, 'F' );
   factory->AddVariable( "jjdr"       ,  "JJ DR"   , ""       	, 'F' );
   factory->AddVariable( "hJet_csv[0]"      ,  "CSV 1"   , ""		, 'F' );
   factory->AddVariable( "hJet_csv[1]"      ,  "CSV 2"   , "" 	, 'F' );
   factory->AddVariable( "hJet_pt[0]"       ,  "J1 PT"   , "GeV"	, 'F' );
   factory->AddVariable( "nofLeptons" ,  "Num Lep" , ""		, 'F' );
   factory->AddVariable( "zeeMass"    ,  "Z Mass"  , "GeV"      , 'F' );
   factory->AddVariable( "hjjzeedPhi" ,  "H-Z DPhi", "rad"      , 'F' );
   factory->AddVariable( "zeePt" ,  "Z PT", "GeV"      , 'F' );
*/

   // Apply additional cuts on the signal and background samples (can be different)
   //H.mass<126 && H.mass>102 &&
TCut mycuts,mycutb;
bool origCuts = false;

mycuts = "";//hJet_pt[0] > 60 && hJet_pt[1] > 20  && hJet_csv[0] > 0.68 && hJet_csv[1] >0.68 && BestCSVPair > 0 && V.mass < 60";
mycutb = "";//hJet_pt[0] > 60 && hJet_pt[1] > 20  && hJet_csv[0] > 0.68 && hJet_csv[1] >0.68 && BestCSVPair > 0 && V.mass < 60";
//mycuts = "hJet_pt[0] > 30 && hJet_pt[1] && V.pt > 150 && H.pt > 150 && ((hJet_csv[0] > 0.9 && hJet_csv[1] >0.5) || (hJet_csv[1] > 0.9 && hJet_csv[0] >0.5)) && hjjzdPhi > 2.70 && V.mass < 105 && V.mass > 75";
//mycutb = "hJet_pt[0] > 30 && hJet_pt[1] && V.pt > 150 && H.pt > 150 && ((hJet_csv[0] > 0.9 && hJet_csv[1] >0.5) || (hJet_csv[1] > 0.9 && hJet_csv[0] >0.5)) && hjjzdPhi > 2.70 && V.mass < 105 && V.mass > 75";


//mycuts = "H.pt > 150 && hJet_csv[0] > 0.5 && hJet_csv[1] > 0.5 && Z.pt > 150"; 
//mycutb = "H.pt > 150 && hJet_csv[0] > 0.5 && hJet_csv[1] > 0.5 && Z.pt > 150"; 



   // Tell the factory how to use the training and testing events
   //
   // If no numbers of events are given, half of the events in the tree are used 
   // for training, and the other half for testing:
   //    factory->PrepareTrainingAndTestTree( mycut, "SplitMode=random:!V" );
   // To also specify the number of testing events, use:
   // factory->PrepareTrainingAndTestTree( mycuts,"NSigTrain=300:NBkgTrain=300:NSigTest=300:NBkgTest=300:SplitMode=Random:!V" );


    //ZJETS
 //NumEvents 
  //ALL
 /*here*/if( trainALL)  factory->PrepareTrainingAndTestTree( mycuts, mycutb,"nTrain_Signal=0:nTrain_Background=0:SplitMode=Random:NormMode=None:!V" );

  
   // ---- Book MVA methods
   //
   // Please lookup the various method configuration options in the corresponding cxx files, eg:
   // src/MethoCuts.cxx, etc, or here: http://tmva.sourceforge.net/optionRef.html
   // it is possible to preset ranges in the option string in which the cut optimisation should be done:
   // "...:CutRangeMin[2]=-1:CutRangeMax[2]=1"...", where [2] is the third input variable

   // Cut optimisation
   if (Use["Cuts"])
      factory->BookMethod( TMVA::Types::kCuts, "Cuts",
                           "!H:!V:FitMethod=MC:EffSel:SampleSize=20.68:VarProp=FSmart" );

   if (Use["CutsD"])
      factory->BookMethod( TMVA::Types::kCuts, "CutsD",
                           "!H:!V:FitMethod=MC:EffSel:SampleSize=20.68:VarProp=FSmart:VarTransform=Decorrelate" );

   if (Use["CutsPCA"])
      factory->BookMethod( TMVA::Types::kCuts, "CutsPCA",
                           "!H:!V:FitMethod=MC:EffSel:SampleSize=20.68:VarProp=FSmart:VarTransform=PCA" );

   if (Use["CutsGA"])
      factory->BookMethod( TMVA::Types::kCuts, "CutsGA",
                           "H:!V:FitMethod=GA:CutRangeMin[0]=-10:CutRangeMax[0]=10:VarProp[1]=FMax:EffSel:Steps=30:Cycles=3:PopSize=400:SC_steps=10:SC_rate=5:SC_factor=0.95" );

   if (Use["CutsSA"])
      factory->BookMethod( TMVA::Types::kCuts, "CutsSA",
                           "!H:!V:FitMethod=SA:EffSel:MaxCalls=150000:KernelTemp=IncAdaptive:InitialTemp=1e+6:MinTemp=1e-6:Eps=1e-10:UseDefaultScale" );

   // Likelihood ("naive Bayes estimator")
   if (Use["Likelihood"])
      factory->BookMethod( TMVA::Types::kLikelihood, "Likelihood",
                           "H:!V:!TransformOutput:PDFInterpol=Spline2:NSmoothSig[0]=20:NSmoothBkg[0]=20:NSmoothBkg[1]=10:NSmooth=1:NAvEvtPerBin=50" );

   // Decorrelated likelihood
   if (Use["LikelihoodD"])
      factory->BookMethod( TMVA::Types::kLikelihood, "LikelihoodD",
                           "!H:!V:TransformOutput:PDFInterpol=Spline2:NSmoothSig[0]=20:NSmoothBkg[0]=20:NSmooth=5:NAvEvtPerBin=50:VarTransform=Decorrelate" );

   // PCA-transformed likelihood
   if (Use["LikelihoodPCA"])
      factory->BookMethod( TMVA::Types::kLikelihood, "LikelihoodPCA",
                           "!H:!V:!TransformOutput:PDFInterpol=Spline2:NSmoothSig[0]=20:NSmoothBkg[0]=20:NSmooth=5:NAvEvtPerBin=50:VarTransform=PCA" ); 

   // Use a kernel density estimator to approximate the PDFs
   if (Use["LikelihoodKDE"])
      factory->BookMethod( TMVA::Types::kLikelihood, "LikelihoodKDE",
                           "!H:!V:!TransformOutput:PDFInterpol=KDE:KDEtype=Gauss:KDEiter=Adaptive:KDEFineFactor=0.3:KDEborder=None:NAvEvtPerBin=50" ); 

   // Use a variable-dependent mix of splines and kernel density estimator
   if (Use["LikelihoodMIX"])
      factory->BookMethod( TMVA::Types::kLikelihood, "LikelihoodMIX",
                           "!H:!V:!TransformOutput:PDFInterpolSig[0]=KDE:PDFInterpolBkg[0]=KDE:PDFInterpolSig[1]=KDE:PDFInterpolBkg[1]=KDE:PDFInterpolSig[2]=Spline2:PDFInterpolBkg[2]=Spline2:PDFInterpolSig[3]=Spline2:PDFInterpolBkg[3]=Spline2:KDEtype=Gauss:KDEiter=Nonadaptive:KDEborder=None:NAvEvtPerBin=50" ); 

   // Test the multi-dimensional probability density estimator
   // here are the options strings for the MinMax and RMS methods, respectively:
   //      "!H:!V:VolumeRangeMode=MinMax:DeltaFrac=0.2:KernelEstimator=Gauss:GaussSigma=0.3" );
   //      "!H:!V:VolumeRangeMode=RMS:DeltaFrac=3:KernelEstimator=Gauss:GaussSigma=0.3" );
   if (Use["PDERS"])
      factory->BookMethod( TMVA::Types::kPDERS, "PDERS",
                           "!H:!V:NormTree=T:VolumeRangeMode=Adaptive:KernelEstimator=Gauss:GaussSigma=0.3:NEventsMin=400:NEventsMax=600" );

   if (Use["PDERSD"])
      factory->BookMethod( TMVA::Types::kPDERS, "PDERSD",
                           "!H:!V:VolumeRangeMode=Adaptive:KernelEstimator=Gauss:GaussSigma=0.3:NEventsMin=400:NEventsMax=600:VarTransform=Decorrelate" );

   if (Use["PDERSPCA"])
      factory->BookMethod( TMVA::Types::kPDERS, "PDERSPCA",
                           "!H:!V:VolumeRangeMode=Adaptive:KernelEstimator=Gauss:GaussSigma=0.3:NEventsMin=400:NEventsMax=600:VarTransform=PCA" );

   // Multi-dimensional likelihood estimator using self-adapting phase-space binning
   if (Use["PDEFoam"])
      factory->BookMethod( TMVA::Types::kPDEFoam, "PDEFoam",
                           "H:!V:SigBgSeparate=F:TailCut=0.001:VolFrac=0.0333:nActiveCells=500:nSampl=2000:nBin=5:Nmin=100:Kernel=None:Compress=T" );

   if (Use["PDEFoamBoost"])
      factory->BookMethod( TMVA::Types::kPDEFoam, "PDEFoamBoost",
                           "!H:!V:Boost_Num=30:Boost_Transform=linear:SigBgSeparate=F:MaxDepth=4:UseYesNoCell=T:DTLogic=MisClassificationError:FillFoamWithOrigWeights=F:TailCut=0:nActiveCells=500:nBin=20:Nmin=400:Kernel=None:Compress=T" );

   // K-Nearest Neighbour classifier (KNN)
   if (Use["KNN"])
      factory->BookMethod( TMVA::Types::kKNN, "KNN",
                           "H:nkNN=20:ScaleFrac=0.8:SigmaFact=1.0:Kernel=Gaus:UseKernel=F:UseWeight=T:!Trim" );

   // H-Matrix (chi2-squared) method
   if (Use["HMatrix"])
      factory->BookMethod( TMVA::Types::kHMatrix, "HMatrix", "!H:!V" );

   // Linear discriminant (same as Fisher discriminant)
   if (Use["LD"])
      factory->BookMethod( TMVA::Types::kLD, "LD", "H:!V:VarTransform=None:CreateMVAPdfs:PDFInterpolMVAPdf=Spline2:NbinsMVAPdf=50:NsmoothMVAPdf=10" );

   // Fisher discriminant (same as LD)
   if (Use["Fisher"])
      factory->BookMethod( TMVA::Types::kFisher, "Fisher", "H:!V:Fisher:CreateMVAPdfs:PDFInterpolMVAPdf=Spline2:NbinsMVAPdf=50:NsmoothMVAPdf=10" );

   // Fisher with Gauss-transformed input variables
   if (Use["FisherG"])
      factory->BookMethod( TMVA::Types::kFisher, "FisherG", "H:!V:VarTransform=Gauss" );

   // Composite classifier: ensemble (tree) of boosted Fisher classifiers
   if (Use["BoostedFisher"])
      factory->BookMethod( TMVA::Types::kFisher, "BoostedFisher", 
                           "H:!V:Boost_Num=20:Boost_Transform=log:Boost_Type=AdaBoost:Boost_AdaBoostBeta=0.2" );

   // Function discrimination analysis (FDA) -- test of various fitters - the recommended one is Minuit (or GA or SA)
   if (Use["FDA_MC"])
      factory->BookMethod( TMVA::Types::kFDA, "FDA_MC",
                           "H:!V:Formula=(0)+(1)*x0+(2)*x1+(3)*x2+(4)*x3:ParRanges=(-1,1);(-10,10);(-10,10);(-10,10);(-10,10):FitMethod=MC:SampleSize=10.68:Sigma=0.1" );

   if (Use["FDA_GA"]) // can also use Simulated Annealing (SA) algorithm (see Cuts_SA options])
      factory->BookMethod( TMVA::Types::kFDA, "FDA_GA",
                           "H:!V:Formula=(0)+(1)*x0+(2)*x1+(3)*x2+(4)*x3:ParRanges=(-1,1);(-10,10);(-10,10);(-10,10);(-10,10):FitMethod=GA:PopSize=300:Cycles=3:Steps=20:Trim=True:SaveBestGen=1" );

   if (Use["FDA_SA"]) // can also use Simulated Annealing (SA) algorithm (see Cuts_SA options])
      factory->BookMethod( TMVA::Types::kFDA, "FDA_SA",
                           "H:!V:Formula=(0)+(1)*x0+(2)*x1+(3)*x2+(4)*x3:ParRanges=(-1,1);(-10,10);(-10,10);(-10,10);(-10,10):FitMethod=SA:MaxCalls=15000:KernelTemp=IncAdaptive:InitialTemp=1e+6:MinTemp=1e-6:Eps=1e-10:UseDefaultScale" );

   if (Use["FDA_MT"])
      factory->BookMethod( TMVA::Types::kFDA, "FDA_MT",
                           "H:!V:Formula=(0)+(1)*x0+(2)*x1+(3)*x2+(4)*x3:ParRanges=(-1,1);(-10,10);(-10,10);(-10,10);(-10,10):FitMethod=MINUIT:ErrorLevel=1:PrintLevel=-1:FitStrategy=2:UseImprove:UseMinos:SetBatch" );

   if (Use["FDA_GAMT"])
      factory->BookMethod( TMVA::Types::kFDA, "FDA_GAMT",
                           "H:!V:Formula=(0)+(1)*x0+(2)*x1+(3)*x2+(4)*x3:ParRanges=(-1,1);(-10,10);(-10,10);(-10,10);(-10,10):FitMethod=GA:Converger=MINUIT:ErrorLevel=1:PrintLevel=-1:FitStrategy=0:!UseImprove:!UseMinos:SetBatch:Cycles=1:PopSize=5:Steps=5:Trim" );

   if (Use["FDA_MCMT"])
      factory->BookMethod( TMVA::Types::kFDA, "FDA_MCMT",
                           "H:!V:Formula=(0)+(1)*x0+(2)*x1+(3)*x2+(4)*x3:ParRanges=(-1,1);(-10,10);(-10,10);(-10,10);(-10,10):FitMethod=MC:Converger=MINUIT:ErrorLevel=1:PrintLevel=-1:FitStrategy=0:!UseImprove:!UseMinos:SetBatch:SampleSize=20" );

   // TMVA ANN: MLP (recommended ANN) -- all ANNs in TMVA are Multilayer Perceptrons
   if (Use["MLP"])
      factory->BookMethod( TMVA::Types::kMLP, "MLP", "H:!V:NeuronType=tanh:VarTransform=N:NCycles=600:HiddenLayers=N+5:TestRate=5:!UseRegulator" );

   if (Use["MLPBFGS"])
      factory->BookMethod( TMVA::Types::kMLP, "MLPBFGS", "H:!V:NeuronType=tanh:VarTransform=N:NCycles=600:HiddenLayers=N+5:TestRate=5:TrainingMethod=BFGS:!UseRegulator" );

   if (Use["MLPBNN"])
      factory->BookMethod( TMVA::Types::kMLP, "MLPBNN", "H:!V:NeuronType=tanh:VarTransform=N:NCycles=600:HiddenLayers=N+5:TestRate=5:TrainingMethod=BFGS:UseRegulator" ); // BFGS training with bayesian regulators

   // CF(Clermont-Ferrand)ANN
   if (Use["CFMlpANN"])
      factory->BookMethod( TMVA::Types::kCFMlpANN, "CFMlpANN", "!H:!V:NCycles=2000:HiddenLayers=N+1,N"  ); // n_cycles:#nodes:#nodes:...  

   // Tmlp(Root)ANN
   if (Use["TMlpANN"])
      factory->BookMethod( TMVA::Types::kTMlpANN, "TMlpANN", "!H:!V:NCycles=200:HiddenLayers=N+1,N:LearningMethod=BFGS:ValidationFraction=0.3"  ); // n_cycles:#nodes:#nodes:...

   // Support Vector Machine
   if (Use["SVM"])
      factory->BookMethod( TMVA::Types::kSVM, "SVM", "Gamma=0.25:Tol=0.001:VarTransform=Norm" );

   // Boosted Decision Trees
   if (Use["BDTG"]) // Gradient Boost
      factory->BookMethod( TMVA::Types::kBDT, "BDTG",
                           "!H:!V:NTrees=1000:BoostType=Grad:Shrinkage=0.10:UseBaggedGrad:GradBaggingFraction=0.5:nCuts=20:NNodesMax=5" );

   if (Use["BDT"])  // Adaptive Boost
      factory->BookMethod( TMVA::Types::kBDT, "BDT",

                       "!H:!V:NTrees=100:nEventsMin=100:MaxDepth=3:BoostType=AdaBoost:AdaBoostBeta=0.3:SeparationType=GiniIndex:nCuts=35:PruneMethod=NoPruning" );
   ///*MJJD*/                     "!H:!V:NTrees=100:nEventsMin=150:MaxDepth=2:BoostType=AdaBoost:AdaBoostBeta=0.3:SeparationType=GiniIndex:nCuts=35:PruneMethod=NoPruning" );
   if (Use["BDTB"]) // Bagging
      factory->BookMethod( TMVA::Types::kBDT, "BDTB",
                           "!H:!V:NTrees=400:BoostType=Bagging:SeparationType=GiniIndex:nCuts=20:PruneMethod=NoPruning" );

   if (Use["BDTD"]) // Decorrelation + Adaptive Boost
      factory->BookMethod( TMVA::Types::kBDT, "BDTD",
                           "!H:!V:NTrees=400:nEventsMin=400:MaxDepth=3:BoostType=AdaBoost:SeparationType=GiniIndex:nCuts=20:PruneMethod=NoPruning:VarTransform=Decorrelate" );

   // RuleFit -- TMVA implementation of Friedman's method
   if (Use["RuleFit"])
      factory->BookMethod( TMVA::Types::kRuleFit, "RuleFit",
                           "H:!V:RuleFitModule=RFTMVA:Model=ModRuleLinear:MinImp=0.001:RuleMinDist=0.001:NTrees=20:fEventsMin=0.01:fEventsMax=0.5:GDTau=-1.0:GDTauPrec=0.01:GDStep=0.01:GDNSteps=10000:GDErrScale=1.02" );

   // For an example of the category classifier usage, see: TMVAClassificationCategory

   // --------------------------------------------------------------------------------------------------

   // ---- Now you can optimize the setting (configuration) of the MVAs using the set of training events

   // factory->OptimizeAllMethods("SigEffAt001","Scan");
   // factory->OptimizeAllMethods("ROCIntegral","GA");

   // --------------------------------------------------------------------------------------------------

   // ---- Now you can tell the factory to train, test, and evaluate the MVAs

   // Train MVAs using the set of training events
   factory->TrainAllMethods();

   // ---- Evaluate all MVAs using the set of test events
   factory->TestAllMethods();

   // ----- Evaluate and compare performance of all configured MVAs
   factory->EvaluateAllMethods();

   // --------------------------------------------------------------

   // Save the output
   outputFile->Close();

   std::cout << "==> Wrote root file: " << outputFile->GetName() << std::endl;
   std::cout << "==> TMVAClassification is done!" << std::endl;

   delete factory;

   // Launch the GUI for the root macros
   if (!gROOT->IsBatch()) TMVAGui( outfileName );
   std::cout << "Trained "<< lnum <<std::endl;
}
