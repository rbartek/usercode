/**
    @brief Example of an analysis code to read tree ntuple
    and create histogram
    Usage:

    \verbatim
    mkdir [directoryName]
    Abb_commonTuples <inputfile> [outputfilename] [directoryName]
    \endverbatim

    @param inputfile Either a ROOT file or an ASCII file containing list of
    ROOT files.

    @param outputfile Name of the ROOT file which contains the histogram.
    Defaulted to 'output.root'

    @author Rachel Wilken <rachel.wilken@cern.ch>

    @date Fri Nov 9 2011

 */

#include "Zbb_commonTuples.h"

using namespace std;

int main(int argc, char** argv) {
	
    // There must be at least one argument, the input file.
    if (argc < 2) {
        std::cout << "Input filename is not specified! Exiting!" << std::endl;
        exit(1);
    }
	
    // First argument: a ROOT file or an ASCII file containing list of
    // ROOT files.
    std::string ifilename(argv[1]);
    std::string ofilename("output.root");
    std::string directory("plots");
	
    if (argc == 3) {
        // If there is a second argument, use it for the output filename.
        ofilename = std::string(argv[2]);
    }
	if (argc == 4) {
        // If there is a third argument, use it for the plots drectory
		ofilename = std::string(argv[2]);
		directory = std::string(argv[3]);
	}
	
    // Create the treeReader object.
    treeReader sample(ifilename, std::string("tree"));
	
    // If the input file(s) doesn't contain any event, exit.
    if (!(sample.readEvent(0))) {
        return 0;
    }
	
    // Set the sumw2 option for all histogram
    TH1::SetDefaultSumw2(true);
	
    // Create the output ROOT file
    TFile ofile(ofilename.c_str(),"RECREATE");
    ofile.cd();
	weight = SetWeight(ifilename);

    // Check how many events and file(s) are analyzed
    long int nEvent     = sample.nEvent();
    size_t   nFile      = sample.nFile();
	
    std::cout << "There are " << nEvent << " events in the chain."
	<< std::endl ;
    std::cout << "There are " << nFile  << " files in the chain."
	<< std::endl ;
    std::cout << "The tree name  is: " << sample.treeName()  << std::endl;
    std::cout << "The tree title is: " << sample.treeTitle() << std::endl;
	
	TTree *TMVA_tree = new TTree("TMVA_tree","Tree for TMVA input");
	TMVA_tree->Branch("nJets",&nJets, "nJets/I");
	TMVA_tree->Branch("CSV0",&CSV0, "CSV0/F");
	TMVA_tree->Branch("CSV1",&CSV1, "CSV1/F");
	TMVA_tree->Branch("Zmass",&Zmass, "Zmass/F");
	TMVA_tree->Branch("Hmass",&Hmass, "Hmass/F");
	TMVA_tree->Branch("DeltaPhiHV",&DeltaPhiHV, "DeltaPhiHV/F");
	TMVA_tree->Branch("Hpt",&Hpt, "Hpt/F");
	TMVA_tree->Branch("Zpt",&Zpt, "Zpt/F");
	TMVA_tree->Branch("mu0pt",&mu0pt, "mu0pt/F");
	TMVA_tree->Branch("Ht",&Ht, "Ht/F");
	TMVA_tree->Branch("EtaStandDev",&EtaStandDev, "EtaStandDev/F");
	TMVA_tree->Branch("UnweightedEta",&UnweightedEta, "UnweightedEta/F");
	TMVA_tree->Branch("EvntShpCircularity",&EvntShpCircularity, "EvntShpCircularity/F");
	TMVA_tree->Branch("alpha_j",&alpha_j, "alpha_j/F");
	TMVA_tree->Branch("qtb1",&qtb1, "qtb1/F");
	TMVA_tree->Branch("nSV",&nSV, "nSV/I");
	TMVA_tree->Branch("Trigweight",&Trigweight, "Trigweight/F");
	TMVA_tree->Branch("DetaJJ",&DetaJJ, "DetaJJ/F");
	TMVA_tree->Branch("jetCHF0",&jetCHF[0], "jetCHF0/F");
	TMVA_tree->Branch("jetCHF1",&jetCHF[1], "jetCHF1/F");
	TMVA_tree->Branch("mu1pt",&muonPt[1], "mu1pt/F");
	TMVA_tree->Branch("muonPFiso0",&muonPFiso[0], "muonPFiso0/F");
	TMVA_tree->Branch("muonPFiso1",&muonPFiso[1], "muonPFiso1/F");
	TMVA_tree->Branch("DphiJJ",&DphiJJ, "DphiJJ/F");
	TMVA_tree->Branch("RMS_eta",&RMS_eta, "RMS_eta/F");
	TMVA_tree->Branch("PtbalZH",&PtbalZH, "PtbalZH/F");
	TMVA_tree->Branch("EventPt",&EventPt, "EventPt/F");
	TMVA_tree->Branch("Angle",&Angle, "Angle/F");
	TMVA_tree->Branch("Centrality",&Centrality, "Centrality/F");
	TMVA_tree->Branch("MET",&MET, "MET/F");
	TMVA_tree->Branch("EvntShpAplanarity",&EvntShpAplanarity, "EvntShpAplanarity/F");
	TMVA_tree->Branch("EvntShpSphericity",&EvntShpSphericity, "EvntShpSphericity/F");
	TMVA_tree->Branch("EvntShpIsotropy",&EvntShpIsotropy, "EvntShpIsotropy/F");


	TTree *BDT_tree = new TTree("BDT_tree","Tree for BDT output");
	BDT_tree->Branch("nJets",&nJets, "nJets/I");
	BDT_tree->Branch("CSV0",&CSV0, "CSV0/F");
	BDT_tree->Branch("CSV1",&CSV1, "CSV1/F");
	BDT_tree->Branch("Zmass",&Zmass, "Zmass/F");
	BDT_tree->Branch("Hmass",&Hmass, "Hmass/F");
	BDT_tree->Branch("DeltaPhiHV",&DeltaPhiHV, "DeltaPhiHV/F");
	BDT_tree->Branch("Hpt",&Hpt, "Hpt/F");
	BDT_tree->Branch("Zpt",&Zpt, "Zpt/F");
	BDT_tree->Branch("mu0pt",&mu0pt, "mu0pt/F");
	BDT_tree->Branch("Ht",&Ht, "Ht/F");
	BDT_tree->Branch("EtaStandDev",&EtaStandDev, "EtaStandDev/F");
	BDT_tree->Branch("UnweightedEta",&UnweightedEta, "UnweightedEta/F");
	BDT_tree->Branch("EvntShpCircularity",&EvntShpCircularity, "EvntShpCircularity/F");
	BDT_tree->Branch("alpha_j",&alpha_j, "alpha_j/F");
	BDT_tree->Branch("qtb1",&qtb1, "qtb1/F");
	BDT_tree->Branch("nSV",&nSV, "nSV/I");
	BDT_tree->Branch("Trigweight",&Trigweight, "Trigweight/F");
	BDT_tree->Branch("DetaJJ",&DetaJJ, "DetaJJ/F");
	BDT_tree->Branch("jetCHF0",&jetCHF[0], "jetCHF[0]/F");
	BDT_tree->Branch("jetCHF1",&jetCHF[1], "jetCHF[1]/F");
	BDT_tree->Branch("mu1pt",&muonPt[1], "mu1pt/F");
	BDT_tree->Branch("muonPFiso0",&muonPFiso[0], "muonPFiso0/F");
	BDT_tree->Branch("muonPFiso1",&muonPFiso[1], "muonPFiso1/F");
	BDT_tree->Branch("DphiJJ",&DphiJJ, "DphiJJ/F");
	BDT_tree->Branch("RMS_eta",&RMS_eta, "RMS_eta/F");
	BDT_tree->Branch("PtbalZH",&PtbalZH, "PtbalZH/F");
	BDT_tree->Branch("EventPt",&EventPt, "EventPt/F");
	BDT_tree->Branch("Angle",&Angle, "Angle/F");
	BDT_tree->Branch("Centrality",&Centrality, "Centrality/F");
	BDT_tree->Branch("MET",&MET, "MET/F");
	BDT_tree->Branch("EvntShpAplanarity",&EvntShpAplanarity, "EvntShpAplanarity/F");
	BDT_tree->Branch("EvntShpSphericity",&EvntShpSphericity, "EvntShpSphericity/F");
	BDT_tree->Branch("EvntShpIsotropy",&EvntShpIsotropy, "EvntShpIsotropy/F");
	

	bool debug = false;	
	
    // Here one can declare histograms
    // In compiled C++, it possible not to use pointer type (TH1F*) and
    // use the non-pointer type (TH1F) of ROOT object.
	
	//Declare Histograms
	
    TH1F hallhJet_pt  ("hallhJet_pt","Pt of all jets in event",		150, 0.0, 150);
	TH1F hCutFlow	("hCutFlow",  "Selection",					8, 0, 8);

	//NotThisCut
	TH1F hPtjj_NotThisCut		("hPtjj_NotThisCut","Pt of two b jets with highest CSV Not This Cut", 150, 0.0, 300);
	TH1F hPtmumu_NotThisCut	("hPtmumu_NotThisCut","Pt of two muons with highest pt Not This Cut", 150, 0.0, 300);
	TH1F hCSV0_NotThisCut		("hCSV0_NotThisCut","Jet with highest CSV Not This Cut",			50, -0.1, 1.5);
	TH1F hCSV1_NotThisCut		("hCSV1_NotThisCut","Jet with second highest CSV Not This Cut",		50, -0.1, 1.5);
	TH1F hdphiVH_NotThisCut	("hdphiVH_NotThisCut","Delta phi between Z and Higgs Not This Cut", 50, -0.1, 4.5);
	TH1F hMmumu_NotThisCut		("hMmumu_NotThisCut",  "Invariant Mass of two muons Not This Cut",	75, 0, 200);
	TH1F hNaj_NotThisCut		("hNaj_NotThisCut",  "Number of Additional Jets Not This Cut",		13, -2.5, 10.5);

	if (debug) std::cout << "all histograms declared " << std::endl;
	
    // Loop over all events.
    // For other methods to access event/navigate through the sample,
    // see the documentation of RootTreeReader class.
	int event =0, Npreselect =0, NpT_jj =0, NpT_Z =0, NCSV0 =0, NCSV1 =0;
	int NdPhi =0, N_Naj =0, N_Mjj =0, Ntree = 0, NTMVAtree =0, NBDTtree =0;
	//	double Zmass = 91.1976;
	
    do {
		event++;
		if (sample.Vtype == 0){
			if (debug) std::cout << "entered event loop " << event << std::endl;
			std::vector< std::pair<size_t,double> > indexedJetPt;
			std::vector<size_t> PtSortedJetIndex;
			std::vector< std::pair<size_t,double> > indexedJetCSV;
			std::vector<size_t> CSVSortedJetIndex;
			std::vector< std::pair<size_t,double> > indexedPt;
			std::vector<size_t> PtSortedIndex;
			
			// Analysis loop.
			// One can access the ntuple leaves directly from sample object
			
			for (int i=0; i < 10; ++i) {
			jetPt[i] = -99.99;
			jetEta[i] = -99.99;
			jetPhi[i] = -99.99;
			jetCSV[i] = -99.99;
			jetCHF[i] = -99.99;
			muonPt[i] = -99.99;
			muonEta[i] = -99.99;
			muonPhi[i] = -99.99;
			muonPFiso[i] = -99.99;
			}
			nJets =0, nSV =-99, nMuons = 0, Na_mu= 0, nPV= -99, MET= -99.99, Naj= -99;
			CSV0 = -1.0, CSV1 = -1.0, Zmass = -99.99, Hmass = -99.99, DeltaPhiHV = -99.99, Hpt = -99.99, Zpt = -99.99;
			mu0pt = -99.99, Ht = -99.99, EtaStandDev = -99.99, UnweightedEta = -99.99, EvntShpCircularity = -99.99;
			alpha_j = -99.99, qtb1 = 0.0, DphiJJ = -99.99, Trigweight = 1.0;
			RMS_eta = -99.99, PtbalZH= -99.99, EventPt= -99.99, Angle= -99.99, Centrality = -99.99;
			qtmu1 = 0.0, alpha_mu = -99.99;
			EvntShpAplanarity = -99.99, EvntShpSphericity = -99.99, EvntShpIsotropy = -99.99;
			DetaJJ = -99.99;

			double StandDevEta[4];
			double AverageEta = -99.99;
			TLorentzVector FirstJet;
			TLorentzVector SecondJet;
			TLorentzVector FirstMuon;
			TLorentzVector SecondMuon;
			TLorentzVector Higgs;
			TLorentzVector ZBoson;
			TLorentzVector HVsystem;
			double plmu1 = 0.0, plmu2 = 0.0 , qtmu2 = 0.0;
			double plb1 = 0.0, plb2 = 0.0 , qtb2 = 0.0;
			
			
			
			for (int k=0;k<sample.nhJets;k++){
				if (debug) cout << "for "<< k << "th pass" << endl;
				indexedJetPt.push_back(std::pair<size_t,double>(k,(double) sample.hJet_pt[k]));
				indexedPt.push_back(std::pair<size_t,double>(k,(double) sample.hJet_pt[k]));
				hallhJet_pt.Fill(sample.hJet_pt[k]); 
				nJets++;
				indexedJetCSV.push_back(std::pair<size_t,double>(k,(double) sample.hJet_csv[k]));
				if (debug)	cout << "CSV discriminator: " << sample.hJet_csv[k] << endl;
				
			}// end k for loop
			for (int a=0;a<sample.naJets;a++){
				hallhJet_pt.Fill(sample.aJet_pt[a]);
				nJets++; 
				jetPt[sample.nhJets+a] = sample.aJet_pt[a];
				jetEta[sample.nhJets+a] = sample.aJet_eta[a];
				jetPhi[sample.nhJets+a] = sample.aJet_phi[a];
				jetCSV[sample.nhJets+a] = sample.aJet_csv[a];
				jetCHF[sample.nhJets+a] = sample.aJet_chf[a];
			}
			indexedPt.push_back(std::pair<size_t,double>(2,(double) sample.vLepton_pt[0]));
			indexedPt.push_back(std::pair<size_t,double>(3,(double) sample.vLepton_pt[1]));
			
			
			if (event == 13) for (size_t i = 0 ; (i != indexedPt.size()) ; ++i) {cout << "Unsorted pt of objects: " << indexedPt[i].second << endl; }
			
			std::sort(indexedJetPt.begin(),indexedJetPt.end(),::IndexedQuantityGreaterThan<double>);
			std::sort(indexedJetCSV.begin(),indexedJetCSV.end(),::IndexedQuantityGreaterThan<double>);
			std::sort(indexedPt.begin(),indexedPt.end(),::IndexedQuantityGreaterThan<double>);
			
			for (size_t i = 0 ; (i != indexedJetPt.size()) ; ++i) {   PtSortedJetIndex.push_back(indexedJetPt[i].first);        }
			for (size_t i = 0 ; (i != indexedJetCSV.size()) ; ++i) {   CSVSortedJetIndex.push_back(indexedJetCSV[i].first);        }
			for (size_t i = 0 ; (i != indexedPt.size()) ; ++i) {   PtSortedIndex.push_back(indexedPt[i].first);        }
			
			if (event == 13) for (size_t i = 0 ; (i != indexedPt.size()) ; ++i) {cout << "Sorted pt of objects: " << indexedPt[i].second << endl; }
			CSV0 = sample.hJet_csv[indexedJetCSV[0].first];
			jetPt[0] = sample.hJet_pt[indexedJetCSV[0].first];
			jetEta[0] = sample.hJet_eta[indexedJetCSV[0].first];
			jetPhi[0] = sample.hJet_phi[indexedJetCSV[0].first];
			jetCSV[0] = sample.hJet_csv[indexedJetCSV[0].first];
			jetCHF[0] = sample.hJet_chf[indexedJetCSV[0].first];	
			if (sample.nhJets > 1) { 
				CSV1 = sample.hJet_csv[indexedJetCSV[1].first];	
				jetPt[1] = sample.hJet_pt[indexedJetCSV[1].first];
				jetEta[1] = sample.hJet_eta[indexedJetCSV[1].first];
				jetPhi[1] = sample.hJet_phi[indexedJetCSV[1].first];
				jetCSV[1] = sample.hJet_csv[indexedJetCSV[1].first];
				jetCHF[1] = sample.hJet_chf[indexedJetCSV[1].first];	
				DetaJJ = sample.hJet_eta[indexedJetCSV[0].first]-sample.hJet_eta[indexedJetCSV[1].first];				
				if (CSV1 > -1) { 
					Hmass = sample.H_mass;
					Hpt = sample.H_pt;
					JetDistributions("allEvts", weight);
				}//end csv non -1 requirement
			}// end njet requirement			
			
			nSV = sample.nSvs;
			nPV = sample.nPVs;
			Naj = sample.naJets;
			MET = sample.MET_sumet;

			Trigweight = sample.weightTrig;
			if (debug) cout << "halfway through filling histos" << endl;
			
			if ((sample.vLepton_type[0] == 13) && (sample.vLepton_type[1] == 13)){
				for (int z = 0; z< sample.nvlep;z++){
				muonPt[z] = sample.vLepton_pt[z];
				muonEta[z] = sample.vLepton_eta[z];
				muonPhi[z] = sample.vLepton_phi[z];
				muonPFiso[z] = sample.vLepton_pfCombRelIso[z];
				nMuons++;
				}
				for (int zz = 0; zz< sample.nalep;zz++){
					if (sample.aLepton_type[zz] == 13){
					muonPt[zz+sample.nvlep] = sample.aLepton_pt[zz];
					muonEta[zz+sample.nvlep] = sample.aLepton_eta[zz];
					muonPhi[zz+sample.nvlep] = sample.aLepton_phi[zz];
					muonPFiso[zz+sample.nvlep] = sample.aLepton_pfCombRelIso[zz];
					nMuons++;
					Na_mu++;
					}}
				if (sample.nvlep > 1) {
					Zmass = sample.V_mass;
					Zpt = sample.V_pt;
				}
				MuonDistributions("allEvts", weight);
				if( (sample.nvlep > 1) && (sample.nhJets > 1)) {
					RMS_eta = sample.vLepton_eta[1]*sample.vLepton_eta[1]+sample.vLepton_eta[0]*sample.vLepton_eta[0]+sample.hJet_eta[indexedJetCSV[0].first]*sample.hJet_eta[indexedJetCSV[0].first]+sample.hJet_eta[indexedJetCSV[1].first]*sample.hJet_eta[indexedJetCSV[1].first];
					RMS_eta = sqrt(RMS_eta/4);
					StandDevEta[0] =sample.vLepton_eta[1];
					StandDevEta[1] =sample.hJet_eta[indexedJetCSV[0].first];
					StandDevEta[2] =sample.vLepton_eta[0];
					StandDevEta[3] =sample.hJet_eta[indexedJetCSV[1].first];
					EtaStandDev = TMath::RMS(4,StandDevEta);
					AverageEta = sample.vLepton_eta[0]+sample.vLepton_eta[1]+sample.hJet_eta[indexedJetCSV[0].first];
					AverageEta = (AverageEta+sample.hJet_eta[indexedJetCSV[1].first])/4;
					UnweightedEta = (sample.vLepton_eta[0]-AverageEta)*(sample.vLepton_eta[0]-AverageEta);
					UnweightedEta = UnweightedEta + (sample.vLepton_eta[1]-AverageEta)*(sample.vLepton_eta[1]-AverageEta);
					UnweightedEta = UnweightedEta + (sample.hJet_eta[indexedJetCSV[0].first]-AverageEta)*(sample.hJet_eta[indexedJetCSV[0].first]-AverageEta);
					UnweightedEta = UnweightedEta + (sample.hJet_eta[indexedJetCSV[1].first]-AverageEta)*(sample.hJet_eta[indexedJetCSV[1].first]-AverageEta);
					
					DeltaPhiHV = sample.HVdPhi;			
					PtbalZH = (Hpt-Zpt);
					mu0pt = sample.vLepton_pt[0];
					
					FirstJet.SetPtEtaPhiM(sample.hJet_pt[indexedJetCSV[0].first],sample.hJet_eta[indexedJetCSV[0].first],sample.hJet_phi[indexedJetCSV[0].first],4.2/1000.0);
					SecondJet.SetPtEtaPhiM(sample.hJet_pt[indexedJetCSV[1].first],sample.hJet_eta[indexedJetCSV[1].first],sample.hJet_phi[indexedJetCSV[1].first],4.2/1000.0);
					FirstMuon.SetPtEtaPhiM(mu0pt,sample.vLepton_eta[0],sample.vLepton_phi[0],105.65836668/1000.0);
					SecondMuon.SetPtEtaPhiM(sample.vLepton_pt[1],sample.vLepton_eta[1],sample.vLepton_phi[1],105.65836668/1000.0);
					Higgs = FirstJet+SecondJet;
					ZBoson = FirstMuon+SecondMuon;
					
					Ht = FirstJet.Pt()+SecondJet.Pt()+FirstMuon.Pt()+SecondMuon.Pt();
					DphiJJ = FirstJet.DeltaPhi(SecondJet);
					
					Centrality = (indexedPt[2].second+indexedPt[3].second)/Ht;
					Angle = Higgs.Angle(ZBoson.Vect());
					HVsystem = 	Higgs+ZBoson;
					EventPt = HVsystem.Pt();
					
					//Armenteros-Podolansky Plots
					TVector3 VectZ = ZBoson.Vect();
					TVector3 Unit_Z = VectZ.Unit();
					plmu1 = Unit_Z.Dot(FirstMuon.Vect());
					plmu2 = Unit_Z.Dot(SecondMuon.Vect());
					alpha_mu = (plmu1-plmu2)/(plmu1+plmu2);
					qtmu1= FirstMuon.Perp(ZBoson.Vect()); 
					qtmu2 = SecondMuon.Perp(ZBoson.Vect());
					
					TVector3 VectH = Higgs.Vect();
					TVector3 Unit_H = VectH.Unit();
					plb1 = Unit_H.Dot(FirstJet.Vect());
					plb2 = Unit_H.Dot(SecondJet.Vect());
					if (FirstJet.Pt()>SecondJet.Pt()){
						alpha_j = (plb1-plb2)/(plb1+plb2);
					}
					if (FirstJet.Pt()<SecondJet.Pt()){
						alpha_j = (plb2-plb1)/(plb1+plb2);
					}
					qtb1= FirstJet.Perp(Higgs.Vect()); 
					qtb2 = SecondJet.Perp(Higgs.Vect());
					TH2FDistributions("allEvts", weight);
					
					//EventShapeVariables
					std::vector<math::XYZVector> particles;
					TVector3 TVect_particle = FirstMuon.Vect();
					math::XYZVector Vparticle1(TVect_particle.X(), TVect_particle.Y(), TVect_particle.Z());
					particles.push_back(Vparticle1);
					TVect_particle = SecondMuon.Vect();
					math::XYZVector Vparticle2(TVect_particle.X(), TVect_particle.Y(), TVect_particle.Z());
					particles.push_back(Vparticle2);
					TVect_particle = FirstJet.Vect();
					math::XYZVector Vparticle3(TVect_particle.X(), TVect_particle.Y(), TVect_particle.Z());
					particles.push_back(Vparticle3);
					TVect_particle = SecondJet.Vect();
					math::XYZVector Vparticle4(TVect_particle.X(), TVect_particle.Y(), TVect_particle.Z());
					particles.push_back(Vparticle4);
					
					EventShapeVariables eventshape(particles);
					EvntShpAplanarity = eventshape.aplanarity();
					EvntShpSphericity = eventshape.sphericity();
					EvntShpCircularity = eventshape.circularity();
					EvntShpIsotropy = eventshape.isotropy();
					
					EventShapeDistributions("allEvts", weight);
				}// two muons and two jets requirement
				
				
				if (debug) cout << "done filling histograms for event: " << event << endl;
				
				//Calculate cut efficiencies
				if ((sample.nhJets > 1) && (sample.nvlep > 1)){
					//all samples are of typle Zmumu
					if (CSV1 > -1 ){
						if ((Hpt > 75) && (Zpt>75) && (CSV0 > 0.244)) { //One b tag
							if (Zmass >70 && Zmass<110) { //Dan's Selection
					Ntree++;
							if (Ntree%2){
							TMVA_tree->Fill();
							NTMVAtree++;}else{
							BDT_tree->Fill();
							NBDTtree++;
							}
							JetDistributions("tree", weight);
							MuonDistributions("tree", weight);
							TH2FDistributions("tree", weight);
							EventShapeDistributions("tree", weight);
							EventDistributions("tree", weight); 
						} //CSV1 > -1
					}// Zmass requirement
				}//Cuts for trees
						if ( (CSV0 > 0.898) && (CSV1>0.5) && (sample.naJets < 1) && ((DeltaPhiHV>=2.90)||(DeltaPhiHV<-2.90)) ){
						/*	if (sample.hJet_pt[indexedJetCSV[0].first] >=20 && sample.hJet_pt[indexedJetCSV[1].first] >=20 ) {
						MuonDistributions("NoPtCut", weight);
							}// end jet pt cuts
							if (sample.vLepton_pt[1] > 20 && (Zmass >75 && Zmass<105) ) {
								JetDistributions("NoPtCut", weight);
							}//end muon pt cuts
						TH2FDistributions("NoPtCut", weight);
						EventShapeDistributions("NoPtCut", weight);
						EventDistributions("NoPtCut", weight);
						}// end non pt cuts */
						if (sample.hJet_pt[indexedJetCSV[0].first] >=20 && sample.hJet_pt[indexedJetCSV[1].first] >=20 && sample.vLepton_pt[1] > 20) {
							if ((Hpt > 100) && (Zpt>100) && (CSV0 > 0.898) && (CSV1>0.5) && ((DeltaPhiHV>=2.90)||(DeltaPhiHV<-2.90)) && (sample.naJets < 1) ) hMmumu_NotThisCut.Fill(Zmass, weight);
							if (Zmass >75 && Zmass<105){
								Npreselect++;
							/*	JetDistributions("preSelect", weight);
								MuonDistributions("preSelect", weight);
								TH2FDistributions("preSelect", weight);
								EventShapeDistributions("preSelect", weight);
								EventDistributions("preSelect", weight);*/
								if ((Zpt>100) && (CSV0 > 0.898) && (CSV1>0.5) && ((DeltaPhiHV>=2.90)||(DeltaPhiHV<-2.90)) && (sample.naJets < 1) ) hPtjj_NotThisCut.Fill(Hpt, weight);
								if(Hpt > 100) {
									NpT_jj++;
									//JetDistributions("pT_jj", weight);
									if ((CSV0 > 0.898) && (CSV1>0.5) && ((DeltaPhiHV>=2.90)||(DeltaPhiHV<-2.90)) && (sample.naJets < 1) ) hPtmumu_NotThisCut.Fill(Zpt, weight);				   
									if(Zpt > 100) {
										NpT_Z++;
										//JetDistributions("pT_Z", weight);
										if ((CSV1>0.5) && ((DeltaPhiHV>=2.90)||(DeltaPhiHV<-2.90)) && (sample.naJets < 1) ) hCSV0_NotThisCut.Fill(CSV0, weight);
										if(CSV0 > 0.898){
											NCSV0++;
											//JetDistributions("CSV0", weight);
											if (((DeltaPhiHV>=2.90)||(DeltaPhiHV<-2.90)) && (sample.naJets < 1) ) hCSV1_NotThisCut.Fill(CSV1, weight);
											if(CSV1>0.5){
												NCSV1++;
												//JetDistributions("CSV1", weight);
												if (sample.naJets < 1) {
												/*	JetDistributions("NoDphiCut", weight);
													MuonDistributions("NoDphiCut", weight);
													TH2FDistributions("NoDphiCut", weight);
													EventShapeDistributions("NoDphiCut", weight);
													EventDistributions("NoDphiCut", weight);*/
												}
												if((DeltaPhiHV>=2.90)||(DeltaPhiHV<-2.90)){
													NdPhi++;
													//JetDistributions("DPhi", weight);
													hNaj_NotThisCut.Fill(sample.naJets);
													if(sample.naJets < 2){
														N_Naj++;
														//JetDistributions("Naj", weight);
													   /* JetDistributions("aftercuts", weight);
														MuonDistributions("aftercuts", weight);
														TH2FDistributions("aftercuts", weight);
														EventShapeDistributions("aftercuts", weight);
														EventDistributions("aftercuts", weight); */
														if((Hmass>=95)&&(Hmass<=125)){
														N_Mjj++; 
													/*	JetDistributions("Mjj", weight);
														MuonDistributions("Mjj", weight);
														TH2FDistributions("Mjj", weight);
														EventShapeDistributions("Mjj", weight);
														EventDistributions("Mjj", weight);*/
														}//higgs mass window
													}//number of additional Jets cut
												}//Delta Phi Cut
											}//CSV1 cut
										}//CSV0 cut
									}// Z pt cut
								}// dijet pt cut
							}// preselection z mass requirement
						}//preselection pt requirement
					}// preselection CSV not -1
				}// preselection number of muons/jets requirement
			}//end requirement muon not electron
			EventDistributions("allEvts", weight);
		}//end requirement Zmumu event
    } while (sample.nextEvent());
	if (debug){	std::cout << "Number of events " << event << endl;}
	std::cout << "BDT tree: " << NBDTtree << endl;
	std::cout << "TMVA tree: " << NTMVAtree << endl;
	std::cout << "PreSelection: " << Npreselect << endl;
	std::cout << "pT_jj: " << NpT_jj << endl;
	std::cout << "pT_Z: " << NpT_Z << endl;
	std::cout << "CSV0: " << NCSV0 << endl;
	std::cout << "CSV1: " << NCSV1 << endl;
	std::cout << "DPhi: " << NdPhi << endl;
	std::cout << "Naj: " << N_Naj << endl;
	std::cout << "Mjj: " << N_Mjj << endl;
	std::cout << "Number of Events Passing the Selection: " << N_Mjj << endl;
	
    // Here one can create canvas, draw and print the histogram.
    TCanvas c1("c1","c1");
    c1.cd();
    c1.SetFillColor(kWhite);
    string suffixps = ".gif";
	
    c1.Clear(); // don't create a new canvas
    hallhJet_pt.Draw();
    c1.Print((directory+"/PtAllJets"+suffixps).c_str());
	
	c1.Clear(); // don't create a new canvas
    hPtjj_NotThisCut.Draw();
    c1.Print((directory+"/Ptjj_NotThisCut"+suffixps).c_str());

	c1.Clear(); // don't create a new canvas
    hPtmumu_NotThisCut.Draw();
    c1.Print((directory+"/Ptmumu_NotThisCut"+suffixps).c_str());

	c1.Clear(); // don't create a new canvas
    hCSV0_NotThisCut.Draw();
    c1.Print((directory+"/CSV0_NotThisCut"+suffixps).c_str());
	
	c1.Clear(); // don't create a new canvas
    hCSV1_NotThisCut.Draw();
    c1.Print((directory+"/CSV1_NotThisCut"+suffixps).c_str());

	c1.Clear(); // don't create a new canvas
	hNaj_NotThisCut.Draw();
	c1.Print((directory+"/Naj_NotThisCut"+suffixps).c_str());
		
	c1.Clear(); // don't create a new canvas
	hMmumu_NotThisCut.Draw();
	c1.Print((directory+"/Mmumu_NotThisCut"+suffixps).c_str());

	c1.Clear(); // don't create a new canvas
	hdphiVH_NotThisCut.Draw();
	c1.Print((directory+"/dphiVH_NotThisCut"+suffixps).c_str());		
			
	c1.Clear(); // don't create a new canvas
	/*	hCutFlow.GetXaxis().SetBinLabel(1,"PreSelection");
	 hCutFlow.GetXaxis().SetBinLabel(2,"p_{T}(jj)");
	 hCutFlow.GetXaxis().SetBinLabel(3,"p_{T}(Z)");
	 hCutFlow.GetXaxis().SetBinLabel(4,"CSV0");
	 hCutFlow.GetXaxis().SetBinLabel(5,"CSV1");
	 hCutFlow.GetXaxis().SetBinLabel(6,"#Delta#phi");
	 hCutFlow.GetXaxis().SetBinLabel(7,"N_{aj}");
	 hCutFlow.GetXaxis().SetBinLabel(8,"M_{jj}");*/
	hCutFlow.SetBinContent(1,Npreselect );
	hCutFlow.SetBinContent(2,NpT_jj );
	hCutFlow.SetBinContent(3,NpT_Z );
	hCutFlow.SetBinContent(4,NCSV0 );
	hCutFlow.SetBinContent(5,NCSV1 );
	hCutFlow.SetBinContent(6,NdPhi );
	hCutFlow.SetBinContent(7,N_Naj );
	hCutFlow.SetBinContent(8,N_Mjj );
	hCutFlow.Draw();
	c1.Print((directory+"/CutFlow"+suffixps).c_str());
		
	TMVA_tree->Write();

	////////
	// write histos to file:
	////////
	for (map<std::string,TH1*>::iterator it=histmap.begin(); it!=histmap.end();it++) {
		(*it).second->Write();
		c1.Clear();
		(*it).second->Draw();
		c1.Print((directory+"/"+(*it).first+suffixps).c_str());
if (debug)		cout << "This is the histmap string: " << (*it).first << endl;
	}
	for (map<string,TH2*>::iterator it2=bidimhistmap.begin(); it2!=bidimhistmap.end();it2++) {
		(*it2).second->Write();
		c1.Clear();
		(*it2).second->Draw();
		c1.Print((directory+"/"+(*it2).first+suffixps).c_str());
		delete (*it2).second;
	}
		
	
    // Write and Close the output file.
    ofile.Write();
    ofile.Close();
	
    return 0;
}


double SetWeight( std::string filename){
	double SampleWeight = 1.0;
	if (findString(filename, "ZH_ZToLL_HToBB_M-115")){ 
	SampleWeight = 0.4107*0.704*0.101*10000/219999.0;
	cout << "found ZH_ZToLL_HToBB_M-115 string" << endl;}
	if (findString(filename, "DYJetsToLL_PtZ")){ SampleWeight = 37.90933333*10000/ 1137280.0;}
	if (findString(filename, "DYJetsToLL_TuneZ2")){ 
	SampleWeight = 3151.864553*10000/ 36228691.0;
	cout << " found DYJetsToLL_TuneZ2 string " << endl;}
	if (findString(filename, "120to170")){ SampleWeight = 115613.7358*10000/ 6127528.0;}
	if (findString(filename, "170to300")){ SampleWeight = 24297.5*10000/ 6220160.0;}
	if (findString(filename, "300to470")){ SampleWeight = 1169.576182*10000/ 6432669.0;}
	if (findString(filename, "470to600")){ SampleWeight = 70.21630282*10000/ 3598283.0;}
	if (findString(filename, "80to120")){ SampleWeight = 823744.5*10000/ 6199951.0;}
	if (findString(filename, "TTJets")){ SampleWeight = 157.5*10000/ 3701947.0;}
	if (findString(filename, "T_TuneZ2_s")){ SampleWeight = 3.19*10000/ 259971.0;}
	if (findString(filename, "T_TuneZ2_t-channel")){ SampleWeight = 41.92*10000/ 3900171.0;}
	if (findString(filename, "T_TuneZ2_tW-channel-DR")){ SampleWeight = 7.87*10000/ 814390.0;}
	if (findString(filename, "T_TuneZ2_tW-channel-DS")){ SampleWeight = 7.87*10000/ 795379.0;}
	if (findString(filename, "Tbar_TuneZ2_s")){ SampleWeight = 1.44*10000/ 137980.0;}
	if (findString(filename, "Tbar_TuneZ2_t-channel")){ SampleWeight = 22.65*10000/ 1944826.0;}
	if (findString(filename, "Tbar_TuneZ2_tW-channel-DR")){ SampleWeight = 7.87*10000/ 809984.0;}
	if (findString(filename, "Tbar_TuneZ2_tW-channel-DS")){ SampleWeight = 7.87*10000/ 785246.0;}
	if (findString(filename, "WJetsToLNu_Pt-100")){ SampleWeight = 732.6421333*10000/ 21955936.0;}
	if (findString(filename, "TestWJetsToLNu_TuneZ2")){ SampleWeight = 10438.0*3.0*10000/ 81345381.0;}
	if (findString(filename, "WW")){ SampleWeight = 47*10000/ 4225916.0;}
	if (findString(filename, "WZ")){ SampleWeight = 18.2*10000/ 18.2*10000/4265243.0;}
	if (findString(filename, "ZJetsToLL")){ SampleWeight = 3048*10000/2750000.0;}
	if (findString(filename, "ZJetsToNuNu")){ SampleWeight = 1.3*25.8*10000/ 5500000.0;}
	if (findString(filename, "ZZ")){ 
	SampleWeight = 7.41*10000/4191045.0;
	cout << "found ZZ string" << endl;
	}
	
	return SampleWeight;
}

bool findString(std::string strToSearch, std::string strPattern){
	size_t found;
	
	bool foundStr = false;
	
	found=strToSearch.find(strPattern);
    if (found!=string::npos)
		foundStr = true;
	
	return foundStr;		
}

void fillhisto(string histname,double val, double hweight, string title, int nbins, double min, double max) {
	if (!histmap[histname]) {
		if (nbins==0) {
			cout << " ERROR: must define histogram " << histname << " properties on first call of fillhisto(...)" << endl;
			exit(1);
		}
		
		histmap[histname]=new TH1F(histname.c_str(),title.c_str(),nbins,min,max);
		
	}
	histmap[histname]->Fill(val, hweight);
}
void fill2Dhisto(string histname,double valx, double valy, double weight, string title, int nxbins, double xmin, double xmax, int nybins, double ymin, double ymax) {
	if (!bidimhistmap[histname]) {
		if (nxbins==0 || nybins==0) {
			//cout << " ERROR: must define bi dim histogram " << histname << " properties on first call of fillbidimhisto(...)" << endl;
			exit(1);
		}
		bidimhistmap[histname]=new TH2F(histname.c_str(),title.c_str(),nxbins,xmin,xmax,nybins,ymin,ymax);
		
	}
	bidimhistmap[histname]->Fill(valx,valy,weight);
}// fill2Dhisto


void JetDistributions(string cut, double ph_weight){
		
		for (int i=0; i != nJets && i < 10; ++i) {
			string jetpT = Form("hPtb%i_", i) + cut;
			string jeteta = Form("hEtab%i_", i) + cut;
				string jetphi = Form("hPhib%i_", i) + cut;
			string jetcsv = Form("hCSV%i_", i) + cut;
			string jetchf = Form("hCHFb%i_", i) + cut;
			
			fillhisto(jetpT, jetPt[i], ph_weight, "jet pT", 150, 0.0, 150);
			fillhisto(jeteta, jetEta[i], ph_weight, "jet #eta", 75, -3, 3.5);
			fillhisto(jetphi, jetPhi[i], ph_weight, "jet #phi", 75, -3.25, 4.5);
			fillhisto(jetcsv, jetCSV[i], ph_weight, "jet CSV", 50, -0.1, 1.5);
			fillhisto(jetchf, jetCHF[i], ph_weight, "charged Hadron Energy Fraction b1", 101, 0.0, 1.2);
		}
		
		string mjj = "hMjj_" + cut;
	string ptjj = "hPtjj_" +cut;
	string detajj = "hdetaJJ_" +cut;
	string dphijj = "hdphiJJ_vect_" +cut;
		
	fillhisto(mjj, Hmass, ph_weight, "Invariant Mass of two Jets", 200, 0, 200);
	fillhisto(ptjj, Hpt, ph_weight, "Pt of two b jets with highest CSV", 150, 0, 300);
	fillhisto(detajj, DetaJJ, ph_weight, "Delta eta between two jets", 101, -3, 3);
	fillhisto(dphijj, DphiJJ, ph_weight, "Delta phi between two jets",  71, -3.5, 4);
		
}// JetDistributions


void MuonDistributions(string cut, double ph_weight){
	
	for (int i=0; i != nMuons && i < 10; ++i) {
		string muonpT = Form("hPtmu%i_", i) + cut;
		string muoneta = Form("hEtamu%i_", i) + cut;
		string muonphi = Form("hPhimu%i_", i) + cut;
		string muonPFiso = Form("hPFRelIsomu%i_", i) + cut;
		
		fillhisto(muonpT, muonPt[i], ph_weight, "muon pT", 200, 0.0, 200);
		fillhisto(muoneta, muonEta[i], ph_weight, "muon #eta", 75, -3, 3.5);
		fillhisto(muonphi, muonPhi[i], ph_weight, "muon #phi", 75, -3.25, 4.5);
		fillhisto(muonPFiso, muonPFiso[i], ph_weight, "PF Rel Iso of muon", 101, 0.0001, 0.2);
	}
	
	string mmumu = "hMmumu_" + cut;
	string ptmumu = "hPtmumu_" +cut;
	
	fillhisto(mmumu, Zmass, ph_weight, "Invariant Mass of two Muons", 200, 0, 200);
	fillhisto(ptmumu, Zpt, ph_weight, "Pt of two muons", 150, 0, 300);
	
}// MuonDistributions	

void TH2FDistributions(string cut, double ph_weight){
		
	string hDphiDetajj = "hDphiDetajj_" + cut;
	string hqtvsalphaJJ = "hqtvsalphaJJ_" +cut;
	string hqtvsalphaZ = "hqtvsalphaZ_" +cut;
	
	fill2Dhisto(hDphiDetajj, DphiJJ, DetaJJ, ph_weight, "#Delta#phi vs #Delta#eta JJ", 75, -5, 5, 75, -3.5, 3.5);
	fill2Dhisto(hqtvsalphaJJ, alpha_j, qtb1, ph_weight, "Armenteros-Podolansky Plot Jets", 70, -2, 2, 75, 0, 75);
	fill2Dhisto(hqtvsalphaZ, alpha_mu, qtmu1, ph_weight,  "Armenteros-Podolansky Plot Z", 70, -2, 2, 75, 0, 75);
	
}// TH2FDistributions		

void EventShapeDistributions(string cut, double ph_weight){
	
	string hPtbalZH = "hPtbalZH_" + cut;
	string hdphiVH = "hdphiVH_" +cut;
	string hRMSeta = "hRMSeta_" +cut;
	string hStaDeveta = "hStaDeveta_" +cut;
	string hUnweightedEta = "hUnweightedEta_" +cut;
	string hHt = "hHt_" +cut;
	string hCentrality = "hCentrality_" +cut;
	string hEventPt = "hEventPt_" +cut;
	string hAngle = "hAngle_" +cut;
	string hSphericity = "hSphericity_" +cut;
	string hAplanarity = "hAplanarity_" +cut;
	string hCircularity = "hCircularity_" +cut;
	string hIsotropy = "hIsotropy_" +cut;

	fillhisto(hPtbalZH, PtbalZH, ph_weight, "Pt balance of Z and H", 101, -75.0, 75);
	fillhisto(hdphiVH, DeltaPhiHV, ph_weight, "Delta phi between Z and Higgs", 50, -0.1, 4.5);
	fillhisto(hRMSeta, RMS_eta, ph_weight, "RMS Eta",		71, 0, 3);
	fillhisto(hStaDeveta, EtaStandDev, ph_weight, "Standard Deviation Eta",		71, 0, 3);
	fillhisto(hUnweightedEta, UnweightedEta, ph_weight, "Unweighted Eta",		101, 0, 10);
	fillhisto(hHt, Ht, ph_weight, "scalar sum of pt of four particles", 101, 0.0, 500);
	fillhisto(hCentrality, Centrality, ph_weight, "Centrality", 71, 0.0, 1.0);
	fillhisto(hEventPt, EventPt, ph_weight, "Pt of HV system", 50, 0.0, 100);
	fillhisto(hAngle, Angle, ph_weight, "Angle between H and Z", 101, -0.1, 4);
	fillhisto(hSphericity, EvntShpSphericity, ph_weight,"EventShapeVariables sphericity", 71, 0.0, 1);
	fillhisto(hAplanarity, EvntShpAplanarity, ph_weight, "EventShapeVariables Aplanarity", 50, 0.0, .5);
	fillhisto(hCircularity, EvntShpCircularity, ph_weight,  "EventShapeVariables circularity", 45, 0.0, 1.2);
	fillhisto(hIsotropy, EvntShpIsotropy, ph_weight,  "EventShapeVariables isotropy", 50, 0.0, 1.2);

}// EventShapeDistributions


void EventDistributions(string cut, double ph_weight){
		
	string hnJets = "hnJets_" + cut;
	string hnMuons = "hnMuons_" +cut;
	string hnSV = "hnSV_" +cut;
	string hnPV = "hnPV_" +cut;
	string hMET = "hMET_" +cut;
	string hNaj = "hNaj_" +cut;
	string Namu_ = "hNamu_" +cut;
	
	fillhisto(hnJets, nJets, ph_weight, "Number of Good Jets", 11, -0.5, 10.5);
	fillhisto(hnMuons, nMuons, ph_weight, "Number of Good Muons", 6, -0.5, 5.5);
	fillhisto(hnSV, nSV, ph_weight, "Number of Secondary Verticies", 6, -0.5, 5.5);
	fillhisto(hnPV, nPV, ph_weight, "Number of Primary Verticies", 21, -0.5, 20.5);
	fillhisto(hMET, MET, ph_weight, "Missing Et",		150, 0.0, 1200);
	fillhisto(hNaj, Naj, ph_weight, "Number of Additional Jets",		13, -2.5, 10.5);
	fillhisto(Namu_, Na_mu, ph_weight, "Number of Additional Muons",		13, -2.5, 10.5);
	
}// EventDistributions	
