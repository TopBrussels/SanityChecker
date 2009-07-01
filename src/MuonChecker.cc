#include "../interface/MuonChecker.h"


//
// constructors and destructor
//
MuonChecker::MuonChecker(const edm::ParameterSet& iConfig,string label)
{
	muonLabel_ = iConfig.getParameter<edm::InputTag>("labelMuons");
	verbose_          =   iConfig.getParameter<bool>( "verbose" );
	NbOfMu           = 0;
	NbOfGlobalMu           = 0;
	label_ = label;

        NbOfMuonId = 13;
        MuonID = new char*[NbOfMuonId];
        MuonID[0] = "All";                             // dummy options - always true
	MuonID[1] = "AllGlobalMuons";                  // checks isGlobalMuon flag
	MuonID[2] = "AllStandAloneMuons";              // checks isStandAloneMuon flag
	MuonID[3] = "AllTrackerMuons";                 // checks isTrackerMuon flag
	MuonID[4] = "TrackerMuonArbitrated";           // resolve ambiguity of sharing segments
	MuonID[5] = "AllArbitrated";                   // all muons with the tracker muon arbitrated
	MuonID[6] = "GlobalMuonPromptTight";           // global muons with tighter fit requirement
	MuonID[7] = "TMLastStationLoose";              // penetration depth loose selector
	MuonID[8] = "TMLastStationTight";              // penetration depth tight selector
	MuonID[9] = "TM2DCompatibilityLoose";          // likelihood based loose selector
	MuonID[10] = "TM2DCompatibilityTight";          // likelihood based tight selector
	MuonID[11] = "TMLastStationOptimizedLowPtLoose";// combination of TMLastStation and TMOneStation
	MuonID[12] = "TMLastStationOptimizedLowPtTight"; // combination of TMLastStation and TMOneStation
						

}

MuonChecker::~MuonChecker()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}

//
// member functions
//

// ------------ method called to for each event  ------------
void
MuonChecker::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{

   // reco::Muon
   edm::Handle<std::vector<reco::Muon> > muonHandle;
   iEvent.getByLabel(muonLabel_,muonHandle);

   //Check if branch is available  
   if (!muonHandle.isValid())
   {
     edm::LogWarning  ("NoDataFound_MC_noMuon") << "--- NoMuonsFound ---";
     return; 
   }
   const std::vector<reco::Muon> & muons = *muonHandle;


   for(std::vector<reco::Muon>::const_iterator muon_iter = muons.begin(); muon_iter!=muons.end(); ++muon_iter)
   {
	NbOfMu++;

	for(Int_t i=0;i<NbOfMuonId;i++)
	{
		if(muon_iter->isGood(reco::Muon::SelectionType(i))) histocontainer_["MuonId"]->Fill(MuonID[i],1);
	}

	if(muon_iter->isGlobalMuon())
	{
		NbOfGlobalMu++;
		histocontainer_["GlobalMuonInnerTrackNbOfValidHits"]  ->Fill(muon_iter->innerTrack()->numberOfValidHits());
		histocontainer_["GlobalMuonInnerTrackNbOfLostHits"]   ->Fill(muon_iter->innerTrack()->numberOfLostHits());
		histocontainer_["GlobalMuonInnerTrackPt"]   	  ->Fill(muon_iter->innerTrack()->pt());

		histocontainer_["GlobalMuonOuterTrackNbOfValidHits"]  ->Fill(muon_iter->outerTrack()->numberOfValidHits());
		histocontainer_["GlobalMuonOuterTrackNbOfLostHits"]   ->Fill(muon_iter->outerTrack()->numberOfLostHits());
		histocontainer_["GlobalMuonOuterTrackPt"]    	  ->Fill(muon_iter->outerTrack()->pt());

		histocontainer_["GlobalMuonGlobalTrackPt"]		  ->Fill(muon_iter->globalTrack()->pt());
		histocontainer_["GlobalMuonGlobalTrackD0"]		  ->Fill(muon_iter->globalTrack()->d0());
		histocontainer_["GlobalMuonGlobalTrackChi2"]	  ->Fill(muon_iter->globalTrack()->normalizedChi2());
	}
	
	histocontainer_["MuonCaloCompatibility"] ->Fill(muon_iter->caloCompatibility());
	histocontainer_["MuonSegmCompatibility"] ->Fill(muon_iter->segmentCompatibility());

	histocontainer_["MuonIsoR03SumPt"] ->Fill(muon_iter->isolationR03().sumPt);
	histocontainer_["MuonIsoR03emEt"]  ->Fill(muon_iter->isolationR03().emEt);
	histocontainer_["MuonIsoR03hadEt"] ->Fill(muon_iter->isolationR03().hadEt);
  	histocontainer_["MuonRelIso"]      ->Fill(muon_iter->pt()/(muon_iter->pt()+muon_iter->isolationR03().sumPt + muon_iter->isolationR03().emEt + muon_iter->isolationR03().hadEt)); 
  
  	//get the veto cone information
  	//  WARNING TO BE RECAL ...
	//histocontainer_["MuonVetoEm"]  ->Fill(muon_iter->ecalIsoDeposit()->candEnergy());
  	//histocontainer_["MuonVetoHad"] ->Fill(muon_iter->hcalIsoDeposit()->candEnergy());

	histocontainer_["MuonCharge"] ->Fill(muon_iter->charge());
   }
   histocontainer_["NofMuons"]  ->Fill(NbOfMu);
   histocontainer_["NofGlobalMuons"]  ->Fill(NbOfGlobalMu);

}


// ------------ method called once each job just before starting event loop  ------------
void 
MuonChecker::beginJob(const edm::EventSetup&)
{
  edm::Service<TFileService> fs;
  if (!fs) throw edm::Exception(edm::errors::Configuration, "TFileService missing from configuration!");

  char dirname[100];
  sprintf(dirname,"Muons_%s",label_.c_str());
  TFileDirectory subDir = fs->mkdir( dirname );
    
  histocontainer_["NofMuons"]      		            =subDir.make<TH1D>("NofMuons","Number of muons",20,0,20);histocontainer_["NofMuons"]->LabelsDeflate("X");
  histocontainer_["NofGlobalMuons"]      		            =subDir.make<TH1D>("NofGlobalMuons","Number of muons",20,0,20);histocontainer_["NofGlobalMuons"]->LabelsDeflate("X");
  histocontainer_["MuonId"]      		            =subDir.make<TH1D>("MuonId","Muon identification algorithms",13,0,13);histocontainer_["MuonId"]->LabelsDeflate("X");
  ///////////////////////////////////////////////////////////////////////////////
  // for all global muons (except the one from the semi-muonic top quark decay) :

  histocontainer_["GlobalMuonInnerTrackNbOfValidHits"]      =subDir.make<TH1D>("GlobalMuonInnerTrackNbOfValidHits","Number of valid hits in silicon fit for global muons",400,0,20);
  histocontainer_["GlobalMuonInnerTrackNbOfLostHits"]	    =subDir.make<TH1D>("GlobalMuonInnerTrackNbOfLostHits","Number of lost hits in silicon fit for global muons",400,0,20);
  histocontainer_["GlobalMuonInnerTrackPt"]	            =subDir.make<TH1D>("GlobalMuonInnerTrackPt","Inner track transverse momentum for global muons",400,0,400);

  histocontainer_["GlobalMuonOuterTrackNbOfValidHits"]      =subDir.make<TH1D>("GlobalMuonOuterTrackNbOfValidHits","Number of valid hits in silicon fit for global muons",400,0,20);
  histocontainer_["GlobalMuonOuterTrackNbOfLostHits"]       =subDir.make<TH1D>("GlobalMuonOuterTrackNbOfLostHits","Number of lost hits in silicon fit for global muons",400,0,20);
  histocontainer_["GlobalMuonOuterTrackPt"]		    =subDir.make<TH1D>("GlobalMuonOuterTrackPt","Outer track transverse momentum for global muons",400,0,400);

  histocontainer_["GlobalMuonGlobalTrackPt"]                =subDir.make<TH1D>("GlobalMuonGlobalTrackPt","Global track transverse momentum for global muons",400,0,400);
  histocontainer_["GlobalMuonGlobalTrackD0"]    	    =subDir.make<TH1D>("GlobalMuonGlobalTrackD0","Global track impact parameter for global muons",400,-0.4,0.4);
  histocontainer_["GlobalMuonGlobalTrackChi2"] 		    =subDir.make<TH1D>("GlobalMuonGlobalTrackChi2","Global track normalized #chi^{2} ",400,0,20);
  
  histocontainer_["MuonCaloCompatibility"]                  =subDir.make<TH1D>("CaloCompatibility","Value of the LR measuring the probability that the muon is calo-compatible",100,0,1);
  histocontainer_["MuonSegmCompatibility"]                  =subDir.make<TH1D>("SegmCompatibility","Value of the LR measuring the probability that the muon is segment-compatible",100,0,1);
  histocontainer_["MuonIsoR03SumPt"]   		            =subDir.make<TH1D>("MuonIsoR03SumPt","Sum of the track transverse momenta in a cone of 0.3 around the muon",200,0,20);
  histocontainer_["MuonIsoR03emEt"]			    =subDir.make<TH1D>("MuonIsoR03emEt","Sum of the electromagnetic transverse energy in a cone of 0.3 around the muon",200,0,20);
  histocontainer_["MuonIsoR03hadEt"]		            =subDir.make<TH1D>("MuonIsoR03hadEt","Sum of the hadronic transverse energy in a cone of 0.3 around the muon",200,0,20);
  histocontainer_["MuonRelIso"] 			    =subDir.make<TH1D>("MuonRelIso","Relative isolation the muon",100,0,1);
  histocontainer_["MuonVetoEm"] 			    =subDir.make<TH1D>("MuonVetoEm","Veto electromagnetic energy deposit in a cone of 0.07",200,0,20);
  histocontainer_["MuonVetoHad"]			    =subDir.make<TH1D>("MuonVetoHad","Veto hadronic energy deposit in a cone of 0.1",200,0,20);
  histocontainer_["MuonCharge"]			   	    =subDir.make<TH1D>("MuonCharge","Charge of the muon",4,-2,2);
  
}
// ------------ method called once each job just after ending the event loop  ------------
void MuonChecker::endJob() {
}
