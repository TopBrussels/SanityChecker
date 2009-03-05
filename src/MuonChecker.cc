 // -*- C++ -*-
//
// Package:    MuonChecker
// Class:      MuonChecker
// 
/**\class MuonChecker MuonChecker.cc UserCode/TopBrussels/MuonChecker/src/MuonChecker.cc

 Description: performs simple sanity checks for the muons

 Implementation:
 
*/
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
//needed for TFileService
#include "PhysicsTools/UtilAlgos/interface/TFileService.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
//needed for MessageLogger
#include "FWCore/MessageLogger/interface/MessageLogger.h"


#include "AnalysisDataFormats/TopObjects/interface/TtGenEvent.h"

#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/Candidate/interface/OverlapChecker.h"

#include "TDirectory.h"
#include "TH1D.h"

//
// class decleration
//

class MuonChecker : public edm::EDAnalyzer {
   public:
      explicit MuonChecker(const edm::ParameterSet&);
      ~MuonChecker();


   private:
      virtual void beginJob(const edm::EventSetup&) ;
      virtual void analyze(const edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;

      // ----------member data ---------------------------

     std::map<std::string,TH1D*> histocontainer_; // simple map to contain all histograms. Histograms are booked in the beginJob() method
     std::map<int,std::string> pdgid_;
  
     edm::InputTag muoLabel_;
     bool verbose_;
  
     int NbOfEvents;
     int NbOfSemiMuEvents;

};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
MuonChecker::MuonChecker(const edm::ParameterSet& iConfig)
{
	muoLabel_ = iConfig.getParameter<edm::InputTag>("muonTag");
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

   // PAT Muons
   edm::Handle<std::vector<pat::Muon> > muonHandle;
   iEvent.getByLabel(muoLabel_,muonHandle);

   //Check if branch is available  
   if (!muonHandle.isValid())
   {
     edm::LogWarning  ("NoMuonsFound") << "My warning message - NoMuonsFound";
     throw cms::Exception("ProductNotFound") <<"MUON collection not found"<<std::endl;
   }
   const std::vector<pat::Muon> & muons = *muonHandle;

   edm::Handle<TtGenEvent> genEvtHandle;
   iEvent.getByLabel ("genEvt",genEvtHandle);

   //Check if branch is available  
   if (!genEvtHandle.isValid())
   {
     edm::LogWarning  ("NoGenEvtFound") << "My warning message - NoGenEvtFound";
   }
   const TtGenEvent & genEvt = *genEvtHandle;

   bool GenLeplMatch = false;
   OverlapChecker MyChecker;
   std::string MatchedMuon = "";

   int NbOfMuonId = 13;
   char **MuonID    = new char* [NbOfMuonId];
   MuonID[0] = new char[100]; MuonID[0] = "All";                              // dummy options - always true
   MuonID[1] = new char[100]; MuonID[1] = "AllGlobalMuons";                   // checks isGlobalMuon flag
   MuonID[2] = new char[100]; MuonID[2] = "AllStandAloneMuons"; 	      // checks isStandAloneMuon flag
   MuonID[3] = new char[100]; MuonID[3] = "AllTrackerMuons";		      // checks isTrackerMuon flag
   MuonID[4] = new char[100]; MuonID[4] = "TrackerMuonArbitrated";	      // resolve ambiguity of sharing segments
   MuonID[5] = new char[100]; MuonID[5] = "AllArbitrated";		      // all muons with the tracker muon arbitrated
   MuonID[6] = new char[100]; MuonID[6] = "GlobalMuonPromptTight";	      // global muons with tighter fit requirements
   MuonID[7] = new char[100]; MuonID[7] = "TMLastStationLoose"; 	      // penetration depth loose selector
   MuonID[8] = new char[100]; MuonID[8] = "TMLastStationTight"; 	      // penetration depth tight selector
   MuonID[9] = new char[100]; MuonID[9] = "TM2DCompatibilityLoose";	      // likelihood based loose selector
   MuonID[10]= new char[100]; MuonID[10]= "TM2DCompatibilityTight";	      // likelihood based tight selector
   MuonID[11]= new char[100]; MuonID[11]= "TMLastStationOptimizedLowPtLoose"; // combination of TMLastStation and TMOneStation
   MuonID[12]= new char[100]; MuonID[12]= "TMLastStationOptimizedLowPtTight";  // combination of TMLastStation and TMOneStation

   for(std::vector<pat::Muon>::const_iterator muon_iter = muons.begin(); muon_iter!=muons.end(); ++muon_iter)
   {
	//if(muon_iter->genLepton() != 0 && muon_iter->genLepton()->numberOfMothers() != 0) histocontainer_["GenMuonMotherPid"]->Fill(muon_iter->genLepton()->mother(0)->pdgId());
	if(muon_iter->genLepton() == 0 || !genEvtHandle.isValid()) GenLeplMatch = false;
	else if( genEvt.isSemiLeptonic(genEvt.kMuon) )             GenLeplMatch	= MyChecker(*muon_iter->genLepton(),*genEvt.singleLepton());

	(GenLeplMatch ? MatchedMuon = "TopMuonMatch_" : MatchedMuon = "");

	for(Int_t i=0;i<NbOfMuonId;i++)
	{
		if(muon_iter->isGood(reco::Muon::SelectionType(i))) histocontainer_["MuonId"]->Fill(MuonID[i],1);
	}

	if(muon_iter->isGlobalMuon())
	{
		histocontainer_[MatchedMuon+"GlobalMuonNbOfValidHits"]  ->Fill(muon_iter->innerTrack()->numberOfValidHits());
		histocontainer_[MatchedMuon+"GlobalMuonNbOfLostHits"]   ->Fill(muon_iter->innerTrack()->numberOfLostHits());
		histocontainer_[MatchedMuon+"GlobalMuonSiliconFitPt"]	->Fill(muon_iter->innerTrack()->pt());

		histocontainer_[MatchedMuon+"GlobalMuonNbOfValidHits"]  ->Fill(muon_iter->outerTrack()->numberOfValidHits());
		histocontainer_[MatchedMuon+"GlobalMuonNbOfLostHits"]   ->Fill(muon_iter->outerTrack()->numberOfLostHits());
		histocontainer_[MatchedMuon+"GlobalMuonSiliconFitPt"]	->Fill(muon_iter->outerTrack()->pt());

		histocontainer_[MatchedMuon+"GlobalMuonGlobalTrackPt"]	->Fill(muon_iter->globalTrack()->pt());
		histocontainer_[MatchedMuon+"GlobalMuonGlobalTrackD0"]	->Fill(muon_iter->globalTrack()->d0());
		histocontainer_[MatchedMuon+"GlobalMuonGlobalTrackChi2"]->Fill(muon_iter->globalTrack()->chi2());
	}
	
	histocontainer_[MatchedMuon+"MuonCaloCompatibility"] ->Fill(muon_iter->caloCompatibility());
	histocontainer_[MatchedMuon+"MuonSegmCompatibility"] ->Fill(muon_iter->segmentCompatibility());

	histocontainer_[MatchedMuon+"MuonIsoR03SumPt"] ->Fill(muon_iter->isolationR03().sumPt);
	histocontainer_[MatchedMuon+"MuonIsoR03emEt"]  ->Fill(muon_iter->isolationR03().emEt);
	histocontainer_[MatchedMuon+"MuonIsoR03hadEt"] ->Fill(muon_iter->isolationR03().hadEt);
  	histocontainer_[MatchedMuon+"MuonRelIso"]      ->Fill(muon_iter->pt()/(muon_iter->pt()+muon_iter->isolationR03().sumPt + muon_iter->isolationR03().emEt + muon_iter->isolationR03().hadEt)); 
  
  	//get the veto cone information
  	histocontainer_[MatchedMuon+"MuonVetoEm"]  ->Fill(muon_iter->ecalIsoDeposit()->candEnergy());
  	histocontainer_[MatchedMuon+"MuonVetoHad"] ->Fill(muon_iter->hcalIsoDeposit()->candEnergy());

   }

   edm::LogWarning  ("NoDataFound") << "My warning message - NoDataFound"; // or  edm::LogPrint    (not formated)
   edm::LogWarning  ("LinkBroken") << "My warning message - LinkBroken"; // or  edm::LogPrint    (not formated)
}


// ------------ method called once each job just before starting event loop  ------------
void 
MuonChecker::beginJob(const edm::EventSetup&)
{
  edm::Service<TFileService> fs;
  if (!fs) throw edm::Exception(edm::errors::Configuration, "TFileService missing from configuration!");

  histocontainer_["GenMuonMotherPid"]      		    =fs->make<TH1D>("GenMuonMotherPid","Mother PDG Id of generated muon matching the leading reconstructed muon",1000,0,1000);
  histocontainer_["GlobalMuonSiliconFitPt"]	 	    =fs->make<TH1D>("GlobalMuonSiliconFitPt","Silicon fit transverse momentum (global muons)",400,0,400);
  histocontainer_["GlobalMuonGlobalTrackPt"]        	    =fs->make<TH1D>("GlobalMuonPt","global muons transverse momentum",400,0,400);
  histocontainer_["GlobalMuonGlobalTrackChi2"]              =fs->make<TH1D>("GlobalMuonGlobalTrackChi2","global muons Chi^2 ",400,0,400);
  histocontainer_["GlobalMuonGlobalTrackD0"]                =fs->make<TH1D>("GlobalMuonGlobalTrackD0","global muon impact parameter (with p_{T}>20GeV/c,|#eta|<2.4)",400,-0.4,0.4);
  histocontainer_["GlobalMuonNbOfValidHits"]		    =fs->make<TH1D>("GlobalMuonNbOfValidHits","Number of valid hits in silicon fit for track muons",400,0,200);
  histocontainer_["GlobalMuonNbOfLostHits"]		    =fs->make<TH1D>("GlobalMuonNbOfLostHits","Number of lost hits in silicon fit for track muons",400,0,200);
  histocontainer_["MuonCaloCompatibility"]                  =fs->make<TH1D>("CaloCompatibility","Value of the LR measuring the probability that the muon is calo-compatible",100,0,1);
  histocontainer_["MuonSegmCompatibility"]                  =fs->make<TH1D>("SegmCompatibility","Value of the LR measuring the probability that the muon is segment-compatible",100,0,1);
  histocontainer_["MuonIsoR03SumPt"]   		            =fs->make<TH1D>("MuonIsoR03SumPt","Sum of the track transverse momenta in a cone of 0.3 around the muon",200,0,20);
  histocontainer_["MuonIsoR03emEt"]			    =fs->make<TH1D>("MuonIsoR03emEt","Sum of the electromagnetic transverse energy in a cone of 0.3 around the muon",200,0,20);
  histocontainer_["MuonIsoR03hadEt"]		            =fs->make<TH1D>("MuonIsoR03hadEt","Sum of the hadronic transverse energy in a cone of 0.3 around the muon",200,0,20);
  histocontainer_["MuonRelIso"] 			    =fs->make<TH1D>("MuonRelIso","Relative isolation the muon",100,0,1);
  histocontainer_["MuonVetoEm"] 			    =fs->make<TH1D>("MuonVetoEm","Veto electromagnetic energy deposit in a cone of 0.07",200,0,20);
  histocontainer_["MuonVetoHad"]			    =fs->make<TH1D>("MuonVetoHad","Veto hadronic energy deposit in a cone of 0.1",200,0,20);
  histocontainer_["TopMuonMatch_SiliconFitPt"]		    =fs->make<TH1D>("TopMuonMatch_SiliconFitPt","silicon fit transverse momentum (global muons) matching the muon from the top decay",400,0,400);
  histocontainer_["TopMuonMatch_GlobalMuonSiliconFitPt"]    =fs->make<TH1D>("TopMuonMatch_GlobalMuonSiliconFitPt","Silicon fit transverse momentum (global muons) matching the muon from the top decay",400,0,400);
  histocontainer_["TopMuonMatch_GlobalMuonGlobalTrackPt"]   =fs->make<TH1D>("TopMuonMatch_GlobalMuonGlobalTrackPt","muon global track pt (muon matching the gen muon from top decay, with p_{T}>20GeV/c,|#eta|<2.4)",400,0,400);
  histocontainer_["TopMuonMatch_GlobalMuonGlobalTrackChi2"] =fs->make<TH1D>("TopMuonMatch_GlobalMuonGlobalTrackChi2","muon global track #chi^{2} (muon matching the gen muon from top decay, with p_{T}>20GeV/c,|#eta|<2.4)",400,0,400);
  histocontainer_["TopMuonMatch_GlobalMuonGlobalTrackD0"]   =fs->make<TH1D>("TopMuonMatch_GlobalMuonGlobalTrackD0","muon global track impact parameter (muon matching the gen muon from top decay, with p_{T}>20GeV/c,|#eta|<2.4)",400,-0.4,0.4);
  histocontainer_["TopMuonMatch_GlobalMuonNbOfValidHits"]   =fs->make<TH1D>("TopMuonMatch_GlobalMuonNbOfValidHits","Number of valid hits in silicon fit for global muons (muon matching the gen muon from top decay",400,0,200);
  histocontainer_["TopMuonMatch_GlobalMuonNbOfLostHits"]    =fs->make<TH1D>("TopMuonMatch_GlobalMuonNbOfLostHits","Number of lost hits in silicon fit for global muons (muon matching the gen muon from top decay",400,0,200);
  histocontainer_["TopMuonMatch_MuonCaloCompatibility"]     =fs->make<TH1D>("TopMuonMatch_MuonCaloCompatibility","Value of the LR measuring the probability that the muon is calo-compatible with a MIP (muon matching the gen muon from top decay), ",100,0,1);
  histocontainer_["TopMuonMatch_MuonSegmCompatibility"]     =fs->make<TH1D>("TopMuonMatch_MuonSegmCompatibility","Value of the LR measuring the probability that the muon is segment-compatible with a MIP (muon matching the gen muon from top decay), ",100,0,1);
  histocontainer_["TopMuonMatch_MuonIsoR03SumPt"]	    =fs->make<TH1D>("TopMuonMatch_MuonIsoR03SumPt","Sum of the track transverse momenta in a cone of 0.3 around the muon",200,0,20);
  histocontainer_["TopMuonMatch_MuonIsoR03emEt"]	    =fs->make<TH1D>("TopMuonMatch_MuonIsoR03emEt","Sum of the electromagnetic transverse energy in a cone of 0.3 around the muon",200,0,20);
  histocontainer_["TopMuonMatch_MuonIsoR03hadEt"]	    =fs->make<TH1D>("TopMuonMatch_MuonIsoR03hadEt","Sum of the hadronic transverse energy in a cone of 0.3 around the muon",200,0,20);
  histocontainer_["TopMuonMatch_MuonRelIso"]		    =fs->make<TH1D>("TopMuonMatch_MuonRelIso","Relative isolation the muon",100,0,1);
  histocontainer_["TopMuonMatch_MuonVetoEm"]		    =fs->make<TH1D>("TopMuonMatch_MuonVetoEm","Veto electromagnetic energy deposit in a cone of 0.07",200,0,20);
  histocontainer_["TopMuonMatch_MuonVetoHad"]		    =fs->make<TH1D>("TopMuonMatch_MuonVetoHad","Veto hadronic energy deposit in a cone of 0.1",200,0,20);
}

// ------------ method called once each job just after ending the event loop  ------------
void 
MuonChecker::endJob()
{
   //use LogError to summarise the error that happen in the execution (by example from warning) (ex: Nof where we cannot access such variable)
   edm::LogError  ("SummaryError") << "My error message \n";    // or  edm::LogProblem  (not formated)
   //use LogInfo to summarise information (ex: pourcentage of events matched ...)
   edm::LogInfo   ("MainResults") << "My LogInfo message \n";  // or  edm::LogVerbatim (not formated)
}

//define this as a plug-in
DEFINE_FWK_MODULE(MuonChecker);
