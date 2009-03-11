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
  
     edm::InputTag muoLabel_;
     bool verbose_;
  
     int NbOfEvents;
     int NbOfSemiMuEvents;
     int NbOfNoGenMu;
     int NbOfMu;
     int MuonsOrigin[12];
     int Muoncharge[2];

};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

   const int NbOfMuonId = 13;
   const char *MuonID[NbOfMuonId] = {
   "All",                             // dummy options - always true
   "AllGlobalMuons",                  // checks isGlobalMuon flag
   "AllStandAloneMuons", 	      // checks isStandAloneMuon flag
   "AllTrackerMuons",		      // checks isTrackerMuon flag
   "TrackerMuonArbitrated",	      // resolve ambiguity of sharing segments
   "AllArbitrated",		      // all muons with the tracker muon arbitrated
   "GlobalMuonPromptTight",	      // global muons with tighter fit requirements
   "TMLastStationLoose", 	      // penetration depth loose selector
   "TMLastStationTight", 	      // penetration depth tight selector
   "TM2DCompatibilityLoose",	      // likelihood based loose selector
   "TM2DCompatibilityTight",	      // likelihood based tight selector
   "TMLastStationOptimizedLowPtLoose",// combination of TMLastStation and TMOneStation
   "TMLastStationOptimizedLowPtTight" // combination of TMLastStation and TMOneStation
   };

   const char *MesonsPdgId[6] = {
   "W boson",                          // PdgId = 24 or -24
   "Light Mesons (I=1)",               // PdgId/100 = 1
   "Light Mesons (I=0)",               // PdgId/100 = 2
   "Strange Mesons",                   // PdgId/100 = 3
   "Charmed Mesons",                   // PdgId/100 = 4
   "Bottom Mesons",                    // PdgId/100 = 5
   };

   const char *BaryonsPdgId[6] = {
   "Other"
   "Light Baryons",                    // PdgId/1000 = 1
   "Light Baryons",                    // PdgId/1000 = 2
   "Strange Baryons",                  // PdgId/1000 = 3
   "Charmed Baryons",                  // PdgId/1000 = 4
   "Bottom Baryons"                    // PdgId/1000 = 5
   };

//
// constructors and destructor
//
MuonChecker::MuonChecker(const edm::ParameterSet& iConfig)
{
	muoLabel_ = iConfig.getParameter<edm::InputTag>("muonTag");
	NbOfEvents       = 0;
	NbOfSemiMuEvents = 0;
	NbOfNoGenMu      = 0;
	NbOfMu           = 0;
	for(int i=0;i<12;i++) MuonsOrigin[i]=0;
	Muoncharge[0]    = 0;
	Muoncharge[1]    = 0;
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
   NbOfEvents++;

   // PAT Muons
   edm::Handle<std::vector<pat::Muon> > muonHandle;
   iEvent.getByLabel(muoLabel_,muonHandle);

   //Check if branch is available  
   if (!muonHandle.isValid())
   {
     edm::LogWarning  ("NoDataFound") << "--- NoMuonsFound ---";
     throw cms::Exception("ProductNotFound") <<"Muon collection not found"<<std::endl;
   }
   const std::vector<pat::Muon> & muons = *muonHandle;

   edm::Handle<TtGenEvent> genEvtHandle;
   iEvent.getByLabel ("genEvt",genEvtHandle);

   //Check if branch is available  
   if (!genEvtHandle.isValid())
   {
     edm::LogWarning  ("NoDataFound") << "--- NoGenEvtFound ---";
   }
   const TtGenEvent & genEvt = *genEvtHandle;

   bool GenLeplMatch = false;
   //OverlapChecker MyChecker;
   std::string MatchedMuon = "";

   for(std::vector<pat::Muon>::const_iterator muon_iter = muons.begin(); muon_iter!=muons.end(); ++muon_iter)
   {
	NbOfMu++;
	GenLeplMatch = false;
	
	if(muon_iter->genParticle() == 0)
	{
		edm::LogWarning ("LinkBroken") << "--- NoGenMuonFound ---";
		NbOfNoGenMu++;
	}
	else if( muon_iter->genParticle()->numberOfMothers() != 0 )
	{
		Int_t PDG = fabs(muon_iter->genParticle()->mother(0)->pdgId());
		if(PDG<100)
		{
			if(PDG == 24) {histocontainer_["GenMuonMotherPid"]->Fill(MesonsPdgId[0],1);MuonsOrigin[0]++;}
			else if(muon_iter->genParticle()->mother(0)->numberOfMothers() != 0 )
			{
				if(fabs(muon_iter->genParticle()->mother(0)->mother(0)->pdgId()) == 24) 
				{
					histocontainer_["GenMuonMotherPid"]->Fill(MesonsPdgId[0],1);MuonsOrigin[0]++;
					if(genEvtHandle.isValid())//->isNonnull())
					{
						if( genEvt.isSemiLeptonic(genEvt.kMuon) && muon_iter->genParticle()->mother(0)->mother(0)->numberOfMothers() != 0 )
						{
							NbOfSemiMuEvents++;
							if(fabs(muon_iter->genParticle()->mother(0)->mother(0)->mother(0)->pdgId()) == 6) GenLeplMatch = true;
							//std::cout<<"Check overlap with the genparticle : "<<GenLeplMatch<<std::endl;
						}
					}
				}
			}
		}
		else if(100<PDG && PDG<600)
		{
			histocontainer_["GenMuonMotherPid"]->Fill(MesonsPdgId[PDG/100],1);
			MuonsOrigin[PDG/100]++;
		}
		else if(999<PDG && PDG<6000)
		{
			histocontainer_["GenMuonMotherPid"]->Fill(BaryonsPdgId[PDG/1000],1);
			MuonsOrigin[(PDG/1000)+6]++;
		}
		else
		{
			histocontainer_["GenMuonMotherPid"]->Fill(BaryonsPdgId[0],1);
			MuonsOrigin[6]++;
		}
	}

	(GenLeplMatch ? MatchedMuon = "TopMuonMatch_" : MatchedMuon = "");

	for(Int_t i=0;i<NbOfMuonId;i++)
	{
		if(muon_iter->isGood(reco::Muon::SelectionType(i))) histocontainer_["MuonId"]->Fill(MuonID[i],1);
	}

	if(muon_iter->isGlobalMuon())
	{
		histocontainer_[MatchedMuon+"GlobalMuonInnerTrackNbOfValidHits"]  ->Fill(muon_iter->innerTrack()->numberOfValidHits());
		histocontainer_[MatchedMuon+"GlobalMuonInnerTrackNbOfLostHits"]   ->Fill(muon_iter->innerTrack()->numberOfLostHits());
		histocontainer_[MatchedMuon+"GlobalMuonInnerTrackPt"]   	  ->Fill(muon_iter->innerTrack()->pt());

		histocontainer_[MatchedMuon+"GlobalMuonOuterTrackNbOfValidHits"]  ->Fill(muon_iter->outerTrack()->numberOfValidHits());
		histocontainer_[MatchedMuon+"GlobalMuonOuterTrackNbOfLostHits"]   ->Fill(muon_iter->outerTrack()->numberOfLostHits());
		histocontainer_[MatchedMuon+"GlobalMuonOuterTrackPt"]    	  ->Fill(muon_iter->outerTrack()->pt());

		histocontainer_[MatchedMuon+"GlobalMuonGlobalTrackPt"]		  ->Fill(muon_iter->globalTrack()->pt());
		histocontainer_[MatchedMuon+"GlobalMuonGlobalTrackD0"]		  ->Fill(muon_iter->globalTrack()->d0());
		histocontainer_[MatchedMuon+"GlobalMuonGlobalTrackChi2"]	  ->Fill(muon_iter->globalTrack()->normalizedChi2());
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

	histocontainer_[MatchedMuon+"MuonCharge"] ->Fill(muon_iter->charge());
	(muon_iter->charge()>0 ? Muoncharge[0]++ : Muoncharge[1]++);
   }
}


// ------------ method called once each job just before starting event loop  ------------
void 
MuonChecker::beginJob(const edm::EventSetup&)
{
  edm::Service<TFileService> fs;
  if (!fs) throw edm::Exception(edm::errors::Configuration, "TFileService missing from configuration!");

  histocontainer_["GenMuonMotherPid"]      		    =fs->make<TH1D>("GenMuonMotherPid","Mother PDG Id of generated muon matching the leading reconstructed muon",12,0,12);histocontainer_["GenMuonMotherPid"]->LabelsDeflate("X");
  histocontainer_["MuonId"]      		            =fs->make<TH1D>("MuonId","Muon identification algorithms",13,0,13);histocontainer_["MuonId"]->LabelsDeflate("X");
  ///////////////////////////////////////////////////////////////////////////////
  // for all global muons (except the one from the semi-muonic top quark decay) :

  histocontainer_["GlobalMuonInnerTrackNbOfValidHits"]      =fs->make<TH1D>("GlobalMuonInnerTrackNbOfValidHits","Number of valid hits in silicon fit for global muons",400,0,20);
  histocontainer_["GlobalMuonInnerTrackNbOfLostHits"]	    =fs->make<TH1D>("GlobalMuonInnerTrackNbOfLostHits","Number of lost hits in silicon fit for global muons",400,0,20);
  histocontainer_["GlobalMuonInnerTrackPt"]	            =fs->make<TH1D>("GlobalMuonInnerTrackPt","Inner track transverse momentum for global muons",400,0,400);

  histocontainer_["GlobalMuonOuterTrackNbOfValidHits"]      =fs->make<TH1D>("GlobalMuonOuterTrackNbOfValidHits","Number of valid hits in silicon fit for global muons",400,0,20);
  histocontainer_["GlobalMuonOuterTrackNbOfLostHits"]       =fs->make<TH1D>("GlobalMuonOuterTrackNbOfLostHits","Number of lost hits in silicon fit for global muons",400,0,20);
  histocontainer_["GlobalMuonOuterTrackPt"]		    =fs->make<TH1D>("GlobalMuonOuterTrackPt","Outer track transverse momentum for global muons",400,0,400);

  histocontainer_["GlobalMuonGlobalTrackPt"]                =fs->make<TH1D>("GlobalMuonGlobalTrackPt","Global track transverse momentum for global muons",400,0,400);
  histocontainer_["GlobalMuonGlobalTrackD0"]    	    =fs->make<TH1D>("GlobalMuonGlobalTrackD0","Global track impact parameter for global muons",400,-0.4,0.4);
  histocontainer_["GlobalMuonGlobalTrackChi2"] 		    =fs->make<TH1D>("GlobalMuonGlobalTrackChi2","Global track normalized #chi^{2} ",400,0,20);
  
  histocontainer_["MuonCaloCompatibility"]                  =fs->make<TH1D>("CaloCompatibility","Value of the LR measuring the probability that the muon is calo-compatible",100,0,1);
  histocontainer_["MuonSegmCompatibility"]                  =fs->make<TH1D>("SegmCompatibility","Value of the LR measuring the probability that the muon is segment-compatible",100,0,1);
  histocontainer_["MuonIsoR03SumPt"]   		            =fs->make<TH1D>("MuonIsoR03SumPt","Sum of the track transverse momenta in a cone of 0.3 around the muon",200,0,20);
  histocontainer_["MuonIsoR03emEt"]			    =fs->make<TH1D>("MuonIsoR03emEt","Sum of the electromagnetic transverse energy in a cone of 0.3 around the muon",200,0,20);
  histocontainer_["MuonIsoR03hadEt"]		            =fs->make<TH1D>("MuonIsoR03hadEt","Sum of the hadronic transverse energy in a cone of 0.3 around the muon",200,0,20);
  histocontainer_["MuonRelIso"] 			    =fs->make<TH1D>("MuonRelIso","Relative isolation the muon",100,0,1);
  histocontainer_["MuonVetoEm"] 			    =fs->make<TH1D>("MuonVetoEm","Veto electromagnetic energy deposit in a cone of 0.07",200,0,20);
  histocontainer_["MuonVetoHad"]			    =fs->make<TH1D>("MuonVetoHad","Veto hadronic energy deposit in a cone of 0.1",200,0,20);
  histocontainer_["MuonCharge"]			   	    =fs->make<TH1D>("MuonCharge","Charge of the muon",4,-2,2);
  
  ///////////////////////////////////////////////////////////////////////////////////////
  // for the global muon that matches the muon coming from the top quark leptonic decay : 

  histocontainer_["TopMuonMatch_GlobalMuonInnerTrackNbOfValidHits"] =fs->make<TH1D>("TopMuonMatch_GlobalMuonInnerTrackNbOfValidHits","Number of valid hits in silicon fit for global muons",400,0,20);
  histocontainer_["TopMuonMatch_GlobalMuonInnerTrackNbOfLostHits"]  =fs->make<TH1D>("TopMuonMatch_GlobalMuonInnerTrackNbOfLostHits","Number of lost hits in silicon fit for global muons",400,0,20);
  histocontainer_["TopMuonMatch_GlobalMuonInnerTrackPt"]	    =fs->make<TH1D>("TopMuonMatch_GlobalMuonInnerTrackPt","Inner track transverse momentum for global muons",400,0,400);

  histocontainer_["TopMuonMatch_GlobalMuonOuterTrackNbOfValidHits"] =fs->make<TH1D>("TopMuonMatch_GlobalMuonOuterTrackNbOfValidHits","Number of valid hits in silicon fit for global muons",400,0,20);
  histocontainer_["TopMuonMatch_GlobalMuonOuterTrackNbOfLostHits"]  =fs->make<TH1D>("TopMuonMatch_GlobalMuonOuterTrackNbOfLostHits","Number of lost hits in silicon fit for global muons",400,0,20);
  histocontainer_["TopMuonMatch_GlobalMuonOuterTrackPt"]	    =fs->make<TH1D>("TopMuonMatch_GlobalMuonOuterTrackPt","Outer track transverse momentum for global muons",400,0,400);

  histocontainer_["TopMuonMatch_GlobalMuonGlobalTrackPt"]           =fs->make<TH1D>("TopMuonMatch_GlobalMuonGlobalTrackPt","Global track transverse momentum for global muons",400,0,400);
  histocontainer_["TopMuonMatch_GlobalMuonGlobalTrackD0"]    	    =fs->make<TH1D>("TopMuonMatch_GlobalMuonGlobalTrackD0","Global track impact parameter for global muons",400,-0.4,0.4);
  histocontainer_["TopMuonMatch_GlobalMuonGlobalTrackChi2"] 	    =fs->make<TH1D>("TopMuonMatch_GlobalMuonGlobalTrackChi2","Global track normalized #chi^{2} ",400,0,20);

  histocontainer_["TopMuonMatch_MuonCaloCompatibility"]    	    =fs->make<TH1D>("TopMuonMatch_MuonCaloCompatibility","Value of the LR measuring the probability that the muon is calo-compatible with a MIP (muon matching the gen muon from top decay), ",100,0,1);
  histocontainer_["TopMuonMatch_MuonSegmCompatibility"]             =fs->make<TH1D>("TopMuonMatch_MuonSegmCompatibility","Value of the LR measuring the probability that the muon is segment-compatible with a MIP (muon matching the gen muon from top decay), ",100,0,1);
  histocontainer_["TopMuonMatch_MuonIsoR03SumPt"]	            =fs->make<TH1D>("TopMuonMatch_MuonIsoR03SumPt","Sum of the track transverse momenta in a cone of 0.3 around the muon",200,0,20);
  histocontainer_["TopMuonMatch_MuonIsoR03emEt"]	            =fs->make<TH1D>("TopMuonMatch_MuonIsoR03emEt","Sum of the electromagnetic transverse energy in a cone of 0.3 around the muon",200,0,20);
  histocontainer_["TopMuonMatch_MuonIsoR03hadEt"]	            =fs->make<TH1D>("TopMuonMatch_MuonIsoR03hadEt","Sum of the hadronic transverse energy in a cone of 0.3 around the muon",200,0,20);
  histocontainer_["TopMuonMatch_MuonRelIso"]		            =fs->make<TH1D>("TopMuonMatch_MuonRelIso","Relative isolation the muon",100,0,1);
  histocontainer_["TopMuonMatch_MuonVetoEm"]		            =fs->make<TH1D>("TopMuonMatch_MuonVetoEm","Veto electromagnetic energy deposit in a cone of 0.07",200,0,20);
  histocontainer_["TopMuonMatch_MuonVetoHad"]		            =fs->make<TH1D>("TopMuonMatch_MuonVetoHad","Veto hadronic energy deposit in a cone of 0.1",200,0,20);
  histocontainer_["TopMuonMatch_MuonCharge"]		            =fs->make<TH1D>("TopMuonMatch_MuonCharge","Charge of the muon matching the gen muon from top decay",4,-2,2);
}

// ------------ method called once each job just after ending the event loop  ------------
void 
MuonChecker::endJob()
{
   //use LogError to summarise the error that happen in the execution (by example from warning) (ex: Nof where we cannot access such variable)
   //edm::LogError  ("SummaryError") << "My error message \n";    // or  edm::LogProblem  (not formated)
   //use LogInfo to summarise information (ex: pourcentage of events matched ...)
   //edm::LogInfo   ("MainResults") << "My LogInfo message \n";  // or  edm::LogVerbatim (not formated)

  edm::LogVerbatim ("MainResults") << " ------------------------------------------";
  edm::LogVerbatim ("MainResults") << " ------------------------------------------";
  edm::LogVerbatim ("MainResults") << " --   Report from MuonSanityChecker     -- ";
  edm::LogVerbatim ("MainResults") << " ------------------------------------------";
  edm::LogVerbatim ("MainResults") << " ------------------------------------------";
  edm::LogVerbatim ("MainResults") << " Nb of processed events :"<<NbOfEvents;
  edm::LogVerbatim ("MainResults") << " Nb of semi-muonic events :"<<NbOfSemiMuEvents;
  edm::LogVerbatim ("MainResults") << " Nb of NoGen muon found : "<<NbOfNoGenMu;
  edm::LogVerbatim ("MainResults") << " Nb of muons found :"<<NbOfMu;
  edm::LogVerbatim ("MainResults") << " -- average : "<<NbOfMu/NbOfSemiMuEvents;
  edm::LogVerbatim ("MainResults") << " -- Composition : ";
  edm::LogVerbatim ("MainResults") << " ---- W bosons : "<<MuonsOrigin[0];
  edm::LogVerbatim ("MainResults") << " ---- Light mesons (I=1) (Pions, rho, etc...): "<<MuonsOrigin[1];
  edm::LogVerbatim ("MainResults") << " ---- Light mesons (I=0) (ita, omega, etc...): "<<MuonsOrigin[2];
  edm::LogVerbatim ("MainResults") << " ---- s-mesons (K^{0/+-}, etc...): "<<MuonsOrigin[3];
  edm::LogVerbatim ("MainResults") << " ---- c-mesons (D^{0/+-}, etc...): "<<MuonsOrigin[4];
  edm::LogVerbatim ("MainResults") << " ---- b-mesons (B^{0/+-}, etc...): "<<MuonsOrigin[5];
  edm::LogVerbatim ("MainResults") << " ---- Light baryons (n, p, etc...): "<<MuonsOrigin[7]+MuonsOrigin[8];
  edm::LogVerbatim ("MainResults") << " ---- s-baryons (lambda, sigma, etc...): "<<MuonsOrigin[9];
  edm::LogVerbatim ("MainResults") << " ---- c-baryons (lambda_c, sigma_c, etc...): "<<MuonsOrigin[10];
  edm::LogVerbatim ("MainResults") << " ---- b-baryons (lambda_b, sigma_b, etc...): "<<MuonsOrigin[11];
  edm::LogVerbatim ("MainResults") << " ---- Others : "<<MuonsOrigin[6];
  edm::LogVerbatim ("MainResults") << " Muon charge : ";
  edm::LogVerbatim ("MainResults") << " -- (+) : "<<Muoncharge[0];
  edm::LogVerbatim ("MainResults") << " -- (-) : "<<Muoncharge[1];
  
}

//define this as a plug-in
DEFINE_FWK_MODULE(MuonChecker);
