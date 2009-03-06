// -*- C++ -*-
//
// Package:    DummyChecker
// Class:      DummyChecker
// 
/**\class DummyChecker DummyChecker.cc UserCode/DummyChecker/src/DummyChecker.cc

 Description: <one line class summary>

 Implementation:
     <Notes on implementation>
*/
//
// Original Author:  local user
//         Created:  Wed Feb 18 16:39:03 CET 2009
// $Id: KinematicsChecker.cc,v 1.2 2009/03/06 14:00:44 jmmaes Exp $
//
//



#include <memory>
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "PhysicsTools/UtilAlgos/interface/TFileService.h"
#include "AnalysisDataFormats/TopObjects/interface/TtGenEvent.h"
#include "TopQuarkAnalysis/TopTools/interface/JetPartonMatching.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/MET.h"

#include "TDirectory.h"
#include "TH1D.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TLorentzVector.h"

//
// class decleration
//

class KinematicsChecker : public edm::EDAnalyzer {
   public:
      explicit KinematicsChecker(const edm::ParameterSet&);
      ~KinematicsChecker();


   private:
      virtual void beginJob(const edm::EventSetup&) ;
      virtual void analyze(const edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;

      edm::InputTag jets_;
      edm::InputTag muons_;
      edm::InputTag mets_;
  int matchingAlgo_;
      bool useMaxDist_;
      bool useDeltaR_;
      double maxDist_;

      //Histograms are booked in the beginJob() method
      std::map<std::string,TDirectory*> TDirectorycontainer_; // simple map to contain all TDirectory.
      //std::map<std::string,TH1D*> TH1Dcontainer_; // simple map to contain all TH1D.
      std::map<std::string,TH1D*> TH1Dcontainer_[4]; // simple map to contain all TH1D.
      std::map<std::string,TH1F*> TH1Fcontainer_; // simple map to contain all TH1F.
      std::map<std::string,TH2F*> TH2Fcontainer_; // simple map to contain all TH2F.

  std::vector< double > jetsAcceptance_, muonsAcceptance_;
  std::vector< double > nJetsAcceptance, nMuonsAcceptance;
  
   
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
KinematicsChecker::KinematicsChecker(const edm::ParameterSet& iConfig)

{
   //now do what ever initialization is needed
  jets_             =   iConfig.getParameter<edm::InputTag>( "jetsName" );
  muons_            =   iConfig.getParameter<edm::InputTag>( "muonsName" );
  mets_             =   iConfig.getParameter<edm::InputTag>( "metsName" );
  matchingAlgo_     =   iConfig.getParameter<int>( "matchingAlgo" );
  useMaxDist_       =   iConfig.getParameter<bool>( "useMaxDist" );
  useDeltaR_        =   iConfig.getParameter<bool>( "useDeltaR" );
  maxDist_          =   iConfig.getParameter<double>( "maxDist" );
  jetsAcceptance_   =  iConfig.getParameter< std::vector< double > >( "jetsAcceptance" );
  muonsAcceptance_  =  iConfig.getParameter< std::vector< double > >( "muonsAcceptance" );
}


KinematicsChecker::~KinematicsChecker()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called to for each event  ------------
void
KinematicsChecker::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;
   using namespace std;

#ifdef THIS_IS_AN_EVENT_EXAMPLE
   Handle<ExampleData> pIn;
   iEvent.getByLabel("example",pIn);
#endif
   
#ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
   ESHandle<SetupData> pSetup;
   iSetup.get<SetupRecord>().get(pSetup);
#endif
  
   Handle< std::vector<pat::Jet> > jets;
   iEvent.getByLabel(jets_, jets);
      
   Handle< std::vector<pat::Muon> > muons;
   iEvent.getByLabel(muons_, muons);

   Handle< std::vector<pat::MET> > mets;
   iEvent.getByLabel(mets_, mets);
   
   //Check if branches are available
   if (!jets.isValid()){
     edm::LogWarning  ("NoJetsFound") << "My warning message - NoJetsFound"; 
     throw cms::Exception("ProductNotFound") <<"jet collection not found"<<std::endl;
   }
   if (!muons.isValid()){
     edm::LogWarning  ("NoMuonsFound") << "My warning message - NoMuonsFound"; 
     throw cms::Exception("ProductNotFound") <<"muon collection not found"<<std::endl;
   }
   if (!mets.isValid()){
     edm::LogWarning  ("NoMetsFound") << "My warning message - NoMetsFound";
     throw cms::Exception("ProductNotFound") <<"MET collection not found"<<std::endl;
   }

   //Create TLorentzVector objets from PAT objects
   vector< reco::Particle::LorentzVector > ObjectP4s[3];
   for( unsigned int i=0;i<jets->size();i++) {ObjectP4s[0].push_back((*jets)[i].p4());}
   for( unsigned int i=0;i<muons->size();i++) {ObjectP4s[1].push_back((*muons)[i].p4());}
   for( unsigned int i=0;i<mets->size();i++) {ObjectP4s[2].push_back((*mets)[i].p4());}

   //Fill Histograms
   for(unsigned int i=0;i<3;i++){
     TH1Dcontainer_[i]["number"]->Fill(ObjectP4s[i].size()); 
     for(unsigned int j=0; j<ObjectP4s[i].size(); j++) TH1Dcontainer_[i]["pt"]->Fill(ObjectP4s[i][j].Pt());  
     for(unsigned int j=0;j<ObjectP4s[i].size();j++) TH1Dcontainer_[i]["eta"]->Fill(ObjectP4s[i][j].Eta()); 
     for(unsigned int j=0;j<ObjectP4s[i].size();j++) TH1Dcontainer_[i]["et"]->Fill(ObjectP4s[i][j].Et());
     for(unsigned int j=0;j<ObjectP4s[i].size();j++) TH1Dcontainer_[i]["energy"]->Fill(ObjectP4s[i][j].E());
     for(unsigned int j=0;j<ObjectP4s[i].size();j++) TH1Dcontainer_[i]["theta"]->Fill(ObjectP4s[i][j].Theta());
     for(unsigned int j=0;j<ObjectP4s[i].size();j++) TH1Dcontainer_[i]["phi"]->Fill(ObjectP4s[i][j].Phi());
   }
   
   //Check acceptance 
   bool bj0=false;
   bool bj1=false;
   bool bm0=false;
   bool bm1=false;
   
   for(unsigned int j=0; j<ObjectP4s[0].size(); j++){
     if(ObjectP4s[0][j].Pt()<jetsAcceptance_[1]) bj1=true;
     if(fabs(ObjectP4s[0][j].Eta())>jetsAcceptance_[0]) bj0=true;
   }
   if(bj0) nJetsAcceptance[0]++;
   if(bj1) nJetsAcceptance[1]++;

   for(unsigned int j=0; j<ObjectP4s[1].size(); j++){
     if(ObjectP4s[1][j].Pt()<muonsAcceptance_[1]) bm1=true;
     if(fabs(ObjectP4s[1][j].Eta())>muonsAcceptance_[0]) bm0=true;
   }
   if(bm0) nMuonsAcceptance[0]++;
   if(bm1) nMuonsAcceptance[1]++;
   

   
   for(unsigned int i=0;i<3;i++){
     ObjectP4s[i].clear();
   }
   
   //matching
   Handle<TtGenEvent> genEvt;
   iEvent.getByLabel ("genEvt",genEvt);
  
   //Check if branch is available  
   if (!genEvt.isValid()){
     edm::LogWarning  ("NoGenEvtFound") << "My warning message - NoGenEvtFound";
   }
   
   if (genEvt.isValid()){
     
     if(genEvt->isSemiLeptonic(genEvt->kMuon)) {
       
       // Matching index : Hadronic Q  = 0, Hadronic Q' = 1, Hadronic b  = 2, Leptonic b  = 3;
       std::vector<const reco::Candidate *> quarks;
       quarks.push_back(genEvt->hadronicDecayQuark());
       quarks.push_back(genEvt->hadronicDecayQuarkBar());
       quarks.push_back(genEvt->hadronicDecayB());
       quarks.push_back(genEvt->leptonicDecayB());
       if(jets->size()>=4) { 

	 JetPartonMatching *GenMatch = new JetPartonMatching(quarks, *jets, matchingAlgo_, useMaxDist_, useDeltaR_, maxDist_);

	 for(unsigned int i=0; i<4; i++){
	   Int_t Idx = GenMatch->getMatchForParton(i,0);
	   if(Idx>=0) ObjectP4s[0].push_back(((*jets)[Idx]).p4());
	 }
	 
	 TH1Dcontainer_[3]["number"]->Fill(ObjectP4s[0].size()); 
	 for(unsigned int j=0; j<ObjectP4s[0].size(); j++) TH1Dcontainer_[3]["pt"]->Fill(ObjectP4s[0][j].Pt());  
	 for(unsigned int j=0;j<ObjectP4s[0].size();j++) TH1Dcontainer_[3]["eta"]->Fill(ObjectP4s[0][j].Eta()); 
	 for(unsigned int j=0;j<ObjectP4s[0].size();j++) TH1Dcontainer_[3]["et"]->Fill(ObjectP4s[0][j].Et());
	 for(unsigned int j=0;j<ObjectP4s[0].size();j++) TH1Dcontainer_[3]["energy"]->Fill(ObjectP4s[0][j].E());
	 for(unsigned int j=0;j<ObjectP4s[0].size();j++) TH1Dcontainer_[3]["theta"]->Fill(ObjectP4s[0][j].Theta());
	 for(unsigned int j=0;j<ObjectP4s[0].size();j++) TH1Dcontainer_[3]["phi"]->Fill(ObjectP4s[0][j].Phi());

       }
     }
   }

 

    //tree levels of message
   //no endl needed 
   //use '\n' to go to next line
   //we can use different category for the same EDAnalyser
   //ex: NoDataFound - LinkBroken - TooMuchDataFound - SummaryError - MainResults
   
   //edm::LogError  ("category") << "My error message";    // or  edm::LogProblem  (not formated)
   //edm::LogWarning  ("category") << "My warning message"; // or  edm::LogPrint    (not formated)
   //edm::LogInfo   ("category") << "My LogInfo message";  // or  edm::LogVerbatim (not formated)

   //use Warning for event by event problem 
   //edm::LogWarning  ("NoDataFound") << "My warning message - NoDataFound"; // or  edm::LogPrint    (not formated)
   //edm::LogWarning  ("LinkBroken") << "My warning message - LinkBroken"; // or  edm::LogPrint    (not formated)

}


// ------------ method called once each job just before starting event loop  ------------
void 
KinematicsChecker::beginJob(const edm::EventSetup&)
{
  edm::Service<TFileService> fs;
  if (!fs) throw edm::Exception(edm::errors::Configuration, "TFileService missing from configuration!");
  
  std::string ObjectNames[4] = {"Jets","Muons","METs","MatchedJets"};

  std::vector< TFileDirectory > subDirs;
  for(unsigned int i=0;i<4;i++) subDirs.push_back(fs->mkdir( ObjectNames[i] ));

  for(unsigned int i=0;i<4;i++){
    TH1Dcontainer_[i]["number"] = subDirs[i].make<TH1D>("number" ,"number of objects",50,0,50);
    TH1Dcontainer_[i]["pt"] = subDirs[i].make<TH1D>("pt" ,"pt",50,0,100);
    TH1Dcontainer_[i]["et"] = subDirs[i].make<TH1D>("et" ,"et",50,0,100);
    TH1Dcontainer_[i]["eta"] = subDirs[i].make<TH1D>("eta" ,"eta",50,-6,6);
    TH1Dcontainer_[i]["energy"] = subDirs[i].make<TH1D>("energy" ,"energy",50,0,200);
    TH1Dcontainer_[i]["theta"] = subDirs[i].make<TH1D>("theta" ,"theta",50,0,3.2);
    TH1Dcontainer_[i]["phi"] = subDirs[i].make<TH1D>("phi" ,"phi",50,-3.2,3.2);
  }
 
  nJetsAcceptance.push_back(0);
  nJetsAcceptance.push_back(0);
  nMuonsAcceptance.push_back(0);
  nMuonsAcceptance.push_back(0);

}

// ------------ method called once each job just after ending the event loop  ------------
void 
KinematicsChecker::endJob() {
   //use LogError to summarise the error that happen in the execution (by example from warning) (ex: Nof where we cannot access such variable)
   //edm::LogError  ("SummaryError") << "My error message \n";    // or  edm::LogProblem  (not formated)
   //use LogInfo to summarise information (ex: pourcentage of events matched ...)
   //edm::LogInfo   ("MainResults") << "My LogInfo message \n";  // or  edm::LogVerbatim (not formated)
  /* edm::LogInfo   ("MainResults") << "Number of events where at least one jet is outside |eta|<" << jetsAcceptance_[0] << ": " << nJetsAcceptance[0] << "\n";
  edm::LogInfo   ("MainResults") << "Number of events where at least one jet has pt<" << jetsAcceptance_[1] << ": " << nJetsAcceptance[1] << "\n";
  edm::LogInfo   ("MainResults") << "Number of events where at least one muon is outside |eta|<" << muonsAcceptance_[0] << ": " << nMuonsAcceptance[0] << "\n";
  edm::LogInfo   ("MainResults") << "Number of events where at least one muon has pt<" << muonsAcceptance_[1] << ": " << nMuonsAcceptance[1] << "\n";*/
  edm::LogError   ("MainError") << "Number of events where at least one jet is outside |eta|<" << jetsAcceptance_[0] << ": " << nJetsAcceptance[0] << "\n";
  edm::LogError   ("MainError") << "Number of events where at least one jet has pt<" << jetsAcceptance_[1] << ": " << nJetsAcceptance[1] << "\n";
  edm::LogError   ("MainError") << "Number of events where at least one muon is outside |eta|<" << muonsAcceptance_[0] << ": " << nMuonsAcceptance[0] << "\n";
  edm::LogError   ("MainError") << "Number of events where at least one muon has pt<" << muonsAcceptance_[1] << ": " << nMuonsAcceptance[1] << "\n";
}

//define this as a plug-in
DEFINE_FWK_MODULE(KinematicsChecker);
