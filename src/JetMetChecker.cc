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
// $Id: DummyChecker.cc,v 1.1 2009/02/19 11:59:06 echabert Exp $
//
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
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/MET.h"

#include "TDirectory.h"
#include "TH1D.h"
#include "TH1F.h"
#include "TH2F.h"

//
// class decleration
//

class JetMetChecker : public edm::EDAnalyzer {
   public:
      explicit JetMetChecker(const edm::ParameterSet&);
      ~JetMetChecker();


   private:
      virtual void beginJob(const edm::EventSetup&) ;
      virtual void analyze(const edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;

      // ----------member data ---------------------------
  edm::InputTag jets_;
  edm::InputTag vertex_;
  edm::InputTag mets_;
  
  //Histograms are booked in the beginJob() method
  std::map<std::string,TDirectory*> TDirectorycontainer_; // simple map to contain all TDirectory.
  std::map<std::string,TH1D*> TH1Dcontainer_; // simple map to contain all TH1D.
  std::map<std::string,TH1F*> TH1Fcontainer_; // simple map to contain all TH1F.
  std::map<std::string,TH2F*> TH2Fcontainer_; // simple map to contain all TH2F.
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
JetMetChecker::JetMetChecker(const edm::ParameterSet& iConfig)

{
   //now do what ever initialization is needed
  jets_          =   iConfig.getParameter<edm::InputTag>( "jetsName" );
  vertex_        =   iConfig.getParameter<edm::InputTag>( "vertexName" );
  mets_          =   iConfig.getParameter<edm::InputTag>( "metsName" );

}


JetMetChecker::~JetMetChecker()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called to for each event  ------------
void
JetMetChecker::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
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
   //Here you handle the collection you want to access
   
   Handle< std::vector<pat::Jet> > jets;
   iEvent.getByLabel(jets_, jets);
 
   Handle< std::vector<reco::Vertex> > vertex;
   iEvent.getByLabel(vertex_, vertex);
   
   Handle< std::vector<pat::MET> > mets;
   iEvent.getByLabel(mets_, mets);
   
   //Check if branches are available
   if (!jets.isValid()){
     edm::LogWarning  ("NoJetsFound") << "My warning message - NoJetsFound"; 
     throw cms::Exception("ProductNotFound") <<"jet collection not found"<<std::endl;
   }
   if (!vertex.isValid()){
     edm::LogWarning  ("NoVertexFound") << "My warning message - NoVertexFound"; 
     throw cms::Exception("ProductNotFound") <<"Primary vertex collection not found"<<std::endl;
   }
   if (!mets.isValid()){
     edm::LogWarning  ("NoMetsFound") << "My warning message - NoMetsFound";
     throw cms::Exception("ProductNotFound") <<"MET collection not found"<<std::endl;
   }

   vector <CaloTowerPtr> jettowers;
 
   for( unsigned int i=0;i<jets->size();i++) {
   
     TH1Dcontainer_["Jetn90"]->Fill((*jets)[i].n90());
     TH1Dcontainer_["JetTowersArea"]->Fill((*jets)[i].towersArea());

  //towers infos
     jettowers = (*jets)[i].getCaloConstituents();
     std::vector <CaloTowerPtr>::const_iterator caloiter;
     for(caloiter=jettowers.begin();caloiter!=jettowers.end();caloiter++){
       //       double caloet=(*caloiter)->et();
       TH1Dcontainer_["JetTwrEt"]->Fill((*caloiter)->et());
       TH1Dcontainer_["JetTwrPhi"]->Fill((*caloiter)->phi());
       TH1Dcontainer_["JetTwrEta"]->Fill((*caloiter)->eta());
     }

     //tracks infos
     edm::RefVector< reco::TrackCollection > tracks;
     tracks = (*jets)[i].associatedTracks();
     for(unsigned int ti=0; ti<tracks.size(); ti++){
       TH1Dcontainer_["JetTrkPt"]->Fill(tracks[ti]->pt());
     }

   }//close for jets

   //  Primary Vertex infos;
   if(vertex->size() > 0){
     for (unsigned int j=0; j< vertex->size(); j++){
       TH1Dcontainer_["JetVtrTrkSize"]->Fill((*vertex)[j].tracksSize());
       reco::Vertex::trackRef_iterator tr ;  
       for ( tr = (*vertex)[j].tracks_begin(); tr !=(*vertex)[j].tracks_end(); tr++){
	 TH1Dcontainer_["JetVrtTrkPt"]->Fill((*tr)->pt());
       }
     }//close for vertex
   }

   //   TH1Dcontainer_["njets"]->Fill(1);
   //   TH1Dcontainer_["njetsInSubDir"]->Fill(2);
   // TH1Dcontainer_["njetsInSubSubDir"]->Fill(3);

   //tree levels of message
   //no endl needed 
   //use '\n' to go to next line
   //we can use different category for the same EDAnalyser
   //ex: NoDataFound - LinkBroken - TooMuchDataFound - SummaryError - MainResults
   
   //edm::LogError  ("category") << "My error message";    // or  edm::LogProblem  (not formated)
   //edm::LogWarning  ("category") << "My warning message"; // or  edm::LogPrint    (not formated)
   //edm::LogInfo   ("category") << "My LogInfo message";  // or  edm::LogVerbatim (not formated)

   //use Warning for event by event problem 
   edm::LogWarning  ("NoDataFound") << "My warning message - NoDataFound"; // or  edm::LogPrint    (not formated)
   edm::LogWarning  ("LinkBroken") << "My warning message - LinkBroken"; // or  edm::LogPrint    (not formated)

}


// ------------ method called once each job just before starting event loop  ------------
void 
JetMetChecker::beginJob(const edm::EventSetup&)
{
  edm::Service<TFileService> fs;
  if (!fs) throw edm::Exception(edm::errors::Configuration, "TFileService missing from configuration!");

   TFileDirectory subDir = fs->mkdir( "PatJets" );
   // TFileDirectory subsubDir = subDir.mkdir( "PatJets" );

  //define the histograms booked
  //TH1D
   //  TH1Dcontainer_["JetTwrEt"] = fs->make<TH1D>("JetTwrEt" ,"jet towers Et ",100,0,1000);// wrote in main directory
  TH1Dcontainer_["Jetn90"] = subDir.make<TH1D>("Jetn90" ,"n90 ",10,0,10);
  TH1Dcontainer_["JetTowersArea"] = subDir.make<TH1D>("JetTowersArea" ," Jet Towers Area ",50,0,1);
  TH1Dcontainer_["JetTwrEt"] = subDir.make<TH1D>("JetTwrEt" ,"jet towers Et ",50,0,100);
  TH1Dcontainer_["JetTwrEta"] = subDir.make<TH1D>("JetTwrEta" ,"jet towers Eta ",50,-6, 6);
  TH1Dcontainer_["JetTwrPhi"] = subDir.make<TH1D>("JetTwrPhi" ,"jet towers Phi ",50,-3.2, 3.2);
  TH1Dcontainer_["JetTrkPt"] = subDir.make<TH1D>("JetTrkPt" ,"jet tracks Pt ",50,0, 100);
  TH1Dcontainer_["JetVtrTrkSize"] = subDir.make<TH1D>("JetVtrTrkSize" ,"Number of tracks associated to Primary vertex ",10,0, 10);
  TH1Dcontainer_["JetVrtTrkPt"] = subDir.make<TH1D>("JetVrtTrkPt" ,"track associated to primavy vertex Pt ",50,0, 100);
  // TH1Dcontainer_["njetsInSubDir"] = subDir.make<TH1D>("njets-2" ,"jet multiplicity for jets with E_{T} > 30 GeV",20,0,20);//wrote in subdirectory
  // TH1Dcontainer_["njetsInSubSubDir"] = subsubDir.make<TH1D>("njets-3" ,"jet multiplicity for jets with E_{T} > 30 GeV",20,0,20);//wrote in subsubdirectory
  //TH1F
  //TH2F
}

// ------------ method called once each job just after ending the event loop  ------------
void 
JetMetChecker::endJob() {
   //use LogError to summarise the error that happen in the execution (by example from warning) (ex: Nof where we cannot access such variable)
   edm::LogError  ("SummaryError") << "My error message \n";    // or  edm::LogProblem  (not formated)
   //use LogInfo to summarise information (ex: pourcentage of events matched ...)
   edm::LogInfo   ("MainResults") << "My LogInfo message \n";  // or  edm::LogVerbatim (not formated)
}

//define this as a plug-in
DEFINE_FWK_MODULE(JetMetChecker);