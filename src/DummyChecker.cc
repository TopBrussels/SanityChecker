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
// $Id$
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

#include "TDirectory.h"
#include "TH1D.h"
#include "TH1F.h"
#include "TH2F.h"

//
// class decleration
//

class DummyChecker : public edm::EDAnalyzer {
   public:
      explicit DummyChecker(const edm::ParameterSet&);
      ~DummyChecker();


   private:
      virtual void beginJob(const edm::EventSetup&) ;
      virtual void analyze(const edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;

      // ----------member data ---------------------------
      edm::InputTag genEventCollectionName_;

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
DummyChecker::DummyChecker(const edm::ParameterSet& iConfig)

{
   //now do what ever initialization is needed
   genEventCollectionName_ = iConfig.getParameter<edm::InputTag>( "genEventCollectionName" );


}


DummyChecker::~DummyChecker()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called to for each event  ------------
void
DummyChecker::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
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
   Handle<TtGenEvent> TtGenEvent_;
   iEvent.getByLabel(genEventCollectionName_,TtGenEvent_);

   TH1Dcontainer_["njets"]->Fill(1);
   TH1Dcontainer_["njetsInSubDir"]->Fill(2);
   TH1Dcontainer_["njetsInSubSubDir"]->Fill(3);

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
DummyChecker::beginJob(const edm::EventSetup&)
{
  edm::Service<TFileService> fs;
  if (!fs) throw edm::Exception(edm::errors::Configuration, "TFileService missing from configuration!");

   TFileDirectory subDir = fs->mkdir( "mySubDirectory" );
   TFileDirectory subsubDir = subDir.mkdir( "mySubSubDirectory" );

  //define the histograms booked
  //TH1D
  TH1Dcontainer_["njets"] = fs->make<TH1D>("njets" ,"jet multiplicity for jets with E_{T} > 30 GeV",20,0,20);// wrote in main directory
  TH1Dcontainer_["njetsInSubDir"] = subDir.make<TH1D>("njets-2" ,"jet multiplicity for jets with E_{T} > 30 GeV",20,0,20);//wrote in subdirectory
  TH1Dcontainer_["njetsInSubSubDir"] = subsubDir.make<TH1D>("njets-3" ,"jet multiplicity for jets with E_{T} > 30 GeV",20,0,20);//wrote in subsubdirectory
  //TH1F
  //TH2F
}

// ------------ method called once each job just after ending the event loop  ------------
void 
DummyChecker::endJob() {
   //use LogError to summarise the error that happen in the execution (by example from warning) (ex: Nof where we cannot access such variable)
   edm::LogError  ("SummaryError") << "My error message \n";    // or  edm::LogProblem  (not formated)
   //use LogInfo to summarise information (ex: pourcentage of events matched ...)
   edm::LogInfo   ("MainResults") << "My LogInfo message \n";  // or  edm::LogVerbatim (not formated)
}

//define this as a plug-in
DEFINE_FWK_MODULE(DummyChecker);
