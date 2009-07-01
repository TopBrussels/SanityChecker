// -*- C++ -*-
//
// Package:    VertexChecker
// Class:      VertexChecker
// 
/**\class VertexChecker VertexChecker.cc UserCode/DummyChecker/src/VertexCheckery.cc

 Description: <one line class summary>

 Implementation:
     <Notes on implementation>
*/
//
// Original Author:  local user
//         Created:  Wed Feb 18 16:39:03 CET 2009
// $Id: VertexChecker.cc,v 1.4 2009/03/20 16:54:52 villella Exp $
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


#include "DataFormats/VertexReco/interface/Vertex.h"

#include "TDirectory.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TH1F.h"
#include "TH2F.h"

//
// class decleration
//

class VertexChecker : public edm::EDAnalyzer {
   public:
      explicit VertexChecker(const edm::ParameterSet&);
      ~VertexChecker();


   private:
      virtual void beginJob(const edm::EventSetup&) ;
      virtual void analyze(const edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;

      // ----------member data ---------------------------
  edm::InputTag vertex_;
   
  //Histograms are booked in the beginJob() method
  std::map<std::string,TDirectory*> TDirectorycontainer_; // simple map to contain all TDirectory.
  std::map<std::string,TH1D*> TH1Dcontainer_;       // simple map to contain all TH1D.
  std::map<std::string,TH2D*> TH2Dcontainer_;       // simple map to contain all TH2D.
  std::map<std::string,TH1F*> TH1Fcontainer_;       // simple map to contain all TH1F.
  std::map<std::string,TH2F*> TH2Fcontainer_;       // simple map to contain all TH2F.
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
VertexChecker::VertexChecker(const edm::ParameterSet& iConfig)

{
   //now do what ever initialization is needed
  vertex_        =   iConfig.getParameter<edm::InputTag>( "vertexName" );

}


VertexChecker::~VertexChecker()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called to for each event  ------------
void
VertexChecker::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
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
   Handle< std::vector<reco::Vertex> > vertex;
   iEvent.getByLabel(vertex_, vertex);
   
   //Check if branches are available
   if (!vertex.isValid()){
     edm::LogWarning  ("LinkBroken") << "My warning message - NoVertexFound"; 
     return;//   throw cms::Exception("ProductNotFound") <<"Primary vertex collection not found"<<std::endl;
   }
   
   double pt2_vtx =0;
   double sumpt2_vtx =0;
   double ptAssTrk = -100;
   int prim_vtxidx = -2;
   double sumpt =0;
   double sumpt_prim =0;
   double pt2_primvtx =  0;
   double sumpt2_primvtx =0;
   //  Primary Vertex infos;
   if(vertex->size() > 0){
     TH1Dcontainer_["VtxNumber"]->Fill(vertex->size());
     for (unsigned int j=0; j< vertex->size(); j++){

      TH1Dcontainer_["Vrt_x"]->Fill((*vertex)[j].x());
      TH1Dcontainer_["Vrt_y"]->Fill((*vertex)[j].y());
      TH1Dcontainer_["Vrt_z"]->Fill((*vertex)[j].z());
      TH2Dcontainer_["Vrt_xy"]->Fill((*vertex)[j].x(),(*vertex)[j].y() );
       TH1Dcontainer_["VrtTrkNumber"]->Fill((*vertex)[j].tracksSize());
       reco::Vertex::trackRef_iterator tr ;  
       for ( tr = (*vertex)[j].tracks_begin(); tr !=(*vertex)[j].tracks_end(); tr++){
	 pt2_vtx =  (*tr)->pt() * (*tr)->pt() ;
	 sumpt2_vtx += pt2_vtx;
	 sumpt += (*tr)->pt();
	 TH1Dcontainer_["VrtTrkPt"]->Fill((*tr)->pt());
       }
       if( sumpt2_vtx > ptAssTrk ){ //defines as "primary" the vertex associated with the highest sum of pt2 of the tracks in the event
	 ptAssTrk = sumpt2_vtx;
	 prim_vtxidx =j; 
      }
       //      if(j !=prim_vtxidx) sec_vtxidx =j;
      }//close for vertex
	 TH1Dcontainer_["PrimaryVrt_x"]->Fill((*vertex)[ prim_vtxidx].x());
	 TH1Dcontainer_["PrimaryVrt_y"]->Fill((*vertex)[ prim_vtxidx].y());
	 TH1Dcontainer_["PrimaryVrt_z"]->Fill((*vertex)[ prim_vtxidx].z());   
	 TH1Dcontainer_["VrtTrkSumPt"]->Fill(sumpt);
	 TH2Dcontainer_["PrimaryVrt_xy"]->Fill((*vertex)[ prim_vtxidx].x(), (*vertex)[ prim_vtxidx].y());
	 TH1Dcontainer_["PrimaryVrtTrkNumber"]->Fill((*vertex)[prim_vtxidx].tracksSize());
       reco::Vertex::trackRef_iterator trp ;  
    for ( trp = (*vertex)[ prim_vtxidx].tracks_begin(); trp !=(*vertex)[ prim_vtxidx].tracks_end(); trp++){
	 sumpt_prim += (*trp)->pt();
	 pt2_primvtx =  (*trp)->pt() * (*trp)->pt() ;
	 sumpt2_primvtx += pt2_primvtx;
    }
	 TH1Dcontainer_["PrimVrtTrkSumPt"]->Fill(sumpt_prim);
	 TH1Dcontainer_["PrimVrtTrkSumPt2"]->Fill(sumpt2_primvtx);
   }//close if
   
 
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
VertexChecker::beginJob(const edm::EventSetup&)
{
  edm::Service<TFileService> fs;
  if (!fs) throw edm::Exception(edm::errors::Configuration, "TFileService missing from configuration!");
  
  TFileDirectory subDir = fs->mkdir( "Vertex" );
  // TFileDirectory subsubDir = subDir.mkdir( "Vertex" );
  
  //define the histograms booked
  //TH1D
  
  TH1Dcontainer_["VtxNumber"] = subDir.make<TH1D>("VtxNumber" ,"Number of vertices ",5,0, 5);
  TH1Dcontainer_["Vrt_x"]= subDir.make<TH1D>("Vtx_x" ,"vertices x ",50,-1, 1);
  TH1Dcontainer_["Vrt_y"]= subDir.make<TH1D>("Vtx_y" ,"vertices y ",50,-1, 1);
  TH1Dcontainer_["Vrt_z"]= subDir.make<TH1D>("Vtx_z" ,"vertices z ",100,-10, 10);
  TH2Dcontainer_["Vrt_xy"]= subDir.make<TH2D>("Vtx_xy" ,"vertices x vs y ",50,-1, 1, 50, -1, 1);
  TH1Dcontainer_["PrimaryVrt_x"]= subDir.make<TH1D>("PrimaryVtx_x" ,"Primary vertex x ",50,-1, 1);
  TH1Dcontainer_["PrimaryVrt_y"]= subDir.make<TH1D>("PrimaryVtx_y" ,"Primary vertex y ",50,-1, 1);
  TH1Dcontainer_["PrimaryVrt_z"]= subDir.make<TH1D>("PrimaryVtx_z" ,"Primary vertex z ",100,-10, 10);
  TH2Dcontainer_["PrimaryVrt_xy"] =  subDir.make<TH2D>("PrimaryVtx_xy" ,"Primary vertex x vs y ",50,-1, 1, 50, -1, 1);
  TH1Dcontainer_["PrimaryVrtTrkNumber"]= subDir.make<TH1D>("PrimaryVrtTrkNumber","Number of tracks associated to primary vertex",10,0, 10);
  TH1Dcontainer_["VrtTrkNumber"] = subDir.make<TH1D>("VrtTrkNumber" ,"Number of tracks associated to vertices ",10,0, 10);
  TH1Dcontainer_["VrtTrkPt"] = subDir.make<TH1D>("VrtTrkPt" ,"track associated to vertices Pt ",100,0, 50);
  TH1Dcontainer_["VrtTrkSumPt"]= subDir.make<TH1D>("VrtTrkSumPt" ," Sum of the pt of the tracks associated to vertices",100,0, 300);
  TH1Dcontainer_["PrimVrtTrkSumPt"]= subDir.make<TH1D>("PrimVrtTrkSumPt" ," Sum of the pt of the tracks associated to Primary vertex",50,0, 300);
  TH1Dcontainer_["PrimVrtTrkSumPt2"]= subDir.make<TH1D>("PrimVrtTrkSumPt2" ,"Sum of the pt^2 of the tracks associated to Primary vertex",50,0, 10000);
}

// ------------ method called once each job just after ending the event loop  ------------
void 
VertexChecker::endJob() {
   //use LogError to summarise the error that happen in the execution (by example from warning) (ex: Nof where we cannot access such variable)
   //edm::LogError  ("SummaryError") << "My error message \n";    // or  edm::LogProblem  (not formated)
   //use LogInfo to summarise information (ex: pourcentage of events matched ...)

  edm::LogVerbatim ("MainResults") << " -------------------------------------------";
  edm::LogVerbatim ("MainResults") << " -------------------------------------------";
  edm::LogVerbatim ("MainResults") << " --     Report from Vertex Checker       -- ";
  edm::LogVerbatim ("MainResults") << " -------------------------------------------";
  edm::LogVerbatim ("MainResults") << " -------------------------------------------";
  edm::LogVerbatim ("MainResults") << "Mean Number of Vertices  = " << TH1Dcontainer_["VtxNumber"]->GetMean() ;  
  edm::LogVerbatim ("MainResults") << "Vertices Mean x = " << TH1Dcontainer_["Vrt_x"]->GetMean() ;  
  edm::LogVerbatim ("MainResults") << "Vertices Mean y = " << TH1Dcontainer_["Vrt_y"]->GetMean() ;  
  edm::LogVerbatim ("MainResults") << "Vertices Mean z = " << TH1Dcontainer_["Vrt_z"]->GetMean() ;  
  edm::LogVerbatim ("MainResults") << " -------------------------------------------";
  edm::LogVerbatim ("MainResults") << " -----    Primary Vertex Info  -------------";
  edm::LogVerbatim ("MainResults") << "Prim Vertex Mean x = " << TH1Dcontainer_["PrimaryVrt_x"]->GetMean() ;  
  edm::LogVerbatim ("MainResults") << "Prim Vertex Mean y = " << TH1Dcontainer_["PrimaryVrt_y"]->GetMean() ;  
  edm::LogVerbatim ("MainResults") << "Prim Vertex Mean z = " << TH1Dcontainer_["PrimaryVrt_z"]->GetMean() ;  

}

//define this as a plug-in
DEFINE_FWK_MODULE(VertexChecker);
