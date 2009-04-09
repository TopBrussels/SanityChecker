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
// $Id: KinematicsChecker.cc,v 1.10 2009/04/08 12:34:47 echabert Exp $
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
#include <vector>
#include <Math/VectorUtil.h>

//
// class decleration
//

/*struct HighestPt{
  bool operator()( reco::Candidate * j1, reco::Candidate * j2 ) const{
  //return j1->p4().Pt() > j2->p4().Pt() ;
  return true;
  }
  };*/

struct Highest{
  bool operator()( double j1, double j2 ) const{
    return j1 > j2 ;
  }
};

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
  std::map<std::string,TH1D*> TH1Dcontainer_[8]; // simple map to contain all TH1D.
  std::map<std::string,TH1F*> TH1Fcontainer_; // simple map to contain all TH1F.
  std::map<std::string,TH2F*> TH2Fcontainer_; // simple map to contain all TH2F.

  std::vector< double > jetsAcceptance_, muonsAcceptance_;
  std::vector< double > nJetsAcceptance, nMuonsAcceptance;
  
  bool verbose_;
  int NbOfEvents;
   
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
  jetsAcceptance_   =   iConfig.getParameter< std::vector< double > >( "jetsAcceptance" );
  muonsAcceptance_  =   iConfig.getParameter< std::vector< double > >( "muonsAcceptance" );
  verbose_          =   iConfig.getParameter<bool>( "verbose" );
  NbOfEvents        = 0;
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
  NbOfEvents++;
  
  Handle< std::vector<pat::Jet> > jets;
  iEvent.getByLabel(jets_, jets);
      
  Handle< std::vector<pat::Muon> > muons;
  iEvent.getByLabel(muons_, muons);

  Handle< std::vector<pat::MET> > mets;
  iEvent.getByLabel(mets_, mets);


  //Check if branches are available
  if (!jets.isValid()){
    edm::LogWarning  ("LinkBroken_noJetsFound") << "My warning message - No Jets Found"; 
    return; //throw cms::Exception("ProductNotFound") <<"jet collection not found"<<std::endl;
  }
  if (!muons.isValid()){
    edm::LogWarning  ("LinkBroken_NoMuonsFound") << "My warning message - No Muons Found"; 
    return; //throw cms::Exception("ProductNotFound") <<"muon collection not found"<<std::endl;
  }
  if (!mets.isValid()){
    edm::LogWarning  ("LinkBroken_NoMetsFound") << "My warning message - No Mets Found";
    return;// throw cms::Exception("ProductNotFound") <<"MET collection not found"<<std::endl;
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
   

  Handle<TtGenEvent> genEvt;
  iEvent.getByLabel ("genEvt",genEvt);
  
  //Check if branch is available  
  if (!genEvt.isValid()){
    edm::LogWarning  ("LinkBroken_noGenEvtFound") << "KinematicsChecker - No GenEvt Found";
  }
  
  if(  jets->size() >= 4){  
    if (genEvt.isValid()){
     
      if(genEvt->isSemiLeptonic(WDecay::kMuon)) {

	//make a copy of the jet collection to be able to drop jets from it
	std::vector<pat::Jet> jets_clone;
	for(unsigned int i=0; i<jets->size(); i++){
	  jets_clone.push_back((*jets)[i]);
	}
      
	if(verbose_){std::cout << " jets->size(): " << jets->size() << std::endl;}
	if(verbose_){std::cout << " jets_clone.size(): " << jets_clone.size() << std::endl;}

	// Make plots only with jets from ttbar events (matched with a certain algo)
	// Matching index : Hadronic Q  = 0, Hadronic Q' = 1, Hadronic b  = 2, Leptonic b  = 3;
	std::vector<const reco::Candidate *> TopQuarks;
	if(genEvt->b() && genEvt->bBar() && genEvt->hadronicDecayQuark() && genEvt->hadronicDecayQuarkBar()){
	  TopQuarks.push_back(genEvt->hadronicDecayQuark());
	  TopQuarks.push_back(genEvt->hadronicDecayQuarkBar());
	  TopQuarks.push_back(genEvt->hadronicDecayB());
	  TopQuarks.push_back(genEvt->leptonicDecayB());
	}
	if(TopQuarks.size()==4) { 

	  JetPartonMatching *GenMatchTopQuarks = new JetPartonMatching(TopQuarks, jets_clone, matchingAlgo_, useMaxDist_, useDeltaR_, maxDist_);

	  for(unsigned int i=0; i<4; i++){
	    Int_t Idx = GenMatchTopQuarks->getMatchForParton(i,0);
	    if(Idx>=0){
	      ObjectP4s[0].push_back((jets_clone[Idx]).p4());
	      //drop jets from collection
	      jets_clone.erase(jets_clone.begin()+Idx);
	    }
	  }
	  
	  if(verbose_){std::cout << " jets_clone.size(): " << jets_clone.size() << std::endl;}

	  TH1Dcontainer_[3]["number"]->Fill(ObjectP4s[0].size()); 
	  for(unsigned int j=0; j<ObjectP4s[0].size(); j++) TH1Dcontainer_[3]["pt"]->Fill(ObjectP4s[0][j].Pt());  
	  for(unsigned int j=0;j<ObjectP4s[0].size();j++) TH1Dcontainer_[3]["eta"]->Fill(ObjectP4s[0][j].Eta()); 
	  for(unsigned int j=0;j<ObjectP4s[0].size();j++) TH1Dcontainer_[3]["et"]->Fill(ObjectP4s[0][j].Et());
	  for(unsigned int j=0;j<ObjectP4s[0].size();j++) TH1Dcontainer_[3]["energy"]->Fill(ObjectP4s[0][j].E());
	  for(unsigned int j=0;j<ObjectP4s[0].size();j++) TH1Dcontainer_[3]["theta"]->Fill(ObjectP4s[0][j].Theta());
	  for(unsigned int j=0;j<ObjectP4s[0].size();j++) TH1Dcontainer_[3]["phi"]->Fill(ObjectP4s[0][j].Phi());

	}
       
	//  Make plots only with jets from ISR (matched with a certain algo)  

	for(unsigned int i=0;i<3;i++){
	  ObjectP4s[i].clear();
	}
   
	std::vector<const reco::Candidate *> ISRquarks;
	std::vector<const reco::Candidate *> ISRquarksSorted;
	std::vector<double> IRSquarksPt;

	unsigned int NumberISR = genEvt->topSisters().size();	
	if(NumberISR>0) {
	  if(verbose_){std::cout << " NumberISR: " <<NumberISR << std::endl;}

	  //Fill vector with ISR quarks
	  for(unsigned int i=0; i<NumberISR; i++){
	    ISRquarks.push_back((genEvt->topSisters())[i]);
	    IRSquarksPt.push_back((genEvt->topSisters())[i]->p4().pt());
	  }
	 
	  //sort vector
	  //std::sort(ISRquarks.begin(),ISRquarks.begin(),HighestPt());
	  //sort by ugly method
	  std::sort(IRSquarksPt.begin(),IRSquarksPt.end(),Highest());
		   
	  for(unsigned int i=0; i<NumberISR; i++){
	    for(unsigned int j=0; j<NumberISR; j++){
	      double diff=fabs(IRSquarksPt[i]-ISRquarks[j]->p4().pt());
	      if(diff<0.001){
		ISRquarksSorted.push_back(ISRquarks[j]);
	      }
	    }
	  }
	 
	  if(verbose_){std::cout << " ISRquarks.size(): " << ISRquarks.size() << std::endl;}
	  if(verbose_){std::cout << " ISRquarksSorted.size(): " << ISRquarksSorted.size() << std::endl;}
	  //Remove quarks if there are more then jets, to have unambiguios jet parton matching
	  if(NumberISR > jets_clone.size()){  
	    edm::LogWarning  ("LinkBroken_ISRJets") << "there are more quarks from initial state radiation than there are jets in the jet collection";//TODO: waht to do with "category"?
	    unsigned int diffSize=ISRquarksSorted.size()-jets_clone.size();
	    if(diffSize>0){
	      for(unsigned int i = 0; i<diffSize; i++ ){
		ISRquarksSorted.pop_back();
	      }
	    }
	  }
	  if(verbose_){std::cout << " ISRquarksSorted.size(): " << ISRquarksSorted.size() << std::endl;}
	  
	  JetPartonMatching GenMatchISRquarks(ISRquarksSorted, jets_clone, matchingAlgo_, useMaxDist_, useDeltaR_, maxDist_); 
	 
	  if(verbose_){std::cout << " jets_clone.size(): " << jets_clone.size() << std::endl;}
	  for(unsigned int i=0; i<ISRquarksSorted.size(); i++){
	    Int_t Idx = GenMatchISRquarks.getMatchForParton(i,0);
	    if(Idx>=0){
	      ObjectP4s[0].push_back((jets_clone[Idx]).p4()); 
	      //drop jets from collection
	      jets_clone.erase(jets_clone.begin()+Idx);
	    }
	  }
	  if(verbose_){std::cout << " jets_clone.size(): " << jets_clone.size() << std::endl;}

	  TH1Dcontainer_[4]["number"]->Fill(ObjectP4s[0].size()); 
	  for(unsigned int j=0; j<ObjectP4s[0].size(); j++) TH1Dcontainer_[4]["pt"]->Fill(ObjectP4s[0][j].Pt());  
	  for(unsigned int j=0;j<ObjectP4s[0].size();j++) TH1Dcontainer_[4]["eta"]->Fill(ObjectP4s[0][j].Eta()); 
	  for(unsigned int j=0;j<ObjectP4s[0].size();j++) TH1Dcontainer_[4]["et"]->Fill(ObjectP4s[0][j].Et());
	  for(unsigned int j=0;j<ObjectP4s[0].size();j++) TH1Dcontainer_[4]["energy"]->Fill(ObjectP4s[0][j].E());
	  for(unsigned int j=0;j<ObjectP4s[0].size();j++) TH1Dcontainer_[4]["theta"]->Fill(ObjectP4s[0][j].Theta());
	  for(unsigned int j=0;j<ObjectP4s[0].size();j++) TH1Dcontainer_[4]["phi"]->Fill(ObjectP4s[0][j].Phi());

	}

	//  Make plots only with jets from topRadiation (matched with a certain algo)  
	for(unsigned int i=0;i<3;i++){
	  ObjectP4s[i].clear();
	}

	std::vector<const reco::Candidate *> TopRadQuarks;
	std::vector<const reco::Candidate *> TopRadQuarksSorted;  
	std::vector<double> TopRadQuarksPt;
       
	unsigned int NumberTopRadLep = genEvt->leptonicDecayTopRadiation().size();		
	unsigned int NumberTopRadHad = genEvt->hadronicDecayTopRadiation().size();		
	unsigned int NumberTopRad = NumberTopRadLep+NumberTopRadHad;	
       
	if(NumberTopRadLep!=0 ) {
	  for(unsigned int i=0;i<NumberTopRadLep; i++){
	    TopRadQuarks.push_back((genEvt->leptonicDecayTopRadiation())[i]);
	    TopRadQuarksPt.push_back((genEvt->leptonicDecayTopRadiation())[i]->p4().pt());
	  }
	}
	if(NumberTopRadHad!=0 ) {
	  for(unsigned int i=0; i<NumberTopRadHad; i++){
	    TopRadQuarks.push_back((genEvt->hadronicDecayTopRadiation())[i]);
	    TopRadQuarksPt.push_back((genEvt->hadronicDecayTopRadiation())[i]->p4().pt());
	  }
	}

	if(NumberTopRad>0){ 
 
	  //sort by ugly method
	  std::sort(TopRadQuarksPt.begin(),TopRadQuarksPt.end(),Highest());
	 
	  for(unsigned int i=0; i<NumberTopRad; i++){
	    for(unsigned int j=0; j<NumberTopRad; j++){
	      double diff=fabs(TopRadQuarksPt[i]-TopRadQuarks[j]->p4().pt());
	      if(diff<0.001){
		TopRadQuarksSorted.push_back(TopRadQuarks[j]);
	      }
	    }
	  }
	  if(verbose_){std::cout << " jets_clone.size(): " << jets_clone.size() << std::endl;} 
	  if(verbose_) std::cout << " TopRadQuarksSorted.size() " << TopRadQuarksSorted.size() << std::endl; 
	  if(verbose_) std::cout << " TopRadQuarks.size() " << TopRadQuarks.size() << std::endl; 
	  //Remove quarks if there are more than jets, to have unambiguios jet parton matching
	  if(NumberTopRad > jets_clone.size()){
	    edm::LogWarning  ("LinkBroken_TopRadJets") << "there are more quarks from top radiation than there are jets in the jet collection";//TODO: waht to do with "category"?
	    unsigned int diffSize=TopRadQuarksSorted.size()-jets_clone.size();
	    if(verbose_){std::cout << " diffSize (TopRadQuarksSorted.size()-jets_clone.size()): " <<  diffSize << std::endl;}
	    for(unsigned int i = 0; i<diffSize; i++ ){
	      TopRadQuarksSorted.pop_back();
	    }
	  }
	
	  JetPartonMatching *GenMatchTopRadQuarks = new JetPartonMatching(TopRadQuarksSorted, jets_clone, matchingAlgo_, useMaxDist_, useDeltaR_, maxDist_); 

	  for(unsigned int i=0; i<TopRadQuarksSorted.size(); i++){
	    Int_t Idx = GenMatchTopRadQuarks->getMatchForParton(i,0);
	    if(Idx>=0) ObjectP4s[0].push_back((jets_clone[Idx]).p4());
	  }

	  TH1Dcontainer_[5]["number"]->Fill(ObjectP4s[0].size()); 
	  for(unsigned int j=0; j<ObjectP4s[0].size(); j++) TH1Dcontainer_[5]["pt"]->Fill(ObjectP4s[0][j].Pt());  
	  for(unsigned int j=0;j<ObjectP4s[0].size();j++) TH1Dcontainer_[5]["eta"]->Fill(ObjectP4s[0][j].Eta()); 
	  for(unsigned int j=0;j<ObjectP4s[0].size();j++) TH1Dcontainer_[5]["et"]->Fill(ObjectP4s[0][j].Et());
	  for(unsigned int j=0;j<ObjectP4s[0].size();j++) TH1Dcontainer_[5]["energy"]->Fill(ObjectP4s[0][j].E());
	  for(unsigned int j=0;j<ObjectP4s[0].size();j++) TH1Dcontainer_[5]["theta"]->Fill(ObjectP4s[0][j].Theta());
	  for(unsigned int j=0;j<ObjectP4s[0].size();j++) TH1Dcontainer_[5]["phi"]->Fill(ObjectP4s[0][j].Phi());


	}
     
	for(unsigned int i=0;i<3;i++){
	  ObjectP4s[i].clear();
	}
       
	//Make plots for matched muon 
	if(muons->size()>0){

	  double ClosestDelR=999.;
	  int ClosestIdx=-1;

	  for(unsigned int i=0; i<muons->size();i++){
	 	
	    if(ROOT::Math::VectorUtil::DeltaR(genEvt->singleLepton()->p4(),(*muons)[i].p4())<ClosestDelR){
	      ClosestDelR=ROOT::Math::VectorUtil::DeltaR(genEvt->singleLepton()->p4(),(*muons)[i].p4());
	      ClosestIdx=i;
	    }
	  }

	  if(ClosestDelR<0.3){
	    ObjectP4s[1].push_back(((*muons)[ClosestIdx].p4()));
	   
	    TH1Dcontainer_[6]["number"]->Fill(ObjectP4s[1].size()); 
	    TH1Dcontainer_[6]["pt"]->Fill(ObjectP4s[1][0].Pt());  
	    TH1Dcontainer_[6]["eta"]->Fill(ObjectP4s[1][0].Eta()); 
	    TH1Dcontainer_[6]["et"]->Fill(ObjectP4s[1][0].Et());
	    TH1Dcontainer_[6]["energy"]->Fill(ObjectP4s[1][0].E());
	    TH1Dcontainer_[6]["theta"]->Fill(ObjectP4s[1][0].Theta());
	    TH1Dcontainer_[6]["phi"]->Fill(ObjectP4s[1][0].Phi());
	   
	    TH1Fcontainer_["resEtMuons"]->Fill((genEvt->singleLepton()->p4().Et()-ObjectP4s[1][0].Et())/ObjectP4s[1][0].Et());
	    TH1Fcontainer_["resPtMuons"]->Fill((genEvt->singleLepton()->p4().Pt()-ObjectP4s[1][0].Pt())/ObjectP4s[1][0].Pt()); 
	    TH1Fcontainer_["resPhiMuons"]->Fill((genEvt->singleLepton()->p4().Phi()-ObjectP4s[1][0].Phi())/ObjectP4s[1][0].Phi()); 
	    TH1Fcontainer_["resThetaMuons"]->Fill((genEvt->singleLepton()->p4().Theta()-ObjectP4s[1][0].Theta())/ObjectP4s[1][0].Theta()); 
	    TH1Fcontainer_["resEtaMuons"]->Fill((genEvt->singleLepton()->p4().Eta()-ObjectP4s[1][0].Eta())/ObjectP4s[1][0].Eta()); 


	  }
	  if(ClosestDelR==999.){ClosestDelR=-1;}

	  TH1Fcontainer_["deltaRMuons"]->Fill(ClosestDelR);

	}
  
	for(unsigned int i=0;i<3;i++){
	  ObjectP4s[i].clear();
	}
       
	if(mets->size()>0){

	  double ClosestDelPhi=999.;
	  ClosestDelPhi=fabs(ROOT::Math::VectorUtil::DeltaPhi(genEvt->singleNeutrino()->p4(),(*mets)[0].p4()));
	 

	  ObjectP4s[2].push_back(((*mets)[0].p4()));
	
	  TH1Dcontainer_[7]["number"]->Fill(ObjectP4s[2].size()); 
	  TH1Dcontainer_[7]["pt"]->Fill(ObjectP4s[2][0].Pt());  
	  TH1Dcontainer_[7]["eta"]->Fill(ObjectP4s[2][0].Eta()); 
	  TH1Dcontainer_[7]["et"]->Fill(ObjectP4s[2][0].Et());
	  TH1Dcontainer_[7]["energy"]->Fill(ObjectP4s[2][0].E());
	  TH1Dcontainer_[7]["theta"]->Fill(ObjectP4s[2][0].Theta());
	  TH1Dcontainer_[7]["phi"]->Fill(ObjectP4s[2][0].Phi());
	 
	  TH1Fcontainer_["resEtMets"]->Fill((genEvt->singleLepton()->p4().Et()-ObjectP4s[2][0].Et())/ObjectP4s[2][0].Et());
	  TH1Fcontainer_["resPtMets"]->Fill((genEvt->singleLepton()->p4().Pt()-ObjectP4s[2][0].Pt())/ObjectP4s[2][0].Pt()); 
	  TH1Fcontainer_["resPhiMets"]->Fill((genEvt->singleLepton()->p4().Phi()-ObjectP4s[2][0].Phi())/ObjectP4s[2][0].Phi()); 
	  TH1Fcontainer_["resPxMets"]->Fill((genEvt->singleLepton()->p4().Px()-ObjectP4s[2][0].Px())/ObjectP4s[2][0].Px()); 
	  TH1Fcontainer_["resPyMets"]->Fill((genEvt->singleLepton()->p4().Py()-ObjectP4s[2][0].Py())/ObjectP4s[2][0].Py()); 
 
	  if(ClosestDelPhi==999.){ClosestDelPhi=-1;}
	  TH1Fcontainer_["deltaPhiMets"]->Fill(ClosestDelPhi);
	 
	}
	 
      }//issemileptonic
    }//genevent.isvalid
  }//more than 4 jets
}


// ------------ method called once each job just before starting event loop  ------------
void 
KinematicsChecker::beginJob(const edm::EventSetup&)
{
  edm::Service<TFileService> fs;
  if (!fs) throw edm::Exception(edm::errors::Configuration, "TFileService missing from configuration!");
 
  std::string ObjectNames[8] = {"Jets","Muons","METs","MatchedTopJets","MatchedISRjets","MatchedTopRadJets","MatchedMuons","MatchedMets"};

  std::vector< TFileDirectory > subDirs;
  for(unsigned int i=0;i<8;i++) subDirs.push_back(fs->mkdir( ObjectNames[i] ));

  for(unsigned int i=0;i<8;i++){  
    TH1Dcontainer_[i]["number"] = subDirs[i].make<TH1D>("number" ,"number of objects",50,0,30);
    TH1Dcontainer_[i]["pt"] = subDirs[i].make<TH1D>("pt" ,"pt",50,0,200);
    TH1Dcontainer_[i]["et"] = subDirs[i].make<TH1D>("et" ,"et",50,0,200);
    TH1Dcontainer_[i]["eta"] = subDirs[i].make<TH1D>("eta" ,"eta",50,-6,6);
    TH1Dcontainer_[i]["energy"] = subDirs[i].make<TH1D>("energy" ,"energy",50,0,400);
    TH1Dcontainer_[i]["theta"] = subDirs[i].make<TH1D>("theta" ,"theta",50,0,3.2);
    TH1Dcontainer_[i]["phi"] = subDirs[i].make<TH1D>("phi" ,"phi",50,-3.2,3.2);
  }
 
  
  TH1Fcontainer_["deltaRMuons"] = subDirs[6].make<TH1F>("deltaRMuons" ,"deltaR of closest reconstructed muon to generated muon",50,-1,5);
  TH1Fcontainer_["resEtMuons"] = subDirs[6].make<TH1F>("resEtMuons" ,"resolution Et (gen-rec)/rec",50,-.2,.5);
  TH1Fcontainer_["resPtMuons"] = subDirs[6].make<TH1F>("resPtMuons" ,"resolution Pt (gen-rec)/rec",50,-.2,.5);
  TH1Fcontainer_["resPhiMuons"] = subDirs[6].make<TH1F>("resPhiMuons" ,"resolution Phi (gen-rec)/rec",50,-.1,.1);
  TH1Fcontainer_["resThetaMuons"] = subDirs[6].make<TH1F>("resThetaMuons" ,"resolution Theta (gen-rec)/rec",50,-.1,.1);
  TH1Fcontainer_["resEtaMuons"] = subDirs[6].make<TH1F>("resEtaMuons" ,"resolution Eta (gen-rec)/rec",50,-.1,.1);


  TH1Fcontainer_["deltaPhiMets"] = subDirs[7].make<TH1F>("deltaPhiMets" ,"deltaPhi of neutrino and missing transverse energy",50,-1,5);
  TH1Fcontainer_["resEtMets"] = subDirs[7].make<TH1F>("resEtMets" ,"resolution Et (gen-rec)/rec",50,-1,3);
  TH1Fcontainer_["resPtMets"] = subDirs[7].make<TH1F>("resPtMets" ,"resolution Pt (gen-rec)/rec",50,-1,3);
  TH1Fcontainer_["resPhiMets"] = subDirs[7].make<TH1F>("resPhiMets" ,"resolution Phi (gen-rec)/rec",50,-4,2);
  TH1Fcontainer_["resPxMets"] = subDirs[7].make<TH1F>("resPxMets" ,"resolution Px (gen-rec)/rec",50,-3,3);
  TH1Fcontainer_["resPyMets"] = subDirs[7].make<TH1F>("resPyMets" ,"resolution Py (gen-rec)/rec",50,-3,3);
 

  nJetsAcceptance.push_back(0);
  nJetsAcceptance.push_back(0);
  nMuonsAcceptance.push_back(0);
  nMuonsAcceptance.push_back(0);

}

// ------------ method called once each job just after ending the event loop  ------------
void 
KinematicsChecker::endJob() {

  edm::LogVerbatim ("SummaryError") << " -------------------------------------------";
  edm::LogVerbatim ("SummaryError") << " -------------------------------------------";
  edm::LogVerbatim ("SummaryError") << " --   Report from Kinematics Checker     -- ";
  edm::LogVerbatim ("SummaryError") << " -------------------------------------------";
  edm::LogVerbatim ("SummaryError") << " -------------------------------------------";
  edm::LogVerbatim ("SummaryError") << " ";
  edm::LogVerbatim ("SummaryError") << " Number of events this module ran on: " << NbOfEvents;

  edm::LogVerbatim ("SummaryError") << " ";
  edm::LogVerbatim ("SummaryError") << " -------------------------------";
  edm::LogVerbatim ("SummaryError") << "  Info on acceptance";
  edm::LogVerbatim ("SummaryError") << " -------------------------------";
  edm::LogVerbatim ("SummaryError") << " ";

  edm::LogProblem   ("SummaryError") << "Number of events where at least one jet is outside |eta|<" << jetsAcceptance_[0] << ": " << nJetsAcceptance[0];
  edm::LogProblem   ("SummaryError") << "Number of events where at least one jet has pt<" << jetsAcceptance_[1] << ": " << nJetsAcceptance[1];
  edm::LogProblem   ("SummaryError") << "Number of events where at least one muon is outside |eta|<" << muonsAcceptance_[0] << ": " << nMuonsAcceptance[0];
  edm::LogProblem   ("SummaryError") << "Number of events where at least one muon has pt<" << muonsAcceptance_[1] << ": " << nMuonsAcceptance[1];
 
  edm::LogVerbatim ("SummaryResults") << " ";
  edm::LogVerbatim ("SummaryResults") << " -------------------------------";
  edm::LogVerbatim ("SummaryResults") << "  Info on resolution of Muon/MET";
  edm::LogVerbatim ("SummaryResults") << " -------------------------------";
  edm::LogVerbatim ("SummaryResults") << " ";
 
  edm::LogVerbatim ("SummaryResults") << "The numbers and plots are made only for semi-muonic ttbar events";
  double a,b;
  a = TH1Fcontainer_["resEtMuons"]->GetMean();
  b = TH1Fcontainer_["resEtMuons"]->GetRMS();
  edm::LogVerbatim ("MainResults") << "MUON:Mean value of Et resolution: " << a*100 << "% +/- " << b*100 << "%";
  a = TH1Fcontainer_["resPtMuons"]->GetMean();
  b = TH1Fcontainer_["resPtMuons"]->GetRMS();
  edm::LogVerbatim ("MainResults") << "MUON:Mean value of Pt resolution: " << a*100 << "% +/- " << b*100 << "%";
  a = TH1Fcontainer_["resPhiMuons"]->GetMean();
  b = TH1Fcontainer_["resPhiMuons"]->GetRMS();
  edm::LogVerbatim ("MainResults") << "MUON:Mean value of Phi resolution: " << a*100 << "% +/- " << b*100 << "%";
  a = TH1Fcontainer_["resThetaMuons"]->GetMean();
  b = TH1Fcontainer_["resThetaMuons"]->GetRMS();
  edm::LogVerbatim ("MainResults") << "MUON:Mean value of Theta resolution: " << a*100 << "% +/- " << b*100 << "%";
  a = TH1Fcontainer_["resEtaMuons"]->GetMean();
  b = TH1Fcontainer_["resEtaMuons"]->GetRMS();
  edm::LogVerbatim ("MainResults") << "MUON:Mean value of Eta resolution: " << a*100 << "% +/- " << b*100 << "%";


  a = TH1Fcontainer_["resEtMets"]->GetMean();
  b = TH1Fcontainer_["resEtMets"]->GetRMS();
  edm::LogVerbatim ("MainResults") << "MET:Mean value of Et resolution: " << a*100 << "% +/- " << b*100 << "%";
  a = TH1Fcontainer_["resPtMets"]->GetMean();
  b = TH1Fcontainer_["resPtMets"]->GetRMS();
  edm::LogVerbatim ("MainResults") << "MET:Mean value of Pt resolution: " << a*100 << "% +/- " << b*100 << "%";
  a = TH1Fcontainer_["resPhiMets"]->GetMean();
  b = TH1Fcontainer_["resPhiMets"]->GetRMS();
  edm::LogVerbatim ("MainResults") << "MET:Mean value of Phi resolution: " << a*100 << "% +/- " << b*100 << "%";
  a = TH1Fcontainer_["resPxMets"]->GetMean();
  b = TH1Fcontainer_["resPxMets"]->GetRMS();
  edm::LogVerbatim ("MainResults") << "MET:Mean value of Px resolution: " << a*100 << "% +/- " << b*100 << "%";
  a = TH1Fcontainer_["resPyMets"]->GetMean();
  b = TH1Fcontainer_["resPyMets"]->GetRMS();
  edm::LogVerbatim ("MainResults") << "MET:Mean value of Py resolution: " << a*100 << "% +/- " << b*100 << "%";
 

}

//define this as a plug-in
DEFINE_FWK_MODULE(KinematicsChecker);
