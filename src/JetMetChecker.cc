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
// $Id: JetMetChecker.cc,v 1.5 2009/03/09 15:01:48 jmmaes Exp $
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

#include "TopQuarkAnalysis/TopTools/interface/JetPartonMatching.h"
#include "AnalysisDataFormats/TopObjects/interface/TtGenEvent.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/MET.h"

#include "TDirectory.h"
#include "TH1D.h"
#include "TH1.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TGraph.h"
#include <string>
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
  int matchingAlgo_;
  bool useMaxDist_;
  bool useDeltaR_;
  double maxDist_;
  
  std::string  JetCorrName;
 
  //Histograms are booked in the beginJob() method
  std::map<std::string,TDirectory*> TDirectorycontainer_; // simple map to contain all TDirectory.
  std::map<std::string,TH1D*> TH1Dcontainer_; // simple map to contain all TH1D.
  std::map<std::string,TH1F*> TH1Fcontainer_; // simple map to contain all TH1F.
  std::map<std::string,TH2F*> TH2Fcontainer_; // simple map to contain all TH2F.
  std::map<std::string,TH1D*> TH1DcontainerForbTagging_[11]; // simple map to contain all TH1D.
  //std::map<std::string,TGraph*> TGraphcontainerForbTagging_[11]; // simple map to contain all TGraph.

  std::string objectNames_[5];
  std::string bTaggerNames_[11]; 
  double lowerRanges_[11];
  double upperRanges_[11];
  int nBins;

  std::vector< std::string > availableTaggers;

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
  matchingAlgo_     =   iConfig.getParameter<int>( "matchingAlgo" );
  useMaxDist_       =   iConfig.getParameter<bool>( "useMaxDist" );
  useDeltaR_        =   iConfig.getParameter<bool>( "useDeltaR" );
  maxDist_          =   iConfig.getParameter<double>( "maxDist" );

  objectNames_[0] = "Inclusive";
  objectNames_[1] = "bOnly";
  objectNames_[2] = "cOnly";
  objectNames_[3] = "gOnly";
  objectNames_[4] = "lOnly";
  
  bTaggerNames_[0] = "trackCountingHighPurBJetTags";
  bTaggerNames_[1] = "trackCountingHighEffBJetTags";
  bTaggerNames_[2] = "softMuonNoIPBJetTags";
  bTaggerNames_[3] = "softMuonBJetTags";
  bTaggerNames_[4] = "softElectronBJetTags";
  bTaggerNames_[5] = "simpleSecondaryVertexBJetTags";
  bTaggerNames_[6] = "jetProbabilityBJetTags";
  bTaggerNames_[7] = "jetBProbabilityBJetTags";
  bTaggerNames_[8] = "impactParameterMVABJetTags";
  bTaggerNames_[9] = "combinedSecondaryVertexMVABJetTags";
  bTaggerNames_[10] = "combinedSecondaryVertexBJetTags";
  
  lowerRanges_[0]=-10;
  lowerRanges_[1]=-10;
  lowerRanges_[2]=0;
  lowerRanges_[3]=0;
  lowerRanges_[4]=0;
  lowerRanges_[5]=0;
  lowerRanges_[6]=0;
  lowerRanges_[7]=0;
  lowerRanges_[8]=0;
  lowerRanges_[9]=0;
  lowerRanges_[10]=0;

  upperRanges_[0]=30;
  upperRanges_[1]=30;
  upperRanges_[2]=1;
  upperRanges_[3]=1;
  upperRanges_[4]=1;
  upperRanges_[5]=8;
  upperRanges_[6]=2.5;
  upperRanges_[7]=8;
  upperRanges_[8]=1;
  upperRanges_[9]=1;
  upperRanges_[10]=1;

  nBins=50;
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
    edm::LogWarning  ("NoJetsFound") << "JetMetCheckerWarning - NoJetsFound"; 
    throw cms::Exception("ProductNotFound") <<"jet collection not found"<<std::endl;
  }
  if (!mets.isValid()){
    edm::LogWarning  ("NoMetsFound") << "JetMetCheckerWarning - NoMetsFound";
    throw cms::Exception("ProductNotFound") <<"MET collection not found"<<std::endl;
  }
  
  vector <CaloTowerPtr> jettowers;
  
  
  for( unsigned int i=0;i<jets->size();i++) {
    
    TH1Dcontainer_["Jetn90"]->Fill((*jets)[i].n90());
    TH1Dcontainer_["JetTowersArea"]->Fill((*jets)[i].towersArea());
    TH1Dcontainer_["emEnergyFraction"]->Fill((*jets)[i].emEnergyFraction());
    TH1Dcontainer_["energyFractionHadronic"]->Fill((*jets)[i].energyFractionHadronic());
    TH1Dcontainer_["maxEInEmTowers"]->Fill((*jets)[i].maxEInEmTowers());
    TH1Dcontainer_["maxEInHadTowers"]->Fill((*jets)[i].maxEInHadTowers());
    
    JetCorrName = (*jets)[i].jetCorrName();   
    
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
  
  
  Handle<TtGenEvent> genEvt;
  iEvent.getByLabel ("genEvt",genEvt);
  
  //Check if branch is available  
  if (!genEvt.isValid()){
    edm::LogWarning  ("NoGenEvtFound") << "JetMetCheckerWarning - NoGenEvtFound";
  }
  
  bool taggerIsAvailable=false;
  bool taggerAlreadyIn=false;
  
  for(unsigned int i=0;i<11;i++){
    
    if(  jets->size() == 0) continue;
    
    //test if b-tagger is available
    for(unsigned int k=0; k<(*jets)[0].getPairDiscri().size(); k++){     
      if( (*jets)[0].getPairDiscri()[k].first != bTaggerNames_[i]) taggerIsAvailable = true;
      if(taggerIsAvailable==true){
	for(unsigned int l=0; l<availableTaggers.size(); l++){
	  if(availableTaggers[l]==bTaggerNames_[i]) taggerAlreadyIn=true;
	}
	if(!taggerAlreadyIn) availableTaggers.push_back(bTaggerNames_[i]);
	taggerAlreadyIn=false;
      }
    }
    if(taggerIsAvailable = false) continue;
    
    //fill b-tag histo's for each jet
    for(std::vector<pat::Jet>::const_iterator iJet = jets->begin(); iJet != jets->end(); ++iJet) {
      
      TH1DcontainerForbTagging_[i]["Inclusive"]->Fill(iJet->bDiscriminator(bTaggerNames_[i]));
      
      if (genEvt.isValid()){
	//std::cout << "iJet->partonFlavour() : " << iJet->partonFlavour() << std::endl;
	if(fabs(iJet->partonFlavour())==5){
	  TH1DcontainerForbTagging_[i]["bjets"]->Fill(iJet->bDiscriminator(bTaggerNames_[i]));
	}
	if(fabs(iJet->partonFlavour())==4){
	  TH1DcontainerForbTagging_[i]["cjets"]->Fill(iJet->bDiscriminator(bTaggerNames_[i]));
	}
	if(fabs(iJet->partonFlavour())==21){
	  TH1DcontainerForbTagging_[i]["gjets"]->Fill(iJet->bDiscriminator(bTaggerNames_[i]));
	} 
	if(fabs(iJet->partonFlavour())==1||fabs(iJet->partonFlavour())==2||fabs(iJet->partonFlavour())==3){
	  TH1DcontainerForbTagging_[i]["ljets"]->Fill(iJet->bDiscriminator(bTaggerNames_[i]));
	}
	if(fabs(iJet->partonFlavour())!=1&&fabs(iJet->partonFlavour())!=2&&fabs(iJet->partonFlavour())!=3&&fabs(iJet->partonFlavour())!=4&&fabs(iJet->partonFlavour())!=5&&fabs(iJet->partonFlavour())!=21){
	  TH1DcontainerForbTagging_[i]["ijets"]->Fill(iJet->bDiscriminator(bTaggerNames_[i]));
	  //std::cout << "partonflavour is not a udscb or g" << std::endl;
	}
      }
    }
    //fill b-tag histo's only for jets fom top quark decay partons
    if(  jets->size() >= 4){  
      if (genEvt.isValid()){
	if(genEvt->isSemiLeptonic(genEvt->kMuon)) {
	  //make a copy of the jet collection to add the matched jets
	  std::vector<pat::Jet> jets_clone;
	  
	  // Make plots only with jets from ttbar events (matched with a certain algo)
	  // Matching index : Hadronic Q  = 0, Hadronic Q' = 1, Hadronic b  = 2, Leptonic b  = 3;
	  std::vector<const reco::Candidate *> TopQuarks;
	  TopQuarks.push_back(genEvt->hadronicDecayQuark());
	  TopQuarks.push_back(genEvt->hadronicDecayQuarkBar());
	  TopQuarks.push_back(genEvt->hadronicDecayB());
	  TopQuarks.push_back(genEvt->leptonicDecayB());
	  if(TopQuarks.size()==4) { 
	    JetPartonMatching *GenMatchTopQuarks = new JetPartonMatching(TopQuarks, *jets, matchingAlgo_, useMaxDist_, useDeltaR_, maxDist_);
	    for(unsigned int l=0; l<4; l++){
	      Int_t Idx = GenMatchTopQuarks->getMatchForParton(l,0);
	      if(Idx>=0){ jets_clone.push_back((*jets)[Idx]); }  
	    }
	    //fill b-tag histo's for each jet
	    for(std::vector<pat::Jet>::const_iterator iJet = jets_clone.begin(); iJet != jets_clone.end(); ++iJet) {
	      TH1DcontainerForbTagging_[i]["InclusiveTtSemiMu"]->Fill(iJet->bDiscriminator(bTaggerNames_[i]));
	      
	      if (genEvt.isValid()){
		//std::cout << "iJet->partonFlavour() : " << iJet->partonFlavour() << std::endl;
		if(fabs(iJet->partonFlavour())==5){
		  TH1DcontainerForbTagging_[i]["bjetsTtSemiMu"]->Fill(iJet->bDiscriminator(bTaggerNames_[i]));
		}
		if(fabs(iJet->partonFlavour())==4){
		  TH1DcontainerForbTagging_[i]["cjetsTtSemiMu"]->Fill(iJet->bDiscriminator(bTaggerNames_[i]));
		}
		if(fabs(iJet->partonFlavour())==21){
		  TH1DcontainerForbTagging_[i]["gjetsTtSemiMu"]->Fill(iJet->bDiscriminator(bTaggerNames_[i]));
		} 
		if(fabs(iJet->partonFlavour())==1||fabs(iJet->partonFlavour())==2||fabs(iJet->partonFlavour())==3){
		  TH1DcontainerForbTagging_[i]["ljetsTtSemiMu"]->Fill(iJet->bDiscriminator(bTaggerNames_[i]));
		}
		if(fabs(iJet->partonFlavour())!=1&&fabs(iJet->partonFlavour())!=2&&fabs(iJet->partonFlavour())!=3&&fabs(iJet->partonFlavour())!=4&&fabs(iJet->partonFlavour())!=5&&fabs(iJet->partonFlavour())!=21){
		   TH1DcontainerForbTagging_[i]["ijetsTtSemiMu"]->Fill(iJet->bDiscriminator(bTaggerNames_[i]));
		  //std::cout << "partonflavour is not a udscb or g" << std::endl;
		}
	      }
	    }
	  }
	}
      }
    }
  }
}


// ------------ method called once each job just before starting event loop  ------------
void 
JetMetChecker::beginJob(const edm::EventSetup&)
{
  edm::Service<TFileService> fs;
  if (!fs) throw edm::Exception(edm::errors::Configuration, "TFileService missing from configuration!");
  
  TFileDirectory subDir = fs->mkdir( "PatJets" );
  // TFileDirectory subsubDir = subDir.mkdir( "PatJets" );
  
  
  //  TH1Dcontainer_["JetTwrEt"] = fs->make<TH1D>("JetTwrEt" ,"jet towers Et ",100,0,1000);// wrote in main directory
  TH1Dcontainer_["Jetn90"] = subDir.make<TH1D>("Jetn90" ,"n90 ",10,0,10);
  TH1Dcontainer_["JetTowersArea"] = subDir.make<TH1D>("JetTowersArea" ," Jet Towers Area ",nBins,0,1);
  TH1Dcontainer_["JetTwrEt"] = subDir.make<TH1D>("JetTwrEt" ,"jet towers Et ",nBins,0,100);
  TH1Dcontainer_["JetTwrEta"] = subDir.make<TH1D>("JetTwrEta" ,"jet towers Eta ",nBins,-6, 6);
  TH1Dcontainer_["JetTwrPhi"] = subDir.make<TH1D>("JetTwrPhi" ,"jet towers Phi ",nBins,-3.2, 3.2);
  TH1Dcontainer_["JetTrkPt"] = subDir.make<TH1D>("JetTrkPt" ,"jet tracks Pt ",nBins,0, 100);
  TH1Dcontainer_["emEnergyFraction"] = subDir.make<TH1D>("emEnergyFraction" ,"jet electomagnetic energy fraction ",nBins,0, 100);
  TH1Dcontainer_["energyFractionHadronic"]= subDir.make<TH1D>("energyFractionHadronic" ,"jet hadronic energy fraction ",nBins,0, 100);
  TH1Dcontainer_["maxEInEmTowers"]= subDir.make<TH1D>("maxEInEmTowers" ,"maximum energy deposited in ECAL towers ",nBins,0, 100);
  TH1Dcontainer_["maxEInHadTowers"]= subDir.make<TH1D>("maxEInHadTowers" ,"maximum energy deposited in HCAL towers ",nBins,0, 100);
  
  
  std::vector< TFileDirectory > subDirsbTagging;
  for(unsigned int i=0;i<11;i++) subDirsbTagging.push_back(fs->mkdir( bTaggerNames_[i] ));
  
  //for more on b-tagging validation: http://nippon.fnal.gov:8888/lpc1/cmsroc/yumiceva/validation/ 
  for(unsigned int i=0;i<11;i++){
    TH1DcontainerForbTagging_[i]["Inclusive"] = subDirsbTagging[i].make<TH1D>("Inclusive" ,"distribution of b discriminant",nBins,lowerRanges_[i],upperRanges_[i]);
    TH1DcontainerForbTagging_[i]["bjets"] = subDirsbTagging[i].make<TH1D>("bjets" ,"distribution of b discriminant",nBins,lowerRanges_[i],upperRanges_[i]);
    TH1DcontainerForbTagging_[i]["cjets"] = subDirsbTagging[i].make<TH1D>("cjets" ,"distribution of b discriminant",nBins,lowerRanges_[i],upperRanges_[i]);
    TH1DcontainerForbTagging_[i]["ljets"] = subDirsbTagging[i].make<TH1D>("ljets" ,"distribution of b discriminant",nBins,lowerRanges_[i],upperRanges_[i]);
    TH1DcontainerForbTagging_[i]["gjets"] = subDirsbTagging[i].make<TH1D>("gjets" ,"distribution of b discriminant",nBins,lowerRanges_[i],upperRanges_[i]);
    TH1DcontainerForbTagging_[i]["ijets"] = subDirsbTagging[i].make<TH1D>("ijets" ,"distribution of b discriminant",nBins,lowerRanges_[i],upperRanges_[i]);

    TH1DcontainerForbTagging_[i]["Effcjets"] = subDirsbTagging[i].make<TH1D>("Effcjets" ,"b tag efficiency versus c mistag rate",1000,0,1);
    TH1DcontainerForbTagging_[i]["Effljets"] = subDirsbTagging[i].make<TH1D>("Effljets" ,"b tag efficiency versus c mistag rate",1000,0,1);
    TH1DcontainerForbTagging_[i]["Effgjets"] = subDirsbTagging[i].make<TH1D>("Effgjets" ,"b tag efficiency versus c mistag rate",1000,0,1);
    TH1DcontainerForbTagging_[i]["Effijets"] = subDirsbTagging[i].make<TH1D>("Effijets" ,"b tag efficiency versus c mistag rate",1000,0,1);
  }
  
  for(unsigned int i=0;i<11;i++){
    TH1DcontainerForbTagging_[i]["InclusiveTtSemiMu"] = subDirsbTagging[i].make<TH1D>("InclusiveTtSemiMu" ,"distribution of b discriminant",nBins,lowerRanges_[i],upperRanges_[i]);
    TH1DcontainerForbTagging_[i]["bjetsTtSemiMu"] = subDirsbTagging[i].make<TH1D>("bjetsTtSemiMu" ,"distribution of b discriminant (TtSemiMu)",nBins,lowerRanges_[i],upperRanges_[i]);
    TH1DcontainerForbTagging_[i]["cjetsTtSemiMu"] = subDirsbTagging[i].make<TH1D>("cjetsTtSemiMu" ,"distribution of b discriminant (TtSemiMu)",nBins,lowerRanges_[i],upperRanges_[i]);
    TH1DcontainerForbTagging_[i]["ljetsTtSemiMu"] = subDirsbTagging[i].make<TH1D>("ljetsTtSemiMu" ,"distribution of b discriminant (TtSemiMu)",nBins,lowerRanges_[i],upperRanges_[i]);
    TH1DcontainerForbTagging_[i]["gjetsTtSemiMu"] = subDirsbTagging[i].make<TH1D>("gjetsTtSemiMu" ,"distribution of b discriminant (TtSemiMu)",nBins,lowerRanges_[i],upperRanges_[i]);
    TH1DcontainerForbTagging_[i]["ijetsTtSemiMu"] = subDirsbTagging[i].make<TH1D>("ijetsTtSemiMu" ,"distribution of b discriminant (TtSemiMu)",nBins,lowerRanges_[i],upperRanges_[i]);

    TH1DcontainerForbTagging_[i]["EffcjetsTtSemiMu"] = subDirsbTagging[i].make<TH1D>("EffcjetsTtSemiMu" ,"b tag efficiency versus c mistag rate (TtSemiMu)",1000,0,1);
    TH1DcontainerForbTagging_[i]["EffljetsTtSemiMu"] = subDirsbTagging[i].make<TH1D>("EffljetsTtSemiMu" ,"b tag efficiency versus c mistag rate (TtSemiMu)",1000,0,1);
    TH1DcontainerForbTagging_[i]["EffgjetsTtSemiMu"] = subDirsbTagging[i].make<TH1D>("EffgjetsTtSemiMu" ,"b tag efficiency versus c mistag rate (TtSemiMu)",1000,0,1);
    TH1DcontainerForbTagging_[i]["EffijetsTtSemiMu"] = subDirsbTagging[i].make<TH1D>("EffijetsTtSemiMu" ,"b tag efficiency versus c mistag rate (TtSemiMu)",1000,0,1);
  }

}

// ------------ method called once each job just after ending the event loop  ------------
void 
JetMetChecker::endJob() {

  //Piece of code to make b-tag efficiency as function of non-b tag efficiency

  //std::cout << "start endjob" << std::endl;

  double xVal=0;
  int intxVal=0;
  double xValAdder=0;

  double cVal=0;
  double cValAdder=0;

  double lVal=0;
  double lValAdder=0;

  double gVal=0;
  double gValAdder=0;
    
  double iVal=0;
  double iValAdder=0;
  
  for(unsigned int i=0;i<11;i++){
    
    xVal=0;
    xValAdder=0;
    intxVal=0;
    cVal=0;
    cValAdder=0;
    lVal=0;
    lValAdder=0;
    gVal=0;
    gValAdder=0;
    iVal=0;
    iValAdder=0;
    
    for(int bin=nBins+1; bin>0; bin--){//the number of bins should be the same, ALLWAYS!
      xValAdder += TH1DcontainerForbTagging_[i]["bjets"]->GetBinContent(bin);
      if(TH1DcontainerForbTagging_[i]["bjets"]->Integral()!=0) xVal=xValAdder/TH1DcontainerForbTagging_[i]["bjets"]->Integral();
      intxVal=int (1000*xVal);

      cValAdder += TH1DcontainerForbTagging_[i]["cjets"]->GetBinContent(bin);
      if(TH1DcontainerForbTagging_[i]["cjets"]->Integral()!=0) cVal=cValAdder/TH1DcontainerForbTagging_[i]["cjets"]->Integral();
      TH1DcontainerForbTagging_[i]["Effcjets"]->SetBinContent(intxVal,cVal);

      lValAdder += TH1DcontainerForbTagging_[i]["ljets"]->GetBinContent(bin);
      if(TH1DcontainerForbTagging_[i]["ljets"]->Integral()!=0) lVal=lValAdder/TH1DcontainerForbTagging_[i]["ljets"]->Integral();
      TH1DcontainerForbTagging_[i]["Effljets"]->SetBinContent(intxVal,lVal);
  
      gValAdder += TH1DcontainerForbTagging_[i]["gjets"]->GetBinContent(bin);
      if(TH1DcontainerForbTagging_[i]["gjets"]->Integral()!=0) gVal=gValAdder/TH1DcontainerForbTagging_[i]["gjets"]->Integral();
      TH1DcontainerForbTagging_[i]["Effgjets"]->SetBinContent(intxVal,gVal);
 
      iValAdder += TH1DcontainerForbTagging_[i]["ijets"]->GetBinContent(bin);
      if(TH1DcontainerForbTagging_[i]["ijets"]->Integral()!=0) iVal=iValAdder/TH1DcontainerForbTagging_[i]["ijets"]->Integral();
      TH1DcontainerForbTagging_[i]["Effijets"]->SetBinContent(intxVal,iVal);
    }  
  
    xVal=0;
    xValAdder=0;
    intxVal=0;
    cVal=0;
    cValAdder=0;
    lVal=0;
    lValAdder=0;
    gVal=0;
    gValAdder=0;
    iVal=0;
    iValAdder=0;
    
    for(int bin=nBins+1; bin>0; bin--){//the number of bins should be the same, ALLWAYS!
      xValAdder += TH1DcontainerForbTagging_[i]["bjetsTtSemiMu"]->GetBinContent(bin);
      if(TH1DcontainerForbTagging_[i]["bjetsTtSemiMu"]->Integral()!=0) xVal=xValAdder/TH1DcontainerForbTagging_[i]["bjetsTtSemiMu"]->Integral();
      intxVal=int (1000*xVal);

      cValAdder += TH1DcontainerForbTagging_[i]["cjetsTtSemiMu"]->GetBinContent(bin);
      if(TH1DcontainerForbTagging_[i]["cjetsTtSemiMu"]->Integral()!=0) cVal=cValAdder/TH1DcontainerForbTagging_[i]["cjetsTtSemiMu"]->Integral();
      TH1DcontainerForbTagging_[i]["EffcjetsTtSemiMu"]->SetBinContent(intxVal,cVal);

      lValAdder += TH1DcontainerForbTagging_[i]["ljetsTtSemiMu"]->GetBinContent(bin);
      if(TH1DcontainerForbTagging_[i]["ljetsTtSemiMu"]->Integral()!=0) lVal=lValAdder/TH1DcontainerForbTagging_[i]["ljetsTtSemiMu"]->Integral();
      TH1DcontainerForbTagging_[i]["EffljetsTtSemiMu"]->SetBinContent(intxVal,lVal);
  
      gValAdder += TH1DcontainerForbTagging_[i]["gjetsTtSemiMu"]->GetBinContent(bin);
      if(TH1DcontainerForbTagging_[i]["gjetsTtSemiMu"]->Integral()!=0) gVal=gValAdder/TH1DcontainerForbTagging_[i]["gjetsTtSemiMu"]->Integral();
      TH1DcontainerForbTagging_[i]["EffgjetsTtSemiMu"]->SetBinContent(intxVal,gVal);

      iValAdder += TH1DcontainerForbTagging_[i]["ijetsTtSemiMu"]->GetBinContent(bin);
      if(TH1DcontainerForbTagging_[i]["ijetsTtSemiMu"]->Integral()!=0) iVal=iValAdder/TH1DcontainerForbTagging_[i]["ijetsTtSemiMu"]->Integral();
      TH1DcontainerForbTagging_[i]["EffijetsTtSemiMu"]->SetBinContent(intxVal,iVal);

    }
  }  

  edm::LogVerbatim ("MainResults") << " -------------------------------------------";
  edm::LogVerbatim ("MainResults") << " -------------------------------------------";
  edm::LogVerbatim ("MainResults") << " --     Report from JetMMET Checker      -- ";
  edm::LogVerbatim ("MainResults") << " -------------------------------------------";
  edm::LogVerbatim ("MainResults") << " -------------------------------------------";

  edm::LogVerbatim ("MainResults") << " ";
  edm::LogVerbatim ("MainResults") << " -------------------------------";
  edm::LogVerbatim ("MainResults") << "  Info from jets";
  edm::LogVerbatim ("MainResults") << " -------------------------------";
  edm::LogVerbatim ("MainResults") << " ";

  edm::LogVerbatim ("MainResults") << "Current level of jet energy corrections = " << JetCorrName ; 

  edm::LogVerbatim ("MainResults") << " ";
  
  edm::LogVerbatim ("MainResults") << " -------------------------------";
  edm::LogVerbatim ("MainResults") << "  Info from b-tagging";
  edm::LogVerbatim ("MainResults") << " -------------------------------";
  edm::LogVerbatim ("MainResults") << " ";
  edm::LogVerbatim ("MainResults") << "Available b-tag algorithms:";

  for(unsigned int l=0; l<availableTaggers.size(); l++){
    edm::LogVerbatim ("MainResults") <<  availableTaggers[l];
  } 
  edm::LogVerbatim ("MainResults") << " ";
  //edm::LogVerbatim ("MainResults") << "b-tag efficiency vs. mistag rate at certain thresholds ";
 
  /*for(unsigned int i=0;i<11;i++){
    edm::LogVerbatim ("MainResults") << bTaggerNames_[i] << ": " <<  TH1DcontainerForbTagging_[i]["Effbjets"]->GetBinContent(500) << " " << TH1DcontainerForbTagging_[i]["Effcjets"]->GetBinContent(500) << " " << TH1DcontainerForbTagging_[i]["Effljets"]->GetBinContent(500) << " " << TH1DcontainerForbTagging_[i]["Effgjets"]->GetBinContent(500) ;
    }*/


}

//define this as a plug-in
DEFINE_FWK_MODULE(JetMetChecker);
