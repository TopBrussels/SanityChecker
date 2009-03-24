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
// $Id: JetMetChecker.cc,v 1.9 2009/03/23 15:11:57 jmmaes Exp $
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
#include "TString.h"
#include <boost/lexical_cast.hpp>
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
  std::map<std::string,TH1D*> TH1DcontainerMatchedJets_[2];
  std::map<std::string,TH1F*> TH1Fcontainer_; // simple map to contain all TH1F.
  std::map<std::string,TH2F*> TH2Fcontainer_; // simple map to contain all TH2F.
  std::map<std::string,TH2D*> TH2Dcontainer_[2]; // simple map to contain all TH2F.
  std::map<std::string,TH1D*> TH1DcontainerForbTagging_[11]; // simple map to contain all TH1D.
  //std::map<std::string,TGraph*> TGraphcontainerForbTagging_[11]; // simple map to contain all TGraph.

  std::string objectNames_[5];
  std::string bTaggerNames_[11]; 
  double lowerRanges_[11];
  double upperRanges_[11];
  int nBins;
  std::string PatJetsNames_[2]; 


  std::vector< std::string > availableTaggers;

};

//
// constants, enums and typedefs
//

//
// static data member definitions
//
std::string UncorrName[4] ={"uncorrALL","uncorrJES","uncorrMUON","uncorrMAXN"};
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

  PatJetsNames_[0] = "PatJetsAll";
  PatJetsNames_[1] = "MatchedTopJets";
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
    return; //  throw cms::Exception("ProductNotFound") <<"jet collection not found"<<std::endl;
  }
  if (!mets.isValid()){
    edm::LogWarning  ("NoMetsFound") << "JetMetCheckerWarning - NoMetsFound";
    return;// throw cms::Exception("ProductNotFound") <<"MET collection not found"<<std::endl;
  }
  
  vector <CaloTowerPtr> jettowers;
  double sumtrkpt;
  double sumtwrpt; 
  for( unsigned int i=0;i<jets->size();i++) {
 

    TH1DcontainerMatchedJets_[0]["Jetn90"]->Fill((*jets)[i].n90());
    TH1DcontainerMatchedJets_[0]["JetTowersArea"]->Fill((*jets)[i].towersArea());
    TH1DcontainerMatchedJets_[0]["emEnergyFraction"]->Fill((*jets)[i].emEnergyFraction());
    TH1DcontainerMatchedJets_[0]["energyFractionHadronic"]->Fill((*jets)[i].energyFractionHadronic());
    TH1DcontainerMatchedJets_[0]["maxEInEmTowers"]->Fill((*jets)[i].maxEInEmTowers());
    TH1DcontainerMatchedJets_[0]["maxEInHadTowers"]->Fill((*jets)[i].maxEInHadTowers());

   
    JetCorrName = (*jets)[i].jetCorrName();   
    
    //towers infos
    sumtwrpt   =0;
    jettowers = (*jets)[i].getCaloConstituents(); 
    std::vector <CaloTowerPtr>::const_iterator caloiter;
    for(caloiter=jettowers.begin();caloiter!=jettowers.end();caloiter++){
      //       double caloet=(*caloiter)->et();
      TH1Dcontainer_["JetTwrEt"]->Fill((*caloiter)->et());
      TH1Dcontainer_["JetTwrPhi"]->Fill((*caloiter)->phi());
      TH1Dcontainer_["JetTwrEta"]->Fill((*caloiter)->eta());
      sumtwrpt   +=(*caloiter)->pt();
    }
     TH1Dcontainer_["JetTwrSumPt"]->Fill( sumtwrpt);
     TH1Dcontainer_["JetdiffTwrSumPt"]->Fill( fabs(sumtwrpt- (*jets)[i].pt()));
 
     double diff =  jettowers.size() - (*jets)[i].n90();
     TH2Dcontainer_[0]["Jetn90vsE-n90"]->Fill((*jets)[i].n90(), diff  );
	

    //tracks infos
    sumtrkpt   =0;
    edm::RefVector< reco::TrackCollection > tracks;
    tracks = (*jets)[i].associatedTracks();
    for(unsigned int ti=0; ti<tracks.size(); ti++){
      TH1Dcontainer_["JetTrkPt"]->Fill(tracks[ti]->pt());
      sumtrkpt+=tracks[ti]->pt();
    }
    TH1Dcontainer_["JetTrkSumPt"]->Fill( sumtrkpt);
    TH1Dcontainer_["JetdiffTrkSumPt"]->Fill( fabs(sumtrkpt- (*jets)[i].pt()));
    
  }
  
   string Corrlevel[6] = { "REL","ABS","EMF","HAD", "UE", "PART" };
   string CorrFlav[4] = {"GLU","UDS","C","B" } ;
   string HistoName = "";
   
  for(std::vector<pat::Jet>::const_iterator jet_iter = jets->begin(); jet_iter!=jets->end(); ++jet_iter)
  {
       if(jet_iter->genJet() == 0) continue;
       TH2Fcontainer_["JetEtaResponse_UpToL2"]         ->Fill(jet_iter->genJet()->eta(),jet_iter->correctedJet(Corrlevel[0]).eta()/jet_iter->genJet()->eta());
       TH2Fcontainer_["JetPtResponse_UpToL3"]	       ->Fill(jet_iter->genJet()->pt(), jet_iter->correctedJet(Corrlevel[1]).pt()/jet_iter->genJet()->pt());
       TH2Fcontainer_["JetPtResponse_UpToL4"]	       ->Fill(jet_iter->genJet()->pt(), jet_iter->correctedJet(Corrlevel[2]).pt()/jet_iter->genJet()->pt());
       
       for(int ii = 5; ii < 8; ii++) //Loop over L5,L6 and L7
       {
               Int_t jj = 0;
               //std::cout<<"Parton flavour : "<<jet_iter->partonFlavour()<<std::endl;
               if(fabs(jet_iter->partonFlavour()) == 21)      jj = 0;
               else if(fabs(jet_iter->partonFlavour()) < 4)   jj = 1;
               else if(fabs(jet_iter->partonFlavour()) == 4)  jj = 2;
               else if(fabs(jet_iter->partonFlavour()) == 5)  jj = 3;
               else std::cout<<"Muchas problemas, capitan!"<<std::endl; //FIXME : Of course!
               
	       // the jet reponse is produced according to the parton flavour : L5b, L6b and L7b for b-jets for instance
	       HistoName  = CorrFlav[jj];
               HistoName += "_JetPtResponse_UpToL";
               HistoName += boost::lexical_cast<std::string>(ii);
               HistoName += "_";
               HistoName += CorrFlav[jj];
               
               TH2Fcontainer_[HistoName]->Fill(jet_iter->genJet()->pt(), jet_iter->correctedJet(Corrlevel[ii-2],CorrFlav[jj]).pt()/jet_iter->genJet()->pt());
       }
  }
  //close for jets
  
  if(mets->size()>0)
  {
  	TH1Dcontainer_["CaloMETInmHF"]->Fill((*mets)[0].CaloMETInmHF());
  	TH1Dcontainer_["CaloMETInpHF"]->Fill((*mets)[0].CaloMETInpHF());
  	TH1Dcontainer_["CaloMETPhiInmHF"]->Fill((*mets)[0].CaloMETPhiInmHF());
  	TH1Dcontainer_["CaloMETPhiInpHF"]->Fill((*mets)[0].CaloMETPhiInpHF());

	std::string histoname;
	Int_t Idx = 0;
	//uncorrALL      //! uncorrect to bare bones
	//uncorrJES,     //! uncorrect for JES only
	//uncorrMUON,    //! uncorrect for MUON only
	for(pat::MET::UncorrectionType i = (*mets)[0].uncorrALL ; i != (*mets)[0].uncorrMAXN ;i = pat::MET::UncorrectionType(i+1))
  	{
  		histoname = "CorEx_";histoname += UncorrName[Idx];
		TH1Dcontainer_[histoname] 	->Fill((*mets)[0].corEx(i));
 		histoname  = "CorEy_";histoname  += UncorrName[Idx];
  		TH1Dcontainer_[histoname]   	->Fill((*mets)[0].corEy(i));
 		histoname  = "CorSumEt_"; histoname  += UncorrName[Idx];
  		TH1Dcontainer_[histoname]	->Fill((*mets)[0].corSumEt(i));
 		histoname  = "uncorrectedPhi_"; histoname  += UncorrName[Idx];
  		TH1Dcontainer_[histoname]	->Fill((*mets)[0].uncorrectedPhi(i));
 		histoname  = "uncorrectedPt_"; histoname  += UncorrName[Idx];
  		TH1Dcontainer_[histoname]	->Fill((*mets)[0].uncorrectedPt(i));
		Idx++;
  	}
  	TH1Dcontainer_["emEtFraction"]		->Fill((*mets)[0].emEtFraction());
  	TH1Dcontainer_["emEtInEB"] 		->Fill((*mets)[0].emEtInEB());
  	TH1Dcontainer_["emEtInEE"] 		->Fill((*mets)[0].emEtInEE());
  	TH1Dcontainer_["emEtInHF"] 	   	->Fill((*mets)[0].emEtInHF());
  	TH1Dcontainer_["etFractionHadronic"]	->Fill((*mets)[0].etFractionHadronic());
  	TH1Dcontainer_["hadEtInHB"] 		->Fill((*mets)[0].hadEtInHB());
  	TH1Dcontainer_["hadEtInHE"] 		->Fill((*mets)[0].hadEtInHE());
  	TH1Dcontainer_["hadEtInHF"] 		->Fill((*mets)[0].hadEtInHF());
  	TH1Dcontainer_["hadEtInHO"] 		->Fill((*mets)[0].hadEtInHO());
  	TH1Dcontainer_["maxEtInEmTowers"]	->Fill((*mets)[0].maxEtInEmTowers());
  	TH1Dcontainer_["maxEtInHadTowers"]	->Fill((*mets)[0].maxEtInHadTowers());
  	TH1Dcontainer_["metSignificance"] 	->Fill((*mets)[0].metSignificance());

  	TH1Dcontainer_["genMET"] 		->Fill((*mets)[0].genMET()->et());
  	TH1Dcontainer_["genMET_eta"] 		->Fill((*mets)[0].genMET()->eta());
  	TH1Dcontainer_["genMET_phi"] 		->Fill((*mets)[0].genMET()->phi());
  	TH1Dcontainer_["genMET_px"] 		->Fill((*mets)[0].genMET()->px());
  	TH1Dcontainer_["genMET_py"] 		->Fill((*mets)[0].genMET()->py());
  	TH1Dcontainer_["genMET_pz"] 		->Fill((*mets)[0].genMET()->pz());
  	TH1Dcontainer_["genMET_emEnergy"] 	->Fill((*mets)[0].genMET()->emEnergy());  
  	TH1Dcontainer_["genMET_hadEnergy"] 	->Fill((*mets)[0].genMET()->hadEnergy());  
  	TH1Dcontainer_["genMET_auxiliaryEnergy"]->Fill((*mets)[0].genMET()->auxiliaryEnergy());
  	TH1Dcontainer_["genMET_mEtSig"] 	->Fill((*mets)[0].genMET()->mEtSig());
  }
  // close for MET
  
  Handle<TtGenEvent> genEvt;
  iEvent.getByLabel ("genEvt",genEvt);
  
  //Check if branch is available  
  if (!genEvt.isValid()){
    edm::LogWarning  ("NoGenEvtFound") << "JetMetCheckerWarning - NoGenEvtFound";
  }

  
  bool taggerIsAvailable=false;
  bool taggerAlreadyIn=false;
  std::vector<pat::Jet> jets_clone;

  
  if(  jets->size() >= 4){  
    if (genEvt.isValid()){
      if(genEvt->isSemiLeptonic(genEvt->kMuon)) {
	//make a copy of the jet collection to add the matched jets
	///////	  std::vector<pat::Jet> jets_clone;
	
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
	  JetPartonMatching *GenMatchTopQuarks = new JetPartonMatching(TopQuarks, *jets, matchingAlgo_, useMaxDist_, useDeltaR_, maxDist_);
	  for(unsigned int l=0; l<4; l++){
	    Int_t Idx = GenMatchTopQuarks->getMatchForParton(l,0);
	    if(Idx>=0){ jets_clone.push_back((*jets)[Idx]); }  
	  }
	}
      }
    }
  }

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
    //fill b-tag histo's only for jets fom top quark decay partons

    //fill b-tag histo's for each jet
    for(std::vector<pat::Jet>::const_iterator iJet = jets_clone.begin(); iJet != jets_clone.end(); ++iJet) {
      TH1DcontainerForbTagging_[i]["InclusiveTtSemiMu"]->Fill(iJet->bDiscriminator(bTaggerNames_[i]));
      
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

  vector <CaloTowerPtr> jettowers_match; 
  for(std::vector<pat::Jet>::const_iterator iJet = jets_clone.begin(); iJet != jets_clone.end(); ++iJet) {
    
    TH1DcontainerMatchedJets_[1]["Jetn90"]->Fill( iJet->n90());
    TH1DcontainerMatchedJets_[1]["JetTowersArea"]->Fill(iJet->towersArea());
    TH1DcontainerMatchedJets_[1]["emEnergyFraction"]->Fill(iJet->emEnergyFraction());
    TH1DcontainerMatchedJets_[1]["energyFractionHadronic"]->Fill(iJet->energyFractionHadronic());
    TH1DcontainerMatchedJets_[1]["maxEInEmTowers"]->Fill(iJet->maxEInEmTowers());
    TH1DcontainerMatchedJets_[1]["maxEInHadTowers"]->Fill(iJet->maxEInHadTowers());
    jettowers_match = iJet->getCaloConstituents(); 
    TH2Dcontainer_[1]["Jetn90vsE-n90"]->Fill(iJet->n90(), jettowers_match.size() - iJet->n90());
  }

}


// ------------ method called once each job just before starting event loop  ------------
void 
JetMetChecker::beginJob(const edm::EventSetup&)
{
  edm::Service<TFileService> fs;
  if (!fs) throw edm::Exception(edm::errors::Configuration, "TFileService missing from configuration!");
 
  std::vector< TFileDirectory > subDirsPatJets;
  for(unsigned int i=0;i<2;i++) subDirsPatJets.push_back(fs->mkdir( PatJetsNames_[i] ));
  
  TFileDirectory subDir = fs->mkdir( "PatJets" );
  // TFileDirectory subsubDir = subDir.mkdir( "PatJets" );
  
  
  //  TH1Dcontainer_["JetTwrEt"] = fs->make<TH1D>("JetTwrEt" ,"jet towers Et ",100,0,1000);// wrote in main directory
  for(unsigned int i=0;i<2;i++){
    TH1DcontainerMatchedJets_[i]["Jetn90"] = subDirsPatJets[i].make<TH1D>("Jetn90" ,"n90 ",10,0,10);
    TH1DcontainerMatchedJets_[i]["JetTowersArea"] = subDirsPatJets[i].make<TH1D>("JetTowersArea" ," Jet Towers Area ",nBins,0,1);
    TH1DcontainerMatchedJets_[i]["emEnergyFraction"] = subDirsPatJets[i].make<TH1D>("emEnergyFraction" ,"jet electomagnetic energy fraction ",nBins,0, 1);
    TH1DcontainerMatchedJets_[i]["energyFractionHadronic"]= subDirsPatJets[i].make<TH1D>("energyFractionHadronic" ,"jet hadronic energy fraction ",nBins,0, 1);
    TH1DcontainerMatchedJets_[i]["maxEInEmTowers"]= subDirsPatJets[i].make<TH1D>("maxEInEmTowers" ,"maximum energy deposited in ECAL towers ",nBins,0, 100);
    TH1DcontainerMatchedJets_[i]["maxEInHadTowers"]= subDirsPatJets[i].make<TH1D>("maxEInHadTowers" ,"maximum energy deposited in HCAL towers ",nBins,0, 100);
    TH2Dcontainer_[i]["Jetn90vsE-n90"]= subDirsPatJets[i].make<TH2D>("Jetn90vsE-n90" ,"n90 vs jet n-n90 ",10,0, 10, 10, 0, 10);
  }
  TH1Dcontainer_["JetTwrEt"] = subDir.make<TH1D>("JetTwrEt" ,"jet towers Et ",nBins,0,100);
  TH1Dcontainer_["JetTwrEta"] = subDir.make<TH1D>("JetTwrEta" ,"jet towers Eta ",nBins,-6, 6);
  TH1Dcontainer_["JetTwrPhi"] = subDir.make<TH1D>("JetTwrPhi" ,"jet towers Phi ",25,-3.2, 3.2);
  TH1Dcontainer_["JetTrkPt"] = subDir.make<TH1D>("JetTrkPt" ,"jet tracks Pt ",nBins,0, 100);
  TH1Dcontainer_["JetTwrSumPt"] = subDir.make<TH1D>("JetTwrSumPt" ,"Sum of Pt towers",nBins,0, 100);
  TH1Dcontainer_["JetdiffTwrSumPt"] = subDir.make<TH1D>("JetdiffTwrSumPt" ,"Diff between Sum of Pt towers and Pt Jet",nBins,0, 50);
  TH1Dcontainer_["JetTrkSumPt"] = subDir.make<TH1D>("JetTrkSumPt" ,"Sum of Pt tracks",nBins,0, 100);
  TH1Dcontainer_["JetdiffTrkSumPt"] = subDir.make<TH1D>("JetdiffTrkSumPt" ,"Diff between Sum of Pt tracks and Pt Jet",nBins,0, 50);
  
  TH2Fcontainer_["JetEtaResponse_UpToL2"]         = subDir.make<TH2F>("JetEtaResponse_UpToL2","",4000,-10,10,500,0,10);
  TH2Fcontainer_["JetPtResponse_UpToL3"]          = subDir.make<TH2F>("JetPtResponse_UpToL3","",4000,0,800,500,0,10);
  TH2Fcontainer_["JetPtResponse_UpToL4"]          = subDir.make<TH2F>("JetPtResponse_UpToL4","",4000,0,800,500,0,10);
  TH2Fcontainer_["GLU_JetPtResponse_UpToL5_GLU"]  = subDir.make<TH2F>("GLU_JetPtResponse_UpToL5_GLU","",4000,0,800,500,0,10);
  TH2Fcontainer_["GLU_JetPtResponse_UpToL6_GLU"]  = subDir.make<TH2F>("GLU_JetPtResponse_UpToL6_GLU","",4000,0,800,500,0,10);
  TH2Fcontainer_["GLU_JetPtResponse_UpToL7_GLU"]  = subDir.make<TH2F>("GLU_JetPtResponse_UpToL7_GLU","",4000,0,800,500,0,10);
  TH2Fcontainer_["UDS_JetPtResponse_UpToL5_UDS"]  = subDir.make<TH2F>("UDS_JetPtResponse_UpToL5_UDS","",4000,0,800,500,0,10);
  TH2Fcontainer_["UDS_JetPtResponse_UpToL6_UDS"]  = subDir.make<TH2F>("UDS_JetPtResponse_UpToL6_UDS","",4000,0,800,500,0,10);
  TH2Fcontainer_["UDS_JetPtResponse_UpToL7_UDS"]  = subDir.make<TH2F>("UDS_JetPtResponse_UpToL7_UDS","",4000,0,800,500,0,10);
  TH2Fcontainer_["C_JetPtResponse_UpToL5_C"]      = subDir.make<TH2F>("C_JetPtResponse_UpToL5_C","",4000,0,800,500,0,10);
  TH2Fcontainer_["C_JetPtResponse_UpToL6_C"]      = subDir.make<TH2F>("C_JetPtResponse_UpToL6_C","",4000,0,800,500,0,10);
  TH2Fcontainer_["C_JetPtResponse_UpToL7_C"]      = subDir.make<TH2F>("C_JetPtResponse_UpToL7_C","",4000,0,800,500,0,10);
  TH2Fcontainer_["B_JetPtResponse_UpToL5_B"]      = subDir.make<TH2F>("B_JetPtResponse_UpToL5_B","",4000,0,800,500,0,10);
  TH2Fcontainer_["B_JetPtResponse_UpToL6_B"]      = subDir.make<TH2F>("B_JetPtResponse_UpToL6_B","",4000,0,800,500,0,10);
  TH2Fcontainer_["B_JetPtResponse_UpToL7_B"]      = subDir.make<TH2F>("B_JetPtResponse_UpToL7_B","",4000,0,800,500,0,10);
  //TH2Fcontainer_[""]= subDir.make<TH2F>("","",XX,XX,XX,XX,XX,XX);
  //TH2Fcontainer_[""]= subDir.make<TH2F>("","",XX,XX,XX,XX,XX,XX);
  
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

  TFileDirectory metsubDir = fs->mkdir( "PatMets" );
  TH1Dcontainer_["CaloMETInmHF"] 		  = metsubDir.make<TH1D>("CaloMETInmHF","ME_{T} in the forward (-) hadronic calorimeter",800,0,400);
  TH1Dcontainer_["CaloMETInpHF"]		  = metsubDir.make<TH1D>("CaloMETInpHF","ME_{T} in the forward (+) hadronic calorimeter",800,0,400);
  TH1Dcontainer_["CaloMETPhiInmHF"]		  = metsubDir.make<TH1D>("CaloMETPhiInmHF","ME_{T} Phi in the forward (-) hadronic calorimeter",400,-4,4);
  TH1Dcontainer_["CaloMETPhiInpHF"]		  = metsubDir.make<TH1D>("CaloMETPhiInpHF","ME_{T} Phi in the forward (+) hadronic calorimeter",400,-4,4);
  std::string histoname, histotitle;
  for(unsigned int i=0;i<3;i++)
  {
 	histoname  = "CorEx_";                  histoname  += UncorrName[i];
 	histotitle = "Uncorrected ME_{x} : ";   histotitle += UncorrName[i];
  	TH1Dcontainer_[histoname]  = metsubDir.make<TH1D>(histoname.c_str(),histotitle.c_str(),400,0,200);
 	histoname  = "CorEy_";                  histoname  += UncorrName[i];
 	histotitle = "Uncorrected ME_{y} : ";   histotitle += UncorrName[i];
  	TH1Dcontainer_[histoname]  = metsubDir.make<TH1D>(histoname.c_str(),histotitle.c_str(),400,0,200);
 	histoname  = "CorSumEt_";               histoname  += UncorrName[i];
 	histotitle = "Uncorrected E_{T} sum : ";histotitle += UncorrName[i];
  	TH1Dcontainer_[histoname]  = metsubDir.make<TH1D>(histoname.c_str(),histotitle.c_str(),800,0,400);
 	histoname  = "uncorrectedPhi_";          histoname  += UncorrName[i];
 	histotitle = "Uncorrected ME_{T} phi : ";histotitle += UncorrName[i];
  	TH1Dcontainer_[histoname] = metsubDir.make<TH1D>(histoname.c_str(),histotitle.c_str(),400,-4,4);
 	histoname  = "uncorrectedPt_";          histoname  += UncorrName[i];
 	histotitle = "Uncorrected ME_{T} pt : ";histotitle += UncorrName[i];
  	TH1Dcontainer_[histoname] = metsubDir.make<TH1D>(histoname.c_str(),histotitle.c_str(),800,0,400);
  }

  TH1Dcontainer_["emEtFraction"] 	  = metsubDir.make<TH1D>("emEtFraction","Event electromagnetic energy fraction",100,0,1);
  TH1Dcontainer_["emEtInEB"] 		  = metsubDir.make<TH1D>("emEtInEB","Event electromagnetic energy in the ECAL barrel",800,0,400);
  TH1Dcontainer_["emEtInEE"] 		  = metsubDir.make<TH1D>("emEtInEE","Event electromagnetic energy in the ECAL end-cap",800,0,400);
  TH1Dcontainer_["emEtInHF"] 		  = metsubDir.make<TH1D>("emEtInHF","Event electromagnetic energy extracted from the forward HCAL",800,0,400);
  TH1Dcontainer_["etFractionHadronic"]    = metsubDir.make<TH1D>("etFractionHadronic","Event hadronic energy fraction",100,0,1);
  TH1Dcontainer_["hadEtInHB"] 		  = metsubDir.make<TH1D>("hadEtInHB","Event hadronic energy in the HCAL barrel",800,0,400);
  TH1Dcontainer_["hadEtInHE"] 		  = metsubDir.make<TH1D>("hadEtInHE","Event hadronic energy in the HCAL end-cap",800,0,400);
  TH1Dcontainer_["hadEtInHF"] 		  = metsubDir.make<TH1D>("hadEtInHF","Event hadronic energy in the forward HCAL",800,0,400);
  TH1Dcontainer_["hadEtInHO"] 		  = metsubDir.make<TH1D>("hadEtInHO","Event hadronic energy in the forward HCAL",800,0,400);
  TH1Dcontainer_["maxEtInEmTowers"] 	  = metsubDir.make<TH1D>("maxEtInEmTowers","Maximum energy deposited in ECAL towers",800,0,400);
  TH1Dcontainer_["maxEtInHadTowers"] 	  = metsubDir.make<TH1D>("maxEtInHadTowers","Maximum energy deposited in HCAL towers",800,0,400);
  TH1Dcontainer_["metSignificance"] 	  = metsubDir.make<TH1D>("metSignificance","Missing transverse energy significance (MET/sqrt(SumEt))",200,-20,20);

  TH1Dcontainer_["genMET"] 		  = metsubDir.make<TH1D>("genMET","Generated missing transverse energy",800,0,400);
  TH1Dcontainer_["genMET_eta"] 		  = metsubDir.make<TH1D>("genMET_eta","Generated missing transverse energy",500,-5,5);
  TH1Dcontainer_["genMET_phi"] 		  = metsubDir.make<TH1D>("genMET_phi","Generated missing  energy phi",400,-4,4);
  TH1Dcontainer_["genMET_px"] 		  = metsubDir.make<TH1D>("genMET_px","Generated missing energy px",800,0,400);
  TH1Dcontainer_["genMET_py"] 		  = metsubDir.make<TH1D>("genMET_py","Generated missing energy py",800,0,400);
  TH1Dcontainer_["genMET_pz"] 		  = metsubDir.make<TH1D>("genMET_pz","Generated missing energy pz",800,0,400);
  TH1Dcontainer_["genMET_emEnergy"] 	  = metsubDir.make<TH1D>("genMET_emEnergy","Energy of electromagnetic particles",800,0,400);
  TH1Dcontainer_["genMET_hadEnergy"] 	  = metsubDir.make<TH1D>("genMET_hadEnergy","Energy of hadronic particles",800,0,400);
  TH1Dcontainer_["genMET_auxiliaryEnergy"]= metsubDir.make<TH1D>("genMET_auxiliaryEnergy","Other energy (undecayed Sigmas etc.)",800,0,400);
  TH1Dcontainer_["genMET_mEtSig"] 	  = metsubDir.make<TH1D>("genMET_mEtSig","Generated missing transverse energy significance (MET/sqrt(SumEt))",200,-20,20);

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
