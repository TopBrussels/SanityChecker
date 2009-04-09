// -*- C++ -*-
//
// Package:    TtGenEventChecker
// Class:      TtGenEventChecker
// 
/**\class TtGenEventChecker TtGenEventChecker.cc UserCode/TtGenEventChecker/src/TtGenEventChecker.cc

 Description: <one line class summary>

 Implementation:
     <Notes on implementation>
*/
//
// Original Author:  local user
//         Created:  Wed Feb 18 16:39:03 CET 2009
// $Id: TtGenEvenChecker.cc,v 1.9 2009/04/08 12:34:47 echabert Exp $
//
//


// system include files
#include <memory>
#include <Math/VectorUtil.h>

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

#include "DataFormats/Candidate/interface/Particle.h"
#include "DataFormats/Math/interface/Vector3D.h"


#include "AnalysisDataFormats/TopObjects/interface/TtGenEvent.h"

#include "TDirectory.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TF2.h"

struct LowestPt{
  bool operator()( reco::Particle::LorentzVector j1, reco::Particle::LorentzVector j2 ) const{
    return j1.pt() < j2.pt() ;
  }
};
struct HighestPt{
  bool operator()( reco::Particle::LorentzVector j1, reco::Particle::LorentzVector j2 ) const{
    return j1.pt() > j2.pt() ;
  }
};
struct HighestEta{
  bool operator()( reco::Particle::LorentzVector j1, reco::Particle::LorentzVector j2 ) const{
    return fabs(j1.eta()) > fabs(j2.eta()) ;
  }
};
struct Lowest{
  bool operator()( float j1, float j2 ) const{
    return j1 < j2 ;
  }
};
struct Highest{
  bool operator()( float j1, float j2 ) const{
    return j1 > j2 ;
  }
};



//
// class decleration
//

class TtGenEventChecker:public
  edm::EDAnalyzer
{
public:
  explicit
  TtGenEventChecker (const edm::ParameterSet &);
   ~
  TtGenEventChecker ();


private:
  virtual void
  beginJob (const edm::EventSetup &);
  virtual void
  analyze (const edm::Event &, const edm::EventSetup &);
  virtual void
  endJob ();

  // ----------member data ---------------------------
  edm::InputTag genEvtCollectionName_;
  bool doit_all;
  bool doit_semileptonic;
  bool doit_semileptonicMuon;

  //Histograms are booked in the beginJob() method
  std::map < std::string, TDirectory * > TDirectorycontainer_;	// simple map to contain all TDirectory.
  std::map < std::string, TH1D * > TH1Dcontainer_;		// simple map to contain all TH1D.
  std::map < std::string, TH2D * > TH2Dcontainer_;		// simple map to contain all TH2D.
  std::map < std::string, TH1F * > TH1Fcontainer_;		// simple map to contain all TH1F.
  std::map < std::string, TH2F * > TH2Fcontainer_;		// simple map to contain all TH2F.
  TF2 fit2LB_;
  TF2 fit2LQ_;
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
TtGenEventChecker::TtGenEventChecker (const edm::ParameterSet & iConfig)
{
  //now do what ever initialization is needed
  genEvtCollectionName_ = iConfig.getParameter < edm::InputTag > ("genEvtCollectionName");
  doit_all = iConfig.getUntrackedParameter <bool>("doit_all",true);
  doit_semileptonic = iConfig.getUntrackedParameter <bool>("doit_semileptonic",false);
  doit_semileptonicMuon = iConfig.getUntrackedParameter <bool>("doit_semileptonicMuon",false);

}


TtGenEventChecker::~TtGenEventChecker ()
{

  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called to for each event  ------------
void
TtGenEventChecker::analyze (const edm::Event & iEvent, const edm::EventSetup & iSetup)
{
  using namespace edm;
  using namespace std;
  
  
  
#ifdef THIS_IS_AN_EVENT_EXAMPLE
  Handle < ExampleData > pIn;
  iEvent.getByLabel ("example", pIn);
#endif
  
#ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
  ESHandle < SetupData > pSetup;
  iSetup.get < SetupRecord > ().get (pSetup);
#endif
  //Here you handle the collection you want to access
  Handle < TtGenEvent > TtGenEvent_;
  iEvent.getByLabel (genEvtCollectionName_, TtGenEvent_);
  TtGenEvent genEvt = *TtGenEvent_;

 
 
  //tree levels of message
  //no endl needed 
  //use '\n' to go to next line
  //we can use different category for the same EDAnalyser
  //ex: NoDataFound - LinkBroken - TooMuchDataFound - SummaryError - MainResults
  
  //edm::LogError  ("category") << "My error message";    // or  edm::LogProblem  (not formated)
  //edm::LogWarning  ("category") << "My warning message"; // or  edm::LogPrint    (not formated)
  //edm::LogVerbatim   ("category") << "My LogVerbatim message";  // or  edm::LogVerbatim (not formated)
  
  //use Warning for event by event problem 
  //edm::LogWarning ("NoDataFound") << "My warning message - NoDataFound";	// or  edm::LogPrint    (not formated)
  //edm::LogWarning ("LinkBroken") << "My warning message - LinkBroken";	// or  edm::LogPrint    (not formated)
  
  bool doit = false;
  if(doit_all) doit=true; 
  if(doit_semileptonic) doit=true;
  if(doit_semileptonicMuon) doit=true;

  // Do only if requirement are fulfilled
  
  TH1Dcontainer_["BRatio"]->Fill(0);
  if(genEvt.isTtBar()){
    //test BR
    TH1Dcontainer_["BRatio"]->Fill(1);
    if(genEvt.isFullHadronic()) TH1Dcontainer_["BRatio"]->Fill(2);
    if(genEvt.isFullLeptonic()) TH1Dcontainer_["BRatio"]->Fill(3);
    if(genEvt.isSemiLeptonic()) TH1Dcontainer_["BRatio"]->Fill(4);
    if(genEvt.isSemiLeptonic(WDecay::kMuon)) TH1Dcontainer_["BRatio"]->Fill(5);
    if(genEvt.isSemiLeptonic(WDecay::kElec)) TH1Dcontainer_["BRatio"]->Fill(6);
    if(genEvt.isSemiLeptonic(WDecay::kTau)) TH1Dcontainer_["BRatio"]->Fill(7);
    if(genEvt.isSemiLeptonic() && genEvt.isSemiLeptonic(WDecay::kNone)) TH1Dcontainer_["BRatio"]->Fill(8);
    
    
    if(!genEvt.top()) edm::LogWarning ("NoDataFound_top") << "No access to top in ttbar event" ;
    if(!genEvt.topBar()) edm::LogWarning ("NoDataFound_top") << "No access to topBar in ttbar event" ;
    if(!genEvt.b()) edm::LogWarning ("NoDataFound_b") << "No access to b-quark in ttbar event" ;
    if(!genEvt.bBar()) edm::LogWarning ("NoDataFound_b") << "No access to bBar-quark in ttbar event" ;
    if(!genEvt.wMinus()) edm::LogWarning ("NoDataFound_w") << "No access to W- in ttbar event" ;
    if(!genEvt.wPlus()) edm::LogWarning ("NoDataFound_w") << "No access to W+ in ttbar event" ;
    
  }

  //embedded main part of the code !
  if(doit){
  
  
 
  TH1Dcontainer_["NofTopsRadiation"]->Fill (genEvt.radiatedGluons(-6).size()+genEvt.radiatedGluons(6).size());
  TH1Dcontainer_["NofISR"]->Fill (genEvt.topSisters ().size());
  TH1Dcontainer_["NofLepton"]->Fill(genEvt.numberOfLeptons(false));
  TH1Dcontainer_["NofLeptonFromW"]->Fill(genEvt.numberOfLeptons(true));
  TH1Dcontainer_["NofB"]->Fill(genEvt.numberOfBQuarks(false));
  TH1Dcontainer_["NofBFromTop"]->Fill(genEvt.numberOfBQuarks(true));

  ///////////////////////////////////////////////// 
  // Mass of particles with different status     //
  ///////////////////////////////////////////////// 

  // particles status 4
  //genEvt.setDefaultStatus (4);
  reco::Particle::LorentzVector ttbar;
  if(genEvt.top () && genEvt.topBar ()){
    ttbar = genEvt.top ()->p4 ()+ genEvt.topBar ()->p4 ();
    TH1Fcontainer_["Mtt"]->Fill(ttbar.mass());
    TH1Fcontainer_["PtTtbar"]->Fill(ttbar.pt());
    TH1Fcontainer_["EtaTtbar"]->Fill(ttbar.eta());
    TH1Fcontainer_["MassTop"]->Fill (genEvt.top ()->p4 ().mass ());
    TH1Fcontainer_["MassTop"]->Fill(genEvt.top()->p4().mass());
    TH1Fcontainer_["MassTop"]->Fill(genEvt.topBar()->p4().mass());
  }
  if(genEvt.wMinus ()) TH1Fcontainer_["MassW"]->Fill (genEvt.wMinus ()->p4 ().mass ());
  if(genEvt.wPlus ()) TH1Fcontainer_["MassW"]->Fill (genEvt.wPlus ()->p4 ().mass ());
  if (genEvt.b ()) TH1Fcontainer_["MassB"]->Fill (genEvt.b ()->p4 ().mass ());
  if (genEvt.bBar ()) TH1Fcontainer_["MassB"]->Fill (genEvt.bBar ()->p4 ().mass ());
  if (genEvt.hadronicDecayQuark ()) TH1Fcontainer_["MasslQuarks"]->Fill (genEvt.hadronicDecayQuark ()->p4 ().mass ());
  if (genEvt.hadronicDecayQuarkBar ())TH1Fcontainer_["MasslQuarks"]->Fill (genEvt.hadronicDecayQuarkBar ()->p4 ().mass ());
  

  ///////////////////////////////////////////////// 
  // For semiLep                                 //
  ///////////////////////////////////////////////// 
  
  if(genEvt.isSemiLeptonic()){
  
  ///////////////////////////////////////////////// 
  // Debug                                       //
  ///////////////////////////////////////////////// 
  
  if(genEvt.numberOfBQuarks()!=2) edm::LogWarning ("NoDataFound_Not2B") << genEvt.numberOfBQuarks()<< " b-quarks coming from Top instead of 1 in semileptonic event" ;
  if(genEvt.numberOfLeptons()!=1) edm::LogWarning ("NoDataFound_Not1LepFromW") << genEvt.numberOfLeptons()<< " leptons coming from W instead of 1 in semileptonic event" ;
  if(!genEvt.singleLepton()) edm::LogWarning ("NoDataFound_Lepton") << "No access to singleLepton in semileptonic event" ;
  if(!genEvt.singleNeutrino()) edm::LogWarning ("NoDataFound_Neutrino") << "No access to singleNeutrino in semileptonic event" ;

  if(!genEvt.hadronicDecayTop()) edm::LogWarning ("NoDataFound_HadTop") << "No access to hadronicDecayTop in semileptonic event" ;
  if(!genEvt.hadronicDecayW()) edm::LogWarning ("NoDataFound_HadW") << "No access to hadronicDecayW in semileptonic event" ;
  if(!genEvt.hadronicDecayB()) edm::LogWarning ("NoDataFound_HadB") << "No access to hadronicDecayB in semileptonic event" ;
  if(!genEvt.hadronicDecayQuark()) edm::LogWarning ("NoDataFound_HadQ") << "No access to hadronicDecayQuark in semileptonic event" ;
  if(!genEvt.hadronicDecayQuarkBar()) edm::LogWarning ("NoDataFound_HadQBar") << "No access to hadronicDecayQuarkBar in semileptonic event" ;
  if(!genEvt.leptonicDecayB()) edm::LogWarning ("NoDataFound_LepB") << "No access to leptonicDecayB in semileptonic event" ;

  if(genEvt.topSisters().size()>0 && !genEvt.topSisters()[0]) edm::LogWarning ("NoAcess_ISR") << " No access to ISR in semileptonic event" ;
  if(genEvt.radiatedGluons(genEvt.leptonicDecayTop()->pdgId()).size()>0 && !genEvt.radiatedGluons(genEvt.leptonicDecayTop()->pdgId())[0]) edm::LogWarning ("NoAccess_LepRad") << " No access to leptonicDecayTopRadiation in semileptonic event"<<endl;
  if(genEvt.radiatedGluons(genEvt.hadronicDecayTop()->pdgId()).size()>0 && !genEvt.radiatedGluons(genEvt.hadronicDecayTop()->pdgId())[0]) edm::LogWarning ("NoAccess_HadRad") << " No access to hadronicDecayTopRadiation in semileptonic event"<<endl;

  
  ///////////////////////////////////////////////// 
  // Kinematics                                  //
  ///////////////////////////////////////////////// 

 
  
   if(genEvt.top()){
    TH1Fcontainer_["TopPt"]->Fill(genEvt.top()->pt());
    TH1Fcontainer_["TopEta"]->Fill(genEvt.top()->eta()); 
   }
   if(genEvt.topBar()){
    TH1Fcontainer_["TopPt"]->Fill(genEvt.topBar()->pt());
    TH1Fcontainer_["TopEta"]->Fill(genEvt.topBar()->eta()); 
   }
   if(genEvt.wMinus()){
    TH1Fcontainer_["WPt"]->Fill(genEvt.wMinus()->pt());
    TH1Fcontainer_["WPt"]->Fill(genEvt.wMinus()->pt());
   }
   if(genEvt.wPlus()){
    TH1Fcontainer_["WPt"]->Fill(genEvt.wPlus()->pt());
    TH1Fcontainer_["WEta"]->Fill(genEvt.wPlus()->eta());
   }
   if(genEvt.b()){
    TH1Fcontainer_["BQuarksPt"]->Fill(genEvt.b()->pt());
    TH1Fcontainer_["BQuarksEta"]->Fill(genEvt.b()->eta());
   }
   if(genEvt.bBar()){
    TH1Fcontainer_["BQuarksPt"]->Fill(genEvt.bBar()->pt());
    TH1Fcontainer_["BQuarksEta"]->Fill(genEvt.bBar()->eta());
   }
   if(genEvt.hadronicDecayQuark()){
    TH1Fcontainer_["lQuarksPt"]->Fill(genEvt.hadronicDecayQuark()->pt()); 
    TH1Fcontainer_["lQuarksEta"]->Fill(genEvt.hadronicDecayQuark()->eta()); 
   }
   if(genEvt.hadronicDecayQuarkBar()){
    TH1Fcontainer_["lQuarksPt"]->Fill(genEvt.hadronicDecayQuarkBar()->pt()); 
    TH1Fcontainer_["lQuarksEta"]->Fill(genEvt.hadronicDecayQuarkBar()->eta()); 
   }
   if(genEvt.singleLepton()){
    TH1Fcontainer_["LeptonPt"]->Fill(genEvt.singleLepton()->pt());
    TH1Fcontainer_["LeptonEta"]->Fill(genEvt.singleLepton()->eta()); 
   }
   if(genEvt.singleNeutrino()){
    TH1Fcontainer_["NeutrinoPt"]->Fill(genEvt.singleNeutrino()->pt()); 
    TH1Fcontainer_["NeutrinoEta"]->Fill(genEvt.singleNeutrino()->eta());
   }

   bool allquarks = false;
   if(genEvt.b() && genEvt.bBar() && genEvt.hadronicDecayQuark() && genEvt.hadronicDecayQuarkBar() && genEvt.hadronicDecayB() && genEvt.leptonicDecayB() ) allquarks = true;
   
   //Declaration of vectors used
   std::vector<reco::Particle::LorentzVector> vec;
   vector<float> DeltaR;
   
   if(allquarks){
    reco::Particle::LorentzVector a,b;
    vec.push_back(genEvt.b()->p4());
    vec.push_back(genEvt.bBar()->p4());
    vec.push_back(genEvt.hadronicDecayQuark()->p4());
    vec.push_back(genEvt.hadronicDecayQuarkBar()->p4());
    std::sort(vec.begin(),vec.end(),LowestPt());
    TH1Fcontainer_["LowestPtQuark"]->Fill(vec[0].pt());
    TH1Fcontainer_["2ndPtQuark"]->Fill(vec[1].pt());
    TH1Fcontainer_["3ndPtQuark"]->Fill(vec[2].pt());
    TH1Fcontainer_["HighestPtQuark"]->Fill(vec[3].pt());
    
    std::sort(vec.begin(),vec.end(),HighestEta());
    TH1Fcontainer_["HighestEtaQuark"]->Fill(abs(vec[0].eta())); 
    //check if methods to sort are correct
   
    TH1Fcontainer_["DeltaRlQuarks"]->Fill(ROOT::Math::VectorUtil::DeltaR(genEvt.hadronicDecayQuark()->p4(),genEvt.hadronicDecayQuarkBar()->p4()));
    DeltaR.push_back(ROOT::Math::VectorUtil::DeltaR(genEvt.hadronicDecayQuark()->p4(),genEvt.hadronicDecayB()->p4()));
    DeltaR.push_back(ROOT::Math::VectorUtil::DeltaR(genEvt.hadronicDecayQuarkBar()->p4(),genEvt.hadronicDecayB()->p4()));
    std::sort(DeltaR.begin(),DeltaR.end(),Lowest());
    TH1Fcontainer_["DeltaRBClosestlQuarksSame"]->Fill(DeltaR[0]);
    DeltaR.push_back(ROOT::Math::VectorUtil::DeltaR(genEvt.hadronicDecayQuark()->p4(),genEvt.leptonicDecayB()->p4()));
    DeltaR.push_back(ROOT::Math::VectorUtil::DeltaR(genEvt.hadronicDecayQuarkBar()->p4(),genEvt.leptonicDecayB()->p4()));
    std::sort(DeltaR.begin(),DeltaR.end(),Lowest());
    TH1Fcontainer_["DeltaRBClosestlQuarks"]->Fill(DeltaR[0]);
    DeltaR.push_back(ROOT::Math::VectorUtil::DeltaR(genEvt.hadronicDecayQuark()->p4(),genEvt.hadronicDecayQuarkBar()->p4()));
    DeltaR.push_back(ROOT::Math::VectorUtil::DeltaR(genEvt.leptonicDecayB()->p4(),genEvt.hadronicDecayB()->p4()));
    std::sort(DeltaR.begin(),DeltaR.end(),Lowest());
    TH1Fcontainer_["MinimalDeltaRQuarks"]->Fill(DeltaR[0]);
    DeltaR.clear();
    if(genEvt.singleLepton()){
     DeltaR.push_back(ROOT::Math::VectorUtil::DeltaR(genEvt.singleLepton()->p4(),genEvt.hadronicDecayQuark()->p4()));
     DeltaR.push_back(ROOT::Math::VectorUtil::DeltaR(genEvt.singleLepton()->p4(),genEvt.hadronicDecayQuarkBar()->p4()));
     DeltaR.push_back(ROOT::Math::VectorUtil::DeltaR(genEvt.singleLepton()->p4(),genEvt.hadronicDecayB()->p4()));
     DeltaR.push_back(ROOT::Math::VectorUtil::DeltaR(genEvt.singleLepton()->p4(),genEvt.leptonicDecayB()->p4()));
     std::sort(DeltaR.begin(),DeltaR.end(),Lowest());
     TH1Fcontainer_["MinimalDeltaRQuarksLepton"]->Fill(DeltaR[0]);
    }
   }
  ///////////////////////////////////////////////// 
  //   Radiations                                //
  ///////////////////////////////////////////////// 
 

   DeltaR.clear(); 
   vector<float> PtRad;
   if(genEvt.topSisters().size()==0){
    TH1Fcontainer_["ISRPt"]->Fill(-1);
   }
   for(unsigned int x=0;x<genEvt.topSisters().size();x++){
     TH1Fcontainer_["ISRPt"]->Fill(genEvt.topSisters()[x]->pt());
     TH1Fcontainer_["ISREta"]->Fill(genEvt.topSisters()[x]->eta());
     PtRad.push_back(genEvt.topSisters()[x]->pt());
     if(allquarks){
      DeltaR.push_back(ROOT::Math::VectorUtil::DeltaR(genEvt.topSisters()[x]->p4(),genEvt.hadronicDecayQuark()->p4()));
      DeltaR.push_back(ROOT::Math::VectorUtil::DeltaR(genEvt.topSisters()[x]->p4(),genEvt.hadronicDecayQuarkBar()->p4()));
      DeltaR.push_back(ROOT::Math::VectorUtil::DeltaR(genEvt.topSisters()[x]->p4(),genEvt.hadronicDecayB()->p4()));
      DeltaR.push_back(ROOT::Math::VectorUtil::DeltaR(genEvt.topSisters()[x]->p4(),genEvt.leptonicDecayB()->p4()));
    }
   }

   std::sort(DeltaR.begin(),DeltaR.end(),Lowest());
   std::sort(PtRad.begin(),PtRad.end(),Highest());
   if(DeltaR.size()>0) TH1Fcontainer_["DeltaRISRQuark"]->Fill(DeltaR[0]);
   std::sort(vec.begin(),vec.end(),HighestPt());
   if(PtRad.size()>0){
     bool found = false;
     for(unsigned int x=0;x<vec.size();x++)
       if(PtRad[0]>vec[x].pt()){
         TH1Dcontainer_["RankHighestISR"]->Fill(x+1);
         found = true;
	 break;
       }
     if(!found) TH1Dcontainer_["RankHighestISR"]->Fill(vec.size()+1);
   }
   else TH1Dcontainer_["RankHighestISR"]->Fill(999);

   std::sort(vec.begin(),vec.end(),LowestPt());
   int nof = 0;
   if(vec.size()>0)
   for(unsigned int x=0;x<PtRad.size();x++)
     if(PtRad[x]>vec[0].pt()) nof++;
   TH1Dcontainer_["NofISRWithHighestPtThanQuarks"]->Fill(nof);
   DeltaR.clear();
   PtRad.clear();
   if(genEvt.isSemiLeptonic() && genEvt.leptonicDecayTop())
   for(unsigned int x=0;x<genEvt.radiatedGluons(genEvt.leptonicDecayTop()->pdgId()).size();x++){
     TH1Fcontainer_["TopRadiationPt"]->Fill(genEvt.radiatedGluons(genEvt.leptonicDecayTop()->pdgId())[x]->pt());
     TH1Fcontainer_["TopRadiationEta"]->Fill(genEvt.radiatedGluons(genEvt.leptonicDecayTop()->pdgId())[x]->eta()); 
     PtRad.push_back(genEvt.radiatedGluons(genEvt.leptonicDecayTop()->pdgId())[x]->pt());
     if(allquarks){
      DeltaR.push_back(ROOT::Math::VectorUtil::DeltaR(genEvt.radiatedGluons(genEvt.leptonicDecayTop()->pdgId())[x]->p4(),genEvt.hadronicDecayQuark()->p4()));
      DeltaR.push_back(ROOT::Math::VectorUtil::DeltaR(genEvt.radiatedGluons(genEvt.leptonicDecayTop()->pdgId())[x]->p4(),genEvt.hadronicDecayQuarkBar()->p4()));
      DeltaR.push_back(ROOT::Math::VectorUtil::DeltaR(genEvt.radiatedGluons(genEvt.leptonicDecayTop()->pdgId())[x]->p4(),genEvt.hadronicDecayB()->p4()));
      DeltaR.push_back(ROOT::Math::VectorUtil::DeltaR(genEvt.radiatedGluons(genEvt.leptonicDecayTop()->pdgId())[x]->p4(),genEvt.leptonicDecayB()->p4()));
     }
   }
   if(genEvt.isSemiLeptonic() && genEvt.hadronicDecayTop())
   for(unsigned int x=0;x<genEvt.radiatedGluons(genEvt.hadronicDecayTop()->pdgId()).size();x++){
     TH1Fcontainer_["TopRadiationPt"]->Fill(genEvt.radiatedGluons(genEvt.hadronicDecayTop()->pdgId())[x]->pt());
     TH1Fcontainer_["TopRadiationEta"]->Fill(genEvt.radiatedGluons(genEvt.hadronicDecayTop()->pdgId())[x]->eta()); 
     PtRad.push_back(genEvt.radiatedGluons(genEvt.hadronicDecayTop()->pdgId())[x]->pt());
     if(allquarks){
      DeltaR.push_back(ROOT::Math::VectorUtil::DeltaR(genEvt.radiatedGluons(genEvt.hadronicDecayTop()->pdgId())[x]->p4(),genEvt.hadronicDecayQuark()->p4()));
      DeltaR.push_back(ROOT::Math::VectorUtil::DeltaR(genEvt.radiatedGluons(genEvt.hadronicDecayTop()->pdgId())[x]->p4(),genEvt.hadronicDecayQuarkBar()->p4()));
      DeltaR.push_back(ROOT::Math::VectorUtil::DeltaR(genEvt.radiatedGluons(genEvt.hadronicDecayTop()->pdgId())[x]->p4(),genEvt.hadronicDecayB()->p4()));
      DeltaR.push_back(ROOT::Math::VectorUtil::DeltaR(genEvt.radiatedGluons(genEvt.hadronicDecayTop()->pdgId())[x]->p4(),genEvt.hadronicDecayB()->p4()));
     }
   }

   std::sort(DeltaR.begin(),DeltaR.end(),Lowest());
   std::sort(PtRad.begin(),PtRad.end(),Highest());
   if(DeltaR.size()>0) TH1Fcontainer_["DeltaRTopRadiationQuark"]->Fill(DeltaR[0]);
 
   std::sort(vec.begin(),vec.end(),HighestPt());
   if(PtRad.size()>0){
     bool found = false;
     for(unsigned int x=0;x<vec.size();x++)
       if(PtRad[0]>vec[x].pt()){
         TH1Dcontainer_["RankHighestTopRadiation"]->Fill(x+1);
         found = true;
	 break;
       }
     if(!found) TH1Dcontainer_["RankHighestTopRadiation"]->Fill(vec.size()+1);
   }
   else TH1Dcontainer_["RankHighestTopRadiation"]->Fill(999);
   nof=0;
   std::sort(vec.begin(),vec.end(),LowestPt());
   for(unsigned int x=0;x<PtRad.size();x++)
     if(vec.size()>0 && PtRad[x]>vec[0].pt()) nof++;
   TH1Dcontainer_["NofTopRadiationWithHighestPtThanQuarks"]->Fill(nof);
   
   vec.clear();
   }
  
  ///////////////////////////////////////////////// 
  // Spin Correlations                           //
  ///////////////////////////////////////////////// 
 
  //Helicity basis
// build CM sytems
  if(genEvt.isSemiLeptonic(WDecay::kMuon)){
    reco::Particle::LorentzVector topsZMF(genEvt.top()->p4()+genEvt.topBar()->p4());
  
  // boost particle 4-vectors to tt ZMF
    if(genEvt.leptonicDecayTop() && genEvt.hadronicDecayTop() && genEvt.singleLepton() && genEvt.hadronicDecayB() && genEvt.hadronicDecayQuark()&& genEvt.hadronicDecayQuarkBar() ){
    
      reco::Particle::LorentzVector tLeptonicZMF(ROOT::Math::VectorUtil::boost(genEvt.leptonicDecayTop()->p4(), topsZMF.BoostToCM()));
      reco::Particle::LorentzVector tHadronicZMF(ROOT::Math::VectorUtil::boost(genEvt.hadronicDecayTop()->p4(), topsZMF.BoostToCM()));
      reco::Particle::LorentzVector lLeptonicZMF(ROOT::Math::VectorUtil::boost(genEvt.singleLepton()->p4(), topsZMF.BoostToCM()));
      reco::Particle::LorentzVector bHadronicZMF(ROOT::Math::VectorUtil::boost(genEvt.hadronicDecayB()->p4(), topsZMF.BoostToCM()));
      reco::Particle::LorentzVector q1HadronicZMF(ROOT::Math::VectorUtil::boost(genEvt.hadronicDecayQuark()->p4(), topsZMF.BoostToCM()));
      reco::Particle::LorentzVector q2HadronicZMF(ROOT::Math::VectorUtil::boost(genEvt.hadronicDecayQuarkBar()->p4(), topsZMF.BoostToCM()));
     
  //--------------------------------------------------------------------------------
  // build spin basis unit vectors
      reco::Particle::Vector leptHelZMF(tLeptonicZMF.Vect().Unit());    
      reco::Particle::Vector hadrHelZMF(tHadronicZMF.Vect().Unit()); // = -leptHelZMF
      
      reco::Particle::Vector q1HadZMF(q1HadronicZMF.Vect().Unit());
      reco::Particle::Vector q2HadZMF(q2HadronicZMF.Vect().Unit());
      reco::Particle::Vector qHadZMF(q1HadronicZMF.energy() < q2HadronicZMF.energy() ? // only lower energy quark used
				     q1HadZMF :
				     q2HadZMF);
      
  // boost 4-vectors to t(bar) rest frames
  reco::Particle::LorentzVector lLeptonicTRest(ROOT::Math::VectorUtil::boost(lLeptonicZMF, tLeptonicZMF.BoostToCM()));
  reco::Particle::LorentzVector bHadronicTRest(ROOT::Math::VectorUtil::boost(bHadronicZMF, tHadronicZMF.BoostToCM()));
  reco::Particle::LorentzVector q1HadronicTRest(ROOT::Math::VectorUtil::boost(q1HadronicZMF, tHadronicZMF.BoostToCM()));
  reco::Particle::LorentzVector q2HadronicTRest(ROOT::Math::VectorUtil::boost(q2HadronicZMF, tHadronicZMF.BoostToCM()));
  reco::Particle::LorentzVector qHadronicTRest(q1HadronicTRest.energy() < q2HadronicTRest.energy() ? // only lower energy quark used
                                               q1HadronicTRest :
                                               q2HadronicTRest);

  // extract particle directions in t(bar) rest frames
  reco::Particle::Vector lDirectionTRest(lLeptonicTRest.Vect().Unit());
  reco::Particle::Vector bDirectionTRest(bHadronicTRest.Vect().Unit());
  reco::Particle::Vector q1DirectionTRest(q1HadronicTRest.Vect().Unit());
  reco::Particle::Vector q2DirectionTRest(q2HadronicTRest.Vect().Unit());
  reco::Particle::Vector qDirectionTRest(qHadronicTRest.Vect().Unit());
  
  TH1Dcontainer_["cosThetaTLHel"]->Fill (leptHelZMF.Dot(lDirectionTRest));
  TH1Dcontainer_["cosThetaTBHel"]->Fill (hadrHelZMF.Dot(bDirectionTRest));
  TH1Dcontainer_["cosThetaTQHel"]->Fill (hadrHelZMF.Dot(qDirectionTRest));
  TH2Dcontainer_["hCosTQCosTLHel"]->Fill ( hadrHelZMF.Dot(qDirectionTRest),leptHelZMF.Dot(lDirectionTRest) );
  TH2Dcontainer_["hCosTQCosTBHel"]->Fill ( hadrHelZMF.Dot(bDirectionTRest),leptHelZMF.Dot(lDirectionTRest) );
    }// close if 
  }//close if semilep
  

  }
  

}


// ------------ method called once each job just before starting event loop  ------------
void
TtGenEventChecker::beginJob (const edm::EventSetup &)
{
  edm::Service < TFileService > fs;
  if (!fs)
    throw edm::Exception (edm::errors::Configuration, "TFileService missing from configuration!");

  TFileDirectory subDirKin = fs->mkdir ("Kinematics");
  TFileDirectory subDirNof = fs->mkdir ("Multiplicity");
  TFileDirectory subDirRad = fs->mkdir ("Radiation");
  TFileDirectory subDirComp = fs->mkdir ("Comp");
  TFileDirectory subDirSpinCorr = fs->mkdir ("SpinCorr");

  //define the histograms booked

  //Branching ratio
  TH1Dcontainer_["BRatio"] = fs->make < TH1D > ("BRatio", "BRatio", 10, 0, 10);
  TH1Dcontainer_["BRatio"]->GetXaxis()->SetBinLabel(1,"All");
  TH1Dcontainer_["BRatio"]->GetXaxis()->SetBinLabel(2,"ttbar");
  TH1Dcontainer_["BRatio"]->GetXaxis()->SetBinLabel(3,"Hadr");
  TH1Dcontainer_["BRatio"]->GetXaxis()->SetBinLabel(4,"Lept");
  TH1Dcontainer_["BRatio"]->GetXaxis()->SetBinLabel(5,"Semi");
  TH1Dcontainer_["BRatio"]->GetXaxis()->SetBinLabel(6,"SemiMu");
  TH1Dcontainer_["BRatio"]->GetXaxis()->SetBinLabel(7,"SemiEl");
  TH1Dcontainer_["BRatio"]->GetXaxis()->SetBinLabel(8,"SemiTau");
  TH1Dcontainer_["BRatio"]->GetXaxis()->SetBinLabel(9,"SemiOther");
  
  //Nof 
  TH1Dcontainer_["NofLepton"] = subDirNof.make < TH1D > ("NofLepton", "Nof Leptons", 10, 0, 10);
  TH1Dcontainer_["NofLeptonFromW"] = subDirNof.make < TH1D > ("NofLeptonFromW", "Nof Leptons from W decay", 10, 0, 10);
  TH1Dcontainer_["NofB"] = subDirNof.make < TH1D > ("NofB", "Nof B quarks ", 10, 0, 10);
  TH1Dcontainer_["NofBFromTop"] = subDirNof.make < TH1D > ("NofBFromTop", "Nof B quarks from Top decay ", 10, 0, 10);
  TH1Dcontainer_["NofTopsRadiation"] = subDirNof.make < TH1D > ("NofTopRadiation", "Nof Top Radiation", 10, 0, 10);
  TH1Dcontainer_["NofISR"] = subDirNof.make < TH1D > ("NofISR", "Nof ISR", 10, 0, 10);

  //Kinematics

  TH1Fcontainer_["TopPt"] = subDirKin.make <TH1F> ("TopPt","Pt of Top Status 4",100,0,500);
  TH1Fcontainer_["WPt"] = subDirKin.make <TH1F> ("WPt","Pt of W Status 4",100,0,500);
  TH1Fcontainer_["BQuarksPt"] = subDirKin.make <TH1F> ("BQuarksPt","Pt of b-quarks Status 4",100,0,500);
  TH1Fcontainer_["lQuarksPt"] = subDirKin.make <TH1F> ("lQuarksPt","Pt of l-quarks Status 4",100,0,500);
  TH1Fcontainer_["LeptonPt"] = subDirKin.make <TH1F> ("LeptonPt","Pt of lepton Status 4",100,0,500);
  TH1Fcontainer_["NeutrinoPt"] = subDirKin.make <TH1F> ("NeutrinoPt","Pt of neutrino Status 4",100,0,500);
  TH1Fcontainer_["TopEta"] = subDirKin.make <TH1F> ("TopEta","Eta of Top Status 4",100,0,5);
  TH1Fcontainer_["WEta"] = subDirKin.make <TH1F> ("WEta","Eta of W Status 4",100,0,5);
  TH1Fcontainer_["BQuarksEta"] = subDirKin.make <TH1F> ("BQuarksEta","Eta of b-quarks Status 4",100,0,5);
  TH1Fcontainer_["lQuarksEta"] = subDirKin.make <TH1F> ("lQuarksEta","Eta of l-quarks Status 4",100,0,5);
  TH1Fcontainer_["LeptonEta"] = subDirKin.make <TH1F> ("LeptonEta","Eta of lepton Status 4",100,0,5);
  TH1Fcontainer_["NeutrinoEta"] = subDirKin.make <TH1F> ("NeutrinoEta","Eta of neutrino Status 4",100,0,5);
  TH1Fcontainer_["LowestPtQuark"] = subDirKin.make <TH1F> ("LowestPtQuark","Pt of lowest Pt quark Status 4",500,0,500);
  TH1Fcontainer_["HighestPtQuark"] = subDirKin.make <TH1F> ("HighestPtQuark","Pt of hightest Pt quark Status 4",500,0,500);
  TH1Fcontainer_["2ndPtQuark"] = subDirKin.make <TH1F> ("2ndPtQuark","Pt of 2nd Pt quark Status 4",500,0,500);
  TH1Fcontainer_["3ndPtQuark"] = subDirKin.make <TH1F> ("3ndPtQuark","Pt of 3nd Pt quark Status 4",500,0,500);
  TH1Fcontainer_["HighestEtaQuark"] = subDirKin.make <TH1F> ("HighestEtaQuark","Eta of hightest Eta quark Status 4",500,0,5);
  TH1Fcontainer_["DeltaRlQuarks"] = subDirKin.make <TH1F> ("DeltaRlQuarks","#Delta(R) between l-Quarks",500,0,5);
  TH1Fcontainer_["DeltaRBClosestlQuarksSame"] = subDirKin.make <TH1F> ("DeltaRBClosestlQuarksSame","#Delta(R) between B and closest l-Quarks (same top) ",500,0,5);
  TH1Fcontainer_["DeltaRBClosestlQuarks"] = subDirKin.make <TH1F> ("DeltaRBClosestlQuarks","#Delta(R) between B and closest l-Quarks ",500,0,5);
  TH1Fcontainer_["MinimalDeltaRQuarks"] = subDirKin.make <TH1F> ("MinimalDeltaRQuarks","Minimal #Delta(R) between quarks ",500,0,5);
  TH1Fcontainer_["MinimalDeltaRQuarksLepton"] = subDirKin.make <TH1F> ("MinimalDeltaRQuarksLepton","Minimal #Delta(R) between quarks and lepton ",500,0,5);
  
  //Radiations
  TH1Dcontainer_["RankHighestISR"] = subDirRad.make < TH1D > ("RankHighestISR", "Rank of the hightest ISR", 10, 0, 10);
  TH1Dcontainer_["RankHighestTopRadiation"] = subDirRad.make < TH1D > ("RankHighestTopRadiation", "Rank of the hightest TopRadiation", 10, 0, 10);
  TH1Fcontainer_["ISRPt"] = subDirRad.make <TH1F> ("ISRPt","Pt of ISR",50,0,100);
  TH1Fcontainer_["TopRadiationPt"] = subDirRad.make <TH1F> ("TopRadiationPt","Pt of Top Radiation",50,0,100);
  TH1Fcontainer_["ISREta"] = subDirRad.make <TH1F> ("ISREta","Eta of ISR",500,0,5);
  TH1Fcontainer_["TopRadiationEta"] = subDirRad.make <TH1F> ("TopRadiationEta","Eta of Top Radiation",50,-5,5);
  TH1Fcontainer_["DeltaRISRQuark"] = subDirRad.make <TH1F> ("DeltaRISRQuark","#Delta(R) between ISR and closest quark",500,0,5);
  TH1Fcontainer_["DeltaRTopRadiationQuark"] = subDirRad.make <TH1F> ("DeltaRTopRadiationQuark","#Delta(R) between Top radiation and closest quark",500,0,5);
  TH1Dcontainer_["NofISRWithHighestPtThanQuarks"] = subDirRad.make < TH1D > ("NofISRWithHighestPtThanQuarks", "Nof ISR with Highest Pt than quarks", 10, 0, 10);
  TH1Dcontainer_["NofTopRadiationWithHighestPtThanQuarks"] = subDirRad.make < TH1D > ("NofTopRadiationWithHighestPtThanQuarks", "Nof Top radiation with Highest Pt than quarks", 10, 0, 10);
   

  //Mtt
  
  TH1Fcontainer_["Mtt"] = subDirComp.make <TH1F>("Mtt","Mttbar status 4",150,0,1500);
  TH1Fcontainer_["PtTtbar"] = subDirComp.make <TH1F>("PtTtbar","PtTtbar status 2",100,0,500);
  TH1Fcontainer_["EtaTtbar"] = subDirComp.make <TH1F>("EtaTtbar","EtaTtbar status 2",100,-5,5);


 
  //masses
  TH1Fcontainer_["MassTop"] = subDirComp.make < TH1F > ("MassTop", "Mass of Top Status 4", 100, 0, 300);
  TH1Fcontainer_["MassW"] = subDirComp.make < TH1F > ("MassW", "Mass of W Status 4", 60, 0, 120);
  TH1Fcontainer_["MassB"] = subDirComp.make < TH1F > ("MassB", "Mass of B Status 4", 240, 0, 120);
  TH1Fcontainer_["MasslQuarks"] = subDirComp.make < TH1F > ("MasslQuarks", "Mass of light quarks Status 4", 120, 0, 30);

  //Spin Correlations
  TH1Dcontainer_["cosThetaTLHel"] = subDirSpinCorr.make < TH1D > ("cosThetaTLHel", "Pseudo-angle between t-quark and lepton", 10, -1, 1);
  TH1Dcontainer_["cosThetaTBHel"] = subDirSpinCorr.make < TH1D > ("cosThetaTBHel", "Pseudo-angle between t-quark and b-quark", 10, -1, 1);
  TH1Dcontainer_["cosThetaTQHel"] = subDirSpinCorr.make < TH1D > ("cosThetaTQHel", "Pseudo-angle between t-quark and low-energy quark", 10, -1, 1);
  TH2Dcontainer_["hCosTQCosTLHel"] = subDirSpinCorr.make < TH2D > ("hCosTQCosTLHel", "Double differential distribution: t-quark and lepton/low-energy quark", 10, -1, 1, 10, -1, 1);
  TH2Dcontainer_["hCosTQCosTBHel"] = subDirSpinCorr.make < TH2D > ("hCosTQCosTBHel", "Double differential distribution: t-quark and b/low-energy quark", 10, -1, 1, 10, -1, 1);
  //fit the double differential distributions
  fit2LB_ = TF2("fit2LB_", "[0]*(1.-[1]*[2]*[3]*x*y)");
  fit2LB_.SetParameter(2,1);
  fit2LB_.FixParameter(2,1);
  fit2LB_.SetParameter(3,-0.41);
  fit2LB_.FixParameter(3,-0.41); 
  fit2LQ_ = TF2("fit2LQ_", "[0]*(1.-[1]*[2]*[3]*x*y)");
  fit2LQ_.SetParameter(2,1);
  fit2LQ_.FixParameter(2,1);
  fit2LQ_.SetParameter(3,0.51);
  fit2LQ_.FixParameter(3,0.51);

}

// ------------ method called once each job just after ending the event loop  ------------
void
TtGenEventChecker::endJob ()
{

  TH2Dcontainer_["hCosTQCosTLHel"]->Fit(&fit2LQ_,"0");
  //  TH2Dcontainer_["hCosTQCosTBHel"]->Fit(&fit2LB_,"0");

  //use LogError to summarise the error that happen in the execution (by example from warning) (ex: Nof where we cannot access such variable)
  //edm::LogError ("SummaryError") << "My error message \n";	// or  edm::LogProblem  (not formated)
  //use LogVerbatim to summarise information (ex: pourcentage of events matched ...)
  //edm::LogVerbatim ("MainResults") << "My LogVerbatim message \n";	// or  edm::LogVerbatim (not formated)

  edm::LogVerbatim ("MainResults") << " -------------------------------------------";
  edm::LogVerbatim ("MainResults") << " -------------------------------------------";
  edm::LogVerbatim ("MainResults") << " -- Report from TtGenEventSanityChecker   -- ";
  edm::LogVerbatim ("MainResults") << " -------------------------------------------";
  edm::LogVerbatim ("MainResults") << " -------------------------------------------";
  
  float a = 0;
  edm::LogVerbatim ("MainResults") << " -------------------------------";
  edm::LogVerbatim ("MainResults") << "  Info from semileptonic events";
  edm::LogVerbatim ("MainResults") << " -------------------------------";
  //take care
  //here born of integral are hard coded and depend on the numberOfbins !!
  a = TH1Fcontainer_["LowestPtQuark"]->Integral(0,15)/ TH1Fcontainer_["LowestPtQuark"]->GetEntries();
  edm::LogVerbatim ("MainResults") << " Quark with Pt < 15 GeV: "<<a*100<<" %";
  a = TH1Fcontainer_["LowestPtQuark"]->Integral(0,20)/ TH1Fcontainer_["LowestPtQuark"]->GetEntries();
  edm::LogVerbatim ("MainResults") << " Quark with Pt < 20 GeV: "<<a*100<<" %";
  a = TH1Fcontainer_["LowestPtQuark"]->Integral(0,25)/ TH1Fcontainer_["LowestPtQuark"]->GetEntries();
  edm::LogVerbatim ("MainResults") << " Quark with Pt < 25 GeV: "<<a*100<<" %";
  a = TH1Fcontainer_["LowestPtQuark"]->Integral(0,30)/ TH1Fcontainer_["LowestPtQuark"]->GetEntries();
  edm::LogVerbatim ("MainResults") << " Quark with Pt < 30 GeV: "<<a*100<<" %";
  a = TH1Fcontainer_["HighestEtaQuark"]->Integral(0,240)/TH1Fcontainer_["HighestEtaQuark"]->GetEntries();
  edm::LogVerbatim ("MainResults") << " Quark with Eta > 2.4: "<<(1-a)*100<<" %";
  a = TH1Fcontainer_["HighestEtaQuark"]->Integral(0,300)/TH1Fcontainer_["HighestEtaQuark"]->GetEntries();
  edm::LogVerbatim ("MainResults") << " Quark with Eta > 3: "<<(1-a)*100<<" %";
  a = TH1Fcontainer_["HighestEtaQuark"]->Integral(0,400)/TH1Fcontainer_["HighestEtaQuark"]->GetEntries();
  edm::LogVerbatim ("MainResults") << " Quark with Eta > 4: "<<(1-a)*100<<" %";
  a = TH1Fcontainer_["DeltaRlQuarks"]->Integral(0,50)/TH1Fcontainer_["DeltaRlQuarks"]->GetEntries();
  edm::LogVerbatim ("MainResults") << " l-Quarks with DeltaR < 0.5: "<<a*100<<" %";
  a = TH1Fcontainer_["DeltaRBClosestlQuarksSame"]->Integral(0,50)/TH1Fcontainer_["DeltaRBClosestlQuarksSame"]->GetEntries();
  edm::LogVerbatim ("MainResults") << " l-Quarks & B (same top) with DeltaR < 0.5: "<<a*100<<" %";
  a = TH1Fcontainer_["DeltaRBClosestlQuarks"]->Integral(0,50)/TH1Fcontainer_["DeltaRBClosestlQuarks"]->GetEntries();
  edm::LogVerbatim ("MainResults") << " l-Quarks & B  with DeltaR < 0.5: "<<a*100<<" %";
  a = TH1Fcontainer_["MinimalDeltaRQuarks"]->Integral(0,50)/TH1Fcontainer_["MinimalDeltaRQuarks"]->GetEntries();
  edm::LogVerbatim ("MainResults") << " minimal DeltaR  between 4 quarks < 0.5: "<<a*100<<" %";
  a = TH1Fcontainer_["MinimalDeltaRQuarks"]->Integral(0,100)/TH1Fcontainer_["MinimalDeltaRQuarks"]->GetEntries();
  edm::LogVerbatim ("MainResults") << " minimal DeltaR  between 4 quarks < 1.0: "<<a*100<<" %";
  a = TH1Fcontainer_["MinimalDeltaRQuarksLepton"]->Integral(0,30)/TH1Fcontainer_["MinimalDeltaRQuarksLepton"]->GetEntries();
  edm::LogVerbatim ("MainResults") << " DeltaR(quarks, muon) < 0.3: "<<a*100<<" %";
  //

  //Radiations
  //take care
  //here born of integral are hard coded and depend on the numberOfbins !!
  a = TH1Dcontainer_["RankHighestISR"]->Integral(0,5)/TH1Dcontainer_["RankHighestISR"]->GetEntries();
  edm::LogVerbatim ("MainResults") << " ISR with Pt > quarks: "<<a*100<<" %";
  a = TH1Dcontainer_["RankHighestTopRadiation"]->Integral(0,5)/TH1Dcontainer_["RankHighestTopRadiation"]->GetEntries();
  edm::LogVerbatim ("MainResults") << " Top Radiation with Pt > quarks: "<<a*100<<" %";
  a = TH1Fcontainer_["ISRPt"]->Integral(0,30)/TH1Fcontainer_["ISRPt"]->GetEntries();
  edm::LogVerbatim ("MainResults") << " ISR  with Pt > 30 GeV: "<<(1-a)*100<<" %";
  a = TH1Fcontainer_["TopRadiationPt"]->Integral(0,30)/TH1Fcontainer_["TopRadiationPt"]->GetEntries();
  edm::LogVerbatim ("MainResults") << " TopRadiation  with Pt > 30 GeV: "<<(1-a)*100<<" %";
  //
 
  a = TH1Dcontainer_["NofISRWithHighestPtThanQuarks"]->Integral(0,1)/TH1Dcontainer_["NofISRWithHighestPtThanQuarks"]->GetEntries();
  edm::LogVerbatim ("MainResults") << " At least one ISR  with Pt > Quarks: "<<(1-a)*100<<" %";
  a = TH1Dcontainer_["NofISRWithHighestPtThanQuarks"]->Integral(0,2)/TH1Dcontainer_["NofISRWithHighestPtThanQuarks"]->GetEntries();
  edm::LogVerbatim ("MainResults") << " At least 2 ISR  with Pt > Quarks: "<<(1-a)*100<<" %";
  a = TH1Dcontainer_["NofISRWithHighestPtThanQuarks"]->Integral(0,3)/TH1Dcontainer_["NofISRWithHighestPtThanQuarks"]->GetEntries();
  edm::LogVerbatim ("MainResults") << " At least 3 ISR  with Pt > Quarks: "<<(1-a)*100<<" %";
  a = TH1Dcontainer_["NofISRWithHighestPtThanQuarks"]->Integral(0,4)/TH1Dcontainer_["NofISRWithHighestPtThanQuarks"]->GetEntries();
  edm::LogVerbatim ("MainResults") << " At least 4 ISR  with Pt > Quarks: "<<(1-a)*100<<" %";
  a = TH1Dcontainer_["NofTopRadiationWithHighestPtThanQuarks"]->Integral(0,1)/TH1Dcontainer_["NofTopRadiationWithHighestPtThanQuarks"]->GetEntries();
  edm::LogVerbatim ("MainResults") << " At least one TopRadiation  with Pt > Quarks: "<<(1-a)*100<<" %";
  a = TH1Dcontainer_["NofTopRadiationWithHighestPtThanQuarks"]->Integral(0,2)/TH1Dcontainer_["NofTopRadiationWithHighestPtThanQuarks"]->GetEntries();
  edm::LogVerbatim ("MainResults") << " At least 2 TopRadiation  with Pt > Quarks: "<<(1-a)*100<<" %";
  a = TH1Dcontainer_["NofTopRadiationWithHighestPtThanQuarks"]->Integral(0,3)/TH1Dcontainer_["NofTopRadiationWithHighestPtThanQuarks"]->GetEntries();
  edm::LogVerbatim ("MainResults") << " At least 3 TopRadiation  with Pt > Quarks: "<<(1-a)*100<<" %";
  a = TH1Dcontainer_["NofTopRadiationWithHighestPtThanQuarks"]->Integral(0,4)/TH1Dcontainer_["NofTopRadiationWithHighestPtThanQuarks"]->GetEntries();
  edm::LogVerbatim ("MainResults") << " At least 4 TopRadiation  with Pt > Quarks: "<<(1-a)*100<<" %";
   
  // Bquarks
  edm::LogVerbatim ("MainResults") << "-------------------- ";
  edm::LogVerbatim ("MainResults") << " Other ";
  edm::LogVerbatim ("MainResults") << "-------------------- ";
  a = TH1Dcontainer_["NofBFromTop"]->Integral(0,2)/TH1Dcontainer_["NofBFromTop"]->GetEntries();
  edm::LogVerbatim ("MainResults") << " Events with 1 b-quarks instead of 2: "<<a*100<<"%";
  if(TH1Dcontainer_["NofTopsRadiation"]->Integral(0,1) == TH1Dcontainer_["NofTopsRadiation"]->GetEntries()) 
    edm::LogVerbatim ("MainResults") << " There is no top radiations in this sample";
  if(TH1Dcontainer_["NofISR"]->Integral(0,1) == TH1Dcontainer_["NofISR"]->GetEntries()) 
    edm::LogVerbatim ("MainResults") << " There is no ISR in this sample (normally present in MadGraph sample)";

  
  //Branching Ratio 
  edm::LogVerbatim ("MainResults") << "-------------------- ";
  edm::LogVerbatim ("MainResults") << " Branching ratio ";
  edm::LogVerbatim ("MainResults") << "-------------------- ";
  a = (float) (TH1Dcontainer_["BRatio"]->GetBinContent(2)/TH1Dcontainer_["BRatio"]->GetBinContent(1));
  edm::LogVerbatim ("MainResults") << "ttbar evt: "<<a*100<<" %";
  a = (float) (TH1Dcontainer_["BRatio"]->GetBinContent(3)/TH1Dcontainer_["BRatio"]->GetBinContent(2));
  edm::LogVerbatim ("MainResults") << "FullHadronic: "<<a*100<<" %";
  a = (float) (TH1Dcontainer_["BRatio"]->GetBinContent(4)/TH1Dcontainer_["BRatio"]->GetBinContent(2));
  edm::LogVerbatim ("MainResults") << "FullLeptonic: "<<a*100<<" %";
  a = (float) (TH1Dcontainer_["BRatio"]->GetBinContent(5)/TH1Dcontainer_["BRatio"]->GetBinContent(2));
  edm::LogVerbatim ("MainResults") << "SemiLeptonic: "<<a*100<<" %";
  a = (float) (TH1Dcontainer_["BRatio"]->GetBinContent(6)/TH1Dcontainer_["BRatio"]->GetBinContent(2));
  edm::LogVerbatim ("MainResults") << "SemiLeptonic Muon: "<<a*100<<" %";
  a = (float) (TH1Dcontainer_["BRatio"]->GetBinContent(7)/TH1Dcontainer_["BRatio"]->GetBinContent(2));
  edm::LogVerbatim ("MainResults") << "SemiLeptonic Electron: "<<a*100<<" %";
  a = (float) (TH1Dcontainer_["BRatio"]->GetBinContent(8)/TH1Dcontainer_["BRatio"]->GetBinContent(2));
  edm::LogVerbatim ("MainResults") << "SemiLeptonic Tau: "<<a*100<<" %";
  a = (float) (TH1Dcontainer_["BRatio"]->GetBinContent(9)/TH1Dcontainer_["BRatio"]->GetBinContent(2));
  edm::LogVerbatim ("MainResults") << "SemiLeptonic Other: "<<a*100<<" %";
  
  edm::LogVerbatim ("MainResults") << "-------------------- ";
  edm::LogVerbatim ("MainResults") << " Spin Correlations ";
  edm::LogVerbatim ("MainResults") << "-------------------- ";
  edm::LogVerbatim ("MainResults") << " CosTQCosTL A =  " <<  fit2LQ_.GetParameter(1);;
  //  edm::LogVerbatim ("MainResults") << " CosTLCosTB A =  " <<  fit2LB_.GetParameter(1);;


}

//define this as a plug-in
DEFINE_FWK_MODULE (TtGenEventChecker);
