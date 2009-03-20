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
#include "TopQuarkAnalysis/TopTools/interface/JetPartonMatching.h"
#include "DataFormats/Candidate/interface/Candidate.h"

#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/MET.h"

#include "DataFormats/Math/interface/deltaR.h"

#include "TDirectory.h"
#include "TH1D.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TF1.h"

#include <Math/VectorUtil.h>

//
// class declaration
//

class TruthReco : public edm::EDAnalyzer {
   public:
      explicit TruthReco(const edm::ParameterSet&);
      ~TruthReco();

   private:
      virtual void beginJob(const edm::EventSetup&) ;
      virtual void analyze(const edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;

      // ----------member data ---------------------------
      edm::InputTag genEventCollectionName_,patMuons_,patJets_,patMET_;
      bool useMatchingFromPAT_;
      int matchingAlgo_;
      bool useMaxDist_,useDeltaR_;
      double maxDist_;
			bool doMuonIsolation_;
			bool verbose_;
			double muonRelIso_,muonECALVetoConeEt_,muonHCALVetoConeEt_,muonPt_,jetPt_,muonEta_,jetEta_,muonMinDR_;
			bool doNonOverlappingJets_;
			double nonOverlappingJetsMinDR_;
			bool calculateJECinbins_;
			std::vector<double> etabinVals_,pTbinVals_;
			int NJets_;
			double minDRunmatchedParton_;
			
			int NbSemiMuRecoEvents,NbSemiMuGenEvents,recoEventsPassed,genObjectsNotNULL,genObjectsNULL,EventNotMatched,ptnrbins,etanrbins;
			int UnmatchedQuarks,MatchedQuarks,MatchedISR,UnmatchedISR,UnmatchedTopRadiation,MatchedTopRadiation,NoTopRadiation,NoISR,NoRadiation,ISRTopRadiation,TopRadiationNoISR;
			int Highest4NoRad,Highest4AtLeast1Rad,Highest4AtLeast2Rad,Highest4AtLeast3Rad,Highest4AtLeast4Rad;
			int UnmatchedQuarkCategoryEta,UnmatchedQuarkCategoryPt,UnmatchedQuarkCategoryDR,UnmatchedQuarkNotCategorized;
			int UnmatchedISRCategoryEta,UnmatchedISRCategoryPt,UnmatchedISRCategoryDR,UnmatchedISRNotCategorized;
			int UnmatchedRadCategoryEta,UnmatchedRadCategoryPt,UnmatchedRadCategoryDR,UnmatchedRadNotCategorized;
			int NrNotSemiMu,NrAnalyzed,NrNot4Jets,NrNot1Mu,NrNot1MET,NrNot4SelJets,NrNot1SelMu,NrNot1SelIsoMu;
			unsigned int NumberQuarks,NumberISR,NumberRadiationLep,NumberRadiationHad,NumberRadiation;	
			int MatchedQuarksRadiation, MatchedQuarksNoRadiation, MatchedQuarksNoISR, MatchedQuarksNoTopRadiation;
			int UnmatchedQuark1NoRadCategoryEta,UnmatchedQuark1NoRadCategoryPt,UnmatchedQuark1NoRadCategoryDR,UnmatchedQuark1NoRadNotCategorized;
			int UnmatchedQuark1CategoryEta,UnmatchedQuark1CategoryPt,UnmatchedQuark1CategoryDR,UnmatchedQuark1NotCategorized;
			int Only1UnmatchedQuarkCategoryEta,Only1UnmatchedQuarkCategoryPt,Only1UnmatchedQuarkCategoryDR,Only1UnmatchedQuarkNotCategorized;			
			int Only1UnmatchedQuarks;
			int	QuarksMatchedUnmatchedISRCategoryEta,QuarksMatchedUnmatchedISRCategoryPt,QuarksMatchedUnmatchedISRCategoryDR,QuarksMatchedUnmatchedISRNotCategorized;
			int QuarksMatchedMatchedISR,QuarksMatchedUnmatchedISR,QuarksMatchedISRTopRadiation,QuarksMatchedMatchedTopRadiation,QuarksMatchedUnmatchedTopRadiation;
			int QuarksMatchedUnmatchedRadCategoryEta,QuarksMatchedUnmatchedRadCategoryPt,QuarksMatchedUnmatchedRadCategoryDR,QuarksMatchedUnmatchedRadNotCategorized;
			int QuarksMatchedNoTopRadiation,QuarksMatchedUnmatchedISRtotal,QuarksMatchedUnmatchedTopRadiationtotal;
			int QuarksMatchedNoISR,QuarksMatchedTopRadiationNoISR;		
			int QuarksMatchedNoRadiation;
			int QuarksMatchedHighest4NoRad,QuarksMatchedHighest4AtLeast1Rad,QuarksMatchedHighest4AtLeast2Rad,QuarksMatchedHighest4AtLeast3Rad,QuarksMatchedHighest4AtLeast4Rad;
			
			//Histograms are booked in the beginJob() method
  		TH1F 	*hExpLightJECPtEtaBin[20][20],*hExpBJECPtEtaBin[20][20];
  		TF1 	*fExpLightJECPtEtaBin[20][20],*fExpBJECPtEtaBin[20][20];
  		TH1F 	*hExpLightJECIncl,*hExpBJECIncl;
  		TF1 	*fExpLightJECIncl,*fExpBJECIncl;
			TH1F 	*LowestPtMatchedJet_Rank,	*LowestPtMatchedJet_Pt,*HighestPtRadMatchedJet_Rank,*HighestPtRadMatchedJet_Pt;
			TH1F 	*mWHad_Gen,*mWLep_Gen,*mtHad_Gen,*mtLep_Gen,*mWHad_Rec,*mWLep_Rec,*mtHad_Rec,*mtLep_Rec; 
			TH1F 	*UnmatchedQuarkMultiplicity,*JetRadiationMultiplicity;	
			TH2F 	*minJetMuVSminQuarkMu,*minJetsVSminQuarks;	
			TH1F 	*Dphi_HadTops,*Dphi_LepTops,*Deta_HadTops,*Deta_LepTops;
			TH1F	*DR_HadTops,*DR_LepTops,*Dphi_Tops_Rec,*Dphi_Tops_Gen,*DR_Tops_Rec,*DR_Tops_Gen;
			TH1F	*QuarksMatchedJetRadiationMultiplicity,*QuarksMatchedHighestPtRadMatchedJet_Rank,*QuarksMatchedHighestPtRadMatchedJet_Pt;
			TH2F 	*LowestPtQuarkMatchedJetPt_vs_HighestPtRadMatchedJetPt;


};


//
// constructors and destructor
//
TruthReco::TruthReco(const edm::ParameterSet& iConfig)
{
	genEventCollectionName_ = iConfig.getParameter<edm::InputTag>( "genEventCollectionName" );
	patMuons_								= iConfig.getParameter<edm::InputTag>( "patMuons" );
	patJets_								= iConfig.getParameter<edm::InputTag>( "patJets" );
	patMET_									= iConfig.getParameter<edm::InputTag>( "patMET" );
	useMatchingFromPAT_			= iConfig.getParameter<bool>( "useMatchingFromPAT" );
  matchingAlgo_						= iConfig.getParameter<int>("matchingAlgo");
  useMaxDist_							=	iConfig.getParameter<bool>("useMaxDist");
  useDeltaR_							=	iConfig.getParameter<bool>("useDeltaR");
  maxDist_								=	iConfig.getParameter<double>("maxDist");
	doMuonIsolation_				=	iConfig.getParameter<bool>("doMuonIsolation");
	muonRelIso_							= iConfig.getParameter<double>( "muonRelIso" );
	muonECALVetoConeEt_			= iConfig.getParameter<double>( "muonECALVetoConeEt" );
	muonHCALVetoConeEt_			= iConfig.getParameter<double>( "muonHCALVetoConeEt" );
	muonPt_									= iConfig.getParameter<double>( "muonPt" );
	jetPt_									= iConfig.getParameter<double>( "jetPt" );
	muonEta_								= iConfig.getParameter<double>( "muonEta" );
	jetEta_									= iConfig.getParameter<double>( "jetEta" );
	muonMinDR_							= iConfig.getParameter<double>( "muonMinDR" );
	doNonOverlappingJets_		= iConfig.getParameter<bool>( "doNonOverlappingJets" );
	nonOverlappingJetsMinDR_= iConfig.getParameter<double>( "nonOverlappingJetsMinDR" );
	calculateJECinbins_			= iConfig.getParameter<bool>( "calculateJECinbins" );
	etabinVals_							= iConfig.getParameter< std::vector<double> > 	("etabinValues");
  pTbinVals_							= iConfig.getParameter< std::vector<double> > 	("pTbinValues");
	NJets_									= iConfig.getParameter<int>( "NJets" );
	minDRunmatchedParton_		=	iConfig.getParameter<double>("minDRunmatchedParton");
	verbose_								= iConfig.getParameter<bool>("verbose");
}


TruthReco::~TruthReco()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called to for each event  ------------
void
TruthReco::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
 //module should be run on semi-leptonic ttbar events
 
  using namespace edm;
  using namespace std;
  using namespace pat;

  //Here you handle the collection you want to access
  // PAT Muons
  edm::Handle<std::vector<pat::Muon> > muonHandle;
  iEvent.getByLabel(patMuons_,muonHandle);
  const std::vector<pat::Muon> & muons = *muonHandle;
  // PAT jets
  edm::Handle<std::vector<pat::Jet> > jetHandle;
  iEvent.getByLabel(patJets_,jetHandle);
  const std::vector<pat::Jet> & jets = *jetHandle;
  // PAT MET
  edm::Handle<std::vector<pat::MET> > metHandle;
  iEvent.getByLabel(patMET_,metHandle);
  const std::vector<pat::MET> & mets = *metHandle;
  // Top GenEvt
  edm::Handle<TtGenEvent> genEvtHandle;
  iEvent.getByLabel ("genEvt",genEvtHandle);
  const TtGenEvent & genEvt = *genEvtHandle;
 
	NrAnalyzed++;

	if(verbose_) cout << "analyzing event number " << NrAnalyzed << endl;
	
	if(genEvt.isSemiLeptonic(genEvt.kMuon)) NbSemiMuGenEvents++;
	if(jets.size()>3 && muons.size()>0 && mets.size()>0) NbSemiMuRecoEvents++;

	if(verbose_) cout << "counters done " << endl;
	
	bool allok = true;
	bool muonselected = true;
	bool jetsselected = true;
	
	if(!genEvt.isSemiLeptonic(genEvt.kMuon)) {
		edm::LogWarning  ("NoDataFound") << "Not a semi-muonic event..."; 
		NrNotSemiMu++;
		allok = false;		
	}
	
	if(jets.size()<4){
		NrNot4Jets++;
		allok = false;		
  }
	if(muons.size()<1){
		NrNot1Mu++;
		allok = false;		
  }
	if(mets.size()<1){
		NrNot1MET++;
		allok = false;		
	}

	if(verbose_) cout << "some checks done " << endl;

  std::vector<pat::Jet> myjetsetagood;
  std::vector<pat::Jet> myselectedjets;
  for(unsigned int i =0;i<jets.size();i++) {if(fabs(jets[i].eta())<jetEta_) myjetsetagood.push_back(jets[i]);}
	for(unsigned int i =0;i<myjetsetagood.size();i++) {if(jets[i].pt()>jetPt_ ) myselectedjets.push_back(myjetsetagood[i]);}		

	if(verbose_) cout << "new jet collection filled " << endl;
	
	//check if at least 4 jets with eta ok and pT exceeding jetPt threshold!
	if(myselectedjets.size()<(unsigned)NJets_ ){ 
		NrNot4SelJets++;
		jetsselected = false;
		allok = false;		
		edm::LogWarning  ("NoDataFound") << "Less then 4 jets with p_{T}>"<<jetPt_<<"GeV and |#eta| <"<<jetEta_<<"\n"; 
	}

	if(verbose_) cout << "check if 4 jets " << endl;

	//check if muon passes eta/pT cut!
	if( jetsselected && (muons.size()<1 || muons[0].pt()<muonPt_ || fabs(muons[0].eta())>muonEta_)){
		NrNot1SelMu++;
		allok = false;			
		muonselected = false;
		edm::LogWarning  ("NoDataFound") << "The leading muon has no p_{T}>"<<muonPt_<<"GeV or |#eta| >"<<muonEta_<<"\n"; 
	}

	if(verbose_) cout << "check muon kinematics " << endl;

	double RelIsoMuon = 0.;
 	if(muons.size()>0){
 		RelIsoMuon = muons[0].pt()/(muons[0].pt()+muons[0].isolationR03().sumPt+muons[0].isolationR03().emEt+muons[0].isolationR03().hadEt);
	}
	if(verbose_) cout << "calculated reliso " << endl;

	//veto cone size is 0.07 in ECAL and 0.1 in HCAL, proposed ET cuts are resp 4 and 6
	const reco::IsoDeposit * ecalIsoDep;
	const reco::IsoDeposit * hcalIsoDep;
 	if(muons.size()>0){
		ecalIsoDep = (muons)[0].ecalIsoDeposit();
		hcalIsoDep = (muons)[0].hcalIsoDeposit();
	}

	if(verbose_) cout << "calculated veto isolation " << endl;
	
	//if muonisolation is required check if muon passes reliso and veto cone energy requirements!
	if(doMuonIsolation_ && jetsselected && muons.size()>0){
	  if(muonselected && (RelIsoMuon<muonRelIso_ || ecalIsoDep->candEnergy()>muonECALVetoConeEt_ || hcalIsoDep->candEnergy()>muonHCALVetoConeEt_)){
			NrNot1SelIsoMu++;
			allok = false;			
			edm::LogWarning  ("NoDataFound") << "The leading muon does not pass the isolation requirements \n"; 
		}
	}

	if(verbose_) cout << "check isolation criteria " << endl;

	
	if(!useMatchingFromPAT_ && allok){ 
 		if(verbose_) cout << "event passed selection cuts" << endl;
	 	// Matching index : Hadronic Q  = 0, Hadronic Q' = 1, Hadronic b  = 2, Leptonic b  = 3;
  	std::vector<const reco::Candidate *> quarks;
		const reco::Candidate * muon;
		const reco::Candidate * neutrino;	
		
		recoEventsPassed++;
		
		if(genEvt.hadronicDecayQuark() != 0 && genEvt.hadronicDecayQuarkBar() != 0  && genEvt.hadronicDecayB() != 0 
		&& genEvt.leptonicDecayB() != 0 && genEvt.singleLepton() != 0 && genEvt.singleNeutrino() != 0){ 		
		
	 		if(verbose_) cout << "event is semi-leptonic and all generated objects are not null" << endl;		
			genObjectsNotNULL++; 

  		quarks.push_back(genEvt.hadronicDecayQuark());
  		quarks.push_back(genEvt.hadronicDecayQuarkBar());
  		quarks.push_back(genEvt.hadronicDecayB());
  		quarks.push_back(genEvt.leptonicDecayB());

	 		muon = genEvt.singleLepton();
  		neutrino = genEvt.singleNeutrino();

			if((deltaR(muon->p4(),muons[0].p4()) > muonMinDR_)){
				EventNotMatched++;
	 			if(verbose_) cout << "the reco muon is not matched with the generated muon within "<< muonMinDR_ << endl;		
				edm::LogWarning  ("NoDataFound") << "The muon is not matched!!! \n"; 
	 		}
   	
			std::vector< double > matchjetpt;
			std::vector<const reco::Candidate *> vectorISR, TopRadiation;
			std::vector<pair<int,pat::Jet> > patJetsMatchedWithRadiation;
			std::vector< int > matchquarkjetidx,matchISRjetidx,matchRadiationjetidx,originalrankjetsISR,originalrankjetsradiation;
			std::vector< int > originalrankjetsISRinput,originalrankjetsradiationinput,rankjetsmatchedtorad;
			std::vector<pat::Jet> myselectedjetsISR,myselectedjetsradiation,myselectedjetsradiationNoISR,FourHighestPtJets;
			std::vector< double > ISRpT,TopRadiationpT,TopRadiationpTNoISR,ptjetsmatchedtorad;
			std::vector<const reco::Candidate *> NfirstpTISR,NfirstpTradiation,NfirstpTradiationNoISR;

			///check number additional ISR partons
			unsigned int NumberISR = genEvt.ISR().size();		
			///check number additional radiation partons
			unsigned int NumberRadiationLep = genEvt.leptonicDecayTopRadiation().size();		
			unsigned int NumberRadiationHad = genEvt.hadronicDecayTopRadiation().size();		
			unsigned int NumberRadiation = NumberRadiationLep+NumberRadiationHad;		
	 		if(verbose_){
				cout << "#quarks "<<quarks.size() << endl;
				cout << "#ISR "<<NumberISR << endl;
				cout << "#radiation "<<NumberRadiation << endl;
			}		
		
			if(NumberISR!=0) {
				for(unsigned int i=0; i<NumberISR; i++){
					vectorISR.push_back((genEvt.ISR())[i]);
				}
			}
			if(NumberRadiationLep!=0 ) {
				for(unsigned int i=0;i<NumberRadiationLep; i++){
					TopRadiation.push_back((genEvt.leptonicDecayTopRadiation())[i]);
				}
			}
			if(NumberRadiationHad!=0 ) {
				for(unsigned int i=0; i<NumberRadiationHad; i++){
					TopRadiation.push_back((genEvt.hadronicDecayTopRadiation())[i]);
				}
			}
		
			double DRquarki[4];
			double DRISRi[20];
			double DRRadi[20];
			double smallestDRquarki;
			double smallestDRISRi;
			double smallestDRRadi;

 			if(verbose_) cout << "will match quarks " << endl;		
			//match the quarks
			JetPartonMatching *GenMatchQuarks = new JetPartonMatching(quarks, myselectedjets, matchingAlgo_, useMaxDist_, useDeltaR_, maxDist_);
			
			//check if all quarks are matched ---> if so, we are not interested in radiation
			if(GenMatchQuarks->getNumberOfUnmatchedPartons(0)== 0){
				MatchedQuarks++;				
 				
				if(NumberISR != 0 && NumberRadiation != 0) MatchedQuarksRadiation++;
				if(NumberISR == 0 && NumberRadiation != 0) MatchedQuarksNoISR++;
				if(NumberISR != 0 && NumberRadiation == 0) MatchedQuarksNoTopRadiation++;
				if(NumberISR == 0 && NumberRadiation == 0) MatchedQuarksNoRadiation++; 				
				if(verbose_) cout << "all quarks are matched"<< endl;		

				for(unsigned int i=0; i<quarks.size(); i++){
					///////////////////// ---start --- calculating expected JEC /////////////////////
					if(verbose_) cout << "start calculating expected jet energy corrections"<< endl;		
					double quarkEnergy   = quarks[i]->energy(); 
    			double jetEnergy     = myselectedjets[GenMatchQuarks->getMatchForParton(i,0)].energy(); 
      		if(i<2)hExpLightJECIncl -> Fill((quarkEnergy-jetEnergy)/jetEnergy);
      		if(i>=2)hExpBJECIncl -> Fill((quarkEnergy-jetEnergy)/jetEnergy);
  		
					if(calculateJECinbins_){
						if(verbose_) cout << "start calculating expected jet energy corrections ... now in bins"<< endl;		
						// find eta and pt bin
      			int etabin  =  0;
      			if(etanrbins > 1){
        			for(unsigned int b=0; b<etabinVals_.size()-1; b++) {
          			if(fabs(myselectedjets[GenMatchQuarks->getMatchForParton(i,0)].eta()) > etabinVals_[b]) etabin = b;
        			}
      			}	 	    
      			int ptbin  =  0;
      			for(unsigned int b=0; b<pTbinVals_.size()-1; b++) {
        			if(myselectedjets[GenMatchQuarks->getMatchForParton(i,0)].pt() > pTbinVals_[b]) ptbin = b;
      			}
      			if(i<2)hExpLightJECPtEtaBin[etabin][ptbin] -> Fill((quarkEnergy-jetEnergy)/jetEnergy);
      			if(i>=2)hExpBJECPtEtaBin[etabin][ptbin] -> Fill((quarkEnergy-jetEnergy)/jetEnergy);
					}
					///////////////////// ---done--- calculating expected JEC /////////////////////
					if(verbose_) cout << "corresponding indexes of jets matched with quarks"<< endl;		
					int IdxQuarkJets = GenMatchQuarks->getMatchForParton(i,0);
					double jetpt = myselectedjets[IdxQuarkJets].pt();
					matchjetpt.push_back(jetpt);
					matchquarkjetidx.push_back(IdxQuarkJets);
				}
				
   			//the highest rank among the indexes of the jets stored (i.e. the ones matching with the partons), 
				//corresponds to the jet with the lowest pt 
   			if(matchquarkjetidx.size()!=0){
   				std::sort(matchquarkjetidx.begin(), matchquarkjetidx.end(), std::greater<int>() );
   				std::sort(matchjetpt.begin(), matchjetpt.end(), std::greater<double>() );
					LowestPtMatchedJet_Rank->Fill(matchquarkjetidx.front()+1);
   				LowestPtMatchedJet_Pt  ->Fill(matchjetpt.back());
   			}

				///////////////////// ---start --- calculating W and top masses /////////////////////
				if(verbose_) cout << "start calculating generated and reconstructed W and top masses"<< endl;		
				math::XYZTLorentzVector WHadGen = quarks[0]->p4()+quarks[1]->p4();
				math::XYZTLorentzVector WLepGen = muon->p4()+neutrino->p4();
				math::XYZTLorentzVector tHadGen = quarks[0]->p4()+quarks[1]->p4()+quarks[2]->p4();
				math::XYZTLorentzVector tLepGen = muon->p4()+neutrino->p4()+quarks[3]->p4();
				math::XYZTLorentzVector WHadRec = myselectedjets[GenMatchQuarks->getMatchForParton(0,0)].p4()+myselectedjets[GenMatchQuarks->getMatchForParton(1,0)].p4();
				math::XYZTLorentzVector WLepRec = muons[0].p4()+mets[0].p4();
				math::XYZTLorentzVector tHadRec = myselectedjets[GenMatchQuarks->getMatchForParton(0,0)].p4()+myselectedjets[GenMatchQuarks->getMatchForParton(1,0)].p4()+myselectedjets[GenMatchQuarks->getMatchForParton(2,0)].p4();
				math::XYZTLorentzVector tLepRec = muons[0].p4()+mets[0].p4()+myselectedjets[GenMatchQuarks->getMatchForParton(3,0)].p4();
				mWHad_Gen -> Fill(WHadGen.mass());
				mWLep_Gen -> Fill(WLepGen.mass());
				mtHad_Gen -> Fill(tHadGen.mass());
				mtLep_Gen -> Fill(tLepGen.mass());
				mWHad_Rec -> Fill(WHadRec.mass());
				mWLep_Rec -> Fill(WLepRec.mass());
				mtHad_Rec -> Fill(tHadRec.mass());
				mtLep_Rec -> Fill(tLepRec.mass());
				Dphi_HadTops -> Fill(ROOT::Math::VectorUtil::DeltaPhi(tHadGen,tHadRec));
				Dphi_LepTops -> Fill(ROOT::Math::VectorUtil::DeltaPhi(tLepGen,tLepRec));
				Deta_HadTops -> Fill(fabs(tHadGen.eta() - tHadRec.eta()));
				Deta_LepTops -> Fill(fabs(tLepGen.eta() - tLepRec.eta()));
				DR_HadTops -> Fill(ROOT::Math::VectorUtil::DeltaR(tHadGen,tHadRec));
				DR_LepTops -> Fill(ROOT::Math::VectorUtil::DeltaR(tHadGen,tHadRec));
				Dphi_Tops_Rec -> Fill(ROOT::Math::VectorUtil::DeltaPhi(tHadRec,tLepRec));
				Dphi_Tops_Gen -> Fill(ROOT::Math::VectorUtil::DeltaPhi(tHadGen,tLepGen));
				DR_Tops_Rec -> Fill(ROOT::Math::VectorUtil::DeltaR(tHadRec,tLepRec));
				DR_Tops_Gen -> Fill(ROOT::Math::VectorUtil::DeltaR(tHadGen,tLepGen));
				///////////////////// ---done --- calculating W and top masses /////////////////////

				///////////////////// ---start --- calculating smallest angle jets/mu and jets /////////////////////
				if(verbose_) cout << "start calculating smallest angle between jets/muon and jets (idem quarks/muon and quarks)"<< endl;		
				double smallestDRjetmu = 1000.; 
				double DRjetmu[4];
				double smallestDRquarkmu = 1000.; 
				double DRquarkmu[4];
				for(int i = 0; i<4; i++){
					DRjetmu[i] = deltaR(myselectedjets[GenMatchQuarks->getMatchForParton(i,0)].p4(),muons[0].p4());
					DRquarkmu[i] = deltaR(quarks[i]->p4(),muon->p4());
	  			if(DRjetmu[i] < smallestDRjetmu) smallestDRjetmu = DRjetmu[i];
	  			if(DRquarkmu[i] < smallestDRquarkmu) smallestDRquarkmu = DRquarkmu[i];
				}
				minJetMuVSminQuarkMu->Fill(smallestDRjetmu,smallestDRquarkmu);

				double smallestDRjets = 1000.; 
				double DRjets[6];
				double smallestDRquarks = 1000.; 
				double DRquarks[6];
				int k = 0;
				for(int i = 0; i<3; i++){
					for(int j = 3; j>i;j--){
						DRjets[k] = deltaR(myselectedjets[GenMatchQuarks->getMatchForParton(i,0)].p4(),myselectedjets[GenMatchQuarks->getMatchForParton(j,0)].p4());
						DRquarks[k] = deltaR(quarks[i]->p4(),quarks[j]->p4());
	  				if(DRjets[k] < smallestDRjets) smallestDRjets = DRjets[k];
	  				if(DRquarks[k] < smallestDRquarks) smallestDRquarks = DRquarks[k];
						k=k+1;
					}
				}
				minJetsVSminQuarks->Fill(smallestDRjets,smallestDRquarks);
				///////////////////// ---done --- calculating smallest angle jets/mu and jets /////////////////////


				/////////////////////////////////
				////////// LOOK AT ISR //////////
				/////////////////////////////////		
				if(verbose_) cout << "check if there is ISR"<< endl;		
				if( NumberISR !=0){
					if(verbose_) cout << "there are partons coming from ISR"<< endl;		
					//make new list of jets to be used as input for matching (remove quark-matched jets from jetlist used as input for quark matching)
					//make this new list of jets, but also keep the original rank of the jet! ---> this rank can be used for plots and numbers
					double arrayISR[100];
					for(unsigned int i=0; i<myselectedjets.size(); i++ ){
						for(unsigned int k =0; k<quarks.size(); k++){
							if((unsigned int) matchquarkjetidx[k] !=i) {
								arrayISR[i]=0;
							}else{
								arrayISR[i]=1;
							}
						}
					}
					for(unsigned int a=0; a<myselectedjets.size();a++){
						if(arrayISR[a]==0){
							myselectedjetsISR.push_back(myselectedjets[a]);
							originalrankjetsISRinput.push_back(a);
						}
					}

					if(verbose_) cout << "check how many partons from ISR compared to the number of unmatched jets"<< endl;							
				  ///////////////////// #ISR > myselectedjetsISR.size() /////////////////////
					if(NumberISR > myselectedjetsISR.size()){
						if(verbose_) cout << "---> more partons than jets: match highest pT partons"<< endl;							
						//match highest pT ISR parton(s) with jets
						for(unsigned int i=0; i<vectorISR.size();i++){ISRpT.push_back(vectorISR[i]->pt());}
						// 1. sort ISRpt on descending pt
						std::sort(ISRpT.begin(),ISRpT.end(),greater<double>());
						// 2. fill vector with the n first pt partons
						for(unsigned int n=0;n<myselectedjetsISR.size();n++){
							for(unsigned int k=0;k<NumberISR;k++){
								if(ISRpT[n]==vectorISR[k]->pt())NfirstpTISR.push_back(vectorISR[k]);
							}
						}
						// 3. do JetPartonMatching using only the highest pT ISR partons						
						if(verbose_) cout << "do matching"<< endl;							
						JetPartonMatching *GenMatchISR = new JetPartonMatching(NfirstpTISR, myselectedjetsISR, matchingAlgo_, useMaxDist_, useDeltaR_, maxDist_);
					
						//check if all partons are matched ---> if so, just store jet ranks, nothing else (no jets left)
						if(GenMatchISR->getNumberOfUnmatchedPartons(0)== 0){
							QuarksMatchedMatchedISR++;				
							if(verbose_) cout << "all partons coming from ISR are matched"<< endl;							
							
							//store the rank of the jet for later use (if parton is unmatched ---> rank is negative)
							if(verbose_) cout << "store the rank of the jets matched with the partons"<< endl;							
 							for(unsigned int i = 0; i < NfirstpTISR.size(); i++){			
								int IdxISRJets = GenMatchISR->getMatchForParton(i,0);
								matchISRjetidx.push_back(IdxISRJets);
								if(IdxISRJets>=0) originalrankjetsISR.push_back(originalrankjetsISRinput[IdxISRJets]);
								else originalrankjetsISR.push_back(IdxISRJets);
							}
							
						//check if unmatched partons ---> if so, categorize, store jet ranks, nothing else (no jets left)
						}else if(GenMatchISR->getNumberOfUnmatchedPartons(0)!= 0){ //unmatched quarks
							QuarksMatchedUnmatchedISR++;
							if(verbose_) cout << "not all partons coming from ISR are matched"<< endl;							
						
							//find parton which is not matched -> check eta/pT/DR ---> count different possibilities
							for(unsigned int i = 0; i<NfirstpTISR.size(); i++){
								if(verbose_) cout << "check which ISR partons are not matched and categorize"<< endl;							
								int Idx = GenMatchISR->getMatchForParton(i,0);
								if(Idx<0){
									QuarksMatchedUnmatchedISRtotal++;
									smallestDRISRi = 1000.;
									for(unsigned int j = 0; j<NfirstpTISR.size(); j++){
										DRISRi[j] = 1000.;
										if(j!=i){
											DRISRi[j] = deltaR(NfirstpTISR[i]->p4(),NfirstpTISR[j]->p4());
	  									if(DRISRi[j] < smallestDRISRi) smallestDRISRi = DRISRi[j];
										}
									}
									if(fabs(NfirstpTISR[i]->eta())>jetEta_){
										QuarksMatchedUnmatchedISRCategoryEta++;
									}else if(NfirstpTISR[i]->pt()<jetPt_){
										QuarksMatchedUnmatchedISRCategoryPt++;
									}else if(smallestDRISRi>minDRunmatchedParton_){
										QuarksMatchedUnmatchedISRCategoryDR++;
									}else{
										QuarksMatchedUnmatchedISRNotCategorized++;
									}
								}
							}
							
 							//store the rank of the jet for later use (if parton is unmatched ---> rank is negative)
							if(verbose_) cout << "store the rank of the jets matched with the partons"<< endl;							
 							for(unsigned int i = 0; i < NfirstpTISR.size(); i++){			
								int IdxISRJets = GenMatchISR->getMatchForParton(i,0);
								matchISRjetidx.push_back(IdxISRJets);
								if(IdxISRJets>=0) originalrankjetsISR.push_back(originalrankjetsISRinput[IdxISRJets]);
								else originalrankjetsISR.push_back(IdxISRJets);
							}
						}	
						
				  ///////////////////// #ISR <= myselectedjetsISR.size() /////////////////////					
          }else if( NumberISR <= myselectedjetsISR.size()){
						if(verbose_) cout << "less or as much partons from ISR compared to the number of jets"<< endl;							
						JetPartonMatching *GenMatchISR = new JetPartonMatching(vectorISR, myselectedjetsISR, matchingAlgo_, useMaxDist_, useDeltaR_, maxDist_);
						
						//check if unmatched ISR  ---> if so, categorize, store jet ranks, look at radiation (still jets left)
						if(GenMatchISR->getNumberOfUnmatchedPartons(0)!= 0){
							if(verbose_) cout << "not all partons are matched"<< endl;							
							QuarksMatchedUnmatchedISR++; 

							//find parton which is not matched -> check eta/pT/DR ---> count different possibilities
							for(unsigned int i = 0; i<NumberISR; i++){
								if(verbose_) cout << "check which ISR partons are not matched and categorize"<< endl;							
								int Idx = GenMatchISR->getMatchForParton(i,0);
								if(Idx<0){
									QuarksMatchedUnmatchedISRtotal++;
									smallestDRISRi = 1000.;
									for(unsigned int j = 0; j<NumberISR; j++){
										DRISRi[j] = 1000.;
										if(j!=i){
											DRISRi[j] = deltaR(vectorISR[i]->p4(),vectorISR[j]->p4());
	  									if(DRISRi[j] < smallestDRISRi) smallestDRISRi = DRISRi[j];
										}
									}
									if(fabs(vectorISR[i]->eta())>jetEta_){
										QuarksMatchedUnmatchedISRCategoryEta++;
									}else if(vectorISR[i]->pt()<jetPt_){
										QuarksMatchedUnmatchedISRCategoryPt++;
									}else if(smallestDRISRi>minDRunmatchedParton_){
										QuarksMatchedUnmatchedISRCategoryDR++;
									}else{
										QuarksMatchedUnmatchedISRNotCategorized++;
									}
								}
							}

				
							//store the rank of the jet for later use (if parton is unmatched ---> rank is negative)
							if(verbose_) cout << "store the rank of the jets matched with the partons"<< endl;							
 							for(unsigned int i = 0; i < NumberISR; i++){			
								int IdxISRJets = GenMatchISR->getMatchForParton(i,0);
								matchISRjetidx.push_back(IdxISRJets);
								if(IdxISRJets >=0) originalrankjetsISR.push_back(originalrankjetsISRinput[IdxISRJets]);
								else originalrankjetsISR.push_back(IdxISRJets);
							}


							//////////////////////////////////////////////////////////////////////////////
							////////// LOOK AT TOPRADIATION - JETS LEFT - ISR NOT ALL MATCHED - //////////
							//////////////////////////////////////////////////////////////////////////////
							if(verbose_) cout << "check if there is radiation from the top quarks"<< endl;		
							if(NumberRadiation !=0){ 
								if(verbose_) cout << "there is radiation from the top quarks"<< endl;		
								QuarksMatchedISRTopRadiation++;
								//make new list of jets to be used as input for matching (remove ISR-matched jets from jetlist used as input for ISR)
								double arrayRad[100];
								for(unsigned int i = 0; i<myselectedjets.size() ; i++){
									for(unsigned int l =0; l<quarks.size(); l++){
										for(unsigned int k =0; k<NumberISR; k++){
											if( (unsigned int) originalrankjetsISR[k]!=i && (unsigned int) matchquarkjetidx[l]!=i){
												arrayRad[i]=0;
											}else{
												arrayRad[i]=1;
											}
										}
									}
								}
								for(unsigned int a=0; a<myselectedjets.size();a++){
									if(arrayRad[a]==0){
										myselectedjetsradiation.push_back(myselectedjets[a]);
										originalrankjetsradiationinput.push_back(a);
									}
								}
								
				  			///////////////////// #top radiation > myselectedjetsradiation.size() /////////////////////
								if(verbose_) cout << "check the number of partons coming from top radiation compared to the number of jets"<< endl;		
								if(NumberRadiation > myselectedjetsradiation.size()){
									if(verbose_) cout << "more partons than jets"<< endl;		
									//match highest pT top radiation parton(s) with jets
									for(unsigned int i=0; i<TopRadiation.size();i++){TopRadiationpT.push_back(TopRadiation[i]->pt());}
									// 1. sort radiationpt on descending pt
									std::sort(TopRadiationpT.begin(),TopRadiationpT.end(),greater<double>());
									// 2. fill vector with the n first pt partons
									for(unsigned int n=0;n<myselectedjetsradiation.size();n++){
										for(unsigned int k=0;k<NumberRadiation;k++){
											if(TopRadiationpT[n]==TopRadiation[k]->pt())NfirstpTradiation.push_back(TopRadiation[k]);
										}
									}
									// 3. do JetPartonMatching using only the highest pT top radiation partons						
									if(verbose_) cout << "do matching"<< endl;		
									JetPartonMatching *GenMatchRadiation = new JetPartonMatching(NfirstpTradiation, myselectedjetsradiation, matchingAlgo_, useMaxDist_, useDeltaR_, maxDist_);

									//check if all partons are matched ---> if so, just store jet ranks, nothing else (no jets left)
									if(GenMatchRadiation->getNumberOfUnmatchedPartons(0)== 0){
										if(verbose_) cout << "all partons are matched"<< endl;		
										QuarksMatchedMatchedTopRadiation++;				
						
										//store the rank of the jet for later use (if parton is unmatched ---> rank is negative)
										if(verbose_) cout << "store the ranks of the jets matched to the partons"<< endl;		
 										for(unsigned int i = 0; i < NfirstpTradiation.size(); i++){			
											int IdxRadJets = GenMatchRadiation->getMatchForParton(i,0);
											matchRadiationjetidx.push_back(IdxRadJets);
											if(IdxRadJets>=0) originalrankjetsradiation.push_back(originalrankjetsradiationinput[IdxRadJets]);
											else originalrankjetsradiation.push_back(IdxRadJets);
										}
							
									//check if unmatched partons ---> if so, categorize, store jet ranks, nothing else (no jets left)
									}else if(GenMatchRadiation->getNumberOfUnmatchedPartons(0)!= 0){ 
										if(verbose_) cout << "not all partons matched"<< endl;		
										QuarksMatchedUnmatchedTopRadiation++;
					
										//find parton which is not matched -> check eta/pT/DR ---> count different possibilities
										for(unsigned int i = 0; i<NfirstpTradiation.size(); i++){
											if(verbose_) cout << "find the partons which are not matched and categorize"<< endl;		
											int Idx = GenMatchRadiation->getMatchForParton(i,0);
											if(Idx<0){
												QuarksMatchedUnmatchedTopRadiationtotal++;
												smallestDRRadi = 1000.;
												for(unsigned int j = 0; j<NfirstpTradiation.size(); j++){
													DRRadi[j] = 1000.;
													if(j!=i){
														DRRadi[j] = deltaR(NfirstpTradiation[i]->p4(),NfirstpTradiation[j]->p4());
	  												if(DRRadi[j] < smallestDRRadi) smallestDRRadi = DRRadi[j];
													}
												}
												if(fabs(NfirstpTradiation[i]->eta())>jetEta_){
													QuarksMatchedUnmatchedRadCategoryEta++;
												}else if(NfirstpTradiation[i]->pt()<jetPt_){
													QuarksMatchedUnmatchedRadCategoryPt++;
												}else if(smallestDRRadi>minDRunmatchedParton_){
													QuarksMatchedUnmatchedRadCategoryDR++;
												}else{
													QuarksMatchedUnmatchedRadNotCategorized++;
												}
											}
										}
							
 										//store the rank of the jet for later use (if parton is unmatched ---> rank is negative)
 										if(verbose_) cout << "store the ranks of the jets matched top the partons"<< endl;		
										for(unsigned int i = 0; i < NfirstpTradiation.size(); i++){			
											int IdxRadJets = GenMatchRadiation->getMatchForParton(i,0);
											matchRadiationjetidx.push_back(IdxRadJets);
											if(IdxRadJets>=0)originalrankjetsradiation.push_back(originalrankjetsradiationinput[IdxRadJets]);
											else originalrankjetsradiation.push_back(IdxRadJets);
										}
									}

				  			///////////////////// #top radiation <= myselectedjetsradiation.size() /////////////////////
          			}else if( NumberRadiation <= myselectedjetsradiation.size()){
 									if(verbose_) cout << "less or equal number of partons compared to the number of jets"<< endl;		
									JetPartonMatching *GenMatchRadiation = new JetPartonMatching(TopRadiation, myselectedjetsradiation, matchingAlgo_, useMaxDist_, useDeltaR_, maxDist_);
									
									//check if unmatched top radiation  ---> if so, categorize, store jet ranks 
									if(GenMatchRadiation->getNumberOfUnmatchedPartons(0)!= 0){
 										if(verbose_) cout << "not all partons matched"<< endl;		
										QuarksMatchedUnmatchedTopRadiation++;

										//find parton which is not matched -> check eta/pT/DR ---> count different possibilities
										for(unsigned int i = 0; i<NumberRadiation; i++){
 											if(verbose_) cout << "check which partons are not matched and categorize"<< endl;		
											int Idx = GenMatchRadiation->getMatchForParton(i,0);
											if(Idx<0){
												QuarksMatchedUnmatchedTopRadiationtotal++;
												smallestDRRadi = 1000.;
												for(unsigned int j = 0; j<NumberRadiation; j++){
													DRRadi[j] = 1000.;
													if(j!=i){
														DRRadi[j] = deltaR(TopRadiation[i]->p4(),TopRadiation[j]->p4());
	  												if(DRRadi[j] < smallestDRRadi) smallestDRRadi = DRRadi[j];
													}
												}
												if(fabs(TopRadiation[i]->eta())>jetEta_){
													QuarksMatchedUnmatchedRadCategoryEta++;
												}else if(TopRadiation[i]->pt()<jetPt_){
													QuarksMatchedUnmatchedRadCategoryPt++;
												}else if(smallestDRRadi>minDRunmatchedParton_){
													QuarksMatchedUnmatchedRadCategoryDR++;
												}else{
													QuarksMatchedUnmatchedRadNotCategorized++;
												}
											}
										}
 							
										//store the rank of the jet for later use (if parton is unmatched ---> rank is negative)
 										if(verbose_) cout << "store the rank of the jets matched with the partons"<< endl;		
 										for(unsigned int i = 0; i < NumberRadiation; i++){			
											int IdxRadJets = GenMatchRadiation->getMatchForParton(i,0);
											matchRadiationjetidx.push_back(IdxRadJets);
											if(IdxRadJets>=0)originalrankjetsradiation.push_back(originalrankjetsradiationinput[IdxRadJets]);
											else originalrankjetsradiation.push_back(IdxRadJets);
										}
							
									//check if all top radiation partons are matched ---> if so, store jet ranks 
									}else if(GenMatchRadiation->getNumberOfUnmatchedPartons(0)== 0){
 										if(verbose_) cout << "all partons matched"<< endl;		
										QuarksMatchedMatchedTopRadiation++;
										
										for(unsigned int i = 0; i < NumberRadiation; i++){			
											int IdxRadJets = GenMatchRadiation->getMatchForParton(i,0);
											matchRadiationjetidx.push_back(IdxRadJets);
											if(IdxRadJets>=0)originalrankjetsradiation.push_back(originalrankjetsradiationinput[IdxRadJets]);
											else originalrankjetsradiation.push_back(IdxRadJets);
										}
									}
								}
								
							}else{
 								if(verbose_) cout << "not radiation from the top quarks"<< endl;		
								QuarksMatchedNoTopRadiation++;
							}
			
						//check if all ISR partons are matched ---> if so, store jet ranks, look at radiation (still jets left)
						}else if(GenMatchISR->getNumberOfUnmatchedPartons(0)== 0){ 
 							if(verbose_) cout << "all ISR partons are matched"<< endl;									
							QuarksMatchedMatchedISR++;

							for(unsigned int i = 0; i < NumberISR; i++){			
								int IdxISRjets = GenMatchISR->getMatchForParton(i,0);
								matchISRjetidx.push_back(IdxISRjets);
								if(IdxISRjets>=0) originalrankjetsISR.push_back(originalrankjetsISRinput[IdxISRjets]);
								else originalrankjetsISR.push_back(IdxISRjets);
							}
				
							///////////////////////////////////////////////////////////////////////////
							////////// LOOK AT TOPRADIATION  - JETS LEFT - ISR ALL MATCHED - //////////
							///////////////////////////////////////////////////////////////////////////
 							if(verbose_) cout << "check if there are partons coming from radiation of top quarks"<< endl;									
							if(NumberRadiation !=0){ 
 								if(verbose_) cout << "there are partons coming from the radiation of top quarks"<< endl;									
								QuarksMatchedISRTopRadiation++;
								//make new list of jets to be used as input for matching (remove ISR-matched jets from jetlist used as input for ISR)
								double arrayRad[100];
								for(unsigned int i = 0; i<myselectedjets.size() ; i++){
									for(unsigned int l =0; l<quarks.size(); l++){
										for(unsigned int k =0; k<NumberISR; k++){
											if( (unsigned int) originalrankjetsISR[k]!=i && (unsigned int) matchquarkjetidx[l]!=i){
												arrayRad[i]=0;
											}else{
												arrayRad[i]=1;
											}
										}
									}
								}
								for(unsigned int a=0; a<myselectedjets.size();a++){
									if(arrayRad[a]==0){
										myselectedjetsradiation.push_back(myselectedjets[a]);
										originalrankjetsradiationinput.push_back(a);
									}
								}
								
 								if(verbose_) cout << "check the number of partons compared to the number of jets"<< endl;									
				  			///////////////////// #top radiation > myselectedjetsradiation.size() /////////////////////
								if(NumberRadiation > myselectedjetsradiation.size()){
 									if(verbose_) cout << "more partons compared to the number of jets"<< endl;									
									//match highest pT top radiation parton(s) with jets
									for(unsigned int i=0; i<TopRadiation.size();i++){TopRadiationpT.push_back(TopRadiation[i]->pt());}
									// 1. sort radiationpt on descending pt
									std::sort(TopRadiationpT.begin(),TopRadiationpT.end(),greater<double>());
									// 2. fill vector with the n first pt partons
									for(unsigned int n=0;n<myselectedjetsradiation.size();n++){
										for(unsigned int k=0;k<NumberRadiation;k++){
											if(TopRadiationpT[n]==TopRadiation[k]->pt())NfirstpTradiation.push_back(TopRadiation[k]);
										}
									}
									// 3. do JetPartonMatching using only the highest pT top radiation partons						
 									if(verbose_) cout << "do matching"<< endl;									
									JetPartonMatching *GenMatchRadiation = new JetPartonMatching(NfirstpTradiation, myselectedjetsradiation, matchingAlgo_, useMaxDist_, useDeltaR_, maxDist_);

									//check if all partons are matched ---> if so, just store jet ranks, nothing else (no jets left)
									if(GenMatchRadiation->getNumberOfUnmatchedPartons(0)== 0){
										QuarksMatchedMatchedTopRadiation++;				
							
 										if(verbose_) cout << "store the rank of the jets matched to the partons"<< endl;									
										//store the rank of the jet for later use (if parton is unmatched ---> rank is negative)
 										for(unsigned int i = 0; i < NfirstpTradiation.size(); i++){			
											int IdxRadJets = GenMatchRadiation->getMatchForParton(i,0);
											matchRadiationjetidx.push_back(IdxRadJets);
											if(IdxRadJets>=0)originalrankjetsradiation.push_back(originalrankjetsradiationinput[IdxRadJets]);
											else originalrankjetsradiation.push_back(IdxRadJets);
										}
							
 									if(verbose_) cout << "check if all partons are matched to the jets"<< endl;									
									//check if unmatched partons ---> if so, categorize, store jet ranks, nothing else (no jets left)
									}else if(GenMatchRadiation->getNumberOfUnmatchedPartons(0)!= 0){ 
 										if(verbose_) cout << "not all partons are matched to the quarks"<< endl;									
										QuarksMatchedUnmatchedTopRadiation++;
						
										//find parton which is not matched -> check eta/pT/DR ---> count different possibilities
										for(unsigned int i = 0; i<NfirstpTradiation.size(); i++){
 											if(verbose_) cout << "check which partons are not matched and categorize"<< endl;									
											int Idx = GenMatchRadiation->getMatchForParton(i,0);
											if(Idx<0){
												QuarksMatchedUnmatchedTopRadiationtotal++;
												smallestDRRadi = 1000.;
												for(unsigned int j = 0; j<NfirstpTradiation.size(); j++){
													DRRadi[j] = 1000.;
													if(j!=i){
														DRRadi[j] = deltaR(NfirstpTradiation[i]->p4(),NfirstpTradiation[j]->p4());
	  												if(DRRadi[j] < smallestDRRadi) smallestDRRadi = DRRadi[j];
													}
												}
												if(fabs(NfirstpTradiation[i]->eta())>jetEta_){
													QuarksMatchedUnmatchedRadCategoryEta++;
												}else if(NfirstpTradiation[i]->pt()<jetPt_){
													QuarksMatchedUnmatchedRadCategoryPt++;
												}else if(smallestDRRadi>minDRunmatchedParton_){
													QuarksMatchedUnmatchedRadCategoryDR++;
												}else{
													QuarksMatchedUnmatchedRadNotCategorized++;
												}
											}
										}
							
 										if(verbose_) cout << "store the rank of the jets matched with the partons"<< endl;									
 										//store the rank of the jet for later use (if parton is unmatched ---> rank is negative)
 										for(unsigned int i = 0; i < NfirstpTradiation.size(); i++){			
											int IdxRadJets = GenMatchRadiation->getMatchForParton(i,0);
											matchRadiationjetidx.push_back(IdxRadJets);
											if(IdxRadJets>=0)originalrankjetsradiation.push_back(originalrankjetsradiationinput[IdxRadJets]);
											else originalrankjetsradiation.push_back(IdxRadJets);
										}
									}

				  			///////////////////// #top radiation <= myselectedjetsradiation.size() /////////////////////
          			}else if( NumberRadiation <= myselectedjetsradiation.size()){
 									if(verbose_) cout << "less or equal number of partons compared to the number of jets"<< endl;									
									JetPartonMatching *GenMatchRadiation = new JetPartonMatching(TopRadiation, myselectedjetsradiation, matchingAlgo_, useMaxDist_, useDeltaR_, maxDist_);
									
									//check if unmatched top radiation  ---> if so, categorize, store jet ranks 
									if(GenMatchRadiation->getNumberOfUnmatchedPartons(0)!= 0){
 										if(verbose_) cout << "not all partons are matched"<< endl;									
										QuarksMatchedUnmatchedTopRadiation++;

										//find parton which is not matched -> check eta/pT/DR ---> count different possibilities
										for(unsigned int i = 0; i<NumberRadiation; i++){
 											if(verbose_) cout << "find partons which are not matched and categorize"<< endl;									
											int Idx = GenMatchRadiation->getMatchForParton(i,0);
											if(Idx<0){
												QuarksMatchedUnmatchedTopRadiationtotal++;
												smallestDRRadi = 1000.;
												for(unsigned int j = 0; j<NumberRadiation; j++){
													DRRadi[j] = 1000.;
													if(j!=i){
														DRRadi[j] = deltaR(TopRadiation[i]->p4(),TopRadiation[j]->p4());
	  												if(DRRadi[j] < smallestDRRadi) smallestDRRadi = DRRadi[j];
													}
												}
												if(fabs(TopRadiation[i]->eta())>jetEta_){
													QuarksMatchedUnmatchedRadCategoryEta++;
												}else if(TopRadiation[i]->pt()<jetPt_){
													QuarksMatchedUnmatchedRadCategoryPt++;
												}else if(smallestDRRadi>minDRunmatchedParton_){
													QuarksMatchedUnmatchedRadCategoryDR++;
												}else{
													QuarksMatchedUnmatchedRadNotCategorized++;
												}
											}
										}
 							
 										if(verbose_) cout << "store the ranks of the jets matched to the partons"<< endl;									
										//store the rank of the jet for later use (if parton is unmatched ---> rank is negative)
 										for(unsigned int i = 0; i < NumberRadiation; i++){			
											int IdxRadJets = GenMatchRadiation->getMatchForParton(i,0);
											matchRadiationjetidx.push_back(IdxRadJets);
											if(IdxRadJets>=0)originalrankjetsradiation.push_back(originalrankjetsradiationinput[IdxRadJets]);
											else originalrankjetsradiation.push_back(IdxRadJets);
										}
							
									//check if all top radiation partons are matched ---> if so, store jet ranks 
									}else if(GenMatchRadiation->getNumberOfUnmatchedPartons(0)== 0){
 										if(verbose_) cout << "all partons are matched"<< endl;									
										QuarksMatchedMatchedTopRadiation++;
										
										for(unsigned int i = 0; i < NumberRadiation; i++){			
											int IdxRadJets = GenMatchRadiation->getMatchForParton(i,0);
											matchRadiationjetidx.push_back(IdxRadJets);
											if(IdxRadJets>=0)originalrankjetsradiation.push_back(originalrankjetsradiationinput[IdxRadJets]);
											else originalrankjetsradiation.push_back(IdxRadJets);
										}
									}
								}
								
							}else{
 								if(verbose_) cout << "no partons from top radiation"<< endl;									
								QuarksMatchedNoTopRadiation++;
							}
						}
					}					
				
				///////////////////////////////////////////////////
				////////// LOOK AT TOPRADIATION - NO ISR //////////
				///////////////////////////////////////////////////	
				//no ISR partons (either not a MadGraph sample or really no ISR happened)	
				}else{ 
 					if(verbose_) cout << "no partons from ISR"<< endl;									
					QuarksMatchedNoISR++;
					
					if(NumberRadiation !=0){ 
 						if(verbose_) cout << "there are partons coming from radiation of top quarks"<< endl;									
						QuarksMatchedTopRadiationNoISR++;
						//make new list of jets to be used as input for matching (remove jets matched with quarks from jetlist (input for	radiation)
						double array[100];
						for(unsigned int i=0; i<myselectedjets.size(); i++ ){
							for(unsigned int k =0; k<quarks.size(); k++){
								if((unsigned int) matchquarkjetidx[k] !=i) {
									array[i]=0;
								}else{
									array[i]=1;
								}
							}
						}
						for(unsigned int a=0; a<myselectedjets.size();a++){
							if(array[a]==0){
								myselectedjetsradiation.push_back(myselectedjets[a]);
								originalrankjetsradiationinput.push_back(a);
							}
						}
						
								
 						if(verbose_) cout << "check the number of partons compared to the number of jets"<< endl;									
						///////////////////// #top radiation > myselectedjetsradiation.size() /////////////////////
						if(NumberRadiation > myselectedjetsradiation.size()){
 							if(verbose_) cout << "more partons than jets"<< endl;									
							//match highest pT top radiation parton(s) with jets
							for(unsigned int i=0; i<TopRadiation.size();i++){TopRadiationpT.push_back(TopRadiation[i]->pt());}
							// 1. sort radiationpt on descending pt
							std::sort(TopRadiationpT.begin(),TopRadiationpT.end(),greater<double>());
							// 2. fill vector with the n first pt partons
							for(unsigned int n=0;n<myselectedjetsradiation.size();n++){
								for(unsigned int k=0;k<NumberRadiation;k++){
									if(TopRadiationpT[n]==TopRadiation[k]->pt())NfirstpTradiation.push_back(TopRadiation[k]);
								}
							}
 							if(verbose_) cout << "do matching"<< endl;									
							// 3. do JetPartonMatching using only the highest pT top radiation partons						
							JetPartonMatching *GenMatchRadiation = new JetPartonMatching(NfirstpTradiation, myselectedjetsradiation, matchingAlgo_, useMaxDist_, useDeltaR_, maxDist_);

							//check if all partons are matched ---> if so, just store jet ranks, nothing else (no jets left)
							if(GenMatchRadiation->getNumberOfUnmatchedPartons(0)== 0){
 								if(verbose_) cout << "all partons are matched"<< endl;									
								QuarksMatchedMatchedTopRadiation++;				
						
 								if(verbose_) cout << "store the ranks of the jets matched with the partons"<< endl;									
								//store the rank of the jet for later use (if parton is unmatched ---> rank is negative)
 								for(unsigned int i = 0; i < NfirstpTradiation.size(); i++){			
									int IdxRadJets = GenMatchRadiation->getMatchForParton(i,0);
									matchRadiationjetidx.push_back(IdxRadJets);
									if(IdxRadJets>=0)originalrankjetsradiation.push_back(originalrankjetsradiationinput[IdxRadJets]);
									else originalrankjetsradiation.push_back(IdxRadJets);
								}
							
							//check if unmatched partons ---> if so, categorize, store jet ranks, nothing else (no jets left)
							}else if(GenMatchRadiation->getNumberOfUnmatchedPartons(0)!= 0){ 
 								if(verbose_) cout << "not all partons are matched"<< endl;									
								QuarksMatchedUnmatchedTopRadiation++;
						
								//find parton which is not matched -> check eta/pT/DR ---> count different possibilities
								for(unsigned int i = 0; i<NfirstpTradiation.size(); i++){
 									if(verbose_) cout << "find partons which are not matched and categorize"<< endl;									
									int Idx = GenMatchRadiation->getMatchForParton(i,0);
									if(Idx<0){
										QuarksMatchedUnmatchedTopRadiationtotal++;
										smallestDRRadi = 1000.;
										for(unsigned int j = 0; j<NfirstpTradiation.size(); j++){
											DRRadi[j] = 1000.;
											if(j!=i){
												DRRadi[j] = deltaR(NfirstpTradiation[i]->p4(),NfirstpTradiation[j]->p4());
	  										if(DRRadi[j] < smallestDRRadi) smallestDRRadi = DRRadi[j];
											}
										}
										if(fabs(NfirstpTradiation[i]->eta())>jetEta_){
											QuarksMatchedUnmatchedRadCategoryEta++;
										}else if(NfirstpTradiation[i]->pt()<jetPt_){
											QuarksMatchedUnmatchedRadCategoryPt++;
										}else if(smallestDRRadi>minDRunmatchedParton_){
											QuarksMatchedUnmatchedRadCategoryDR++;
										}else{
											QuarksMatchedUnmatchedRadNotCategorized++;
										}
									}	
								}	
							
 								if(verbose_) cout << "store the ranks of the jets which are matched with the partons"<< endl;									
 								//store the rank of the jet for later use (if parton is unmatched ---> rank is negative)
 								for(unsigned int i = 0; i < NfirstpTradiation.size(); i++){			
									int IdxRadJets = GenMatchRadiation->getMatchForParton(i,0);
									matchRadiationjetidx.push_back(IdxRadJets);
									if(IdxRadJets>=0)originalrankjetsradiation.push_back(originalrankjetsradiationinput[IdxRadJets]);
									else originalrankjetsradiation.push_back(IdxRadJets);
								}
							}	

				  	///////////////////// #top radiation <= myselectedjetsradiation.size() /////////////////////
          	}else if( NumberRadiation <= myselectedjetsradiation.size()){
 							if(verbose_) cout << "less or equal number of partons compared to the number of jets"<< endl;									
							JetPartonMatching *GenMatchRadiation = new JetPartonMatching(TopRadiation, myselectedjetsradiation, matchingAlgo_, useMaxDist_, useDeltaR_, maxDist_);
							
							//check if unmatched top radiation  ---> if so, categorize, store jet ranks 
							if(GenMatchRadiation->getNumberOfUnmatchedPartons(0)!= 0){
 								if(verbose_) cout << "there are partons which are not matched"<< endl;									
								QuarksMatchedUnmatchedTopRadiation++;

								//find parton which is not matched -> check eta/pT/DR ---> count different possibilities
								for(unsigned int i = 0; i<NumberRadiation; i++){
 									if(verbose_) cout << "check which partons not matched and categorize"<< endl;									
									int Idx = GenMatchRadiation->getMatchForParton(i,0);
									if(Idx<0){
										QuarksMatchedUnmatchedTopRadiationtotal++;
										smallestDRRadi = 1000.;
										for(unsigned int j = 0; j<NumberRadiation; j++){
											DRRadi[j] = 1000.;
											if(j!=i){
												DRRadi[j] = deltaR(TopRadiation[i]->p4(),TopRadiation[j]->p4());
	  										if(DRRadi[j] < smallestDRRadi) smallestDRRadi = DRRadi[j];
											}
										}
										if(fabs(TopRadiation[i]->eta())>jetEta_){
											QuarksMatchedUnmatchedRadCategoryEta++;
										}else if(TopRadiation[i]->pt()<jetPt_){
											QuarksMatchedUnmatchedRadCategoryPt++;
										}else if(smallestDRRadi>minDRunmatchedParton_){
											QuarksMatchedUnmatchedRadCategoryDR++;
										}else{
											QuarksMatchedUnmatchedRadNotCategorized++;
										}
									}
								}
 							
 								if(verbose_) cout << "store the ranks of the jets matched with the partons"<< endl;									
								//store the rank of the jet for later use (if parton is unmatched ---> rank is negative)
 								for(unsigned int i = 0; i < NumberRadiation; i++){			
									int IdxRadJets = GenMatchRadiation->getMatchForParton(i,0);
									matchRadiationjetidx.push_back(IdxRadJets);
									if(IdxRadJets>=0)originalrankjetsradiation.push_back(originalrankjetsradiationinput[IdxRadJets]);
									else originalrankjetsradiation.push_back(IdxRadJets);
								}
						
							//check if all top radiation partons are matched ---> if so, store jet ranks
							}else if(GenMatchRadiation->getNumberOfUnmatchedPartons(0)== 0){
 								if(verbose_) cout << "all partons are matched"<< endl;									
								QuarksMatchedMatchedTopRadiation++;
								
								for(unsigned int i = 0; i < NumberRadiation; i++){			
									int IdxRadJets = GenMatchRadiation->getMatchForParton(i,0);
									matchRadiationjetidx.push_back(IdxRadJets);
									if(IdxRadJets>=0)originalrankjetsradiation.push_back(originalrankjetsradiationinput[IdxRadJets]);
									else originalrankjetsradiation.push_back(IdxRadJets);
								}
							}
						}
					}else{
 						if(verbose_) cout << "no radiation from the top quarks"<< endl;									
						QuarksMatchedNoRadiation++;
					}
				}

				////////////////////////////////////////////////////////////
				////////// MAKE PLOTS WHEN ALL QUARKS ARE MATCHED //////////
				////////////////////////////////////////////////////////////	
			
				bool Radiation = false;
 				if(verbose_) cout << "store the 4 highest pT jets"<< endl;									
				for(unsigned int j=0; j<myselectedjets.size();j++){ 
					//make vector with 4 highest pT jets+rank
					if(j<4) {
						FourHighestPtJets.push_back(myselectedjets[j]);
					}
				}
				
				//keep jets + rank which are matched with radiation!!!
				double arraynew[100];
				if(NumberISR !=0 && NumberRadiation !=0){
	 				if(verbose_) cout << "keep the jets matched to radiation (ISR and top) and their ranks"<< endl;									
					for(unsigned int j=0; j<myselectedjets.size();j++){ 
						for(unsigned int l =0; l<originalrankjetsradiation.size(); l++){
							if(j == (unsigned int) originalrankjetsradiation[l]){
								Radiation = true;
								arraynew[j]=0;
								l=originalrankjetsradiation.size();		
							}else{arraynew[j]=1;}
						}
						if(arraynew[j] != 0){
							for(unsigned int k =0; k<originalrankjetsISR.size(); k++){
								if(j == (unsigned int) originalrankjetsISR[k]){ 
									Radiation = true;
									arraynew[j]=0;
									k=originalrankjetsISR.size();					
								}else{arraynew[j]=1;}
							}
						}
					}
				}else if((NumberISR ==0) && (NumberRadiation !=0)){
	 				if(verbose_) cout << "keep the jets matched to radiation (top) and their ranks"<< endl;									
					for(unsigned int j=0; j<myselectedjets.size();j++){ 
						for(unsigned int l =0; l<originalrankjetsradiation.size(); l++){
							if(j == (unsigned int) originalrankjetsradiation[l]){
								Radiation = true;
								arraynew[j]=0;
								l=originalrankjetsradiation.size();
							}else{arraynew[j]=1;}
						}
					}
				}else if((NumberISR !=0) && (NumberRadiation ==0)){
	 				if(verbose_) cout << "keep the jets matched to radiation (ISR) and their ranks"<< endl;									
					for(unsigned int j=0; j<myselectedjets.size();j++){ 
						for(unsigned int k =0; k<originalrankjetsISR.size(); k++){
							if(j == (unsigned int) originalrankjetsISR[k] ){
								Radiation = true;
								arraynew[j]=0;
								k=originalrankjetsISR.size();
							}else{arraynew[j]=1;}
						}
					}
				}
				
				bool jetmatchedwithradiation = false;
				if(Radiation){
 					if(verbose_) cout << "store the jets and their pt"<< endl;									
					for(unsigned int j=0; j<myselectedjets.size();j++){ 
						if(arraynew[j]==0){
							jetmatchedwithradiation = true;
							pair<int,pat::Jet> pairjetsrad(j,myselectedjets[j]);
							patJetsMatchedWithRadiation.push_back(pairjetsrad);
						}
					}	
				}
			
				if(Radiation && jetmatchedwithradiation){
 					int countradmultiplicity = 0;
					countradmultiplicity = patJetsMatchedWithRadiation.size();
					QuarksMatchedJetRadiationMultiplicity -> Fill(countradmultiplicity);
					
					double arraynewnew[4];
					bool notgood = false;
					for(unsigned int i = 0; i<FourHighestPtJets.size();i++){		
						for(unsigned int j = 0; j<patJetsMatchedWithRadiation.size();j++){
							if(i == (unsigned int) patJetsMatchedWithRadiation[j].first){
								notgood = true;
								arraynewnew[i]=0;
								//if among the 4 highest pT jets, some are matched with radiation:
								//store the jet rank
								//store the jet pT
								rankjetsmatchedtorad.push_back(patJetsMatchedWithRadiation[j].first);
								ptjetsmatchedtorad.push_back((patJetsMatchedWithRadiation[j].second).pt());
								j = patJetsMatchedWithRadiation.size();
							}else{arraynewnew[i]=1;}
						}
					}
					
					int count = 0;
					if(notgood){
						//highest pT jet has per definition lowest rank!
						//sort first on rank --> make plot of lowest
						//sort on pT -->make plot of highest
						std::sort(rankjetsmatchedtorad.begin(), rankjetsmatchedtorad.end(), std::greater<int>() );
   					std::sort(ptjetsmatchedtorad.begin(), ptjetsmatchedtorad.end(), std::greater<double>() );
						QuarksMatchedHighestPtRadMatchedJet_Rank->Fill(rankjetsmatchedtorad.back()+1);
   					QuarksMatchedHighestPtRadMatchedJet_Pt  ->Fill(ptjetsmatchedtorad.front());
						std::sort(matchjetpt.begin(), matchjetpt.end(), std::greater<double>() );
						LowestPtQuarkMatchedJetPt_vs_HighestPtRadMatchedJetPt  ->Fill(matchjetpt.back(),ptjetsmatchedtorad.front()); 						
						
						if(verbose_) cout << "some of the jets are matched with radiation, count how many"<< endl;									
						for(unsigned int i = 0; i<FourHighestPtJets.size(); i++){
							if(arraynewnew[i]==0){
								count++;
							}
						}
					}
					
					if(count == 0)QuarksMatchedHighest4NoRad++;
					if(count >= 1)QuarksMatchedHighest4AtLeast1Rad++;
					if(count >= 2)QuarksMatchedHighest4AtLeast2Rad++;
					if(count >= 3)QuarksMatchedHighest4AtLeast3Rad++;
					if(count >= 4)QuarksMatchedHighest4AtLeast4Rad++;

				}else{
					edm::LogInfo("NoDataFound") << "all quarks matched, no radiation, so no related plots are filled \n";
				}








			
			//check if unmatched quarks ---> if so, categorize, look at radiation
			}else{ //unmatched quarks
				if(verbose_) cout << "not all quarks are matched"<< endl;		
				UnmatchedQuarks++;

				//find parton which is not matched -> check eta/pT/DR ---> count different possibilities
				int unmatchedquarkmultiplicity = 0;
				for(unsigned int i = 0; i<quarks.size(); i++){
					if(verbose_) cout << "check which quarks not matched and categorize unmatched quarks"<< endl;		
					int Idx = GenMatchQuarks->getMatchForParton(i,0);
					if(Idx<0){
						unmatchedquarkmultiplicity++;
						smallestDRquarki = 1000.;
						for(unsigned int j = 0; j<quarks.size(); j++){
							DRquarki[j] = 1000.;
							if(j!=i){
								DRquarki[j] = deltaR(quarks[i]->p4(),quarks[j]->p4());
	  						if(DRquarki[j] < smallestDRquarki) smallestDRquarki = DRquarki[j];
							}
						}
						if(fabs(quarks[i]->eta())>jetEta_){
							UnmatchedQuarkCategoryEta++;
						}else if(quarks[i]->pt()<jetPt_){
							UnmatchedQuarkCategoryPt++;
						}else if(smallestDRquarki>minDRunmatchedParton_){
							UnmatchedQuarkCategoryDR++;
						}else{
							UnmatchedQuarkNotCategorized++;
						}
					}
				}
				
				UnmatchedQuarkMultiplicity->Fill(unmatchedquarkmultiplicity);
				//find parton which is not matched -> check eta/pT/DR ---> count different possibilities
				int countnumberunmatched = 0;
				bool only1 = true;
				for(unsigned int i = 0; i<quarks.size(); i++){
					if(verbose_) cout << "check which quarks not matched and categorize unmatched quarks"<< endl;		
					int Idx = GenMatchQuarks->getMatchForParton(i,0);
					if(Idx<0)countnumberunmatched++;
					if(only1 && Idx<0){
						only1 = false;
						smallestDRquarki = 1000.;
						for(unsigned int j = 0; j<quarks.size(); j++){
							DRquarki[j] = 1000.;
							if(j!=i){
								DRquarki[j] = deltaR(quarks[i]->p4(),quarks[j]->p4());
	  						if(DRquarki[j] < smallestDRquarki) smallestDRquarki = DRquarki[j];
							}
						}
						if(fabs(quarks[i]->eta())>jetEta_){
							UnmatchedQuark1CategoryEta++;
						}else if(quarks[i]->pt()<jetPt_){
							UnmatchedQuark1CategoryPt++;
						}else if(smallestDRquarki>minDRunmatchedParton_){
							UnmatchedQuark1CategoryDR++;
						}else{
							UnmatchedQuark1NotCategorized++;
						}
					}
				}
				if(countnumberunmatched == 1){
					Only1UnmatchedQuarks++;
					for(unsigned int i = 0; i<quarks.size(); i++){
						if(verbose_) cout << "check which quarks not matched and categorize unmatched quarks"<< endl;		
						int Idx = GenMatchQuarks->getMatchForParton(i,0);
						if(Idx<0){
							smallestDRquarki = 1000.;
							for(unsigned int j = 0; j<quarks.size(); j++){
								DRquarki[j] = 1000.;
								if(j!=i){
									DRquarki[j] = deltaR(quarks[i]->p4(),quarks[j]->p4());
	  							if(DRquarki[j] < smallestDRquarki) smallestDRquarki = DRquarki[j];
								}
							}
							if(fabs(quarks[i]->eta())>jetEta_){
								Only1UnmatchedQuarkCategoryEta++;
							}else if(quarks[i]->pt()<jetPt_){
								Only1UnmatchedQuarkCategoryPt++;
							}else if(smallestDRquarki>minDRunmatchedParton_){
								Only1UnmatchedQuarkCategoryDR++;
							}else{
								Only1UnmatchedQuarkNotCategorized++;
							}
						}
					}
				}
				
				//store the rank of the jet (if parton is unmatched ---> rank is negative)
				if(verbose_) cout << "store the rank of the jets which are matched"<< endl;		
				for(unsigned int i = 0; i < quarks.size(); i++){			
					int IdxQuarkJets = GenMatchQuarks->getMatchForParton(i,0);
					matchquarkjetidx.push_back(IdxQuarkJets); 
				}
				
				/////////////////////////////////
				////////// LOOK AT ISR //////////
				/////////////////////////////////		
				if(verbose_) cout << "check if there is ISR"<< endl;		
				if( NumberISR !=0){
					if(verbose_) cout << "there are partons coming from ISR"<< endl;		
					//make new list of jets to be used as input for matching (remove quark-matched jets from jetlist used as input for quark matching)
					//make this new list of jets, but also keep the original rank of the jet! ---> this rank can be used for plots and numbers
					double arrayISR[100];
					for(unsigned int i=0; i<myselectedjets.size(); i++ ){
						for(unsigned int k =0; k<quarks.size(); k++){
							if((unsigned int) matchquarkjetidx[k] !=i) {
								arrayISR[i]=0;
							}else{
								arrayISR[i]=1;
							}
						}
					}
					for(unsigned int a=0; a<myselectedjets.size();a++){
						if(arrayISR[a]==0){
							myselectedjetsISR.push_back(myselectedjets[a]);
							originalrankjetsISRinput.push_back(a);
						}
					}

					if(verbose_) cout << "check how many partons from ISR compared to the number of unmatched jets"<< endl;							
				  ///////////////////// #ISR > myselectedjetsISR.size() /////////////////////
					if(NumberISR > myselectedjetsISR.size()){
						if(verbose_) cout << "---> more partons than jets: match highest pT partons"<< endl;							
						//match highest pT ISR parton(s) with jets
						for(unsigned int i=0; i<vectorISR.size();i++){ISRpT.push_back(vectorISR[i]->pt());}
						// 1. sort ISRpt on descending pt
						std::sort(ISRpT.begin(),ISRpT.end(),greater<double>());
						// 2. fill vector with the n first pt partons
						for(unsigned int n=0;n<myselectedjetsISR.size();n++){
							for(unsigned int k=0;k<NumberISR;k++){
								if(ISRpT[n]==vectorISR[k]->pt())NfirstpTISR.push_back(vectorISR[k]);
							}
						}
						// 3. do JetPartonMatching using only the highest pT ISR partons						
						if(verbose_) cout << "do matching"<< endl;							
						JetPartonMatching *GenMatchISR = new JetPartonMatching(NfirstpTISR, myselectedjetsISR, matchingAlgo_, useMaxDist_, useDeltaR_, maxDist_);
					
						//check if all partons are matched ---> if so, just store jet ranks, nothing else (no jets left)
						if(GenMatchISR->getNumberOfUnmatchedPartons(0)== 0){
							MatchedISR++;				
							if(verbose_) cout << "all partons coming from ISR are matched"<< endl;							
							
							//store the rank of the jet for later use (if parton is unmatched ---> rank is negative)
							if(verbose_) cout << "store the rank of the jets matched with the partons"<< endl;							
 							for(unsigned int i = 0; i < NfirstpTISR.size(); i++){			
								int IdxISRJets = GenMatchISR->getMatchForParton(i,0);
								matchISRjetidx.push_back(IdxISRJets);
								if(IdxISRJets>=0) originalrankjetsISR.push_back(originalrankjetsISRinput[IdxISRJets]);
								else originalrankjetsISR.push_back(IdxISRJets);
							}
							
						//check if unmatched partons ---> if so, categorize, store jet ranks, nothing else (no jets left)
						}else if(GenMatchISR->getNumberOfUnmatchedPartons(0)!= 0){ //unmatched quarks
							UnmatchedISR++;
							if(verbose_) cout << "not all partons coming from ISR are matched"<< endl;							
						
							//find parton which is not matched -> check eta/pT/DR ---> count different possibilities
							for(unsigned int i = 0; i<NfirstpTISR.size(); i++){
								if(verbose_) cout << "check which ISR partons are not matched and categorize"<< endl;							
								int Idx = GenMatchISR->getMatchForParton(i,0);
								if(Idx<0){
									smallestDRISRi = 1000.;
									for(unsigned int j = 0; j<NfirstpTISR.size(); j++){
										DRISRi[j] = 1000.;
										if(j!=i){
											DRISRi[j] = deltaR(NfirstpTISR[i]->p4(),NfirstpTISR[j]->p4());
	  									if(DRISRi[j] < smallestDRISRi) smallestDRISRi = DRISRi[j];
										}
									}
									if(fabs(NfirstpTISR[i]->eta())>jetEta_){
										UnmatchedISRCategoryEta++;
									}else if(NfirstpTISR[i]->pt()<jetPt_){
										UnmatchedISRCategoryPt++;
									}else if(smallestDRISRi>minDRunmatchedParton_){
										UnmatchedISRCategoryDR++;
									}else{
										UnmatchedISRNotCategorized++;
									}
								}
							}
							
 							//store the rank of the jet for later use (if parton is unmatched ---> rank is negative)
							if(verbose_) cout << "store the rank of the jets matched with the partons"<< endl;							
 							for(unsigned int i = 0; i < NfirstpTISR.size(); i++){			
								int IdxISRJets = GenMatchISR->getMatchForParton(i,0);
								matchISRjetidx.push_back(IdxISRJets);
								if(IdxISRJets>=0) originalrankjetsISR.push_back(originalrankjetsISRinput[IdxISRJets]);
								else originalrankjetsISR.push_back(IdxISRJets);
							}
						}	
						
				  ///////////////////// #ISR <= myselectedjetsISR.size() /////////////////////					
          }else if( NumberISR <= myselectedjetsISR.size()){
						if(verbose_) cout << "less or as much partons from ISR compared to the number of jets"<< endl;							
						JetPartonMatching *GenMatchISR = new JetPartonMatching(vectorISR, myselectedjetsISR, matchingAlgo_, useMaxDist_, useDeltaR_, maxDist_);
						
						//check if unmatched ISR  ---> if so, categorize, store jet ranks, look at radiation (still jets left)
						if(GenMatchISR->getNumberOfUnmatchedPartons(0)!= 0){
							if(verbose_) cout << "not all partons are matched"<< endl;							
							UnmatchedISR++; 

							//find parton which is not matched -> check eta/pT/DR ---> count different possibilities
							for(unsigned int i = 0; i<NumberISR; i++){
								if(verbose_) cout << "check which ISR partons are not matched and categorize"<< endl;							
								int Idx = GenMatchISR->getMatchForParton(i,0);
								if(Idx<0){
									smallestDRISRi = 1000.;
									for(unsigned int j = 0; j<NumberISR; j++){
										DRISRi[j] = 1000.;
										if(j!=i){
											DRISRi[j] = deltaR(vectorISR[i]->p4(),vectorISR[j]->p4());
	  									if(DRISRi[j] < smallestDRISRi) smallestDRISRi = DRISRi[j];
										}
									}
									if(fabs(vectorISR[i]->eta())>jetEta_){
										UnmatchedISRCategoryEta++;
									}else if(vectorISR[i]->pt()<jetPt_){
										UnmatchedISRCategoryPt++;
									}else if(smallestDRISRi>minDRunmatchedParton_){
										UnmatchedISRCategoryDR++;
									}else{
										UnmatchedISRNotCategorized++;
									}
								}
							}

				
							//store the rank of the jet for later use (if parton is unmatched ---> rank is negative)
							if(verbose_) cout << "store the rank of the jets matched with the partons"<< endl;							
 							for(unsigned int i = 0; i < NumberISR; i++){			
								int IdxISRJets = GenMatchISR->getMatchForParton(i,0);
								matchISRjetidx.push_back(IdxISRJets);
								if(IdxISRJets >=0) originalrankjetsISR.push_back(originalrankjetsISRinput[IdxISRJets]);
								else originalrankjetsISR.push_back(IdxISRJets);
							}


							//////////////////////////////////////////////////////////////////////////////
							////////// LOOK AT TOPRADIATION - JETS LEFT - ISR NOT ALL MATCHED - //////////
							//////////////////////////////////////////////////////////////////////////////
							if(verbose_) cout << "check if there is radiation from the top quarks"<< endl;		
							if(NumberRadiation !=0){ 
								if(verbose_) cout << "there is radiation from the top quarks"<< endl;		
								ISRTopRadiation++;
								//make new list of jets to be used as input for matching (remove ISR-matched jets from jetlist used as input for ISR)
								double arrayRad[100];
								for(unsigned int i = 0; i<myselectedjets.size() ; i++){
									for(unsigned int l =0; l<quarks.size(); l++){
										for(unsigned int k =0; k<NumberISR; k++){
											if( (unsigned int) originalrankjetsISR[k]!=i && (unsigned int) matchquarkjetidx[l]!=i){
												arrayRad[i]=0;
											}else{
												arrayRad[i]=1;
											}
										}
									}
								}
								for(unsigned int a=0; a<myselectedjets.size();a++){
									if(arrayRad[a]==0){
										myselectedjetsradiation.push_back(myselectedjets[a]);
										originalrankjetsradiationinput.push_back(a);
									}
								}
								
				  			///////////////////// #top radiation > myselectedjetsradiation.size() /////////////////////
								if(verbose_) cout << "check the number of partons coming from top radiation compared to the number of jets"<< endl;		
								if(NumberRadiation > myselectedjetsradiation.size()){
									if(verbose_) cout << "more partons than jets"<< endl;		
									//match highest pT top radiation parton(s) with jets
									for(unsigned int i=0; i<TopRadiation.size();i++){TopRadiationpT.push_back(TopRadiation[i]->pt());}
									// 1. sort radiationpt on descending pt
									std::sort(TopRadiationpT.begin(),TopRadiationpT.end(),greater<double>());
									// 2. fill vector with the n first pt partons
									for(unsigned int n=0;n<myselectedjetsradiation.size();n++){
										for(unsigned int k=0;k<NumberRadiation;k++){
											if(TopRadiationpT[n]==TopRadiation[k]->pt())NfirstpTradiation.push_back(TopRadiation[k]);
										}
									}
									// 3. do JetPartonMatching using only the highest pT top radiation partons						
									if(verbose_) cout << "do matching"<< endl;		
									JetPartonMatching *GenMatchRadiation = new JetPartonMatching(NfirstpTradiation, myselectedjetsradiation, matchingAlgo_, useMaxDist_, useDeltaR_, maxDist_);

									//check if all partons are matched ---> if so, just store jet ranks, nothing else (no jets left)
									if(GenMatchRadiation->getNumberOfUnmatchedPartons(0)== 0){
										if(verbose_) cout << "all partons are matched"<< endl;		
										MatchedTopRadiation++;				
						
										//store the rank of the jet for later use (if parton is unmatched ---> rank is negative)
										if(verbose_) cout << "store the ranks of the jets matched to the partons"<< endl;		
 										for(unsigned int i = 0; i < NfirstpTradiation.size(); i++){			
											int IdxRadJets = GenMatchRadiation->getMatchForParton(i,0);
											matchRadiationjetidx.push_back(IdxRadJets);
											if(IdxRadJets>=0) originalrankjetsradiation.push_back(originalrankjetsradiationinput[IdxRadJets]);
											else originalrankjetsradiation.push_back(IdxRadJets);
										}
							
									//check if unmatched partons ---> if so, categorize, store jet ranks, nothing else (no jets left)
									}else if(GenMatchRadiation->getNumberOfUnmatchedPartons(0)!= 0){ 
										if(verbose_) cout << "not all partons matched"<< endl;		
										UnmatchedTopRadiation++;
					
										//find parton which is not matched -> check eta/pT/DR ---> count different possibilities
										for(unsigned int i = 0; i<NfirstpTradiation.size(); i++){
											if(verbose_) cout << "find the partons which are not matched and categorize"<< endl;		
											int Idx = GenMatchRadiation->getMatchForParton(i,0);
											if(Idx<0){
												smallestDRRadi = 1000.;
												for(unsigned int j = 0; j<NfirstpTradiation.size(); j++){
													DRRadi[j] = 1000.;
													if(j!=i){
														DRRadi[j] = deltaR(NfirstpTradiation[i]->p4(),NfirstpTradiation[j]->p4());
	  												if(DRRadi[j] < smallestDRRadi) smallestDRRadi = DRRadi[j];
													}
												}
												if(fabs(NfirstpTradiation[i]->eta())>jetEta_){
													UnmatchedRadCategoryEta++;
												}else if(NfirstpTradiation[i]->pt()<jetPt_){
													UnmatchedRadCategoryPt++;
												}else if(smallestDRRadi>minDRunmatchedParton_){
													UnmatchedRadCategoryDR++;
												}else{
													UnmatchedRadNotCategorized++;
												}
											}
										}
							
 										//store the rank of the jet for later use (if parton is unmatched ---> rank is negative)
 										if(verbose_) cout << "store the ranks of the jets matched top the partons"<< endl;		
										for(unsigned int i = 0; i < NfirstpTradiation.size(); i++){			
											int IdxRadJets = GenMatchRadiation->getMatchForParton(i,0);
											matchRadiationjetidx.push_back(IdxRadJets);
											if(IdxRadJets>=0)originalrankjetsradiation.push_back(originalrankjetsradiationinput[IdxRadJets]);
											else originalrankjetsradiation.push_back(IdxRadJets);
										}
									}

				  			///////////////////// #top radiation <= myselectedjetsradiation.size() /////////////////////
          			}else if( NumberRadiation <= myselectedjetsradiation.size()){
 									if(verbose_) cout << "less or equal number of partons compared to the number of jets"<< endl;		
									JetPartonMatching *GenMatchRadiation = new JetPartonMatching(TopRadiation, myselectedjetsradiation, matchingAlgo_, useMaxDist_, useDeltaR_, maxDist_);
									
									//check if unmatched top radiation  ---> if so, categorize, store jet ranks 
									if(GenMatchRadiation->getNumberOfUnmatchedPartons(0)!= 0){
 										if(verbose_) cout << "not all partons matched"<< endl;		
										UnmatchedTopRadiation++;

										//find parton which is not matched -> check eta/pT/DR ---> count different possibilities
										for(unsigned int i = 0; i<NumberRadiation; i++){
 											if(verbose_) cout << "check which partons are not matched and categorize"<< endl;		
											int Idx = GenMatchRadiation->getMatchForParton(i,0);
											if(Idx<0){
												smallestDRRadi = 1000.;
												for(unsigned int j = 0; j<NumberRadiation; j++){
													DRRadi[j] = 1000.;
													if(j!=i){
														DRRadi[j] = deltaR(TopRadiation[i]->p4(),TopRadiation[j]->p4());
	  												if(DRRadi[j] < smallestDRRadi) smallestDRRadi = DRRadi[j];
													}
												}
												if(fabs(TopRadiation[i]->eta())>jetEta_){
													UnmatchedRadCategoryEta++;
												}else if(TopRadiation[i]->pt()<jetPt_){
													UnmatchedRadCategoryPt++;
												}else if(smallestDRRadi>minDRunmatchedParton_){
													UnmatchedRadCategoryDR++;
												}else{
													UnmatchedRadNotCategorized++;
												}
											}
										}
 							
										//store the rank of the jet for later use (if parton is unmatched ---> rank is negative)
 										if(verbose_) cout << "store the rank of the jets matched with the partons"<< endl;		
 										for(unsigned int i = 0; i < NumberRadiation; i++){			
											int IdxRadJets = GenMatchRadiation->getMatchForParton(i,0);
											matchRadiationjetidx.push_back(IdxRadJets);
											if(IdxRadJets>=0)originalrankjetsradiation.push_back(originalrankjetsradiationinput[IdxRadJets]);
											else originalrankjetsradiation.push_back(IdxRadJets);
										}
							
									//check if all top radiation partons are matched ---> if so, store jet ranks 
									}else if(GenMatchRadiation->getNumberOfUnmatchedPartons(0)== 0){
 										if(verbose_) cout << "all partons matched"<< endl;		
										MatchedTopRadiation++;
										
										for(unsigned int i = 0; i < NumberRadiation; i++){			
											int IdxRadJets = GenMatchRadiation->getMatchForParton(i,0);
											matchRadiationjetidx.push_back(IdxRadJets);
											if(IdxRadJets>=0)originalrankjetsradiation.push_back(originalrankjetsradiationinput[IdxRadJets]);
											else originalrankjetsradiation.push_back(IdxRadJets);
										}
									}
								}
								
							}else{
 								if(verbose_) cout << "not radiation from the top quarks"<< endl;		
								NoTopRadiation++;
							}
			
						//check if all ISR partons are matched ---> if so, store jet ranks, look at radiation (still jets left)
						}else if(GenMatchISR->getNumberOfUnmatchedPartons(0)== 0){ 
 							if(verbose_) cout << "all ISR partons are matched"<< endl;									
							MatchedISR++;

							for(unsigned int i = 0; i < NumberISR; i++){			
								int IdxISRjets = GenMatchISR->getMatchForParton(i,0);
								matchISRjetidx.push_back(IdxISRjets);
								if(IdxISRjets>=0) originalrankjetsISR.push_back(originalrankjetsISRinput[IdxISRjets]);
								else originalrankjetsISR.push_back(IdxISRjets);
							}
				
							///////////////////////////////////////////////////////////////////////////
							////////// LOOK AT TOPRADIATION  - JETS LEFT - ISR ALL MATCHED - //////////
							///////////////////////////////////////////////////////////////////////////
 							if(verbose_) cout << "check if there are partons coming from radiation of top quarks"<< endl;									
							if(NumberRadiation !=0){ 
 								if(verbose_) cout << "there are partons coming from the radiation of top quarks"<< endl;									
								ISRTopRadiation++;
								//make new list of jets to be used as input for matching (remove ISR-matched jets from jetlist used as input for ISR)
								double arrayRad[100];
								for(unsigned int i = 0; i<myselectedjets.size() ; i++){
									for(unsigned int l =0; l<quarks.size(); l++){
										for(unsigned int k =0; k<NumberISR; k++){
											if( (unsigned int) originalrankjetsISR[k]!=i && (unsigned int) matchquarkjetidx[l]!=i){
												arrayRad[i]=0;
											}else{
												arrayRad[i]=1;
											}
										}
									}
								}
								for(unsigned int a=0; a<myselectedjets.size();a++){
									if(arrayRad[a]==0){
										myselectedjetsradiation.push_back(myselectedjets[a]);
										originalrankjetsradiationinput.push_back(a);
									}
								}
								
 								if(verbose_) cout << "check the number of partons compared to the number of jets"<< endl;									
				  			///////////////////// #top radiation > myselectedjetsradiation.size() /////////////////////
								if(NumberRadiation > myselectedjetsradiation.size()){
 									if(verbose_) cout << "more partons compared to the number of jets"<< endl;									
									//match highest pT top radiation parton(s) with jets
									for(unsigned int i=0; i<TopRadiation.size();i++){TopRadiationpT.push_back(TopRadiation[i]->pt());}
									// 1. sort radiationpt on descending pt
									std::sort(TopRadiationpT.begin(),TopRadiationpT.end(),greater<double>());
									// 2. fill vector with the n first pt partons
									for(unsigned int n=0;n<myselectedjetsradiation.size();n++){
										for(unsigned int k=0;k<NumberRadiation;k++){
											if(TopRadiationpT[n]==TopRadiation[k]->pt())NfirstpTradiation.push_back(TopRadiation[k]);
										}
									}
									// 3. do JetPartonMatching using only the highest pT top radiation partons						
 									if(verbose_) cout << "do matching"<< endl;									
									JetPartonMatching *GenMatchRadiation = new JetPartonMatching(NfirstpTradiation, myselectedjetsradiation, matchingAlgo_, useMaxDist_, useDeltaR_, maxDist_);

									//check if all partons are matched ---> if so, just store jet ranks, nothing else (no jets left)
									if(GenMatchRadiation->getNumberOfUnmatchedPartons(0)== 0){
										MatchedTopRadiation++;				
							
 										if(verbose_) cout << "store the rank of the jets matched to the partons"<< endl;									
										//store the rank of the jet for later use (if parton is unmatched ---> rank is negative)
 										for(unsigned int i = 0; i < NfirstpTradiation.size(); i++){			
											int IdxRadJets = GenMatchRadiation->getMatchForParton(i,0);
											matchRadiationjetidx.push_back(IdxRadJets);
											if(IdxRadJets>=0)originalrankjetsradiation.push_back(originalrankjetsradiationinput[IdxRadJets]);
											else originalrankjetsradiation.push_back(IdxRadJets);
										}
							
 									if(verbose_) cout << "check if all partons are matched to the jets"<< endl;									
									//check if unmatched partons ---> if so, categorize, store jet ranks, nothing else (no jets left)
									}else if(GenMatchRadiation->getNumberOfUnmatchedPartons(0)!= 0){ 
 										if(verbose_) cout << "not all partons are matched to the quarks"<< endl;									
										UnmatchedTopRadiation++;
						
										//find parton which is not matched -> check eta/pT/DR ---> count different possibilities
										for(unsigned int i = 0; i<NfirstpTradiation.size(); i++){
 											if(verbose_) cout << "check which partons are not matched and categorize"<< endl;									
											int Idx = GenMatchRadiation->getMatchForParton(i,0);
											if(Idx<0){
												smallestDRRadi = 1000.;
												for(unsigned int j = 0; j<NfirstpTradiation.size(); j++){
													DRRadi[j] = 1000.;
													if(j!=i){
														DRRadi[j] = deltaR(NfirstpTradiation[i]->p4(),NfirstpTradiation[j]->p4());
	  												if(DRRadi[j] < smallestDRRadi) smallestDRRadi = DRRadi[j];
													}
												}
												if(fabs(NfirstpTradiation[i]->eta())>jetEta_){
													UnmatchedRadCategoryEta++;
												}else if(NfirstpTradiation[i]->pt()<jetPt_){
													UnmatchedRadCategoryPt++;
												}else if(smallestDRRadi>minDRunmatchedParton_){
													UnmatchedRadCategoryDR++;
												}else{
													UnmatchedRadNotCategorized++;
												}
											}
										}
							
 										if(verbose_) cout << "store the rank of the jets matched with the partons"<< endl;									
 										//store the rank of the jet for later use (if parton is unmatched ---> rank is negative)
 										for(unsigned int i = 0; i < NfirstpTradiation.size(); i++){			
											int IdxRadJets = GenMatchRadiation->getMatchForParton(i,0);
											matchRadiationjetidx.push_back(IdxRadJets);
											if(IdxRadJets>=0)originalrankjetsradiation.push_back(originalrankjetsradiationinput[IdxRadJets]);
											else originalrankjetsradiation.push_back(IdxRadJets);
										}
									}

				  			///////////////////// #top radiation <= myselectedjetsradiation.size() /////////////////////
          			}else if( NumberRadiation <= myselectedjetsradiation.size()){
 									if(verbose_) cout << "less or equal number of partons compared to the number of jets"<< endl;									
									JetPartonMatching *GenMatchRadiation = new JetPartonMatching(TopRadiation, myselectedjetsradiation, matchingAlgo_, useMaxDist_, useDeltaR_, maxDist_);
									
									//check if unmatched top radiation  ---> if so, categorize, store jet ranks 
									if(GenMatchRadiation->getNumberOfUnmatchedPartons(0)!= 0){
 										if(verbose_) cout << "not all partons are matched"<< endl;									
										UnmatchedTopRadiation++;

										//find parton which is not matched -> check eta/pT/DR ---> count different possibilities
										for(unsigned int i = 0; i<NumberRadiation; i++){
 											if(verbose_) cout << "find partons which are not matched and categorize"<< endl;									
											int Idx = GenMatchRadiation->getMatchForParton(i,0);
											if(Idx<0){
												smallestDRRadi = 1000.;
												for(unsigned int j = 0; j<NumberRadiation; j++){
													DRRadi[j] = 1000.;
													if(j!=i){
														DRRadi[j] = deltaR(TopRadiation[i]->p4(),TopRadiation[j]->p4());
	  												if(DRRadi[j] < smallestDRRadi) smallestDRRadi = DRRadi[j];
													}
												}
												if(fabs(TopRadiation[i]->eta())>jetEta_){
													UnmatchedRadCategoryEta++;
												}else if(TopRadiation[i]->pt()<jetPt_){
													UnmatchedRadCategoryPt++;
												}else if(smallestDRRadi>minDRunmatchedParton_){
													UnmatchedRadCategoryDR++;
												}else{
													UnmatchedRadNotCategorized++;
												}
											}
										}
 							
 										if(verbose_) cout << "store the ranks of the jets matched to the partons"<< endl;									
										//store the rank of the jet for later use (if parton is unmatched ---> rank is negative)
 										for(unsigned int i = 0; i < NumberRadiation; i++){			
											int IdxRadJets = GenMatchRadiation->getMatchForParton(i,0);
											matchRadiationjetidx.push_back(IdxRadJets);
											if(IdxRadJets>=0)originalrankjetsradiation.push_back(originalrankjetsradiationinput[IdxRadJets]);
											else originalrankjetsradiation.push_back(IdxRadJets);
										}
							
									//check if all top radiation partons are matched ---> if so, store jet ranks 
									}else if(GenMatchRadiation->getNumberOfUnmatchedPartons(0)== 0){
 										if(verbose_) cout << "all partons are matched"<< endl;									
										MatchedTopRadiation++;
										
										for(unsigned int i = 0; i < NumberRadiation; i++){			
											int IdxRadJets = GenMatchRadiation->getMatchForParton(i,0);
											matchRadiationjetidx.push_back(IdxRadJets);
											if(IdxRadJets>=0)originalrankjetsradiation.push_back(originalrankjetsradiationinput[IdxRadJets]);
											else originalrankjetsradiation.push_back(IdxRadJets);
										}
									}
								}
								
							}else{
 								if(verbose_) cout << "no partons from top radiation"<< endl;									
								NoTopRadiation++;
							}
						}
					}					
				
				///////////////////////////////////////////////////
				////////// LOOK AT TOPRADIATION - NO ISR //////////
				///////////////////////////////////////////////////	
				//no ISR partons (either not a MadGraph sample or really no ISR happened)	
				}else{ 
 					if(verbose_) cout << "no partons from ISR"<< endl;									
					NoISR++;
					
					if(NumberRadiation !=0){ 
 						if(verbose_) cout << "there are partons coming from radiation of top quarks"<< endl;									
						TopRadiationNoISR++;
						//make new list of jets to be used as input for matching (remove jets matched with quarks from jetlist (input for	radiation)
						double array[100];
						for(unsigned int i=0; i<myselectedjets.size(); i++ ){
							for(unsigned int k =0; k<quarks.size(); k++){
								if((unsigned int) matchquarkjetidx[k] !=i) {
									array[i]=0;
								}else{
									array[i]=1;
								}
							}
						}
						for(unsigned int a=0; a<myselectedjets.size();a++){
							if(array[a]==0){
								myselectedjetsradiation.push_back(myselectedjets[a]);
								originalrankjetsradiationinput.push_back(a);
							}
						}
						
								
 						if(verbose_) cout << "check the number of partons compared to the number of jets"<< endl;									
						///////////////////// #top radiation > myselectedjetsradiation.size() /////////////////////
						if(NumberRadiation > myselectedjetsradiation.size()){
 							if(verbose_) cout << "more partons than jets"<< endl;									
							//match highest pT top radiation parton(s) with jets
							for(unsigned int i=0; i<TopRadiation.size();i++){TopRadiationpT.push_back(TopRadiation[i]->pt());}
							// 1. sort radiationpt on descending pt
							std::sort(TopRadiationpT.begin(),TopRadiationpT.end(),greater<double>());
							// 2. fill vector with the n first pt partons
							for(unsigned int n=0;n<myselectedjetsradiation.size();n++){
								for(unsigned int k=0;k<NumberRadiation;k++){
									if(TopRadiationpT[n]==TopRadiation[k]->pt())NfirstpTradiation.push_back(TopRadiation[k]);
								}
							}
 							if(verbose_) cout << "do matching"<< endl;									
							// 3. do JetPartonMatching using only the highest pT top radiation partons						
							JetPartonMatching *GenMatchRadiation = new JetPartonMatching(NfirstpTradiation, myselectedjetsradiation, matchingAlgo_, useMaxDist_, useDeltaR_, maxDist_);

							//check if all partons are matched ---> if so, just store jet ranks, nothing else (no jets left)
							if(GenMatchRadiation->getNumberOfUnmatchedPartons(0)== 0){
 								if(verbose_) cout << "all partons are matched"<< endl;									
								MatchedTopRadiation++;				
						
 								if(verbose_) cout << "store the ranks of the jets matched with the partons"<< endl;									
								//store the rank of the jet for later use (if parton is unmatched ---> rank is negative)
 								for(unsigned int i = 0; i < NfirstpTradiation.size(); i++){			
									int IdxRadJets = GenMatchRadiation->getMatchForParton(i,0);
									matchRadiationjetidx.push_back(IdxRadJets);
									if(IdxRadJets>=0)originalrankjetsradiation.push_back(originalrankjetsradiationinput[IdxRadJets]);
									else originalrankjetsradiation.push_back(IdxRadJets);
								}
							
							//check if unmatched partons ---> if so, categorize, store jet ranks, nothing else (no jets left)
							}else if(GenMatchRadiation->getNumberOfUnmatchedPartons(0)!= 0){ 
 								if(verbose_) cout << "not all partons are matched"<< endl;									
								UnmatchedTopRadiation++;
						
								//find parton which is not matched -> check eta/pT/DR ---> count different possibilities
								for(unsigned int i = 0; i<NfirstpTradiation.size(); i++){
 									if(verbose_) cout << "find partons which are not matched and categorize"<< endl;									
									int Idx = GenMatchRadiation->getMatchForParton(i,0);
									if(Idx<0){
										smallestDRRadi = 1000.;
										for(unsigned int j = 0; j<NfirstpTradiation.size(); j++){
											DRRadi[j] = 1000.;
											if(j!=i){
												DRRadi[j] = deltaR(NfirstpTradiation[i]->p4(),NfirstpTradiation[j]->p4());
	  										if(DRRadi[j] < smallestDRRadi) smallestDRRadi = DRRadi[j];
											}
										}
										if(fabs(NfirstpTradiation[i]->eta())>jetEta_){
											UnmatchedRadCategoryEta++;
										}else if(NfirstpTradiation[i]->pt()<jetPt_){
											UnmatchedRadCategoryPt++;
										}else if(smallestDRRadi>minDRunmatchedParton_){
											UnmatchedRadCategoryDR++;
										}else{
											UnmatchedRadNotCategorized++;
										}
									}	
								}	
							
 								if(verbose_) cout << "store the ranks of the jets which are matched with the partons"<< endl;									
 								//store the rank of the jet for later use (if parton is unmatched ---> rank is negative)
 								for(unsigned int i = 0; i < NfirstpTradiation.size(); i++){			
									int IdxRadJets = GenMatchRadiation->getMatchForParton(i,0);
									matchRadiationjetidx.push_back(IdxRadJets);
									if(IdxRadJets>=0)originalrankjetsradiation.push_back(originalrankjetsradiationinput[IdxRadJets]);
									else originalrankjetsradiation.push_back(IdxRadJets);
								}
							}	

				  	///////////////////// #top radiation <= myselectedjetsradiation.size() /////////////////////
          	}else if( NumberRadiation <= myselectedjetsradiation.size()){
 							if(verbose_) cout << "less or equal number of partons compared to the number of jets"<< endl;									
							JetPartonMatching *GenMatchRadiation = new JetPartonMatching(TopRadiation, myselectedjetsradiation, matchingAlgo_, useMaxDist_, useDeltaR_, maxDist_);
							
							//check if unmatched top radiation  ---> if so, categorize, store jet ranks 
							if(GenMatchRadiation->getNumberOfUnmatchedPartons(0)!= 0){
 								if(verbose_) cout << "there are partons which are not matched"<< endl;									
								UnmatchedTopRadiation++;

								//find parton which is not matched -> check eta/pT/DR ---> count different possibilities
								for(unsigned int i = 0; i<NumberRadiation; i++){
 									if(verbose_) cout << "check which partons not matched and categorize"<< endl;									
									int Idx = GenMatchRadiation->getMatchForParton(i,0);
									if(Idx<0){
										smallestDRRadi = 1000.;
										for(unsigned int j = 0; j<NumberRadiation; j++){
											DRRadi[j] = 1000.;
											if(j!=i){
												DRRadi[j] = deltaR(TopRadiation[i]->p4(),TopRadiation[j]->p4());
	  										if(DRRadi[j] < smallestDRRadi) smallestDRRadi = DRRadi[j];
											}
										}
										if(fabs(TopRadiation[i]->eta())>jetEta_){
											UnmatchedRadCategoryEta++;
										}else if(TopRadiation[i]->pt()<jetPt_){
											UnmatchedRadCategoryPt++;
										}else if(smallestDRRadi>minDRunmatchedParton_){
											UnmatchedRadCategoryDR++;
										}else{
											UnmatchedRadNotCategorized++;
										}
									}
								}
 							
 								if(verbose_) cout << "store the ranks of the jets matched with the partons"<< endl;									
								//store the rank of the jet for later use (if parton is unmatched ---> rank is negative)
 								for(unsigned int i = 0; i < NumberRadiation; i++){			
									int IdxRadJets = GenMatchRadiation->getMatchForParton(i,0);
									matchRadiationjetidx.push_back(IdxRadJets);
									if(IdxRadJets>=0)originalrankjetsradiation.push_back(originalrankjetsradiationinput[IdxRadJets]);
									else originalrankjetsradiation.push_back(IdxRadJets);
								}
						
							//check if all top radiation partons are matched ---> if so, store jet ranks
							}else if(GenMatchRadiation->getNumberOfUnmatchedPartons(0)== 0){
 								if(verbose_) cout << "all partons are matched"<< endl;									
								MatchedTopRadiation++;
								
								for(unsigned int i = 0; i < NumberRadiation; i++){			
									int IdxRadJets = GenMatchRadiation->getMatchForParton(i,0);
									matchRadiationjetidx.push_back(IdxRadJets);
									if(IdxRadJets>=0)originalrankjetsradiation.push_back(originalrankjetsradiationinput[IdxRadJets]);
									else originalrankjetsradiation.push_back(IdxRadJets);
								}
							}
						}
					}else{
 						if(verbose_) cout << "no radiation from the top quarks"<< endl;									
						NoRadiation++;
						//find parton which is not matched -> check eta/pT/DR ---> count different possibilities
						bool only1 = true;
						for(unsigned int i = 0; i<quarks.size(); i++){
							if(verbose_) cout << "check which quarks not matched and categorize unmatched quarks"<< endl;		
							int Idx = GenMatchQuarks->getMatchForParton(i,0);
							if(only1 && Idx<0){
								only1 = false;
								smallestDRquarki = 1000.;
								for(unsigned int j = 0; j<quarks.size(); j++){
									DRquarki[j] = 1000.;
									if(j!=i){
										DRquarki[j] = deltaR(quarks[i]->p4(),quarks[j]->p4());
	  								if(DRquarki[j] < smallestDRquarki) smallestDRquarki = DRquarki[j];
									}
								}
								if(fabs(quarks[i]->eta())>jetEta_){
									UnmatchedQuark1NoRadCategoryEta++;
								}else if(quarks[i]->pt()<jetPt_){
									UnmatchedQuark1NoRadCategoryPt++;
								}else if(smallestDRquarki>minDRunmatchedParton_){
									UnmatchedQuark1NoRadCategoryDR++;
								}else{
									UnmatchedQuark1NoRadNotCategorized++;
								}
							}
						}
					}
				}

				////////////////////////////////////////////////////////////////
				////////// MAKE PLOTS WHEN ALL QUARKS ARE NOT MATCHED //////////
				////////////////////////////////////////////////////////////////	
			
				bool Radiation = false;
 				if(verbose_) cout << "store the 4 highest pT jets"<< endl;									
				for(unsigned int j=0; j<myselectedjets.size();j++){ 
					//make vector with 4 highest pT jets+rank
					if(j<4) {
						FourHighestPtJets.push_back(myselectedjets[j]);
					}
				}
				
				//keep jets + rank which are matched with radiation!!!
				double arraynew[100];
				if(NumberISR !=0 && NumberRadiation !=0){
	 				if(verbose_) cout << "keep the jets matched to radiation (ISR and top) and their ranks"<< endl;									
					for(unsigned int j=0; j<myselectedjets.size();j++){ 
						for(unsigned int l =0; l<originalrankjetsradiation.size(); l++){
							if(j == (unsigned int) originalrankjetsradiation[l]){
								Radiation = true;
								arraynew[j]=0;
								l=originalrankjetsradiation.size();		
							}else{arraynew[j]=1;}
						}
						if(arraynew[j] != 0){
							for(unsigned int k =0; k<originalrankjetsISR.size(); k++){
								if(j == (unsigned int) originalrankjetsISR[k]){ 
									Radiation = true;
									arraynew[j]=0;
									k=originalrankjetsISR.size();					
								}else{arraynew[j]=1;}
							}
						}
					}
				}else if((NumberISR ==0) && (NumberRadiation !=0)){
	 				if(verbose_) cout << "keep the jets matched to radiation (top) and their ranks"<< endl;									
					for(unsigned int j=0; j<myselectedjets.size();j++){ 
						for(unsigned int l =0; l<originalrankjetsradiation.size(); l++){
							if(j == (unsigned int) originalrankjetsradiation[l]){
								Radiation = true;
								arraynew[j]=0;
								l=originalrankjetsradiation.size();
							}else{arraynew[j]=1;}
						}
					}
				}else if((NumberISR !=0) && (NumberRadiation ==0)){
	 				if(verbose_) cout << "keep the jets matched to radiation (ISR) and their ranks"<< endl;									
					for(unsigned int j=0; j<myselectedjets.size();j++){ 
						for(unsigned int k =0; k<originalrankjetsISR.size(); k++){
							if(j == (unsigned int) originalrankjetsISR[k] ){
								Radiation = true;
								arraynew[j]=0;
								k=originalrankjetsISR.size();
							}else{arraynew[j]=1;}
						}
					}
				}
				
				bool jetmatchedwithradiation = false;
				if(Radiation){
 					if(verbose_) cout << "store the jets and their pt"<< endl;									
					for(unsigned int j=0; j<myselectedjets.size();j++){ 
						if(arraynew[j]==0){
							jetmatchedwithradiation = true;
							pair<int,pat::Jet> pairjetsrad(j,myselectedjets[j]);
							patJetsMatchedWithRadiation.push_back(pairjetsrad);
						}
					}	
				}
			
				if(Radiation && jetmatchedwithradiation){
 					int countradmultiplicity = 0;
					countradmultiplicity = patJetsMatchedWithRadiation.size();
					JetRadiationMultiplicity -> Fill(countradmultiplicity);
					
					double arraynewnew[4];
					bool notgood = false;
					for(unsigned int i = 0; i<FourHighestPtJets.size();i++){		
						for(unsigned int j = 0; j<patJetsMatchedWithRadiation.size();j++){
							if(i == (unsigned int) patJetsMatchedWithRadiation[j].first){
								notgood = true;
								arraynewnew[i]=0;
								//if among the 4 highest pT jets, some are matched with radiation:
								//store the jet rank
								//store the jet pT
								rankjetsmatchedtorad.push_back(patJetsMatchedWithRadiation[j].first);
								ptjetsmatchedtorad.push_back((patJetsMatchedWithRadiation[j].second).pt());
								j = patJetsMatchedWithRadiation.size();
							}else{arraynewnew[i]=1;}
						}
					}
					
					int count = 0;
					if(notgood){
						//highest pT jet has per definition lowest rank!
						//sort first on rank --> make plot of lowest
						//sort on pT -->make plot of highest
						std::sort(rankjetsmatchedtorad.begin(), rankjetsmatchedtorad.end(), std::greater<int>() );
   					std::sort(ptjetsmatchedtorad.begin(), ptjetsmatchedtorad.end(), std::greater<double>() );
						HighestPtRadMatchedJet_Rank->Fill(rankjetsmatchedtorad.back()+1);
   					HighestPtRadMatchedJet_Pt  ->Fill(ptjetsmatchedtorad.front());
 						if(verbose_) cout << "some of the jets are matched with radiation, count how many"<< endl;									
						for(unsigned int i = 0; i<FourHighestPtJets.size(); i++){
							if(arraynewnew[i]==0){
								count++;
							}
						}
					}
					
					if(count == 0)Highest4NoRad++;
					if(count >= 1)Highest4AtLeast1Rad++;
					if(count >= 2)Highest4AtLeast2Rad++;
					if(count >= 3)Highest4AtLeast3Rad++;
					if(count >= 4)Highest4AtLeast4Rad++;

				}else{
					edm::LogInfo("NoDataFound") << "no radiation, so no related plots are filled \n";
				}
			}
		}else{
			genObjectsNULL++;
		}

	}else if(useMatchingFromPAT_ && allok){
		edm::LogError ("NoDataFound") << "useMatchingFromPAT_ = true, but unfortunately this part of the code is not implemented \n";
	}
}


// ------------ method called once each job just before starting event loop  ------------
void 
TruthReco::beginJob(const edm::EventSetup&)
{
	edm::Service<TFileService> fs;
  if (!fs) throw edm::Exception(edm::errors::Configuration, "TFileService missing from configuration!");
	
	NbSemiMuRecoEvents = 0;
  NbSemiMuGenEvents = 0;
	genObjectsNotNULL = 0;
	recoEventsPassed = 0;
	genObjectsNULL = 0;
	EventNotMatched = 0;
	UnmatchedQuarks = 0;
	MatchedQuarks = 0;
	MatchedISR = 0;
	UnmatchedISR = 0;
	UnmatchedTopRadiation = 0;
	MatchedTopRadiation = 0;
	NoTopRadiation = 0;
	NoISR = 0;
	NoRadiation = 0;
	ISRTopRadiation = 0;
	TopRadiationNoISR = 0;
	NrNotSemiMu = 0;
	NrAnalyzed = 0;
	NrNot4Jets = 0;
	NrNot1Mu = 0;
	NrNot1MET = 0;
	NrNot4SelJets = 0;
	NrNot1SelMu = 0;
	NrNot1SelIsoMu = 0;
	Highest4NoRad = 0;	
	Highest4AtLeast1Rad = 0;
	Highest4AtLeast2Rad = 0;
	Highest4AtLeast3Rad = 0;
	Highest4AtLeast4Rad = 0;
	UnmatchedQuarkCategoryEta = 0;
	UnmatchedQuarkCategoryPt = 0;
	UnmatchedQuarkCategoryDR = 0;
	UnmatchedQuarkNotCategorized = 0;	
	UnmatchedISRCategoryEta = 0;
	UnmatchedISRCategoryPt = 0;
	UnmatchedISRCategoryDR = 0;
	UnmatchedISRNotCategorized = 0;
	UnmatchedRadCategoryEta = 0;
	UnmatchedRadCategoryPt = 0;
	UnmatchedRadCategoryDR = 0;
	UnmatchedRadNotCategorized = 0;
	NumberQuarks = 0;
	NumberISR = 0; 
	NumberRadiationLep = 0;
	NumberRadiationHad = 0;
	NumberRadiation = 0;	
	MatchedQuarksRadiation = 0;
	MatchedQuarksNoRadiation = 0;
	MatchedQuarksNoISR = 0;
	MatchedQuarksNoTopRadiation = 0;
	UnmatchedQuark1NoRadCategoryEta = 0;
	UnmatchedQuark1NoRadCategoryPt = 0;
	UnmatchedQuark1NoRadCategoryDR = 0;
	UnmatchedQuark1NotCategorized = 0;
 	UnmatchedQuark1CategoryEta = 0;
	UnmatchedQuark1CategoryPt = 0;
	UnmatchedQuark1CategoryDR = 0;
	UnmatchedQuark1NotCategorized = 0;
	Only1UnmatchedQuarks = 0;
	Only1UnmatchedQuarkCategoryEta = 0;
	Only1UnmatchedQuarkCategoryPt = 0;
	Only1UnmatchedQuarkCategoryDR = 0;
	Only1UnmatchedQuarkNotCategorized = 0;			
	QuarksMatchedUnmatchedISRtotal = 0;
	QuarksMatchedUnmatchedISRCategoryEta = 0;
	QuarksMatchedUnmatchedISRCategoryPt = 0;
	QuarksMatchedUnmatchedISRCategoryDR = 0;
	QuarksMatchedUnmatchedISRNotCategorized = 0;
	QuarksMatchedUnmatchedISR = 0;
	QuarksMatchedMatchedISR = 0;
	QuarksMatchedISRTopRadiation = 0;
	QuarksMatchedMatchedTopRadiation = 0;
	QuarksMatchedUnmatchedTopRadiation = 0;
	QuarksMatchedUnmatchedRadCategoryEta = 0;
	QuarksMatchedUnmatchedRadCategoryPt = 0;
	QuarksMatchedUnmatchedRadCategoryDR = 0;
	QuarksMatchedUnmatchedRadNotCategorized = 0;
	QuarksMatchedNoTopRadiation = 0;
	QuarksMatchedNoISR = 0;
	QuarksMatchedTopRadiationNoISR = 0;		
	QuarksMatchedNoRadiation = 0;
	QuarksMatchedHighest4NoRad = 0;
	QuarksMatchedHighest4AtLeast1Rad = 0;
	QuarksMatchedHighest4AtLeast2Rad = 0;
	QuarksMatchedHighest4AtLeast3Rad = 0;
	QuarksMatchedHighest4AtLeast4Rad = 0;
		
	TFileDirectory subDir = fs->mkdir( "ExpectedJECinbins" );

	//histograms
	ptnrbins        = pTbinVals_.size()-1;
  double *ptbins  = new double[pTbinVals_.size()];
  for(unsigned int b=0; b<pTbinVals_.size(); b++)  ptbins[b]  = pTbinVals_[b];
	etanrbins  = etabinVals_.size()-1;
  double *etabins    = new double[etabinVals_.size()];
  for(unsigned int b=0; b<etabinVals_.size(); b++) etabins[b] = etabinVals_[b];
  
	
	for(Int_t etab=0; etab<etanrbins; etab++) {	
    for(Int_t ptb=0; ptb<ptnrbins; ptb++) {
      TString obsNamelight = "ExpectedLightJEC_etabin"; obsNamelight += etab; obsNamelight += "_ptbin"; obsNamelight += ptb;
      TString obsNameb = "ExpectedBJEC_etabin"; obsNameb += etab; obsNameb += "_ptbin"; obsNameb += ptb;
			hExpLightJECPtEtaBin[etab][ptb] = subDir.make<TH1F>(obsNamelight,obsNamelight,500,-5.,5.);
			hExpBJECPtEtaBin[etab][ptb] 		= subDir.make<TH1F>(obsNameb,obsNameb,500,-5.,5.);
      fExpLightJECPtEtaBin[etab][ptb] = subDir.make<TF1>("F_"+obsNamelight,"gaus");
      fExpBJECPtEtaBin[etab][ptb] 		= subDir.make<TF1>("F_"+obsNameb,"gaus");
    }
  }
  
	hExpLightJECIncl				= fs->make<TH1F>("hExpLightJECIncl","hExpLightJECIncl",500,-5.,5.);
  hExpBJECIncl						= fs->make<TH1F>("hExpBJECIncl","hExpBJECIncl",500,-5.,5.);
  fExpLightJECIncl				= fs->make<TF1>("F_ExpLightJECIncl","gaus");
  fExpBJECIncl						= fs->make<TF1>("F_ExpLightJECIncl","gaus");

	LowestPtMatchedJet_Rank = fs->make<TH1F>("LowestPtMatchedJet_Rank","LowestPtMatchedJet_Rank",20,0.,20.);
	LowestPtMatchedJet_Pt 	= fs->make<TH1F>("LowestPtMatchedJet_Pt","LowestPtMatchedJet_Pt",500,0.,500.);
	mWHad_Gen 							= fs->make<TH1F>("mWHad_Gen","mWHad_Gen",500,0.,500.);
	mWLep_Gen 							= fs->make<TH1F>("mWLep_Gen","mWLep_Gen",500,0.,500.);
	mtHad_Gen 							= fs->make<TH1F>("mtHad_Gen","mtHad_Gen",1000,0.,1000.);
	mtLep_Gen 							= fs->make<TH1F>("mtLep_Gen","mtLep_Gen",1000,0.,500.);	
	mWHad_Rec 							= fs->make<TH1F>("mWHad_Rec","mWHad_Rec",500,0.,500.);
	mWLep_Rec 							= fs->make<TH1F>("mWLep_Rec","mWLep_Rec",500,0.,500.);
	mtHad_Rec 							= fs->make<TH1F>("mtHad_Rec","mtHad_Rec",1000,0.,1000.);
	mtLep_Rec 							= fs->make<TH1F>("mtLep_Rec","mtLep_Rec",1000,0.,1000.);	
	minJetMuVSminQuarkMu 		= fs->make<TH2F>("minJetMuVSminQuarkMu","minJetMuVSminQuarkMu",500,0.,5.,500,0.,5.);			 	
	minJetsVSminQuarks 			= fs->make<TH2F>("minJetsVSminQuarks","minJetsVSminQuarks",500,0.,5.,500,0.,5.);			 	
	HighestPtRadMatchedJet_Rank = fs->make<TH1F>("HighestPtRadMatchedJet_Rank","HighestPtRadMatchedJet_Rank",20,0.,20.);
	HighestPtRadMatchedJet_Pt 	= fs->make<TH1F>("HighestPtRadMatchedJet_Pt","HighestPtRadMatchedJet_Pt",500,0.,500.);
	UnmatchedQuarkMultiplicity	=	fs->make<TH1F>("UnmatchedQuarkMultiplicity","UnmatchedQuarkMultiplicity",5,0.,5.);
	JetRadiationMultiplicity		=	fs->make<TH1F>("JetRadiationMultiplicity","JetRadiationMultiplicity",20,0.,20.);
	Dphi_HadTops						= fs->make<TH1F>("Dphi_HadTops","Dphi_HadTops",700,0.,3.5);
	Dphi_LepTops						= fs->make<TH1F>("Dphi_LepTops","Dphi_LepTops",700,0.,3.5);
	Deta_HadTops						= fs->make<TH1F>("Deta_HadTops","Deta_HadTops",1000,0.,5.);
	Deta_LepTops						= fs->make<TH1F>("Deta_LepTops","Deta_LepTops",1000,0.,5.);
	DR_HadTops							= fs->make<TH1F>("DR_HadTops","DR_HadTops",1000,0.,5.);
	DR_LepTops							= fs->make<TH1F>("DR_LepTops","DR_LepTops",1000,0.,5.);
	Dphi_Tops_Rec						= fs->make<TH1F>("Dphi_Tops_Rec","Dphi_Tops_Rec",700,0.,3.5);
	Dphi_Tops_Gen						= fs->make<TH1F>("Dphi_Tops_Gen","Dphi_Tops_Gen",700,0.,3.5);
	DR_Tops_Rec							= fs->make<TH1F>("DR_Tops_Rec","DR_Tops_Rec",1000,0.,5.);
	DR_Tops_Gen							= fs->make<TH1F>("DR_Tops_Gen","DR_Tops_Gen",1000,0.,5.);
	QuarksMatchedJetRadiationMultiplicity	=	fs->make<TH1F>("QuarksMatchedJetRadiationMultiplicity","QuarksMatchedJetRadiationMultiplicity",20,0.,20.);
	QuarksMatchedHighestPtRadMatchedJet_Rank = fs->make<TH1F>("QuarksMatchedHighestPtRadMatchedJet_Rank","QuarksMatchedHighestPtRadMatchedJet_Rank",20,0.,20.);
	QuarksMatchedHighestPtRadMatchedJet_Pt 	= fs->make<TH1F>("QuarksMatchedHighestPtRadMatchedJet_Pt","QuarksMatchedHighestPtRadMatchedJet_Pt",500,0.,500.);
	LowestPtQuarkMatchedJetPt_vs_HighestPtRadMatchedJetPt	= fs->make<TH2F>("LowestPtQuarkMatchedJetPt_vs_HighestPtRadMatchedJetPt","LowestPtQuarkMatchedJetPt_vs_HighestPtRadMatchedJetPt",500,0.,500.,500,0.,500.);
  
	delete [] etabins; 
  delete [] ptbins; 


}

// ------------ method called once each job just after ending the event loop  ------------
void 
TruthReco::endJob() {

	for(int etab=0; etab<etanrbins; etab++) {	
    for(int ptb=0; ptb<ptnrbins; ptb++) {
      double maxcontent = 0.;
			int maxbin = 0;
			for(int nb=1; nb<hExpLightJECPtEtaBin[etab][ptb]->GetNbinsX(); nb ++){
	  		if (hExpLightJECPtEtaBin[etab][ptb]->GetBinContent(nb)>maxcontent) {
	    		maxcontent = hExpLightJECPtEtaBin[etab][ptb]->GetBinContent(nb);
	    		maxbin = nb;
	 			}
			}
			int range = (int)(hExpLightJECPtEtaBin[etab][ptb]->GetNbinsX()/6); //in order that ~1/3 of X-axis range is fitted
  		fExpLightJECPtEtaBin[etab][ptb] -> SetRange(hExpLightJECPtEtaBin[etab][ptb]->GetBinCenter(maxbin-range),
																									hExpLightJECPtEtaBin[etab][ptb]->GetBinCenter(maxbin+range));
			fExpLightJECPtEtaBin[etab][ptb] -> SetParameters(hExpLightJECPtEtaBin[etab][ptb] -> GetMaximum(),
	                                             				hExpLightJECPtEtaBin[etab][ptb] -> GetMean(),
	 					     																			hExpLightJECPtEtaBin[etab][ptb] -> GetRMS());
			hExpLightJECPtEtaBin[etab][ptb] -> Fit(fExpLightJECPtEtaBin[etab][ptb]->GetName(),"RQ");
		}
  }

	for(int etab=0; etab<etanrbins; etab++) {	
    for(int ptb=0; ptb<ptnrbins; ptb++) {
      double maxcontent = 0.;
			int maxbin = 0;
			for(int nb=1; nb<hExpBJECPtEtaBin[etab][ptb]->GetNbinsX(); nb ++){
	  		if (hExpBJECPtEtaBin[etab][ptb]->GetBinContent(nb)>maxcontent) {
	    		maxcontent = hExpBJECPtEtaBin[etab][ptb]->GetBinContent(nb);
	    		maxbin = nb;
	  		}
			}
			int range2 = (int)(hExpBJECPtEtaBin[etab][ptb]->GetNbinsX()/6); //in order that ~1/3 of X-axis range is fitted
  		fExpBJECPtEtaBin[etab][ptb] -> SetRange(hExpBJECPtEtaBin[etab][ptb]->GetBinCenter(maxbin-range2),
																							hExpBJECPtEtaBin[etab][ptb]->GetBinCenter(maxbin+range2));
			fExpBJECPtEtaBin[etab][ptb] -> SetParameters(hExpBJECPtEtaBin[etab][ptb] -> GetMaximum(),
	                                            		hExpBJECPtEtaBin[etab][ptb] -> GetMean(),
	 				     																		hExpBJECPtEtaBin[etab][ptb] -> GetRMS());
			hExpBJECPtEtaBin[etab][ptb] -> Fit(fExpBJECPtEtaBin[etab][ptb]->GetName(),"RQ");
		}
  }

 	double maxcontent = 0.;
	int maxbin = 0;
	for(int nb=1; nb<hExpLightJECIncl->GetNbinsX(); nb ++){
		if (hExpLightJECIncl->GetBinContent(nb)>maxcontent) {
	 		maxcontent = hExpLightJECIncl->GetBinContent(nb);
	 		maxbin = nb;
		}
	}
	int rangeIncl = (int)(hExpLightJECIncl->GetNbinsX()/6); //in order that ~1/3 of X-axis range is fitted
  fExpLightJECIncl -> SetRange(hExpLightJECIncl->GetBinCenter(maxbin-rangeIncl),
																hExpLightJECIncl->GetBinCenter(maxbin+rangeIncl));
	fExpLightJECIncl -> SetParameters(hExpLightJECIncl -> GetMaximum(),
	                                  hExpLightJECIncl -> GetMean(),
	 					     										hExpLightJECIncl -> GetRMS());
	hExpLightJECIncl -> Fit(fExpLightJECIncl->GetName(),"RQ");

 	maxcontent = 0.;
	maxbin = 0;

	for(int nb=1; nb<hExpBJECIncl->GetNbinsX(); nb ++){
		if (hExpBJECIncl->GetBinContent(nb)>maxcontent) {
	 		maxcontent = hExpBJECIncl->GetBinContent(nb);
	 		maxbin = nb;
		}
	}
  int rangeIncl2 = (int)(hExpBJECIncl->GetNbinsX()/6); //in order that ~1/3 of X-axis range is fitted
	fExpBJECIncl -> SetRange(hExpBJECIncl->GetBinCenter(maxbin-rangeIncl2),
														hExpBJECIncl->GetBinCenter(maxbin+rangeIncl2));
	fExpBJECIncl -> SetParameters(hExpBJECIncl -> GetMaximum(),
	                              hExpBJECIncl -> GetMean(),
	 					     								hExpBJECIncl -> GetRMS());
	hExpBJECIncl -> Fit(fExpBJECIncl->GetName(),"RQ");

 double eff1 = ((double)NbSemiMuGenEvents/(double)NrAnalyzed) * 100;
 double eff2 = ((double)NbSemiMuRecoEvents/(double)NrAnalyzed) * 100;
 double eff3 = ((double)NrNot4SelJets/(double)NrAnalyzed) * 100;
 double eff4 = ((double)NrNot1SelMu/(double)NrAnalyzed) * 100;
 double eff5 = ((double)NrNot1SelIsoMu/(double)NrAnalyzed) * 100;
 double eff6 = ((double)NrNot1MET/(double)NrAnalyzed) * 100;
 double eff6b = ((double)EventNotMatched/(double)NrAnalyzed) * 100;
 double eff6c = ((double)recoEventsPassed/(double)NrAnalyzed) * 100;
 double eff7 = ((double)NrNotSemiMu/(double)NrAnalyzed) * 100;
 double eff8 = ((double)MatchedQuarks/(double)recoEventsPassed) * 100;
 double eff9 = ((double)MatchedQuarksRadiation/(double)MatchedQuarks) * 100;
 double eff10 = ((double)MatchedQuarksNoISR/(double)MatchedQuarks) * 100;
 double eff11 = ((double)MatchedQuarksNoTopRadiation/(double)MatchedQuarks) * 100;
 double eff12 = ((double)MatchedQuarksNoRadiation/(double)MatchedQuarks) * 100;
 double eff13 = ((double)QuarksMatchedMatchedISR/(double)MatchedQuarks) * 100;
 double eff14 = ((double)QuarksMatchedUnmatchedISR/(double)MatchedQuarks) * 100;
 double eff15 = ((double)QuarksMatchedMatchedTopRadiation/(double)MatchedQuarks) * 100;
 double eff16 = ((double)QuarksMatchedUnmatchedTopRadiation/(double)MatchedQuarks) * 100;
 double eff17 = ((double)QuarksMatchedNoTopRadiation/(double)MatchedQuarks) * 100;
 double eff18 = ((double)QuarksMatchedISRTopRadiation/(double)MatchedQuarks) * 100;
 double eff19 = ((double)QuarksMatchedTopRadiationNoISR/(double)MatchedQuarks) * 100;
 double eff20 = ((double)QuarksMatchedNoRadiation/(double)MatchedQuarks) * 100;
 double eff21 = ((double)QuarksMatchedHighest4NoRad/(double)MatchedQuarks) * 100;
 double eff22 = ((double)QuarksMatchedHighest4AtLeast1Rad/(double)MatchedQuarks) * 100;
 double eff23 = ((double)QuarksMatchedHighest4AtLeast2Rad/(double)MatchedQuarks) * 100;
 double eff24 = ((double)QuarksMatchedHighest4AtLeast3Rad/(double)MatchedQuarks) * 100;
 double eff25 = ((double)QuarksMatchedHighest4AtLeast4Rad/(double)MatchedQuarks) * 100;
 double eff26 = ((double)UnmatchedQuarks/(double)recoEventsPassed) * 100;
 double eff27 = ((double)UnmatchedQuark1CategoryEta/(double)UnmatchedQuarks) * 100;
 double eff28 = ((double)UnmatchedQuark1CategoryPt/(double)UnmatchedQuarks) * 100;
 double eff29 = ((double)UnmatchedQuark1CategoryDR/(double)UnmatchedQuarks) * 100;
 double eff30 = ((double)UnmatchedQuark1NotCategorized/(double)UnmatchedQuarks) * 100;
 double eff31 = ((double)Only1UnmatchedQuarks/(double)recoEventsPassed) * 100;
 double eff32 = ((double)Only1UnmatchedQuarkCategoryEta/(double)Only1UnmatchedQuarks) * 100;
 double eff33 = ((double)Only1UnmatchedQuarkCategoryPt/(double)Only1UnmatchedQuarks) * 100;
 double eff34 = ((double)Only1UnmatchedQuarkCategoryDR/(double)Only1UnmatchedQuarks) * 100;
 double eff35 = ((double)Only1UnmatchedQuarkNotCategorized/(double)Only1UnmatchedQuarks) * 100;
 double eff36 = ((double)MatchedISR/(double)UnmatchedQuarks) * 100;
 double eff37 = ((double)UnmatchedISR/(double)UnmatchedQuarks) * 100;
 double eff38 = ((double)MatchedTopRadiation/(double)UnmatchedQuarks) * 100;
 double eff39 = ((double)UnmatchedTopRadiation/(double)UnmatchedQuarks) * 100;
 double eff40 = ((double)NoTopRadiation/(double)UnmatchedQuarks) * 100;
 double eff41 = ((double)ISRTopRadiation/(double)UnmatchedQuarks) * 100;
 double eff42 = ((double)TopRadiationNoISR/(double)UnmatchedQuarks) * 100;
 double eff43 = ((double)NoRadiation/(double)UnmatchedQuarks) * 100;
 double eff44 = ((double)UnmatchedQuark1NoRadCategoryEta/(double)NoRadiation) * 100;
 double eff45 = ((double)UnmatchedQuark1NoRadCategoryPt/(double)NoRadiation) * 100;
 double eff46 = ((double)UnmatchedQuark1NoRadCategoryDR/(double)NoRadiation) * 100;
 double eff47 = ((double)UnmatchedQuark1NoRadNotCategorized/(double)NoRadiation) * 100;
 double eff48 = ((double)Highest4NoRad/(double)UnmatchedQuarks) * 100;
 double eff49 = ((double)Highest4AtLeast1Rad/(double)UnmatchedQuarks) * 100;
 double eff50 = ((double)Highest4AtLeast2Rad/(double)UnmatchedQuarks) * 100;
 double eff51 = ((double)Highest4AtLeast3Rad/(double)UnmatchedQuarks) * 100;
 double eff52 = ((double)Highest4AtLeast4Rad/(double)UnmatchedQuarks) * 100;

  edm::LogVerbatim ("MainResults") << " \n\n";
  edm::LogVerbatim ("MainResults") << " ------------------------------------";
  edm::LogVerbatim ("MainResults") << " ------------------------------------";
  edm::LogVerbatim ("MainResults") << " -- Report from TruthReco module   --";
  edm::LogVerbatim ("MainResults") << " ------------------------------------";
  edm::LogVerbatim ("MainResults") << " ------------------------------------";

	edm::LogVerbatim  ("MainResults") << "  (1) " <<NrAnalyzed << " number of analyzed events ";
	edm::LogVerbatim  ("MainResults") << "  (2) " <<NbSemiMuGenEvents << " number of semi-mu genevents before any cuts ";
	edm::LogVerbatim  ("MainResults") << "  (3) (2)/(1) " <<eff1 << "% semi-mu genevents ";
	edm::LogVerbatim  ("MainResults") << "  (4) " <<NbSemiMuRecoEvents << " number of semi-mu reco'ed events (>=4jets,>=1muon, >0MET) before any cuts ";
	edm::LogVerbatim  ("MainResults") << "  (5) (4)/(1) " <<eff2 << "% semi-mu reco'ed events ";
	edm::LogVerbatim  ("MainResults") << "  (6) (6)/(1) " <<eff3 << "% or "<<NrNot4SelJets <<" events with less then 4 jets passing cuts"; 
	edm::LogVerbatim  ("MainResults") << "  (7) (7)/(1) " <<eff4 << "% or "<<NrNot1SelMu <<" events with muon not passing kinematic cuts "; 
	edm::LogVerbatim  ("MainResults") << "  (8) (8)/(1) " <<eff5 << "% or "<<NrNot1SelIsoMu <<" events with muon not passing cuts "; 
	edm::LogVerbatim  ("MainResults")	<< "  (9) (9)/(1) " <<eff6 << "% or "<<NrNot1MET <<" events with less then 1 MET "; 
	edm::LogVerbatim  ("MainResults") << "  (10) (10)/(1) " <<eff6b << "% or "<<EventNotMatched << " events where the muon is not matched ";
	edm::LogVerbatim  ("MainResults") << "  (11) (11)/(1) " <<eff6c << "% or "<<recoEventsPassed << " events which pass event selection";
	edm::LogVerbatim  ("MainResults") << "  (12) " <<genObjectsNotNULL << " events for which all the genObjects are there \n";

  edm::LogVerbatim ("MainResults") << " ----------------------------";
  edm::LogVerbatim ("MainResults") << " -- Inclusive expected JEC --";
  edm::LogVerbatim ("MainResults") << " ----------------------------";

	edm::LogVerbatim  ("MainResults") << "  (13) " << hExpLightJECIncl -> GetMean()<< " mean light jets ";
	edm::LogVerbatim  ("MainResults") << "  (14) " << hExpLightJECIncl -> GetRMS()<< " RMS light jets \n";
	edm::LogVerbatim  ("MainResults") << "  (15) " << hExpBJECIncl -> GetMean() << " mean b jets ";
	edm::LogVerbatim  ("MainResults") << "  (16) " << hExpBJECIncl -> GetRMS()<< " RMS b jets \n";

  edm::LogVerbatim ("MainResults") << " --------------------------";
  edm::LogVerbatim ("MainResults") << " -- reconstructed masses --";
  edm::LogVerbatim ("MainResults") << " --------------------------";

	edm::LogVerbatim  ("MainResults") << "  (17) " << mWHad_Rec -> GetMean()<< " mean hadronic W mass ";
	edm::LogVerbatim  ("MainResults") << "  (18) " << mWHad_Rec -> GetRMS()<< " RMS hadronic W mass \n";
	edm::LogVerbatim  ("MainResults") << "  (19) " << mtHad_Rec -> GetMean() << " mean hadronic top masss ";
	edm::LogVerbatim  ("MainResults") << "  (20) " << mtHad_Rec -> GetRMS()<< " RMS hadronic top mass \n";
	edm::LogVerbatim  ("MainResults") << "  (21) " << mWLep_Rec -> GetMean()<< " mean leptonic W mass ";
	edm::LogVerbatim  ("MainResults") << "  (22) " << mWLep_Rec -> GetRMS()<< " RMS leptonic W mass \n";
	edm::LogVerbatim  ("MainResults") << "  (23) " << mtLep_Rec -> GetMean() << " mean leptonic top masss ";
	edm::LogVerbatim  ("MainResults") << "  (24) " << mtLep_Rec -> GetRMS()<< " RMS leptonic top mass \n";


  edm::LogVerbatim ("MainResults") << " -------------------------";
  edm::LogVerbatim ("MainResults") << " -- Results on matching --";
  edm::LogVerbatim ("MainResults") << " -------------------------";

	
	edm::LogVerbatim  ("MainResults") << "  (25) (25)/(11) "<< eff8<< "% or " <<MatchedQuarks << " events where all quarks are matched";	
	edm::LogVerbatim  ("MainResults") << "  (26) \t" <<"of which (26)/(25) " <<eff9<< "% or " <<MatchedQuarksRadiation << " events where there is radiation";	
	edm::LogVerbatim  ("MainResults") << "  (27) \t" <<"of which (27)/(25) " <<eff10<< "% or " <<MatchedQuarksNoISR << " events where there is only top radiation";	
	edm::LogVerbatim  ("MainResults") << "  (28) \t" <<"of which (28)/(25) " <<eff11<< "% or " <<MatchedQuarksNoTopRadiation << " events where there is only ISR";	
	edm::LogVerbatim  ("MainResults") << "  (29) \t" <<"of which (29)/(25) " <<eff12<< "% or " <<MatchedQuarksNoRadiation << " events where there is both ISR and	top radiation";	
	edm::LogVerbatim  ("MainResults") << "  (30) \t\t(30)/(25) "<< eff13<< "% or "  <<QuarksMatchedMatchedISR << " events where all ISR partons are matched";	
	edm::LogVerbatim  ("MainResults") << "  (31a) \t\t(31a)/(25) "<< eff14<< "% or "  <<QuarksMatchedUnmatchedISR << " events where not all ISR partons are matched (there can be more than 1	unmatched parton)";
	edm::LogVerbatim  ("MainResults") << "  (31b) \t\t "<< QuarksMatchedUnmatchedISRtotal << " total	unmatched ISR partons";
	edm::LogVerbatim  ("MainResults") << "  (32) \t\t\t" <<"of which " << QuarksMatchedUnmatchedISRCategoryEta << " because of the eta requirement";
	edm::LogVerbatim  ("MainResults") << "  (33) \t\t\t" <<"of which " << QuarksMatchedUnmatchedISRCategoryPt << " because of the p_{T} requirement";
	edm::LogVerbatim  ("MainResults") << "  (34) \t\t\t" <<"of which " << QuarksMatchedUnmatchedISRCategoryDR << " because of the DR requirement";
	edm::LogVerbatim  ("MainResults") << "  (35) \t\t\t" <<"of which " << QuarksMatchedUnmatchedISRNotCategorized << " not categorized\n";	
	edm::LogVerbatim  ("MainResults") << "  (36) \t\t(36)/(25) "<< eff15<< "% or " << QuarksMatchedMatchedTopRadiation << " events where all top-radiated partons are matched";
	edm::LogVerbatim  ("MainResults") << "  (37) \t\t(37)/(25) "<< eff16<< "% or " << QuarksMatchedUnmatchedTopRadiation << " events where not all top-radiated partons are matched (there can be more than 1 unmatched parton)";
	edm::LogVerbatim  ("MainResults") << "  (37b) \t\t"<< QuarksMatchedUnmatchedTopRadiationtotal << "	total unmatched TopRadiation partons";
	edm::LogVerbatim  ("MainResults") << "  (38) \t\t\t" <<"of which " << QuarksMatchedUnmatchedRadCategoryEta << " because of the eta requirement";
	edm::LogVerbatim  ("MainResults") << "  (39) \t\t\t" <<"of which " << QuarksMatchedUnmatchedRadCategoryPt << " because of the p_{T} requirement";
	edm::LogVerbatim  ("MainResults") << "  (40) \t\t\t" <<"of which " << QuarksMatchedUnmatchedRadCategoryDR << " because of the DR requirement";
	edm::LogVerbatim  ("MainResults") << "  (41) \t\t\t" <<"of which " << QuarksMatchedUnmatchedRadNotCategorized << " not categorized";		
	edm::LogVerbatim  ("MainResults") << "  events with matched quarks";
	edm::LogVerbatim  ("MainResults") << "  first 4 numbers should be the same as (26)-(30)";
	edm::LogVerbatim  ("MainResults") << "  (42) \t (42)/(25) " <<eff17<< "% or " <<QuarksMatchedNoTopRadiation << " events with ISR but without top radiation";
	edm::LogVerbatim  ("MainResults") << "  (43) \t (43)/(25) " <<eff18<< "% or " <<QuarksMatchedISRTopRadiation << " events with ISR and top radiation ";
	edm::LogVerbatim  ("MainResults") << "  (44) \t (44)/(25) " <<eff19<< "% or " <<QuarksMatchedTopRadiationNoISR << " events with top radiation but without ISR ";
	edm::LogVerbatim  ("MainResults") << "  (45) \t (45)/(25) " <<eff20<< "% or " <<QuarksMatchedNoRadiation << " events without ISR or top radiation";
	edm::LogVerbatim  ("MainResults") << "  if there is radiation and more then 4 jets: ";
	edm::LogVerbatim  ("MainResults") << "  (46) \t (46)/(25) " <<eff21<< "% or " <<QuarksMatchedHighest4NoRad << " events in which the 4 highest pT jets are not matched with partons from	radiation";
	edm::LogVerbatim  ("MainResults") << "  (47) \t (47)/(25) " <<eff22<< "% or " <<QuarksMatchedHighest4AtLeast1Rad << " events in which at least 1 of the 4 highest pT jets are matched with partons from	radiation";
	edm::LogVerbatim  ("MainResults") << "  (48) \t (48)/(25) " <<eff23<< "% or " <<QuarksMatchedHighest4AtLeast2Rad << " events in which at least 2 of the 4 highest pT jets are matched with partons from	radiation";
	edm::LogVerbatim  ("MainResults") << "  (49) \t (49)/(25) " <<eff24<< "% or " <<QuarksMatchedHighest4AtLeast3Rad << " events in which at least 3 of the 4 highest pT jets are matched with partons from	radiation";
	edm::LogVerbatim  ("MainResults") << "  (50) \t (50)/(25) " <<eff25<< "% or " <<QuarksMatchedHighest4AtLeast4Rad << " events in which 4 of the 4 highest pT jets	are matched with partons from	radiation\n\n";

	edm::LogVerbatim  ("MainResults") << "  (51) (51)/(11) " <<eff26<< "% or " << UnmatchedQuarks << " events where not all quarks are matched (there can be more than	1 unmatched parton)";
	edm::LogVerbatim  ("MainResults") << "  (52) \t" <<"of which " << UnmatchedQuarkCategoryEta << " because of the eta requirement";
	edm::LogVerbatim  ("MainResults") << "  (53) \t" <<"of which " << UnmatchedQuarkCategoryPt << " because of the p_{T} requirement";
	edm::LogVerbatim  ("MainResults") << "  (54) \t" <<"of which " << UnmatchedQuarkCategoryDR << " because of the DR requirement";
	edm::LogVerbatim  ("MainResults") << "  (55) \t" <<"of which " << UnmatchedQuarkNotCategorized << " not categorized \n";

	edm::LogVerbatim  ("MainResults") << "  (56) \t(56)/(51) " <<"of which " << eff27<< "% or " << UnmatchedQuark1CategoryEta << " because of the eta requirement (look at 1 of unmatched	quarks)";
	edm::LogVerbatim  ("MainResults") << "  (57) \t(57)/(51) " <<"of which " << eff28<< "% or " << UnmatchedQuark1CategoryPt << " because of the p_{T} requirement (look at 1 of unmatched	quarks)";
	edm::LogVerbatim  ("MainResults") << "  (58) \t(58)/(51) " <<"of which " << eff29<< "% or " << UnmatchedQuark1CategoryDR << " because of the DR requirement (look at 1 of unmatched	quarks)";
	edm::LogVerbatim  ("MainResults") << "  (59) \t(59)/(51) " <<"of which " << eff30<< "% or " << UnmatchedQuark1NotCategorized << " not categorized (look at 1 of	unmatched	quarks) \n";
	edm::LogVerbatim  ("MainResults") << "  (60) (60)/(11) " << eff31<< "% or "  <<Only1UnmatchedQuarks << " events with exactly 1 unmatched quark";
	edm::LogVerbatim  ("MainResults") << "  (61) \t(61)/(60) " <<"of which " << eff32<< "% or "  <<Only1UnmatchedQuarkCategoryEta << " because of the eta requirement";
	edm::LogVerbatim  ("MainResults") << "  (62) \t(62)/(60) " <<"of which " << eff33<< "% or "  <<Only1UnmatchedQuarkCategoryPt << " because of the p_{T} requirement";
	edm::LogVerbatim  ("MainResults") << "  (63) \t(63)/(60) " <<"of which " << eff34<< "% or "  <<Only1UnmatchedQuarkCategoryDR << " because of the DR requirement";
	edm::LogVerbatim  ("MainResults") << "  (64) \t(64)/(60) " <<"of which " << eff35<< "% or "  <<Only1UnmatchedQuarkNotCategorized	<< " not categorized\n";
	edm::LogVerbatim  ("MainResults") << "  (65) \t\t(65)/(51) " << eff36<< "% or " << MatchedISR << " events where all ISR partons are matched";	
	edm::LogVerbatim  ("MainResults") << "  (66) \t\t(66)/(51) " << eff37<< "% or " <<  UnmatchedISR << " events where not all ISR partons are matched (there can be more than 1	unmatched parton)";
	edm::LogVerbatim  ("MainResults") << "  (67) \t\t\t" <<"of which " << UnmatchedISRCategoryEta << " because of the eta requirement";
	edm::LogVerbatim  ("MainResults") << "  (68) \t\t\t" <<"of which " << UnmatchedISRCategoryPt << " because of the p_{T} requirement";
	edm::LogVerbatim  ("MainResults") << "  (69) \t\t\t" <<"of which " << UnmatchedISRCategoryDR << " because of the DR requirement";
	edm::LogVerbatim  ("MainResults") << "  (70) \t\t\t" <<"of which " << UnmatchedISRNotCategorized << " not categorized\n";	
	edm::LogVerbatim  ("MainResults") << "  (71) \t\t(71)/(51) " <<eff38<< "% or " <<  MatchedTopRadiation << " events where all top-radiated partons are matched";
	edm::LogVerbatim  ("MainResults") << "  (72) \t\t(72)/(51) " <<eff39<< "% or " <<  UnmatchedTopRadiation << " events where not all top-radiated partons are matched (there can be more than 1 unmatched parton)";
	edm::LogVerbatim  ("MainResults") << "  (73) \t\t\t" <<"of which " << UnmatchedRadCategoryEta << " because of the eta requirement";
	edm::LogVerbatim  ("MainResults") << "  (74) \t\t\t" <<"of which " << UnmatchedRadCategoryPt << " because of the p_{T} requirement";
	edm::LogVerbatim  ("MainResults") << "  (75) \t\t\t" <<"of which " << UnmatchedRadCategoryDR << " because of the DR requirement";
	edm::LogVerbatim  ("MainResults") << "  (76) \t\t\t" <<"of which " << UnmatchedRadNotCategorized << " not categorized";		
	edm::LogVerbatim  ("MainResults") << "  events with unmatched quarks";
	edm::LogVerbatim  ("MainResults") << "  (77) \t(77)/(51) " <<eff40<< "% or " <<  NoTopRadiation << " events with ISR but without top radiation";
	edm::LogVerbatim  ("MainResults") << "  (78) \t(78)/(51) " <<eff41<< "% or " <<  ISRTopRadiation << " events with ISR and top radiation ";
	edm::LogVerbatim  ("MainResults") << "  (79) \t(79)/(51) " <<eff42<< "% or " <<  TopRadiationNoISR << " events with top radiation but without ISR ";
	edm::LogVerbatim  ("MainResults") << "  (80) \t(80)/(51) " <<eff43<< "% or " <<  NoRadiation << " events without ISR or top radiation";
	edm::LogVerbatim  ("MainResults") << "  (81) \t\t of which (81)/(80) " <<eff44<< "% or " <<  UnmatchedQuark1NoRadCategoryEta << " because of eta requirement (look	at 1 of the unmatched quarks)";
	edm::LogVerbatim  ("MainResults") << "  (82) \t\t of which (82)/(80) " <<eff45<< "% or " <<  UnmatchedQuark1NoRadCategoryPt << " because of p_{T} requirement (look	at 1 of the unmatched quarks)";
	edm::LogVerbatim  ("MainResults") << "  (83) \t\t of which (83)/(80) " <<eff46<< "% or " <<  UnmatchedQuark1NoRadCategoryDR << " because of DR requirement (look	at 1 of the unmatched quarks)";
	edm::LogVerbatim  ("MainResults") << "  (84) \t\t of which (84)/(80) " <<eff47<< "% or " <<  UnmatchedQuark1NoRadNotCategorized << " not categorized (look	at 1 of the unmatched quarks)";
	edm::LogVerbatim  ("MainResults") << "  (85) \t (85)/(51)" <<eff48<< "% or " << Highest4NoRad << " events in which the 4 highest pT jets are not matched with partons from	radiation";
	edm::LogVerbatim  ("MainResults") << "  (86) \t (86)/(51)" <<eff49<< "% or " << Highest4AtLeast1Rad << " events in which at least 1 of the 4 highest pT jets are matched with partons from	radiation";
	edm::LogVerbatim  ("MainResults") << "  (87) \t (87)/(51)" <<eff50<< "% or " << Highest4AtLeast2Rad << " events in which at least 2 of the 4 highest pT jets are matched with partons from	radiation";
	edm::LogVerbatim  ("MainResults") << "  (88) \t (88)/(51)" <<eff51<< "% or " << Highest4AtLeast3Rad << " events in which at least 3 of the 4 highest pT jets are matched with partons from	radiation";
	edm::LogVerbatim  ("MainResults") << "  (89) \t (89)/(51)" <<eff52<< "% or " << Highest4AtLeast4Rad << " events in which 4 of the 4 highest pT jets are matched with partons from	radiation\n";

  edm::LogProblem ("SummaryError") << " \n\n";
  edm::LogProblem ("SummaryError") << " ------------------------------------";
  edm::LogProblem ("SummaryError") << " ------------------------------------";
  edm::LogProblem ("SummaryError") << " -- Report from TruthReco module   --";
  edm::LogProblem ("SummaryError") << " ------------------------------------";
  edm::LogProblem ("SummaryError") << " ------------------------------------";

	edm::LogProblem ("SummaryError")  << "  " << eff7<< "% or "<<NrNotSemiMu <<" events no ttsemimu generated event "; 
	edm::LogProblem ("SummaryError")  << "  " << genObjectsNULL << " events where at least 1 of the objects in the semi-mu genevent is zero ";	




}


//define this as a plug-in
DEFINE_FWK_MODULE(TruthReco);
