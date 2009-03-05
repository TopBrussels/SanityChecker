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
			double muonRelIso_,muonECALVetoConeEt_,muonHCALVetoConeEt_,muonPt_,jetPt_,muonEta_,jetEta_,muonMinDR_;
			bool doNonOverlappingJets_;
			double nonOverlappingJetsMinDR_;
			bool calculateJECinbins_;
			std::vector<double> etabinVals_,pTbinVals_;
			int NJets_;
			double minDRunmatchedParton_;
			
			int NbSemiMuRecoEvents,NbSemiMuGenEvents,genObjectsNotNULL,genObjectsNULL,FourQuarksNotMatched,EventNotMatched,ptnrbins,etanrbins;
			int UnmatchedQuarks,MatchedQuarks,MatchedISR,UnmatchedISR,UnmatchedTopRadiation,MatchedTopRadiation,NoTopRadiation,NoISR,NoRadiation;
			int Highest4NoRad,Highest4AtLeast1Rad,Highest4AtLeast2Rad,Highest4AtLeast3Rad,Highest4AtLeast4Rad;
			int UnmatchedQuarkCategoryEta,UnmatchedQuarkCategoryPt,UnmatchedQuarkCategoryDR,UnmatchedQuarkNotCategorized;
			int UnmatchedISRCategoryEta,UnmatchedISRCategoryPt,UnmatchedISRCategoryDR,UnmatchedISRNotCategorized;
			int UnmatchedRadCategoryEta,UnmatchedRadCategoryPt,UnmatchedRadCategoryDR,UnmatchedRadNotCategorized;
			int NrNotSemiMu,NrAnalyzed,NrNot4Jets,NrNot1Mu,NrNot1MET,NrNot4SelJets,NrNot1SelMu,NrNot1SelIsoMu;
			unsigned int NumberQuarks,NumberISR,NumberRadiationLep,NumberRadiationHad,NumberRadiation;	
			
			//Histograms are booked in the beginJob() method
  		TH1F 	*hExpLightJECPtEtaBin[20][20],*hExpBJECPtEtaBin[20][20];
  		TF1 	*fExpLightJECPtEtaBin[20][20],*fExpBJECPtEtaBin[20][20];
  		TH1F 	*hExpLightJECIncl,*hExpBJECIncl;
  		TF1 	*fExpLightJECIncl,*fExpBJECIncl;
			TH1F 	*LowestPtMatchedJet_Rank,	*LowestPtMatchedJet_Pt,*HighestPtRadMatchedJet_Rank,*HighestPtRadMatchedJet_Pt;
			TH1F 	*mWHad_Gen,*mWLep_Gen,*mtHad_Gen,*mtLep_Gen,*mWHad_Rec,*mWLep_Rec,*mtHad_Rec,*mtLep_Rec; 	
			TH2F 	*minJetMuVSminQuarkMu,*minJetsVSminQuarks;	
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

//cout <<NrAnalyzed << endl;

	
	if(genEvt.isSemiLeptonic(genEvt.kMuon)) NbSemiMuGenEvents++;
	if(jets.size()>3 && muons.size()>0 && mets.size()>0) NbSemiMuRecoEvents++;
	
	bool allok = true;
	bool muonselected = true;
	
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
	
  std::vector<pat::Jet> myjetsetagood;
  std::vector<pat::Jet> myselectedjets;
  for(unsigned int i =0;i<jets.size();i++) {if(fabs(jets[i].eta())<jetEta_) myjetsetagood.push_back(jets[i]);}
	for(unsigned int i =0;i<myjetsetagood.size();i++) {if(jets[i].pt()>jetPt_ ) myselectedjets.push_back(myjetsetagood[i]);}		

	
	//check if at least 4 jets with eta ok and pT exceeding jetPt threshold!
	if(myselectedjets.size()<(unsigned)NJets_){ 
		edm::LogWarning  ("NoDataFound") << "Less then 4 jets with p_{T}>"<<jetPt_<<"GeV and |#eta| <"<<jetEta_<<"\n"; 
		NrNot4SelJets++;
		allok = false;		
	}

	//check if muon passes eta/pT cut!
	if(muons.size()<1 || muons[0].pt()<muonPt_ || fabs(muons[0].eta())>muonEta_){
		edm::LogWarning  ("NoDataFound") << "The leading muon has no p_{T}>"<<muonPt_<<"GeV or |#eta| >"<<muonEta_<<"\n"; 
		NrNot1SelMu++;
		allok = false;			
		muonselected = false;
	}

	double RelIsoMuon = 0.;
	if(muons.size()>0){
		RelIsoMuon = muons[0].pt()/(muons[0].pt()+muons[0].isolationR03().sumPt+muons[0].isolationR03().emEt+muons[0].isolationR03().hadEt);
	
		//veto cone size is 0.07 in ECAL and 0.1 in HCAL, proposed ET cuts are resp 4 and 6
		const reco::IsoDeposit * ecalIsoDep = (muons)[0].ecalIsoDeposit();
		const reco::IsoDeposit * hcalIsoDep = (muons)[0].hcalIsoDeposit();

		//if muonisolation is required check if muon passes reliso and veto cone energy requirements!
		if(doMuonIsolation_ && muonselected && (RelIsoMuon<muonRelIso_ || ecalIsoDep->candEnergy()>muonECALVetoConeEt_ || hcalIsoDep->candEnergy()>muonHCALVetoConeEt_)){
			edm::LogWarning  ("NoDataFound") << "The leading muon does not pass the isolation requirements \n"; 
			NrNot1SelIsoMu++;
			allok = false;			
		}
	
	}
	if(!useMatchingFromPAT_ && allok){ 
//cout << "entered main part"<< endl;
  	// Matching index : Hadronic Q  = 0, Hadronic Q' = 1, Hadronic b  = 2, Leptonic b  = 3;
  	std::vector<const reco::Candidate *> quarks;
		const reco::Candidate * muon;
		const reco::Candidate * neutrino;	

		if(genEvt.hadronicDecayQuark() != 0 && genEvt.hadronicDecayQuarkBar() != 0  && genEvt.hadronicDecayB() != 0 
		&& genEvt.leptonicDecayB() != 0 && genEvt.singleLepton() != 0 && genEvt.singleNeutrino() != 0){ 		
		
		genObjectsNotNULL++; 
//cout << genObjectsNotNULL<< endl;

  	quarks.push_back(genEvt.hadronicDecayQuark());
  	quarks.push_back(genEvt.hadronicDecayQuarkBar());
  	quarks.push_back(genEvt.hadronicDecayB());
  	quarks.push_back(genEvt.leptonicDecayB());

	 	muon = genEvt.singleLepton();
  	neutrino = genEvt.singleNeutrino();

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////// MATCH JET WITH QUARKS (NO RADIATION)///////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////		
		JetPartonMatching *GenMatch = new JetPartonMatching(quarks, myselectedjets, matchingAlgo_, useMaxDist_, useDeltaR_, maxDist_);
		//by construction of JetPartonMatching 1 jet can not match 2 different partons			
		//now avoid that 2 jets are matched with the same parton: there should be no unmatched partons
//cout << "matched with quarks"<< endl;
		if(GenMatch->getNumberOfUnmatchedPartons(0)!=0){
			edm::LogWarning  ("NoDataFound") << "One of the quarks is not matched!!! \n"; 
			throw cms::Exception("NotFound") << "One of the quarks is not matched!!!"; 
			FourQuarksNotMatched++;
   	}
		if((deltaR(muon->p4(),muons[0].p4()) > muonMinDR_)){
			edm::LogWarning  ("NoDataFound") << "The muon is not matched!!! \n"; 
			throw cms::Exception("NotFound") << "The muon is not matched!!!"; 
			EventNotMatched++;
	 	}
   	std::vector< double > matchjetpt;
   	std::vector< int > matchjetidx;
//cout << "lala"<< endl;

		for(unsigned int i=0; i<quarks.size(); i++){
			///////////////////// ---start --- calculating expected JEC /////////////////////
			double quarkEnergy   = quarks[i]->energy(); 
    	double jetEnergy     = myselectedjets[GenMatch->getMatchForParton(i,0)].energy(); 
      if(i<2)hExpLightJECIncl -> Fill((quarkEnergy-jetEnergy)/jetEnergy);
      if(i>=2)hExpBJECIncl -> Fill((quarkEnergy-jetEnergy)/jetEnergy);
  		
			if(calculateJECinbins_){
				// find eta and pt bin
      	int etabin  =  0;
      	if(etanrbins > 1){
        	for(unsigned int b=0; b<etabinVals_.size()-1; b++) {
          	if(fabs(myselectedjets[GenMatch->getMatchForParton(i,0)].eta()) > etabinVals_[b]) etabin = b;
        	}
      	}     
      	int ptbin  =  0;
      	for(unsigned int b=0; b<pTbinVals_.size()-1; b++) {
        	if(myselectedjets[GenMatch->getMatchForParton(i,0)].pt() > pTbinVals_[b]) ptbin = b;
      	}
      	if(i<2)hExpLightJECPtEtaBin[etabin][ptbin] -> Fill((quarkEnergy-jetEnergy)/jetEnergy);
      	if(i>=2)hExpBJECPtEtaBin[etabin][ptbin] -> Fill((quarkEnergy-jetEnergy)/jetEnergy);
			}
			///////////////////// ---done--- calculating expected JEC /////////////////////
//	cout << "lala1"<< endl;
		
			int Idx = GenMatch->getMatchForParton(i,0);
			double jetpt = myselectedjets[Idx].pt();
			matchjetpt.push_back(jetpt);
			matchjetidx.push_back(Idx);

		}

   	//the highest rank among the indexes of the jets stored (i.e. the ones matching with the partons), 
		//corresponds to the jet with the lowest pt 
   	if(matchjetidx.size()!=0){
   		std::sort(matchjetidx.begin(), matchjetidx.end(), std::greater<int>() );
   		std::sort(matchjetpt.begin(), matchjetpt.end(), std::greater<double>() );
			LowestPtMatchedJet_Rank->Fill(matchjetidx.front());
   		LowestPtMatchedJet_Pt  ->Fill(matchjetpt.back());
   	}
//	cout << "lala2"<< endl;

		///////////////////// ---start --- calculating W and top masses /////////////////////
		math::XYZTLorentzVector WHadGen = quarks[0]->p4()+quarks[1]->p4();
		math::XYZTLorentzVector WLepGen = muon->p4()+neutrino->p4();
		math::XYZTLorentzVector tHadGen = quarks[0]->p4()+quarks[1]->p4()+quarks[2]->p4();
		math::XYZTLorentzVector tLepGen = muon->p4()+neutrino->p4()+quarks[3]->p4();
		math::XYZTLorentzVector WHadRec = myselectedjets[GenMatch->getMatchForParton(0,0)].p4()+myselectedjets[GenMatch->getMatchForParton(1,0)].p4();
		math::XYZTLorentzVector WLepRec = muons[0].p4()+mets[0].p4();
		math::XYZTLorentzVector tHadRec = myselectedjets[GenMatch->getMatchForParton(0,0)].p4()+myselectedjets[GenMatch->getMatchForParton(1,0)].p4()+myselectedjets[GenMatch->getMatchForParton(2,0)].p4();
		math::XYZTLorentzVector tLepRec = muons[0].p4()+mets[0].p4()+myselectedjets[GenMatch->getMatchForParton(3,0)].p4();
		mWHad_Gen -> Fill(WHadGen.mass());
		mWLep_Gen -> Fill(WLepGen.mass());
		mtHad_Gen -> Fill(tHadGen.mass());
		mtLep_Gen -> Fill(tLepGen.mass());
		mWHad_Rec -> Fill(WHadRec.mass());
		mWLep_Rec -> Fill(WLepRec.mass());
		mtHad_Rec -> Fill(tHadRec.mass());
		mtLep_Rec -> Fill(tLepRec.mass());
		///////////////////// ---done --- calculating W and top masses /////////////////////

//	cout << "lala3"<< endl;
		///////////////////// ---start --- calculating smallest angle jets/mu and jets /////////////////////
		double smallestDRjetmu = 1000.; 
		double DRjetmu[4];
		double smallestDRquarkmu = 1000.; 
		double DRquarkmu[4];
		for(int i = 0; i<4; i++){
			DRjetmu[i] = deltaR(myselectedjets[GenMatch->getMatchForParton(i,0)].p4(),muons[0].p4());
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
				DRjets[k] = deltaR(myselectedjets[GenMatch->getMatchForParton(i,0)].p4(),myselectedjets[GenMatch->getMatchForParton(j,0)].p4());
				DRquarks[k] = deltaR(quarks[i]->p4(),quarks[j]->p4());
	  		if(DRjets[k] < smallestDRjets) smallestDRjets = DRjets[k];
	  		if(DRquarks[k] < smallestDRquarks) smallestDRquarks = DRquarks[k];
				k=k+1;
			}
		}
		minJetsVSminQuarks->Fill(smallestDRjets,smallestDRquarks);
		///////////////////// ---done --- calculating smallest angle jets/mu and jets /////////////////////
//	cout << "lala4"<< endl;
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////// END MATCH JET WITH QUARKS (NO RADIATION)///////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////		

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////// INCLUDE RADIATION ////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

		std::vector<const reco::Candidate *> vectorISR, TopRadiation;
		std::vector<pair<int,pat::Jet> > patJetsMatchedWithRadiation;
		std::vector< int > matchquarkjetidx,matchISRjetidx,matchRadiationjetidx,originalrankjetsISR,originalrankjetsradiation;
		std::vector< int > originalrankjetsISRinput,originalrankjetsradiationinput,rankjetsmatchedtorad;
		std::vector<pat::Jet> myselectedjetsISR,myselectedjetsradiation,myselectedjetsradiationNoISR,FourHighestPtJets;
		std::vector< double > ISRpT,TopRadiationpT,TopRadiationpTNoISR,ptjetsmatchedtorad;
		std::vector<const reco::Candidate *> NfirstpTISR,NfirstpTradiation,NfirstpTradiationNoISR;
		
		///check number quarks
		unsigned int NumberQuarks = quarks.size();		
		///check number additional ISR partons
		unsigned int NumberISR = genEvt.ISR().size();		
		///check number additional radiation partons
		unsigned int NumberRadiationLep = genEvt.leptonicDecayTopRadiation().size();		
		unsigned int NumberRadiationHad = genEvt.hadronicDecayTopRadiation().size();		
		unsigned int NumberRadiation = NumberRadiationLep+NumberRadiationHad;		
		
		if(NumberISR!=0) {
//	cout << "lala5"<< endl;
			for(unsigned int i=NumberQuarks; i<(NumberISR+NumberQuarks); i++){
				vectorISR.push_back((genEvt.ISR())[i]);
			}
		}
		if(NumberRadiationLep!=0 ) {
			for(unsigned int i=(NumberISR+NumberQuarks);i<(NumberISR+NumberQuarks+NumberRadiationLep); i++){
				TopRadiation.push_back((genEvt.leptonicDecayTopRadiation())[i]);
			}
		}
		if(NumberRadiationHad!=0 ) {
			for(unsigned int i=(NumberISR+NumberQuarks+NumberRadiationLep); i<(NumberISR+NumberQuarks+NumberRadiation); i++){
				TopRadiation.push_back((genEvt.hadronicDecayTopRadiation())[i]);
			}
		}
		
		double DRquarki[4];
		double DRISRi[20];
		double DRRadi[20];
		double smallestDRquarki;
		double smallestDRISRi;
		double smallestDRRadi;
		if(myselectedjets.size()>=4){
//	cout << "lala6"<< endl;
			
			//match the quarks
			JetPartonMatching *GenMatchQuarks = new JetPartonMatching(quarks, myselectedjets, matchingAlgo_, useMaxDist_, useDeltaR_, maxDist_);
			
			//check if all quarks are matched ---> if so, we are not interested in radiation
			if(GenMatchQuarks->getNumberOfUnmatchedPartons(0)== 0){
				MatchedQuarks++;				
			
			//check if unmatched quarks ---> if so, categorize, look at radiation
			}else if(GenMatchQuarks->getNumberOfUnmatchedPartons(0)!= 0){ //unmatched quarks
				UnmatchedQuarks++;
				
				//find parton which is not matched -> check eta/pT/DR ---> count different possibilities
				for(unsigned int i = 0; i<quarks.size(); i++){
					int Idx = GenMatchQuarks->getMatchForParton(i,0);
					if(Idx<0){
						smallestDRquarki = 0;
						for(unsigned int j = 0; j<quarks.size(); i++){
							if(j!=i){
								DRquarki[j] = deltaR(quarks[i]->p4(),quarks[j]->p4());
	  						if(DRquarki[j] < smallestDRquarki) smallestDRquarki = DRquarki[j];
							}
						}
						if(fabs(quarks[i]->eta())>jetEta_){
							UnmatchedQuarkCategoryEta++;
						}else if(quarks[i]->pt()<jetPt_){
							UnmatchedQuarkCategoryPt++;
						}else if(smallestDRquarki<minDRunmatchedParton_){
							UnmatchedQuarkCategoryDR++;
						}else{
							UnmatchedQuarkNotCategorized++;
						}
					}
				}
				for(unsigned int i = 0; i < quarks.size(); i++){			
					int IdxQuarkJets = GenMatchQuarks->getMatchForParton(i,0);
					matchquarkjetidx.push_back(IdxQuarkJets); //store the rank of the jet (if parton is unmatched ---> rank is negative)
				}
				
				/////////////////////////////////
				////////// LOOK AT ISR //////////
				/////////////////////////////////		
				if( NumberISR !=0){
//		cout << "lala7"<< endl;
				//make new list of jets to be used as input for matching (remove quark-matched jets from jetlist used as input for quark matching)
					//make this new list of jets, but also keep the original rank of the jet! ---> this rank can be used for plots and numbers
					for(unsigned int i = 0; i<myselectedjets.size() ; i++){
						for(unsigned int k =0; k<quarks.size(); k++){
							if(i != (unsigned int) matchquarkjetidx[k]) {
								myselectedjetsISR.push_back(myselectedjets[i]);
								originalrankjetsISRinput.push_back(i);
							}
						}
					}
					
				  ///////////////////// #ISR > myselectedjetsISR.size() /////////////////////
					if(NumberISR > myselectedjetsISR.size()){
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
						JetPartonMatching *GenMatchISR = new JetPartonMatching(NfirstpTISR, myselectedjetsISR, matchingAlgo_, useMaxDist_, useDeltaR_, maxDist_);
					
						//check if all partons are matched ---> if so, just store jet ranks, nothing else (no jets left)
						if(GenMatchISR->getNumberOfUnmatchedPartons(0)== 0){
							MatchedISR++;				
							
							//store the rank of the jet for later use (if parton is unmatched ---> rank is negative)
 							for(unsigned int i = 0; i < NfirstpTISR.size(); i++){			
								int IdxISRJets = GenMatchISR->getMatchForParton(i,0);
								matchISRjetidx.push_back(IdxISRJets);
								originalrankjetsISR.push_back(originalrankjetsISRinput[IdxISRJets]);
							}
							
						//check if unmatched partons ---> if so, categorize, store jet ranks, nothing else (no jets left)
						}else if(GenMatchISR->getNumberOfUnmatchedPartons(0)!= 0){ //unmatched quarks
							UnmatchedISR++;
						
							//find parton which is not matched -> check eta/pT/DR ---> count different possibilities
							for(unsigned int i = 0; i<NfirstpTISR.size(); i++){
								int Idx = GenMatchISR->getMatchForParton(i,0);
								if(Idx<0){
									smallestDRISRi = 0;
									for(unsigned int j = 0; j<NfirstpTISR.size(); i++){
										if(j!=i){
											DRISRi[j] = deltaR(NfirstpTISR[i]->p4(),NfirstpTISR[j]->p4());
	  									if(DRISRi[j] < smallestDRISRi) smallestDRISRi = DRISRi[j];
										}
									}
									if(fabs(NfirstpTISR[i]->eta())>jetEta_){
										UnmatchedISRCategoryEta++;
									}else if(NfirstpTISR[i]->pt()<jetPt_){
										UnmatchedISRCategoryPt++;
									}else if(smallestDRISRi<minDRunmatchedParton_){
										UnmatchedISRCategoryDR++;
									}else{
										UnmatchedISRNotCategorized++;
									}
								}
							}
							
 							//store the rank of the jet for later use (if parton is unmatched ---> rank is negative)
 							for(unsigned int i = 0; i < NfirstpTISR.size(); i++){			
								int IdxISRJets = GenMatchISR->getMatchForParton(i,0);
								matchISRjetidx.push_back(IdxISRJets);
								originalrankjetsISR.push_back(originalrankjetsISRinput[IdxISRJets]);
							}
						}	
						
				  ///////////////////// #ISR <= myselectedjetsISR.size() /////////////////////					
          }else if( NumberISR <= myselectedjetsISR.size()){
						JetPartonMatching *GenMatchISR = new JetPartonMatching(vectorISR, myselectedjetsISR, matchingAlgo_, useMaxDist_, useDeltaR_, maxDist_);
						
						//check if unmatched ISR  ---> if so, categorize, store jet ranks, look at radiation (still jets left)
						if(GenMatchISR->getNumberOfUnmatchedPartons(0)!= 0){
							UnmatchedISR++; 

							//find parton which is not matched -> check eta/pT/DR ---> count different possibilities
							for(unsigned int i = 0; i<NumberISR; i++){
								int Idx = GenMatchISR->getMatchForParton(i,0);
								if(Idx<0){
									smallestDRISRi = 0;
									for(unsigned int j = 0; j<NumberISR; i++){
										if(j!=i){
											DRISRi[j] = deltaR(vectorISR[i]->p4(),vectorISR[j]->p4());
	  									if(DRISRi[j] < smallestDRISRi) smallestDRISRi = DRISRi[j];
										}
									}
									if(fabs(vectorISR[i]->eta())>jetEta_){
										UnmatchedISRCategoryEta++;
									}else if(vectorISR[i]->pt()<jetPt_){
										UnmatchedISRCategoryPt++;
									}else if(smallestDRISRi<minDRunmatchedParton_){
										UnmatchedISRCategoryDR++;
									}else{
										UnmatchedISRNotCategorized++;
									}
								}
							}
 							
							//store the rank of the jet for later use (if parton is unmatched ---> rank is negative)
 							for(unsigned int i = 0; i < NumberISR; i++){			
								int IdxISRJets = GenMatchISR->getMatchForParton(i,0);
								matchISRjetidx.push_back(IdxISRJets);
								originalrankjetsISR.push_back(originalrankjetsISRinput[IdxISRJets]);
							}

							//////////////////////////////////////////////////////////////////////////////
							////////// LOOK AT TOPRADIATION - JETS LEFT - ISR NOT ALL MATCHED - //////////
							//////////////////////////////////////////////////////////////////////////////
							if(NumberRadiation !=0){ 
								//make new list of jets to be used as input for matching (remove ISR-matched jets from jetlist used as input for ISR)
								for(unsigned int i = 0; i<myselectedjets.size() ; i++){
									for(unsigned int l =0; l<quarks.size(); l++){
										for(unsigned int k =0; k<NumberISR; k++){
											if(i != (unsigned int) originalrankjetsISR[k] && i != (unsigned int) matchquarkjetidx[l]){
												myselectedjetsradiation.push_back(myselectedjets[i]);
												originalrankjetsradiationinput.push_back(i);
											}
										}
									}
								}
								
				  			///////////////////// #top radiation > myselectedjetsradiation.size() /////////////////////
								if(NumberRadiation > myselectedjetsradiation.size()){
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
									JetPartonMatching *GenMatchRadiation = new JetPartonMatching(NfirstpTradiation, myselectedjetsradiation, matchingAlgo_, useMaxDist_, useDeltaR_, maxDist_);

									//check if all partons are matched ---> if so, just store jet ranks, nothing else (no jets left)
									if(GenMatchRadiation->getNumberOfUnmatchedPartons(0)== 0){
										MatchedTopRadiation++;				
							
										//store the rank of the jet for later use (if parton is unmatched ---> rank is negative)
 										for(unsigned int i = 0; i < NfirstpTradiation.size(); i++){			
											int IdxRadJets = GenMatchRadiation->getMatchForParton(i,0);
											matchRadiationjetidx.push_back(IdxRadJets);
											originalrankjetsradiation.push_back(originalrankjetsradiationinput[IdxRadJets]);
										}
							
									//check if unmatched partons ---> if so, categorize, store jet ranks, nothing else (no jets left)
									}else if(GenMatchRadiation->getNumberOfUnmatchedPartons(0)!= 0){ 
										UnmatchedTopRadiation++;
						
										//find parton which is not matched -> check eta/pT/DR ---> count different possibilities
										for(unsigned int i = 0; i<NfirstpTradiation.size(); i++){
											int Idx = GenMatchRadiation->getMatchForParton(i,0);
											if(Idx<0){
												smallestDRRadi = 1000;
												for(unsigned int j = 0; j<NfirstpTradiation.size(); i++){
													if(j!=i){
														DRRadi[j] = deltaR(NfirstpTradiation[i]->p4(),NfirstpTradiation[j]->p4());
	  												if(DRRadi[j] < smallestDRRadi) smallestDRRadi = DRRadi[j];
													}
												}
												if(fabs(NfirstpTradiation[i]->eta())>jetEta_){
													UnmatchedRadCategoryEta++;
												}else if(NfirstpTradiation[i]->pt()<jetPt_){
													UnmatchedRadCategoryPt++;
												}else if(smallestDRRadi<minDRunmatchedParton_){
													UnmatchedRadCategoryDR++;
												}else{
													UnmatchedRadNotCategorized++;
												}
											}
										}
							
 										//store the rank of the jet for later use (if parton is unmatched ---> rank is negative)
 										for(unsigned int i = 0; i < NfirstpTradiation.size(); i++){			
											int IdxRadJets = GenMatchRadiation->getMatchForParton(i,0);
											matchRadiationjetidx.push_back(IdxRadJets);
											originalrankjetsradiation.push_back(originalrankjetsradiationinput[IdxRadJets]);
										}
									}

				  			///////////////////// #top radiation <= myselectedjetsradiation.size() /////////////////////
          			}else if( NumberRadiation <= myselectedjetsradiation.size()){
									JetPartonMatching *GenMatchRadiation = new JetPartonMatching(TopRadiation, myselectedjetsradiation, matchingAlgo_, useMaxDist_, useDeltaR_, maxDist_);
									
									//check if unmatched top radiation  ---> if so, categorize, store jet ranks 
									if(GenMatchRadiation->getNumberOfUnmatchedPartons(0)!= 0){
										UnmatchedTopRadiation++;

										//find parton which is not matched -> check eta/pT/DR ---> count different possibilities
										for(unsigned int i = 0; i<NumberRadiation; i++){
											int Idx = GenMatchRadiation->getMatchForParton(i,0);
											if(Idx<0){
												smallestDRRadi = 0;
												for(unsigned int j = 0; j<NumberRadiation; i++){
													if(j!=i){
														DRRadi[j] = deltaR(TopRadiation[i]->p4(),TopRadiation[j]->p4());
	  												if(DRRadi[j] < smallestDRRadi) smallestDRRadi = DRRadi[j];
													}
												}
												if(fabs(TopRadiation[i]->eta())>jetEta_){
													UnmatchedRadCategoryEta++;
												}else if(TopRadiation[i]->pt()<jetPt_){
													UnmatchedRadCategoryPt++;
												}else if(smallestDRRadi<minDRunmatchedParton_){
													UnmatchedRadCategoryDR++;
												}else{
													UnmatchedRadNotCategorized++;
												}
											}
										}
 							
										//store the rank of the jet for later use (if parton is unmatched ---> rank is negative)
 										for(unsigned int i = 0; i < NumberRadiation; i++){			
											int IdxRadJets = GenMatchRadiation->getMatchForParton(i,0);
											matchRadiationjetidx.push_back(IdxRadJets);
											originalrankjetsradiation.push_back(originalrankjetsradiationinput[IdxRadJets]);
										}
							
									//check if all top radiation partons are matched ---> if so, store jet ranks 
									}else if(GenMatchRadiation->getNumberOfUnmatchedPartons(0)== 0){
										MatchedTopRadiation++;
										
										for(unsigned int i = 0; i < NumberRadiation; i++){			
											int IdxRadJets = GenMatchRadiation->getMatchForParton(i,0);
											matchRadiationjetidx.push_back(IdxRadJets);
											originalrankjetsradiation.push_back(originalrankjetsradiationinput[IdxRadJets]);
										}
									}
								}
								
							}else{
								NoTopRadiation++;
							}
			
						//check if all ISR partons are matched ---> if so, store jet ranks, look at radiation (still jets left)
						}else if(GenMatchISR->getNumberOfUnmatchedPartons(0)== 0){ 
							MatchedISR++;

							for(unsigned int i = 0; i < NumberISR; i++){			
								int IdxISRjets = GenMatchISR->getMatchForParton(i,0);
								matchISRjetidx.push_back(IdxISRjets);
								originalrankjetsISR.push_back(originalrankjetsISRinput[IdxISRjets]);
							}
				
							///////////////////////////////////////////////////////////////////////////
							////////// LOOK AT TOPRADIATION  - JETS LEFT - ISR ALL MATCHED - //////////
							///////////////////////////////////////////////////////////////////////////
							if(NumberRadiation !=0){ 
//	cout << "lala8"<< endl;
								//make new list of jets to be used as input for matching (remove ISR-matched jets from jetlist used as input for ISR)
								for(unsigned int i = 0; i<myselectedjets.size() ; i++){
									for(unsigned int l =0; l<quarks.size(); l++){
										for(unsigned int k =0; k<NumberISR; k++){
											if(i != (unsigned int) originalrankjetsISR[k] && i != (unsigned int) matchquarkjetidx[l]){
												myselectedjetsradiation.push_back(myselectedjets[i]);
												originalrankjetsradiationinput.push_back(i);
											}
										}
									}
								}
								
				  			///////////////////// #top radiation > myselectedjetsradiation.size() /////////////////////
								if(NumberRadiation > myselectedjetsradiation.size()){
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
									JetPartonMatching *GenMatchRadiation = new JetPartonMatching(NfirstpTradiation, myselectedjetsradiation, matchingAlgo_, useMaxDist_, useDeltaR_, maxDist_);

									//check if all partons are matched ---> if so, just store jet ranks, nothing else (no jets left)
									if(GenMatchRadiation->getNumberOfUnmatchedPartons(0)== 0){
										MatchedTopRadiation++;				
							
										//store the rank of the jet for later use (if parton is unmatched ---> rank is negative)
 										for(unsigned int i = 0; i < NfirstpTradiation.size(); i++){			
											int IdxRadJets = GenMatchRadiation->getMatchForParton(i,0);
											matchRadiationjetidx.push_back(IdxRadJets);
											originalrankjetsradiation.push_back(originalrankjetsradiationinput[IdxRadJets]);
										}
							
									//check if unmatched partons ---> if so, categorize, store jet ranks, nothing else (no jets left)
									}else if(GenMatchRadiation->getNumberOfUnmatchedPartons(0)!= 0){ 
										UnmatchedTopRadiation++;
						
										//find parton which is not matched -> check eta/pT/DR ---> count different possibilities
										for(unsigned int i = 0; i<NfirstpTradiation.size(); i++){
											int Idx = GenMatchRadiation->getMatchForParton(i,0);
											if(Idx<0){
												smallestDRRadi = 1000;
												for(unsigned int j = 0; j<NfirstpTradiation.size(); i++){
													if(j!=i){
														DRRadi[j] = deltaR(NfirstpTradiation[i]->p4(),NfirstpTradiation[j]->p4());
	  												if(DRRadi[j] < smallestDRRadi) smallestDRRadi = DRRadi[j];
													}
												}
												if(fabs(NfirstpTradiation[i]->eta())>jetEta_){
													UnmatchedRadCategoryEta++;
												}else if(NfirstpTradiation[i]->pt()<jetPt_){
													UnmatchedRadCategoryPt++;
												}else if(smallestDRRadi<minDRunmatchedParton_){
													UnmatchedRadCategoryDR++;
												}else{
													UnmatchedRadNotCategorized++;
												}
											}
										}
							
 										//store the rank of the jet for later use (if parton is unmatched ---> rank is negative)
 										for(unsigned int i = 0; i < NfirstpTradiation.size(); i++){			
											int IdxRadJets = GenMatchRadiation->getMatchForParton(i,0);
											matchRadiationjetidx.push_back(IdxRadJets);
											originalrankjetsradiation.push_back(originalrankjetsradiationinput[IdxRadJets]);
										}
									}

				  			///////////////////// #top radiation <= myselectedjetsradiation.size() /////////////////////
          			}else if( NumberRadiation <= myselectedjetsradiation.size()){
									JetPartonMatching *GenMatchRadiation = new JetPartonMatching(TopRadiation, myselectedjetsradiation, matchingAlgo_, useMaxDist_, useDeltaR_, maxDist_);
									
									//check if unmatched top radiation  ---> if so, categorize, store jet ranks 
									if(GenMatchRadiation->getNumberOfUnmatchedPartons(0)!= 0){
										UnmatchedTopRadiation++;

										//find parton which is not matched -> check eta/pT/DR ---> count different possibilities
										for(unsigned int i = 0; i<NumberRadiation; i++){
											int Idx = GenMatchRadiation->getMatchForParton(i,0);
											if(Idx<0){
												smallestDRRadi = 0;
												for(unsigned int j = 0; j<NumberRadiation; i++){
													if(j!=i){
														DRRadi[j] = deltaR(TopRadiation[i]->p4(),TopRadiation[j]->p4());
	  												if(DRRadi[j] < smallestDRRadi) smallestDRRadi = DRRadi[j];
													}
												}
												if(fabs(TopRadiation[i]->eta())>jetEta_){
													UnmatchedRadCategoryEta++;
												}else if(TopRadiation[i]->pt()<jetPt_){
													UnmatchedRadCategoryPt++;
												}else if(smallestDRRadi<minDRunmatchedParton_){
													UnmatchedRadCategoryDR++;
												}else{
													UnmatchedRadNotCategorized++;
												}
											}
										}
 							
										//store the rank of the jet for later use (if parton is unmatched ---> rank is negative)
 										for(unsigned int i = 0; i < NumberRadiation; i++){			
											int IdxRadJets = GenMatchRadiation->getMatchForParton(i,0);
											matchRadiationjetidx.push_back(IdxRadJets);
											originalrankjetsradiation.push_back(originalrankjetsradiationinput[IdxRadJets]);
										}
							
									//check if all top radiation partons are matched ---> if so, store jet ranks 
									}else if(GenMatchRadiation->getNumberOfUnmatchedPartons(0)== 0){
										MatchedTopRadiation++;
										
										for(unsigned int i = 0; i < NumberRadiation; i++){			
											int IdxRadJets = GenMatchRadiation->getMatchForParton(i,0);
											matchRadiationjetidx.push_back(IdxRadJets);
											originalrankjetsradiation.push_back(originalrankjetsradiationinput[IdxRadJets]);
										}
									}
								}
								
							}else{
								NoTopRadiation++;
							}
						}					
					}					
				
				///////////////////////////////////////////////////
				////////// LOOK AT TOPRADIATION - NO ISR //////////
				///////////////////////////////////////////////////	
				//no ISR partons (either not a MadGraph sample or really no ISR happened)	
				}else{ 
//	cout << "lala9"<< endl;
					NoISR++;
					
					if(NumberRadiation !=0){ 
						//make new list of jets to be used as input for matching (remove ISR-matched jets from jetlist used as input for ISR)
						for(unsigned int i = 0; i<myselectedjets.size() ; i++){
							for(unsigned int l =0; l<quarks.size(); l++){
								if(i != (unsigned int) matchquarkjetidx[l]){
									myselectedjetsradiation.push_back(myselectedjets[i]);
									originalrankjetsradiationinput.push_back(i);
								}
							}
						}
								
						///////////////////// #top radiation > myselectedjetsradiation.size() /////////////////////
						if(NumberRadiation > myselectedjetsradiation.size()){
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
							JetPartonMatching *GenMatchRadiation = new JetPartonMatching(NfirstpTradiation, myselectedjetsradiation, matchingAlgo_, useMaxDist_, useDeltaR_, maxDist_);

							//check if all partons are matched ---> if so, just store jet ranks, nothing else (no jets left)
							if(GenMatchRadiation->getNumberOfUnmatchedPartons(0)== 0){
								MatchedTopRadiation++;				
						
								//store the rank of the jet for later use (if parton is unmatched ---> rank is negative)
 								for(unsigned int i = 0; i < NfirstpTradiation.size(); i++){			
									int IdxRadJets = GenMatchRadiation->getMatchForParton(i,0);
									matchRadiationjetidx.push_back(IdxRadJets);
									originalrankjetsradiation.push_back(originalrankjetsradiationinput[IdxRadJets]);
								}
							
							//check if unmatched partons ---> if so, categorize, store jet ranks, nothing else (no jets left)
							}else if(GenMatchRadiation->getNumberOfUnmatchedPartons(0)!= 0){ 
								UnmatchedTopRadiation++;
						
								//find parton which is not matched -> check eta/pT/DR ---> count different possibilities
								for(unsigned int i = 0; i<NfirstpTradiation.size(); i++){
									int Idx = GenMatchRadiation->getMatchForParton(i,0);
									if(Idx<0){
										smallestDRRadi = 1000;
										for(unsigned int j = 0; j<NfirstpTradiation.size(); i++){
											if(j!=i){
												DRRadi[j] = deltaR(NfirstpTradiation[i]->p4(),NfirstpTradiation[j]->p4());
	  										if(DRRadi[j] < smallestDRRadi) smallestDRRadi = DRRadi[j];
											}
										}
										if(fabs(NfirstpTradiation[i]->eta())>jetEta_){
											UnmatchedRadCategoryEta++;
										}else if(NfirstpTradiation[i]->pt()<jetPt_){
											UnmatchedRadCategoryPt++;
										}else if(smallestDRRadi<minDRunmatchedParton_){
											UnmatchedRadCategoryDR++;
										}else{
											UnmatchedRadNotCategorized++;
										}
									}	
								}	
							
 								//store the rank of the jet for later use (if parton is unmatched ---> rank is negative)
 								for(unsigned int i = 0; i < NfirstpTradiation.size(); i++){			
									int IdxRadJets = GenMatchRadiation->getMatchForParton(i,0);
									matchRadiationjetidx.push_back(IdxRadJets);
									originalrankjetsradiation.push_back(originalrankjetsradiationinput[IdxRadJets]);
								}
							}	

				  	///////////////////// #top radiation <= myselectedjetsradiation.size() /////////////////////
          	}else if( NumberRadiation <= myselectedjetsradiation.size()){
							JetPartonMatching *GenMatchRadiation = new JetPartonMatching(TopRadiation, myselectedjetsradiation, matchingAlgo_, useMaxDist_, useDeltaR_, maxDist_);
							
							//check if unmatched top radiation  ---> if so, categorize, store jet ranks 
							if(GenMatchRadiation->getNumberOfUnmatchedPartons(0)!= 0){
								UnmatchedTopRadiation++;

								//find parton which is not matched -> check eta/pT/DR ---> count different possibilities
								for(unsigned int i = 0; i<NumberRadiation; i++){
									int Idx = GenMatchRadiation->getMatchForParton(i,0);
									if(Idx<0){
										smallestDRRadi = 0;
										for(unsigned int j = 0; j<NumberRadiation; i++){
											if(j!=i){
												DRRadi[j] = deltaR(TopRadiation[i]->p4(),TopRadiation[j]->p4());
	  										if(DRRadi[j] < smallestDRRadi) smallestDRRadi = DRRadi[j];
											}
										}
										if(fabs(TopRadiation[i]->eta())>jetEta_){
											UnmatchedRadCategoryEta++;
										}else if(TopRadiation[i]->pt()<jetPt_){
											UnmatchedRadCategoryPt++;
										}else if(smallestDRRadi<minDRunmatchedParton_){
											UnmatchedRadCategoryDR++;
										}else{
											UnmatchedRadNotCategorized++;
										}
									}
								}
 							
								//store the rank of the jet for later use (if parton is unmatched ---> rank is negative)
 								for(unsigned int i = 0; i < NumberRadiation; i++){			
									int IdxRadJets = GenMatchRadiation->getMatchForParton(i,0);
									matchRadiationjetidx.push_back(IdxRadJets);
									originalrankjetsradiation.push_back(originalrankjetsradiationinput[IdxRadJets]);
								}
						
							//check if all top radiation partons are matched ---> if so, store jet ranks
							}else if(GenMatchRadiation->getNumberOfUnmatchedPartons(0)== 0){
								MatchedTopRadiation++;
								
								for(unsigned int i = 0; i < NumberRadiation; i++){			
									int IdxRadJets = GenMatchRadiation->getMatchForParton(i,0);
									matchRadiationjetidx.push_back(IdxRadJets);
									originalrankjetsradiation.push_back(originalrankjetsradiationinput[IdxRadJets]);
								}
							}
						}
					}else{
						NoTopRadiation++;
					}
				}


				////////////////////////////////////////////////////////////////
				////////// MAKE PLOTS WHEN ALL QUARKS ARE NOT MATCHED //////////
				////////////////////////////////////////////////////////////////	
//	cout << "lala10"<< endl;
			
				bool Radiation = false;
				for(unsigned int j=0; j<myselectedjets.size();j++){ 
					//make vector with 4 highest pT jets+rank
					if(j<4) {
						FourHighestPtJets.push_back(myselectedjets[j]);
					}
//	cout << myselectedjets.size()<< endl;

//	cout << NumberISR<< endl;
//	cout << NumberRadiation<< endl;

					//keep jets + rank which are matched with radiation!!!
					for(unsigned int l =0; l<NumberRadiation; l++){
						for(unsigned int k =0; k<NumberISR; k++){
							if(NumberISR !=0 && NumberRadiation !=0){
//	cout << "lala10a"<< endl;
								Radiation = true;
//	cout << "lala10b " << originalrankjetsISR[k] <<  " " << originalrankjetsradiation[l]<< endl;
								if(j == (unsigned int) originalrankjetsISR[k] || j == (unsigned int) originalrankjetsradiation[l]){
//	cout << "lala10bbis"<< endl;
									pair<int,pat::Jet> pairjetsrad(j,myselectedjets[j]);
//	cout << "lala10c"<< endl;
									patJetsMatchedWithRadiation.push_back(pairjetsrad);
//	cout << "lala10d"<< endl;
									rankjetsmatchedtorad.push_back(j);
//	cout << "lala10e"<< endl;
									ptjetsmatchedtorad.push_back(myselectedjets[j].pt());
//	cout << "lala10f"<< endl;
								}
//	cout << "lala10g"<< endl;
							}else if(NumberISR ==0 && NumberRadiation !=0){
								Radiation = true;
								if(j == (unsigned int) originalrankjetsradiation[l]){
									pair<int,pat::Jet> pairjetsrad(j,myselectedjets[j]);
									patJetsMatchedWithRadiation.push_back(pairjetsrad);
									rankjetsmatchedtorad.push_back(j);
									ptjetsmatchedtorad.push_back(myselectedjets[j].pt());
								}
							}else if(NumberISR !=0 && NumberRadiation ==0){
								Radiation = true;
								if(j == (unsigned int) originalrankjetsISR[k] ){
									pair<int,pat::Jet> pairjetsrad(j,myselectedjets[j]);
									patJetsMatchedWithRadiation.push_back(pairjetsrad);
									rankjetsmatchedtorad.push_back(j);
									ptjetsmatchedtorad.push_back(myselectedjets[j].pt());
								}
							}
						}
					}

				}
			
				if(Radiation){
//	cout << "lala11"<< endl;
					//highest pT jet has per definition lowest rank!
					//sort first on rank --> make plot of lowest
					//sort on pT -->make plot of highest
					std::sort(rankjetsmatchedtorad.begin(), rankjetsmatchedtorad.end(), std::greater<int>() );
   				std::sort(ptjetsmatchedtorad.begin(), ptjetsmatchedtorad.end(), std::greater<double>() );
					HighestPtRadMatchedJet_Rank->Fill(rankjetsmatchedtorad.back());
   				HighestPtRadMatchedJet_Pt  ->Fill(ptjetsmatchedtorad.front());

					int count = 0;
					for(unsigned int i = 0; i<FourHighestPtJets.size();i++){		
						for(unsigned int j = 0; j<patJetsMatchedWithRadiation.size();i++){
							if(i == (unsigned int) patJetsMatchedWithRadiation[j].first){
								count++;
							}
						}
					}
					if(count = 0)Highest4NoRad++;
					if(count >= 1)Highest4AtLeast1Rad++;
					if(count >= 2)Highest4AtLeast2Rad++;
					if(count >= 3)Highest4AtLeast3Rad++;
					if(count >= 4)Highest4AtLeast4Rad++;

				}else{
					edm::LogInfo("NoDataFound") << "no radiation, so no related plots are filled \n";
				}
			}

		}
		}else{
			genObjectsNULL++;
		}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////// END INCLUDE RADIATION ////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//	cout << "lala12"<< endl;

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
	genObjectsNULL = 0;
	FourQuarksNotMatched = 0;
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
			hExpLightJECPtEtaBin[etab][ptb] = fs->make<TH1F>(obsNamelight,obsNamelight,500,-5.,5.);
			hExpBJECPtEtaBin[etab][ptb] 		= fs->make<TH1F>(obsNameb,obsNameb,500,-5.,5.);
      fExpLightJECPtEtaBin[etab][ptb] = fs->make<TF1>("F_"+obsNamelight,"gaus");
      fExpBJECPtEtaBin[etab][ptb] 		= fs->make<TF1>("F_"+obsNameb,"gaus");
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

  delete [] etabins; 
  delete [] ptbins; 


}

// ------------ method called once each job just after ending the event loop  ------------
void 
TruthReco::endJob() {
	
	edm::LogInfo  ("MainResults") << NrAnalyzed << " number of analyzed events \n";
	edm::LogInfo  ("MainResults") << NbSemiMuGenEvents << " number of semi-mu genevents before any cuts \n";
	edm::LogInfo  ("MainResults") << NbSemiMuGenEvents/NrAnalyzed<< "% semi-mu genevents \n";
	edm::LogInfo  ("MainResults") << NrNot4SelJets/NrAnalyzed << "% events with less then 4 jets passing cuts\n"; 
	edm::LogInfo  ("MainResults") << NrNot1SelIsoMu/NrAnalyzed << "% events with less then 1 muon passing cuts \n"; 
	edm::LogInfo  ("MainResults")	<< NrNot1MET/NrAnalyzed << "% events with less then 1 MET \n"; 
	edm::LogInfo  ("MainResults") << NbSemiMuRecoEvents << " number of semi-mu reco'ed events (>=4jets,>=1muon, >0MET) before any cuts \n";
	edm::LogInfo  ("MainResults") << NbSemiMuRecoEvents/NrAnalyzed<< "% semi-mu reco'ed events \n";
	edm::LogInfo  ("MainResults") << EventNotMatched << " events where the muon is not matched \n \n";

	edm::LogInfo  ("MainResults") << NumberQuarks << " quarks \n";
	edm::LogInfo  ("MainResults") << NumberISR << " partons from ISR \n";
	edm::LogInfo  ("MainResults") << NumberRadiation << " partons from top radiation \n\n";
	
	edm::LogInfo  ("MainResults") << UnmatchedQuarks << " events where not all quarks are matched\n";
	edm::LogInfo  ("MainResults") << "of which " << UnmatchedQuarkCategoryEta << " because of the eta requirement\n";
	edm::LogInfo  ("MainResults") << "of which " << UnmatchedQuarkCategoryPt << " because of the p_{T} requirement\n";
	edm::LogInfo  ("MainResults") << "of which " << UnmatchedQuarkCategoryDR << " because of the DR requirement\n";
	edm::LogInfo  ("MainResults") << "of which " << UnmatchedQuarkNotCategorized << " not categorized\n";
	edm::LogInfo  ("MainResults") << MatchedQuarks << " events where all quarks are matched\n\n";	
	edm::LogInfo  ("MainResults") << UnmatchedISR << " events where not all ISR partons are matched\n";
	edm::LogInfo  ("MainResults") << "of which " << UnmatchedISRCategoryEta << " because of the eta requirement\n";
	edm::LogInfo  ("MainResults") << "of which " << UnmatchedISRCategoryPt << " because of the p_{T} requirement\n";
	edm::LogInfo  ("MainResults") << "of which " << UnmatchedISRCategoryDR << " because of the DR requirement\n";
	edm::LogInfo  ("MainResults") << "of which " << UnmatchedISRNotCategorized << " not categorized\n";	
	edm::LogInfo  ("MainResults") << MatchedISR << " events where all ISR partons are matched\n\n";	
	edm::LogInfo  ("MainResults") << UnmatchedTopRadiation << " events where not all top-radiated partons are matched\n";
	edm::LogInfo  ("MainResults") << "of which " << UnmatchedRadCategoryEta << " because of the eta requirement\n";
	edm::LogInfo  ("MainResults") << "of which " << UnmatchedRadCategoryPt << " because of the p_{T} requirement\n";
	edm::LogInfo  ("MainResults") << "of which " << UnmatchedRadCategoryDR << " because of the DR requirement\n";
	edm::LogInfo  ("MainResults") << "of which " << UnmatchedRadNotCategorized << " not categorized\n";		
	edm::LogInfo  ("MainResults") << MatchedTopRadiation << " events where all top-radiated partons are matched\n";
	edm::LogInfo  ("MainResults") << NoTopRadiation << " events without partons from top radiation\n";
	edm::LogInfo  ("MainResults") << NoISR << " events without partons from ISR\n";
	edm::LogInfo  ("MainResults") << NoRadiation << " events without partons coming from radiation\n \n";

	edm::LogInfo  ("MainResults") << Highest4NoRad << " events in which the 4 highest pT jets are not matched with partons from	radiation\n";
	edm::LogInfo  ("MainResults") << Highest4AtLeast1Rad << " events in which at least 1 of the 4 highest pT jets are matched with partons from	radiation\n";
	edm::LogInfo  ("MainResults") << Highest4AtLeast2Rad << " events in which at least 2 of the 4 highest pT jets are matched with partons from	radiation\n";
	edm::LogInfo  ("MainResults") << Highest4AtLeast3Rad << " events in which at least 3 of the 4 highest pT jets are matched with partons from	radiation\n";
	edm::LogInfo  ("MainResults") << Highest4AtLeast4Rad << " events in which 4 of the 4 highest pT jets are matched with partons from	radiation\n";


	if(NrAnalyzed !=0)edm::LogError ("SummaryError") << NrNotSemiMu/NrAnalyzed << "% events no ttsemimu generated event \n"; 
	edm::LogError ("SummaryError") << genObjectsNULL << " events where at least 1 of the objects in the semi-mu genevent is zero \n";	
	edm::LogError ("SummaryError") << FourQuarksNotMatched << " events where at least 1 of the quarks is not matched \n";

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

}


//define this as a plug-in
DEFINE_FWK_MODULE(TruthReco);
