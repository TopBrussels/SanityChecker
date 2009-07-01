#include "../interface/KinematicsChecker.h"

using namespace std;

//
// constructors and destructor
//
KinematicsChecker::KinematicsChecker(const edm::ParameterSet& iConfig, string label)
{
  //now do what ever initialization is needed
  jets_             =   iConfig.getParameter<edm::InputTag>( "labelJets" );
  mets_             =   iConfig.getParameter<edm::InputTag>( "labelMETs" );
  muons_            =   iConfig.getParameter<edm::InputTag>( "labelMuons" );
  electrons_        =   iConfig.getParameter<edm::InputTag>( "labelElectrons" );
  verbose_          =   iConfig.getParameter<bool>( "verbose" );
  NbOfEvents        = 0;
  
  label_            = label;

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
  
  Handle< std::vector<reco::CaloJet> > jets;
  iEvent.getByLabel(jets_, jets);
      
  Handle< std::vector<reco::CaloMET> > mets;
  iEvent.getByLabel(mets_, mets);

  Handle< std::vector<reco::Muon> > muons;
  iEvent.getByLabel(muons_, muons);

  Handle< std::vector<reco::GsfElectron> > electrons;
  iEvent.getByLabel(electrons_, electrons);


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

  //Create TLorentzVector objets from reco objects
  vector< reco::Particle::LorentzVector > ObjectP4s[4];
  for( unsigned int i=0;i<jets->size();i++) {ObjectP4s[0].push_back((*jets)[i].p4());}
  for( unsigned int i=0;i<mets->size();i++) {ObjectP4s[1].push_back((*mets)[i].p4());}
  for( unsigned int i=0;i<muons->size();i++) {ObjectP4s[2].push_back((*muons)[i].p4());}
  for( unsigned int i=0;i<electrons->size();i++) {ObjectP4s[3].push_back((*electrons)[i].p4());}

  //Fill Histograms
  for(unsigned int i=0;i<4;i++){
    TH1Dcontainer_[i]["number"]->Fill(ObjectP4s[i].size()); 
    for(unsigned int j=0; j<ObjectP4s[i].size(); j++) TH1Dcontainer_[i]["pt"]->Fill(ObjectP4s[i][j].Pt());  
    for(unsigned int j=0;j<ObjectP4s[i].size();j++) TH1Dcontainer_[i]["eta"]->Fill(ObjectP4s[i][j].Eta()); 
    for(unsigned int j=0;j<ObjectP4s[i].size();j++) TH1Dcontainer_[i]["et"]->Fill(ObjectP4s[i][j].Et());
    for(unsigned int j=0;j<ObjectP4s[i].size();j++) TH1Dcontainer_[i]["energy"]->Fill(ObjectP4s[i][j].E());
    for(unsigned int j=0;j<ObjectP4s[i].size();j++) TH1Dcontainer_[i]["theta"]->Fill(ObjectP4s[i][j].Theta());
    for(unsigned int j=0;j<ObjectP4s[i].size();j++) TH1Dcontainer_[i]["phi"]->Fill(ObjectP4s[i][j].Phi());
  }
   
   
  for(unsigned int i=0;i<4;i++){
    ObjectP4s[i].clear();
  }
   

}


// ------------ method called once each job just before starting event loop  ------------
void 
KinematicsChecker::beginJob(const edm::EventSetup&)
{
  edm::Service<TFileService> fs;
  if (!fs) throw edm::Exception(edm::errors::Configuration, "TFileService missing from configuration!");
 
  std::string ObjectNames[4] = {"CaloJets","CaloMETs","Muons","Electrons"};

  std::vector< TFileDirectory > subDirs;
  //for(unsigned int i=0;i<4;i++) subDirs.push_back(fs->mkdir( ObjectNames[i] ));
  for(unsigned int i=0;i<4;i++) subDirs.push_back(fs->mkdir( ObjectNames[i]+"_"+label_ ));

  for(unsigned int i=0;i<4;i++){  
    TH1Dcontainer_[i]["number"] = subDirs[i].make<TH1D>("number" ,"number of objects",50,0,30);
    TH1Dcontainer_[i]["pt"] = subDirs[i].make<TH1D>("pt" ,"pt",50,0,200);
    TH1Dcontainer_[i]["et"] = subDirs[i].make<TH1D>("et" ,"et",50,0,200);
    TH1Dcontainer_[i]["eta"] = subDirs[i].make<TH1D>("eta" ,"eta",50,-6,6);
    TH1Dcontainer_[i]["energy"] = subDirs[i].make<TH1D>("energy" ,"energy",50,0,400);
    TH1Dcontainer_[i]["theta"] = subDirs[i].make<TH1D>("theta" ,"theta",50,0,3.2);
    TH1Dcontainer_[i]["phi"] = subDirs[i].make<TH1D>("phi" ,"phi",50,-3.2,3.2);
  }
 

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

 

}

