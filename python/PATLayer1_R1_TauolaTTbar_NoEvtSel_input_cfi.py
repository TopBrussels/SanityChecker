import FWCore.ParameterSet.Config as cms

readFiles = cms.untracked.vstring(
    'dcap://maite.iihe.ac.be/pnfs/iihe/cms/store/user/ghammad/CMSSW223/Common/PATLayer1/TauolaTTbar_NoEvtSel/PATLayer1_10.root',
    'dcap://maite.iihe.ac.be/pnfs/iihe/cms/store/user/ghammad/CMSSW223/Common/PATLayer1/TauolaTTbar_NoEvtSel/PATLayer1_11.root',
    'dcap://maite.iihe.ac.be/pnfs/iihe/cms/store/user/ghammad/CMSSW223/Common/PATLayer1/TauolaTTbar_NoEvtSel/PATLayer1_12.root',
    'dcap://maite.iihe.ac.be/pnfs/iihe/cms/store/user/ghammad/CMSSW223/Common/PATLayer1/TauolaTTbar_NoEvtSel/PATLayer1_13.root',
    'dcap://maite.iihe.ac.be/pnfs/iihe/cms/store/user/ghammad/CMSSW223/Common/PATLayer1/TauolaTTbar_NoEvtSel/PATLayer1_14.root',
    'dcap://maite.iihe.ac.be/pnfs/iihe/cms/store/user/ghammad/CMSSW223/Common/PATLayer1/TauolaTTbar_NoEvtSel/PATLayer1_15.root',
    'dcap://maite.iihe.ac.be/pnfs/iihe/cms/store/user/ghammad/CMSSW223/Common/PATLayer1/TauolaTTbar_NoEvtSel/PATLayer1_16.root',
    'dcap://maite.iihe.ac.be/pnfs/iihe/cms/store/user/ghammad/CMSSW223/Common/PATLayer1/TauolaTTbar_NoEvtSel/PATLayer1_17.root',
    'dcap://maite.iihe.ac.be/pnfs/iihe/cms/store/user/ghammad/CMSSW223/Common/PATLayer1/TauolaTTbar_NoEvtSel/PATLayer1_18.root',
    'dcap://maite.iihe.ac.be/pnfs/iihe/cms/store/user/ghammad/CMSSW223/Common/PATLayer1/TauolaTTbar_NoEvtSel/PATLayer1_19.root',
    'dcap://maite.iihe.ac.be/pnfs/iihe/cms/store/user/ghammad/CMSSW223/Common/PATLayer1/TauolaTTbar_NoEvtSel/PATLayer1_1.root',
    'dcap://maite.iihe.ac.be/pnfs/iihe/cms/store/user/ghammad/CMSSW223/Common/PATLayer1/TauolaTTbar_NoEvtSel/PATLayer1_20.root',
    'dcap://maite.iihe.ac.be/pnfs/iihe/cms/store/user/ghammad/CMSSW223/Common/PATLayer1/TauolaTTbar_NoEvtSel/PATLayer1_21.root',
    'dcap://maite.iihe.ac.be/pnfs/iihe/cms/store/user/ghammad/CMSSW223/Common/PATLayer1/TauolaTTbar_NoEvtSel/PATLayer1_22.root',
    'dcap://maite.iihe.ac.be/pnfs/iihe/cms/store/user/ghammad/CMSSW223/Common/PATLayer1/TauolaTTbar_NoEvtSel/PATLayer1_23.root',
    'dcap://maite.iihe.ac.be/pnfs/iihe/cms/store/user/ghammad/CMSSW223/Common/PATLayer1/TauolaTTbar_NoEvtSel/PATLayer1_24.root',
    'dcap://maite.iihe.ac.be/pnfs/iihe/cms/store/user/ghammad/CMSSW223/Common/PATLayer1/TauolaTTbar_NoEvtSel/PATLayer1_25.root',
    'dcap://maite.iihe.ac.be/pnfs/iihe/cms/store/user/ghammad/CMSSW223/Common/PATLayer1/TauolaTTbar_NoEvtSel/PATLayer1_26.root',
    'dcap://maite.iihe.ac.be/pnfs/iihe/cms/store/user/ghammad/CMSSW223/Common/PATLayer1/TauolaTTbar_NoEvtSel/PATLayer1_27.root',
    'dcap://maite.iihe.ac.be/pnfs/iihe/cms/store/user/ghammad/CMSSW223/Common/PATLayer1/TauolaTTbar_NoEvtSel/PATLayer1_28.root',
    'dcap://maite.iihe.ac.be/pnfs/iihe/cms/store/user/ghammad/CMSSW223/Common/PATLayer1/TauolaTTbar_NoEvtSel/PATLayer1_29.root',
    'dcap://maite.iihe.ac.be/pnfs/iihe/cms/store/user/ghammad/CMSSW223/Common/PATLayer1/TauolaTTbar_NoEvtSel/PATLayer1_2.root',
    'dcap://maite.iihe.ac.be/pnfs/iihe/cms/store/user/ghammad/CMSSW223/Common/PATLayer1/TauolaTTbar_NoEvtSel/PATLayer1_30.root',
    'dcap://maite.iihe.ac.be/pnfs/iihe/cms/store/user/ghammad/CMSSW223/Common/PATLayer1/TauolaTTbar_NoEvtSel/PATLayer1_31.root',
    'dcap://maite.iihe.ac.be/pnfs/iihe/cms/store/user/ghammad/CMSSW223/Common/PATLayer1/TauolaTTbar_NoEvtSel/PATLayer1_3.root',
    'dcap://maite.iihe.ac.be/pnfs/iihe/cms/store/user/ghammad/CMSSW223/Common/PATLayer1/TauolaTTbar_NoEvtSel/PATLayer1_4.root',
    'dcap://maite.iihe.ac.be/pnfs/iihe/cms/store/user/ghammad/CMSSW223/Common/PATLayer1/TauolaTTbar_NoEvtSel/PATLayer1_5.root',
    'dcap://maite.iihe.ac.be/pnfs/iihe/cms/store/user/ghammad/CMSSW223/Common/PATLayer1/TauolaTTbar_NoEvtSel/PATLayer1_6.root',
    'dcap://maite.iihe.ac.be/pnfs/iihe/cms/store/user/ghammad/CMSSW223/Common/PATLayer1/TauolaTTbar_NoEvtSel/PATLayer1_7.root',
    'dcap://maite.iihe.ac.be/pnfs/iihe/cms/store/user/ghammad/CMSSW223/Common/PATLayer1/TauolaTTbar_NoEvtSel/PATLayer1_8.root',
    'dcap://maite.iihe.ac.be/pnfs/iihe/cms/store/user/ghammad/CMSSW223/Common/PATLayer1/TauolaTTbar_NoEvtSel/PATLayer1_9.root'
)
secFiles = cms.untracked.vstring() 
source = cms.Source ("PoolSource",duplicateCheckMode = cms.untracked.string('noDuplicateCheck') ,fileNames = readFiles#,
)
