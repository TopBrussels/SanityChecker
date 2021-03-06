import FWCore.ParameterSet.Config as cms

readFiles = cms.untracked.vstring(
    'dcap://maite.iihe.ac.be/pnfs/iihe/cms/store/user/pvmulder/CMSSW223/Common/PATLayer1/QCD250to500-madgraph/PATLayer1_601.root',
'dcap://maite.iihe.ac.be/pnfs/iihe/cms/store/user/pvmulder/CMSSW223/Common/PATLayer1/QCD250to500-madgraph/PATLayer1_602.root',
'dcap://maite.iihe.ac.be/pnfs/iihe/cms/store/user/pvmulder/CMSSW223/Common/PATLayer1/QCD250to500-madgraph/PATLayer1_603.root',
'dcap://maite.iihe.ac.be/pnfs/iihe/cms/store/user/pvmulder/CMSSW223/Common/PATLayer1/QCD250to500-madgraph/PATLayer1_604.root',
'dcap://maite.iihe.ac.be/pnfs/iihe/cms/store/user/pvmulder/CMSSW223/Common/PATLayer1/QCD250to500-madgraph/PATLayer1_605.root',
'dcap://maite.iihe.ac.be/pnfs/iihe/cms/store/user/pvmulder/CMSSW223/Common/PATLayer1/QCD250to500-madgraph/PATLayer1_606.root',
'dcap://maite.iihe.ac.be/pnfs/iihe/cms/store/user/pvmulder/CMSSW223/Common/PATLayer1/QCD250to500-madgraph/PATLayer1_607.root',
'dcap://maite.iihe.ac.be/pnfs/iihe/cms/store/user/pvmulder/CMSSW223/Common/PATLayer1/QCD250to500-madgraph/PATLayer1_610.root',
'dcap://maite.iihe.ac.be/pnfs/iihe/cms/store/user/pvmulder/CMSSW223/Common/PATLayer1/QCD250to500-madgraph/PATLayer1_611.root',
'dcap://maite.iihe.ac.be/pnfs/iihe/cms/store/user/pvmulder/CMSSW223/Common/PATLayer1/QCD250to500-madgraph/PATLayer1_612.root',
'dcap://maite.iihe.ac.be/pnfs/iihe/cms/store/user/pvmulder/CMSSW223/Common/PATLayer1/QCD250to500-madgraph/PATLayer1_613.root',
'dcap://maite.iihe.ac.be/pnfs/iihe/cms/store/user/pvmulder/CMSSW223/Common/PATLayer1/QCD250to500-madgraph/PATLayer1_614.root',
'dcap://maite.iihe.ac.be/pnfs/iihe/cms/store/user/pvmulder/CMSSW223/Common/PATLayer1/QCD250to500-madgraph/PATLayer1_615.root',
'dcap://maite.iihe.ac.be/pnfs/iihe/cms/store/user/pvmulder/CMSSW223/Common/PATLayer1/QCD250to500-madgraph/PATLayer1_617.root',
'dcap://maite.iihe.ac.be/pnfs/iihe/cms/store/user/pvmulder/CMSSW223/Common/PATLayer1/QCD250to500-madgraph/PATLayer1_618.root',
'dcap://maite.iihe.ac.be/pnfs/iihe/cms/store/user/pvmulder/CMSSW223/Common/PATLayer1/QCD250to500-madgraph/PATLayer1_619.root',
'dcap://maite.iihe.ac.be/pnfs/iihe/cms/store/user/pvmulder/CMSSW223/Common/PATLayer1/QCD250to500-madgraph/PATLayer1_620.root',
'dcap://maite.iihe.ac.be/pnfs/iihe/cms/store/user/pvmulder/CMSSW223/Common/PATLayer1/QCD250to500-madgraph/PATLayer1_621.root',
'dcap://maite.iihe.ac.be/pnfs/iihe/cms/store/user/pvmulder/CMSSW223/Common/PATLayer1/QCD250to500-madgraph/PATLayer1_622.root',
'dcap://maite.iihe.ac.be/pnfs/iihe/cms/store/user/pvmulder/CMSSW223/Common/PATLayer1/QCD250to500-madgraph/PATLayer1_625.root',
'dcap://maite.iihe.ac.be/pnfs/iihe/cms/store/user/pvmulder/CMSSW223/Common/PATLayer1/QCD250to500-madgraph/PATLayer1_626.root',
'dcap://maite.iihe.ac.be/pnfs/iihe/cms/store/user/pvmulder/CMSSW223/Common/PATLayer1/QCD250to500-madgraph/PATLayer1_627.root',
'dcap://maite.iihe.ac.be/pnfs/iihe/cms/store/user/pvmulder/CMSSW223/Common/PATLayer1/QCD250to500-madgraph/PATLayer1_628.root',
'dcap://maite.iihe.ac.be/pnfs/iihe/cms/store/user/pvmulder/CMSSW223/Common/PATLayer1/QCD250to500-madgraph/PATLayer1_629.root',
'dcap://maite.iihe.ac.be/pnfs/iihe/cms/store/user/pvmulder/CMSSW223/Common/PATLayer1/QCD250to500-madgraph/PATLayer1_630.root',
'dcap://maite.iihe.ac.be/pnfs/iihe/cms/store/user/pvmulder/CMSSW223/Common/PATLayer1/QCD250to500-madgraph/PATLayer1_631.root',
'dcap://maite.iihe.ac.be/pnfs/iihe/cms/store/user/pvmulder/CMSSW223/Common/PATLayer1/QCD250to500-madgraph/PATLayer1_632.root',
'dcap://maite.iihe.ac.be/pnfs/iihe/cms/store/user/pvmulder/CMSSW223/Common/PATLayer1/QCD250to500-madgraph/PATLayer1_633.root',
'dcap://maite.iihe.ac.be/pnfs/iihe/cms/store/user/pvmulder/CMSSW223/Common/PATLayer1/QCD250to500-madgraph/PATLayer1_636.root',
'dcap://maite.iihe.ac.be/pnfs/iihe/cms/store/user/pvmulder/CMSSW223/Common/PATLayer1/QCD250to500-madgraph/PATLayer1_637.root',
'dcap://maite.iihe.ac.be/pnfs/iihe/cms/store/user/pvmulder/CMSSW223/Common/PATLayer1/QCD250to500-madgraph/PATLayer1_638.root',
'dcap://maite.iihe.ac.be/pnfs/iihe/cms/store/user/pvmulder/CMSSW223/Common/PATLayer1/QCD250to500-madgraph/PATLayer1_640.root',
'dcap://maite.iihe.ac.be/pnfs/iihe/cms/store/user/pvmulder/CMSSW223/Common/PATLayer1/QCD250to500-madgraph/PATLayer1_641.root',
'dcap://maite.iihe.ac.be/pnfs/iihe/cms/store/user/pvmulder/CMSSW223/Common/PATLayer1/QCD250to500-madgraph/PATLayer1_642.root',
'dcap://maite.iihe.ac.be/pnfs/iihe/cms/store/user/pvmulder/CMSSW223/Common/PATLayer1/QCD250to500-madgraph/PATLayer1_643.root',
'dcap://maite.iihe.ac.be/pnfs/iihe/cms/store/user/pvmulder/CMSSW223/Common/PATLayer1/QCD250to500-madgraph/PATLayer1_644.root',
'dcap://maite.iihe.ac.be/pnfs/iihe/cms/store/user/pvmulder/CMSSW223/Common/PATLayer1/QCD250to500-madgraph/PATLayer1_647.root',
'dcap://maite.iihe.ac.be/pnfs/iihe/cms/store/user/pvmulder/CMSSW223/Common/PATLayer1/QCD250to500-madgraph/PATLayer1_650.root',
'dcap://maite.iihe.ac.be/pnfs/iihe/cms/store/user/pvmulder/CMSSW223/Common/PATLayer1/QCD250to500-madgraph/PATLayer1_651.root',
'dcap://maite.iihe.ac.be/pnfs/iihe/cms/store/user/pvmulder/CMSSW223/Common/PATLayer1/QCD250to500-madgraph/PATLayer1_653.root',
'dcap://maite.iihe.ac.be/pnfs/iihe/cms/store/user/pvmulder/CMSSW223/Common/PATLayer1/QCD250to500-madgraph/PATLayer1_654.root',
'dcap://maite.iihe.ac.be/pnfs/iihe/cms/store/user/pvmulder/CMSSW223/Common/PATLayer1/QCD250to500-madgraph/PATLayer1_655.root',
'dcap://maite.iihe.ac.be/pnfs/iihe/cms/store/user/pvmulder/CMSSW223/Common/PATLayer1/QCD250to500-madgraph/PATLayer1_656.root',
'dcap://maite.iihe.ac.be/pnfs/iihe/cms/store/user/pvmulder/CMSSW223/Common/PATLayer1/QCD250to500-madgraph/PATLayer1_657.root',
'dcap://maite.iihe.ac.be/pnfs/iihe/cms/store/user/pvmulder/CMSSW223/Common/PATLayer1/QCD250to500-madgraph/PATLayer1_658.root',
'dcap://maite.iihe.ac.be/pnfs/iihe/cms/store/user/pvmulder/CMSSW223/Common/PATLayer1/QCD250to500-madgraph/PATLayer1_659.root',
'dcap://maite.iihe.ac.be/pnfs/iihe/cms/store/user/pvmulder/CMSSW223/Common/PATLayer1/QCD250to500-madgraph/PATLayer1_660.root',
'dcap://maite.iihe.ac.be/pnfs/iihe/cms/store/user/pvmulder/CMSSW223/Common/PATLayer1/QCD250to500-madgraph/PATLayer1_661.root',
'dcap://maite.iihe.ac.be/pnfs/iihe/cms/store/user/pvmulder/CMSSW223/Common/PATLayer1/QCD250to500-madgraph/PATLayer1_662.root',
'dcap://maite.iihe.ac.be/pnfs/iihe/cms/store/user/pvmulder/CMSSW223/Common/PATLayer1/QCD250to500-madgraph/PATLayer1_663.root',
'dcap://maite.iihe.ac.be/pnfs/iihe/cms/store/user/pvmulder/CMSSW223/Common/PATLayer1/QCD250to500-madgraph/PATLayer1_664.root',
'dcap://maite.iihe.ac.be/pnfs/iihe/cms/store/user/pvmulder/CMSSW223/Common/PATLayer1/QCD250to500-madgraph/PATLayer1_665.root',
'dcap://maite.iihe.ac.be/pnfs/iihe/cms/store/user/pvmulder/CMSSW223/Common/PATLayer1/QCD250to500-madgraph/PATLayer1_666.root',
'dcap://maite.iihe.ac.be/pnfs/iihe/cms/store/user/pvmulder/CMSSW223/Common/PATLayer1/QCD250to500-madgraph/PATLayer1_667.root',
'dcap://maite.iihe.ac.be/pnfs/iihe/cms/store/user/pvmulder/CMSSW223/Common/PATLayer1/QCD250to500-madgraph/PATLayer1_670.root',
'dcap://maite.iihe.ac.be/pnfs/iihe/cms/store/user/pvmulder/CMSSW223/Common/PATLayer1/QCD250to500-madgraph/PATLayer1_672.root',
'dcap://maite.iihe.ac.be/pnfs/iihe/cms/store/user/pvmulder/CMSSW223/Common/PATLayer1/QCD250to500-madgraph/PATLayer1_673.root',
'dcap://maite.iihe.ac.be/pnfs/iihe/cms/store/user/pvmulder/CMSSW223/Common/PATLayer1/QCD250to500-madgraph/PATLayer1_674.root',
'dcap://maite.iihe.ac.be/pnfs/iihe/cms/store/user/pvmulder/CMSSW223/Common/PATLayer1/QCD250to500-madgraph/PATLayer1_675.root',
'dcap://maite.iihe.ac.be/pnfs/iihe/cms/store/user/pvmulder/CMSSW223/Common/PATLayer1/QCD250to500-madgraph/PATLayer1_676.root',
'dcap://maite.iihe.ac.be/pnfs/iihe/cms/store/user/pvmulder/CMSSW223/Common/PATLayer1/QCD250to500-madgraph/PATLayer1_677.root',
'dcap://maite.iihe.ac.be/pnfs/iihe/cms/store/user/pvmulder/CMSSW223/Common/PATLayer1/QCD250to500-madgraph/PATLayer1_678.root',
'dcap://maite.iihe.ac.be/pnfs/iihe/cms/store/user/pvmulder/CMSSW223/Common/PATLayer1/QCD250to500-madgraph/PATLayer1_679.root',
'dcap://maite.iihe.ac.be/pnfs/iihe/cms/store/user/pvmulder/CMSSW223/Common/PATLayer1/QCD250to500-madgraph/PATLayer1_680.root',
'dcap://maite.iihe.ac.be/pnfs/iihe/cms/store/user/pvmulder/CMSSW223/Common/PATLayer1/QCD250to500-madgraph/PATLayer1_681.root',
'dcap://maite.iihe.ac.be/pnfs/iihe/cms/store/user/pvmulder/CMSSW223/Common/PATLayer1/QCD250to500-madgraph/PATLayer1_682.root',
'dcap://maite.iihe.ac.be/pnfs/iihe/cms/store/user/pvmulder/CMSSW223/Common/PATLayer1/QCD250to500-madgraph/PATLayer1_683.root',
'dcap://maite.iihe.ac.be/pnfs/iihe/cms/store/user/pvmulder/CMSSW223/Common/PATLayer1/QCD250to500-madgraph/PATLayer1_684.root',
'dcap://maite.iihe.ac.be/pnfs/iihe/cms/store/user/pvmulder/CMSSW223/Common/PATLayer1/QCD250to500-madgraph/PATLayer1_685.root',
'dcap://maite.iihe.ac.be/pnfs/iihe/cms/store/user/pvmulder/CMSSW223/Common/PATLayer1/QCD250to500-madgraph/PATLayer1_686.root',
'dcap://maite.iihe.ac.be/pnfs/iihe/cms/store/user/pvmulder/CMSSW223/Common/PATLayer1/QCD250to500-madgraph/PATLayer1_687.root',
'dcap://maite.iihe.ac.be/pnfs/iihe/cms/store/user/pvmulder/CMSSW223/Common/PATLayer1/QCD250to500-madgraph/PATLayer1_688.root',
'dcap://maite.iihe.ac.be/pnfs/iihe/cms/store/user/pvmulder/CMSSW223/Common/PATLayer1/QCD250to500-madgraph/PATLayer1_689.root',
'dcap://maite.iihe.ac.be/pnfs/iihe/cms/store/user/pvmulder/CMSSW223/Common/PATLayer1/QCD250to500-madgraph/PATLayer1_690.root',
'dcap://maite.iihe.ac.be/pnfs/iihe/cms/store/user/pvmulder/CMSSW223/Common/PATLayer1/QCD250to500-madgraph/PATLayer1_691.root',
'dcap://maite.iihe.ac.be/pnfs/iihe/cms/store/user/pvmulder/CMSSW223/Common/PATLayer1/QCD250to500-madgraph/PATLayer1_692.root',
'dcap://maite.iihe.ac.be/pnfs/iihe/cms/store/user/pvmulder/CMSSW223/Common/PATLayer1/QCD250to500-madgraph/PATLayer1_693.root',
'dcap://maite.iihe.ac.be/pnfs/iihe/cms/store/user/pvmulder/CMSSW223/Common/PATLayer1/QCD250to500-madgraph/PATLayer1_694.root',
'dcap://maite.iihe.ac.be/pnfs/iihe/cms/store/user/pvmulder/CMSSW223/Common/PATLayer1/QCD250to500-madgraph/PATLayer1_695.root',
'dcap://maite.iihe.ac.be/pnfs/iihe/cms/store/user/pvmulder/CMSSW223/Common/PATLayer1/QCD250to500-madgraph/PATLayer1_696.root',
'dcap://maite.iihe.ac.be/pnfs/iihe/cms/store/user/pvmulder/CMSSW223/Common/PATLayer1/QCD250to500-madgraph/PATLayer1_697.root',
'dcap://maite.iihe.ac.be/pnfs/iihe/cms/store/user/pvmulder/CMSSW223/Common/PATLayer1/QCD250to500-madgraph/PATLayer1_698.root',
'dcap://maite.iihe.ac.be/pnfs/iihe/cms/store/user/pvmulder/CMSSW223/Common/PATLayer1/QCD250to500-madgraph/PATLayer1_699.root',
'dcap://maite.iihe.ac.be/pnfs/iihe/cms/store/user/pvmulder/CMSSW223/Common/PATLayer1/QCD250to500-madgraph/PATLayer1_700.root',
'dcap://maite.iihe.ac.be/pnfs/iihe/cms/store/user/pvmulder/CMSSW223/Common/PATLayer1/QCD250to500-madgraph/PATLayer1_701.root',
'dcap://maite.iihe.ac.be/pnfs/iihe/cms/store/user/pvmulder/CMSSW223/Common/PATLayer1/QCD250to500-madgraph/PATLayer1_702.root',
'dcap://maite.iihe.ac.be/pnfs/iihe/cms/store/user/pvmulder/CMSSW223/Common/PATLayer1/QCD250to500-madgraph/PATLayer1_703.root',
'dcap://maite.iihe.ac.be/pnfs/iihe/cms/store/user/pvmulder/CMSSW223/Common/PATLayer1/QCD250to500-madgraph/PATLayer1_704.root',
'dcap://maite.iihe.ac.be/pnfs/iihe/cms/store/user/pvmulder/CMSSW223/Common/PATLayer1/QCD250to500-madgraph/PATLayer1_705.root',
'dcap://maite.iihe.ac.be/pnfs/iihe/cms/store/user/pvmulder/CMSSW223/Common/PATLayer1/QCD250to500-madgraph/PATLayer1_706.root',
'dcap://maite.iihe.ac.be/pnfs/iihe/cms/store/user/pvmulder/CMSSW223/Common/PATLayer1/QCD250to500-madgraph/PATLayer1_707.root',
'dcap://maite.iihe.ac.be/pnfs/iihe/cms/store/user/pvmulder/CMSSW223/Common/PATLayer1/QCD250to500-madgraph/PATLayer1_708.root',
'dcap://maite.iihe.ac.be/pnfs/iihe/cms/store/user/pvmulder/CMSSW223/Common/PATLayer1/QCD250to500-madgraph/PATLayer1_709.root',
'dcap://maite.iihe.ac.be/pnfs/iihe/cms/store/user/pvmulder/CMSSW223/Common/PATLayer1/QCD250to500-madgraph/PATLayer1_710.root',
'dcap://maite.iihe.ac.be/pnfs/iihe/cms/store/user/pvmulder/CMSSW223/Common/PATLayer1/QCD250to500-madgraph/PATLayer1_711.root',
'dcap://maite.iihe.ac.be/pnfs/iihe/cms/store/user/pvmulder/CMSSW223/Common/PATLayer1/QCD250to500-madgraph/PATLayer1_712.root',
'dcap://maite.iihe.ac.be/pnfs/iihe/cms/store/user/pvmulder/CMSSW223/Common/PATLayer1/QCD250to500-madgraph/PATLayer1_713.root',
'dcap://maite.iihe.ac.be/pnfs/iihe/cms/store/user/pvmulder/CMSSW223/Common/PATLayer1/QCD250to500-madgraph/PATLayer1_717.root',
'dcap://maite.iihe.ac.be/pnfs/iihe/cms/store/user/pvmulder/CMSSW223/Common/PATLayer1/QCD250to500-madgraph/PATLayer1_719.root',
'dcap://maite.iihe.ac.be/pnfs/iihe/cms/store/user/pvmulder/CMSSW223/Common/PATLayer1/QCD250to500-madgraph/PATLayer1_721.root',
'dcap://maite.iihe.ac.be/pnfs/iihe/cms/store/user/pvmulder/CMSSW223/Common/PATLayer1/QCD250to500-madgraph/PATLayer1_722.root',
'dcap://maite.iihe.ac.be/pnfs/iihe/cms/store/user/pvmulder/CMSSW223/Common/PATLayer1/QCD250to500-madgraph/PATLayer1_724.root',
'dcap://maite.iihe.ac.be/pnfs/iihe/cms/store/user/pvmulder/CMSSW223/Common/PATLayer1/QCD250to500-madgraph/PATLayer1_725.root',
'dcap://maite.iihe.ac.be/pnfs/iihe/cms/store/user/pvmulder/CMSSW223/Common/PATLayer1/QCD250to500-madgraph/PATLayer1_726.root',
'dcap://maite.iihe.ac.be/pnfs/iihe/cms/store/user/pvmulder/CMSSW223/Common/PATLayer1/QCD250to500-madgraph/PATLayer1_728.root',
'dcap://maite.iihe.ac.be/pnfs/iihe/cms/store/user/pvmulder/CMSSW223/Common/PATLayer1/QCD250to500-madgraph/PATLayer1_729.root',
'dcap://maite.iihe.ac.be/pnfs/iihe/cms/store/user/pvmulder/CMSSW223/Common/PATLayer1/QCD250to500-madgraph/PATLayer1_731.root',
'dcap://maite.iihe.ac.be/pnfs/iihe/cms/store/user/pvmulder/CMSSW223/Common/PATLayer1/QCD250to500-madgraph/PATLayer1_732.root',
'dcap://maite.iihe.ac.be/pnfs/iihe/cms/store/user/pvmulder/CMSSW223/Common/PATLayer1/QCD250to500-madgraph/PATLayer1_733.root',
'dcap://maite.iihe.ac.be/pnfs/iihe/cms/store/user/pvmulder/CMSSW223/Common/PATLayer1/QCD250to500-madgraph/PATLayer1_734.root',
'dcap://maite.iihe.ac.be/pnfs/iihe/cms/store/user/pvmulder/CMSSW223/Common/PATLayer1/QCD250to500-madgraph/PATLayer1_735.root',
'dcap://maite.iihe.ac.be/pnfs/iihe/cms/store/user/pvmulder/CMSSW223/Common/PATLayer1/QCD250to500-madgraph/PATLayer1_736.root',
'dcap://maite.iihe.ac.be/pnfs/iihe/cms/store/user/pvmulder/CMSSW223/Common/PATLayer1/QCD250to500-madgraph/PATLayer1_737.root',
'dcap://maite.iihe.ac.be/pnfs/iihe/cms/store/user/pvmulder/CMSSW223/Common/PATLayer1/QCD250to500-madgraph/PATLayer1_738.root',
'dcap://maite.iihe.ac.be/pnfs/iihe/cms/store/user/pvmulder/CMSSW223/Common/PATLayer1/QCD250to500-madgraph/PATLayer1_739.root',
'dcap://maite.iihe.ac.be/pnfs/iihe/cms/store/user/pvmulder/CMSSW223/Common/PATLayer1/QCD250to500-madgraph/PATLayer1_740.root',
'dcap://maite.iihe.ac.be/pnfs/iihe/cms/store/user/pvmulder/CMSSW223/Common/PATLayer1/QCD250to500-madgraph/PATLayer1_741.root',
'dcap://maite.iihe.ac.be/pnfs/iihe/cms/store/user/pvmulder/CMSSW223/Common/PATLayer1/QCD250to500-madgraph/PATLayer1_742.root',
'dcap://maite.iihe.ac.be/pnfs/iihe/cms/store/user/pvmulder/CMSSW223/Common/PATLayer1/QCD250to500-madgraph/PATLayer1_743.root',
'dcap://maite.iihe.ac.be/pnfs/iihe/cms/store/user/pvmulder/CMSSW223/Common/PATLayer1/QCD250to500-madgraph/PATLayer1_745.root',
'dcap://maite.iihe.ac.be/pnfs/iihe/cms/store/user/pvmulder/CMSSW223/Common/PATLayer1/QCD250to500-madgraph/PATLayer1_746.root',
'dcap://maite.iihe.ac.be/pnfs/iihe/cms/store/user/pvmulder/CMSSW223/Common/PATLayer1/QCD250to500-madgraph/PATLayer1_747.root',
'dcap://maite.iihe.ac.be/pnfs/iihe/cms/store/user/pvmulder/CMSSW223/Common/PATLayer1/QCD250to500-madgraph/PATLayer1_748.root',
'dcap://maite.iihe.ac.be/pnfs/iihe/cms/store/user/pvmulder/CMSSW223/Common/PATLayer1/QCD250to500-madgraph/PATLayer1_749.root',
'dcap://maite.iihe.ac.be/pnfs/iihe/cms/store/user/pvmulder/CMSSW223/Common/PATLayer1/QCD250to500-madgraph/PATLayer1_750.root',
'dcap://maite.iihe.ac.be/pnfs/iihe/cms/store/user/pvmulder/CMSSW223/Common/PATLayer1/QCD250to500-madgraph/PATLayer1_751.root',
'dcap://maite.iihe.ac.be/pnfs/iihe/cms/store/user/pvmulder/CMSSW223/Common/PATLayer1/QCD250to500-madgraph/PATLayer1_752.root',
'dcap://maite.iihe.ac.be/pnfs/iihe/cms/store/user/pvmulder/CMSSW223/Common/PATLayer1/QCD250to500-madgraph/PATLayer1_754.root',
'dcap://maite.iihe.ac.be/pnfs/iihe/cms/store/user/pvmulder/CMSSW223/Common/PATLayer1/QCD250to500-madgraph/PATLayer1_755.root',
'dcap://maite.iihe.ac.be/pnfs/iihe/cms/store/user/pvmulder/CMSSW223/Common/PATLayer1/QCD250to500-madgraph/PATLayer1_757.root',
'dcap://maite.iihe.ac.be/pnfs/iihe/cms/store/user/pvmulder/CMSSW223/Common/PATLayer1/QCD250to500-madgraph/PATLayer1_759.root',
'dcap://maite.iihe.ac.be/pnfs/iihe/cms/store/user/pvmulder/CMSSW223/Common/PATLayer1/QCD250to500-madgraph/PATLayer1_760.root',
'dcap://maite.iihe.ac.be/pnfs/iihe/cms/store/user/pvmulder/CMSSW223/Common/PATLayer1/QCD250to500-madgraph/PATLayer1_763.root',
'dcap://maite.iihe.ac.be/pnfs/iihe/cms/store/user/pvmulder/CMSSW223/Common/PATLayer1/QCD250to500-madgraph/PATLayer1_764.root',
'dcap://maite.iihe.ac.be/pnfs/iihe/cms/store/user/pvmulder/CMSSW223/Common/PATLayer1/QCD250to500-madgraph/PATLayer1_765.root',
'dcap://maite.iihe.ac.be/pnfs/iihe/cms/store/user/pvmulder/CMSSW223/Common/PATLayer1/QCD250to500-madgraph/PATLayer1_768.root',
'dcap://maite.iihe.ac.be/pnfs/iihe/cms/store/user/pvmulder/CMSSW223/Common/PATLayer1/QCD250to500-madgraph/PATLayer1_769.root',
'dcap://maite.iihe.ac.be/pnfs/iihe/cms/store/user/pvmulder/CMSSW223/Common/PATLayer1/QCD250to500-madgraph/PATLayer1_770.root',
'dcap://maite.iihe.ac.be/pnfs/iihe/cms/store/user/pvmulder/CMSSW223/Common/PATLayer1/QCD250to500-madgraph/PATLayer1_772.root',
'dcap://maite.iihe.ac.be/pnfs/iihe/cms/store/user/pvmulder/CMSSW223/Common/PATLayer1/QCD250to500-madgraph/PATLayer1_773.root',
'dcap://maite.iihe.ac.be/pnfs/iihe/cms/store/user/pvmulder/CMSSW223/Common/PATLayer1/QCD250to500-madgraph/PATLayer1_774.root',
'dcap://maite.iihe.ac.be/pnfs/iihe/cms/store/user/pvmulder/CMSSW223/Common/PATLayer1/QCD250to500-madgraph/PATLayer1_775.root',
'dcap://maite.iihe.ac.be/pnfs/iihe/cms/store/user/pvmulder/CMSSW223/Common/PATLayer1/QCD250to500-madgraph/PATLayer1_776.root',
'dcap://maite.iihe.ac.be/pnfs/iihe/cms/store/user/pvmulder/CMSSW223/Common/PATLayer1/QCD250to500-madgraph/PATLayer1_777.root',
'dcap://maite.iihe.ac.be/pnfs/iihe/cms/store/user/pvmulder/CMSSW223/Common/PATLayer1/QCD250to500-madgraph/PATLayer1_778.root',
'dcap://maite.iihe.ac.be/pnfs/iihe/cms/store/user/pvmulder/CMSSW223/Common/PATLayer1/QCD250to500-madgraph/PATLayer1_780.root',
'dcap://maite.iihe.ac.be/pnfs/iihe/cms/store/user/pvmulder/CMSSW223/Common/PATLayer1/QCD250to500-madgraph/PATLayer1_781.root',
'dcap://maite.iihe.ac.be/pnfs/iihe/cms/store/user/pvmulder/CMSSW223/Common/PATLayer1/QCD250to500-madgraph/PATLayer1_782.root',
'dcap://maite.iihe.ac.be/pnfs/iihe/cms/store/user/pvmulder/CMSSW223/Common/PATLayer1/QCD250to500-madgraph/PATLayer1_783.root',
'dcap://maite.iihe.ac.be/pnfs/iihe/cms/store/user/pvmulder/CMSSW223/Common/PATLayer1/QCD250to500-madgraph/PATLayer1_784.root',
'dcap://maite.iihe.ac.be/pnfs/iihe/cms/store/user/pvmulder/CMSSW223/Common/PATLayer1/QCD250to500-madgraph/PATLayer1_785.root',
'dcap://maite.iihe.ac.be/pnfs/iihe/cms/store/user/pvmulder/CMSSW223/Common/PATLayer1/QCD250to500-madgraph/PATLayer1_786.root',
'dcap://maite.iihe.ac.be/pnfs/iihe/cms/store/user/pvmulder/CMSSW223/Common/PATLayer1/QCD250to500-madgraph/PATLayer1_787.root',
'dcap://maite.iihe.ac.be/pnfs/iihe/cms/store/user/pvmulder/CMSSW223/Common/PATLayer1/QCD250to500-madgraph/PATLayer1_788.root',
'dcap://maite.iihe.ac.be/pnfs/iihe/cms/store/user/pvmulder/CMSSW223/Common/PATLayer1/QCD250to500-madgraph/PATLayer1_789.root',
'dcap://maite.iihe.ac.be/pnfs/iihe/cms/store/user/pvmulder/CMSSW223/Common/PATLayer1/QCD250to500-madgraph/PATLayer1_790.root',
'dcap://maite.iihe.ac.be/pnfs/iihe/cms/store/user/pvmulder/CMSSW223/Common/PATLayer1/QCD250to500-madgraph/PATLayer1_791.root',
'dcap://maite.iihe.ac.be/pnfs/iihe/cms/store/user/pvmulder/CMSSW223/Common/PATLayer1/QCD250to500-madgraph/PATLayer1_792.root',
'dcap://maite.iihe.ac.be/pnfs/iihe/cms/store/user/pvmulder/CMSSW223/Common/PATLayer1/QCD250to500-madgraph/PATLayer1_793.root',
'dcap://maite.iihe.ac.be/pnfs/iihe/cms/store/user/pvmulder/CMSSW223/Common/PATLayer1/QCD250to500-madgraph/PATLayer1_794.root',
'dcap://maite.iihe.ac.be/pnfs/iihe/cms/store/user/pvmulder/CMSSW223/Common/PATLayer1/QCD250to500-madgraph/PATLayer1_795.root',
'dcap://maite.iihe.ac.be/pnfs/iihe/cms/store/user/pvmulder/CMSSW223/Common/PATLayer1/QCD250to500-madgraph/PATLayer1_796.root',
'dcap://maite.iihe.ac.be/pnfs/iihe/cms/store/user/pvmulder/CMSSW223/Common/PATLayer1/QCD250to500-madgraph/PATLayer1_797.root',
'dcap://maite.iihe.ac.be/pnfs/iihe/cms/store/user/pvmulder/CMSSW223/Common/PATLayer1/QCD250to500-madgraph/PATLayer1_798.root'
)
secFiles = cms.untracked.vstring() 
source = cms.Source ("PoolSource",duplicateCheckMode = cms.untracked.string('noDuplicateCheck') ,fileNames = readFiles#,
)
