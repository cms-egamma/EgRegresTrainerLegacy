{

 gROOT->ProcessLine(".L makeResPlots.C+");
 TTree* regTestTree = HistFuncs::makeChain("een_analyzer/ClusterTree","/eos/cms/store/group/phys_egamma/EgRegression/SCReg/allBar1.root"); 
 regTestTree->AddFriend("regCorr = een_analyzer/ClusterTree","results/reg_sc_stdVarsNoWidthNoSigma_stdCutsAllEvts_ntrees1500_applied.root");
}
