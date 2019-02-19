{

 gROOT->ProcessLine(".L makeResPlots.C+");
 TTree* regTestTree = HistFuncs::makeChain("een_analyzer/ClusterTree","/eos/cms/store/group/phys_egamma/EgRegression/SCReg/allBar1.root"); 
 regTestTree->AddFriend("regCorr = een_analyzer/ClusterTree","results/reg_sc_stdVarsNoWidthNoSigma_stdCutsAllEvts_ntrees1500_applied.root");
}

/*now example plans in the root session
auto hists = makeHists(regTestTree,{-3.0,-2.5,-2.,-1.6,-1.566,-1.4442,-1.1,-0.7,0.,0.7,1.1,1.4442,1.566,1.6,2.,2.5},150,0,1.5,{"sc.rawEnergy/mc.energy:sc.seedEta","sc.corrEnergy74X/mc.energy:sc.seedEta","regCorr.mean*sc.rawEnergy/mc.energy:sc.seedEta"},"mc.energy>0 && sc.sigmaIEtaIEta>0 && mc.dR<0.1 && mc.pt>20 && mc.pt<60");
compareRes({hists[0],"raw energy"},{hists[1],"current energy"},{hists[2],"new energy"},6); //6 is the bin number, adjust as you like
*/
