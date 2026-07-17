 {
   //binning
   std::vector<double> puBins = {0,5,10,15,20,25,30,35,40,45,50,55,60,65,70};
   std::vector<double> resBins = {0.005, 0.0075, 0.01, 0.0125, 0.015, 0.0175, 0.02, 0.0225, 0.025, 0.0275, 0.03, 0.0325, 0.035, 0.0375, 0.04, 0.0425, 0.045, 0.0475, 0.05, 0.055, 0.06, 0.065, 0.07, 0.075, 0.08, 0.09, 0.1,0.12,0.2,0.4,0.5};
   std::vector<double> etaBins = {0.,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.3,1.4,1.4442,1.566,1.7,1.8,1.9,2.,2.25,2.5};//,2.75,3.0}
   std::vector<double> etaBinsPho = {0.,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.3,1.4,1.4442,1.566,1.7,1.8,1.9,2.,2.25,2.5,2.75,3.0};
   std::vector<double> etaBinsEB = {0.,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.3,1.4,1.4442};//,2.75,3.0}
   std::vector<double> etaBinsGap = {0,0.25,0.70,1.0,1.4442};//,2.75,3.0}
   std::vector<double> etaBinsScout = {0.,0.4,0.8,1.2,1.4442,1.566,2.,2.5};//,2.75,3.0}
   std::vector<double> etaBinsEBEE = {0.,1.4442,1.566,2.5};
   std::vector<double> etBins = {5,15,30,50,100,150,300};
   std::vector<double> etBinsTest = {25,35,50,100,150,300};
   std::vector<double> etBinsScout = {5,10,15,20,30,50,100};
 
   
   std::vector<double> etBinsPho = {10,20,30,50,100,150,300};
   std::vector<double> etBinsSC = {25,40,50,60};
   
   //std::vector<double> etBins = {10,20,30,50,75,100,125,150,200,250,300};
   std::vector<double> etBinsVs = {5,10,15,20,25,30,35,40,45,50,60,65,70,75,80,85,90,95,100,120,140,160,180,200,220,240,260,280,300};
   std::vector<double> etBinsVsPho = {20,25,30,35,40,45,50,60,65,70,75,80,85,90,95,100,120,140,160,180,200,220,240,260,280,300};
   
   std::vector<double> etBinsMedium = {300,400,500,600,700,800,900,1000};
   std::vector<double> etBinsEleHigh = {1000,1100,1200,1300,1400,1500};
   std::vector<double> etBinsPhoHigh = {1500,1600,1700,1800,1900,2000,2100,2200,2300,2400,2500,2600,2700,2800,2900,3000};
   std::vector<double> etaBinsSimple = {0.,0.7,1.1,1.4442,1.566,2.,2.25,2.5};
   
   //suppressing noisy fits
   RooMsgService::instance().setGlobalKillBelow(RooFit::FATAL); 
   RooMsgService::instance().setSilentMode(true);
   gErrorIgnoreLevel = kError;
  
   //trees
   //with 102X regression applied by default
   //TTree* regTreeScouting = HistFuncs::makeChain("egScoutRegTree","resultsEleV1/regEleEcalScout2023_IdealIC_IdealTraining_stdVar_stdCuts_ntrees1500_applied.root",1,1,1);
   //TTree* regTreeScoutingFriend = HistFuncs::makeChain("egScoutRegTreeFriend","resultsEleV1/regEleEcalTrkScout2023_RealIC_stdVar_stdCuts_ntrees1500_applied.root",1,1,1);

   TTree* regTree = HistFuncs::makeChain("egRegTree","/mercury/data3/harper/workexpfiles/DoubleElectron_FlatPT-1t500_skimmedntup.root",1,1,1);
   TTree* regTreeFriend = HistFuncs::makeChain("egRegTreeFriend","resultsEle2/regEleEcalRetrain2026_stdVar_stdCuts_ntrees1500_applied.root",1,1,1);
   TTree* regTreeFriendOld = HistFuncs::makeChain("egRegTreeFriend","resultsEle1/regEleEcalRetrain2026_stdVar_stdCuts_ntrees1500_applied.root",1,1,1);
   TTree* regTreeFriendMin = HistFuncs::makeChain("egRegTreeFriend","resultsEle2/regEleEcalRetrain2026_Minimimal_stdVar_stdCuts_ntrees1500_applied.root",1,1,1);
   //TTree* regTreeFriendMinEta = HistFuncs::makeChain("egRegTreeFriend","resultsEle2/regEleEcalRetrain2026_StdPosEtaZQ_stdVar_stdCuts_ntrees1500_applied.root",1,1,1);
   //TTree* regTreeFriendMinEta = HistFuncs::makeChain("egRegTreeFriend","resultsEle2/regEleEcalRetrain2026_FullPosEtaZQ_15kTree_stdVar_stdCuts_ntrees15000_applied.root",1,1,1);
   TTree* regTreeFriendGapTrain = HistFuncs::makeChain("egRegTreeFriend","resultsEle2/regEleEcalRetrain2026_FullPosEtaZQ_EtaGapNoPhiGap1Clus_stdVar_stdCuts_ntrees15000_applied.root",1,1,1);
   TTree* regTreeFriendFullEtaZQ = HistFuncs::makeChain("egRegTreeFriend","resultsEle2/regEleEcalRetrain2026_FullPosEtaZQ_stdVar_stdCuts_ntrees1500_applied.root",1,1,1);
   regTree->AddFriend(regTreeFriend);
   regTree->AddFriend(regTreeFriendGapTrain,"fullEtaZQEtaGap");
   regTree->AddFriend(regTreeFriendFullEtaZQ,"fullEtaZQ");
  
   /*************************************
   #now as an example do the following, R
   #note the second tree argument is for when I was comparing to a different sample, 
   #ie 102X 2018, now we just set it to null

   ResPlotter resPlotter
   resPlotter.makeHists({regTreeEleReal2018V52018Reg,nullptr},"Real IC, 1.566< |#eta|<2.5","mc.energy>0 && sc.rawEnergy>0 && ssFrac.sigmaIEtaIEta>0 && mc.dR<0.1 && ele.et>0 && eventnr%5>=3","mc.pt","sc.seedEta",etBins,etaBins)
   resPlotter.printFits({3,5,6},"plots/regresFitsThreeComp")

   #or compare two variables 
   resPlotter.printFits({3,6},"plots/regresFitsTwoComp")
  
   ************************************/
   
  ResPlotter resPlotter;
   //resPlotter.makeHists({regTree,nullptr},"z_{vtx}>0, |#eta|>0, q>0, #eta/#phi gap=1,0, #clus=1","mc.energy>0 && sc.rawEnergy>0 && ssFrac.sigmaIEtaIEta>0 && mc.dR<0.1 && ele.et>0 && eventnr%2==1 && sc.scEta>0 && eleExtra.q>0 && eleExtra.z0>0 && sc.etaGapCode==-1 && sc.phiGapCode==0 && sc.numberOfClusters==1","mc.pt","sc.seedEta",etBinsTest,etaBinsGap)
   resPlotter.makeHists({regTree,nullptr},"z_{vtx}>0, |#eta|>0, q>0","mc.energy>0 && sc.rawEnergy>0 && ssFrac.sigmaIEtaIEta>0 && mc.dR<0.1 && ele.et>0 && eventnr%2==1 && sc.scEta>0 && eleExtra.q>0 && eleExtra.z0>0","mc.pt","sc.seedEta",etBinsTest,etaBinsGap)
}
