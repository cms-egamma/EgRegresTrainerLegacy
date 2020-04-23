{
  //binning
  std::vector<double> ptOneBin{15,300};
  std::vector<double> puBins = {0,5,10,15,20,25,30,35,40,45,50,55,60,65,70};
  std::vector<double> eBins= {0,3500};
  std::vector<double> puOneBin = {0,70};
  std::vector<double> resBins = {0.005, 0.0075, 0.01, 0.0125, 0.015, 0.0175, 0.02, 0.0225, 0.025, 0.0275, 0.03, 0.0325, 0.035, 0.0375, 0.04, 0.0425, 0.045, 0.0475, 0.05, 0.055, 0.06, 0.065, 0.07, 0.075, 0.08, 0.09, 0.1,0.12,0.2,0.4,0.5};
  std::vector<double> etaBinsSC = {0.,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.3,1.4,1.4442,1.566,1.7,1.8,1.9,2.,2.25,2.5,2.75,3.0};
  std::vector<double> etaBinsPho = {0.,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.3,1.4,1.4442,1.566,1.7,1.8,1.9,2.,2.25,2.5,2.75,3.0};
  std::vector<double> etBins = {15,30,50,100,150,300};
  std::vector<double> etBinsPho = {10,20,30,50,100,150,300};
  std::vector<double> etBinsSC = {25,40,50,60};
  std::vector<double> etBinsCoarse = {1,20,60,100,120,160,200,220,260,300};
  std::vector<double> etBinSingle = {1,300};
  std::vector<double> etBinsHE = {5,300,1000,1500,3000,4000,5000}; 
  std::vector<double> etBinsHEFine = {5,300,400,500,600,700,800,900,1000,1100,1200,1300,1400,1500,1600,1700,1800,1900,2000,2100,2200,2300,2400,2500,2600,2700,2800,2900,3000,3100,3200,3300,3400,3500,3600,3700,3800,3900,4000,4100,4200,4300,4400,4500,4600,4700,4800,4900,5000}; 
  std::vector<double> etBinsOneBin = {5,5000};

  //mcEnergy Bins
  std::vector<double> mcEBins = {1,500,1000,1500,2000,3500};
  std::vector<double> mcEBinsOneBin = {1,3500};

  //Eta Bins
  std::vector<double> etaBinsEBEE = {0.0,1.442,1.566,2.5};
  std::vector<double> etaBins3 = {-3,-2.75,-2.5,-2.25,-2,-1.8,-1.566,-1.4442,-1.25,-1.1,-0.7,0.,0.7,1.1,1.25,1.4442,1.566,1.8,2.,2.25,2.5,2.75,3.0};
  std::vector<double> etaBins2p5 = {-2.5,-2.25,-2,-1.8,-1.566,-1.4442,-1.25,-1.1,-0.7,0.,0.7,1.1,1.25,1.4442,1.566,1.8,2.,2.25,2.5};
  std::vector<double> absEtaBins2p5 = {0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.3,1.4,1.4442,1.566,1.7,1.8,1.9,2.0,2.25,2.5};
  std::vector<double> etaES = {-2.6,-1.65,1.65,2.6};//seperate preshower and other bins

  std::vector<double> etBinsVs = {5,10,15,20,25,30,35,40,45,50,60,65,70,75,80,85,90,95,100,120,140,160,180,200,220,240,260,280,300};
  std::vector<double> etBinsVsPho = {20,25,30,35,40,45,50,60,65,70,75,80,85,90,95,100,120,140,160,180,200,220,240,260,280,300};
  
  std::vector<double> etBinsMedium = {300,400,500,600,700,800,900,1000};
  std::vector<double> etBinsEleHigh = {1000,1100,1200,1300,1400,1500};
  std::vector<double> etBinsPhoHigh = {1500,1600,1700,1800,1900,2000,2100,2200,2300,2400,2500,2600,2700,2800,2900,3000};
  std::vector<double> etaBinsSimple = {0.,0.7,1.1,1.4442,1.566,2.,2.25,2.5,3.0};
  
  //suppressing noisy fits
  RooMsgService::instance().setGlobalKillBelow(RooFit::FATAL); 
  RooMsgService::instance().setSilentMode(true);
  gErrorIgnoreLevel = kError;
 
  //Directories for making trees
  std::string resultsDirectory = "/home/hep/wrtabb/Egamma/results/";
  std::string inputDirectory = "/home/hep/wrtabb/Egamma/input_trees/2016UL/";

  //2016UL Electrons 
  std::string resultsEle = resultsDirectory + "2016UL/";
  std::string step1InputName = "DoubleElectron_FlatPt-1To300_2016ConditionsFlatPU0to70ECALGT_105X_realistic_IdealEcalIC_v2-v2.root";
  std::string step1Name = "regEleEcal2016UL_IdealIC_IdealTraining_stdVar_stdCuts_ntrees1500_applied.root"; 
  std::string step2Name = "regEleEcal2016UL_RealIC_IdealTraining_stdVar_stdCuts_ntrees1500_applied.root";
  std::string step3Name = "regEleEcal2016UL_RealIC_RealTraining_stdVar_stdCuts_ntrees1500_applied.root";
  std::string step4Name = "regEleEcalTrk2016UL_RealIC_stdVar_stdCuts_ntrees1500_applied.root";

 //The result of training on step 4 only includes egRegTreeFriend
 //egRegTree comes from the output of the prevoius step
 TTree*treeEleStep4 = HistFuncs::makeChain("egRegTree",resultsEle+step3Name,1,1,1);
 TTree*treeEleStep4Friend = HistFuncs::makeChain("egRegTreeFriend",resultsEle+step4Name,1,1,1);
 treeEleStep4->AddFriend(treeEleStep4Friend); 
 
 TTree*treeEleStep3 = HistFuncs::makeChain("egRegTree",resultsEle+step3Name,1,1,1);
 TTree*treeEleStep3Friend = HistFuncs::makeChain("egRegTreeFriend",resultsEle+step3Name,1,1,1);
 treeEleStep3->AddFriend(treeEleStep3Friend); 

 //THe result of training step 1 only includes egRegTreeFriend
 //egRegTree comes from the original sample that training is applied to
 TTree*treeEleStep1 = HistFuncs::makeChain("egRegTree",inputDirectory+step1InputName,1,1,1);
 TTree*treeEleStep1FriendEle = HistFuncs::makeChain("egRegTreeFriend",resultsEle+step1Name,1,1,1);
 treeEleStep1->AddFriend(treeEleStep1FriendEle); 



   /*************************************
   #now as an example do the following, 
   #note the second tree argument is for when I was comparing to a different sample, 
   #ie 102X 2018, now we just set it to null

   ResPlotter resPlotter
   resPlotter.makeHists({regTreeEleReal2018V52018Reg,nullptr},"Real IC, 1.566< |#eta|<2.5","mc.energy>0 && sc.rawEnergy>0 && ssFrac.sigmaIEtaIEta>0 && mc.dR<0.1 && ele.et>0 && eventnr%5>=3","mc.pt","sc.seedEta",etBins,etaBins)
   resPlotter.printFits({3,5,6},"plots/regresFitsThreeComp")

   #or compare two variables 
   resPlotter.printFits({3,6},"plots/regresFitsTwoComp")
  
   ************************************/
}
