
void plot(bool etaPlot, bool puPlot, bool ePlot, bool dcbFit);

void plot2016Ele()
{
 //-----parameters-----//
 //Which plots to plot
 bool etaPlot = true;
 bool puPlot = true;
 bool ePlot = true;
 //Double crystal ball fit option (default is Cruijff)
 bool dcbFit = true;

 gSystem->Exec("gmake RegressionTrainerExe -j 8");
 gSystem->Exec("gmake RegressionApplierExe -j 8");

 gROOT->ProcessLine("gROOT->SetBatch(true)");
 gROOT->ProcessLine(".x rootScripts/setupResPlotter.c");

 plot(etaPlot,puPlot,ePlot,dcbFit);
}

void plot(bool etaPlot, bool puPlot, bool ePlot,bool dcbFit)
{
 gROOT->ProcessLine("ResPlotter res");
 if(dcbFit) gROOT->ProcessLine("res.setFitType(ResFitter::FitType::DCB)");
 //Eta
 if(etaPlot){
  gROOT->ProcessLine("res.makeHists({treeEleStep3,treeEleStep4},\"\",\"mc.energy>0 && sc.rawEnergy>0 && ssFrac.sigmaIEtaIEta>0 && mc.dR<0.1 && ele.et>0 && eventnr%5>=3\",\"mc.pt\",\"sc.seedEta\",etBins,etaBins2p5)");
  gROOT->ProcessLine("res.printFits({1,2,3},\"../plots/2016UL/Eta_\")");
 }

 //Pileup
 if(puPlot){
  gROOT->ProcessLine("res.makeHists({treeEleStep3,treeEleStep4},\"Barrel\",\"abs(sc.seedEta)<1.442 && mc.energy>0 && sc.rawEnergy>0 && ssFrac.sigmaIEtaIEta>0 && mc.dR<0.1 && ele.et>0 && eventnr%5>=3\",\"mc.pt\",\"nrVert\",ptOneBin,puBins)");
  gROOT->ProcessLine("res.printFits({1,2,3},\"../plots/2016UL/PU_EB_\")");

  gROOT->ProcessLine("res.makeHists({treeEleStep3,treeEleStep4},\"Endcap\",\"abs(sc.seedEta)>1.566 && mc.energy>0 && sc.rawEnergy>0 && ssFrac.sigmaIEtaIEta>0 && mc.dR<0.1 && ele.et>0 && eventnr%5>=3\",\"mc.pt\",\"nrVert\",ptOneBin,puBins)");
  gROOT->ProcessLine("res.printFits({1,2,3},\"../plots/2016UL/PU_EE_\")");
 }

 //Et
 if(ePlot){
  gROOT->ProcessLine("res.makeHists({treeEleStep3,treeEleStep4},\"Barrel\",\"abs(sc.seedEta)<1.442 && mc.energy>0 && sc.rawEnergy>0 && ssFrac.sigmaIEtaIEta>0 && mc.dR<0.1 && ele.et>0 && eventnr%5>=3\",\"mc.pt\",\"mc.pt\",ptOneBin,etBins)");
  gROOT->ProcessLine("res.printFits({1,2,3},\"../plots/2016UL/Et_EB_\")");

  gROOT->ProcessLine("res.makeHists({treeEleStep3,treeEleStep4},\"Endcap\",\"abs(sc.seedEta)>1.566 && mc.energy>0 && sc.rawEnergy>0 && ssFrac.sigmaIEtaIEta>0 && mc.dR<0.1 && ele.et>0 && eventnr%5>=3\",\"mc.pt\",\"mc.pt\",ptOneBin,etBins)");
  gROOT->ProcessLine("res.printFits({1,2,3},\"../plots/2016UL/Et_EE_\")");
 }
}
