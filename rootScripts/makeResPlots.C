#include "Utility/DetIdTools.hh"

#include "resStudies.C"


namespace resStudies{
  std::string etStr = "20 < E_{T}^{gen} < 60 GeV";
  std::string eTypeStr = "";
  float fitMin = 0.6;
  float fitMax = 1.2;
}

bool isNextToRingBoundary(int iX,int iY,float eta)
{
  int detId = DetIdTools::makeEcalEndcapId(iX,iY,eta<0 ? -1 : 1); 
  return DetIdTools::isNextToRingBoundary(detId);
}
 

TGraph* makeResVsPlot(TH2* hist2D,float fitMin,float fitMax,bool isMean)
{
  std::vector<float> xValues;
  std::vector<float> xErrs;
  std::vector<float> yValues;
  std::vector<float> yErrs;
  for(int binNr=1;binNr<=hist2D->GetNbinsX();binNr++){
    auto cbFit = getProjWithFit(hist2D,binNr,fitMin,fitMax);
    xValues.push_back(hist2D->GetXaxis()->GetBinLowEdge(binNr)+hist2D->GetXaxis()->GetBinWidth(binNr)/2.);
    xErrs.push_back(hist2D->GetXaxis()->GetBinWidth(binNr)/2.);
    yValues.push_back(isMean ? cbFit.scale : cbFit.sigma);
    yErrs.push_back(isMean ? cbFit.scaleErr : cbFit.sigmaErr);
  }
  auto graph = new TGraphErrors(xValues.size(),xValues.data(),yValues.data(),xErrs.data(),yErrs.data());
  return graph;
}

TGraph* compResVsPlot(TH2* newCorrHist2D,TH2* oldCorrHist2D,TH2* rawHist2D,float fitMin,float fitMax,bool isMean)
{
  auto newCorrGraph = makeResVsPlot(newCorrHist2D,fitMin,fitMax,isMean);
  auto oldCorrGraph = makeResVsPlot(oldCorrHist2D,fitMin,fitMax,isMean);
  auto rawGraph = makeResVsPlot(rawHist2D,fitMin,fitMax,isMean);
  AnaFuncs::setHistAttributes(rawGraph,kAzure+8,2,8,kAzure+8);
  AnaFuncs::setHistAttributes(oldCorrGraph,kOrange+7,2,22,kOrange+7);
  AnaFuncs::setHistAttributes(newCorrGraph,kBlue+2,2,23,kBlue+2);

  newCorrGraph->Draw("AP");
  oldCorrGraph->Draw("P");
  rawGraph->Draw("P");
  auto leg = HistFuncs::makeLegend({{rawGraph,"raw energy"},{oldCorrGraph,"corr energy 74X"},{newCorrGraph,"corr energy 94X"}});
  leg->Draw();
  return newCorrGraph;
}



std::vector<TH2*> makeHists(TTree* tree,int nrBinsX,float xmin,float xmax,int nrBinsY,float ymin,float ymax,const std::vector<std::string>& vars,const std::string& cuts)
{
  std::vector<TH2*> hists;
  for(auto& var : vars){
    hists.push_back(HistFuncs::makeHist(tree,nrBinsX,xmin,xmax,nrBinsY,ymin,ymax,var,cuts));
  }
  return hists;
}

std::vector<TH2*> makeHists(TTree* tree,const std::vector<double>& xbins,int nrBinsY,float ymin,float ymax,const std::vector<std::string>& vars,const std::string& cuts)
{
  std::vector<TH2*> hists;
  for(auto& var : vars){
    hists.push_back(HistFuncs::makeHist(tree,xbins,nrBinsY,ymin,ymax,var,cuts));
  }
  return hists;
}

TH1* compareRes(const std::pair<TH2*,std::string>& newCorrHist2D,
		const std::pair<TH2*,std::string>& oldCorrHist2D,
		const std::pair<TH2*,std::string>& rawHist2D,int binNr)
{
  gStyle->SetOptFit(0);
  auto make1DHist = [](TH2* hist2D,const std::string& baseName,int binNr){ 
    std::ostringstream histName;
    histName <<baseName<<hist2D->GetXaxis()->GetBinLowEdge(binNr)<<"To"<<hist2D->GetXaxis()->GetBinUpEdge(binNr)<<"GeV";
    TH1* hist = (TH1*) hist2D->ProjectionY(histName.str().c_str(),binNr,binNr)->Clone((histName.str()+"Tmp").c_str());
    hist->SetDirectory(0);
    return hist;
  }; 
  
  TH1* newCorrHist = make1DHist(newCorrHist2D.first,"newCorrHist",binNr);
  TH1* oldCorrHist = make1DHist(oldCorrHist2D.first,"oldCorrHist",binNr);
  TH1* rawHist = make1DHist(rawHist2D.first,"rawHist",binNr);
  rawHist->GetXaxis()->SetNdivisions(520);
  rawHist->GetXaxis()->SetRangeUser(0.5,1.4);
  
  float max = std::max(rawHist->GetMaximum(),std::max(newCorrHist->GetMaximum(),oldCorrHist->GetMaximum()));
  rawHist->GetYaxis()->SetRangeUser(0.5,max*3);
    
  
  auto newCorrFitParam = makeResFit(newCorrHist,resStudies::fitMin,resStudies::fitMax);
  auto oldCorrFitParam = makeResFit(oldCorrHist,resStudies::fitMin,resStudies::fitMax);
  auto rawFitParam = makeResFit(rawHist,resStudies::fitMin,resStudies::fitMax);
  newCorrFitParam.legName = newCorrHist2D.second;
  oldCorrFitParam.legName = oldCorrHist2D.second;
  rawFitParam.legName = rawHist2D.second;
  
  plotResHist(rawFitParam,oldCorrFitParam,newCorrFitParam);
    
  auto meanNewLabel = HistFuncs::makeLabel(newCorrFitParam.resLabel(),0.149,0.677,0.402,0.744);
  auto meanOldLabel = HistFuncs::makeLabel(oldCorrFitParam.resLabel(),0.149,0.740,0.402,0.810);
  auto meanRawLabel = HistFuncs::makeLabel(rawFitParam.resLabel(),0.149,0.811,0.402,0.880);
  meanOldLabel->SetTextColor(kOrange+7);
  meanNewLabel->SetTextColor(kBlue+2);
  meanRawLabel->SetTextColor(kAzure+8);
  meanNewLabel->Draw();
  meanOldLabel->Draw();
  meanRawLabel->Draw();
  std::ostringstream binStr;
  binStr<<std::fixed<<std::setprecision(2)<<newCorrHist2D.first->GetXaxis()->GetBinLowEdge(binNr)<<" < #eta < "<<newCorrHist2D.first->GetXaxis()->GetBinLowEdge(binNr+1);

  //  auto labelEta = HistFuncs::makeLabel(binStr.str(),0.685,0.811,0.938,0.88);
  auto labelEta = HistFuncs::makeLabel(binStr.str(),0.149,0.503484,0.402,0.573171);
  labelEta->Draw();

  auto labelEt = HistFuncs::makeLabel(resStudies::etStr,0.149,0.58885,0.402,0.658537);
  labelEt->Draw();

  auto labelEType = HistFuncs::makeLabel(resStudies::eTypeStr,0.149,0.433798,0.402,0.503484);
  labelEType->Draw();

  // rawHist->GetXaxis()->SetNdivisions(510);  
  //rawHist->GetXaxis()->SetRangeUser(0.5,1.4);
  auto c1 = static_cast<TCanvas*>(gROOT->FindObject("c1"));
  c1->Update();
  auto leg = HistFuncs::getFromCanvas<TLegend>(c1,"TLegend")[0];
  //HistFuncs::XYCoord(0.14922,0.5,0.469933,0.674216).setNDC(leg)
  HistFuncs::XYCoord(0.150334,0.254355,0.471047,0.440767).setNDC(leg);
  leg->Draw();

  std::cout <<newCorrHist2D.second<<" "<<newCorrHist->Integral()<<std::endl;
  std::cout <<oldCorrHist2D.second<<" "<<oldCorrHist->Integral()<<std::endl;
  std::cout <<rawHist2D.second<<" "<<rawHist->Integral()<<std::endl;
  

  return rawHist;
}

TH1* compareRes(TH2* newCorrHist2D,TH2* oldCorrHist2D,TH2* rawHist2D,int binNr)
{
  return compareRes({newCorrHist2D,""},{oldCorrHist2D,""},{rawHist2D,""},binNr);
}


void printAllResPlots(TH2* newCorrHist2D,TH2* oldCorrHist2D,TH2* rawHist2D,const std::string& baseName)
{
  auto label = HistFuncs::makeLabel("20 < E^{gen}_{T} < 60 GeV",0.693207,0.759582,0.945991,0.829268);

  for(int binNr=1;binNr<=newCorrHist2D->GetNbinsX();binNr++){
    compareRes(newCorrHist2D,oldCorrHist2D,rawHist2D,binNr);
    std::string histName = baseName+"_binNr"+std::to_string(binNr)+"_eta"+AnaFuncs::convertToTTreeStr(newCorrHist2D->GetXaxis()->GetBinLowEdge(binNr))+"To"+AnaFuncs::convertToTTreeStr(newCorrHist2D->GetXaxis()->GetBinLowEdge(binNr+1));
    label->Draw();
    HistFuncs::print(histName,"c1");
  }
}
