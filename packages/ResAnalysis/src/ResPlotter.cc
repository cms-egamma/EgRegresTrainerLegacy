#include "ResAnalysis/ResPlotter.hh"

#include "Utility/DetIdTools.hh"
#include "Utility/LogErr.hh"
#include "Utility/HistFuncs.hh"

#include "TRandom3.h"
#include "TH2D.h"
#include "TTree.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TStyle.h"
#include "RooPlot.h"
#include "RooHist.h"
#include "TROOT.h"


void ResPlotter::Config::setDefaults()
{
  nrResBins = 300;
  resMin = 0.;
  resMax = 1.5;
  fitMin = 0.5;
  fitMax = 1.3;
  fitMinHigh = 0.8;
  fitMaxHigh = 1.1;
  fitHighThres = 50;

  normalise = true;

  binLabelPrecision = 3;
  divideMeanBySigma = true;

  std::vector<std::pair<std::string,std::string> > varsTree1 = {
    {"sc.rawEnergy/mc.energy","raw energy"},
    {"sc.corrEnergy/mc.energy","74X correction"},
    {"sc.corrEnergyAlt/mc.energy","2018 UL correction"},
    {"eleAltEnergy1.ecal/mc.energy","80X ecal"},
    {"eleAltEnergy1.ecalTrk/mc.energy","80X ecal-trk"},
    {"phoAltEnergy1.ecal/mc.energy","80X pho"},
    {"ele.ecalEnergy/mc.energy","2018UL ecal"}, 
    {"ele.energy/mc.energy","2018UL ecal-trk"},
    {"pho.energy/mc.energy","2018UL pho"}
  };

  std::vector<std::pair<std::string,std::string> > varsTree2 = {
    {"sc.rawEnergy/mc.energy","raw energy, 102X"},
    {"sc.corrEnergy/mc.energy","74X corr, 102X"},
    {"ele.ecalEnergy/mc.energy","80X ecal, 102X"}, 
    {"ele.energy/mc.energy","80X ecal-trk, 102X"},
    {"pho.energy/mc.energy","80X pho, 102X"}
  };
  vars.clear();
  vars.push_back(varsTree1);
  vars.push_back(varsTree2);
}

void ResPlotter::VarNameData::autoFill()
{
  if(name=="mc.pt"){
    plotname = "E_{T}^{gen}";
    filename = "Et";
    unit = "GeV";
  }else if(name=="sc.seedEta"){
    plotname = "#eta";
    filename = "Eta";
    unit = "";
  }else if(name=="abs(sc.seedEta)"){
    plotname = "|#eta|";
    filename = "AbsEta";
    unit = "";
  }else if(name=="nrPUIntTrue"){
    plotname = "<#PU>";
    filename = "NrTruePU";
    unit = "";
  }else{
    plotname = name;
    filename = name;
    unit = "";
  }
}

void ResPlotter::makeHists(std::vector<TTree*> trees,const std::string& label,const std::string& cuts,
			   const std::string& vsVar1,const std::string& vsVar2,
			   const std::vector<double>& vsVar1Bins,const std::vector<double>& vsVar2Bins)
{
  
  if(trees.size()!=cfg_.vars.size()){
    LogErr<<" error trees size "<<trees.size()<<" does not equal vars size "<<cfg_.vars.size()<<std::endl;
    return;
  }

  histsVec_.clear(); //nasty memory leak here fix it
  histsVec_.resize(vsVar1Bins.size()-1);
  vsVar1Bins_ = vsVar1Bins;
  vsVar2Bins_ = vsVar2Bins;
  vsVar1_ = vsVar1;
  vsVar2_ = vsVar2;
  label_ = label;

  for(size_t treeNr=0;treeNr<trees.size();treeNr++){
    if(trees[treeNr]==nullptr) continue;

    auto treeHists = makeHists(trees[treeNr],cfg_.vars[treeNr],cuts);
    for(size_t vsVar1BinNr=0;vsVar1BinNr<histsVec_.size();vsVar1BinNr++){
      for(auto& treeHist : treeHists[vsVar1BinNr]){
	histsVec_[vsVar1BinNr].emplace_back(std::move(treeHist));
      }
    }
  }

  if(cfg_.normalise) normaliseHists();
  
}

std::vector<std::vector<std::pair<TH2*,std::string> > > 
ResPlotter::makeHists(TTree* tree,const std::vector<std::pair<std::string,std::string> >& vars,
		      const std::string& cuts)const			    
{

  std::vector<std::vector<std::pair<TH2*,std::string> > > outHistsVec(vsVar1Bins_.size()-1);
  for(auto& hists : outHistsVec){
    for(const auto& var : vars){
      TH2* hist = new TH2D("hist",(";"+vsVar2_.plotname).c_str(),vsVar2Bins_.size()-1,&vsVar2Bins_[0],cfg_.nrResBins,cfg_.resMin,cfg_.resMax);
      hist->Sumw2();
      hist->SetDirectory(0);
      hists.push_back({hist,var.second});
    }
    
  }
  std::string varStr = vsVar1_.name+":"+vsVar2_.name;
  for(auto& var : vars){
    varStr+=":"+var.first;
  }
  auto vsVar1BinNr = [this](float vsVar)->size_t{
    for(size_t binNr=0;binNr+1<vsVar1Bins_.size();binNr++){
      if(vsVar>=vsVar1Bins_[binNr] && vsVar<vsVar1Bins_[binNr+1]) return binNr;
    }
    return vsVar1Bins_.size();
  };
  auto treeData = HistFuncs::readTree(tree,varStr,cuts);
  for(auto entry : treeData){
    size_t binNr=vsVar1BinNr(entry[0]);
    if(binNr<outHistsVec.size()){
      for(size_t histNr=0;histNr<outHistsVec[binNr].size();histNr++){
	outHistsVec[binNr][histNr].first->Fill(entry[1],entry[histNr+2]);
      }
    }
  }
  return outHistsVec;
}


void ResPlotter::printFits(const std::vector<int>& histNrs,const std::string& baseOutName)const
{
  bool twoComp = histNrs.size()==2;  
  if(histNrs.size()!=2 && histNrs.size()!=3){
    LogErr << "Error, number of selected histograms must be either 2 or 3, not "<<histNrs.size()<<std::endl;
    return;
  }
  auto vsVar1Label = HistFuncs::makeLabel("",0.657016,0.30662,0.9098,0.374564);
  auto idealLabel = HistFuncs::makeLabel(label_,0.159106,0.798303,0.41364,0.864361);
   
  for(size_t vsVar1BinNr=0;vsVar1BinNr<histsVec_.size();vsVar1BinNr++){

    float fitMin=0;
    float fitMax=0;

    if(vsVar1Bins_[vsVar1BinNr]<cfg_.fitHighThres){
      fitMin = cfg_.fitMin;
      fitMax = cfg_.fitMax;
    }else{
      fitMin = cfg_.fitMinHigh;
      fitMax = cfg_.fitMaxHigh;     
    }

    std::ostringstream vsVar1LabelStr;
    vsVar1LabelStr <<vsVar1Bins_[vsVar1BinNr]<<" < "<<vsVar1_.plotname<<" < "<<vsVar1Bins_[vsVar1BinNr+1];
    if(!vsVar1_.unit.empty()) vsVar1LabelStr<<" "<<vsVar1_.unit<<std::endl;
    vsVar1Label->SetLabel(vsVar1LabelStr.str().c_str());
    
    std::string nameSuffix = vsVar1_.filename+"BinNr"+std::to_string(vsVar1BinNr)+"_"+vsVar1_.filename+""+AnaFuncs::convertToTTreeStr(vsVar1Bins_[vsVar1BinNr])+"To"+AnaFuncs::convertToTTreeStr(vsVar1Bins_[vsVar1BinNr+1]);
    std::string outName = baseOutName+"_"+nameSuffix;
    std::string outNameMean = baseOutName+"Mean_"+nameSuffix;
    std::string outNameSigma = baseOutName+"Sigma_"+nameSuffix;
    
    
    auto hists = histsVec_[vsVar1BinNr];
  
    //now do all the fits
    std::vector<ResFitter::ParamsVsVar> fitParams;
    for(int histNr : histNrs){
      std::pair<TH2*,std::string>& histPair = hists[histNr];
      fitParams.push_back(resFitter_.makeFitVsVar(histPair.first,fitMin,fitMax,histPair.second));
    }
  
    //now we plot the fit params vs the variable of interets
    auto graphSigma = plotFitParamsVsVarComp(fitParams,ResFitter::ValType::Sigma,cfg_.divideMeanBySigma);
    if(twoComp) formatTwoComp(graphSigma,vsVar1Label,idealLabel);
    else formatThreeComp(graphSigma,vsVar1Label,idealLabel);
    if(!baseOutName.empty()) HistFuncs::print(outNameSigma,"c1",true);
    
    auto graphMean = plotFitParamsVsVarComp(fitParams,ResFitter::ValType::Mean);
    if(twoComp) formatTwoComp(graphMean,vsVar1Label,idealLabel,true);
    else formatThreeComp(graphMean,vsVar1Label,idealLabel,true);
    if(!baseOutName.empty()) HistFuncs::print(outNameMean,"c1",true);

    //now we plot the individual fits to check for errors
    if(!baseOutName.empty()){
      if(twoComp){
	//need to reset canvas to normal non-ratio
	auto c1 = static_cast<TCanvas*>(gROOT->FindObject("c1"));
	delete c1;
      }
      printResComps(fitParams,outName,{0.55,1.2},vsVar1LabelStr.str());
    }   
  }  
}

void ResPlotter::printResComps(const std::vector<ResFitter::ParamsVsVar>& fitParamsVsVars,
			       const std::string& baseName,const std::pair<float,float>& plotRange,
			       const std::string& regionStr)const
{
  gStyle->SetOptFit(0);

  for(size_t binNr=0;binNr<fitParamsVsVars[0].params().size();binNr++){
    std::vector<ResFitter::Param> fitParams;
    for(auto& fitParamsVsVar : fitParamsVsVars){
      fitParams.push_back(fitParamsVsVar.params()[binNr]);
    }
      
    plotResComp(fitParams,plotRange);
    double binMin = fitParamsVsVars[0].binLowEdges()[binNr];
    double binMax = fitParamsVsVars[0].binLowEdges()[binNr+1];
    

    std::ostringstream binStr;
    binStr<<std::fixed<<std::setprecision(cfg_.binLabelPrecision)<<binMin<<" #leq "<<vsVar2_.plotname<<" < "<<binMax;
  
    auto binLabel = HistFuncs::makeLabel(binStr.str(),0.149,0.503484,0.402,0.573171);
    binLabel->Draw();
    
    auto regionLabel = HistFuncs::makeLabel(regionStr,0.149,0.58885,0.402,0.658537);
    regionLabel->Draw();
    
    auto mainLabel = HistFuncs::makeLabel(label_,0.149,0.433798,0.402,0.503484);
    mainLabel->Draw();
    
    auto c1 = static_cast<TCanvas*>(gROOT->FindObject("c1"));
    c1->Update();
    auto leg = HistFuncs::getFromCanvas<TLegend>(c1,"TLegend")[0];
    HistFuncs::XYCoord(0.150334,0.254355,0.471047,0.440767).setNDC(leg);
    leg->Draw();

    std::string histName = baseName+"_binNr"+std::to_string(binNr)+"_"+vsVar2_.filename+
      AnaFuncs::convertToTTreeStr(binMin)+"To"+AnaFuncs::convertToTTreeStr(binMax);
    HistFuncs::print(histName,"c1",true);
  }
}

TGraph* ResPlotter::plotFitParamsVsVarComp(const std::vector<ResFitter::ParamsVsVar>& fits,
					   ResFitter::ValType valType,
					   bool divideSigmaByMean)const
{
  if(fits.size()!=2 && fits.size()!=3){
    std::cout <<"Error: comparisions are only valid for 2 or 3 quanities, not "<<fits.size()<<std::endl;
    return nullptr;
  }

  std::vector<TGraph*> graphs;
  for(const auto& fit : fits) graphs.push_back(fit.makeGraph(valType,divideSigmaByMean));

  for(size_t graphNr=0;graphNr<graphs.size();graphNr++) setStyle(graphs[graphNr],graphNr);

  double min = 0.; //actually set to 0 for now
  double max = 0;
  for(int pointNr=0;pointNr<graphs[0]->GetN();pointNr++){
    for(auto graph : graphs){
      min = std::min(graph->GetY()[pointNr],min);
      max = std::max(graph->GetY()[pointNr],max);
    }
  }

  std::vector<TPad*> pads;
  if(fits.size()==2){
    auto c1 = HistFuncs::makeRatioCanvas("c1");
    pads = HistFuncs::getFromCanvas<TPad>(c1,"TPad");
  }

  std::vector<std::pair<TGraph*,std::string> > legEntries;
  for(size_t graphNr=0;graphNr<graphs.size();graphNr++){
    if(graphNr==0){
      graphs[graphNr]->GetYaxis()->SetRangeUser(min,max*1.1);
      graphs[graphNr]->Draw("AP");
    }else{
      graphs[graphNr]->Draw("P");
    }
    legEntries.push_back({graphs[graphNr],fits[graphNr].legName()});
  }
  
  auto leg = HistFuncs::makeLegend(legEntries);
  leg->Draw();

  if(fits.size()==2){
    pads[1]->cd();
    auto ratioGraph = makeRatio(graphs[1],graphs[0]);
    ratioGraph->SetTitle(";;ratio");
    ratioGraph->Draw("AP");
    ratioGraph->GetXaxis()->SetLabelSize(0.1);
    ratioGraph->GetXaxis()->SetTitleSize(0.1);
    ratioGraph->GetYaxis()->SetLabelSize(0.1);
    ratioGraph->GetYaxis()->SetTitleSize(0.11);
    ratioGraph->GetYaxis()->SetTitleOffset(0.60); 
    ratioGraph->GetYaxis()->SetTickLength(0.04);
    ratioGraph->GetXaxis()->SetTickLength(0.06);
    pads[0]->cd();
    return ratioGraph;
  }else{
    return graphs[0];
  }
}

RooPlot* ResPlotter::plotResComp(std::vector<ResFitter::Param>& fitParams,
				 const std::pair<float,float>& xRange)const
{
  double max = 0.;
  std::vector<std::pair<TGraph*,std::string>> legEntries;
  for(size_t fitNr=0;fitNr<fitParams.size();fitNr++){
    auto& fitParam = fitParams[fitNr];
    TGraph* graph = static_cast<TGraph*>(fitParam.plot->getHist("h_res"));
    setStyle(graph,fitNr);
    fitParam.plot->getCurve("model_Norm[res]")->SetLineColor(getColour(fitNr));
    fitParam.plot->remove("model_paramBox");
    fitParam.plot->SetTitle("");
    
    max = std::max(fitParam.plot->GetMaximum(),max);
    if(fitNr==0){
      if(xRange.first!=xRange.second){
	fitParam.plot->GetXaxis()->SetRangeUser(xRange.first,xRange.second);
      }
      fitParam.plot->Draw();

      auto c1 = static_cast<TCanvas*>(gROOT->FindObject("c1"));
      // c1->SetLogy(1);
      c1->SetGridx(1);
      c1->SetGridy(1);

    }else{
      fitParam.plot->Draw("SAME");
    }
    auto meanLabel = HistFuncs::makeLabel(fitParams[fitNr].resLabel(),
					  0.149,0.811-fitNr*0.07,
					  0.402,0.880-fitNr*0.07);
    meanLabel->SetTextColor(getColour(fitNr));
    meanLabel->Draw();
    legEntries.push_back({graph,fitParam.legName});
  }
  fitParams[0].plot->SetMaximum(max*1.05);
  auto leg = HistFuncs::makeLegend(legEntries,0.167038,0.56446,0.488864,0.735192);
  leg->SetFillStyle(0);
  leg->Draw();

  return fitParams[0].plot;
}

void ResPlotter::formatTwoComp(TGraph* graph,TPaveLabel* vsVar1Label,TPaveLabel* infoLabel,bool isMean)const
{ 
  float vsVar2Max = vsVar2Bins_.back();
   
  auto c1 = static_cast<TCanvas*>(gROOT->FindObject("c1"));
  auto pads = HistFuncs::getFromCanvas<TPad>(c1,"TPad");
  for(auto &pad : pads){
    pad->SetGridx();
    pad->SetGridy();
    pad->Update();
  }
  pads[0]->cd();
  auto topGraph = HistFuncs::getFromCanvas<TGraph>(pads[0],"TGraphErrors")[0];
  topGraph->GetXaxis()->SetRangeUser(0.,vsVar2Max);
  if(isMean) topGraph->GetYaxis()->SetRangeUser(0.9,1.05);
  topGraph->SetTitle("");
  auto leg = HistFuncs::getFromCanvas<TLegend>(pads[0],"TLegend")[0];
  leg->SetTextSize(0.0440388);
  //      HistFuncs::XYCoord(0.644766,0.155052,0.997773,0.301394).setNDC(leg);
  HistFuncs::XYCoord(0.57386,0.155137,0.927253,0.301265).setNDC(leg);
  leg->SetFillStyle(0);
  leg->Draw();
  vsVar1Label->Draw();
  infoLabel->Draw();
  
  graph->SetTitle((";"+vsVar2_.axisLabel()).c_str());
  graph->GetXaxis()->SetRangeUser(0,vsVar2Max);
  graph->SetMarkerStyle(8);
  graph->GetYaxis()->SetRangeUser(0.8,1.1);
}

void ResPlotter::formatThreeComp(TGraph* graph,TPaveLabel* vsVar1Label,TPaveLabel* infoLabel,bool isMean)const 
{
  float vsVar2Max = vsVar2Bins_.back();
  
  auto c1 = static_cast<TCanvas*>(gROOT->FindObject("c1"));
  c1->SetGridx();
  c1->SetGridy();
  c1->Update();
  auto leg = HistFuncs::getFromCanvas<TLegend>(c1,"TLegend")[0];
  //HistFuncs::XYCoord(0.644766,0.155052,0.997773,0.301394).setNDC(leg);   
  HistFuncs::XYCoord(0.57386,0.155137,0.927253,0.301265).setNDC(leg);
  leg->SetFillStyle(0);
  leg->Draw();
  
  vsVar1Label->Draw();
  infoLabel->Draw();
  graph->SetTitle((";"+vsVar2_.axisLabel()).c_str());
  graph->GetXaxis()->SetRangeUser(0,vsVar2Max); 
  if(isMean) graph->GetYaxis()->SetRangeUser(0.9,1.05);
}



TGraph* ResPlotter::makeRatio(TGraph* numer,TGraph* denom)
{
  std::vector<float> xValues;
  std::vector<float> xErrs;
  std::vector<float> yValues;
  std::vector<float> yErrs;
  if(numer->GetN()!=denom->GetN()){
    LogErr <<"makeRatio: graphs have different entries "<<numer->GetN()<<" "<<denom->GetN()<<std::endl;
    return nullptr;
  }
  for(int pointNr=0;pointNr<numer->GetN();pointNr++){
    xValues.push_back(denom->GetX()[pointNr]);
    xErrs.push_back(denom->GetEX()[pointNr]);
    float numerVal = numer->GetY()[pointNr];
    float numerErr = numer->GetEY()[pointNr];
    float denomVal = denom->GetY()[pointNr];
    float denomErr = denom->GetEY()[pointNr];
    auto sq = [](float x){return x*x;};
    float val = denomVal !=0 ? numerVal/denomVal : 0;
    float err = 0;
    if(denomVal!=0) err+= sq(denomErr/denomVal);
    if(numerVal!=0) err+= sq(numerErr/numerVal);
    err = std::sqrt(err)*val;
    
    yValues.push_back(val);
    yErrs.push_back(err);
  }
  return new TGraphErrors(xValues.size(),xValues.data(),yValues.data(),xErrs.data(),yErrs.data());
}

int ResPlotter::getColour(unsigned int colourNr)
{
  switch(colourNr%3){
  case 0:
    return kAzure+8;
  case 1:
    return kOrange+7;
  case 2:
    return kBlue+2;
  }
  return 0;
}

int ResPlotter::getMarkerStyle(unsigned int markerNr)
{
  switch(markerNr%3){
  case 0:
    return 8;
  case 1:
    return 22;
  case 2:
    return 23;
  }
  return 0;
}
    
 
//urgh this painful 
//ideally we would do this all at hist creation time but we dont know what to normalise to
//as it depends on the number of events passing selection the tree and we dont have an easy access
//so we just normalise hists after the fact
//note due to our rather odd vector layout (it organically grew) a single histogram is consists 
//of multiple 2D histograms
void ResPlotter::normaliseHists()
{
  size_t nrHists = !histsVec_.empty() ? histsVec_[0].size() : 0;
  float minIntegral = std::numeric_limits<float>::max();
  std::vector<float> histIntegrals;
  for(size_t histNr=0;histNr<nrHists;histNr++){
    float histIntegral = 0;
    //now we sum over all vsVar1 bins
    for(const auto& histVec : histsVec_){
      TH2* hist = histVec[histNr].first;
      histIntegral+=hist->Integral(0,hist->GetNbinsX()+1,0,hist->GetNbinsY());
    }
    minIntegral = std::min(minIntegral,histIntegral);
    histIntegrals.push_back(histIntegral);
  }
  
  for(size_t histNr=0;histNr<nrHists;histNr++){
    for(const auto& histVec : histsVec_){
      histVec[histNr].first->Scale(minIntegral/histIntegrals[histNr]);
    }      
  }
  
}
