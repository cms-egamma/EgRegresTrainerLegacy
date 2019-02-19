#include "Utility/HistFuncs.hh"
#include "RegresTrainer/CruijffPdf.h"


#include "RooAbsPdf.h"
#include "RooAddPdf.h"
#include "RooDataHist.h"
#include "RooArgList.h"
#include "RooBreitWigner.h"
#include "RooCBShape.h"
#include "RooDataSet.h"
#include "RooExponential.h"
#include "RooFFTConvPdf.h"
#include "RooGaussian.h"
#include "RooPlot.h"
#include "RooRealVar.h"
#include "RooWorkspace.h"
#include "RooVoigtian.h"
#include "RooHist.h"


#include "TCanvas.h"
#include "TROOT.h"
#include "TStopwatch.h"
#include "TStyle.h"
#include "TLatex.h"
#include "TTree.h"
#include "TGraphErrors.h"
#include "TSystem.h"
#include "TColor.h"
#include "TFitResult.h"

#include "Math/PdfFuncMathCore.h"

#include "RooGlobalFunc.h"

#include <sstream>
#include <iomanip>
size_t gPrecision=3;
namespace resFit{
  enum class FitType{CB,Cruijff};
  FitType fitType = FitType::CB;
}

struct ResFitParam {
  enum class FitType{CB,Cruijff};
  float scale,scaleErr;
  float sigma,sigmaErr;
  float sigmaR,sigmaRErr;
  float sigmaL,sigmaLErr;
  FitType fitType;
  RooPlot* plot;
  std::string legName;
  void fill(const RooRealVar& iScale,const RooRealVar& iSigma,RooPlot* iPlot){
    scale = iScale.getValV();
    scaleErr = iScale.getError();
    sigma = iSigma.getValV();
    sigmaErr = iSigma.getError();
    fitType=FitType::CB;
    plot = iPlot;
  }
  void fill(const RooRealVar& iScale,const RooRealVar& iSigmaL,const RooRealVar& iSigmaR,RooPlot* iPlot){
    scale = iScale.getValV();
    scaleErr = iScale.getError();
    sigmaL = iSigmaL.getValV();
    sigmaLErr = iSigmaL.getError();
    sigmaR = iSigmaR.getValV();
    sigmaRErr = iSigmaR.getError();
    sigma = 0;
    sigmaErr = 0;
    
    fitType=FitType::Cruijff;
    plot = iPlot;
  }
  std::string scaleLabel(){
    std::ostringstream retVal;
    retVal<<std::fixed<<std::setprecision(gPrecision)<<"mean: "<<scale<<" #pm "<<scaleErr<<"^{stat.}";
    return retVal.str();
  }
  std::string scaleAndResCBLabel(){
    std::ostringstream retVal;
    retVal<<std::fixed<<std::setprecision(gPrecision)<<"CB mean: "<<scale<<" CB #sigma: "<<sigma;
    return retVal.str();
  }
  std::string scaleAndResCruijffLabel(){
    std::ostringstream retVal;
    retVal<<std::fixed<<std::setprecision(gPrecision)<<" peak: "<<scale<<" #sigma_{L}: "<<std::endl<<sigmaL<<" #sigma_{R}: "<<sigmaR;
    return retVal.str();
  }
  std::string scaleAndResLabel(){
    if(fitType==FitType::Cruijff) return scaleAndResCruijffLabel();
    else return scaleAndResCBLabel();
  }
  std::string resCBLabel(){
    std::ostringstream retVal;
    retVal<<std::fixed<<std::setprecision(gPrecision)<<"#sigma_{CB}: "<<sigma<<" #pm "<<sigmaErr<<"^{stat.}";
    return retVal.str();
  }
  std::string resCruijffLabel(){
    std::ostringstream retVal;
    retVal<<std::fixed<<std::setprecision(gPrecision)<<"#sigma_{L}: "<<sigmaL<<" #pm "<<sigmaLErr<<"^{stat.}"<<" #sigma_{R}: "<<sigmaR<<" #pm "<<sigmaRErr<<"^{stat.}";
    return retVal.str();
  }
  std::string resLabel(){
    if(fitType==FitType::Cruijff) return resCruijffLabel();
    else return resCBLabel();
  }
  
};

struct CBFitParam {
  float scale,scaleErr;
  float sigma,sigmaErr;
  TH1* hist;
  TF1* func;
  void fill(TH1* iHist,TF1* iFunc){
    scale = iFunc->GetParameter(1);
    scaleErr = iFunc->GetParError(1);
    sigma = iFunc->GetParameter(2);
    sigmaErr = iFunc->GetParError(2);
    hist = iHist;
    func = iFunc;
    
  }  
  std::string scaleLabel(){
    std::ostringstream retVal;
    retVal<<std::fixed<<std::setprecision(gPrecision)<<"#sigma_{mean} "<<scale<<" #pm "<<scaleErr<<"^{stat.}";
    return retVal.str();
  }
  std::string scaleAndResLabel(){
    std::ostringstream retVal;
    retVal<<std::fixed<<std::setprecision(gPrecision)<<"CB mean: "<<scale<<" CB #sigma: "<<sigma;
    return retVal.str();
  }
  std::string resLabel(){
    std::ostringstream retVal;
    retVal<<std::fixed<<std::setprecision(gPrecision)<<"#sigma_{CB} "<<sigma<<" #pm "<<sigmaErr<<"^{stat.}";
    return retVal.str();
  }
};

class CBFitFunc {
public:
  double operator()(double* x,double *p){
    //p0 = norm
    //p1 = mean 
    //p2 = sigma
    //p3 = alpha
    //p4 = n
    return p[0]*ROOT::Math::crystalball_function(*x,p[3],p[4],p[2],p[1]);
  }
};
CBFitParam makeResFitSimple(TH1* hist,float fitMin,float fitMax);
ResFitParam makeResFit(TH1* hist,float fitMin,float fitMax);

namespace cbfit {
  float fitMin=0.8;
  float fitMax=1.1;
  
}

CBFitParam getProjWithFitSimple(TH2* hist2D,int binNr,float fitMin,float fitMax)
{
  std::ostringstream histName;
  histName <<"resHist"<<hist2D->GetXaxis()->GetBinLowEdge(binNr)<<"To"<<hist2D->GetXaxis()->GetBinUpEdge(binNr)<<"GeV";
  TH1* hist = (TH1*) hist2D->ProjectionY(histName.str().c_str(),binNr,binNr)->Clone((histName.str()+"Tmp").c_str());
  hist->SetDirectory(0);
  auto fitRes = makeResFitSimple(hist,cbfit::fitMin,cbfit::fitMax);
  return fitRes;
}

ResFitParam getProjWithFit(TH2* hist2D,int binNr,float fitMin,float fitMax)
{
  std::ostringstream histName;
  histName <<"resHist"<<hist2D->GetXaxis()->GetBinLowEdge(binNr)<<"To"<<hist2D->GetXaxis()->GetBinUpEdge(binNr)<<"GeV";
  TH1* hist = (TH1*) hist2D->ProjectionY(histName.str().c_str(),binNr,binNr)->Clone((histName.str()+"Tmp").c_str());
  hist->SetDirectory(0);
  auto fitRes = makeResFit(hist,fitMin,fitMax);
  return fitRes;
}

void plotResHist(CBFitParam& fitStd,CBFitParam& fitAlt,CBFitParam& fitSC)
{
  //std::vector<std::string> legNames = {"current","proposed","SC energy"};
  std::vector<std::string> legNames = {"raw energy","corr energy 74X","corr energy 94X"};
  

  AnaFuncs::setHistAttributes(fitStd.hist,kAzure+8,2,8,kAzure+8);
  AnaFuncs::setHistAttributes(fitAlt.hist,kOrange+7,2,22,kOrange+7);
  AnaFuncs::setHistAttributes(fitSC.hist,kBlue+2,2,23,kBlue+2);
  fitSC.hist->GetFunction(fitSC.func->GetName())->SetLineColor(kBlue+2);
  //  fitSC.hist->GetFunction(fitSC.func->GetName())->SetLineStyle(2);
  fitAlt.hist->GetFunction(fitAlt.func->GetName())->SetLineColor(kOrange+7);
  //fitAlt.hist->GetFunction(fitAlt.func->GetName())->SetLineStyle(2);
  fitStd.hist->GetFunction(fitStd.func->GetName())->SetLineColor(kAzure+8);
  //fitStd.hist->GetFunction(fitStd.func->GetName())->SetLineStyle(2);
  auto c1 = static_cast<TCanvas*>(gROOT->FindObject("c1"));
  c1->SetLogy(1);
  c1->SetGridx(1);
  // c1->SetGridy(1);

  auto addOverFlow =[](TH1* hist,float xmax){
    int binNr = AnaFuncs::getBinNr(hist,xmax);
    if(xmax==hist->GetBinLowEdge(binNr)) binNr--;
    double intErr=0;
    double intVal = hist->IntegralAndError(binNr,hist->GetNbinsX()+1,intErr);
    hist->SetBinContent(binNr,intVal);
    hist->SetBinError(binNr,intErr);
    for(int i=binNr+1;i<=hist->GetNbinsX()+1;i++){
      hist->SetBinContent(i,0);
      hist->SetBinError(i,0);
    }
  };
  addOverFlow(fitStd.hist,1.5);
  addOverFlow(fitAlt.hist,1.5);
  addOverFlow(fitSC.hist,1.5);

  fitStd.hist->Draw("EP");
  fitStd.hist->GetXaxis()->SetRangeUser(0,1.5);
  fitStd.hist->SetTitle(";E^{reco}/E^{gen};#events");
  fitStd.hist->GetXaxis()->SetTitleOffset(0.9);
  fitAlt.hist->Draw("SAME EP");
  fitSC.hist->Draw("SAME EP");
  
  auto leg = HistFuncs::makeLegend({{fitStd.hist,legNames[0]},{fitAlt.hist,legNames[1]},{fitSC.hist,legNames[2]}},{0.167038,0.56446,0.488864,0.735192});
  leg->SetFillStyle(0);
  leg->Draw();
}

void plotResHist(ResFitParam& fitStd,ResFitParam& fitAlt,ResFitParam& fitSC)
{

  //std::vector<std::string> legNames = {"current","proposed","SC energy"};

  TGraph * histStd = static_cast<TGraph*>(fitStd.plot->getHist("h_res"));
  TGraph * histAlt = static_cast<TGraph*>(fitAlt.plot->getHist("h_res"));
  TGraph * histSC = static_cast<TGraph*>(fitSC.plot->getHist("h_res"));
  std::cout <<"hists "<<histStd<<" "<<histAlt<<" "<<histSC<<std::endl;
  AnaFuncs::setHistAttributes(histStd,kAzure+8,2,8,kAzure+8);
  AnaFuncs::setHistAttributes(histAlt,kOrange+7,2,22,kOrange+7);
  AnaFuncs::setHistAttributes(histSC,kBlue+2,2,23,kBlue+2);
  fitSC.plot->getCurve("model_Norm[res]")->SetLineColor(kBlue+2);
  fitAlt.plot->getCurve("model_Norm[res]")->SetLineColor(kOrange+7);
  fitStd.plot->getCurve("model_Norm[res]")->SetLineColor(kAzure+8);
  //  fitSC.hist->GetFunction(fitSC.func->GetName())->SetLineColor(kBlue+2);
  //  fitSC.hist->GetFunction(fitSC.func->GetName())->SetLineStyle(2);
  //fitAlt.hist->GetFunction(fitAlt.func->GetName())->SetLineColor(kOrange+7);
  //fitAlt.hist->GetFunction(fitAlt.func->GetName())->SetLineStyle(2);
  //fitStd.hist->GetFunction(fitStd.func->GetName())->SetLineColor(kAzure+8);
  //fitStd.hist->GetFunction(fitStd.func->GetName())->SetLineStyle(2);aram.plot->SetTitle("")

  double max = std::max(fitStd.plot->GetMaximum(),fitAlt.plot->GetMaximum());
  max = std::max(fitSC.plot->GetMaximum(),max);
  fitStd.plot->SetMaximum(max);
  fitStd.plot->remove("model_paramBox");
  fitAlt.plot->remove("model_paramBox");
  fitSC.plot->remove("model_paramBox");
  fitStd.plot->SetTitle("");
  fitStd.plot->Draw();
  auto c1 = static_cast<TCanvas*>(gROOT->FindObject("c1"));
// c1->SetLogy(1);
  c1->SetGridx(1);
  c1->SetGridy(1);



  // auto addOverFlow =[](TH1* hist,float xmax){
  //   int binNr = AnaFuncs::getBinNr(hist,xmax);
  //   if(xmax==hist->GetBinLowEdge(binNr)) binNr--;
  //   double intErr=0;
  //   double intVal = hist->IntegralAndError(binNr,hist->GetNbinsX()+1,intErr);
  //   hist->SetBinContent(binNr,intVal);
  //   hist->SetBinError(binNr,intErr);
  //   for(int i=binNr+1;i<=hist->GetNbinsX()+1;i++){
  //     hist->SetBinContent(i,0);
  //     hist->SetBinError(i,0);
  //   }
  // };
  // addOverFlow(fitStd.hist,1.5);
  // addOverFlow(fitAlt.hist,1.5);
  // addOverFlow(fitSC.hist,1.5);

  // fitStd.hist->Draw("EP");
  // fitStd.hist->GetXaxis()->SetRangeUser(0,1.5);
  // fitStd.hist->SetTitle(";E^{reco}/E^{gen};#events");
  // fitStd.hist->GetXaxis()->SetTitleOffset(0.9);
  // fitAlt.hist->Draw("SAME EP");
  // fitSC.hist->Draw("SAME EP");
  
  fitStd.plot->Draw();
  fitAlt.plot->Draw("SAME");
  fitSC.plot->Draw("SAME");

  auto leg = HistFuncs::makeLegend({{fitStd.plot->getHist("h_res"),fitStd.legName},{fitAlt.plot->getHist("h_res"),fitAlt.legName},{fitSC.plot->getHist("h_res"),fitSC.legName}},0.167038,0.56446,0.488864,0.735192);
  leg->SetFillStyle(0);
  leg->Draw();
}


void makeResPlots(TTree* tree,int nrBins,float xmin,float xmax,std::string cuts,int region,const std::string& dir)
{
  if(region==0) cuts+" && ele.region==0"; 
  if(region==1) cuts+" && ele.region==1"; 

  TH2* histStd2D = new TH2D("histStd2D","hist2D",nrBins,xmin,xmax,300,0,3);
  tree->Draw("ele.phoNrgy/truthNrgy:eleTruthEt>>histStd2D",cuts.c_str(),"goff");
  TH2* histAlt2D = new TH2D("histAlt2D","hist2D",nrBins,xmin,xmax,300,0,3);
  tree->Draw("ele.phoNrgyAlt/truthNrgy:eleTruthEt>>histAlt2D",cuts.c_str(),"goff");
  TH2* histSC2D = new TH2D("histSC2D","hist2D",nrBins,xmin,xmax,300,0,3);
  tree->Draw("ele.clusNrgy/truthNrgy:eleTruthEt>>histSC2D",cuts.c_str(),"goff");
  histStd2D->SetDirectory(0);
  histAlt2D->SetDirectory(0);
  histSC2D->SetDirectory(0);
  gSystem->mkdir(dir.c_str(),true);
  
  for(int binNr=1;binNr<nrBins+1;binNr++){
    
    auto fitResStd = getProjWithFit(histStd2D,binNr,cbfit::fitMin,cbfit::fitMax);
    auto fitResAlt = getProjWithFit(histAlt2D,binNr,cbfit::fitMin,cbfit::fitMax);
    auto fitResSC = getProjWithFit(histSC2D,binNr,cbfit::fitMin,cbfit::fitMax);
    
    plotResHist(fitResStd,fitResAlt,fitResSC);
    if(region==0 || region==1){
      //      std::string regionLabelStr = region==0 ? "barrel electrons" : "endcap electrons";
      std::string regionLabelStr = region==0 ? "barrel photons" : "endcap photons";
      auto labelRegion = HistFuncs::makeLabel(regionLabelStr,0.178174,0.820557,0.430958,0.888502);
      labelRegion->Draw();
    }
    std::ostringstream etLabelStr;
    etLabelStr <<histStd2D->GetXaxis()->GetBinLowEdge(binNr)<<" < E_{T} < "<<histStd2D->GetXaxis()->GetBinUpEdge(binNr)<<" GeV";
    
    auto labelEt = HistFuncs::makeLabel(etLabelStr.str(),0.174833,0.740418,0.427617,0.808362);
    labelEt->Draw();
    
    std::ostringstream outputName;
    outputName<<"energyResPho"<<histStd2D->GetXaxis()->GetBinLowEdge(binNr)<<"To"<<histStd2D->GetXaxis()->GetBinUpEdge(binNr)<<"GeV";
    if(region==0) outputName<<"EB";
    if(region==1) outputName<<"EE";
    
    HistFuncs::print(dir+"/"+outputName.str(),"c1");
    
  }
}


std::pair<TGraph*,TGraph*> makeRes(TTree* tree,int nrBins,float xmin,float xmax,const std::string& resVar,const std::string& cuts)
{
  
  TH2* hist2D = new TH2D("hist2D","hist2D",nrBins,xmin,xmax,300,0,3);
  tree->Draw((resVar+":eleTruthEt>>hist2D").c_str(),cuts.c_str(),"goff");
  hist2D->SetDirectory(0);
  std::vector<double> binCentres,scales,scaleErrs,sigmas,sigmaErrs;
  
  for(int binNr=1;binNr<hist2D->GetNbinsX()+1;binNr++){
    std::ostringstream histName;
    histName <<"resHist"<<hist2D->GetXaxis()->GetBinLowEdge(binNr)<<"To"<<hist2D->GetXaxis()->GetBinUpEdge(binNr)<<"GeV";
    TH1* hist = (TH1*) hist2D->ProjectionY(histName.str().c_str(),binNr,binNr)->Clone((histName.str()+"Tmp").c_str());
    hist->SetDirectory(0);
    auto fitRes = makeResFitSimple(hist,cbfit::fitMin,cbfit::fitMax);
    binCentres.push_back(hist2D->GetXaxis()->GetBinCenter(binNr));
    scales.push_back(fitRes.scale);
    scaleErrs.push_back(fitRes.scaleErr);
    sigmas.push_back(fitRes.sigma);
    sigmaErrs.push_back(fitRes.sigmaErr);
  }

  TGraphErrors* scaleGraph = new TGraphErrors(binCentres.size(),&binCentres[0],&scales[0],
					      nullptr,&scaleErrs[0]);
  TGraphErrors* sigmaGraph = new TGraphErrors(binCentres.size(),&binCentres[0],&sigmas[0],
					      nullptr,&sigmaErrs[0]);
					      
  

  return {scaleGraph,sigmaGraph};

}




CBFitParam makeResFitSimple(TH1* hist,float fitMin,float fitMax)
{
  
  TF1* fitFunc = new TF1("fitFunc",new CBFitFunc,0,3,5,"CBFitFunc");
  fitFunc->SetParameters(hist->Integral(),1,0.02,1,1);
  fitFunc->SetParLimits(1,0.9,1.1);
  fitFunc->SetParLimits(2,0.,0.2);
  fitFunc->SetParLimits(3,0.,100);
  fitFunc->SetParLimits(4,0.,100);
  int status = 0;
  int tryNr=0;
  while(status==0 && tryNr<100){
    fitFunc->SetParameter(2,fitFunc->GetParameter(2));
    hist->Fit(fitFunc,"S","goff",0.9,1.1);
    hist->Fit(fitFunc,"S","goff",0.8,1.2);
    hist->Fit(fitFunc,"S","goff",0.95,1.05);
    hist->Fit(fitFunc,"S","goff",fitMin,fitMax);
    status = hist->Fit(fitFunc,"S","goff",fitMin,fitMax)->IsValid();
    tryNr++;
  }
  std::cout <<"status "<<status<<std::endl;
  CBFitParam fitParam;
  fitParam.fill(hist,fitFunc);
  return fitParam;
}



ResFitParam makeResCBFit(TH1* hist,float xmin,float xmax)
{
  RooRealVar  res("res","E^{reco}/E^{gen}", xmin,xmax,"");
  res.setBins(10000,"cache") ;
  res.setMin("cache",xmin) ;
  res.setMax("cache",xmax) ;

  RooRealVar  nsig("N_{S}", "#signal events", 90000, 0, 100000000.);
  RooRealVar  cbSigma("#sigma_{CB}","CB Width", 1.5, 0.0, 10,"");
  RooRealVar mean( "#DeltaE", "mean_{cb}", 1. ,0.5,1.5,"");
  //RooRealVar alpha( "alpha_{cb}", "alpha_{cb}", 1.2 ,0,10);
  //  RooRealVar n( "n_{cb}", "n_{cb}", 0.81 ,0,20);
  RooRealVar alpha( "alpha_{cb}", "alpha_{cb}", 1.2 ,0,20);
  RooRealVar n( "n_{cb}", "n_{cb}", 0.81 ,0,40);
  RooCBShape cb( "cb", "cb",res, mean, cbSigma, alpha, n );
  

  RooDataHist data("res","E^{reco}/E^{gen}",res,hist);
   
  RooAddPdf      model("model", "model", RooArgList(cb), RooArgList(nsig));
  //auto& model = cb;
  model.fitTo(data,RooFit::FitOptions("mh"),RooFit::Optimize(0),RooFit::Timer(1));
  model.fitTo(data,RooFit::FitOptions("mh"),RooFit::Optimize(0),RooFit::Timer(1));

  RooPlot* plot = res.frame(RooFit::Range(xmin,xmax),RooFit::Bins(100));

  data.plotOn(plot,RooFit::MarkerSize(1.0));
  model.plotOn(plot); 


  // model.paramOn(plot,RooFit::Format("NEU", RooFit::AutoPrecision(2)),RooFit::ShowConstants(true),RooFit::Layout(0.6,0.95,0.8),RooFit::Parameters(RooArgSet(mean,cbSigma))); 
  model.paramOn(plot,RooFit::Format("NEU", RooFit::AutoPrecision(2)),RooFit::ShowConstants(true),RooFit::Layout(0.6,0.95,0.8));//,RooFit::Parameters(RooArgSet(mean,cbSigma))); 

  ResFitParam fitParam;
  fitParam.fill(mean,cbSigma,plot);
  return fitParam;
}

ResFitParam makeResCruijffFit(TH1* hist,float xmin,float xmax)
{
  RooRealVar  res("res","E^{reco}/E^{gen}", xmin,xmax,"");
  res.setBins(10000,"cache") ;
  res.setMin("cache",xmin) ;
  res.setMax("cache",xmax) ;

  RooRealVar  nsig("N_{S}", "#signal events", 90000, 0, 100000000.);
  RooRealVar mean( "#DeltaE", "mean_{cb}", 1. ,0.5,1.5,"");
  RooRealVar sigmaL("#sigma_{L}","#sigma_{L}", 0.02, 0.0, 0.5);
  RooRealVar sigmaR("#sigma_{R}","#sigma_{R}", 0.02, 0.0, 0.5);
  RooRealVar alphaL( "alpha_{L}", "alpha_{L}", 0.1 ,0,2);
  RooRealVar alphaR( "alpha_{R}", "alpha_{R}", 0.1,0,2);
  CruijffPdf cruijff("cruijff","cruijff",res,mean,sigmaL,sigmaR,alphaL,alphaR);

  RooDataHist data("res","E^{reco}/E^{gen}",res,hist);
   
  RooAddPdf      model("model", "model", RooArgList(cruijff), RooArgList(nsig));
  //auto& model = cb;
  model.fitTo(data,RooFit::FitOptions("mh"),RooFit::Optimize(0),RooFit::Timer(1));
  model.fitTo(data,RooFit::FitOptions("mh"),RooFit::Optimize(0),RooFit::Timer(1));

  RooPlot* plot = res.frame(RooFit::Range(xmin,xmax),RooFit::Bins(100));

  data.plotOn(plot,RooFit::MarkerSize(1.0));
  model.plotOn(plot); 


  // model.paramOn(plot,RooFit::Format("NEU", RooFit::AutoPrecision(2)),RooFit::ShowConstants(true),RooFit::Layout(0.6,0.95,0.8),RooFit::Parameters(RooArgSet(mean,cbSigma))); 
  model.paramOn(plot,RooFit::Format("NEU", RooFit::AutoPrecision(2)),RooFit::ShowConstants(true),RooFit::Layout(0.6,0.95,0.8));//,RooFit::Parameters(RooArgSet(mean,cbSigma))); 

  ResFitParam fitParam;
  fitParam.fill(mean,sigmaL,sigmaR,plot);
  return fitParam;
}



ResFitParam makeResFit(TH1* hist,float xmin,float xmax)
{
  if(resFit::fitType == resFit::FitType::CB){
    return makeResCBFit(hist,xmin,xmax);
  }else{
    return makeResCruijffFit(hist,xmin,xmax);
  }
}



