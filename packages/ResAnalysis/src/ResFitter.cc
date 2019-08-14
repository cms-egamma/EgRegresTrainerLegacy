#include "ResAnalysis/ResFitter.hh"
#include "RegresTrainer/CruijffPdf.h"
#include "GBRLikelihood/RooDoubleCBFast.h"
#include "Utility/HistFuncs.hh"

#include "RooAddPdf.h"
#include "RooConstVar.h"
#include "RooDataHist.h"
#include "RooArgList.h"
#include "RooCBShape.h"
#include "RooPlot.h"
#include "RooRealVar.h"
#include "RooFitResult.h"
#include "RooWorkspace.h"
#include "RooHist.h"

#include "TH1.h"
#include "TH2.h"
#include "TGraphErrors.h"

ResFitter::Param ResFitter::makeFit(TH1* hist,float xmin,float xmax,const std::string& fitVarName)const
{
  switch(fitType_){
  case FitType::CB: return makeCBFit(hist,xmin,xmax,fitVarName);
  case FitType::DCB: return makeDCBFit(hist,xmin,xmax,fitVarName);
  case FitType::Cruijff: return makeCruijffFit(hist,xmin,xmax,fitVarName);
  }
  
  return ResFitter::Param();
}

ResFitter::Param ResFitter::makeFit(TH2* hist2D,int binNr,float xmin,float xmax,const std::string& fitVarName)const
{
  std::ostringstream histName;
  histName <<"resHist"<<hist2D->GetXaxis()->GetBinLowEdge(binNr)<<"To"<<hist2D->GetXaxis()->GetBinUpEdge(binNr)<<"GeV";
  TH1* hist = static_cast<TH1*>(hist2D->ProjectionY(histName.str().c_str(),binNr,binNr)->Clone((histName.str()+"Tmp").c_str()));
  hist->SetDirectory(0);
  return makeFit(hist,xmin,xmax,fitVarName);
}

ResFitter::Param ResFitter::makeCBFit(TH1* hist,float xmin,float xmax,const std::string& fitVarName)const
{
  RooRealVar  res("res","E^{reco}/E^{gen}", xmin,xmax,"");
  res.setBins(10000,"cache") ;
  res.setMin("cache",xmin) ;
  res.setMax("cache",xmax) ;

  RooRealVar  nsig("N_{S}", "#signal events", 90000, 0, 100000000.);
  RooRealVar  cbSigma("#sigma_{CB}","CB Width", 1.5, 0.0, 10,"");
  RooRealVar mean( "#DeltaE", "mean_{cb}", 1. ,0.5,1.5,"");
  RooRealVar alpha( "alpha_{cb}", "alpha_{cb}", 1.2 ,0,20);
  RooRealVar n( "n_{cb}", "n_{cb}", 0.81 ,0,40);
  RooCBShape cb( "cb", "cb",res, mean, cbSigma, alpha, n );
  

  RooDataHist data("res","E^{reco}/E^{gen}",res,hist);
   
  RooAddPdf      model("model", "model", RooArgList(cb), RooArgList(nsig));
  model.fitTo(data,RooFit::FitOptions("mh"),RooFit::Optimize(0),RooFit::Timer(1));
  auto fitRes = model.fitTo(data,RooFit::FitOptions("mh"),RooFit::Optimize(0),RooFit::Timer(1));
  std::cout <<" fit status "<<fitRes->status()<<std::endl;
  RooPlot* plot = res.frame(RooFit::Range(xmin,xmax),RooFit::Bins(100));

  data.plotOn(plot,RooFit::MarkerSize(1.0));
  model.plotOn(plot); 

  model.paramOn(plot,RooFit::Format("NEU", RooFit::AutoPrecision(2)),RooFit::ShowConstants(true),RooFit::Layout(0.6,0.95,0.8));

  ResFitter::Param fitParam;
  fitParam.fill(mean,cbSigma,FitType::CB,plot,fitVarName);
  return fitParam;
}


ResFitter::Param ResFitter::makeDCBFit(TH1* hist,float xmin,float xmax,const std::string& fitVarName)const
{
  RooRealVar  res("res","E^{reco}/E^{gen}", xmin,xmax,"");
  res.setBins(10000,"cache") ;
  res.setMin("cache",xmin) ;
  res.setMax("cache",xmax) ;

  float nrInRange = AnaFuncs::getHistIntegral(hist,xmin,xmax);

  RooRealVar  nsig("N_{S}", "#signal events", nrInRange,nrInRange*0.9,nrInRange*1.1);
  RooRealVar mean( "#DeltaE", "mean_{cb}", 1. ,0.5,1.5,""); 
  if(fixMeanDCB_) mean.setRange(1,1);
  RooRealVar cbSigma("#sigma_{CB}","CB Width", 0.05, 0.0002, 0.5,"");
  RooRealVar alpha1( "alpha_{1}", "alpha_{1}", 1.2 ,0,20);
  //RooRealVar n1( "n_{1}", "n_{1}", 3 ,0,40);
  RooRealVar n1( "n_{1}", "n_{1}", 2 ,1.01,5000.);
  RooRealVar alpha2( "alpha_{2}", "alpha_{2}", 1.2 ,0,20);
  //  RooRealVar n2( "n_{2}", "n_{2}", 0.81 ,0,40);
  RooRealVar n2( "n_{2}", "n_{2}", 2 ,1.01,5000.);

  RooDoubleCBFast dcb = fixAlphaDCB_ ? 
    RooDoubleCBFast("dcb","dcb",res,mean,cbSigma,RooFit::RooConst(1.),n1,RooFit::RooConst(2.),n2):
    RooDoubleCBFast("dcb","dcb",res,mean,cbSigma,alpha1,n1,alpha2,n2);

  RooDataHist data("res","E^{reco}/E^{gen}",res,hist);
   
  RooAddPdf      model("model", "model", RooArgList(dcb), RooArgList(nsig));
  //auto& model = cb;
  model.fitTo(data,RooFit::FitOptions("mh"),RooFit::Optimize(0),RooFit::Timer(1));
  model.fitTo(data,RooFit::FitOptions("mh"),RooFit::Optimize(0),RooFit::Timer(1));

  RooPlot* plot = res.frame(RooFit::Range(xmin,xmax),RooFit::Bins(100));

  data.plotOn(plot,RooFit::MarkerSize(1.0));
  model.plotOn(plot); 

  model.paramOn(plot,RooFit::Format("NEU", RooFit::AutoPrecision(2)),RooFit::ShowConstants(true),RooFit::Layout(0.6,0.95,0.8));

  ResFitter::Param fitParam;
  fitParam.fill(mean,cbSigma,FitType::DCB,plot,fitVarName);
  return fitParam;
}

ResFitter::Param ResFitter::makeCruijffFit(TH1* hist,float xmin,float xmax,const std::string& fitVarName)const
{
  RooRealVar  res("res","E^{reco}/E^{gen}", xmin,xmax,"");
  res.setBins(10000,"cache") ;
  res.setMin("cache",xmin) ;
  res.setMax("cache",xmax) ;
  float nrEntries = AnaFuncs::getHistIntegral(hist,xmin,xmax);
  RooRealVar  nsig("N_{S}", "#signal events", nrEntries, nrEntries*0.5, nrEntries*2);
  RooRealVar mean( "#DeltaE", "mean_{cb}", 1. ,0.8,1.2,"");
  RooRealVar sigmaL("#sigma_{L}","#sigma_{L}", 0.02, 0.001, 0.5);
  RooRealVar sigmaR("#sigma_{R}","#sigma_{R}", 0.02, 0.001, 0.5);
  RooRealVar alphaL( "alpha_{L}", "alpha_{L}", 0.1,0,2);
  RooRealVar alphaR( "alpha_{R}", "alpha_{R}", 0.1,0,2);
  CruijffPdf cruijff("cruijff","cruijff",res,mean,sigmaL,sigmaR,alphaL,alphaR);

  RooDataHist data("res","E^{reco}/E^{gen}",res,hist);
   
  RooAddPdf model("model", "model", RooArgList(cruijff), RooArgList(nsig));
  model.fitTo(data,RooFit::FitOptions("mh"),RooFit::Optimize(0),RooFit::Timer(1),RooFit::PrintEvalErrors(-1));
  auto fitRes = model.fitTo(data,RooFit::Optimize(1),RooFit::Timer(0),RooFit::PrintEvalErrors(-1),RooFit::Save(1));
  if( (sigmaL.getValV()<=0.001 || sigmaR.getValV()<=0.001) || fitRes->edm()>10){
    std::cout <<" fit status "<<fitRes->status()<<" sigma "<<sigmaL.getValV()<<" "<<sigmaR.getValV()<<" nrsig "<<nsig.getValV()<<" edm "<<fitRes->edm()<<std::endl;
    double maxSigma = std::max(sigmaL.getValV(),sigmaR.getValV());
    sigmaL.setVal(maxSigma);
    sigmaR.setVal(maxSigma);
    fitRes = model.fitTo(data,RooFit::Optimize(1),RooFit::Timer(0),RooFit::PrintEvalErrors(-1),RooFit::Save(1));
    std::cout <<"trying again "<<fitRes->status()<<" sigma "<<sigmaL.getValV()<<" "<<sigmaR.getValV()<<" nrsig "<<nsig.getValV()<<" edm "<<fitRes->edm()<<std::endl;
  }
  RooPlot* plot = res.frame(RooFit::Range(xmin,xmax),RooFit::Bins(100));

  data.plotOn(plot,RooFit::MarkerSize(1.0));
  model.plotOn(plot); 

  model.paramOn(plot,RooFit::Format("NEU", RooFit::AutoPrecision(2)),RooFit::ShowConstants(true),RooFit::Layout(0.6,0.95,0.8));

  if(alphaL.getValV()>=2 || alphaR.getValV()>=2) std::cout <<"alphaL "<<alphaL.getValV()<<" alphaR "<<alphaR.getValV()<<std::endl;

  ResFitter::Param fitParam;
  fitParam.fill(mean,sigmaL,sigmaR,plot,fitVarName);
 
  //  std::cout <<" res NLL "<<fitRes->minNll()<<" edm "<<fitRes->edm()<<std::endl;

  return fitParam;
}

ResFitter::ParamsVsVar ResFitter::makeFitVsVar(TH2* hist2D,float fitMin,float fitMax,const std::string& fitVarName)const
{
  std::vector<Param> fitParams;
  std::vector<double> binLowEdges;
  for(int binNr=1;binNr<=hist2D->GetNbinsX();binNr++){
    fitParams.push_back(makeFit(hist2D,binNr,fitMin,fitMax,fitVarName));
    binLowEdges.push_back(hist2D->GetXaxis()->GetBinLowEdge(binNr));
  }
  binLowEdges.push_back(hist2D->GetXaxis()->GetBinLowEdge(hist2D->GetNbinsX()+1));
  return ParamsVsVar(std::move(fitParams),std::move(binLowEdges));
}

ResFitter::Param::Param():
  scale(0.),scaleErr(0.),
  sigma(0.),sigmaErr(0.),
  sigmaR(0.),sigmaRErr(0.),
  sigmaL(0.),sigmaLErr(0.),
  plot(nullptr),
  legName(),
  fitType(FitType::CB),
  precision(3)
{

}

void ResFitter::Param::fill(const RooRealVar& iScale,const RooRealVar& iSigma,ResFitter::FitType iFitType,
			    RooPlot* iPlot,const std::string& iLegName){
  scale = iScale.getValV();
  scaleErr = iScale.getError();
  sigma = iSigma.getValV();
  sigmaErr = iSigma.getError();
  fitType=iFitType;
  plot = iPlot;
  legName = iLegName;
}

void ResFitter::Param::fill(const RooRealVar& iScale,const RooRealVar& iSigmaL,const RooRealVar& iSigmaR,
			    RooPlot* iPlot,const std::string& iLegName){
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
  legName = iLegName;
}

TGraph* ResFitter::ParamsVsVar::makeGraph(ValType valType,bool divideSigmaByMean)const
{
  std::vector<float> xValues;
  std::vector<float> xErrs;
  std::vector<float> yValues;
  std::vector<float> yErrs;
  std::string yAxisTitle="";

  //params doesnt have under/overflows
  for(size_t binNr=0;binNr<params_.size();binNr++){
    Param fitParam = params_[binNr];
    FitType fitType = fitParam.fitType;
    if(binNr+1<params_.size() && fitType!=params_[binNr+1].fitType){
      LogErr << "Warning parameters have different fit types. This should be impossible and should be fixed"<<std::endl;
    }
    float binMin = binLowEdges_[binNr];
    float binMax = binLowEdges_[binNr+1];
    float binHalfWidth = (binMax-binMin)/2.;
    xValues.push_back(binMin+binHalfWidth);
    xErrs.push_back(binHalfWidth);

    std::pair<float,float> yVal = {0.,0.};
   
    switch(valType){
    case ValType::Mean:
      yVal = {fitParam.scale,fitParam.scaleErr};
      yAxisTitle = fitType!=FitType::Cruijff ? "CB mean" : "Cruijff mean";
      break;
    case ValType::Sigma:
      if (fitType!=FitType::Cruijff){
	yVal = {fitParam.sigma,fitParam.sigmaErr};
	yAxisTitle = "CB #sigma";
      }else{
	auto sqrtOfSqSum = [](float x,float y){return std::sqrt(x*x+y*y);};
	yVal = {(fitParam.sigmaL+fitParam.sigmaR)/2.,sqrtOfSqSum(fitParam.sigmaLErr,fitParam.sigmaRErr)/2.};
	yAxisTitle = "Cruijff #sigma_{ave}";
      }
      break;
    case ValType::SigmaL:
      yVal = {fitParam.sigmaL,fitParam.sigmaLErr};
      yAxisTitle = "Cruiff #sigma_{L}";
      break;
    case ValType::SigmaR:
      yVal = {fitParam.sigmaR,fitParam.sigmaRErr};
      yAxisTitle = "Cruiff #sigma_{R}";
      break;
    };
    if(valType==ValType::Mean || !divideSigmaByMean){
      yValues.push_back(yVal.first);
      yErrs.push_back(yVal.second);
    }else{
      yAxisTitle+=" / #mu";
      yValues.push_back(yVal.first/fitParam.scale);
      yErrs.push_back(yVal.second/fitParam.scale);
    }
  }
  auto graph = new TGraphErrors(xValues.size(),xValues.data(),yValues.data(),xErrs.data(),yErrs.data());
  graph->GetYaxis()->SetTitle(yAxisTitle.c_str());
  return graph;
}
