#include "Utility/HistFuncs.hh"
#include "Utility/AnaFuncs.hh"
#include "Utility/TempFuncs.hh"
#include "Utility/LogErr.hh"
#include "Utility/MathFuncs.hh"


#include "Math/QuantFuncMathCore.h"
#include "TLegend.h"
#include "TProfile.h"
#include "TCanvas.h"
#include "TROOT.h"
#include "TGraphErrors.h"
#include "TEventList.h"

#include <boost/algorithm/string/replace.hpp>
#include <boost/algorithm/string.hpp>

int HistFuncs::nrBinsXForEffHist_ = 20;
float HistFuncs::xMinForEffHist_ = 25;
float HistFuncs::xMaxForEffHist_ = 1025;

bool HistFuncs::plotXErrOnAsymEff_=true;

//should really be two functions
void HistFuncs::plotHistProjections(TH2* inputHist,int nrProjections,int maxNrToDisplay,const char* title,bool scale)
{
  int nrBins = inputHist->GetNbinsX();
  int nrBinsPerProj = nrBins/nrProjections;

  std::vector<TH1*> hists;
  std::vector<std::pair<double,double> > binRanges;
  for(int projNr=0;projNr<nrProjections;projNr++){
    int minBinNr = 1+nrBinsPerProj*projNr;
    int maxBinNr = nrBinsPerProj*(projNr+1);
    if(maxBinNr>nrBins) maxBinNr = nrBins;
    
    TString histName(inputHist->GetName());
    histName+="_proj";
    histName+=projNr;
    
    TH1* tempHist = inputHist->ProjectionY("temp",minBinNr,maxBinNr);
    
    TH1* histProj = (TH1*) tempHist->Clone(histName.Data());
    histProj->SetDirectory(0);

    hists.push_back(histProj);
    binRanges.push_back(std::make_pair<double,double>(inputHist->GetBinLowEdge(minBinNr),inputHist->GetBinLowEdge(maxBinNr+1)));

  }

  //now plot it
  hists[0]->Draw("HIST");
  hists[0]->SetTitle(title);
  hists[0]->GetXaxis()->SetTitleSize(0.047);
  hists[0]->GetYaxis()->SetTitleSize(0.047);
  hists[0]->GetXaxis()->SetTitleOffset(0.9);
  hists[0]->GetYaxis()->SetTitleOffset(1.3);
  TLegend *leg = new TLegend(0.3,0.3,0.5,0.5);
  for(size_t histNr=0;histNr<hists.size() && histNr<(unsigned)maxNrToDisplay;histNr++){
    int markerNr = histNr%2==0 ? 22 : 23;
    AnaFuncs::setHistAttributes(hists[histNr],getColourNr(histNr),2,markerNr,getColourNr(histNr));
    hists[histNr]->Sumw2();
    int lastHistNr = hists.size() < (unsigned) maxNrToDisplay ? hists.size() : maxNrToDisplay;
    if(scale) hists[histNr]->Scale(hists[lastHistNr-1]->GetEntries()/hists[histNr]->GetEntries());

    hists[histNr]->Draw("HIST SAME");
    std::ostringstream legText;
    legText << binRanges[histNr].first <<" < E_{T} < "<<binRanges[histNr].second<<" GeV";
    leg->AddEntry(hists[histNr],legText.str().c_str(),"L");
    
  }
  leg->SetFillColor(0);
  leg->SetBorderSize(0);
  leg->Draw();
}

int HistFuncs::getColourNr(int indx)
{
  const int nrColours=6;
  const int colourNrs[nrColours]={1,4,2,16,7,6};
  int colourNr = indx%nrColours; //so it repeats after nrColours
  return colourNrs[colourNr];
}

int HistFuncs::getMarkerNr(int indx)
{
  const int nrMarkers=4;
  const int markerNrs[nrMarkers]={8,4,22,23};
  int markerNr = indx%nrMarkers; //so it repeats after nrMarkers
  return markerNrs[markerNr];
}





TChain* HistFuncs::makeChain(const std::string& chainName,std::string filelist,int nrJobs,int jobNr,int verbose)
{
  return makeChain(chainName,std::vector<std::string>{filelist},nrJobs,jobNr,verbose);
}

TChain* HistFuncs::makeChain(const std::string& chainName,std::vector<std::string> filelists,int nrJobs,int jobNr,int verbose)
{
  std::vector<std::string> filenames;
  AnaFuncs::readFilelist(filelists,filenames,nrJobs,jobNr,verbose); 
  TChain* chain = new TChain(chainName.c_str(),"chain");
  for(size_t i=0;i<filenames.size();i++){
    chain->Add(filenames[i].c_str());
  }
  return chain;
}

TH1* HistFuncs::getHist(const std::string& histName,const std::string& filename)
{
  TFile* file = TFile::Open(filename.c_str(),"READ");
  TH1* hist = (TH1*) file->Get(histName.c_str());
  if(hist) hist->SetDirectory(0);
  delete file;
  return hist;
}

TH1* HistFuncs::makeEffHistFromTree(TTree* tree,int nrBins,float xMin,float xMax,const std::string& var,const std::string& sampleCuts,const std::string& cuts)
{
  TH1* pass = new TH1D("passTemp","pass",nrBins,xMin,xMax);
  TH1* all = new TH1D("allTemp","all",nrBins,xMin,xMax);
  pass->Sumw2();
  all->Sumw2();
  
  float nrPass = tree->Draw((var+">>passTemp").c_str(),(sampleCuts+" && "+cuts).c_str(),"goff");
  float nrAll = tree->Draw((var+">>allTemp").c_str(),(sampleCuts).c_str(),"goff");
  
  all->SetDirectory(0);
  pass->SetDirectory(0);

  pass->Divide(pass,all,1,1,"B");
  std::cout <<"nrPass "<<nrPass<< " / "<<nrAll<<std::endl;

  delete all;

  pass->SetTitle((";"+vsVarAxisLabel(var)+";Efficiency").c_str());
  
  return pass;
 
}

TH1* HistFuncs::makeEffHistBothEles(TTree* tree,int nrBins,float xMin,float xMax,const std::string& var,const std::string& sampleCuts,const std::string& cuts)
{
  
  TH1* all = makeHistBothEles(tree,nrBins,xMin,xMax,var,sampleCuts);
  TH1* pass = makeHistBothEles(tree,nrBins,xMin,xMax,var,sampleCuts+" && "+cuts);

  float nrPass = pass->GetEntries();
  float nrAll = all->GetEntries();
  
  all->SetDirectory(0);
  pass->SetDirectory(0);

  pass->Divide(pass,all,1,1,"B");
  std::cout <<"nrPass "<<nrPass<< " / "<<nrAll<<std::endl;

  delete all;

  pass->SetTitle((";"+vsVarAxisLabel(var)+";Efficiency").c_str());
  
  return pass;
 
}

TH1* HistFuncs::makeEffHistBothEles(TTree* tree,int nrBins,float xMin,float xMax,const std::string& var,const std::string& sampleCuts,const std::string& cuts,const std::string& weight)
{
  
  TH1* all = makeHistBothEles(tree,nrBins,xMin,xMax,var,"("+sampleCuts+")*"+weight);
  TH1* pass = makeHistBothEles(tree,nrBins,xMin,xMax,var,"("+sampleCuts+" && "+cuts+")*"+weight);
  
  float nrPass = pass->GetEntries();
  float nrAll = all->GetEntries();
  
  all->SetDirectory(0);
  pass->SetDirectory(0);

  pass->Divide(pass,all,1,1,"B");
  std::cout <<"nrPass "<<nrPass<< " / "<<nrAll<<std::endl;

  delete all;

  pass->SetTitle((";"+vsVarAxisLabel(var)+";Efficiency").c_str());
  
  return pass;
 
}


TH2* HistFuncs::makeEffHistFromTree(TTree* tree,int nrBinsX,float xMin,float xMax,int nrBinsY,float yMin,float yMax,const std::string& varX,const std::string& varY,const std::string& sampleCuts,const std::string& cuts)
{
  TH2* pass = new TH2D("passTemp","pass",nrBinsX,xMin,xMax,nrBinsY,yMin,yMax);
  TH2* all = new TH2D("allTemp","all",nrBinsX,xMin,xMax,nrBinsY,yMin,yMax);
  pass->Sumw2();
  all->Sumw2();
  
  float nrPass = tree->Draw((varY+":"+varX+">>passTemp").c_str(),(sampleCuts+" && "+cuts).c_str(),"goff");
  float nrAll = tree->Draw((varY+":"+varX+">>allTemp").c_str(),(sampleCuts).c_str(),"goff");
  
  all->SetDirectory(0);
  pass->SetDirectory(0);

  pass->Divide(pass,all,1,1,"B");
  pass->SetTitle((";"+vsVarAxisLabel(varX)+";"+vsVarAxisLabel(varY)).c_str());
  std::cout <<"nrPass "<<nrPass<< " / "<<nrAll<<std::endl;

  delete all;
  
  return pass;

}


TH2* HistFuncs::makeEffHistFromTree(TTree* tree,int nrBinsX,float xMin,float xMax,const std::vector<double>& yBins,const std::string& varX,const std::string& varY,const std::string& sampleCuts,const std::string& cuts)
{
  TH2* pass = new TH2D("passTemp","pass",nrBinsX,xMin,xMax,yBins.size()-1,&yBins[0]);
  TH2* all = new TH2D("allTemp","all",nrBinsX,xMin,xMax,yBins.size()-1,&yBins[0]);
  pass->Sumw2();
  all->Sumw2();
  
  float nrPass = tree->Draw((varY+":"+varX+">>passTemp").c_str(),(sampleCuts+" && "+cuts).c_str(),"goff");
  float nrAll = tree->Draw((varY+":"+varX+">>allTemp").c_str(),(sampleCuts).c_str(),"goff");
  
  all->SetDirectory(0);
  pass->SetDirectory(0);

  pass->Divide(pass,all,1,1,"B");
  pass->SetTitle((";"+vsVarAxisLabel(varX)+";"+vsVarAxisLabel(varY)).c_str());
  std::cout <<"nrPass "<<nrPass<< " / "<<nrAll<<std::endl;

  delete all;
  
  return pass;
 
}

TH1* HistFuncs::makeEffHistFromTree1DProj(TTree* tree,int nrBinsX,float xMin,float xMax,int nrBinsY,float yMin,float yMax,const std::string& varX,const std::string& varY,const std::string& sampleCuts,const std::string& cuts)
{
  TH2* hist2D = makeEffHistFromTree(tree,nrBinsX,xMin,xMax,nrBinsY,yMin,yMax,varX,varY,sampleCuts,cuts);
  const std::string& varName = getNiceName(varY);
  const std::string& varUnits = getUnits(varY);
  TLegend* leg = new TLegend(0.3,0.4,0.5,0.6);
  leg->SetBorderSize(0);
  
  TH1* axisHist=nullptr;
  for(int binNr=1;binNr<=nrBinsY;binNr++){
    TH1* hist=hist2D->ProjectionX(("histBinNr"+AnaFuncs::convertToStr(binNr)).c_str(),binNr,binNr);
    hist->SetDirectory(0);
    
    int colour = getColourNr(binNr);
    int marker = getMarkerNr(binNr);
    AnaFuncs::setHistAttributes(hist,colour,1,marker,colour);
    if(!axisHist){
      axisHist = hist;
      hist->Draw("EP");
    }else hist->Draw("SAME EP");
    
    std::ostringstream label;
    label<<AnaFuncs::getBinLowEdge(nrBinsY,yMin,yMax,binNr)<<" #leq "<<varName<<" < "<<AnaFuncs::getBinLowEdge(nrBinsY,yMin,yMax,binNr+1)<<" "<<varUnits;
    leg->AddEntry(hist,label.str().c_str(),"LP");
    

  }
  leg->Draw();
  delete hist2D;
  
  return axisHist;
  
}


TH1* HistFuncs::makeEffHistFromTree1DProj(TTree* tree,int nrBinsX,float xMin,float xMax,const std::vector<double>& yBins,const std::string& varX,const std::string& varY,const std::string& sampleCuts,const std::string& cuts)
{
  TH2* hist2D = makeEffHistFromTree(tree,nrBinsX,xMin,xMax,yBins,varX,varY,sampleCuts,cuts);
  const std::string& varName = getNiceName(varY);
  const std::string& varUnits = getUnits(varY);
  TLegend* leg = new TLegend(0.3,0.4,0.5,0.6);
  leg->SetBorderSize(0);
  
  TH1* axisHist=nullptr;
  for(size_t binNr=0;binNr+1<yBins.size();binNr++){
    TH1* hist=hist2D->ProjectionX(("histBinNr"+AnaFuncs::convertToStr(binNr+1)).c_str(),binNr+1,binNr+1);
    hist->SetDirectory(0);
    
    int colour = getColourNr(binNr);
    int marker = getMarkerNr(binNr);
    AnaFuncs::setHistAttributes(hist,colour,1,marker,colour);
    if(!axisHist){
      axisHist = hist;
      hist->Draw("EP");
    }else hist->Draw("SAME EP");
    
    std::ostringstream label;
    label<<yBins[binNr]<<" #leq "<<varName<<" < "<<yBins[binNr+1]<<" "<<varUnits;
    leg->AddEntry(hist,label.str().c_str(),"LP");
    

  }
  leg->Draw();
  delete hist2D;
  
  return axisHist;
  
}

TH1* HistFuncs::makeIntEffHistFromTree(TTree* tree,int nrBins,float xMin,float xMax,const std::string& var,const std::string& sampleCuts,const std::string& cuts,bool intIsGreatThan)
{
  TH1* pass = new TH1D("passTemp","pass",nrBins,xMin,xMax);
  TH1* all = new TH1D("allTemp","all",nrBins,xMin,xMax);
  pass->Sumw2();
  all->Sumw2();
  
  float nrPass = tree->Draw((var+">>passTemp").c_str(),(sampleCuts+" && "+cuts).c_str(),"goff");
  float nrAll = tree->Draw((var+">>allTemp").c_str(),(sampleCuts).c_str(),"goff");
  
  all->SetDirectory(0);
  pass->SetDirectory(0);

  TH1* intPass = makeCHist(pass,intIsGreatThan);
  TH1* intAll = makeCHist(all,intIsGreatThan);

  


  intPass->Divide(intPass,intAll,1,1,"B");
  std::cout <<"nrPass "<<nrPass<< " / "<<nrAll<<std::endl;

  delete all;
  delete pass;
  delete intAll;
  
  return intPass;
 
}


TH1* HistFuncs::makeEffHistFromTree(TTree* tree,int nrBins,float xMin,float xMax,const std::string& var,const std::string& sampleCuts,const std::string& cuts,const std::string& weight)
{
  TH1* pass = new TH1D("pass","pass",nrBins,xMin,xMax);
  TH1* all = new TH1D("all","all",nrBins,xMin,xMax);
  pass->Sumw2();
  all->Sumw2();
  
  float nrPass = tree->Draw((var+">>pass").c_str(),("("+sampleCuts+" && "+cuts+")*"+weight).c_str(),"goff");
  float nrAll = tree->Draw((var+">>all").c_str(),("("+sampleCuts+")*"+weight).c_str(),"goff");
  
  all->SetDirectory(0);
  pass->SetDirectory(0);

  pass->Divide(pass,all,1,1,"B");
  
  std::cout <<"nrPass "<<nrPass<< " / "<<nrAll<<std::endl;
  pass->SetTitle((";"+vsVarAxisLabel(var)+";Efficiency").c_str());

  delete all;

  return pass;
 
}


TH1* HistFuncs::makeEffHistFromTree(TTree* tree,std::vector<float> bins,const std::string& var,const std::string& sampleCuts,const std::string& cuts,const std::string& weight)
{
  TH1* pass = new TH1D("pass","pass",bins.size()-1,&bins[0]);
  TH1* all = new TH1D("all","all",bins.size()-1,&bins[0]);
  pass->Sumw2();
  all->Sumw2();
  
  float nrPass = tree->Draw((var+">>pass").c_str(),("("+sampleCuts+" && "+cuts+")*"+weight).c_str(),"goff");
  float nrAll = tree->Draw((var+">>all").c_str(),("("+sampleCuts+")*"+weight).c_str(),"goff");
  
  all->SetDirectory(0);
  pass->SetDirectory(0);

  pass->Divide(pass,all,1,1,"B");
  
  std::cout <<"nrPass "<<nrPass<< " / "<<nrAll<<std::endl;

  pass->SetTitle((";"+vsVarAxisLabel(var)+";Efficiency").c_str());
  delete all;

  return pass;
 
}

TGraph* HistFuncs::makeEffHistFromTreeAsymErr(TTree* tree,int nrBins,float xMin,float xMax,const std::string& var,const std::string& sampleCuts,const std::string& cuts)
{
  TH1* pass = new TH1D("pass","pass",nrBins,xMin,xMax);
  TH1* all = new TH1D("all","all",nrBins,xMin,xMax);
  pass->Sumw2();
  all->Sumw2();
  
 float nrPass = tree->Draw((var+">>pass").c_str(),(sampleCuts+" && "+cuts).c_str(),"goff");
  float nrAll = tree->Draw((var+">>all").c_str(),(sampleCuts).c_str(),"goff");
  
  all->SetDirectory(0);
  pass->SetDirectory(0);

  TGraphAsymmErrors* graph = new TGraphAsymmErrors(pass,all,"B");
  if(!plotXErrOnAsymEff_){
    for(int pointNr=0;pointNr<graph->GetN();pointNr++){
      graph->SetPointEXhigh(pointNr,0);
      graph->SetPointEXlow(pointNr,0);
    }
  }

  delete all;
  delete pass;
  
  std::cout <<"nrPass "<<nrPass<< " / "<<nrAll<<std::endl;

  graph->SetTitle((";"+vsVarAxisLabel(var)+";Efficiency").c_str());
  return graph;
 
}


TGraph* HistFuncs::makeEffHistFromTreeAsymErr(TTree* tree,std::vector<float> bins,const std::string& var,const std::string& sampleCuts,const std::string& cuts)
{
  TH1* pass = new TH1D("pass","pass",bins.size()-1,&bins[0]);
  TH1* all = new TH1D("all","all",bins.size()-1,&bins[0]);
  pass->Sumw2();
  all->Sumw2();
  
 float nrPass = tree->Draw((var+">>pass").c_str(),(sampleCuts+" && "+cuts).c_str(),"goff");
  float nrAll = tree->Draw((var+">>all").c_str(),(sampleCuts).c_str(),"goff");
  
  all->SetDirectory(0);
  pass->SetDirectory(0);

  TGraphAsymmErrors* graph = new TGraphAsymmErrors(pass,all,"B");

  if(!plotXErrOnAsymEff_){
    for(int pointNr=0;pointNr<graph->GetN();pointNr++){
      graph->SetPointEXhigh(pointNr,0);
      graph->SetPointEXlow(pointNr,0);
    }
  }

 
  delete all;
  delete pass;
  
  std::cout <<"nrPass "<<nrPass<< " / "<<nrAll<<std::endl;
  graph->SetTitle((";"+vsVarAxisLabel(var)+";Efficiency").c_str());
  return graph;
 
}


void HistFuncs::getRebinedBins(const TH1* hist,std::vector<float>& binLowEdges,float minInEachBin)
{
  binLowEdges.clear();

  int nrBins = hist->GetNbinsX();
  binLowEdges.push_back(hist->GetBinLowEdge(nrBins+1));
			
  float binTot=0.;

  // std::cout <<"nr bins "<<hist->GetNbinsX();
 

  for(int binNr=nrBins;binNr>1;binNr--){
    float nrEvts = hist->GetBinContent(binNr);
    
    binTot+=nrEvts;
    //   std::cout <<"binTot "<<binTot<<std::endl;
    if(binTot>minInEachBin){
      binTot=0.;
      binLowEdges.push_back(hist->GetBinLowEdge(binNr));
    }
    
  }
   
  binLowEdges.push_back(hist->GetBinLowEdge(1));
  std::sort(binLowEdges.begin(),binLowEdges.end());
}


double HistFuncs::getValueAtXLinear(const TH1* hist,double x)
{
  int binNr = AnaFuncs::getBinNr(hist,x);
  if(binNr<=1 || binNr>=hist->GetNbinsX()) return AnaFuncs::getBinContent(hist,x);
  float binMean = hist->GetBinLowEdge(binNr)+hist->GetBinWidth(binNr)/2.;
  if(x==binMean) return hist->GetBinContent(binNr); //the easy trival case
  else if(x<binMean){ //go towards lower bin
    if(binNr==1) return hist->GetBinContent(binNr); //for now we dont do this for the first bin
    float lowerBinMean = hist->GetBinLowEdge(binNr-1)+hist->GetBinWidth(binNr-1)/2.;
    float dMean = binMean-lowerBinMean;
    float dMeanX = binMean -x;
    return dMeanX/dMean*hist->GetBinContent(binNr-1) + (1-dMeanX/dMean)*hist->GetBinContent(binNr);
  }else{ //go higher
    if(binNr==hist->GetNbinsX()) return hist->GetBinContent(binNr); //for now we dont do this for the last bin
    float upperBinMean = hist->GetBinLowEdge(binNr+1)+hist->GetBinWidth(binNr+1)/2.;
    float dMean = fabs(binMean-upperBinMean);
    float dMeanX = fabs(binMean -x);
    return dMeanX/dMean*hist->GetBinContent(binNr+1) + (1-dMeanX/dMean)*hist->GetBinContent(binNr);
  }
}

void HistFuncs::normBinsByBinWidth(TH1* hist,float widthToNormTo)
{
  for(int binNr=1;binNr<=hist->GetNbinsX();binNr++){
    hist->SetBinError(binNr,hist->GetBinError(binNr)*widthToNormTo/hist->GetBinWidth(binNr));
    hist->SetBinContent(binNr,hist->GetBinContent(binNr)*widthToNormTo/hist->GetBinWidth(binNr));
 
  }
}

TH1* HistFuncs::makeDataMinusBkgPlot(const TH1* data,const TH1* bkg,float minBkgExpectInBin,bool ratio)
{
  std::vector<float> binLowEdges;
  getRebinedBins(bkg,binLowEdges,minBkgExpectInBin);
 

  for(size_t i=0;i<binLowEdges.size();i++) std::cout <<"bin low edge "<<binLowEdges[i]<<std::endl;
  
  TH1* resultHist = new TH1D("data","data",binLowEdges.size()-1,&binLowEdges[0]);
  resultHist->SetDirectory(0);

  for(size_t binNr=0;binNr<binLowEdges.size()-1;binNr++){
    double nrData,dataErr;
    double nrBkg,bkgErr;
    AnaFuncs::getHistIntegral(data,binLowEdges[binNr],binLowEdges[binNr+1],nrData,dataErr);
    AnaFuncs::getHistIntegral(bkg,binLowEdges[binNr],binLowEdges[binNr+1],nrBkg,bkgErr);
 
    std::cout <<" bin low edges "<<binLowEdges[binNr]<<" "<<binLowEdges[binNr+1]<<" nrData "<<nrData<<" nrBkg "<<nrBkg<<std::endl;
    if(ratio){
      if(nrBkg!=0){
	resultHist->SetBinContent(binNr+1,nrData/nrBkg -1);
	// float err = sqrt(nrData)/nrBkg;
	// double dataLow =  (N==0) ? 0  : (ROOT::Math::gamma_quantile(alpha/2,N,1.));
	// double dataHigh =  ROOT::Math::gamma_quantile_c(alpha/2,N+1,1) ;
	float err = sqrt(1/nrData+bkgErr*bkgErr/(nrBkg*nrBkg))*nrData/nrBkg;
	resultHist->SetBinError(binNr+1,err);
      }else{
	resultHist->SetBinContent(binNr+1,0);
	resultHist->SetBinError(binNr+1,0.00001);
      }
    }else{
       resultHist->SetBinContent(binNr+1,nrData-nrBkg);
       //float err = sqrt(nrData);
       float err = sqrt(nrData+bkgErr*bkgErr);
       resultHist->SetBinError(binNr+1,err);

    }
  
  }
  
  return resultHist;

}

TGraph* HistFuncs::makeDataMinusBkgPlotAsymErr(const TH1* data,const TH1* bkg,float minBkgExpectInBin)
{
  std::vector<float> binLowEdges;
  getRebinedBins(bkg,binLowEdges,minBkgExpectInBin);
 

  for(size_t i=0;i<binLowEdges.size();i++) std::cout <<"bin low edge "<<binLowEdges[i]<<std::endl;
  
  std::vector<double> xPoint,yPoint,xErrLow,xErrHigh,yErrLow,yErrHigh;

  for(size_t binNr=0;binNr<binLowEdges.size()-1;binNr++){
    double nrData,dataErr;
    double nrBkg,bkgErr;
    AnaFuncs::getHistIntegral(data,binLowEdges[binNr],binLowEdges[binNr+1],nrData,dataErr);
    AnaFuncs::getHistIntegral(bkg,binLowEdges[binNr],binLowEdges[binNr+1],nrBkg,bkgErr);
 
    const double alpha = (1 - 0.6827)/2;
    const double beta  = (1 - 0.6827)/2;
    double dataLowBound = 0.5*ROOT::Math::chisquared_quantile_c(1-alpha, 2*nrData);
    double dataHighBound = 0.5*ROOT::Math::chisquared_quantile_c(beta, 2*(nrData+1));
   
    auto dataBkgErrComb = [](double data,double dataErr,double bkg,double bkgErr)->double{
      if(bkg==0 || data==0) return 0.;
      else return std::sqrt(dataErr*dataErr/data/data + bkgErr*bkgErr/bkg/bkg)*data/bkg;
    };
    std::cout <<"bin lowEdge "<<binLowEdges[binNr]<<" bkg "<<nrBkg<<" bkgErr "<<bkgErr<<std::endl;
    double totYErrLow = dataBkgErrComb(nrData,nrData-dataLowBound,nrBkg,bkgErr);
    double totYErrHigh = dataBkgErrComb(nrData,dataHighBound-nrData,nrBkg,bkgErr);
 
    xPoint.push_back((binLowEdges[binNr+1]+binLowEdges[binNr])/2);
    yPoint.push_back((nrData-nrBkg)/nrBkg);
    xErrLow.push_back(xPoint.back()-binLowEdges[binNr]);
    xErrHigh.push_back(binLowEdges[binNr+1]-xPoint.back());
    yErrLow.push_back(totYErrLow);
    yErrHigh.push_back(totYErrHigh);
  }

  //  TGraphAsymmErrors* resultGraph = new TGraphAsymmErrors(xPoint.size(),&xPoint[0],&yPoint[0],0,0,&yErrLow[0],&yErrHigh[0]);
  TGraphAsymmErrors* resultGraph = new TGraphAsymmErrors(xPoint.size(),&xPoint[0],&yPoint[0],&xErrLow[0],&xErrHigh[0],&yErrLow[0],&yErrHigh[0]);
  if(!plotXErrOnAsymEff_){
    for(int pointNr=0;pointNr<resultGraph->GetN();pointNr++){
      resultGraph->SetPointEXhigh(pointNr,0);
      resultGraph->SetPointEXlow(pointNr,0);
    }
  }

  return resultGraph;

}



TH1* HistFuncs::makeCHist(const TH1* hist,bool intIsGreatThan)
{
  TH1* cHist = (TH1*) hist->Clone("cHist");
  cHist->SetDirectory(0);
  
  int maxBin = hist->GetNbinsX()+1;

  for(int binNr=0;binNr<=hist->GetNbinsX();binNr++){
    double binErr=0;
    float nrEntries = intIsGreatThan ? hist->IntegralAndError(binNr,maxBin,binErr) : hist->IntegralAndError(0,binNr,binErr);
    cHist->SetBinContent(binNr,nrEntries);
    cHist->SetBinError(binNr,binErr);
  }
  return cHist;
    
}

TGraph* HistFuncs::makeEffVsRejCurve(const TH1* sig,const TH1* bkg,bool sigGreaterThan,bool bkgGreaterThan)
{
  std::vector<float> sigPoints;
  std::vector<float> bkgPoints;

  // std::cout <<"about to make curve"<<std::endl;

  int nrBins = sig->GetNbinsX();

  float totSig = sig->Integral(0,sig->GetNbinsX()+1);
  float totBkg = bkg->Integral(0,bkg->GetNbinsX()+1);

  std::cout <<"totSig "<<totSig<<" tot bkg "<<totBkg<<std::endl;

  // std::cout <<"loop starting"<<std::endl;
  for(int binNr=1;binNr<=sig->GetNbinsX()+1;binNr++){
    float sigEff;
    if(sigGreaterThan) sigEff = sig->Integral(binNr,nrBins+1)/totSig;
    else sigEff = sig->Integral(0,binNr-1)/totSig;
    
    float bkgEff;
    if(bkgGreaterThan)  bkgEff = bkg->Integral(binNr,nrBins+1)/totBkg;
    else bkgEff = bkg->Integral(0,binNr-1)/totBkg;

    sigPoints.push_back(sigEff);
    bkgPoints.push_back(bkgEff);   
  }

  // std::cout <<"make curve"<<std::endl;

  TGraph* effVsRejCurve = new TGraph(sigPoints.size(),&bkgPoints[0],&sigPoints[0]);
  return effVsRejCurve;

}

void HistFuncs::print(const std::string& fileName,const std::string& canvasName,bool reduced)
{
  TCanvas* canvas = (TCanvas*) gROOT->FindObject(canvasName.c_str());
  //canvas->RedrawAxis();
  
  std::string outputName(fileName);
  
  std::string outputNameEps = outputName + ".eps";
  std::string outputNameGif = outputName + ".gif";
  std::string outputNameC = outputName + ".C";
  std::string outputNamePdf = outputName + ".pdf";
  std::string outputNamePng = outputName + ".png";
  if(!reduced){
    canvas->Print(outputNameEps.c_str());
    canvas->Print(outputNameGif.c_str());
    canvas->Print(outputNamePdf.c_str());
  }
  canvas->Print(outputNamePng.c_str());
  canvas->Print(outputNameC.c_str());
}


TGraph* HistFuncs::makeEffVsRejCurve2D(const TH2* sig,const TH2* bkg,bool sigGreaterThan,bool bkgGreaterThan)
{
  TH1* temp = sig->ProjectionY();
  TH1* sig1DHist = (TH1*) temp->Clone("temp");
  sig1DHist->SetDirectory(0);
  
  temp = bkg->ProjectionY("temp2");
  TH1* bkg1DHist = (TH1*) temp->Clone("bkg1D");
  bkg1DHist->SetDirectory(0);

  //std::cout <<"made hists"<<std::endl;

  TGraph* curve = makeEffVsRejCurve(sig,bkg,sigGreaterThan,bkgGreaterThan);
  delete bkg1DHist;
  delete sig1DHist;

  return curve;
}

TH1* HistFuncs::getProjection(const char* name,const char* filename,int startBin,int endBin)
{
  TFile* file = new TFile(filename,"READ");
  TH2* hist = (TH2*) file->Get(name);
  if(hist==NULL) std::cout <<"error "<<name<<" no found"<<std::endl;
  TH1* temp = hist->ProjectionY("tempProj",startBin,endBin);
  std::string newName(name);
  newName+="_1DProj";
  TH1* clone = (TH1*) temp->Clone(newName.c_str());
  clone->SetDirectory(0);
  delete file;
  return clone;
  
}

TGraph* HistFuncs::makeBestEffVsRejCurve(const TH2* sig,const TH2* bkg,bool sigGreaterThan,bool bkgGreaterThan)
{

 

  const int nrDecPlaces = 3;
  const double step  = 1./MathFuncs::power(10.,nrDecPlaces);
  std::vector<std::pair<float,float> > sigBkgEffs;
  for(int i=0;i<=MathFuncs::power(10,nrDecPlaces);i++){
    sigBkgEffs.push_back(std::make_pair(i*step,999.));
  }

  

  // std::cout <<"about to make curve"<<std::endl;

  int nrXBins = sig->GetNbinsX();
  int nrYBins = sig->GetNbinsY();

  float totSig = sig->Integral(0,sig->GetNbinsX()+1,0,sig->GetNbinsY()+1);
  float totBkg = bkg->Integral(0,bkg->GetNbinsX()+1,0,bkg->GetNbinsY()+1);

  // std::cout <<"loop starting"<<std::endl;
  for(int xBinNr=1;xBinNr<=sig->GetNbinsX()+1;xBinNr++){
    for(int yBinNr=1;yBinNr<=sig->GetNbinsY()+1;yBinNr++){
      float sigEff =  sigGreaterThan ? sig->Integral(xBinNr,nrXBins+1,yBinNr,nrYBins+1)/totSig : sig->Integral(0,xBinNr-1,0,yBinNr-1)/totSig;   
      float bkgEff = bkgGreaterThan ?  bkg->Integral(xBinNr,nrXBins+1,0,nrYBins+1)/totBkg : bkg->Integral(0,xBinNr-1,0,yBinNr-1)/totBkg;
 
      std::vector<std::pair<float,float> >::iterator effEntry = std::upper_bound(sigBkgEffs.begin(),sigBkgEffs.end(),std::make_pair(sigEff,-1),TempFuncs::PairComp<float,float,std::less<float> >());
    
      //if(sigEff!=0) std::cout <<"sig "<<sigEff<<std::endl;

      //bounds check
      if(effEntry==sigBkgEffs.end() && effEntry==sigBkgEffs.begin()){
       	LogErr <<" warning "<<sigEff<<" not found "<<std::endl;
	continue;
      }
      --effEntry; //to get into it so binLowEdge<= X < nextBinLowEdge
      
      if(sigEff<effEntry->first || sigEff>=(effEntry+1)->first){
	LogErr<<" eff "<<sigEff<<" bin lowEdge "<<effEntry->first<<" bin highEdge "<<(effEntry==sigBkgEffs.end()-1 ? 999 : (effEntry+1)->first)<<std::endl;
	continue;
      }
      if(effEntry->second>bkgEff) effEntry->second = bkgEff; //pick the smallest bkg eff for a given signal eff
    }
  }

  // std::cout <<"make curve"<<std::endl;

  std::vector<float> sigPoints;
  std::vector<float> bkgPoints;
  for(size_t i=0;i<sigBkgEffs.size();i++){
    if(sigBkgEffs[i].second<10){
      sigPoints.push_back(sigBkgEffs[i].first);
      bkgPoints.push_back(sigBkgEffs[i].second);
    }
  }

  TGraph* effVsRejCurve = new TGraph(sigPoints.size(),&bkgPoints[0],&sigPoints[0]);
  return effVsRejCurve;

}


TH1* HistFuncs::makeCutValueForFixedEffHist(TTree* tree,int nrBins,float xMin,float xMax,const std::string& vsVar,const std::string& sampleCuts,const std::string& var,float eff)
{
  TH1* hist = new TH1D("histTemp","pass",nrBins,xMin,xMax);
  hist->Sumw2();
  
  tree->SetEstimate(tree->GetEntries());  
  size_t nrPass = tree->Draw((var+":"+vsVar).c_str(),sampleCuts.c_str(),"goff");
  nrPass = tree->GetSelectedRows() % tree->GetEstimate();
  double* varArray = tree->GetVal(0);
  double* vsVarArray = tree->GetVal(1);
  
  //  std::cout <<"here "<<varArray<<" "<<vsVarArray<<std::endl;

  std::vector<std::vector<double> > varBinned(nrBins+2); //+2 for the under/over flow
  for(size_t entryNr=0;entryNr<nrPass;entryNr++){
    size_t binNr=AnaFuncs::getBinNr(nrBins,xMin,xMax,vsVarArray[entryNr]);
    varBinned[binNr].push_back(varArray[entryNr]);
  }
  for(size_t binNr=0;binNr<varBinned.size();binNr++){
    //std::cout <<"binNr "<<binNr<<" "<<varBinned.size()<<" "<<varBinned[binNr].size()<<std::endl;
    std::sort(varBinned[binNr].begin(),varBinned[binNr].end());
    if(varBinned[binNr].size()>10){
      double cutValForEff= MathFuncs::getWeightedListValue(varBinned[binNr],eff*varBinned[binNr].size());
      hist->SetBinContent(binNr,cutValForEff);
      float errUpNr = std::min(eff+0.005,1.)*varBinned[binNr].size();
      float errDownNr =  std::max(eff-0.005,0.)*varBinned[binNr].size();
      float errUp = MathFuncs::getWeightedListValue(varBinned[binNr],errUpNr) - cutValForEff;
      float errDown = cutValForEff - MathFuncs::getWeightedListValue(varBinned[binNr],errDownNr);
      float err = std::max(errUp,errDown);
      //std::cout <<"err "<<err<<" up "<<errUp<<" down "<<errDown<<" errDownNr "<<errDownNr<<" errUpNr "<<errUpNr<<std::endl;
      hist->SetBinError(binNr,err);
    }else{
      hist->SetBinContent(binNr,-1);
      hist->SetBinError(binNr,0.0001);
    }
  }
     
  hist->SetDirectory(0);
  
  return hist;

}

std::vector<TH1*> HistFuncs::makeCutValueForFixedEffHist(TTree* tree,int nrBins,float xMin,float xMax,const std::string& vsVar,const std::string& sampleCuts,const std::string& var,const std::vector<float>& effs)
{
  std::vector<TH1*> hists;
  for(size_t i=0;i<effs.size();i++){
    hists.push_back(new TH1D("histTemp","pass",nrBins,xMin,xMax));
    hists.back()->Sumw2();
    hists.back()->SetDirectory(0);
  }
  
  tree->SetEstimate(tree->GetEntries());  
  size_t nrPass = tree->Draw((var+":"+vsVar).c_str(),sampleCuts.c_str(),"goff");
  nrPass = tree->GetSelectedRows() % tree->GetEstimate();
  double* varArray = tree->GetVal(0);
  double* vsVarArray = tree->GetVal(1);
  
  std::vector<std::vector<double> > varBinned(nrBins+2); //+2 for the under/over flow
  for(size_t entryNr=0;entryNr<nrPass;entryNr++){
    size_t binNr=AnaFuncs::getBinNr(nrBins,xMin,xMax,vsVarArray[entryNr]);
    varBinned[binNr].push_back(varArray[entryNr]);
  }
  for(size_t binNr=0;binNr<varBinned.size();binNr++){
    std::sort(varBinned[binNr].begin(),varBinned[binNr].end());
    for(size_t effNr=0;effNr<effs.size();effNr++){
      TH1* hist = hists[effNr];	
      const float eff = effs[effNr];
      
      if(varBinned[binNr].size()>10){

	
	double cutValForEff= MathFuncs::getWeightedListValue(varBinned[binNr],eff*varBinned[binNr].size());
	hist->SetBinContent(binNr,cutValForEff);
	float errUpNr = std::min(eff+0.005,1.)*varBinned[binNr].size();
	float errDownNr =  std::max(eff-0.005,0.)*varBinned[binNr].size();
	float errUp = MathFuncs::getWeightedListValue(varBinned[binNr],errUpNr) - cutValForEff;
	float errDown = cutValForEff - MathFuncs::getWeightedListValue(varBinned[binNr],errDownNr);
	float err = std::max(errUp,errDown);
	//std::cout <<"err "<<err<<" up "<<errUp<<" down "<<errDown<<" errDownNr "<<errDownNr<<" errUpNr "<<errUpNr<<std::endl;
	hist->SetBinError(binNr,err);
     
      }else{
	hist->SetBinContent(binNr,-1);
	hist->SetBinError(binNr,0.0001);
      }
    }//end eff loop
  }
 
  
  return hists;

}

TH1* HistFuncs::makeBkgEffForFixedSigEffHist(TTree* sigTree,TTree* bkgTree,int nrBins,float xMin,float xMax,const std::string& vsVarSig,const std::string& vsVarBkg,const std::string& sampleCutsSig,const std::string& sampleCutsBkg,const std::string& var,float eff)
{
  TH1* hist = new TH1D("histTemp","pass",nrBins,xMin,xMax);
  hist->Sumw2();
  hist->SetDirectory(0);
  
  TH1* cutValueHist = makeCutValueForFixedEffHist(sigTree,nrBins,xMin,xMax,vsVarSig,sampleCutsSig,var,eff);
  
  
  bkgTree->SetEstimate(bkgTree->GetEntries());  
  size_t nrPass = bkgTree->Draw((var+":"+vsVarBkg).c_str(),sampleCutsBkg.c_str(),"goff");
  nrPass = bkgTree->GetSelectedRows() % bkgTree->GetEstimate();
  double* varArray = bkgTree->GetVal(0);
  double* vsVarArray = bkgTree->GetVal(1);
  double* weightArray = bkgTree->GetW();
  //  std::vector<std::vector<double> > varBinned(nrBins+2); //+2 for the under/over flow
  std::vector<std::vector<std::pair<double,double > > > varBinned(nrBins+2); //+2 for the under/over flow
  for(size_t entryNr=0;entryNr<nrPass;entryNr++){
    size_t binNr=AnaFuncs::getBinNr(nrBins,xMin,xMax,vsVarArray[entryNr]);
    varBinned[binNr].push_back(std::pair<double,double>(varArray[entryNr],weightArray[entryNr]));
  }
  for(size_t binNr=0;binNr<varBinned.size();binNr++){
    //    std::cout <<"binNr "<<binNr<<" "<<varBinned.size()<<" "<<varBinned[binNr].size()<<std::endl;
    std::sort(varBinned[binNr].begin(),varBinned[binNr].end());
    float cutValue = cutValueHist->GetBinContent(binNr);
    
   
    float nrPass=0;
    float nrPassW2=0;
    float nrTot=0;
    float nrTotW2=0;

 
    for(size_t i=0;i<varBinned[binNr].size();i++){
      if(varBinned[binNr][i].first<=cutValue){
	nrPass+=varBinned[binNr][i].second;
	nrPassW2+=varBinned[binNr][i].second*varBinned[binNr][i].second;
      }
      nrTot+=varBinned[binNr][i].second;
      nrTotW2+=varBinned[binNr][i].second*varBinned[binNr][i].second;
    }

    std::cout <<"bin "<<binNr<<" cutValue "<<cutValue<<" nrPass "<<nrPass<<" nrTot "<<nrTot<<std::endl;

    if(nrTot!=0){
      std::pair<float,float> effAndErr = MathFuncs::calEffAndErr(nrPass,nrPassW2,nrTot,nrTotW2);
      
      hist->SetBinContent(binNr,effAndErr.first);
      hist->SetBinError(binNr,effAndErr.second);
    }else{
        
      hist->SetBinContent(binNr,0);
      hist->SetBinError(binNr,0.001);
    }
  }
     
  hist->SetDirectory(0);
  
  return hist;

}
std::pair<float,float> HistFuncs::calTotBkgEffForFixedSigEffHist(TTree* sigTree,TTree* bkgTree,int nrBins,float xMin,float xMax,const std::string& vsVarSig,const std::string& vsVarBkg,const std::string& sampleCutsSig,const std::string& sampleCutsBkg,const std::string& var,float eff)
{
  
  TH1* cutValueHist = makeCutValueForFixedEffHist(sigTree,nrBins,xMin,xMax,vsVarSig,sampleCutsSig,var,eff);
  
  
  bkgTree->SetEstimate(bkgTree->GetEntries());  
  size_t nrPassSample = bkgTree->Draw((var+":"+vsVarBkg).c_str(),sampleCutsBkg.c_str(),"goff");
  nrPassSample = bkgTree->GetSelectedRows() % bkgTree->GetEstimate();
  double* varArray = bkgTree->GetVal(0);
  double* vsVarArray = bkgTree->GetVal(1);
   std::vector<std::vector<double> > varBinned(nrBins+2); //+2 for the under/over flow
  for(size_t entryNr=0;entryNr<nrPassSample;entryNr++){
    size_t binNr=AnaFuncs::getBinNr(nrBins,xMin,xMax,vsVarArray[entryNr]);
    varBinned[binNr].push_back(varArray[entryNr]);
  }
  int nrTot=0;
  int nrPass=0;
  for(size_t binNr=0;binNr<varBinned.size();binNr++){
    // std::cout <<"binNr "<<binNr<<" "<<varBinned.size()<<" "<<varBinned[binNr].size()<<std::endl;
    std::sort(varBinned[binNr].begin(),varBinned[binNr].end());
    float cutValue = cutValueHist->GetBinContent(binNr);
    
    nrTot+= varBinned[binNr].size();
    for(size_t i=0;i<varBinned[binNr].size();i++){
      if(varBinned[binNr][i]<=cutValue) nrPass++;
    }

    std::cout <<"bin "<<binNr<<" cutValue "<<cutValue<<" nrPass "<<nrPass<<" nrTot "<<nrTot<<std::endl;

   
  }
  if(nrTot!=0) return MathFuncs::calEffAndErr(nrPass,nrTot);
  else return std::pair<float,float>(0,0);


}

TGraph* HistFuncs::makeBkgEffForFixedSigEffHist(TTree* sigTree,TTree* bkgTree,int nrBins,float xMin,float xMax,const std::string& vsVarSig,const std::string& vsVarBkg,const std::string& sampleCutsSig,const std::string& sampleCutsBkg,const std::string& var,const std::vector<float>& effs)
{
  std::vector<TH1*> cutValueHist = makeCutValueForFixedEffHist(sigTree,nrBins,xMin,xMax,vsVarSig,sampleCutsSig,var,effs);
  
  std::vector<float> bkgEffs(effs.size(),0.);
  std::vector<float> bkgEffErrs(effs.size(),0.);
  bkgTree->SetEstimate(bkgTree->GetEntries());  
  size_t nrPass = bkgTree->Draw((var+":"+vsVarBkg).c_str(),sampleCutsBkg.c_str(),"goff");
  nrPass = bkgTree->GetSelectedRows() % bkgTree->GetEstimate();
  double* varArray = bkgTree->GetVal(0);
  double* vsVarArray = bkgTree->GetVal(1);
  double* weightArray = bkgTree->GetW();
  std::vector<std::vector<std::pair<double,double > > > varBinned(nrBins+2); //+2 for the under/over flow
  for(size_t entryNr=0;entryNr<nrPass;entryNr++){
    size_t binNr=AnaFuncs::getBinNr(nrBins,xMin,xMax,vsVarArray[entryNr]);
    varBinned[binNr].push_back(std::pair<double,double>(varArray[entryNr],weightArray[entryNr]));
  }

  for(size_t effNr=0;effNr<effs.size();effNr++){
    
    float nrPass=0;
    float nrPassW2=0;
    float nrTot=0;
    float nrTotW2=0;
    for(size_t binNr=1;binNr+1<varBinned.size();binNr++){
      std::sort(varBinned[binNr].begin(),varBinned[binNr].end());
      float cutValue = cutValueHist[effNr]->GetBinContent(binNr);
     
      for(size_t i=0;i<varBinned[binNr].size();i++){ //loop over bkg data values
	if(varBinned[binNr][i].first<=cutValue){
	  nrPass+=varBinned[binNr][i].second;
	  nrPassW2+=varBinned[binNr][i].second*varBinned[binNr][i].second;
	}
	nrTot+=varBinned[binNr][i].second;
	nrTotW2+=varBinned[binNr][i].second*varBinned[binNr][i].second;
      }
    }
    
    std::cout <<"eff "<<effs[effNr]<<" nrPass "<<nrPass<<" nrTot "<<nrTot<<std::endl;

    if(nrTot!=0){
      std::pair<float,float> effAndErr = MathFuncs::calEffAndErr(nrPass,nrPassW2,nrTot,nrTotW2);
      bkgEffs[effNr]=effAndErr.first;
      bkgEffErrs[effNr]=effAndErr.second;
    }
  }
     
  TGraphErrors* graph = new TGraphErrors(effs.size(),&effs[0],&bkgEffs[0],0,&bkgEffErrs[0]);
  
  return graph;

}
			   

TH1* HistFuncs::compTwoVars(TTree* tree,int nrBins,float xmin,float xmax,const std::string& var1,const std::string& var2,const std::string& cuts,const std::string& var1LegName,const std::string& var2LegName)
{

  TH1* var1Hist = makeHist(tree,nrBins,xmin,xmax,var1,cuts);
  TH1* var2Hist = makeHist(tree,nrBins,xmin,xmax,var2,cuts); 
  
  AnaFuncs::setHistAttributes(var1Hist,4,1,8,4);
  AnaFuncs::setHistAttributes(var2Hist,2,1,4,2);

  var1Hist->Draw("EP");
  var2Hist->Draw("SAME EP");

  TLegend* leg = new TLegend(0.3,0.4,0.5,0.6);
  leg->AddEntry(var1Hist,var1LegName.c_str(),"LP");
  leg->AddEntry(var2Hist,var2LegName.c_str(),"LP");
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->Draw();

  return var1Hist;
}

TH1* HistFuncs::makeHist(TTree* tree,int nrBins,float xmin,float xmax,const std::string& var,const std::string& cuts)
{
  TH1* hist = new TH1D("var1Hist","temp",nrBins,xmin,xmax);
  hist->Sumw2();

  tree->Draw((var+">>var1Hist").c_str(),cuts.c_str(),"goff");
  hist->SetDirectory(0);
  hist->SetTitle(("; "+vsVarAxisLabel(var)+";").c_str());
  return hist;
}
TH1* HistFuncs::makeHist(TTree* tree,const HistFuncs::Bins1D& bins,const std::string& var,const std::string& cuts)
{
  TH1* hist=bins.makeHist("var1Hist");
  hist->Sumw2();

  tree->Draw((var+">>var1Hist").c_str(),cuts.c_str(),"goff");
  hist->SetDirectory(0);
  hist->SetTitle(("; "+vsVarAxisLabel(var)+";").c_str());
  return hist;
}

TH2* HistFuncs::makeHist(TTree* tree,int nrBinsX,float xmin,float xmax,int nrBinsY,float ymin,float ymax,const std::string& var,const std::string& cuts)
{
  TH2* hist = new TH2D("var1Hist","temp",nrBinsX,xmin,xmax,nrBinsY,ymin,ymax);
  hist->Sumw2();
  tree->Draw((var+">>var1Hist").c_str(),cuts.c_str(),"goff");
  hist->SetDirectory(0);
  return hist;
}

TH2* HistFuncs::makeHist(TTree* tree,const std::vector<double>& binLowEdges,int nrBinsY,float ymin,float ymax,const std::string& var,const std::string& cuts)
{
  TH2* hist = new TH2D("var1Hist","temp",binLowEdges.size()-1,binLowEdges.data(),nrBinsY,ymin,ymax);
  hist->Sumw2();
  tree->Draw((var+">>var1Hist").c_str(),cuts.c_str(),"goff");
  hist->SetDirectory(0);
  return hist;
}

TH1* HistFuncs::makeHistBothEles(TTree* tree,int nrBins,float xmin,float xmax,const std::string& var,const std::string& cuts)
{  
  std::string var1 = boost::algorithm::replace_all_copy(var,"{1}","1");
  var1 = boost::algorithm::replace_all_copy(var1,"{2}","2");
  std::string cuts1 = boost::algorithm::replace_all_copy(cuts,"{1}","1");
  cuts1 = boost::algorithm::replace_all_copy(cuts1,"{2}","2");


  std::string var2 = boost::algorithm::replace_all_copy(var,"{1}","2");
  var2 = boost::algorithm::replace_all_copy(var2,"{2}","1");
  std::string cuts2 = boost::algorithm::replace_all_copy(cuts,"{1}","2");
  cuts2 = boost::algorithm::replace_all_copy(cuts2,"{2}","1");
  std::cout <<"1: "<<var1<<" : "<<cuts1<<std::endl;
  std::cout <<"2: "<<var2<<" : "<<cuts2<<std::endl;
  
  TH1* totHist1 = HistFuncs::makeHist(tree,nrBins,xmin,xmax,var1,cuts1);
  TH1* totHist2 = HistFuncs::makeHist(tree,nrBins,xmin,xmax,var2,cuts2);
   
  totHist1->Add(totHist2);
 
  delete totHist2; 

  return totHist1;
}

TH2* HistFuncs::makeHistBothEles(TTree* tree,int nrBinsX,float xmin,float xmax,int nrBinsY,float ymin,float ymax,const std::string& var,const std::string& cuts)
{  
  std::string var1 = boost::algorithm::replace_all_copy(var,"{1}","1");
  var1 = boost::algorithm::replace_all_copy(var1,"{2}","2");
  std::string cuts1 = boost::algorithm::replace_all_copy(cuts,"{1}","1");
  cuts1 = boost::algorithm::replace_all_copy(cuts1,"{2}","2");


  std::string var2 = boost::algorithm::replace_all_copy(var,"{1}","2");
  var2 = boost::algorithm::replace_all_copy(var2,"{2}","1");
  std::string cuts2 = boost::algorithm::replace_all_copy(cuts,"{1}","2");
  cuts2 = boost::algorithm::replace_all_copy(cuts2,"{2}","1");
  std::cout <<"1: "<<var1<<" : "<<cuts1<<std::endl;
  std::cout <<"2: "<<var2<<" : "<<cuts2<<std::endl;
  
  TH2* totHist1 = HistFuncs::makeHist(tree,nrBinsX,xmin,xmax,nrBinsY,ymin,ymax,var1,cuts1);
  TH2* totHist2 = HistFuncs::makeHist(tree,nrBinsX,xmin,xmax,nrBinsY,ymin,ymax,var2,cuts2);
   
  totHist1->Add(totHist2);
 
  delete totHist2; 

  return totHist1;
}

TProfile* HistFuncs::makeHistBothElesProf(TTree* tree,int nrBinsX,float xmin,float xmax,const std::string& var,const std::string& cuts)
{  
  std::string var1 = boost::algorithm::replace_all_copy(var,"{1}","1");
  var1 = boost::algorithm::replace_all_copy(var1,"{2}","2");
  std::string cuts1 = boost::algorithm::replace_all_copy(cuts,"{1}","1");
  cuts1 = boost::algorithm::replace_all_copy(cuts1,"{2}","2");


  std::string var2 = boost::algorithm::replace_all_copy(var,"{1}","2");
  var2 = boost::algorithm::replace_all_copy(var2,"{2}","1");
  std::string cuts2 = boost::algorithm::replace_all_copy(cuts,"{1}","2");
  cuts2 = boost::algorithm::replace_all_copy(cuts2,"{2}","1");
  std::cout <<"1: "<<var1<<" : "<<cuts1<<std::endl;
  std::cout <<"2: "<<var2<<" : "<<cuts2<<std::endl;
  
  TProfile* hist = new TProfile("var1Hist","temp",nrBinsX,xmin,xmax);
  hist->Sumw2();
  tree->Draw((var1+">>var1Hist").c_str(),cuts1.c_str(),"goff");
  tree->Draw((var2+">>+var1Hist").c_str(),cuts2.c_str(),"goff");
  hist->SetDirectory(0);
  return hist;
}



TH1* HistFuncs::compTwoVars(TTree* tree1,TTree* tree2,int nrBins,float xmin,float xmax,const std::string& var1,const std::string& var2,const std::string& cuts1,const std::string& cuts2,const std::string& var1LegName,const std::string& var2LegName,bool norm)
{
  TH1* var1Hist = makeHist(tree1,nrBins,xmin,xmax,var1,cuts1);
  TH1* var2Hist = makeHist(tree2,nrBins,xmin,xmax,var2,cuts2); 


  if(norm && var1Hist->Integral(0,nrBins+1)!=0) var1Hist->Scale(var2Hist->Integral(0,nrBins+1)/var1Hist->Integral(0,nrBins+1));
  
  AnaFuncs::setHistAttributes(var1Hist,4,1,8,4);
  AnaFuncs::setHistAttributes(var2Hist,2,1,4,2);

  var1Hist->Draw("EP");
  var2Hist->Draw("SAME EP");

  TLegend* leg = new TLegend(0.3,0.4,0.5,0.6);
  leg->AddEntry(var1Hist,var1LegName.c_str(),"LP");
  leg->AddEntry(var2Hist,var2LegName.c_str(),"LP");
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->Draw();

  return var1Hist;
}

TH1* HistFuncs::compTwoCuts(TTree* tree,int nrBins,float xmin,float xmax,const std::string& var,const std::string& baseCuts,const std::string& cuts1,const std::string& cuts2,const std::string& var1LegName,const std::string& var2LegName)
{
  TH1* var1Hist = new TH1D("var1Hist","temp",nrBins,xmin,xmax);
  var1Hist->Sumw2();
  tree->Draw((var+">>var1Hist").c_str(),(baseCuts+cuts1).c_str());
  TH1* var2Hist = new TH1D("var2Hist","temp",nrBins,xmin,xmax);
  var2Hist->Sumw2();
  tree->Draw((var+">>var2Hist").c_str(),(baseCuts+cuts2).c_str());
  
  var1Hist->SetDirectory(0);
  var2Hist->SetDirectory(0);
  
  float nrEntries1= var1Hist->Integral();
  float nrEntries2= var2Hist->Integral();
  float nrEntriesMin = std::min(nrEntries1,nrEntries2);
  var1Hist->Scale(nrEntriesMin/var1Hist->Integral());
  var2Hist->Scale(nrEntriesMin/var2Hist->Integral());
  

  AnaFuncs::setHistAttributes(var1Hist,4,1,8,4);
  AnaFuncs::setHistAttributes(var2Hist,2,1,4,2);

  var1Hist->Draw("HISTE");
  var2Hist->Draw("SAME HISTE");

  TLegend* leg = new TLegend(0.3,0.4,0.5,0.6);
  leg->AddEntry(var1Hist,var1LegName.c_str(),"LP");
  leg->AddEntry(var2Hist,var2LegName.c_str(),"LP");
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->Draw();

  return var1Hist;
}

TH1* HistFuncs::compTwoEffs(TTree* tree,int nrBins,float xmin,float xmax,const std::string& var,const std::string& baseCuts,const std::string& cuts1,const std::string& cuts2,const std::string& var1LegName,const std::string& var2LegName,const std::string& weight)
{
  TH1* var1Hist = HistFuncs::makeEffHistFromTree(tree,nrBins,xmin,xmax,var,baseCuts,cuts1,weight);
  TH1* var2Hist = HistFuncs::makeEffHistFromTree(tree,nrBins,xmin,xmax,var,baseCuts,cuts2,weight);
  


  AnaFuncs::setHistAttributes(var1Hist,4,1,8,4);
  AnaFuncs::setHistAttributes(var2Hist,2,1,4,2);

  var1Hist->Draw("EP");
  var2Hist->Draw("SAME EP");

  TLegend* leg = new TLegend(0.3,0.4,0.5,0.6);
  leg->AddEntry(var1Hist,var1LegName.c_str(),"LP");
  leg->AddEntry(var2Hist,var2LegName.c_str(),"LP");
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->Draw();

  return var1Hist;
}


TH1* HistFuncs::makeCutValueForFixedEffHist(TH2* hist,float eff)
{
  TAxis* xAxis = hist->GetXaxis();
  std::vector<float> binLowEdges;
  for(int binNr=1;binNr<=xAxis->GetNbins()+1;binNr++) binLowEdges.push_back(xAxis->GetBinLowEdge(binNr));
  TH1* outHist = new TH1D("cutValHist","hist",binLowEdges.size()-1,&binLowEdges[0]); 
  outHist->SetDirectory(0);
  for(int xBinNr=1;xBinNr<=hist->GetNbinsX();xBinNr++){
    std::vector<int> vec;
    for(int yBinNr=1;yBinNr<=hist->GetNbinsY();yBinNr++){
      TAxis* yAxis = hist->GetYaxis();
      int binContent = hist->GetBinContent(xBinNr,yBinNr);
      for(int i=0;i<binContent;i++) vec.push_back(yAxis->GetBinLowEdge(yBinNr));
    }
    std::sort(vec.begin(),vec.end());
    std::pair<float,int> cutVal = MathFuncs::estQuantile(vec,eff);
    outHist->SetBinContent(xBinNr,cutVal.first);
    outHist->SetBinError(xBinNr,cutVal.second);
    if(cutVal.second<0.5) outHist->SetBinError(xBinNr,0.5);
  }
  return outHist;
}

void HistFuncs::fillVecFromTree(TTree* tree,const std::string& var,const std::string& cuts,std::vector<float>& vec)
{
  
  tree->SetEstimate(tree->GetEntries()+1);
  tree->Draw(var.c_str(),cuts.c_str(),"goff");
  size_t nrPass = tree->GetSelectedRows() % tree->GetEstimate()+1;
  double* varArray = tree->GetVal(0);
  vec.clear();
  vec.reserve(nrPass);
  for(size_t entryNr=0;entryNr<nrPass;entryNr++){
    vec.push_back(varArray[entryNr]);
  }
}

TH1* HistFuncs::plotTwoHists(TH1* hist1,TH1* hist2,const std::string& leg1Str,const std::string& leg2Str)
{
  AnaFuncs::setHistAttributes(hist1,4,1,8,4);
  AnaFuncs::setHistAttributes(hist2,2,1,4,2);
  
  hist1->Draw();
  hist2->Draw("SAME EP");


  TLegend* leg = new TLegend(0.138393,0.176573,0.417411,0.34965);
  leg->SetBorderSize(0);
  leg->AddEntry(hist1,leg1Str.c_str(),"LP");
  leg->AddEntry(hist2,leg2Str.c_str(),"LP");
  leg->Draw();
  
  return hist1;
}

TGraph* HistFuncs::plotTwoHists(TGraph* hist1,TGraph* hist2,const std::string& leg1Str,const std::string& leg2Str)
{
  AnaFuncs::setHistAttributes(hist1,4,1,8,4);
  AnaFuncs::setHistAttributes(hist2,2,1,4,2);
  
  hist1->Draw("AP");
  hist2->Draw("P");

  TLegend* leg = new TLegend(0.138393,0.176573,0.417411,0.34965);
  leg->SetBorderSize(0);
  leg->AddEntry(hist1,leg1Str.c_str(),"LP");
  leg->AddEntry(hist2,leg2Str.c_str(),"LP");
  leg->Draw();
  
  return hist1;
}

void HistFuncs::setStdAxisTextSize(TH1* hist)
{
  auto setTextSize = [](TAxis* axis){
    axis->SetTitleSize(0.045);
    axis->SetLabelSize(0.045);
  };
  hist->GetYaxis()->SetTitleOffset(1.5);
  setTextSize(hist->GetXaxis());
  setTextSize(hist->GetYaxis());

}


TH1* HistFuncs::makeRatio(TH1* numerHist,TH1* denomHist,const std::string& option)
{
  TH1* hist = static_cast<TH1*>(numerHist->Clone("ratioHist"));
  hist->Divide(numerHist,denomHist,1,1,option.c_str());
  return hist;
}

TH1* HistFuncs::plotWithRatio(HistFuncs::HistOpts numer,HistFuncs::HistOpts denom,
			      const std::string& divOpt)
{
  TCanvas* c1 =new TCanvas("c1", "c1",900,750);
  c1->cd();
  TPad* spectrumPad = new TPad("spectrumPad", "newpad",0.01,0.30,0.99,0.99);
  spectrumPad->Draw(); 
  spectrumPad->cd();
  std::string xAxisLabel = denom.hist->GetXaxis()->GetTitle();
  denom.hist->GetXaxis()->SetTitle();
  denom.hist->Draw(denom.histDrawOpt.c_str());
  numer.hist->Draw((numer.histDrawOpt+" SAME").c_str());
  TLegend* leg = makeLegend<0>({{numer.hist,numer.legEntry},{denom.hist,denom.legEntry}});
  leg->Draw();
  
  c1->cd();
  TPad* ratioPad = new TPad("ratioPad", "newpad",0.01,0.01,0.99,0.33);
  ratioPad->Draw();
  ratioPad->cd();
  ratioPad->SetTopMargin(0.05);
  ratioPad->SetBottomMargin(0.3);
  //    ratioPad->SetRightMargin(0.1);
  ratioPad->SetFillStyle(0);
  
  TH1* ratioHist = makeRatio(numer.hist,denom.hist,divOpt);
  AnaFuncs::setHistAttributes(ratioHist,1,1,8,1);
  ratioHist->SetTitle("");
  //  ratioHist->GetXaxis()->SetLabelSize(ratioHist->GetXaxis()->GetLabelSize()*(0.99-0.33)/0.33);
  ratioHist->GetXaxis()->SetLabelSize(0.1);
  ratioHist->GetXaxis()->SetTitleSize(0.1);
  ratioHist->GetXaxis()->SetTitle(xAxisLabel.c_str());
  ratioHist->GetYaxis()->SetLabelSize(0.1);
  ratioHist->GetYaxis()->SetTitleSize(0.1);
  ratioHist->GetYaxis()->SetTitleOffset(0.65); 
  ratioHist->GetYaxis()->SetTitle("ratio");   
  
  ratioHist->Draw();
  spectrumPad->cd();
  return denom.hist;
}

TH1* HistFuncs::plotWithRatio(HistFuncs::HistOpts numer1,HistFuncs::HistOpts numer2,HistFuncs::HistOpts denom,
			      const std::string& divOpt)
{
  TCanvas* c1 =new TCanvas("c1", "c1",900,750);
  c1->cd();
  TPad* spectrumPad = new TPad("spectrumPad", "newpad",0.01,0.30,0.99,0.99);
  spectrumPad->Draw(); 
  spectrumPad->cd();
  std::string xAxisLabel = denom.hist->GetXaxis()->GetTitle();
  denom.hist->GetXaxis()->SetTitle();
  denom.hist->Draw(denom.histDrawOpt.c_str());
  numer1.hist->Draw((numer1.histDrawOpt+" SAME").c_str());
  numer2.hist->Draw((numer2.histDrawOpt+" SAME").c_str());
  TLegend* leg = makeLegend<2>({{numer1.hist,numer1.legEntry},{numer2.hist,numer2.legEntry},{denom.hist,denom.legEntry}});
  leg->Draw();
  
  c1->cd();
  TPad* ratioPad = new TPad("ratioPad", "newpad",0.01,0.01,0.99,0.33);
  ratioPad->Draw();
  ratioPad->cd();
  ratioPad->SetTopMargin(0.05);
  ratioPad->SetBottomMargin(0.3);
  //    ratioPad->SetRightMargin(0.1);
  ratioPad->SetFillStyle(0);
  
  TH1* ratio1Hist = makeRatio(numer1.hist,denom.hist,divOpt);
  auto copyHistStyle =[](TH1* refHist,TH1* newHist){
    AnaFuncs::setHistAttributes(newHist,refHist->GetLineColor(),refHist->GetLineWidth(),refHist->GetMarkerStyle(),refHist->GetLineColor());
  };
  copyHistStyle(numer1.hist,ratio1Hist);
  ratio1Hist->SetTitle("");
  //  ratioHist->GetXaxis()->SetLabelSize(ratioHist->GetXaxis()->GetLabelSize()*(0.99-0.33)/0.33);
  ratio1Hist->GetXaxis()->SetLabelSize(0.1);
  ratio1Hist->GetXaxis()->SetTitleSize(0.1);
  ratio1Hist->GetXaxis()->SetTitle(xAxisLabel.c_str());
  ratio1Hist->GetYaxis()->SetLabelSize(0.1);
  ratio1Hist->GetYaxis()->SetTitleSize(0.1);
  ratio1Hist->GetYaxis()->SetTitleOffset(0.65); 
  ratio1Hist->GetYaxis()->SetTitle("ratio");   
  
  ratio1Hist->Draw();

  TH1* ratio2Hist = makeRatio(numer2.hist,denom.hist,divOpt);

  copyHistStyle(numer2.hist,ratio2Hist);
  ratio2Hist->SetTitle("");
  ratio2Hist->Draw("SAME");

  spectrumPad->cd();
  return denom.hist;
}
TH1* HistFuncs::plotWithRatio(std::vector<HistFuncs::HistOpts> numers,HistFuncs::HistOpts denom,
			      const std::string& divOpt)
{
  TCanvas* c1 =new TCanvas("c1", "c1",900,750);
  c1->cd();
  TPad* spectrumPad = new TPad("spectrumPad", "newpad",0.01,0.30,0.99,0.99);
  spectrumPad->Draw(); 
  spectrumPad->cd();
  std::string xAxisLabel = denom.hist->GetXaxis()->GetTitle();
  denom.hist->GetXaxis()->SetTitle();
  denom.hist->Draw(denom.histDrawOpt.c_str());
  std::vector<std::pair<TH1*,std::string> >legEntryVec;
  for(auto& numer : numers){
    numer.hist->Draw((numer.histDrawOpt+" SAME").c_str());
    legEntryVec.push_back({numer.hist,numer.legEntry});
  }
  legEntryVec.push_back({denom.hist,denom.legEntry});
  TLegend* leg = makeLegend<2>(legEntryVec);
  leg->Draw();
  
  c1->cd();
  TPad* ratioPad = new TPad("ratioPad", "newpad",0.01,0.01,0.99,0.33);
  ratioPad->Draw();
  ratioPad->cd();
  ratioPad->SetTopMargin(0.05);
  ratioPad->SetBottomMargin(0.3);
  //    ratioPad->SetRightMargin(0.1);
  ratioPad->SetFillStyle(0);
  
  bool firstHist=true;
  for(auto& numer : numers){
    TH1* ratioHist = makeRatio(numer.hist,denom.hist,divOpt);
    auto copyHistStyle =[](TH1* refHist,TH1* newHist){
      AnaFuncs::setHistAttributes(newHist,refHist->GetLineColor(),refHist->GetLineWidth(),refHist->GetMarkerStyle(),refHist->GetLineColor());
    };
    copyHistStyle(numer.hist,ratioHist);
    ratioHist->SetTitle("");
    //  ratioHist->GetXaxis()->SetLabelSize(ratioHist->GetXaxis()->GetLabelSize()*(0.99-0.33)/0.33);
    ratioHist->GetXaxis()->SetLabelSize(0.1);
    ratioHist->GetXaxis()->SetTitleSize(0.1);
    ratioHist->GetXaxis()->SetTitle(xAxisLabel.c_str());
    ratioHist->GetYaxis()->SetLabelSize(0.1);
    ratioHist->GetYaxis()->SetTitleSize(0.1);
    ratioHist->GetYaxis()->SetTitleOffset(0.65); 
    ratioHist->GetYaxis()->SetTitle("ratio");   
   
    if(firstHist){
      ratioHist->Draw();
      firstHist=false;
    }else ratioHist->Draw("SAME");
    
  }
  spectrumPad->cd();
  return denom.hist;
}
TH1* HistFuncs::plotWithRatio(HistFuncs::HistOpts numer,std::vector<HistFuncs::HistOpts> denoms,
			      const std::string& divOpt)
{
  TCanvas* c1 =new TCanvas("c1", "c1",900,750);
  c1->cd();
  TPad* spectrumPad = new TPad("spectrumPad", "newpad",0.01,0.30,0.99,0.99);
  spectrumPad->Draw(); 
  spectrumPad->cd();
  std::vector<std::pair<TH1*,std::string> >legEntryVec;
  legEntryVec.push_back({numer.hist,numer.legEntry});
  bool firstHist=true;
  std::string xAxisLabel;
  for(auto& denom : denoms){
    if(firstHist){
      denom.hist->Draw((denom.histDrawOpt).c_str());
      xAxisLabel = denom.hist->GetXaxis()->GetTitle();
      firstHist=false;
    }else denom.hist->Draw((denom.histDrawOpt+" SAME").c_str());
    denom.hist->GetXaxis()->SetTitle();
    legEntryVec.push_back({denom.hist,denom.legEntry});
  }
  numer.hist->Draw((numer.histDrawOpt+" SAME").c_str());

 
  TLegend* leg = makeLegend<0>(legEntryVec);
  leg->Draw();
  
  c1->cd();
  TPad* ratioPad = new TPad("ratioPad", "newpad",0.01,0.01,0.99,0.33);
  ratioPad->Draw();
  ratioPad->cd();
  ratioPad->SetTopMargin(0.05);
  ratioPad->SetBottomMargin(0.3);
  //    ratioPad->SetRightMargin(0.1);
  ratioPad->SetFillStyle(0);
  
  firstHist=true;
  for(auto& denom : denoms){
    TH1* ratioHist = makeRatio(numer.hist,denom.hist,divOpt);
    auto copyHistStyle =[](TH1* refHist,TH1* newHist){
      AnaFuncs::setHistAttributes(newHist,refHist->GetLineColor(),refHist->GetLineWidth(),refHist->GetMarkerStyle(),refHist->GetLineColor());
    };
    copyHistStyle(denom.hist,ratioHist);
    ratioHist->SetTitle("");
    //  ratioHist->GetXaxis()->SetLabelSize(ratioHist->GetXaxis()->GetLabelSize()*(0.99-0.33)/0.33);
    ratioHist->GetXaxis()->SetLabelSize(0.1);
    ratioHist->GetXaxis()->SetTitleSize(0.1);
    ratioHist->GetXaxis()->SetTitle(xAxisLabel.c_str());
    ratioHist->GetYaxis()->SetLabelSize(0.1);
    ratioHist->GetYaxis()->SetTitleSize(0.1);
    ratioHist->GetYaxis()->SetTitleOffset(0.65); 
    ratioHist->GetYaxis()->SetTitle("ratio");   
   
    if(firstHist){
      ratioHist->Draw();
      firstHist=false;
    }else ratioHist->Draw("SAME");
    
  }
  spectrumPad->cd();
  if(!denoms.empty()) return denoms[0].hist;
  else return nullptr;
}



TH1* HistFuncs::adjustForRatioCanvas(TH1* ratioHist,const std::string& title)
{
  if(!title.empty()) ratioHist->SetTitle(title.c_str());
  ratioHist->GetXaxis()->SetLabelSize(0.1);
  ratioHist->GetXaxis()->SetTitleSize(0.1);
  ratioHist->GetYaxis()->SetLabelSize(0.1);
  ratioHist->GetYaxis()->SetTitleSize(0.11);
  ratioHist->GetYaxis()->SetTitleOffset(0.60); 
  ratioHist->GetYaxis()->SetTickLength(0.04);
  ratioHist->GetXaxis()->SetTickLength(0.06);
  
  return ratioHist;
}

TCanvas* HistFuncs::makeRatioCanvas(const std::string& name)
{
  TCanvas* c1 =new TCanvas("c1", "c1",900,750);
  c1->cd();
  TPad* spectrumPad = new TPad("spectrumPad", "newpad",0.01,0.30,0.99,0.99);
  spectrumPad->Draw(); 
  spectrumPad->cd();
  c1->cd();
  TPad* ratioPad = new TPad("ratioPad", "newpad",0.01,0.01,0.99,0.33);
  ratioPad->Draw();
  ratioPad->cd();
  ratioPad->SetTopMargin(0.05);
  ratioPad->SetBottomMargin(0.3);
  //    ratioPad->SetRightMargin(0.1);
  ratioPad->SetFillStyle(0);
  spectrumPad->cd();
  return c1;
}

void HistFuncs::dumpEvtList(TTree* tree,const std::string& cuts,const std::string& outputFile)
{
  tree->SetEstimate(tree->GetEntries()+1);

  std::string cuts1 = boost::algorithm::replace_all_copy(cuts,"{1}","1");
  cuts1 = boost::algorithm::replace_all_copy(cuts1,"{2}","2");
  std::string cuts2 = boost::algorithm::replace_all_copy(cuts,"{1}","2");
  cuts2 = boost::algorithm::replace_all_copy(cuts2,"{2}","1");

  const std::string finalCuts = "("+cuts1+") || ("+cuts2+")";

  tree->Draw("runnr:lumiSec:eventnr",finalCuts.c_str(),"goff");
  double* runnrArray = tree->GetV1();
  double* lumiSecArray = tree->GetV2();
  double* eventNrArray = tree->GetV3();
  size_t nrEvents = std::max(tree->GetSelectedRows() % tree->GetEstimate(),Long64_t(0));

  std::set<std::string> eventSet;
  for(size_t eventNr=0;eventNr<nrEvents;eventNr++){
    std::ostringstream eventStr;
    eventStr<<runnrArray[eventNr]<<" "<<lumiSecArray[eventNr]<<" "<<static_cast<unsigned int>(eventNrArray[eventNr])<<std::endl;
    eventSet.insert(eventStr.str());
  }

  std::ofstream file(outputFile);
  file <<"# run lum event"<<std::endl;
  file <<"# cuts: "<<finalCuts<<std::endl;
  for(const auto& eventStr : eventSet){
    file << eventStr;
    //    file <<runnrArray[eventNr]<<" "<<lumiSecArray[eventNr]<<" "<<static_cast<unsigned int>(eventNrArray[eventNr])<<std::endl;
  }
  file.close();
}


void HistFuncs::readTree(TTree* tree,const std::string& vars,const std::string& cuts,std::vector<std::array<double,4>>& output)
{
  {
    std::vector<std::array<double,4>> temp;
    output.swap(temp);
  }

  tree->SetEstimate(tree->GetEntries()+2);
  tree->Draw((vars).c_str(),cuts.c_str(),"goff");
  size_t nrEntries = tree->GetSelectedRows() % tree->GetEstimate();
  output.reserve(nrEntries);
  for(size_t entryNr=0;entryNr<nrEntries;entryNr++){
    output.emplace_back(std::array<double,4>{tree->GetV1()[entryNr],tree->GetV2()[entryNr],
			                    tree->GetV3()[entryNr],tree->GetW()[entryNr]});
  }
}

std::vector<std::vector<float> > HistFuncs::readTree(TTree* tree,const std::string& vars,const std::string& cuts)
{
  std::vector<std::vector<float> > output;

  //format of vars is "var1:var2:var3:...:varN" so number of variables  is number of : +1
  const size_t nrVars = std::count(vars.begin(),vars.end(),':')+1;
  tree->SetEstimate(tree->GetEntries()+2);
  tree->Draw(vars.c_str(),cuts.c_str(),"goff");
  const size_t nrEntries = tree->GetSelectedRows() % tree->GetEstimate();
  output.reserve(nrEntries);
  for(size_t entryNr=0;entryNr<nrEntries;entryNr++){
    std::vector<float> varsVec(nrVars,0.);
    for(size_t varNr=0;varNr<nrVars;varNr++){
      varsVec[varNr] = tree->GetVal(varNr)[entryNr];
    }
    output.emplace_back(varsVec);			
  }
  return output;
}

std::vector<std::vector<float> > HistFuncs::readTreeBothEles(TTree* tree,const std::string& vars,const std::string& cuts)
{
  std::vector<std::vector<float> > output;

  std::string vars1 = boost::algorithm::replace_all_copy(vars,"{1}","1");
  vars1 = boost::algorithm::replace_all_copy(vars1,"{2}","2");
  std::string cuts1 = boost::algorithm::replace_all_copy(cuts,"{1}","1");
  cuts1 = boost::algorithm::replace_all_copy(cuts1,"{2}","2");


  std::string vars2 = boost::algorithm::replace_all_copy(vars,"{1}","2");
  vars2 = boost::algorithm::replace_all_copy(vars2,"{2}","1");
  std::string cuts2 = boost::algorithm::replace_all_copy(cuts,"{1}","2");
  cuts2 = boost::algorithm::replace_all_copy(cuts2,"{2}","1");
  std::cout <<"1: "<<vars1<<" : "<<cuts1<<std::endl;
  std::cout <<"2: "<<vars2<<" : "<<cuts2<<std::endl;

  auto output1 = readTree(tree,vars1,cuts1);
  auto output2 = readTree(tree,vars2,cuts2);
  output1.insert(output1.end(),output2.begin(),output2.end());
  return output1;
}



std::string HistFuncs::vsVarAxisLabel(const std::string& vsVar)
{

  std::string label = getNiceName(vsVar);
  std::string units = getUnits(vsVar);
  if(!units.empty()) label+=" ["+units+"]";
  return label;

}


std::string HistFuncs::getNiceName(const std::string& var)
{

  
  //first we could have added variables together so split them into individual vars
  // std::string subVars;
  // boost::split(subVars,boost::is_any_of("*+-/"));
  // for(subVar : splitVars){
    
  //screw it, I'm just going to hardcode the patterns I want to remove
  std::string varClean =AnaFuncs::replaceSubStr(var,"ele.","");
  varClean =AnaFuncs::replaceSubStr(varClean,"ele1.","");
  varClean =AnaFuncs::replaceSubStr(varClean,"ele2.","");
  varClean =AnaFuncs::replaceSubStr(varClean,"probeHLT.","");
 

  if(varClean=="eleTruthDetEta") return "#eta^{gen}_{z=0}";
  else if(varClean=="abs(eleTruthDetEta)") return "|#eta^{gen}_{z=0}|";
  else if(varClean=="eleTruthEta") return "#eta^{gen}";
  else if(varClean=="eleTruthEt") return "E_{T}^{gen}";
  else if(varClean=="eleTruthPhi") return "#phi^{gen}";
  else if(varClean=="et/eleTruthEt") return "E_{T}^{reco}/E_{T}^{gen}";
  else if(varClean=="et") return "E_{T}";
  else if(varClean=="eta") return "#eta";
  else if(varClean=="detEta") return "#eta_{SC}"; 
  else if(varClean=="abs(detEta)") return "|#eta_{SC}|";
  else if(varClean=="phi") return "#phi";
  else if(varClean=="detPhi") return "#phi_{SC}";
  else if(varClean=="dEtaIn") return "#Delta#eta_{in}";
  else if(varClean=="dEtaInSeed") return "#Delta#eta_{in}^{seed}";
  else if(varClean=="dPhiIn") return "#Delta#phi_{in}";
  else if(varClean=="epIn") return "E/p_{in}";
  else if(varClean=="sigmaIEtaIEta") return "#sigma_{i#etai#eta}";
  else if(varClean=="sigmaIEtaIEta/0.01745") return "#sigma_{i#etai#eta}^{RAW}"; 
  else if(varClean=="sigmaIEtaIEta/0.0447") return "#sigma_{i#etai#eta}^{RAW}"; 
  else if(varClean=="e2x5Over5x5") return "E^{2x5}/E^{5x5}";
  else if(varClean=="e1x5Over5x5") return "E^{1x5}/E^{5x5}";
  else if(varClean=="hadem") return "H/E";
  else if(varClean=="hadem*clusNrgy") return "Hcal Energy";
  else if(varClean=="hademDepth1BC*clusNrgy") return "Hcal Single Tower \"Depth1\" Energy";
  else if(varClean=="hademDepth2BC*clusNrgy") return "Hcal Single Tower \"Depth2\" Energy";
  else if(varClean=="towerHadD1Nrgy") return "Hcal Single Tower \"Depth1\" Energy (Hit E>0.8 GeV)";
  else if(varClean=="towerHadD2Nrgy") return "Hcal Single Tower \"Depth2\" Energy (Hit E>0.8 GeV)";
  else if(varClean=="isolHadDepth1") return "HCAL \"Depth1\" Isol";
  else if(varClean=="isolHadDepth2") return "HCAL \"Depth2\" Isol";
  else if(varClean=="isolHadDepth1+isolHadDepth2") return "HCAL Isol";
  else if(varClean=="isolEcalClus") return "PF ECAL Cluster Isol";
  else if(varClean=="isolHcalClus") return "PF HCAL Cluster Isol";
  else if(varClean=="isolEm") return "ECAL Isol";
  else if(varClean=="isolPtTrks") return "Trk Isol";
  else if(varClean=="isolCharged") return "Charged Isol";
  else if(varClean=="isolPhoton") return "Photon Isol";
  else if(varClean=="isolNeutral") return "Neutral Isol";
  else if(varClean=="nrMissHits") return "# Miss Hits";
  else if(varClean=="dxy") return "dxy";
  else if(varClean=="nrClus") return "# sub clusters";
  else if(varClean=="nrTruePUInt") return "# true PU interactions";
  else if(varClean=="nrVert") return "# vertices";
  else if(varClean=="rho") return "#rho";
  else if(varClean=="r9") return "R9";
  else if(varClean=="bremFrac") return "(p_{in} - p_{out}) / p_{in}";

  else return var;

}
std::string HistFuncs::getUnits(const std::string& var)
{

  //screw it, I'm just going to hardcode the patterns I want to remove
  std::string varClean =AnaFuncs::replaceSubStr(var,"ele.","");
  varClean =AnaFuncs::replaceSubStr(varClean,"ele1.","");
  varClean =AnaFuncs::replaceSubStr(varClean,"ele2.","");
  varClean =AnaFuncs::replaceSubStr(varClean,"probeHLT.","");
 
  if(varClean=="eleTruthEt") return "GeV";
  else if(varClean=="et") return "GeV";
  else if(varClean=="eleTruthPhi") return "rad";
  else if(varClean=="dPhiIn") return "rad";
  else if(varClean=="phi") return "rad";
  else if(varClean=="isolHadDepth1") return "GeV";
  else if(varClean=="isolHadDepth2") return "GeV";
  else if(varClean=="isolHadDepth1+isolHadDepth2") return "GeV"; 
  else if(varClean=="hadem*clusNrgy") return "GeV";  
  else if(varClean=="hademDepth1BC*clusNrgy") return "GeV";
  else if(varClean=="hademDepth2BC*clusNrgy") return "GeV";
  else if(varClean=="towerHadD1Nrgy") return "GeV";
  else if(varClean=="towerHadD2Nrgy") return "GeV";
  else if(varClean=="isolEcalClus") return "GeV";
  else if(varClean=="isolHcalClus") return "GeV";
  else if(varClean=="isolEm") return "GeV";
  else if(varClean=="isolPtTrks") return "GeV";
  else if(varClean=="isolCharged") return "GeV";
  else if(varClean=="isolPhoton") return "GeV";
  else if(varClean=="isolNeutral") return "GeV";
  else if(varClean=="dxy") return "cm";
  else if(varClean=="rho") return "GeV";
 
  else return "";

}

TPaveLabel* HistFuncs::makeLabel(const std::string& labelTxt,float x1,float y1,float x2,float y2)
{
  TPaveLabel* label = new TPaveLabel(x1,y1,x2,y2,
				     labelTxt.c_str(),"brNDC");
  label->SetBorderSize(0);
  label->SetTextFont(42);
  label->SetFillStyle(0);
  label->SetTextAlign(12);
  label->SetTextSize(0.657895);
  return label;


}



TLegend* HistFuncs::makeLegend(const std::vector<std::pair<TGraph*,std::string>>& legEntries,float x1,float y1,float x2,float y2)
{
  TLegend* leg = new TLegend(x1,y1,x2,y2);
  for(auto& entry : legEntries){
    leg->AddEntry(entry.first,entry.second.c_str(),"LP");
  }
  leg->SetBorderSize(0);
  return leg;
}


bool HistFuncs::validExpression(const std::string& expression,TTree* tree,bool verbose)
{
  bool valid =true;
  std::vector<std::string> expressSplit;
  boost::split(expressSplit,expression,boost::is_any_of("?*/-+=&|!(),;><~% "));
  for(const auto& var : expressSplit){
    if(var.empty() || AnaFuncs::isNumber(var)) continue;
    std::vector<std::string> varSplit;
    boost::split(varSplit,var,boost::is_any_of("."));
    if(!varSplit.empty()){
      TBranch*  branch = tree->GetBranch(varSplit[0].c_str());
      for(size_t branchNr=1;branchNr+1<varSplit.size();branchNr++){
	if(branch) branch = branch->FindBranch(varSplit[branchNr].c_str());
      }
      if(varSplit.size()>1){
	TLeaf* leaf = branch ? branch->GetLeaf(varSplit.back().c_str()) : nullptr;
	if(!leaf){
	  if(verbose){
	    std::cout <<"var: \""<<var<<"\" in \""<<expression<<"\" NOT found"<<std::endl;
	    valid = false;
	  }else return false; //as not printing out each fail might as well return now
	}
      }else if(varSplit.size()==1){
	if(!branch){
	  if(verbose){
	    std::cout <<"var: \""<<var<<"\" in \""<<expression<<"\" NOT found"<<std::endl;
	    valid = false;
	  }else return false; //as not printing out each fail might as well return now
	}
      }
    }
  }
  return valid;
}
