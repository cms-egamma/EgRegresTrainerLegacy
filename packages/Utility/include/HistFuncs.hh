#ifndef HISTFUNCS
#define HISTFUNCS

//static class for root histograming functions that root is too crap to interprate (ie almost all of them)

#include "Utility/AnaFuncs.hh" //just for hist style
#include "Utility/LogErr.hh"

#include "TH1.h"
#include "TH2.h"
#include "TFile.h"
#include "TTree.h"
#include "TProfile.h"
#include "TChain.h"
#include "TGraph.h"
#include "TGraphAsymmErrors.h"
#include "TPaveLabel.h"
#include "TLegend.h"
#include "TCanvas.h"

#include <sstream>
#include <iostream>
#include <string>
#include <vector>
#include <array>


class HistFuncs { 

public:
  class Bins1D{
  private:
    int nrBins_;
    float xmin_;
    float xmax_;
    std::vector<float> binLowEdges_;
  public:
    Bins1D(int nrBins,float xmin,float xmax):nrBins_(nrBins),xmin_(xmin),xmax_(xmax){}
    Bins1D(const std::vector<float>& binLowEdges):binLowEdges_(binLowEdges){}
    Bins1D(const std::vector<double>& binLowEdges){
      for(auto binLowEdge : binLowEdges) binLowEdges_.push_back(binLowEdge);
    }
    
    TH1* makeHist(const std::string& name,const std::string& title="")const{
      if(binLowEdges_.empty()) return new TH1D(name.c_str(),title.c_str(),nrBins_,xmin_,xmax_);
      else return new TH1D(name.c_str(),title.c_str(),binLowEdges_.size()-1,&binLowEdges_[0]);
    }
  };

  struct XYCoord {
    float x1,y1,x2,y2;
    XYCoord():x1(0.485491),y1(0.437063),x2(0.873884),y2(0.536713){}
    XYCoord(float iX1,float iY1,float iX2,float iY2):x1(iX1),y1(iY1),x2(iX2),y2(iY2){}
    template<typename T> void setNDC(T* obj)const{obj->SetX1NDC(x1);obj->SetX2NDC(x2);obj->SetY1NDC(y1);obj->SetY2NDC(y2);}
  };

  struct HistOpts {
  public:
    TH1* hist;
    std::string legEntry;
    std::string legDrawOpt;
    std::string histDrawOpt;
    HistOpts(TH1* iHist,const std::string& iLegEntry,
	    const std::string& iLegDrawOpt,const std::string& iHistDrawOpt=""):
      hist(iHist),legEntry(iLegEntry),legDrawOpt(iLegDrawOpt),histDrawOpt(iHistDrawOpt){}
    
  };

  struct HistStyle{
    int colour;
    int lineWidth;
    int markerStyle;
    HistStyle(int iColour,int iLineWidth,int iMarkerStyle):
      colour(iColour),lineWidth(iLineWidth),markerStyle(iMarkerStyle){}

    void operator()(TH1* hist)const{AnaFuncs::setHistAttributes(hist,colour,lineWidth,markerStyle,colour);}
    void operator()(TGraph* hist)const{AnaFuncs::setHistAttributes(hist,colour,lineWidth,markerStyle,colour);}
  };

private:
  HistFuncs(){}
  virtual ~HistFuncs(){}

  static int nrBinsXForEffHist_;
  static float xMinForEffHist_;
  static float xMaxForEffHist_;

  static bool plotXErrOnAsymEff_;
public:
  static void setPlotXErrOnAsymEff(bool val){plotXErrOnAsymEff_=val;}

  static void plotHistProjections(TH2* inputHist,int nrProjections,int maxNrToDisplay,const char* title="",bool scale=false);

  static int getColourNr(int indx); //this has a list of colours so when I'm looping over histograms, I can supply the index to this function and get a nice colour out
 static int getMarkerNr(int indx); //this has a list of markers so when I'm looping over histograms, I can supply the index to this function and get a nice marker out
  static TChain* makeChain(const std::string& chainName,std::string filelist,int nrJobs=1,int jobNr=1,int verbose=2);
  static TChain* makeChain(const std::string& chainName,std::vector<std::string> filelist,int nrJobs=1,int jobNr=1,int verbose=2);
  static TH1* getHist(const std::string& histName,const std::string& filename);

  static TH1* makeEffHistFromTree(TTree* tree,int nrBins,float xMin,float xMax,const std::string& var,const std::string& sampleCuts,const std::string& cuts);
  static TH1* makeEffHistBothEles(TTree* tree,int nrBins,float xMin,float xMax,const std::string& var,const std::string& sampleCuts,const std::string& cuts);
  static TH1* makeEffHistBothEles(TTree* tree,int nrBins,float xMin,float xMax,const std::string& var,const std::string& sampleCuts,const std::string& cuts,const std::string& weight);
  static TH2* makeEffHistFromTree(TTree* tree,int nrBinsX,float xMin,float xMax,int nrBinsY,float yMin,float yMax,const std::string& varX,const std::string& varY,const std::string& sampleCuts,const std::string& cuts);
  static TH2* makeEffHistFromTree(TTree* tree,int nrBinsX,float xMin,float xMax,const std::vector<double>& yBins,const std::string& varX,const std::string& varY,const std::string& sampleCuts,const std::string& cuts);

  static TH1* makeIntEffHistFromTree(TTree* tree,int nrBins,float xMin,float xMax,const std::string& var,const std::string& sampleCuts,const std::string& cuts,bool intIsGreatThan);
  static TH1* makeEffHistFromTree(TTree* tree,int nrBins,float xMin,float xMax,const std::string& var,const std::string& sampleCuts,const std::string& cuts,const std::string& weight);

  static TH1* makeEffHistFromTree(TTree* tree,std::vector<float> bins,const std::string& var,const std::string& sampleCuts,const std::string& cuts,const std::string& weight="1");
  static TGraph* makeEffHistFromTreeAsymErr(TTree* tree,int nrBins,float xMin,float xMax,const std::string& var,const std::string& sampleCuts,const std::string& cuts);
  static TGraph* makeEffHistFromTreeAsymErr(TTree* tree,std::vector<float> bins,const std::string& var,const std::string& sampleCuts,const std::string& cuts);
  static TH1* makeEffHistFromTree(TTree* tree,const std::string& var,const std::string& sampleCuts,const std::string& cuts){
    return makeEffHistFromTree(tree,nrBinsXForEffHist_,xMinForEffHist_,xMaxForEffHist_,var,sampleCuts,cuts);
  }
 
  static TH1* makeEffHistFromTree1DProj(TTree* tree,int nrBinsX,float xMin,float xMax,int nrBinsY,float yMin,float yMax,const std::string& varX,const std::string& varY,const std::string& sampleCuts,const std::string& cuts);
  static TH1* makeEffHistFromTree1DProj(TTree* tree,int nrBinsX,float xMin,float xMax,const std::vector<double>& yBins,const std::string& varX,const std::string& varY,const std::string& sampleCuts,const std::string& cuts);

  static void setStdAxisTextSize(TH1* hist);

  static void setBinsForEffHist(int nrBins,float xMin,float xMax){
    nrBinsXForEffHist_=nrBins;xMinForEffHist_=xMin;xMaxForEffHist_=xMax;
  }
  
  //makes a histogram of the required cut value to obtain a desired efficiency
  static TH1* makeCutValueForFixedEffHist(TTree* tree,int nrBins,float xMin,float xMax,const std::string& vsVar,const std::string& sampleCuts,const std::string& var,float eff);
  static std::vector<TH1*> makeCutValueForFixedEffHist(TTree* tree,int nrBins,float xMin,float xMax,const std::string& vsVar,const std::string& sampleCuts,const std::string& var,const std::vector<float>& eff);
  static std::pair<float,float> calTotBkgEffForFixedSigEffHist(TTree* sigTree,TTree* bkgTree,int nrBins,float xMin,float xMax,const std::string& vsVarSig,const std::string& vsVarBkg,const std::string& sampleCutsSig,const std::string& sampleCutsBkg,const std::string& var,float eff);
  static void print(const std::string& fileName,const std::string& canvasName="c1",bool reduced=false);

  static  double getValueAtXLinear(const TH1* hist,double x);

  static TH1* makeCutValueForFixedEffHist(TH2* hist,float eff);
  static TH1* makeBkgEffForFixedSigEffHist(TTree* sigTree,TTree* bkgTree,int nrBins,float xMin,float xMax,const std::string& vsVarSig,const std::string& vsVarBkg,const std::string& sampleCutsSig,const std::string& sampleCutBkg,const std::string& var,float eff);
  static TGraph* makeBkgEffForFixedSigEffHist(TTree* sigTree,TTree* bkgTree,int nrBins,float xMin,float xMax,const std::string& vsVarSig,const std::string& vsVarBkg,const std::string& sampleCutsSig,const std::string& sampleCutsBkg,const std::string& var,const std::vector<float>& eff);
  static TGraph* makeEffVsRejCurve(const TH1* sig,const TH1* bkg,bool sigGreaterThan,bool bkgGreaterThan);
  static TGraph* makeEffVsRejCurve2D(const TH2* sig,const TH2* bkg,bool sigGreaterThan,bool bkgGreaterThan);
  static TH1* getProjection(const char* name,const char* filename,int startBin=-1,int endBin=-1);

  static TGraph* makeBestEffVsRejCurve(const TH2* sig,const TH2* bkg,bool sigGreaterThan,bool bkgGreaterThan);

  static TH1* makeDataMinusBkgPlot(const TH1* data,const TH1* bkg,float minBkgExpectInBin=3,bool ratio=true);
  static TGraph* makeDataMinusBkgPlotAsymErr(const TH1* data,const TH1* bkg,float minBkgExpectInBin);
  static TH1* makeCHist(const TH1* hist,bool intIsGreatThan=true);
  static void getRebinedBins(const TH1* hist,std::vector<float>& binLowEdges,float minInEachBin);
  static void normBinsByBinWidth(TH1* hist,float widthToNormTo);

  static void fillVecFromTree(TTree* tree,const std::string& var,const std::string& cuts,std::vector<float>& vec);
  static void readTree(TTree* tree,const std::string& vars,const std::string& cuts,std::vector<std::array<double,4>>& output);
  template <typename T>
  static std::vector<T> readTreeEntry(TTree* tree,const std::string& vars,const std::string& cuts,long entryNr);
  template<std::size_t N>
  static std::vector<std::array<double,N+1> > readTree(TTree* tree,const std::string& vars,const std::string& cuts);
  static std::vector<std::vector<float> > readTree(TTree* tree,const std::string& vars,const std::string& cuts);
  static std::vector<std::vector<float> > readTreeBothEles(TTree* tree,const std::string& vars,const std::string& cuts);
  static TH1* compTwoVars(TTree* tree,int nrBins,float min,float xmax,const std::string& var1,const std::string& var2,const std::string& cuts,const std::string& var1LegName,const std::string& var2LegName);
  static TH1* compTwoVars(TTree* tree1,TTree* tree2,int nrBins,float min,float xmax,const std::string& var1,const std::string& var2,const std::string& cuts1,const std::string& cuts2,const std::string& var1LegName,const std::string& var2LegName,bool norm);
  static TH1* compTwoCuts(TTree* tree,int nrBins,float min,float xmax,const std::string& var,const std::string& baseCuts,const std::string& cuts1,const std::string& cuts2,const std::string& var1LegName,const std::string& var2LegName);
  static TH1* compTwoEffs(TTree* tree,int nrBins,float min,float xmax,const std::string& var,const std::string& baseCuts,const std::string& cuts1,const std::string& cuts2,const std::string& var1LegName,const std::string& var2LegName,const std::string& weight);
  static TH1* makeHist(TTree*,int nrBins,float min,float xmax,const std::string& var,const std::string& cuts);
  static TH1* makeHist(TTree*,const Bins1D& bins,const std::string& var,const std::string& cuts);
  static TH2* makeHist(TTree* tree,int nrBinsX,float xmin,float xmax,int nrBinsY,float ymin,float ymax,const std::string& var,const std::string& cuts);
  static TH2* makeHist(TTree* tree,const std::vector<double>& binLowEdges,int nrBinsY,float ymin,float ymax,const std::string& var,const std::string& cuts);
  static TH1* makeHistBothEles(TTree*,int nrBins,float min,float xmax,const std::string& var,const std::string& cuts);
  static TH2* makeHistBothEles(TTree* tree,int nrBinsX,float xmin,float xmax,int nrBinsY,float ymin,float ymax,const std::string& var,const std::string& cuts);
  static TProfile* makeHistBothElesProf(TTree* tree,int nrBinsX,float xmin,float xmax,const std::string& var,const std::string& cuts);
  static TH1* plotTwoHists(TH1* hist1,TH1* hist2,const std::string& leg1Str,const std::string& leg2Str);
  static TGraph* plotTwoHists(TGraph* hist1,TGraph* hist2,const std::string& leg1Str,const std::string& leg2Str);

  static TH1* makeRatio(TH1* numerHist,TH1* denomHist,const std::string& option);
  //TH1* hist,std::string legEntry,std::string histDrawOpt;
  static TH1* plotWithRatio(HistOpts numer,HistOpts denom,
			    const std::string& divOpt);
  static TH1* plotWithRatio(HistOpts numer1,HistOpts numer2,HistOpts denom,
			    const std::string& divOpt);
  static TH1* plotWithRatio(std::vector<HistOpts> numers,HistOpts denom,
			    const std::string& divOpt);
  static TH1* plotWithRatio(HistOpts denom,std::vector<HistOpts> denoms,
			    const std::string& divOpt);
  static TH1* adjustForRatioCanvas(TH1* ratioHist,const std::string& title="");
  static TCanvas* makeRatioCanvas(const std::string& name="c1");

  static void dumpEvtList(TTree* tree,const std::string& cuts,const std::string& outputFile);

  static std::string vsVarAxisLabel(const std::string& vsVar);
  static std::string getNiceName(const std::string& var);
  static std::string getUnits(const std::string& var);
  

  static TPaveLabel* makeLabel(const std::string& labelTxt,float x1,float y1,float x2,float y2);
  static TPaveLabel* makeLabel(const std::string& labelTxt){return makeLabel(labelTxt,0.63029,0.787456,0.884187,0.853659);}
  
  template<int mode=0>
  static TLegend* makeLegend(const std::vector<std::pair<TH1*,std::string>>& legEntries,const XYCoord& xy=XYCoord());
  template<typename T>
  static TLegend* makeLegend(const std::vector<T*>& hists,const std::vector<std::string>& legEntries,const XYCoord& xy=XYCoord());
  static TLegend* makeLegend(const std::vector<std::pair<TGraph*,std::string>>& legEntries,float x1,float y1,float x2,float y2);
  static TLegend* makeLegend(const std::vector<std::pair<TGraph*,std::string>>& legEntries){
    return makeLegend(legEntries,0.485491,0.437063,0.873884,0.536713);
  }
  template<typename T,typename Container=TCanvas>
  static std::vector<T*> getFromCanvas(Container* c1,const std::string& className,int verbose=0);
  template<typename T>
  static void setCoord(T* obj,float x1,float y1,float x2,float y2);
  template<typename T>
  static void printCoord(const T* obj);
  
  template<typename T>
  static T* get(const std::string& objName,const std::string& filename);
 

  static bool validExpression(const std::string& expression,TTree* tree,bool verbose=true);

  ClassDef(HistFuncs,1)
}; 
template<typename T>
T* HistFuncs::get(const std::string& objName,const std::string& filename)
{
  TFile* file = TFile::Open(filename.c_str(),"READ");
  auto obj = static_cast<T*>( file->Get(objName.c_str()));
  if(obj) obj->SetDirectory(0); //need to sort this
  delete file;
  return obj;
}
template<typename T>
void HistFuncs::setCoord(T* obj,float x1,float y1,float x2,float y2)
{
  obj->SetX1NDC(x1);obj->SetX2NDC(x2);obj->SetY1NDC(y1);obj->SetY2NDC(y2);
}
  
template<typename T>
void HistFuncs::printCoord(const T* obj)
{
  std::cout <<","<<obj->GetX1NDC()<<","<<obj->GetY1NDC()<<","<<obj->GetX2NDC()<<","<<obj->GetY2NDC()<<std::endl;
}
  
template<int mode>
TLegend* HistFuncs::makeLegend(const std::vector<std::pair<TH1*,std::string>>& legEntries,const XYCoord& coord)
{
  TLegend* leg = new TLegend(coord.x1,coord.y1,coord.x2,coord.y2);
  for(auto& entry : legEntries){
    std::string drawOpt = "LP";
    if(mode==1 && &entry!=&legEntries.front()) drawOpt="F";
    if(mode==2 && &entry==&legEntries.back()) drawOpt="F";
    leg->AddEntry(entry.first,entry.second.c_str(),drawOpt.c_str());
  }
  leg->SetBorderSize(0);
  return leg;
}

template<typename T>
TLegend* HistFuncs::makeLegend(const std::vector<T*>& hists,const std::vector<std::string>& legEntries,const XYCoord& coord)
{
  TLegend* leg = new TLegend(coord.x1,coord.y1,coord.x2,coord.y2);
  for(size_t histNr=0;histNr<hists.size();histNr++){
    leg->AddEntry(hists[histNr],legEntries[histNr].c_str(),"LP");
  }
  leg->SetBorderSize(0);
  return leg;
}

template<typename T,typename Container>
std::vector<T*> HistFuncs::getFromCanvas(Container* c1,const std::string& className,int verbose)
{
  std::vector<T*> result;
  TIter next(c1->GetListOfPrimitives());
  while (TObject* tObj=next()) {
    if(verbose>0) std::cout <<tObj->ClassName()<<std::endl;
    if(tObj->ClassName()==className){
      T* obj  = static_cast<T*>(tObj);
      result.push_back(obj);
    }
  }
  if(verbose>0) std::cout <<"nr objs "<<result.size()<<std::endl;
  return result;
}
template<typename T>
std::vector<T> HistFuncs::readTreeEntry(TTree* tree,const std::string& vars,const std::string& cuts,long entryNr)
{
  //format of vars is "var1:var2:var3:...:varN" so number of variables  is number of : +1
  const size_t nrVars = std::count(vars.begin(),vars.end(),':')+1;
  std::vector<T> output(nrVars,0.);

  auto gErrorIgnoreLevelOld = gErrorIgnoreLevel;
  gErrorIgnoreLevel = kError;
  tree->Draw(vars.c_str(),(cuts).c_str(),"goff",1,entryNr);
  gErrorIgnoreLevel = gErrorIgnoreLevelOld;
  const size_t nrEntries = tree->GetSelectedRows() % (tree->GetEstimate() + 1);
  //a formula compile error gives -1 for selected rows so guard against
  if(tree->GetSelectedRows()>=0 && nrEntries>=1){
    if(nrEntries>1){
      LogErr<<" Error multiple entries found, should be exactly 0 or 1 entries, returning first"<<nrEntries<<std::endl;
    }
    for(size_t varNr=0;varNr<nrVars;varNr++){
      output[varNr] = tree->GetVal(varNr)[0];
    } 
  }
  return output;
}
template<std::size_t N>
std::vector<std::array<double,N+1> > 
HistFuncs::readTree(TTree* tree,const std::string& vars,const std::string& cuts)
{
  std::cout <<"reading tree "<<vars<<" with cuts "<<cuts<<std::endl;

  typedef  Double_t* (TTree::*TTreeFunc)();
  std::vector<TTreeFunc> funcs={&TTree::GetV1,&TTree::GetV2,&TTree::GetV3};

  std::vector<std::array<double,N+1> > output;
  tree->SetEstimate(tree->GetEntries()+2);
  tree->Draw((vars).c_str(),cuts.c_str(),"goff");
  size_t nrEntries = tree->GetSelectedRows() % tree->GetEstimate();
  output.resize(nrEntries);
  for(size_t entryNr=0;entryNr<nrEntries;entryNr++){
    for(size_t varNr=0;varNr<N;varNr++){ //last entry of N is the weight var
      output[entryNr][varNr]=(tree->*funcs[varNr])()[entryNr];
    }
    output[entryNr][N]=tree->GetW()[entryNr];
  }
  return output;
}

      
#endif
