#include "RegresTester/ResFitter.hh"
#include "Utility/AnaFuncs.hh"

#include <vector>
#include <string>
#include <iostream>

class TH2;
class TTree;
class TGraph;
class TPaveLabel;
  
class RegValidator {
public:
  
  struct Config {
    int nrResBins;
    float resMin;
    float resMax;

    float fitMin;
    float fitMax;
    float fitHighThres;
    float fitMinHigh;
    float fitMaxHigh;

    int binLabelPrecision;
    // std::string binVar;
    // std::string binVarOutname;
    bool divideMeanBySigma;    

    std::vector<std::vector<std::pair<std::string,std::string> > > vars;
    Config(){setDefaults();}
    void setDefaults();
  };

  struct VarNameData {
    std::string name;
    std::string plotname;
    std::string filename;
    std::string unit;
    VarNameData(std::string name=""):name(std::move(name)){
      autoFill();
    }
    VarNameData& operator=(const std::string& rhs){
      name = rhs;
      autoFill();
      return *this;
    }
    void autoFill();
    std::string axisLabel()const{
      if(unit.empty()) return plotname;
      else return plotname+" "+unit;
    }
  };

private:
  Config cfg_;
  
  ResFitter resFitter_;

  //internal histogram cache
  std::vector<std::vector<std::pair<TH2*,std::string> > > histsVec_;
  VarNameData vsVar1_;
  VarNameData vsVar2_;
  std::vector<double> vsVar1Bins_;
  std::vector<double> vsVar2Bins_;
  std::string label_;


public:
  Config& cfg(){return cfg_;}
  const Config& cfg()const{return cfg_;}

  void makeHists(std::vector<TTree*> trees,const std::string& label,const std::string& cuts,const std::string& vsVar1,const std::string& vsVar2,const std::vector<double>& vsVar1Bins,const std::vector<double>& vsVar2Bins);
  std::vector<std::vector<std::pair<TH2*,std::string> > >  
  makeHists(TTree* tree,
	    const std::vector<std::pair<std::string,std::string> >& vars,
	    const std::string& cuts)const;
  void printResHists(const std::vector<int>& histNrs,const std::string& baseOutName="")const;
  void printLabels()const{
    if(!histsVec_.empty()){
      for(size_t histNr=0;histNr<histsVec_[0].size();histNr++){
	std::cout <<histNr<<" "<<histsVec_[0][histNr].second<<std::endl;
      }
    }
  }

private:
  void formatTwoComp(TGraph* graph,TPaveLabel* vsVar1Label,TPaveLabel* infoLabel,bool isMean=false)const;
  void formatThreeComp(TGraph* graph,TPaveLabel* vsVar1Label,TPaveLabel* infoLabel,bool isMean=false)const;
  
  TGraph* plotFitParamsVsVarComp(const std::vector<ResFitter::ParamsVsVar>& fits,
				 ResFitter::ValType valType,
				 bool divideSigmaByMean=false)const;
  RooPlot* plotResComp(std::vector<ResFitter::Param>& fitParams,
		       const std::pair<float,float>& xRange={0,0})const;

  void printResComps(const std::vector<ResFitter::ParamsVsVar>& fitParams,
		     const std::string& baseName,const std::pair<float,float>& plotRange,
		     const std::string& regionStr)const;
  
  static TGraph* makeRatio(TGraph* numer,TGraph* denom);
  static int getColour(unsigned int colourNr);
  static int getMarkerStyle(unsigned int markerNr);
  template<typename T> static void setStyle(T* obj,int objNr){
    int colour = getColour(objNr);
    int markerStyle = getMarkerStyle(objNr);
    AnaFuncs::setHistAttributes(obj,colour,2,markerStyle,colour);
  }
    
};


