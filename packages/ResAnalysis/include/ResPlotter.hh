#include "ResAnalysis/ResFitter.hh"
#include "Utility/AnaFuncs.hh"

#include <vector>
#include <string>
#include <iostream>

class TH2;
class TTree;
class TGraph;
class TPaveLabel;
  

/************************************************************************************************************
// This is class to make all the resolution plots necessary for validation of a new energy correction
// It is a port of code that grew organically under great time pressure so has some questionable 
// design decisions and has not really been properly brought up to standard. It is also feature incomplete
//
// flow:
//   1) the class first creates histograms from a tree and can bin the resolution data in up to 2 dimensions
//      this step is reasonably slow and could do with speed ups
//      the resolution variables are defined by cfg_.vars and it is set up for sensible defaults
//      
//      The histograms are made by: 
//      makeHists(std::vector<TTree*> trees,const std::string& label,const std::string& cuts,
//                const std::string& vsVar1,const std::string& vsVar2,
//                const std::vector<double>& vsVar1Bins,const std::vector<double>& vsVar2Bins);
//
//      * trees : a vector of the input trees which matches the size of cfg_.vars which specifies
//        the resolution variables for each tree
//      * label : the label to put on the plots, eg RealIC, EB 
//      * cuts  : the cuts to apply to trees when making the histograms 
//      * vsVar1, vsVar2:  the variables to bin the resolution in
//      * vsVar1Bins, vsVar2Bins:  the binning of the variables
//
//   2) next the fits can be done and plotted by
//      printFits(const std::vector<int>& histNrs,const std::string& baseOutName="")const
//     
//      * histNrs : the histogram / resolution variables you wish to plot (can be either 2 or 3 entries)
//        use printLabels() to see the options for this
//      * baseOutName: the prefix of the outputed plots, you are responsible for ensuring any directories
//        already exist
//
//      this step can be done multiple times for different variables
//
************************************************************************************************************/


class ResPlotter {
public:
  
  struct Config {
    
    //the binning of the resolutin variables
    int nrResBins;
    float resMin;
    float resMax;

    //fit ranges
    float fitMin;
    float fitMax;
    float fitHighThres;
    float fitMinHigh;
    float fitMaxHigh;

    bool normalise;

    //plotting options
    int binLabelPrecision;
    bool divideMeanBySigma;    

    //the resolution variables to plot for each tree
    //first vector is indexed to trees, the payload is the variables to plot for that tree
    //format is {"name in tree","name for plotting"}
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
  //so this is one of the things I was aluding to when mentioned it organically grew suboptimally
  //the first vector is binning in vsVar1, 
  //the second vector is the individual resolution variables with the 2D hist being the 
  //the resolution variable vs vsVar2 
  //as you might have noticed vsVar1 was hacked onto this setup
  //to be fixed
  std::vector<std::vector<std::pair<TH2*,std::string> > > histsVec_;
  VarNameData vsVar1_;
  VarNameData vsVar2_;
  std::vector<double> vsVar1Bins_;
  std::vector<double> vsVar2Bins_;
  std::string label_;

public:

  Config& cfg(){return cfg_;}
  const Config& cfg()const{return cfg_;}

  //size of trees must equal size of cfg_.vars (this is enforced), however you can put a nullptr into skip that
  //entry eg {tree,nullptr} will skip cfg_vars[1] entry
  void makeHists(std::vector<TTree*> trees,const std::string& label,const std::string& cuts,const std::string& vsVar1,const std::string& vsVar2,const std::vector<double>& vsVar1Bins,const std::vector<double>& vsVar2Bins);
  void printFits(const std::vector<int>& histNrs,const std::string& baseOutName="")const;
  void printLabels()const{
    if(!histsVec_.empty()){
      for(size_t histNr=0;histNr<histsVec_[0].size();histNr++){
	std::cout <<histNr<<" "<<histsVec_[0][histNr].second<<std::endl;
      }
    }
  }
  void setFitType(ResFitter::FitType fitType){
    if(fitType!=ResFitter::FitType::CB) resFitter_.setFitType(fitType);
    else std::cout <<"current CB fitting is bugged and therefore disabled"<<std::endl; 
  }

private:
  
  //makes the histograms for a single tree
  std::vector<std::vector<std::pair<TH2*,std::string> > >  
  makeHists(TTree* tree,
	    const std::vector<std::pair<std::string,std::string> >& vars,
	    const std::string& cuts)const;	    

  //prints all the plots comparing different resolutions for a given bin
  void printResComps(const std::vector<ResFitter::ParamsVsVar>& fitParams,
		     const std::string& baseName,const std::pair<float,float>& plotRange,
		     const std::string& regionStr)const;
  
  //makes a plot comparing the given fit parameter vs a variable
  TGraph* plotFitParamsVsVarComp(const std::vector<ResFitter::ParamsVsVar>& fits,
				 ResFitter::ValType valType,
				 bool divideSigmaByMean=false)const;
  //makes a plot comparing different resolutions for a given bin
  RooPlot* plotResComp(std::vector<ResFitter::Param>& fitParams,
		       const std::pair<float,float>& xRange={0,0})const;

  //formatting and utility functions
  void formatTwoComp(TGraph* graph,TPaveLabel* vsVar1Label,TPaveLabel* infoLabel,bool isMean=false)const;
  void formatThreeComp(TGraph* graph,TPaveLabel* vsVar1Label,TPaveLabel* infoLabel,bool isMean=false)const;
  static TGraph* makeRatio(TGraph* numer,TGraph* denom);
  static int getColour(unsigned int colourNr);
  static int getMarkerStyle(unsigned int markerNr);
  template<typename T> static void setStyle(T* obj,int objNr){
    int colour = getColour(objNr);
    int markerStyle = getMarkerStyle(objNr);
    AnaFuncs::setHistAttributes(obj,colour,2,markerStyle,colour);
  }
  void normaliseHists();
};


