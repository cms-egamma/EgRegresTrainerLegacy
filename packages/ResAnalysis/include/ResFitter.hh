
#include <sstream>
#include <iomanip>
#include <vector>

class RooPlot;
class RooRealVar;
class TH1;
class TH2;
class TGraph;

class ResFitter { 
public:
  enum class FitType{CB,Cruijff,DCB};
  enum class ValType{Mean,Sigma,SigmaL,SigmaR};
  
  struct Param {
    float scale,scaleErr;
    float sigma,sigmaErr;
    float sigmaR,sigmaRErr;
    float sigmaL,sigmaLErr;
    RooPlot* plot; //we do not own this
    std::string legName;
    FitType fitType;
    int precision;
    Param();
    Param(FitType fitType,int precision=3):fitType(fitType),precision(precision){}
    
    void fill(const RooRealVar& iScale,const RooRealVar& iSigma,FitType iFitType,RooPlot* iPlot,const std::string& iLegName="");
    void fill(const RooRealVar& iScale,const RooRealVar& iSigmaL,const RooRealVar& iSigmaR,
	      RooPlot* iPlot,const std::string& iLegName="");

    std::string scaleLabel(){
      std::ostringstream retVal;
      retVal<<std::fixed<<std::setprecision(precision)<<"mean: "<<scale<<" #pm "<<scaleErr<<"^{stat.}";
      return retVal.str();
    }
    std::string scaleAndResCBLabel(){
      std::ostringstream retVal;
      retVal<<std::fixed<<std::setprecision(precision)<<"CB mean: "<<scale<<" CB #sigma: "<<sigma;
      return retVal.str();
    }
    std::string scaleAndResCruijffLabel(){
      std::ostringstream retVal;
      retVal<<std::fixed<<std::setprecision(precision)<<" peak: "<<scale<<" #sigma_{L}: "<<std::endl<<sigmaL<<" #sigma_{R}: "<<sigmaR;
      return retVal.str();
    }
    std::string scaleAndResLabel(){
      if(fitType==FitType::Cruijff) return scaleAndResCruijffLabel();
      else return scaleAndResCBLabel();
    }
    std::string resCBLabel(){
      std::ostringstream retVal;
      retVal<<std::fixed<<std::setprecision(precision)<<"#sigma_{CB}: "<<sigma<<" #pm "<<sigmaErr<<"^{stat.}";
      return retVal.str();
    }
    std::string resCruijffLabel(){
      std::ostringstream retVal;
      retVal<<std::fixed<<std::setprecision(precision)<<"#sigma_{L}: "<<sigmaL<<" #pm "<<sigmaLErr<<"^{stat.}"<<" #sigma_{R}: "<<sigmaR<<" #pm "<<sigmaRErr<<"^{stat.}";
      return retVal.str();
    }
    std::string resLabel(){
      if(fitType==FitType::Cruijff) return resCruijffLabel();
      else return resCBLabel();
    }
  };

  class ParamsVsVar {
  private:
    std::vector<Param> params_; //without under/overflows
    std::vector<double> binLowEdges_;//size = params size +1
  public:
    ParamsVsVar(std::vector<Param> params,std::vector<double> binLowEdges):
      params_(std::move(params)),binLowEdges_(std::move(binLowEdges)){}
    TGraph* makeGraph(ValType valType,bool divideSigmaByMean=false)const;   
    std::vector<Param>& params(){return params_;}
    const std::vector<Param>& params()const{return params_;}
    std::vector<double> binLowEdges(){return binLowEdges_;}
    const std::vector<double> binLowEdges()const{return binLowEdges_;}
    std::string legName()const{return !params_.empty() ? params_[0].legName : "";}
  };

private:
  FitType fitType_;
  bool fixAlphaDCB_;
  bool fixMeanDCB_;
 
public:
  ResFitter():fitType_(FitType::Cruijff),fixAlphaDCB_(false),fixMeanDCB_(false){}
  Param makeFit(TH1* hist,float xmin,float xmax,const std::string& fitVarName="")const;
  Param makeFit(TH2* hist2D,int binNr,float xmin,float xmax,const std::string& fitVarName="")const;
  Param makeCBFit(TH1* hist,float xmin,float xmax,const std::string& fitVarName="")const;
  Param makeDCBFit(TH1* hist,float xmin,float xmax,const std::string& fitVarName="")const;
  Param makeCruijffFit(TH1* hist,float xmin,float xmax,const std::string& fitVarName="")const;
  ParamsVsVar makeFitVsVar(TH2* hist2D,float fitMin,float fitMax,const std::string& fitVarName="")const;

  void setFitType(FitType fitType){fitType_ = fitType;}
  void setFixDCB(bool fixMeanDCB,bool fixAlphaDCB){
    fixMeanDCB_ = fixMeanDCB;
    fixAlphaDCB_ = fixAlphaDCB;
  }
  
};
