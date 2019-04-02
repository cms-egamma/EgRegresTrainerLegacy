#ifndef ANAFUNCS
#define ANAFUNCS

#include "TH1D.h"
#include "TH2.h"
#include "TF1.h"
#include "TTree.h"
#include "TLorentzVector.h"
#include "TGraph.h"

#include <ctime>
#include <cmath>
#include <string>
#include <vector>
#include <algorithm>
#include <iostream>
#include <fstream>
#include <list>
#include <sstream>


//Author: Sam Harper (RAL)

//a misc set of common analysis functions that I use that dont really fit elsewhere
//this class is not supposed to instantanced, ever.
//perhaps a namespace might be better



class AnaFuncs{
public:
  bool isNextToRingBoundary(int detId);//hack shouldnt be here;
  bool isNextToRingBoundary(int ix,int iy);

public:
  //First public classes in this namespace/class

  //lite class to handle program timing
  class Timer{
  private:
    time_t _startTime;
  public:
    Timer();
    ~Timer(){}
    void reset();
    void printTimeElapsed()const;
    friend std::ostream &operator <<(std::ostream& output,const Timer &timer);
    std::ostream &print(std::ostream& output)const;
  };

  //small class to find duplicate events
  class DupEventFinder{
  private:
    std::vector<int> _runList,_eventList; //storing the previous event and run numbers
    int _nrDupEvents; //number of duplicate events
  public:
    DupEventFinder():_nrDupEvents(0){}
    ~DupEventFinder(){}
    void checkEvent(int runnr,int eventnr);
  };

  //small class to check runnr against a goodrun list
  class GoodrunChecker{
  private:
    std::vector<int> _goodrunList; //the list of all the good runs
    //both tempory internal state parameters, hence mutable
    mutable int _lastRunnr;
    mutable bool _lastRunGood; //the parameters of the last run checked
  public:
    GoodrunChecker(const char* goodrunFilename="");
    ~GoodrunChecker(){}
    bool isGoodRun(int runnr)const;
    void setGoodrunList(const char* goodrunFilename="");
  };


 
  class GoodLumiChecker{
    
  public:
    struct GoodLumiData {
      int runnr;
      std::vector<std::pair<int,int> > allowedLumis;
      bool operator<(const GoodLumiData& rhs)const{return runnr<rhs.runnr;}
      bool operator<(const int rhsRunnr)const{return runnr<rhsRunnr;}
      bool operator==(const int rhsRunnr)const{return runnr==rhsRunnr;}
    };

    class GoodLumiDataComp {
    public:
      bool operator()(const GoodLumiData& lhs,const GoodLumiData& rhs)const{return lhs.runnr<rhs.runnr;}
      bool operator()(int lhs,const GoodLumiData& rhs)const{return lhs<rhs.runnr;}
      bool operator()(int lhs,int rhs)const{return lhs<rhs;}
      bool operator()(const GoodLumiData& lhs,int rhs)const{return lhs.runnr<rhs;}
    };
  private:
    std::vector<GoodLumiData> goodLumis_; //sorted by runnr, each run number can only enter once
    mutable int indexOfLastRun_; //temporary cache hence mutable
  public:
    GoodLumiChecker(const std::string& goodLumiList="");
    ~GoodLumiChecker(){}
    bool isGoodLumiSec(int runnr,int lumiSec)const;
    void setGoodLumiList(const std::string& goodLumiList);
    void clear(){goodLumis_.clear();indexOfLastRun_=-1;}
  private:
    void addLumiRange_(int runnr,int lowerLumi,int upperLumi);
  };
  
  class RunList{
  private:
    std::vector<int> _runList;
  public:
    RunList(){}
    ~RunList(){}
    void addRunnr(int runnr){if(!checkRunnr(runnr))_runList.push_back(runnr);}
    bool checkRunnr(int runnr)const;
    int size()const{return _runList.size();}
    int operator[](int i)const{return _runList[i];}
    void clear(){_runList.clear();}
  };
  class Blinder{
  private:
    int _runnrs[3];
    double _masses[3];
  public:
    Blinder(){};
    ~Blinder(){}
    bool blind(int regionCode,int runnr,double mass=99999.)const;
    bool blind(int runnr,double mass=99999.)const;
    void setBlindPara(int regionCode,int runnr,double mass);
    void setBlindPara(int runnr,double mass);
  };

  class LumiData{
  private:
    //run, lumi of that run
    std::vector<std::pair<int,float> > data_;
 
    class RunComp {
    public:
      bool operator()(const std::pair<int,float>& lhs,const std::pair<int,float>& rhs)const{return lhs.first<rhs.first;}
      bool operator()(int lhs,const std::pair<int,float> rhs)const{return lhs<rhs.first;}
      bool operator()(int lhs,int rhs)const{return lhs<rhs;}
      bool operator()(const std::pair<int,float>& lhs,int rhs)const{return lhs.first<rhs;}
    }; 

  public:
    LumiData(){}
    ~LumiData(){}
    
    void readLumiData(const std::string& filename); 

    size_t size()const{return data_.size();}
    std::pair<int,float> runLumi(size_t index){return data_[index];}
    
    float lumi(int runnr)const;
    void clear(){data_.clear();}
  };

  class PUReWeighter{
  private:
    std::vector<float> weights_;
    int nrBins_;
    float xmin_,xmax_;
  public:
    PUReWeighter(){}
    ~PUReWeighter(){}
    
    bool loadWeights();
    bool loadWeights(const std::string& dataFilename,const std::string& mcFilename);
    float weight(int nrPUInt);
  };

 
private:
  static std::string _eleDataContents;
  AnaFuncs(){} //to make sure no instantancing is posible
  virtual ~AnaFuncs(){} //to make sure no instantancing is posible
 
  static const Blinder _blinder;

public:

  static std::string floatToStr(float val);
  static void splitStrings(const char* charString,std::vector<std::string>& array,const char* seperator=":");
  static int getIndexInList(const char* var,const char* listOfVar,const char* seperator=":");
  static std::string replaceSubStr(const std::string& orgString,const std::string& subStrToReplace,const std::string& replacement);
  static std::string makeStringOfBinsLowEdges(int nrBins,float xmin,float xmax);
  static void makeVecFromInputString(std::vector<int>& vec,const std::string& inputStr);
  static void makeVecFromInputString(std::vector<float>& vec,const std::string& inputStr);
  static void makeVecFromInputString(std::vector<double>& vec,const std::string& inputStr);  
  
  static std::vector<float> makeVecFromInputStringF(const std::string& inputStr){std::vector<float> retVal;makeVecFromInputString(retVal,inputStr);return retVal;} 
  static std::vector<double> makeVecFromInputStringD(const std::string& inputStr){std::vector<double> retVal;makeVecFromInputString(retVal,inputStr);return retVal;}
  static void copyArrayToVec(const double* array,size_t nrEntries,std::vector<double>& vec);
  static void readFilelistFromPattern(const std::string& filelistPattern,std::vector<std::string> &filenames);
  static void readFilelistFromFile(const std::string& fileListName,std::vector<std::string> &filenames);
  static void readFilelist(std::string fileListName,std::vector<std::string> &filenames,int nrJobs=1,int jobNr=1,int verbose=2); 
  static void readFilelist(std::vector<std::string> fileListName,std::vector<std::string> &filenames,int nrJobs=1,int jobNr=1,int verbose=2);
  static void splitFilelist(int nrJobs,int jobNr,std::vector<std::string>& filenames);
  static void splitFilelistMixed(int nrJobs,int jobNr,std::vector<std::string>& filenames);//files are maxmally mixed between jobs
  static void splitFilelistConsecutive(int nrJobs,int jobNr,std::vector<std::string>& filenames);//each job has files that were consecutive to each other
  static int nrFilesInJob(int nrFiles,int nrJobs,int jobNr);
  static int nrFilesInPreviousJobs(int nrFiles,int nrJobs,int jobNr);
  static void addJobNrToFilename(char* filename,int jobNr);
  static void addJobNrToFilename(std::string& filename,int jobNr);
 
  static int getBinNr(const TH1* theHist,double x); //returns the bin in which x is
  static int getBinNr(const std::vector<float>& bins,double x);
  static int getBinNr(int nrBins,double min,double max,double x);
  static int getXBinNr(const TH2* theHist,double x); //returns the xbin in which x is
  static int getYBinNr(const TH2* theHist,double y); //returns the ybin in which y is
  static float getDistInBin(int nrBins,float min,float max,float val); //returns the percentage of the bin the xval is along in the bin

  static void convert1DHistToVec(const TH1* hist,std::vector<float>& vec);

  static float getBinLowEdge(const std::vector<float>& binLowEdges,size_t binNr);
  static float getBinHighEdge(const std::vector<float>& binLowEdges,size_t binNr);
  static float getBinLowEdge(const int nrBins,const float xmin,const float xmax,const int binNr){return nrBins!=0 ? (xmax-xmin)/nrBins*(binNr-1)+xmin : -999;}

  static bool hasNonZeroAND(int setOfBits1,int setOfBits2){return (setOfBits1&setOfBits2)!=0x0;}

  static double getBinContent(const TH1* hist,double x){return hist->GetBinContent(getBinNr(hist,x));} //returns the content of the bin in which x is 
  static double getBinContent(const TH2* hist,double x,double y){
    return hist->GetBinContent(getXBinNr(hist,x),getYBinNr(hist,y));
  }  
  static double getNrInRange(const TH1* theHist,double min,double max);
  static double getErrInRange(const TH1* theHist,double min,double max);
  static void setHistAttributes(TH1* theHist,int lineColour=-1,int lineWidth=-1,int markerStyle=-1,int markerColour=-1);
  static void setHistAttributes(TGraph* theHist,int lineColour=-1,int lineWidth=-1,int markerStyle=-1,int markerColour=-1);
 
  static void getHistIntegral(const TH1* theHist,double xMin,double xMax,double& nrEvents,double& nrEventsErr);  
  static double getHistIntegral(const TH1* theHist,double xMin,double xMax);
  static std::pair<double,double> getHistIntAndErr(const TH1* theHist,double xMin,double xMax){
    std::pair<double,double> integral;
    getHistIntegral(theHist,xMin,xMax,integral.first,integral.second);
    return integral;
  }
  static void setHistAxes(TH1* hist);

  static TH1* getHistFromFile(const std::string& histName,const std::string& filename);

  static bool isBlind(int regionCode,int runnr,double mass=999999.){return _blinder.blind(regionCode,runnr,mass);}
  static bool isBlind(int runnr,double mass=999999.){return _blinder.blind(runnr,mass);}

  //sets a branch address and enables it at the same time
  static void setBranchAddress(TTree *tree,const char* branchName,void *address);

  static void getDatasetID(const char *baseFilename,std::string &datasetID);

  static void boost_replace_all(std::string& orgString,const std::string& pattern,const std::string& subTxt);

  static std::string convertToTTreeStr(int val); //for -ve vals, string is M instead 
  static std::string convertToTTreeStr(float val); //for -ve vals, string is M instead and . replaced by p
  static std::string convertToTTreeStr(double val); //for -ve vals, string is M instead and . replaced by p
  static std::string convertToStr(int val){std::ostringstream retVal;retVal<<val;return retVal.str();}
  static std::string getenv(const std::string& var){
    auto val = std::getenv(var.c_str()); return val ? std::string(val) : std::string();
  }
  //static double calChi2(const TH1* hist1,const TH1* hist2,double minRange,double maxRange,double minBinContent);
  static bool isNumber(const std::string& str);


 private:
  static void _fillLogFacLookup(int n);
 //  template <class type1,class type2> static type2 getMapValue(std::map<type1><type2> theMap,type1 theKey);

  ClassDef(AnaFuncs,1)

};


#endif
