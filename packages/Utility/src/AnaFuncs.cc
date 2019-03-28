#include "Utility/AnaFuncs.hh"
#include "Utility/LogErr.hh"
#include "Utility/TempFuncs.hh"

#include <limits>
#include <cstring>
#include <sstream>
#include <ios>

#include "Utility/DetIdTools.hh" //hack which doesnt belong here

#include "TFile.h"

#include <boost/algorithm/string/replace.hpp>

ClassImp(AnaFuncs)


const AnaFuncs::Blinder AnaFuncs::_blinder;


//hack which doesnt belong here
bool AnaFuncs::isNextToRingBoundary(int detId)
{
  return DetIdTools::isNextToRingBoundary(detId);
}
bool AnaFuncs::isNextToRingBoundary(int ix,int iy)
{
  int detId = DetIdTools::makeEcalEndcapId(ix,iy,1);
  return DetIdTools::isNextToRingBoundary(detId);

}

AnaFuncs::Timer::Timer():
  _startTime(time(NULL))
{
  //no operations required
}


void AnaFuncs::Timer::reset()
{
  _startTime = time(NULL);
}

void AnaFuncs::Timer::printTimeElapsed()const
{
  std::cout << "Program run time = "<<*this<<" (hh:mm:ss)"<<std::endl;
}

std::ostream& operator <<(std::ostream& output,const AnaFuncs::Timer &timer)
{
  return timer.print(output);
}

std::ostream& AnaFuncs::Timer::print(std::ostream& output)const
{
  int timeTaken = time(NULL) - _startTime;
  int nrSecs = timeTaken%60;
  int nrMins = ((timeTaken-nrSecs)/60)%60;
  int nrHours = (timeTaken-nrSecs-nrMins*60)/3600;
  if(nrHours<10) output <<"0";
  output <<nrHours<<":";
  if(nrMins<10) output <<"0";
  output <<nrMins<<":";
  if(nrSecs<10) output <<"0";
  output <<nrSecs;

  return output;
}

void AnaFuncs::setHistAxes(TH1* hist)
{
  hist->GetXaxis()->SetTitleSize(0.047);
  hist->GetXaxis()->SetTitleOffset(1.0);
  hist->GetYaxis()->SetTitleSize(0.047);
  hist->GetYaxis()->SetTitleOffset(1.3);
}

void AnaFuncs::getHistIntegral(const TH1* theHist,double xMin,double xMax,double& nrEvents,double& nrEventsErr)
{
  nrEvents=0.;
  nrEventsErr=0.;
  int lowBinNr=getBinNr(theHist,xMin);
  int highBinNr=getBinNr(theHist,xMax);
  if(theHist->GetBinLowEdge(highBinNr)==xMax) highBinNr--; //dont include bin as it only has entries above xMax
  for(int binNr=lowBinNr;binNr<=highBinNr;binNr++){
    nrEvents+=theHist->GetBinContent(binNr);
    nrEventsErr+=theHist->GetBinError(binNr)*theHist->GetBinError(binNr);
  }
  nrEventsErr=sqrt(nrEventsErr);
}

TH1* AnaFuncs::getHistFromFile(const std::string& histName,const std::string& filename)
{
  TFile file(filename.c_str(),"READ");
  TH1* hist = static_cast<TH1*>(file.Get(histName.c_str()));
  hist->SetDirectory(0);
  return hist;
}

double AnaFuncs::getHistIntegral(const TH1* theHist,double xMin,double xMax)
{
  
  double nrEvents=0.;
  double nrEventsErr=0.;
  getHistIntegral(theHist,xMin,xMax,nrEvents,nrEventsErr);
  return nrEvents;
}


int AnaFuncs::getBinNr(const TH1* theHist,double x)
{
  //nasty logic but works perfectly
 int binNr = 0;
  while(binNr<=theHist->GetNbinsX() && x>=theHist->GetBinLowEdge(binNr+1)){
    binNr++;
  }
  return binNr;
}

//binLowEdges.size = GetNbinX()+1
int AnaFuncs::getBinNr(const std::vector<float>& binLowEdges,double x)
{
  if(binLowEdges.empty()) return 0;
  
  if(x<binLowEdges[0]) return 0; //underflow
  if(x>=binLowEdges.back()) return binLowEdges.size(); //overflow 
  
  for(size_t index=0;index<binLowEdges.size()-1;index++){
    if(x>=binLowEdges[index] && x<binLowEdges[index+1]) return index+1;
  }
  LogErr << " Error logic fail, shouldnt ever get here"<<std::endl;
  return 0;
  
}

int AnaFuncs::getBinNr(int nrBins,double min,double max,double x)
{
  //uses root logic
  const float step = (max-min)/nrBins;

  //first over, underflow
  if(x<min) return 0; 
  if(x>=max) return nrBins+1;
  
  for(int binNr=1;binNr<=nrBins;binNr++){
    float lowEdge = step*(binNr-1)+min;
    float highEdge = step*(binNr)+min;
    if(x>=lowEdge && x<highEdge) return binNr;
  }
  
  LogErr << " Error logic fail, shouldnt ever get here"<<std::endl;
  
  return -1;
}

int AnaFuncs::getXBinNr(const TH2* theHist,double x)
{
  //nasty logic but works perfectly
  const TAxis* xAxis = theHist->GetXaxis(); 
  int binNr = 0;
  while(binNr<=xAxis->GetNbins() && x>=xAxis->GetBinLowEdge(binNr+1)){
    binNr++;
  }
  return binNr;
}

int AnaFuncs::getYBinNr(const TH2* theHist,double y)
{
  //nasty logic but works perfectly
  const TAxis* yAxis = theHist->GetYaxis(); 
  int binNr = 0;
  while(binNr<=yAxis->GetNbins() && y>=yAxis->GetBinLowEdge(binNr+1)){
    binNr++;
  }
  return binNr;
}

float AnaFuncs::getDistInBin(int nrBins,float min,float max,float val)
{
  int binNr= getBinNr(nrBins,min,max,val);
  float stepSize = (max-min)/nrBins;
  float binLowEdge = (binNr-1)*stepSize+min;
  return (val-binLowEdge)/stepSize;
  
}

void AnaFuncs::convert1DHistToVec(const TH1* hist,std::vector<float>& vec)
{
  vec.clear();
  for(int binNr=0;binNr<=hist->GetNbinsX()+1;binNr++){
    vec.push_back(hist->GetBinContent(binNr));
  }
}

//logic not checked
float AnaFuncs::getBinLowEdge(const std::vector<float>& binLowEdges,size_t binNr)
{
  if(binNr==0) return std::numeric_limits<float>::min();
  else if(binLowEdges.empty() || binNr>binLowEdges.size()) return 0.;
  else return binLowEdges[binNr-1];
}
//logic not checked
float AnaFuncs::getBinHighEdge(const std::vector<float>& binLowEdges,size_t binNr)
{
  
  if(binLowEdges.empty() || binNr>binLowEdges.size()) return 0.;
  else if(binNr==binLowEdges.size()) return std::numeric_limits<float>::max();
  else return binLowEdges[binNr];
}

double AnaFuncs::getNrInRange(const TH1* theHist,double min,double max)
{
  int minBin = getBinNr(theHist,min);
  int maxBin = getBinNr(theHist,max);
  double binContent = 0.;
  for(int binNr=minBin;binNr<maxBin;binNr++){
    binContent+=theHist->GetBinContent(binNr);
  }
  return binContent;
}

double AnaFuncs::getErrInRange(const TH1* theHist,double min,double max)
{
  int minBin = getBinNr(theHist,min);
  int maxBin = getBinNr(theHist,max);
  double binErr = 0.;
  for(int binNr=minBin;binNr<maxBin;binNr++){
    binErr+=theHist->GetBinError(binNr)*theHist->GetBinError(binNr);
  }
  return sqrt(binErr);
}


            
void AnaFuncs::setHistAttributes(TH1* theHist,int lineColour,int lineWidth,int markerStyle,int markerColour)
{
  if(lineColour!=-1) theHist->SetLineColor(lineColour);
  if(lineWidth!=-1) theHist->SetLineWidth(lineWidth);
  if(markerStyle!=-1) theHist->SetMarkerStyle(markerStyle);
  if(markerColour!=-1) theHist->SetMarkerColor(markerColour);
}

void AnaFuncs::setHistAttributes(TGraph* theHist,int lineColour,int lineWidth,int markerStyle,int markerColour)
{
  if(lineColour!=-1) theHist->SetLineColor(lineColour);
  if(lineWidth!=-1) theHist->SetLineWidth(lineWidth);
  if(markerStyle!=-1) theHist->SetMarkerStyle(markerStyle);
  if(markerColour!=-1) theHist->SetMarkerColor(markerColour);
}

// inline double AnaFuncs::gaus(double x,double mean,double sigma)
// {
//   double gausNorm = 1./sqrt(2*PI()*sigma*sigma); 
//   double gausExp = -.5*(x-mean)*(x-mean)/(sigma*sigma);
  
//   return gausNorm*exp(gausExp);
// }

void AnaFuncs::DupEventFinder::checkEvent(int runnr,int eventnr)
{
  int runnrIndex = -1;
    for(int i=0;i<(int)_runList.size();i++){
      if(_runList[i]==runnr){
	runnrIndex = i;
	break;
      }
    }
    if(runnrIndex == -1){
      _runList.push_back(runnr);
      runnrIndex = _runList.size()-1;
    }
    int eventnrIndex = -1;
    for(int i=0;i<(int)_eventList.size();i++){
      if(_eventList[i]==eventnr){
	std::cout <<"error, duplicate event: runnr ="<<runnr<<" eventnr ="<<eventnr<<std::endl;
	_nrDupEvents++;
	break;
      }
    }
    if(eventnrIndex == -1) _eventList.push_back(eventnr);
}

AnaFuncs::GoodrunChecker::GoodrunChecker(const char* goodrunFilename):
  _lastRunnr(-1),
  _lastRunGood(false)
{

  setGoodrunList(goodrunFilename);
}

void AnaFuncs::GoodrunChecker::setGoodrunList(const char* goodrunFilename)
{
  _goodrunList.clear();

  if(strcmp(goodrunFilename,"null")!=0){
    std::ifstream goodrunFile(goodrunFilename);
    if(goodrunFile.bad()){
      std::cout <<"Good runlist not found"<<std::endl;
    } else if(goodrunFile.eof()){
      std::cout <<"Good runlist found but empty"<<std::endl;
    }

    
    int lastRunnr = -1;
    if(!goodrunFile.bad()){
      while(!goodrunFile.eof()){
	char tempBuffer[64];
	goodrunFile.getline(tempBuffer,64);
	char* runnrChar = strtok(tempBuffer," ");
	if(runnrChar!=NULL){
	  int runnr = atoi(runnrChar);
	  if(runnr < lastRunnr){
	    std::cout <<"AnaFuncs::GoodrunChecker::setGoodrunlist Error runlist is not in order this run "<<runnr<<" last runnr "<<lastRunnr<<std::endl;
	  }else lastRunnr = runnr;
	  //if(lastRunnr!=runnr) _goodrunList.push_back(runnr);
	  _goodrunList.push_back(runnr);
	}
      }
    }

    //  _goodrunList.sort();
  }//end of file name check
}

bool AnaFuncs::GoodrunChecker::isGoodRun(int runnr)const
{
  if(_goodrunList.size()==0) return true;
  if(runnr!=_lastRunnr){
    _lastRunGood = std::binary_search(_goodrunList.begin(),_goodrunList.end(),runnr);
    // _lastRunGood =false;
    // for(unsigned i=0;i<_goodrunList.size();i++){
    //  if(runnr==_goodrunList[i]){
    //	_lastRunGood = true;
    //	break;
    //  }
    //}
  }
  return _lastRunGood;
}

AnaFuncs::GoodLumiChecker::GoodLumiChecker(const std::string& goodLumiList)
{
  setGoodLumiList(goodLumiList);

}


void AnaFuncs::GoodLumiChecker::setGoodLumiList(const std::string& goodLumiList)
{
  clear();

  if(!goodLumiList.empty()){
    std::ifstream goodLumiFile(goodLumiList.c_str());
    if(goodLumiFile.bad()){
      std::cout <<"Good lumi list not found"<<std::endl;
    } else if(goodLumiFile.eof()){
      std::cout <<"Good lumi list found but empty"<<std::endl;
    }

    if(!goodLumiFile.bad()){
      while(!goodLumiFile.eof()){
	char tempBuffer[128];
	goodLumiFile.getline(tempBuffer,128);
	
	char* runAndLumiSec = strtok(tempBuffer,"'");
	runAndLumiSec = strtok(NULL,"'");
	if(runAndLumiSec!=NULL){
	  char* runLumi1 = strtok(runAndLumiSec,"-");
	  char* runLumi2 = strtok(NULL,"-");
	  
	  char* run1Char = strtok(runLumi1,":");
	  char* lumi1Char = strtok(NULL,":");
	  char* run2Char = strtok(runLumi2,":");
	  char* lumi2Char = strtok(NULL,":");
	  int run1 = atoi(run1Char);
	  int run2 = atoi(run2Char);
	  int lumi1 = atoi(lumi1Char);
	  int lumi2 = atoi(lumi2Char);

	  if(run1!=run2) std::cout <<"AnaFuncs::GoodLumiChecker::setGoodLumiList: Error run number missmatch "<<run1<<" : "<<run2<<std::endl;
	  
	  // std::cout <<" \'"<<run1<<":"<<lumi1<<"-"<<run2<<":"<<lumi2<<"\',"<<std::endl;
	  addLumiRange_(run1,lumi1,lumi2);
	  
   
	  

	}//end runAndLumiSec is valid
      }//end of file
      std::sort(goodLumis_.begin(),goodLumis_.end());
    }//bad file check
    
    
  }//empty file name check
}

void  AnaFuncs::GoodLumiChecker::addLumiRange_(int runnr,int lowerLumi,int upperLumi)
{
  std::vector<GoodLumiData>::iterator lumiData = std::find(goodLumis_.begin(),goodLumis_.end(),runnr);
  if(lumiData==goodLumis_.end()){ //that run is new, add it
    GoodLumiData data;
    data.runnr = runnr;
    goodLumis_.push_back(data);
    lumiData = goodLumis_.end(); //last element
    --lumiData;
  }
  
  lumiData->allowedLumis.push_back(std::make_pair(lowerLumi,upperLumi));


}

bool AnaFuncs::GoodLumiChecker::isGoodLumiSec(int runnr,int lumiSec)const
{
  std::pair<std::vector<GoodLumiData>::const_iterator,std::vector<GoodLumiData>::const_iterator> lumiData;
  lumiData = std::equal_range(goodLumis_.begin(),goodLumis_.end(),runnr,GoodLumiDataComp());
  // std::cout <<lumiData.first<<" "<<lumiData.second<<" begin "<<goodLumis_.begin()<<" "<<goodLumis_.end()<<std::endl;
  if(lumiData.second-lumiData.first==1){ //found it, now check is lumi sec is good
    // std::cout <<"found run "<<runnr<<std::endl;
    for(size_t lumiNr=0;lumiNr<lumiData.first->allowedLumis.size();lumiNr++){
      // std::cout <<"lumiNr : "<<lumiNr<<" secs "<<lumiData.first->allowedLumis[lumiNr].first<<" - "<<lumiData.first->allowedLumis[lumiNr].second<<" lumi sec "<<lumiSec<<std::endl;
      if(lumiSec>=lumiData.first->allowedLumis[lumiNr].first && lumiSec<=lumiData.first->allowedLumis[lumiNr].second) return true;
    }
    return false;//good run, lumi not in allowed ranges  though
  }else if(lumiData.second-lumiData.first>1){
    std::cout <<"AnaFuncs::GoodLumiChecker:: warning "<<runnr<<" has "<<lumiData.second-lumiData.first<<" entries in lumi list "<<std::endl;
    return false; //lets call it bad...
  }else return false; //not found, bad run

}

bool AnaFuncs::RunList::checkRunnr(int runnr)const
{
  bool found=false;
  for(int i=0;i<(int)_runList.size();i++){
    if(runnr == _runList[i]){
      found = true;
      break;
    }
  }
  return found;
}


void AnaFuncs::LumiData::readLumiData(const std::string& filename)
{
  clear();
  if(!filename.empty()){
    std::ifstream file(filename.c_str());
    if(file.bad()){
      std::cout<<" file "<<filename<<" not found"<<std::endl;
    }else if(file.eof()){
      std::cout<<" file "<<filename<<" found but empty"<<std::endl;
    }
    
    if(!file.bad()){
      while(!file.eof()){

	char tempBuffer[128]; 


	file.getline(tempBuffer,128);
	
	std::pair<int,float> runLumi;

	char* runnrChar = strtok(tempBuffer,"|");
	char* lumiChar = strtok(NULL,"|");
	if(runnrChar!=NULL && lumiChar!=NULL){
	  runLumi.first = atoi(runnrChar);
	  runLumi.second = atof(lumiChar);
	  std::cout <<" run "<<runLumi.first<<" lumi "<<runLumi.second<<std::endl;
	  data_.push_back(runLumi);
	}
      }
    }
  }
  std::sort(data_.begin(),data_.end(),RunComp());
}


float AnaFuncs::LumiData::lumi(int runnr)const
{
  std::pair<std::vector<std::pair<int,float> >::const_iterator,
    std::vector<std::pair<int,float> >::const_iterator> lumiData;
  
  lumiData = std::equal_range(data_.begin(),data_.end(),runnr,RunComp());
  if(lumiData.second-lumiData.first==1) return lumiData.first->second;
  else return 0.;

}

bool AnaFuncs::PUReWeighter::loadWeights()
{
  return loadWeights("input/dataPUBin_22JanReReco.root","input/mcPUDistS10.root");
}

bool AnaFuncs::PUReWeighter::loadWeights(const std::string& dataFilename,const std::string& mcFilename)
{
  weights_.clear();
  TFile* mcFile = new TFile(mcFilename.c_str(),"READ");
  TFile* dataFile = new TFile(dataFilename.c_str(),"READ");
  
  TH1* mcHist = (TH1*) mcFile->Get("mcPUDist");
  TH1* dataHist = (TH1*) dataFile->Get("pileup");
  
  dataHist->Scale(1./dataHist->Integral());
  mcHist->Scale(1./mcHist->Integral());
  nrBins_=mcHist->GetNbinsX();
  xmin_=mcHist->GetBinLowEdge(1);
  xmax_=mcHist->GetBinLowEdge(nrBins_+1);
  for(int binNr=0;binNr<=mcHist->GetNbinsX()+1;binNr++){
    float weight = mcHist->GetBinContent(binNr)!=0 ?  dataHist->GetBinContent(binNr)/mcHist->GetBinContent(binNr) : -1;
    weights_.push_back(weight);
  }
  delete mcFile;
  delete dataFile;
  
  return true;
}

float AnaFuncs::PUReWeighter::weight(int nrPUInts)
{
  size_t binNr = static_cast<size_t>(AnaFuncs::getBinNr(nrBins_,xmin_,xmax_,nrPUInts));
  if(binNr<weights_.size()) return weights_[binNr];
  else return -1;
}
				     

//get the dataset id or MC process
//by my file naming convention, this will be just before the _ntuples part in the filename
//eg dataset bewk0d will be named as bewk0d_ntuples_$CDFSOFTVERSION_$ADDITIONALINFO
void AnaFuncs::getDatasetID(const char* baseFilename,std::string &datasetID)
{
  datasetID = baseFilename;
  int nameBegin = datasetID.find_last_of("/");
  datasetID.erase(0,nameBegin+1); //erasing the directory path from the file name
  int datasetIDEnd = datasetID.find("_");
  datasetID.erase(datasetIDEnd,datasetID.length()-1);//erasing everything but the dataset id from the filename
  if(datasetID=="bewk0d" || datasetID=="zewk0d" || datasetID=="bhel0h" || datasetID=="bpel0h" || datasetID=="bhel0i" || datasetID=="bhel0d" || datasetID=="bhel0x" || datasetID=="bhelbi" || datasetID=="bhel" || datasetID=="bhelah" || datasetID=="bhelai") datasetID = "Data";
    
  datasetID.replace(0,1,1,toupper(datasetID[0]));//captilizing the first letter
}


std::string AnaFuncs::floatToStr(float val)
{
  std::ostringstream valStr;
  valStr<<val;
  return valStr.str();
}

void AnaFuncs::splitStrings(const char* charString,std::vector<std::string>& array,const char* seperator)
{
  array.clear();
  if(strlen(charString)>=10240){
    std::cout <<"AnaFuncs::splitString : Error string \""<<charString<<"\" is larger than 1024 charactors"<<std::endl;
    return;
  }
  char localString[10240];
  strcpy(localString,charString);
  char* chanId = strtok(localString,seperator);
  while(chanId!=NULL){
    array.push_back(chanId);
    chanId = strtok(NULL,seperator);
  }
}

std::string AnaFuncs::makeStringOfBinsLowEdges(int nrBins,float xmin,float xmax)
{
  std::ostringstream returnStr;
  returnStr<<xmin;
  for(int binNr=1;binNr<=nrBins;binNr++){
    returnStr<<":"<<(xmax-xmin)/nrBins*binNr+xmin;
  }
  return returnStr.str();
}

void AnaFuncs::makeVecFromInputString(std::vector<int>& vec,const std::string& inputStr)
{
  vec.clear();
  std::vector<std::string> strVec;
  splitStrings(inputStr.c_str(),strVec);
  
  
  for(size_t i=0;i<strVec.size();i++){
    vec.push_back(atoi(strVec[i].c_str()));
  }
}
void AnaFuncs::makeVecFromInputString(std::vector<float>& vec,const std::string& inputStr)
{
  TempFuncs::makeVec(inputStr,vec);
}
 

void AnaFuncs::makeVecFromInputString(std::vector<double>& vec,const std::string& inputStr)
{
  TempFuncs::makeVec(inputStr,vec);
}
  
void AnaFuncs::copyArrayToVec(const double* array,size_t nrEntries,std::vector<double>& vec)
{
  vec.clear();
  vec.reserve(nrEntries);
  for(size_t i=0;i<nrEntries;i++){
    // std::cout <<"pushing back "<<array[i]<<std::endl;
    vec.push_back(array[i]);
  }
}


int AnaFuncs::getIndexInList(const char* var,const char* listOfVar,const char* seperator)
{
  char *localListOfVar= new char[strlen(listOfVar)+1];
  strcpy(localListOfVar,listOfVar);
  
  int index=0;
  char* varAtIndex = strtok(localListOfVar,seperator);
  while(varAtIndex!=NULL){
    if(strcmp(varAtIndex,var)==0) return index;//found the varible in the list
    index++;
    varAtIndex = strtok(NULL,seperator);
  }
  std::cout <<"AnaFuncs::getIndexInList : Warning "<<var<<" was not found in "<<listOfVar<<" (delimiter = \""<<seperator<<"\")"<<std::endl;
  return 0;
}



void AnaFuncs::Blinder::setBlindPara(int regionCode,int runnr,double mass)
{
  if(regionCode==-1) setBlindPara(runnr,mass);
  else if(regionCode>=0 && regionCode<=2){
    _runnrs[regionCode]=runnr;
    _masses[regionCode]=mass;
  }else{
    std::cout <<"AnaFuncs::Blinder::setBlindPara : Error : regionCode "<<regionCode<<" not valid"<<std::endl;
  }
}

void AnaFuncs::Blinder::setBlindPara(int runnr,double mass)
{
  for(int i=0;i<3;i++){
    _runnrs[i] =runnr;
    _masses[i] =mass;
  }
}

bool AnaFuncs::Blinder::blind(int regionCode,int runnr,double mass)const
{
  if(regionCode==-1) return blind(runnr,mass);
  else if(regionCode>=0 && regionCode<=2){
    return runnr>_runnrs[regionCode] && mass>_masses[regionCode];
  }else{
    std::cout <<"AnaFuncs::Blinder::blind : Error : regionCode "<<regionCode<<" not valid"<<std::endl;
    return true; //well better to be blind than not
  }
}

bool AnaFuncs::Blinder::blind(int runnr,double mass)const
{
  bool result = false;
  for(int i=0;i<3;i++){
    result = result || ( runnr>_runnrs[i] && mass>_masses[i]);
  }
  return result;
}


void AnaFuncs::setBranchAddress(TTree *tree,const char* branchName,void *address)
{
  tree->SetBranchStatus(branchName,1);
  tree->SetBranchAddress(branchName,address);
}

void AnaFuncs::readFilelistFromPattern(const std::string& filelistPattern,std::vector<std::string> &filenames)
{
  char filename[512];
  auto file = popen(("ls "+filelistPattern).c_str(),"r");
  while(file && !std::feof(file)){
    if(std::fgets(filename,sizeof(filename),file)==NULL) break; //otherwise will read the last entry twice
    std::string filenameStr(filename);
    filenameStr.pop_back(); //removing end of line character
    filenames.push_back(filenameStr);
  }
  std::fclose(file);
}

void AnaFuncs::readFilelistFromFile(const std::string& fileListName,std::vector<std::string> &filenames)
{
  std::ifstream fileList(fileListName.c_str());
  char filename[512];
  while(!fileList.eof() && fileList.is_open()){
    fileList.getline(filename,511);
    if(strcmp(filename,"")!=0 && filename[0]!='#'){
      //std::cout <<"filename ="<<filename<<std::endl;
      filenames.push_back(filename);
    }
  }
  fileList.close();
}

void AnaFuncs::readFilelist(std::string fileListName,std::vector<std::string> &filenames,int nrJobs,int jobNr,int verbose)
{
  readFilelist(std::vector<std::string>{fileListName},filenames,nrJobs,jobNr,verbose);
}

void AnaFuncs::readFilelist(std::vector<std::string> fileListNames,std::vector<std::string> &filenames,int nrJobs,int jobNr,int verbose)  
{
  for(const auto& fileListName : fileListNames){
    const char* extension = strstr(fileListName.c_str(),".root");
    if(fileListName.find("*")!=std::string::npos || fileListName.find("?")!=std::string::npos  ){ 
      readFilelistFromPattern(fileListName,filenames);
    }else if(extension!=NULL && strlen(extension)==5){ //file is a root file, process it
      filenames.push_back(fileListName);
    }else{
      readFilelistFromFile(fileListName,filenames);
    }
  }
  
  for(auto& filename : filenames){
    boost::replace_all(filename,"dcap://heplnx209.pp.rl.ac.uk/","root://dcap.pp.rl.ac.uk:1094/");
  }

  if(nrJobs!=1) splitFilelist(nrJobs,jobNr,filenames);
  if(verbose>=2) for(unsigned i=0;i<filenames.size();i++) std::cout <<"filename = "<<filenames[i].c_str()<<std::endl;
  if(verbose) std::cout <<"nr files: "<<filenames.size()<<std::endl; 
}

void AnaFuncs::splitFilelistConsecutive(int nrJobs,int jobNr,std::vector<std::string>& filenames)
{
  if(jobNr>nrJobs || nrJobs<1){
    std::vector<std::string>().swap(filenames);
    return;
  }
  if(jobNr==0) jobNr=nrJobs; //setting 0th job to be the same as last job

  int nrFiles = filenames.size();
  //last job is the remaining files
  //  int nrFilesPerJob=0;
  // if(nrFiles%nrJobs!=0) nrFilesPerJob = nrFiles/(nrJobs-1); //intentionally throwing away remander
  // else nrFilesPerJob =  nrFiles/nrJobs; 
  
  std::vector<std::string>::iterator start,end;
  start = filenames.begin() + nrFilesInPreviousJobs(nrFiles,nrJobs,jobNr);
  end = start+nrFilesInJob(nrFiles,nrJobs,jobNr); // as it should be just past the indx we want
  if(jobNr==nrJobs && end!=filenames.end()) std::cout <<"AnaFuncs::splitFilelist : Warning : last job does not equal end of filenames end"<<std::endl;
  std::vector<std::string> jobFilenames(start,end);
  filenames.swap(jobFilenames);
}

void AnaFuncs::splitFilelist(int nrJobs,int jobNr,std::vector<std::string>& filenames)
{
  splitFilelistMixed(nrJobs,jobNr,filenames);
}

void AnaFuncs::splitFilelistMixed(int nrJobs,int jobNr,std::vector<std::string>& filenames)
{
  if(jobNr>nrJobs || nrJobs<1){
    std::vector<std::string>().swap(filenames);
    return;
  }
  
  std::vector<std::string> outFilenamesTmp;
  for(size_t fileNr=0;fileNr<filenames.size();fileNr++){
    if( (fileNr+jobNr)%nrJobs==0){
      outFilenamesTmp.push_back(filenames[fileNr]);
    }
  }
  filenames.swap(outFilenamesTmp);
}
   
 

int AnaFuncs::nrFilesInJob(int nrFiles,int nrJobs,int jobNr)
{
  int nrFilesLeftOver = nrFiles%nrJobs;
  int nrFilesInJob = nrFiles/nrJobs; //intentionally throwing away remainder
  if(jobNr<= nrFilesLeftOver) nrFilesInJob++;
  return nrFilesInJob;
}

int AnaFuncs::nrFilesInPreviousJobs(int nrFiles,int nrJobs,int jobNr)
{
  int totNrFiles =0;
  for(int job = jobNr-1;job>0;job--){
    totNrFiles+=nrFilesInJob(nrFiles,nrJobs,job);
  }
  return totNrFiles;
}


//god this is horrible but efficiency doesnt matter as its run once per job
void AnaFuncs::addJobNrToFilename(char* filename,int jobNr)
{
  std::string filenameStr(filename);
  addJobNrToFilename(filenameStr,jobNr);
  std::strcpy(filename,filenameStr.c_str());
}

void AnaFuncs::addJobNrToFilename(std::string& filename,int jobNr)
{ 
  std::string::size_type endLoc = filename.find(".root");
  if(endLoc!=std::string::npos){
    filename.erase(endLoc);
    std::ostringstream filenameStream;
    filenameStream << filename<<"_"<<jobNr<<".root";
    
    filename = filenameStream.str();
  }
}
    
//because root doesnt like it so hiding it here..
void AnaFuncs::boost_replace_all(std::string& orgString,const std::string& pattern,const std::string& subTxt)
{
  std::cout <<"function disabled (root *really* didnt like it)"<<std::endl;
  //  boost::algorithm::replace_all(orgString,pattern,subTxt);

}

std::string AnaFuncs::replaceSubStr(const std::string& orgStr,const std::string& subStrToReplace,const std::string& replacement)
{
  std::string newStr(orgStr);

  size_t subStrIndex = newStr.find(subStrToReplace);
  while(subStrIndex!=std::string::npos){
    newStr.erase(subStrIndex,subStrToReplace.size());
    newStr.insert(subStrIndex,replacement);
    subStrIndex = newStr.find(subStrToReplace,subStrIndex);
    
  }
  return newStr;


}


std::string AnaFuncs::convertToTTreeStr(int val) //for -ve vals, string is M instead so -28 is M28
{
  std::ostringstream valStr;
  if(val<0) valStr<<"M";
  valStr<<abs(val);
  return valStr.str();
}


std::string AnaFuncs::convertToTTreeStr(float val) //for -ve vals, string is M instead so -28 is M28
{
  std::ostringstream valOStr;
  valOStr<<val;
  std::string valStr = valOStr.str();
  auto decPoint = valStr.find(".");
  if(decPoint!=std::string::npos){
    valStr.replace(decPoint,1,"p");
  }
  if(val<0) valStr.replace(0,1,"M");
  return valStr;
}

std::string AnaFuncs::convertToTTreeStr(double val) //for -ve vals, string is M instead so -28 is M28
{
  std::ostringstream valOStr;
  valOStr<<val;
  std::string valStr = valOStr.str();
  auto decPoint = valStr.find(".");
  if(decPoint!=std::string::npos){
    valStr.replace(decPoint,1,"p");
  }
  if(val<0) valStr.replace(0,1,"M");
  return valStr;
}

//does not handle scientific notation
bool AnaFuncs::isNumber(const std::string& str)
{
  // if(str.empty()) return false;
  // else{
  //   size_t start = 0;
  //   if(str[0]=='-') start=1;
  //   int nrDecPoints = 0;
  //   for(size_t i=start;i<str.size();i++){
  //     if(str[i]=='.') nrDecPoints++;
  //     else if(!::isdigit(str[i])) return false;
  //     if(nrDecPoints>1) return false;
  //   }
  //   return true;
  // }
  std::stringstream ss(str);
  float f;
  //  return bool(ss>>f);
  return !((ss >> std::noskipws >> f).rdstate() ^ std::ios_base::eofbit);
}
