#include "Utility/GoodLumiChecker.hh"

#include <cstring>
#include <fstream>
#include <iostream>
#include <algorithm>

GoodLumiChecker::GoodLumiChecker(const std::string& goodLumiList)
{
  setGoodLumiList(goodLumiList);

}


void GoodLumiChecker::setGoodLumiList(const std::string& goodLumiList)
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

	  if(run1!=run2) std::cout <<"GoodLumiChecker::setGoodLumiList: Error run number missmatch "<<run1<<" : "<<run2<<std::endl;
	  
	  // std::cout <<" \'"<<run1<<":"<<lumi1<<"-"<<run2<<":"<<lumi2<<"\',"<<std::endl;
	  addLumiRange_(run1,lumi1,lumi2);
	  
	}//end runAndLumiSec is valid
      }//end of file
      std::sort(goodLumis_.begin(),goodLumis_.end());
    }//bad file check
    
    
  }//empty file name check
}

void  GoodLumiChecker::addLumiRange_(int runnr,int lowerLumi,int upperLumi)
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

bool GoodLumiChecker::isGoodLumiSec(int runnr,int lumiSec)const
{
  std::pair<std::vector<GoodLumiData>::const_iterator,std::vector<GoodLumiData>::const_iterator> lumiData;
  lumiData = std::equal_range(goodLumis_.begin(),goodLumis_.end(),runnr,GoodLumiDataComp());
  if(lumiData.second-lumiData.first==1){ //found it, now check is lumi sec is good
    //just brute force it, run over all lumi range and see if its in an allowed range
    for(size_t lumiNr=0;lumiNr<lumiData.first->allowedLumis.size();lumiNr++){
      if(lumiSec>=lumiData.first->allowedLumis[lumiNr].first && lumiSec<=lumiData.first->allowedLumis[lumiNr].second) return true;
    }
    return false;//good run, lumi not in allowed ranges  though
  }else if(lumiData.second-lumiData.first>1){
    std::cout <<"GoodLumiChecker:: warning "<<runnr<<" has "<<lumiData.second-lumiData.first<<" entries in lumi list "<<std::endl;
    return false; //lets call it bad...
  }else return false; //not found, bad run

}
