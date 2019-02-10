#include "Utility/EventListComp.hh"

//ClassImp(EventListComp)

#include <algorithm>
#include <cstdlib>
#include <cstring>

void EventListComp::compareLists(const char* filenameA1,const char* filenameB2)
{
  readList_(filenameA1,eventListA1_);
  readList_(filenameB2,eventListB2_);

  int nrNotFoundInB2=0,nrNotFoundInA1=0;
  //check all entries in list A1 are in list B2 and flag any that are not
  for(size_t entryNr=0;entryNr<eventListA1_.size();entryNr++){
    if(!std::binary_search(eventListB2_.begin(),eventListB2_.end(),eventListA1_[entryNr])){
      std::cout <<"event "<<filenameA1<<" not found in "<<filenameB2<<" "<<eventListA1_[entryNr]<<std::endl;
      nrNotFoundInB2++;
    }
  }
  std::cout <<nrNotFoundInB2<<"/"<<eventListA1_.size()<<" events in "<<filenameA1<<" not in "<<filenameB2<<std::endl;
  //now check all entries in list B2 are in list A1
  for(size_t entryNr=0;entryNr<eventListB2_.size();entryNr++){
    if(!std::binary_search(eventListA1_.begin(),eventListA1_.end(),eventListB2_[entryNr])){
      std::cout <<"event "<<filenameB2<<" not found in "<<filenameA1<<" "<<eventListB2_[entryNr]<<std::endl;
      nrNotFoundInA1++;
    }
  }  
  std::cout <<nrNotFoundInA1<<"/"<<eventListB2_.size()<<" events in "<<filenameB2<<" not in "<<filenameA1<<std::endl;

}

void EventListComp::readList_(const char* filename,std::vector<EventListCompData>& eventList)
{
  eventList.clear();
  
  std::ifstream file(filename);
  if(file.bad()) std::cout <<"file "<<filename<<" not found"<<std::endl;
  else if(file.eof()) std::cout <<"file "<<filename<<" found but empty"<<std::endl;
  else std::cout <<"Opened file "<<filename<<std::endl;

  if(!file.bad()){
    while(!file.eof()){
      char tempBuffer[128];
      file.getline(tempBuffer,128);
      char* runnrChar = strtok(tempBuffer," ");
      if(runnrChar!=NULL){
	int runnr = atoi(runnrChar);
	char* lumiChar = strtok(NULL," ");
	int lumiSec = atoi(lumiChar);
	char* eventnrChar = strtok(NULL," ");
	int eventnr = atoi(eventnrChar);
	char* userVarChar = strtok(NULL," ");
	double userVar=0.;
	if(userVarChar!=NULL) userVar = atof(userVarChar);
	EventListCompData data;
	data.runnr= runnr;
	data.lumiSec = lumiSec;
	data.eventnr= eventnr;
	data.userVar= userVar;
	//	std::cout <<"read "<<data<<std::endl;
	eventList.push_back(data);
      }//end of empty line check
    } //end end of file check
  }//end file good check
  std::sort(eventList.begin(),eventList.end());
  // // std::cout <<"sorted list"<<std::endl;
  // for(size_t i=0;i<eventList.size();i++) std::cout << eventList[i]<<std::endl;
}


bool EventListComp::EventListCompData::operator<(const EventListCompData& rhs)const
{
  if(runnr<rhs.runnr) return true;
  else if(runnr>rhs.runnr) return false;
  else{ //run numbers are equal, lumi nr for tye break    
    if(lumiSec<rhs.lumiSec) return true;
    else if(lumiSec>rhs.lumiSec) return false;
    else{
      if(eventnr<rhs.eventnr) return true;
      else if(eventnr>rhs.eventnr) return false;
      else { //equal eventnrs, user var tie break
	if(userVar<rhs.userVar) return true;
	else return false;
      }
    }
  }
}
  

std::ostream& operator <<(std::ostream& output,const EventListComp::EventListCompData &data)
{
  output << "runnr "<<data.runnr<<" lumiSec "<<data.lumiSec<<" eventnr "<<data.eventnr<<" userVar "<<data.userVar;
  return output;
}
