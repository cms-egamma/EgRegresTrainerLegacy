#include "Utility/EvtIndexLUT.hh"

#include "Utility/LogErr.hh"

bool operator<(const int lhs,const EvtIndexLUT::RunData& rhs)
{
  return lhs<rhs.runnr_;
}

bool operator<(const int lhs,const EvtIndexLUT::LumiData& rhs)
{
  return lhs<rhs.lumiSec_;
}

bool operator<(const int lhs,const EvtIndexLUT::EvtData& rhs)
{
  return lhs<rhs.eventNr;
}

void EvtIndexLUT::add(int runnr,int lumiSec,int eventNr,const std::string& filename,int fileIndex)
{
  isSorted_=false;
 
  nrEvtsStored_++;
 
  EvtData evtData(runnr,lumiSec,eventNr,convertToFileId_(filename),fileIndex);

  if(lastRunW_ &&lastRunW_->runnr()==runnr){
    lastRunW_->add(evtData);
  }else{
    std::vector<RunData>::iterator runIt = std::find(data_.begin(),data_.end(),runnr);
    if(runIt!=data_.end()) lastRunW_=&*runIt;
    else{
      data_.push_back(RunData(runnr));
      lastRunW_= &data_.back();
    }
    lastRunW_->add(evtData);
  } 
}

EvtIndexLUT::LumiData* EvtIndexLUT::RunData::add(const EvtData& evtData)
{
  //  bool print=evtData.lumiSec==128 && evtData.runnr==207454
  nrEvtsStored_++;
  std::vector<LumiData>::iterator lumiIt = std::find(data_.begin(),data_.end(),evtData.lumiSec);
  if(lumiIt==data_.end()){ 
    data_.push_back(LumiData(evtData.lumiSec));
    lumiIt = data_.end()-1; //evil?
  }
  lumiIt->add(evtData);
  return &*lumiIt;
}
  
void EvtIndexLUT::LumiData::sort()
{
  std::sort(data_.begin(),data_.end());
  
  //this collesses multiple occurences of the same run/lumi/eventnr into the same entry
  std::vector<EvtData> tempData;
  if(!data_.empty()) tempData.push_back(data_[0]);
  for(size_t index=1;index<data_.size();index++){
    if(data_[index].eventNr==tempData.back().eventNr){
      tempData.back()+=data_[index];
    }else{
      tempData.push_back(data_[index]);
    }
  }
  data_.swap(tempData);

}

EvtIndexLUT::EvtData& EvtIndexLUT::EvtData::operator+=(EvtIndexLUT::EvtData& rhs)
{
  fileId.insert(fileId.end(),rhs.fileId.begin(),rhs.fileId.end());
  fileIndex.insert(fileIndex.end(),rhs.fileIndex.begin(),rhs.fileIndex.end());
  return *this;
}

void EvtIndexLUT::sort()
{
  std::sort(data_.begin(),data_.end());
  std::for_each(data_.begin(),data_.end(),std::mem_fun_ref(&RunData::sort));
  //std::sort(filenames_.begin(),filenames_.end());
  isSorted_=true;
  lastRunW_=0;
  //lastLumiW_=0;
  lastRunR_=0;
  lastLumiR_=0;
  
}

EvtIndexLUT::EvtIndex EvtIndexLUT::operator[](size_t index)const
{
  if(index>=nrEvtsStored()){
    LogErr<<" Warning, index "<<index<<" out of bounds ("<<nrEvtsStored()<<")"<<std::endl;
    return EvtIndex();
  }
  size_t sum=0;
  for(size_t runNr=0;runNr<data_.size();runNr++){
    sum+=data_[runNr].nrEvtsStored();
    //std::cout <<"sum "<<sum<<" index "<<index<<" nrEvtsStored "<<nrEvtsStored()<<std::endl;
    if(index<sum){
      const EvtData& evtData =  data_[runNr][index-(sum-data_[runNr].nrEvtsStored())];
      if(!evtData.fileId.empty()){
	EvtIndex evtIndex(evtData.runnr,evtData.lumiSec,evtData.eventNr,convertToFilename_(evtData.fileId[0]),evtData.fileIndex[0]);
	return evtIndex;
      }else{
	LogErr<<" Warning, somehow an entry with no indices stored "<<std::endl;
	return EvtIndex();
      }
    }
  }
  LogErr<<" Warning, logic fail, it should be impossible to actually get here "<<index<<" "<<nrEvtsStored()<<std::endl;
  return EvtIndex();
}

const EvtIndexLUT::EvtData& EvtIndexLUT::RunData::operator[](size_t index)const
{
  
  if(index>=nrEvtsStored()){
    LogErr<<" Error, index "<<index<<" out of bounds ("<<nrEvtsStored()<<"), this will result in memory coruption so terminate"<<std::endl;
    exit(1);
  }

  size_t sum=0;
  for(size_t lumiNr=0;lumiNr<data_.size();lumiNr++){
    sum+=data_[lumiNr].nrEvtsStored();
    if(index<sum) return data_[lumiNr][index-(sum-data_[lumiNr].nrEvtsStored())];  
  }
  LogErr<<" Warning, logic fail, it should be impossible to actually get here "<<data_.size()<<" index "<<index<<" sum "<<sum<<" nrEvtsStored "<<nrEvtsStored()<<" runnr "<<runnr_<<std::endl;
  exit(1);
}

const EvtIndexLUT::EvtData& EvtIndexLUT::LumiData::operator[](size_t index)const
{
  if(index>=nrEvtsStored()){
    LogErr<<" Error, index "<<index<<" out of bounds ("<<nrEvtsStored()<<"), this will result in memory coruption"<<std::endl;
    exit(1);
  }
  return data_[index];
}


EvtIndexLUT::EvtIndex EvtIndexLUT::get(int runnr,int lumiSec,int eventNr)const
{
  if(!isSorted_){
     LogErr<<" Warnining, LUT is not sorted and is const so events can not be looked up "<<std::endl;
     return EvtIndex();
  }

  const EvtData* evtData=getEvtData_(runnr,lumiSec,eventNr);
  
  if(evtData){ 
    if(!evtData->fileId.empty()){
      EvtIndex evtIndex(evtData->runnr,evtData->lumiSec,evtData->eventNr,convertToFilename_(evtData->fileId[0]),evtData->fileIndex[0]);
      return evtIndex;
    }else{
      LogErr<<" Warnining, event "<<runnr<<" "<<lumiSec<<" "<<eventNr<<" found but empty index"<<std::endl;
      return EvtIndex();
    }
  }else{
    if(verbose_==1) LogErr<<" Warnining, event "<<runnr<<" "<<lumiSec<<" "<<eventNr<<" not found"<<std::endl;
    return EvtIndex();
  }

}

std::vector< std::pair<std::string,int> > EvtIndexLUT::getFilenameAndIndexMulti(int runnr,int lumiSec,int eventNr)const
{
  std::vector<std::pair<std::string,int > > retVal;
  const EvtData* evtData=getEvtData_(runnr,lumiSec,eventNr);

  if(evtData){
    for(size_t entryNr=0;entryNr<evtData->fileId.size();entryNr++){
      retVal.push_back(std::pair<std::string,int>(convertToFilename_(evtData->fileId[entryNr]),evtData->fileIndex[entryNr]));
    }
    if(retVal.empty()){
      LogErr<<" Warnining, event "<<runnr<<" "<<lumiSec<<" "<<eventNr<<" found but has zero entries"<<std::endl;
      
      retVal.push_back(std::pair<std::string,int>("",-1));
    }
  }else{
    LogErr<<" Warnining, event "<<runnr<<" "<<lumiSec<<" "<<eventNr<<" not found"<<std::endl;
    
    retVal.push_back(std::pair<std::string,int>("",-1));
  }
  return retVal;
}

const EvtIndexLUT::EvtData* EvtIndexLUT::getEvtData_(int runnr,int lumiSec,int eventNr)const
{

  if(lastRunR_ && lastRunR_->runnr()==runnr){
    if(!lastLumiR_ || lastLumiR_->lumiSec()!=lumiSec) lastLumiR_ = lastRunR_->getLumi(lumiSec);
  }else{
    lastRunR_ = TempFuncs::findSingleSorted<std::vector<RunData>,RunData,int>(data_,runnr);
    lastLumiR_ = lastRunR_ ? lastRunR_->getLumi(lumiSec) : 0;
  }

  if(lastLumiR_) return lastLumiR_->getEvtData(eventNr);
  else return 0;
}

size_t EvtIndexLUT::convertToFileId_(const std::string& filename)
{
//   if(isSorted_){
//     std::vector<std::string>::iterator It;
//     std::pair<It,It> result = std::equal_range(filenames_.begin(),filenames_.end(),filename);
    
//     if(result.second-result.first==0) filenames_.insert(result.first,filename);
//     if(result.second-result.first>1){
//       LogErr<<" warning, "<<result.second-result.first<<" copies of "<<filename<<" in file id index, bad things are going to happen"<<std::endl;
//     }
//     return result.first-filenames_.begin();
//   }

//   }else{
  std::vector<std::string>::iterator fileIt= std::find(filenames_.begin(),filenames_.end(),filename);
  if(fileIt!=filenames_.end()) return fileIt-filenames_.begin();
  else{
    filenames_.push_back(filename);
    return filenames_.size()-1;
  }
  
}

const std::string& EvtIndexLUT::convertToFilename_(size_t fileId)const
{
  if(fileId>=filenames_.size()) return filenames_[0]; //the null file id
  else return filenames_[fileId];
}
