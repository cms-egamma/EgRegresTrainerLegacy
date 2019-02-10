#ifndef UTILITY_EVTLUT_HH
#define UTILITY_EVTLUT_HH

//a lookup table which stores the filename and index  in that file of a given event

#include "Utility/TempFuncs.hh"
#include "Utility/LogErr.hh"
#include <vector>
#include <string>
#include <algorithm>
#include <numeric>

template<typename Payload>
class EvtLUT {
  ClassDef(EvtLUT,1)
  //this is the information we wish to store, the class is designed solely to store and retrieve this info
  //aim: store order 1-10 million events
  //use case: fill LUT and then read back. Expect filled once and then read back a lot
   

public:

  struct EvtData {
    int runnr; //not strickly necessary but makes things so much damn easier
    int lumiSec; //not strickly necessary but makes things so much damn easier
    int eventNr; //the key, only this is used for sorting (EvtData of different runs and lumisecs should not be currently stored in the same container)
    Payload payload;

    EvtData(int iRunnr,int iLumiSec,int iEventNr,const Payload& iPayload):
      runnr(iRunnr),lumiSec(iLumiSec),eventNr(iEventNr),payload(iPayload){}
    EvtData():runnr(0),lumiSec(0),eventNr(0),payload(){}
    bool operator<(const EvtData& rhs)const{return eventNr<rhs.eventNr;}
    bool operator<(const int rhs)const{return eventNr<rhs;}
    friend bool operator<(const int lhs,const typename EvtLUT<Payload>::EvtData& rhs){return lhs<rhs.eventNr;}
    EvtData& operator+=(EvtData& rhs);
  };

  class LumiData {
  private:
    int lumiSec_;
    std::vector<EvtData> data_;
  public:
    LumiData(const int iLumiSec=0):lumiSec_(iLumiSec){}
    void add(const EvtData& evtData){data_.push_back(evtData);}
    bool operator==(const LumiData& rhs)const{return lumiSec_==rhs.lumiSec_;}
    bool operator==(const int rhs)const{return lumiSec_==rhs;}
    bool operator<(const LumiData& rhs)const{return lumiSec_<rhs.lumiSec_;}
    bool operator<(const int rhs)const{return lumiSec_<rhs;}
    friend bool operator<(const int lhs,const LumiData& rhs){return lhs<rhs.lumiSec();}
    const EvtData& operator[](size_t index)const;
    int lumiSec()const{return lumiSec_;}
    size_t nrEvtsStored()const{return data_.size();}
    void sort();
    const EvtData* getEvtData(int eventNr)const{return TempFuncs::findSingleSorted<std::vector<EvtData>,EvtData,int>(data_,eventNr);}
  };
  
  class RunData {
  private:
    int runnr_;
    size_t nrEvtsStored_; //this is invalid when we collapse entries together in the sort so right now this is ignored
    std::vector<LumiData> data_;
  public:
    RunData(const int iRunnr=0):runnr_(iRunnr),nrEvtsStored_(0){}
    LumiData* add(const EvtData& evtData);
    bool operator==(const RunData& rhs)const{return runnr_==rhs.runnr_;}
    bool operator==(const int rhs)const{return runnr_==rhs;}
    bool operator<(const RunData& rhs)const{return runnr_<rhs.runnr_;} 
    bool operator<(const int rhs)const{return runnr_<rhs;}
    friend bool operator<(const int lhs,const typename EvtLUT<Payload>::RunData& rhs){return lhs<rhs.runnr();}

    int runnr()const{return runnr_;}
    //    size_t nrEvtsStored()const{return nrEvtsStored_;}
    
    size_t nrEvtsStored()const{return std::accumulate(data_.begin(),data_.end(),0,
						      [](size_t lhs,const LumiData& rhs){
							return lhs+rhs.nrEvtsStored();});
    }
      
    const EvtData& operator[](size_t index)const;
    
    void sort(){std::sort(data_.begin(),data_.end());std::for_each(data_.begin(),data_.end(),std::mem_fun_ref(&LumiData::sort));}
    const LumiData* getLumi(int lumiSec)const{return TempFuncs::findSingleSorted<std::vector<LumiData>,LumiData,int>(data_,lumiSec);}
  };
private:
  std::vector<RunData> data_; //will be sorted before read back

  bool isSorted_; //checks if sorted, will auto sort if necessary
  //temporary cache to store the last run and lumi sec (assume we're running over events  more or  less in order)
  size_t nrEvtsStored_;//this is invalid when we collapse entries together in the sort so right now this is ignored

  mutable RunData* lastRunW_; //!
  // mutable LumiData* lastLumiW_; //!
  mutable const RunData* lastRunR_; //!
  mutable const LumiData* lastLumiR_; //! 
  mutable int verbose_; //!


public:
  EvtLUT():nrEvtsStored_(0),lastRunW_(0),lastRunR_(0),lastLumiR_(0),verbose_(0){}
  virtual ~EvtLUT(){}
  
  void add(const EvtData& evt);
  void sort();

  void setVerboseLvl(int val)const{verbose_=val;}
 
  //size_t nrEvtsStored()const{return nrEvtsStored_;}
  size_t nrEvtsStored()const{return std::accumulate(data_.begin(),data_.end(),0,
						    [](size_t lhs,const RunData& rhs){
						      return lhs+rhs.nrEvtsStored();});
  }
  const EvtData& operator[](size_t index)const;
  const EvtData* get(int runnr,int lumiSec,int eventNr)const;
  const EvtData* get(int runnr,int lumiSec,int eventNr){sort();return static_cast<const EvtLUT*>(this)->get(runnr,lumiSec,eventNr);}
  
  const EvtData* getEvtData_(int runnr,int lumiSec,int eventNr)const;
};

templateClassImp(EvtLUT)


template<typename T>
void EvtLUT<T>::add(const EvtData& evtData)
{
  isSorted_=false;
 
  nrEvtsStored_++;

  if(lastRunW_ &&lastRunW_->runnr()==evtData.runnr){
    lastRunW_->add(evtData);
  }else{
    auto runIt = std::find(data_.begin(),data_.end(),evtData.runnr);
    if(runIt!=data_.end()) lastRunW_=&*runIt;
    else{
      data_.push_back(RunData(evtData.runnr));
      lastRunW_= &data_.back();
    }
    lastRunW_->add(evtData);
  } 
}

template<typename T>
typename EvtLUT<T>::LumiData* EvtLUT<T>::RunData::add(const EvtData& evtData)
{
  //  bool print=evtData.lumiSec==128 && evtData.runnr==207454
  nrEvtsStored_++;
  if(!data_.empty() && data_.back().lumiSec()==evtData.lumiSec){
    data_.back().add(evtData);
    return &data_.back();
  }else{
    auto lumiIt = std::find(data_.begin(),data_.end(),evtData.lumiSec);
    if(lumiIt==data_.end()){
      data_.push_back(LumiData(evtData.lumiSec));
      lumiIt = data_.end()-1; //evil?
    }
    lumiIt->add(evtData);
    return &*lumiIt;
  }
}
 
template<typename T>
void EvtLUT<T>::LumiData::sort()
{
  std::sort(data_.begin(),data_.end());
  
  //this collesses multiple occurences of the same run/lumi/eventnr into the same entry
  // std::vector<EvtData> tempData;
  // if(!data_.empty()) tempData.push_back(data_[0]);
  // for(size_t index=1;index<data_.size();index++){
  //   if(data_[index].eventNr==tempData.back().eventNr){
  //     tempData.back()+=data_[index];
  //   }else{
  //     tempData.push_back(data_[index]);
  //   }
  // }
  //  data_.swap(tempData);

}

template<typename T>
void EvtLUT<T>::sort()
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

template<typename T>
const typename EvtLUT<T>::EvtData& EvtLUT<T>::operator[](size_t index)const
{
  if(index>=nrEvtsStored()){
    LogErr<<" Warning, index "<<index<<" out of bounds ("<<nrEvtsStored()<<")"<<std::endl;
    exit(0);
  }
  size_t sum=0;
  for(size_t runNr=0;runNr<data_.size();runNr++){
    sum+=data_[runNr].nrEvtsStored();
    //std::cout <<"sum "<<sum<<" index "<<index<<" nrEvtsStored "<<nrEvtsStored()<<std::endl;
    if(index<sum){
      const EvtData& evtData =  data_[runNr][index-(sum-data_[runNr].nrEvtsStored())];
      return evtData;
    }
  }
  LogErr<<" Warning, logic fail, it should be impossible to actually get here "<<index<<" "<<nrEvtsStored()<<std::endl;
  exit(0);
}

template<typename T>
const typename EvtLUT<T>::EvtData& EvtLUT<T>::RunData::operator[](size_t index)const
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

template<typename T>
const typename EvtLUT<T>::EvtData& EvtLUT<T>::LumiData::operator[](size_t index)const
{
  if(index>=nrEvtsStored()){
    LogErr<<" Error, index "<<index<<" out of bounds ("<<nrEvtsStored()<<"), this will result in memory coruption"<<std::endl;
    exit(1);
  }
  return data_[index];
}

template<typename T>
const typename EvtLUT<T>::EvtData* EvtLUT<T>::get(int runnr,int lumiSec,int eventNr)const
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



#endif
