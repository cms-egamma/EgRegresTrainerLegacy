#ifndef EVENTINDEXLUT_HH
#define EVENTINDEXLUT_HH

//a lookup table which stores the filename and index  in that file of a given event

#include "Utility/TempFuncs.hh"

#include <vector>
#include <string>
#include <algorithm>
#include <numeric>

class EvtIndexLUT {
public:
  
  //this is the information we wish to store, the class is designed solely to store and retrieve this info
  //aim: store order 1-10 million events
  //use case: fill LUT and then read back. Expect filled once and then read back a lot
  struct EvtIndex{
    int runnr;
    int lumiSec;
    int eventNr;
    std::string filename;
    int  fileIndex;

    EvtIndex(int iRunnr,int iLumiSec,int iEventNr,const std::string& iFilename,int iFileIndex):
      runnr(iRunnr),lumiSec(iLumiSec),eventNr(iEventNr),filename(iFilename),fileIndex(iFileIndex){}
    EvtIndex():runnr(0),lumiSec(0),eventNr(0),filename(),fileIndex(){}
  };

public:
  //these structs internally store the required data, split for faster reading. Filename is not stored per event but in a table as filename is large and will be max ~500-2000 so no point saving it for every event

  struct EvtData {
    int runnr; //not strickly necessary but makes things so much damn easier
    int lumiSec; //not strickly necessary but makes things so much damn easier
    int eventNr; //the key, only this is used for sorting (EvtData of different runs and lumisecs should not be currently stored in the same container)
    std::vector<size_t> fileId;
    std::vector<int> fileIndex;
    EvtData(int iRunnr,int iLumiSec,int iEventNr,size_t iFileId,int iFileIndex):
      runnr(iRunnr),lumiSec(iLumiSec),eventNr(iEventNr),fileId(1,iFileId),fileIndex(1,iFileIndex){}
    EvtData():runnr(0),lumiSec(0),eventNr(0),fileId(),fileIndex(){}
    bool operator<(const EvtData& rhs)const{return eventNr<rhs.eventNr;}
    bool operator<(const int rhs)const{return eventNr<rhs;}
    friend bool operator<(const int lhs,const EvtData& rhs);
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
    friend bool operator<(const int lhs,const LumiData& rhs);
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
    friend bool operator<(const int lhs,const RunData& rhs);

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
  friend bool operator<(const int lhs,const RunData& rhs);
  friend bool operator<(const int lhs,const LumiData& rhs);
  friend bool operator<(const int lhs,const EvtData& rhs);
private:
  std::vector<RunData> data_; //will be sorted before read back
  std::vector<std::string> filenames_; //will never be sorted, this makes it slow to look up fileId but fast to look up filename as its at index fileId, when writing, we will be dealing with 1-2 files while when reading hundreds so I think we are fine here
  bool isSorted_; //checks if sorted, will auto sort if necessary
  //temporary cache to store the last run and lumi sec (assume we're running over events  more or  less in order)
  size_t nrEvtsStored_;//this is invalid when we collapse entries together in the sort so right now this is ignored

  mutable RunData* lastRunW_; //!
  // mutable LumiData* lastLumiW_; //!
  mutable const RunData* lastRunR_; //!
  mutable const LumiData* lastLumiR_; //! 
  mutable int verbose_; //!


public:
  EvtIndexLUT():nrEvtsStored_(0),lastRunW_(0),lastRunR_(0),lastLumiR_(0),verbose_(0){filenames_.push_back("");}
  ~EvtIndexLUT(){}
  
  void add(const EvtIndex& evt){add(evt.runnr,evt.lumiSec,evt.eventNr,evt.filename,evt.fileIndex);}
  void add(int runnr,int lumiSec,int eventNr,const std::string& filename,int fileIndex);
  void sort();

  void setVerboseLvl(int val)const{verbose_=val;}

  size_t nrFiles()const{return filenames_.size();}
  const std::string& filename(size_t index)const{return filenames_[index];}
 
  //size_t nrEvtsStored()const{return nrEvtsStored_;}
  size_t nrEvtsStored()const{return std::accumulate(data_.begin(),data_.end(),0,
						    [](size_t lhs,const RunData& rhs){
						      return lhs+rhs.nrEvtsStored();});
  }
  EvtIndex operator[](size_t index)const;
  EvtIndex get(int runnr,int lumiSec,int eventNr)const;
  EvtIndex get(int runnr,int lumiSec,int eventNr){sort();return static_cast<const EvtIndexLUT*>(this)->get(runnr,lumiSec,eventNr);}

  std::vector<std::pair<std::string,int> > getFilenameAndIndexMulti(int runnr,int lumiSec,int eventNr)const;
  std::vector<std::pair<std::string,int> > getFilenameAndIndexMulti(int runnr,int lumiSec,int eventNr){sort();return static_cast<const EvtIndexLUT*>(this)->getFilenameAndIndexMulti(runnr,lumiSec,eventNr);}
  std::pair<std::string,int> getFilenameAndIndex(int runnr,int lumiSec,int eventNr)const{return getFilenameAndIndexMulti(runnr,lumiSec,eventNr)[0];}
  std::pair<std::string,int> getFilenameAndIndex(int runnr,int lumiSec,int eventNr){sort();return static_cast<const EvtIndexLUT*>(this)->getFilenameAndIndex(runnr,lumiSec,eventNr);}
private:
  size_t convertToFileId_(const std::string& filename);
  const std::string& convertToFilename_(size_t fileId)const;
  const EvtData* getEvtData_(int runnr,int lumiSec,int eventNr)const;
};

#endif
