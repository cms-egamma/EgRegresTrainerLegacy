#ifndef EVTLIST_HH
#define EVTLIST_HH

#include <vector>
#include <algorithm>
#include <iostream>
#include <limits>
class EvtList {
public:
  struct LumiData {
    int lumiSec;
    std::vector<int> events;

    LumiData(int iLumiSec):lumiSec(iLumiSec){}
    void sort(){std::sort(events.begin(),events.end());}
    bool operator<(const LumiData& rhs)const{return lumiSec<rhs.lumiSec;}
    bool operator<(int rhs)const{return lumiSec<rhs;} 
    bool operator==(const LumiData& rhs)const{return lumiSec==rhs.lumiSec;}
    bool operator==(int rhs)const{return lumiSec==rhs;}
    friend bool operator<(int lhs,const LumiData& rhs);
  };
  
  struct RunData {
    int runnr;
    std::vector<LumiData> lumis;
    
    RunData(int iRunnr):runnr(iRunnr){}
    void sort();
    bool operator<(const RunData& rhs)const{return runnr<rhs.runnr;}
    bool operator<(int rhs)const{return runnr<rhs;}
    bool operator==(const RunData& rhs)const{return runnr==rhs.runnr;}
    bool operator==(int rhs)const{return runnr==rhs;}
    friend bool operator<(int lhs,const RunData& rhs);
  };
  
private:
  std::vector<RunData> runs_;
  
  mutable size_t lastRunIndex_;
  mutable size_t lastLumiIndex_;

public:


  EvtList():lastRunIndex_(std::numeric_limits<size_t>::max()),lastLumiIndex_(std::numeric_limits<size_t>::max()){}
  ~EvtList(){}

  void add(int runnr,int lumiSec,int eventnr);
  void sort();
  void read(std::istream& stream);
  
  bool hasEvt(int runnr,int lumiSec,int eventnr)const;

private:
  const RunData* lastRun_()const{return lastRunIndex_<runs_.size() ? &runs_[lastRunIndex_] : 0;}
  RunData* lastRun_(){return lastRunIndex_<runs_.size() ? &runs_[lastRunIndex_] : 0;}
  
  const LumiData* lastLumi_()const{return lastRun_() && lastLumiIndex_<lastRun_()->lumis.size() ? &(lastRun_()->lumis[lastLumiIndex_]) : 0;}
  LumiData* lastLumi_(){return lastRun_() && lastLumiIndex_<lastRun_()->lumis.size() ? &(lastRun_()->lumis[lastLumiIndex_]) : 0;}


};



#endif
