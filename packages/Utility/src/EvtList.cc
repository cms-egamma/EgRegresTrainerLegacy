
#include "Utility/EvtList.hh"

#include <iostream>
#include <sstream>

template <typename TCont,typename TVal> const size_t findSingleSortedIndex(const TCont& container,const TVal& val)
{
  // typename TCont::const_iterator;
  std::pair<typename TCont::const_iterator,typename TCont::const_iterator> result = std::equal_range(container.begin(),container.end(),val);
  int nrFound = std::distance(result.first,result.second);
  if(nrFound==1) return std::distance(container.begin(),result.first);
  else if(nrFound>=1) {
    std::cout <<"warning "<<nrFound<<" instead of 1 for "<<val<<std::endl;
    return std::distance(container.begin(),result.first);
  }
  return container.size();
}
template <typename TCont,typename TVal> const size_t findSingleIndex(const TCont& container,const TVal& val)
{
  //typename TCont::const_iterator;
  typename TCont::const_iterator result = std::find(container.begin(),container.end(),val);
  return std::distance(container.begin(),result);
}
bool operator<(int lhs,const EvtList::RunData& rhs){return lhs<rhs.runnr;}
bool operator<(int lhs,const EvtList::LumiData& rhs){return lhs<rhs.lumiSec;}
bool operator==(int lhs,const EvtList::RunData& rhs){return lhs==rhs.runnr;}
bool operator==(int lhs,const EvtList::LumiData& rhs){return lhs==rhs.lumiSec;}

bool EvtList::hasEvt(int runnr,int lumiSec,int eventnr)const
{ 
  //std::cout <<"runnr "<<runnr<<" lumiSec "<<lumiSec<<" eventnr "<<eventnr<<" lastRun "<<lastRun_()<<" "<<lastLumi_()<<" "<<lastRunIndex_<<" "<<lastLumiIndex_<<std::endl;
  if(lastRun_() && lastRun_()->runnr==runnr){
    if(!lastLumi_() || lastLumi_()->lumiSec!=lumiSec){
      lastLumiIndex_ =  findSingleSortedIndex(lastRun_()->lumis,lumiSec);
    }
  }else{
    
    lastRunIndex_ = findSingleSortedIndex(runs_,runnr);
    lastLumiIndex_ = lastRun_() ? findSingleSortedIndex(lastRun_()->lumis,lumiSec) : std::numeric_limits<size_t>::max(); 
    //  std::cout <<" now testing runnr "<<runnr<<" lumiSec "<<lumiSec<<" eventnr "<<eventnr<<" lastRun "<<lastRun_()<<" "<<lastLumi_()<<" "<<lastRunIndex_<<" "<<lastLumiIndex_<<std::endl;
  }
  
  // std::cout <<"runnr "<<runnr<<" lumiSec "<<lumiSec<<" eventnr "<<eventnr<<" lastRun "<<lastRun_()<<" "<<lastLumi_()<<" pass "<<findSingleSortedIndex(lastLumi_()->events,eventnr)<<" "<<lastLumi_()->events.size()<<std::endl;

  if(!lastRun_() || !lastLumi_()) return false;
  else return findSingleSortedIndex(lastLumi_()->events,eventnr)<lastLumi_()->events.size();
  
}

void EvtList::add(int runnr,int lumiSec,int eventnr)
{
  if(!lastRun_() || lastRun_()->runnr!=runnr){ 
    lastLumiIndex_=std::numeric_limits<size_t>::max(); //lumisec is from a different run, is invalid
    //std::cout <<"last run index "<<lastRunIndex_<<" runnr "<<runnr<<std::endl;
    lastRunIndex_ = findSingleIndex(runs_,runnr);
    //std::cout <<"updated "<<lastRunIndex_<<std::endl;
    if(!lastRun_()){  
      // std::cout <<"adding run "<<runnr<<" "<<lastRun_()<<" lastRunIndex "<<lastRunIndex_<<std::endl;
      lastRunIndex_ =runs_.size();
      runs_.push_back(runnr);

    }
  }
  
  
  if(!lastLumi_() || lastLumi_()->lumiSec!=lumiSec){
    lastLumiIndex_ =  findSingleIndex(lastRun_()->lumis,lumiSec);
  } 
  
  if(!lastLumi_()){
    lastLumiIndex_ = lastRun_()->lumis.size();
    lastRun_()->lumis.push_back(lumiSec);
    
  }
 
  lastLumi_()->events.push_back(eventnr);
}
      
void EvtList::sort()
{
  lastRunIndex_=std::numeric_limits<size_t>::max();
  lastLumiIndex_=std::numeric_limits<size_t>::max();
  std::sort(runs_.begin(),runs_.end());
  std::for_each(runs_.begin(),runs_.end(),std::mem_fun_ref(&EvtList::RunData::sort));
  
}

void EvtList::RunData::sort()
{
  std::sort(lumis.begin(),lumis.end());
  std::for_each(lumis.begin(),lumis.end(),std::mem_fun_ref(&EvtList::LumiData::sort));
}
void EvtList::read(std::istream& stream)
{ 
  std::string line;
  //  std::cout <<"reading "<<std::endl;
  while(std::getline(stream,line)){
    // std::cout <<"line "<<line<<std::endl;
    line.erase( std::find( line.begin(), line.end(), '#' ), line.end() ); //ignore comments
    std::istringstream lineStream(line);
    int runnr,lumiSec,eventnr;
    while(lineStream >> runnr >> lumiSec >> eventnr ){
      //std::cout <<"pusing back "<<runnr<<" "<<lumiSec<<" "<<eventnr<<std::endl;
      add(runnr,lumiSec,eventnr);
    }
  }
  sort();  
}
