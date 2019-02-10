#include "Utility/BXPUInfo.hh"
#include "Utility/LogErr.hh"

ClassImp(BXPUInfo)

float BXPUInfo::exptNrPUInt(int runnr,int lumiSec,int bx)const
{
  if(bx>kNrBX){
    LogErr<<" for run "<<runnr<<" lumi "<<lumiSec<<" bx "<<bx<<" is greater than max bx of "<<kNrBX<<std::endl;
    return 0;
  }
  auto bxData = getBXData(runnr,lumiSec);
  if(bxData) return bxData->get(bx);
  else return 0;
}

std::pair<int,int> BXPUInfo::nrInTrain(int runnr,int lumiSec,int bx)const
{
  if(bx>kNrBX){
    LogErr<<" for run "<<runnr<<" lumi "<<lumiSec<<" bx "<<bx<<" is greater than max bx of "<<kNrBX<<std::endl;
    return {-1,-1};
  }
  auto bxData = getBXData(runnr,lumiSec);
  if(bxData && bxData->get(bx)>0){
    size_t bxStartTrain = bx;
    while( bxStartTrain>=1 && bxData->get(bxStartTrain-1)>0) bxStartTrain--;
    size_t bxEndTrain = bx;
    while( bxEndTrain+1<bxData->size() && bxData->get(bxEndTrain+1)>0) bxEndTrain++;
    
    return {bx-bxStartTrain,bxEndTrain-bx};
  }else return {-1,-1};
}

int BXPUInfo::nrBXFromStartOfTrain(int runnr,int lumiSec,int bx)const
{
  return nrInTrain(runnr,lumiSec,bx).first;
}

int BXPUInfo::nrBXFromEndOfTrain(int runnr,int lumiSec,int bx)const
{
  return nrInTrain(runnr,lumiSec,bx).second;
}

const BXPUInfo::BXData* BXPUInfo::getBXData(int runnr,int lumiSec)const
{
  auto res = runData_.find(runnr);
  if(res==runData_.end()) return nullptr;
  else return res->second.getBXData(lumiSec);
}
  

void BXPUInfo::addBXData(int runnr,int lumiSec,BXData bxData)
{
  auto res = runData_.find(runnr);
  if(res==runData_.end()) res = runData_.insert({runnr,RunData(runnr)}).first;
  res->second.addBXData(lumiSec,std::move(bxData));
}

void BXPUInfo::addBXData(int runnr,int lumiSec,float liveTimeFrac,const std::vector<float>& bxDataVec)
{
  if(bxDataVec.size()!=kNrBX+1){
    LogErr<<" bxDataVec must have exactly "<<kNrBX+1<<" entries, has "<<bxDataVec.size()<<std::endl;
    return; 
  }
  BXData bxData(liveTimeFrac);
  for(size_t entryNr=0;entryNr<bxData.size();entryNr++){
    //    bxData[entryNr] = std::floor(bxDataVec[entryNr]+0.5);//rounding up when converting to int
    bxData.set(entryNr,bxDataVec[entryNr]);
  }
  addBXData(runnr,lumiSec,std::move(bxData));
}


const BXPUInfo::BXData* BXPUInfo::RunData::getBXData(int lumiSec)const
{
  auto res = lumiData_.find(lumiSec);
  if(res==lumiData_.end()) return nullptr;
  else return &(res->second);
}


void BXPUInfo::RunData::addBXData(int lumiSec,BXData&& bxData)
{
  auto res = lumiData_.emplace(lumiSec,bxData);
  if(!res.second){
    LogErr<<" warning run "<<runnr_<<" lumi sec "<<lumiSec<<" already existed"<<std::endl;
  }
}


// const BXPUInfo::RunData* BXPUInfo::getRunData(int runnr)
// {
//   auto res = std::equal_range(runData_.begin(),runData_.end(),RunData(runnr));
//   size_t nrFound = std::distance(res.first,res.second);
//   if(nrFound==0) return nullptr;
//   else if(nrFound>1){
//     LogErr<<" warning, "<<nrFound<<" entries found for run "<<runnr<<std::endl;
//   }
//   return res.first;
// }

// const BXData* BXPUInfo::getBXData(int runnr,int lumiSec)
// {
//   auto runData = getRunData(runnr);
//   if(runData){
//     auto lumiData = runData->getLumiData(runnr,lumiSec);
//     if(lumiData) return &lumiData->bxData();
//   }
//   return nullptr;
// }
  
// const BXPUInfo::LumiData* BXPUInfo::RunData::getLumiData(int runnr,int lumiSec)
// {
//   auto res = std::equal_range(lumiData_.begin(),lumiData_.end(),LumiData(runnr,lumiSec));
//   size_t nrFound = std::distance(res.first,res.second);
//   if(nrFound==0) return nullptr;
//   else if(nrFound>1){
//     LogErr<<" warning, "<<nrFound<<" entries found for lumi "<<lumiSec<<std::endl;
//   }
//   return res.first;
// }
 
