#ifndef UTILITY_BXPUINFO_HH
#define UTILITY_BXPUINFO_HH

//simple class to give per bx expected pileup  reading in from a previously generated file
//Author: Sam Harper (RAL), Aug 2018

#include "Utility/MiniFloatConverter.hh"
#include "TROOT.h"
#include <unordered_map>
#include <array>

class BXPUInfo {
public:
  static constexpr int kNrBX = 3564;
  static constexpr float kMinBiasXsec = 69200.;
  // using BXData = std::array<float,kNrBX+1>;
  class BXData {
  private:
    std::array<uint16_t,kNrBX+1> data_;
    float liveTimeFrac_;
  public:
    BXData(float liveTimeFrac=0.):liveTimeFrac_(liveTimeFrac){}
    BXData(std::array<float,kNrBX+1> data,float liveTimeFrac):
      liveTimeFrac_(liveTimeFrac){
      for(int bx=0;bx<=kNrBX;bx++) set(bx,data[bx]);
    }

    BXData(BXData&& rhs):data_(std::move(rhs.data_)),liveTimeFrac_(rhs.liveTimeFrac_){}
    BXData(const BXData& rhs)=default;
    BXData& operator=(const BXData& rhs)=default;
    BXData& operator=(BXData&& rhs){
      if(this!=&rhs){
	data_=std::move(rhs.data_);
	liveTimeFrac_=rhs.liveTimeFrac_;
      }
      return *this;
    }

    void setLiveTimeFrac(float val){liveTimeFrac_=val;}

    float get(int bx)const{return MiniFloatConverter::float16to32(data_[bx]);}
    void set(int bx,float val){data_[bx] = MiniFloatConverter::float32to16(val);}
    size_t size()const{return data_.size();}
    float liveTimeFrac()const{return liveTimeFrac_;}
  };

  class RunData{
  private:
    int runnr_;
    std::unordered_map<int,BXData> lumiData_;

  public:
    RunData(int runnr=0):runnr_(runnr){}
    const BXData* getBXData(int lumiSec)const;
    void addBXData(int lumiSec,BXData&& bxData);
  };

private: 
  std::unordered_map<int,RunData> runData_;
  
public:
  BXPUInfo(){}
  virtual ~BXPUInfo(){}
  float exptNrPUInt(int runnr,int lumiSec,int bx)const;
  std::pair<int,int> nrInTrain(int runnr,int lumiSec,int bx)const;
  int nrBXFromStartOfTrain(int runnr,int lumiSec,int bx)const;
  int nrBXFromEndOfTrain(int runnr,int lumiSec,int bx)const;
  const BXData* getBXData(int runnr,int lumiSec)const;
  void addBXData(int runnr,int lumiSec,BXData bxData);
  void addBXData(int runnr,int lumiSec,float liveTimeFrac,const std::vector<float>& bxData); //for easier python access

  //  bool isBXNthInTrain(int runnr,int lumiSec,int bx)const;

  ClassDef(BXPUInfo,3)

};

#endif
