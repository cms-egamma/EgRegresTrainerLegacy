#ifndef ECALCODES
#define ECALCODES


#include "Utility/ComCodes.hh"

#include <cstring>
#include <map>
#include <string>
#include <iostream>

class EcalCodes { //class to handle the ecalcodes 
public:
 
  enum EcalCode{
          kGood=0x1,                   // channel ok, the energy and time measurement are reliable
          kPoorReco=0x2,                 // the energy is available from the UncalibRecHit, but approximate (bad shape, large chi2)
          kOutOfTime=0x4,                // the energy is available from the UncalibRecHit (sync reco), but the event is out of time
          kFaultyHardware=0x8,           // The energy is available from the UncalibRecHit, channel is faulty at some hardware level (e.g. noisy)
          kNoisy=0x10,                    // the channel is very noisy
          kPoorCalib=0x20,                // the energy is available from the UncalibRecHit, but the calibration of the channel is poor
          kSaturated=0x40,                // saturated channel (recovery not tried)
          kLeadingEdgeRecovered=0x80,     // saturated channel: energy estimated from the leading edge before saturation
          kNeighboursRecovered=0x100,      // saturated/isolated dead: energy estimated from neighbours
          kTowerRecovered=0x200,           // channel in TT with no data link, info retrieved from Trigger Primitive
          kDead=0x400,                     // channel is dead and any recovery fails
          kKilled=0x800,                   // MC only flag: the channel is killed in the real detector
          kTPSaturated=0x1000,              // the channel is in a region with saturated TP
          kL1SpikeFlag=0x2000,              // the channel is in a region with TP with sFGVB = 0
          kWeird=0x4000,                    // the signal is believed to originate from an anomalous deposit (spike) 
          kDiWeird=0x8000,                  // the signal is anomalous, and neighbors another anomalous signal  
          kHasSwitchToGain6=0x10000,         // at least one data frame is in G6
          kHasSwitchToGain1=0x20000,         // at least one data frame is in G1
                                     //
          kUnknown=0x40000                   // to ease the interface with functions returning flags. 
  };
private:
  static ComCodes codes_;

private:
  EcalCodes(){} //not going to allow instainitiation
  ~EcalCodes(){}

public:
  static int getCode(const char *descript){return codes_.getCode(descript);}
  static void getCodeName(int code,std::string& id){return codes_.getCodeName(code,id);}
  static std::string getCodeName(int code){return codes_.getCodeName(code);}

private:
  static ComCodes setCodes_();
  
};

#endif
