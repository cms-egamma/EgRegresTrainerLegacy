#include "Utility/EcalCodes.hh"


ComCodes EcalCodes::codes_(EcalCodes::setCodes_());

//horribly inefficient I know but its done once
ComCodes EcalCodes::setCodes_()
{
  ComCodes codes;
  codes.setCode("kGood",kGood);
  codes.setCode("kPoorReco",kPoorReco);
  codes.setCode("kOutOfTime",kOutOfTime);
  codes.setCode("kFaultyHardware",kFaultyHardware);
  codes.setCode("kNoisy",kNoisy);
  codes.setCode("kPoorCalib",kPoorCalib);
  codes.setCode("kSaturated",kSaturated);
  codes.setCode("kLeadingEdgeRecovered",kLeadingEdgeRecovered);
  codes.setCode("kNeighboursRecovered",kNeighboursRecovered);
  codes.setCode("kTowerRecovered",kTowerRecovered);
  codes.setCode("kDead",kDead);
  codes.setCode("kKilled",kKilled);
  codes.setCode("kTPSaturated",kTPSaturated);
  codes.setCode("kL1SpikeFlag",kL1SpikeFlag);
  codes.setCode("kWeird",kWeird);
  codes.setCode("kDiWeird",kDiWeird);
  codes.setCode("kHasSwitchToGain6",kHasSwitchToGain6);
  codes.setCode("kHasSwitchToGain1",kHasSwitchToGain1);
  codes.setCode("kUnknown",kUnknown);
  
  codes.sort();
  return codes;
}

