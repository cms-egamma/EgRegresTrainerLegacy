#include "Utility/DetIdTools.hh"

#include "Utility/LogErr.hh"

#include<cstdlib>
#include<cmath>
#include<cstring>
#include<fstream>

#include <functional>
#include <boost/bind.hpp>
#include <algorithm>
#include <sstream>

std::vector<std::pair<int,int> > DetIdTools::eeDetIdToTowerId_;
const std::vector<std::vector<int> > DetIdTools::towerToRecHits_ =DetIdTools::makeTowerToRecHitsHashTable_();
const std::vector<int> DetIdTools::dummyTowerHits_;
DetIdTools::EcalNavigator::EcalNavigator(int startId):
  startId_(startId)
{ 
  if(DetIdTools::isEcalBarrel(startId)) isBarrel_=true;
  else if(DetIdTools::isEcalEndcap(startId)) isBarrel_=false;
  else{
    isBarrel_=false;
    LogErr <<" Error "<<startId<<" is not ecal barrel or endcap, navigation is undefined"<<std::endl;
    startId_=0;
  }
 
  initPosCoords_();
 
}

void DetIdTools::EcalNavigator::initPosCoords_()
{
  if(isBarrel_) {
    currEtaOrX_= iEtaBarrel(startId_);
    currPhiOrY_= iPhiBarrel(startId_);
    zSide_ = 0;
  }else{
    currEtaOrX_ = iXEndcap(startId_);
    currPhiOrY_ = iYEndcap(startId_);
    zSide_ = zEndcap(startId_); 
  }
}

int DetIdTools::EcalNavigator::goNorth(int nrSteps)
{
  currPhiOrY_+=nrSteps; 
  while(isBarrel_ && currPhiOrY_>360) currPhiOrY_-=360;
  while(isBarrel_ && currPhiOrY_<=0) currPhiOrY_+=360;
  return curPos();
}

int DetIdTools::EcalNavigator::goEast(int nrSteps)
{
  int oldEtaOrX_ = currEtaOrX_;
  currEtaOrX_+=nrSteps;
  //little fix in barrel as there is no crystal at eta =0;
  if(isBarrel_ && oldEtaOrX_*currEtaOrX_<=0){ //sign change
    if(nrSteps>0) currEtaOrX_++;
    else if(nrSteps<0) currEtaOrX_--;
  }
  
  return curPos();
}

int DetIdTools::EcalNavigator::curPos()const
{
  return getDetId_(currEtaOrX_,currPhiOrY_);
}

int DetIdTools::EcalNavigator::getIdAtPos(int nrStepsEast,int nrStepsNorth)const   
{
  int posEtaOrX = currEtaOrX_+nrStepsEast;
  if(isBarrel_ && posEtaOrX*currEtaOrX_<=0){ //little fix for barrel having no crystal at 0, going directly from -1 to 1
    if(nrStepsEast>0) posEtaOrX++;
    else if(nrStepsEast<0) posEtaOrX--;
  }
  int posPhiOrY = currPhiOrY_+nrStepsNorth;
  while(isBarrel_ && posPhiOrY>360) posPhiOrY-=360;
  while(isBarrel_ && posPhiOrY<=0) posPhiOrY+=360;
  return getDetId_(posEtaOrX,posPhiOrY);

}

int DetIdTools::EcalNavigator::getDetId_(int etaOrX,int phiOrY)const
{
  if(isBarrel_){
    if(isValidEcalBarrelId(etaOrX,phiOrY)){
      return makeEcalBarrelId(etaOrX,phiOrY);
    }else return 0;
  }else{
    // std::cout <<"endcap id for x "<<etaOrX<<" ("<<currEtaOrX_<<")"<<" "<<phiOrY<<" ("<<currPhiOrY_<<")";
    if(isValidEcalEndcapId(etaOrX,phiOrY,zSide_)){
      //std::cout <<" is valid, id is "<<makeEcalEndcapId(etaOrX,phiOrY,zSide_)<<std::endl;
      return makeEcalEndcapId(etaOrX,phiOrY,zSide_);
    }else{
      //std::cout <<"is not valid "<<std::endl;
      return 0;
    }
  }
}

DetIdTools::CaloNavigator::CaloNavigator(int startId,bool l1Navigator):
  startId_(startId),
  l1Navigator_(l1Navigator)
{
  if((!l1Navigator && isValidCaloId(startId)) || (l1Navigator_ && isValidL1CaloId(startId))){
    startId_=startId;
    initPosCoords_();
  }else {
    startId_=0;
    LogErr <<" Error "<<startId<<" is not valid calo id, navigation is undefined"<<std::endl;
  }
}

int DetIdTools::CaloNavigator::moveInEta(int nrSteps)
{
  int newEta = offsetInEta(currEta_,nrSteps);
  
  if(changedPhiSeg(newEta,currEta_)){
    if(abs(newEta)>kIEtaPhiSegChange){ //going from 72->36 segments, even numbers are now invalid, so sub one if thats the case
      if(currPhi_%2==0) currPhi_--;
    }//going the other way is fine
  }
  currEta_=newEta;
  
  return curPos();
}

int DetIdTools::CaloNavigator::offsetInEta(int startEta,int nrSteps)
{
  int newEta =startEta+=nrSteps;
  if(newEta*startEta<=0) { //sign change
    if(nrSteps>0) newEta++; 
    else if(nrSteps<0) newEta--;
  }
  return newEta;
}


int DetIdTools::CaloNavigator::moveInPhi(int nrSteps)
{
  //okay so we need to 
  if(abs(currEta_)>kIEtaPhiSegChange && !l1Navigator_) nrSteps*=2; //even segments are now invalid, we move 2 at time
  
  
  int newPhi = currPhi_+=nrSteps;
  while(newPhi<=0) newPhi+= kNrPhiSeg;
  while(newPhi>kNrPhiSeg) newPhi-=kNrPhiSeg;
  
  currPhi_=newPhi;

  return curPos();
}


int DetIdTools::CaloNavigator::offsetInL1Phi(int startPhi,int nrSteps)
{
  
  int newPhi = startPhi+=nrSteps;
  while(newPhi<=0) newPhi+= kNrPhiSeg;
  while(newPhi>kNrPhiSeg) newPhi-=kNrPhiSeg;
  return newPhi;
 
}


int DetIdTools::CaloNavigator::getIdAtPos(int nrStepsEta,int nrStepsPhi)const
{
  int newEta = currEta_+nrStepsEta;
  if(newEta*currEta_<=0) { //sign change
    if(nrStepsEta>0) newEta++;
    else if(nrStepsEta<0) newEta--; 
  }
  int newPhi = currPhi_;
  if(changedPhiSeg(newEta,currEta_)){
    if(abs(newEta)>kIEtaPhiSegChange){ //going from 72->36 segments, even numbers are now invalid, so sub one if thats the case
      if(newPhi%2==0) newPhi--;
    }//going the other way is fine
  }
  
  if(abs(newEta)>kIEtaPhiSegChange && !l1Navigator_) nrStepsPhi*=2; //even segments are now invalid, we move 2 at time  
  newPhi+=nrStepsPhi;
  while(newPhi<=0) newPhi+= kNrPhiSeg;
  while(newPhi>kNrPhiSeg) newPhi-=kNrPhiSeg;

  return l1Navigator_ ? DetIdTools::makeL1CaloDetId(newEta,newPhi) : DetIdTools::makeCaloDetId(newEta,newPhi);
  
}  
void DetIdTools::CaloNavigator::initPosCoords_()
{
  currEta_= iEtaCalo(startId_);
  currPhi_= iPhiCalo(startId_);
}

//magic numbers stolen from CMSSW EEDetId
const int DetIdTools::nBegin_[kIXMax_] = { 41, 41, 41, 36, 36, 26, 26, 26, 21, 21, 21, 21, 21, 16, 16, 14, 14, 14, 14, 14, 9, 9, 9, 9, 9, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 4, 4, 4, 4, 4, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 4, 4, 4, 4, 4, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 9, 9, 9, 9, 9, 14, 14, 14, 14, 14, 16, 16, 21, 21, 21, 21, 21, 26, 26, 26, 36, 36, 41, 41, 41 };
const int DetIdTools::nIntegral_[kIXMax_] = { 0, 20, 40, 60, 90, 120, 170, 220, 270, 330, 390, 450, 510, 570, 640, 710, 784, 858, 932, 1006, 1080, 1164, 1248, 1332, 1416, 1500, 1590, 1680, 1770, 1860, 1950, 2040, 2130, 2220, 2310, 2400, 2494, 2588, 2682, 2776, 2870, 2970, 3070, 3170, 3270, 3370, 3470, 3570, 3670, 3770, 3870, 3970, 4070, 4170, 4270, 4370, 4470, 4570, 4670, 4770, 4870, 4964, 5058, 5152, 5246, 5340, 5430, 5520, 5610, 5700, 5790, 5880, 5970, 6060, 6150, 6240, 6324, 6408, 6492, 6576, 6660, 6734, 6808, 6882, 6956, 7030, 7100, 7170, 7230, 7290, 7350, 7410, 7470, 7520, 7570, 7620, 7650, 7680, 7700, 7720 };

const std::vector<int> DetIdTools::ebFastHashTable_(DetIdTools::makeEBFastHashTable_());
const std::vector<int> DetIdTools::eeFastHashTable_(DetIdTools::makeEEFastHashTable_());
const std::unordered_map<int,int> DetIdTools::hcalFastHashTable_(DetIdTools::makeHcalFastHashTable_());



DetIdTools::DetIdTools()
{

}

int DetIdTools::makeEcalBarrelId(int iEta,int iPhi)
{
 
  int iEtaAbs = std::abs(iEta);
  int detId=(kEcalCode | kBarrelCode);
  detId |= (iEtaAbs<<9);
  detId |= iPhi;
  if(iEta>0) detId |= 0x10000;
  return detId;

}

std::ostream& DetIdTools::printEBDetId(int detId,std::ostream& out)
{
  if(!isEcalBarrel(detId)) out <<"not a EBDetId";
  else out <<" iEta "<<DetIdTools::iEtaBarrel(detId)<<" iPhi "<<DetIdTools::iPhiBarrel(detId);
  return out;
}

int DetIdTools::dIEtaBarrel(int iEta1,int iEta2)
{
  
  int dEta = iEta1-iEta2;
  if(iEta1*iEta2<0) {//-ve to +ve transistion and no crystal at zero
    if(dEta<0) dEta++;
    else dEta--;
  }
  return dEta;
}

int DetIdTools::dIPhiBarrel(int iPhi1,int iPhi2)
{

  int dPhi = iPhi1-iPhi2;

  if(dPhi>180) dPhi-=360;
  else if(dPhi<-180) dPhi+=360;
  
  return dPhi;

}

int DetIdTools::makeEcalEndcapId(int ix,int iy,int iz)
{
  int detId=(kEcalCode | kEndcapCode);
  detId |= (iy&0x7F);
  detId |= ((ix&0x7F)<<7);
  detId |= ((iz>0)?(0x4000):(0));
  return detId;
}

// int DetIdTools::getHashEcal(int detId)
// {
//   if(isEcal(detId)){
//     if(isBarrel(detId)) return getHashEcalBarrel(detId);
//     else if(isEndcap(detId)) return getHashEcalEndcap(detId)+kNrCrysBarrel;
//   }
//   //not ecal barrel or endcap
//   std::cout <<"DetIdTools::getHashEcal, detId "<<std::hex<<detId<<std::dec<<" is invalid"<<std::endl;
//   return -1;
// }


// int DetIdTools::getHashEcalBarrel(int detId)
// {
//   int etaBandNr =  kMaxIEtaBarrel + (positiveZBarrel(detId) ? iEtaAbsBarrel(detId)-1 : -iEtaAbsBarrel(detId)); // 0 - 189 starting at lowest eta (~-1.5) to highest
//   return etaBandNr* kMaxIPhiBarrel + iPhiBarrel(detId)-1;
// }


// int DetIdTools::getHashEcalEndcap(int detId)
// {
//   return iYEndcap(detId) - nBegin_[iXEndcap(detId)-1] + nIntegral_[iXEndcap(detId) -1 ] + (positiveZEndcap(detId) ? kICrFee_ : 0);
// }
int DetIdTools::makeHcalDetIdOld(int iEta,int iPhi,int depth)
{
  int subDetCode=0;
  if(abs(iEta)>=17 || (abs(iEta)==16 && depth==3) ) subDetCode = kEndcapCode;
  else subDetCode=kBarrelCode;
  return makeHcalDetIdOld(subDetCode,iEta,iPhi,depth);
}

int DetIdTools::makeHcalDetIdOld(int subDetCode,int iEta,int iPhi,int depth)
{
  int detId=0x0;
  detId |= kHcalCode;
  detId |= subDetCode;
  detId |= (iPhi&0x7F);
  if(iEta>0) detId |= (0x2000 | (iEta<<7) );
  else detId |= ((-iEta)<<7);
  detId |= ((depth&0x7)<<14);
  return detId;
}

int DetIdTools::makeHcalDetId(int iEta,int iPhi,int depth)
{
  int subDetCode=0;
  if(abs(iEta)>=17 || (abs(iEta)==16 && depth>=3) ) subDetCode = kEndcapCode;
  else subDetCode=kBarrelCode;
  return makeHcalDetId(subDetCode,iEta,iPhi,depth);
}

int DetIdTools::makeHcalDetId(int subDetCode,int iEta,int iPhi,int depth)
{
  int detId=0x0;
  if(isValidHcalId(subDetCode,iEta,iPhi,depth)){
    detId |= kHcalCode;
    detId |= subDetCode;
    detId |= (kHcalIdFormat2);
    detId |= ((depth&kHcalDepthMask2)<<kHcalDepthOffset2);
    if(iEta>0) detId |= kHcalZsideMask2 | (iEta<<kHcalEtaOffset2);
    else detId |= (-iEta<<kHcalEtaOffset2);
    detId |=  (iPhi&kHcalPhiMask2);
  }
  return detId;
}
void DetIdTools::printHcalDetId(int detId)
{
  std::cout <<"detId "<<std::hex<<detId<<std::dec<<" isHcal "<<isHcal(detId)<<" isBarrel "<<isBarrel(detId)<<" isEndcap "<<isEndcap(detId)<<" iEta "<<iEtaHcal(detId)<<" iPhi "<<iPhiHcal(detId)<<" depth "<<depthHcal(detId)<<std::endl;
}

int DetIdTools::makeCaloDetId(int iEta,int iPhi)
{
  if(!isValidCaloId(iEta,iPhi)) return 0;
     
  const int subDetCode=1<< kSubDetOffset;
  
  int detId=0x0;
  detId |= kCaloCode;
  detId |= subDetCode;
  detId |= (iPhi&0x7F);
  if(iEta>0) detId |= (0x2000 | ((iEta&0x3f)<<7) );
 else detId |= (((-iEta)&0x3f)<<7);
  
  return detId;
}

int DetIdTools::makeL1CaloDetId(int iEta,int iPhi)
{
  if(!isValidL1CaloId(iEta,iPhi)) return 0;
     
  const int subDetCode=1<< kSubDetOffset;
  
  int detId=0x0;
  detId |= kCaloCode;
  detId |= subDetCode;
  detId |= (iPhi&0x7F);
  if(iEta>0) detId |= (0x2000 | ((iEta&0x3f)<<7) );
 else detId |= (((-iEta)&0x3f)<<7);
  
  return detId;
}

// int DetIdTools::getHashHcal(int detId)
// {
//   //if(!DetIdTools::isHcal(detId)){
//     //std::cout <<"DetIdTools::getHashHcal: Warning det id "<<std::hex<<detId<<std::dec<<" is not in the hcal"<<std::endl;
//     // return -1;
//     //}

//   //removing the detector and sub detector info
//   int index = detId & ~(kDetMask | kSubDetMask);
//   //if(index<0 || index >= static_cast<int>(hcalFastHashTable_.size())){
//     //std::cout <<"DetIdTools::getHashHcal: Warning det id "<<std::hex<<detId<<std::dec<<" gives index "<<index<<" which is not valid (max="<<hcalFastHashTable_.size()<<")"<<std::endl;
//     // }
//   return hcalFastHashTable_[index];
// }


// int DetIdTools::getHashEcal(int detId)
// {
//   if(DetIdTools::isEcalBarrel(detId)) return getHashEcalBarrel(detId);
//   else if(DetIdTools::isEcalEndcap(detId)) return getHashEcalEndcap(detId)+  kNrEcalCellsBarrel;
//   else{
//     std::cout <<"DetIdTools::getHashEcal: Warning det id "<<std::hex<<detId<<std::dec<<" is not ecal barrel or endcap"<<std::endl;
//     return -1;
//   }
// }

 
// int DetIdTools::getHashEcalBarrel(int detId)
// {
//   //if(!DetIdTools::isEcalBarrel(detId)){
//     //std::cout <<"DetIdTools::getHashEcalBarrel: Warning det id "<<std::hex<<detId<<std::dec<<" is not ecal barrel"<<std::endl;
//     //return -1;
//     //}
//   //removing the detector and sub detector info
//   int index = detId & ~(kDetMask | kSubDetMask);
//   //if(index<0 || index >= static_cast<int>(ebFastHashTable_.size())){
//   // std::cout <<"DetIdTools::getHashEcalBarrel: Warning det id "<<std::hex<<detId<<std::dec<<" gives index "<<index<<" which is not valid (max="<<ebFastHashTable_.size()<<")"<<std::endl;
//   //}

//   return ebFastHashTable_[index];
// }

// int DetIdTools::getHashEcalEndcap(int detId)
// {
//   //if(!DetIdTools::isEcalEndcap(detId)){
//   // std::cout <<"DetIdTools::getHashEcalEndcap: Warning det id "<<std::hex<<detId<<std::dec<<" is not ecal endcap"<<std::endl;
//   // return -1;
//   //}
//   //removing the detector and sub detector info
//   int index = detId & ~(kDetMask | kSubDetMask);
//   //if(index<0 || index >= static_cast<int>(eeFastHashTable_.size())){
//   // std::cout <<"DetIdTools::getHashEcalEndcap: Warning det id "<<std::hex<<detId<<std::dec<<" gives index "<<index<<" which is not valid (max="<<eeFastHashTable_.size()<<")"<<std::endl;
//   //}

//   return eeFastHashTable_[index];
// }

//ported from HcalDetId.hash_index in CMSSW
int DetIdTools::calHashHcalLegacy(int detId)
{
  if(!DetIdTools::isHcal(detId)){
    std::cout <<"DetIdTools::calHashHcalLegacy: Warning det id "<<std::hex<<detId<<std::dec<<" is not in the hcal"<<std::endl;
    return -1;
  }
  
  int index = -1;

  int HBhalf = kNrHcalCellsBarrel/2;
  int HEhalf = kNrHcalCellsEndcap/2;
  
  int iPhi = DetIdTools::iPhiHcal(detId);
  int iEtaAbs = DetIdTools::iEtaAbsHcal(detId);
  int zSide = DetIdTools::zSideHcal(detId);
  int depth = DetIdTools::depthHcal(detId);

  // HB valid DetIds: phi=1-72,eta=1-14,depth=1; phi=1-72,eta=15-16,depth=1-2
  if (DetIdTools::isHcalBarrel(detId)) {
    if (iEtaAbs < 16)   index = (iPhi - 1)*18 + (iEtaAbs - 1) + (depth - 1);
    if (iEtaAbs == 16)  index = (iPhi - 1)*18 + iEtaAbs + (depth - 1);
    
    if (zSide == -1) index += HBhalf;
    // HE valid DetIds: phi=1-72,eta=16-17,depth=1; phi=1-72,eta=18-20,depth=1-2; 
  //                  phi=1-71(in steps of 2),eta=21-26,depth=1-2; phi=1-71(in steps of 2),eta=27-28,depth=1-3
  //                  phi=1-71(in steps of 2),eta=29,depth=1-2
    
  } else if (DetIdTools::isHcalEndcap(detId)) {
    if (iEtaAbs == 16 || iEtaAbs == 17)  index = (iPhi - 1)*8 + (iPhi/2)*20 + (iEtaAbs - 16);
    if (iEtaAbs >= 18 && iEtaAbs <= 20)  index = (iPhi - 1)*8 + (iPhi/2)*20 + 2  + 2*(iEtaAbs-18) + (depth - 1);
    if (iEtaAbs >= 21 && iEtaAbs <= 26)  index = (iPhi - 1)*8 + (iPhi/2)*20 + 8  + 2*(iEtaAbs-21) + (depth - 1);
    if (iEtaAbs >= 27 && iEtaAbs <= 28)  index = (iPhi - 1)*8 + (iPhi/2)*20 + 20 + 3*(iEtaAbs-27) + (depth - 1);
    if (iEtaAbs == 29)                     index = (iPhi - 1)*8 + (iPhi/2)*20 + 26 + 2*(iEtaAbs-29) + (depth - 1);
    
    index += 2*HBhalf;
    if (zSide == -1) index += HEhalf;
  }
  
  return index;
}

int DetIdTools::calHashHcal(int detId)
{
  if(!DetIdTools::isHcal(detId)){
    std::cout <<"DetIdTools::calHashHcal: Warning det id "<<std::hex<<detId<<std::dec<<" is not in the hcal"<<std::endl;
    return -1;
  }
  
  int index = -1;

  int HBhalf = kNrHcalCellsBarrel/2;
  int HEhalf = kNrHcalCellsEndcap/2;
  
  int iPhi = DetIdTools::iPhiHcal(detId);
  int iEtaAbs = DetIdTools::iEtaAbsHcal(detId);
  int zSide = DetIdTools::zSideHcal(detId);
  int depth = DetIdTools::depthHcal(detId);

  // HB valid DetIds: phi=1-72,eta=1-14,depth=1; phi=1-72,eta=15-16,depth=1-2
  if (DetIdTools::isHcalBarrel(detId)) {
    if (iEtaAbs < 16)   index = (iPhi - 1)*18 + (iEtaAbs - 1) + (depth - 1);
    if (iEtaAbs == 16)  index = (iPhi - 1)*18 + iEtaAbs + (depth - 1);
    
    if (zSide == -1) index += HBhalf;
  } else if (DetIdTools::isHcalEndcap(detId)) {
    //we're going to change it to be based on eta for the endcap as its easier (as #phi segements changes in eta)
    //we're also going to be compatible with phase 0 and phase 1
    //therefore eta=16 depth=4 (PhaseI) and eta=16, depth=3 (Phase0) give the same hash
    //likewise eta=17 depth=2 (PhaseI) and eta=17, depth=1 (Phase0) give the same hash
    if (iEtaAbs == 16) index = (iPhi -1);
    if (iEtaAbs == 17){
      int effectiveDepth = depth;
      if(depth>1) effectiveDepth--;
      index = kHcalIPhiMax+(iPhi-1)*2+(effectiveDepth-1);
    }
    if (iEtaAbs == 18)  index = kHcalIPhiMax*3 + (iPhi-1)*5 + (depth-1);
    if (iEtaAbs >= 19 && iEtaAbs <= 20)  index = kHcalIPhiMax*(8 + 6*(iEtaAbs-19)) + (iPhi-1)*6 + (depth-1);
    if (iEtaAbs >= 21 && iEtaAbs <= 25)  index = kHcalIPhiMax*20 + kHcalIPhiMax/2*6*(iEtaAbs-21) + (iPhi/2)*6 + (depth-1);
    if (iEtaAbs >= 26 && iEtaAbs <= 28)  index = kHcalIPhiMax*20 + kHcalIPhiMax/2*6*5 + kHcalIPhiMax/2*7*(iEtaAbs-26) + (iPhi/2)*7 + (depth-1);
    if (iEtaAbs == 29)                   index = kHcalIPhiMax*35 + kHcalIPhiMax/2*21 + (iPhi/2)*3 + (depth-1);
    
    index += 2*HBhalf;
    if (zSide == -1) index += HEhalf;
  }
  
  return index;
}

int DetIdTools::calHashCalo(int detId)
{
  if(!isValidCaloId(detId)) return -1;
  int iPhi = DetIdTools::iPhiCalo(detId);
  int iEtaAbs = DetIdTools::iEtaAbsCalo(detId);
  int zSide = DetIdTools::zSideCalo(detId);
  
  //std::cout <<"iEta "<<iEtaAbs<<" iPhi "<<iPhi<<std::endl;

  int index=0;
  if(iEtaAbs<=20) index=(iEtaAbs-1)*72+iPhi-1;
  else index = 20*72+(iEtaAbs-20)*36+(iPhi+1)/2-1; //only odd iPhi = 1, 3, 5 are now valid, add 1 and divide by 2 to get the compact phi 
  
  if(zSide>0) index+=kNrCaloTowers/2;
  return index;
}

//you know, using the info in SHCaloGeom, I could achieve this faster and in half the lines in one loop
//but it somehow doesnt feel right going there
//if any entry has zero, its a wildcard, to select all cells in tower 20, on both sides at depth 2, do
//getMatchingIdsHcal(20,0,0,2,ids)
void DetIdTools::getMatchingIdsHcal(int etaAbs,int phi,int side,int depth,std::vector<int>& ids)
{
  for(int etaAbsNr=kHcalIEtaAbsMin;etaAbsNr<=kHcalIEtaAbsMax;etaAbsNr++){
    for(int phiNr=kHcalIPhiMin;phiNr<kHcalIPhiMax;phiNr++){
      for(int depthNr=kHcalDepthMin;depthNr<=kHcalDepthMax;depthNr++){
	for(int sideNr=-1;sideNr<=1;sideNr+=2){
	  if(isValidHcalId(etaAbsNr*sideNr,phiNr,depthNr)){
	    if( (etaAbsNr==etaAbs || etaAbs==0) && //eta match
		(phiNr==phi || phi==0) && //phi match
		(depthNr==depth || depth==0) && //depth match
		(sideNr==side || side==0) ) { //side match
	      int subDet=0;
	      if(etaAbsNr<=15 || (etaAbsNr==16 && depthNr<3)) subDet=kBarrelCode;
	      else subDet=kEndcapCode;
	      ids.push_back(makeHcalDetId(subDet,etaAbsNr*sideNr,phiNr,depthNr));
	    }//end of cell match check
	  }//end of valid hcal det id check
	}//end of side loop
      }//end of depth loop
    }//end of phi loop
  }//end of eta loop
}

int DetIdTools::calHashL1Calo(int iEta,int iPhi)
{
  if(abs(iEta)<kL1CaloIEtaAbsMin || abs(iEta)>kL1CaloIEtaAbsMax || 
     iPhi<kL1CaloIPhiMin || iPhi>kL1CaloIPhiMax) return -1; //invalid id

  int etaIndex=iEta+kL1CaloIEtaAbsMax;
  if(iEta>0) etaIndex--; // theres is no 0 in ieta
  int phiIndex = iPhi-1; //going from 1-72 to 0-71
  return etaIndex*kL1CaloIPhiMax+phiIndex;

}

bool DetIdTools::isNextToBarrelPhiGap(int detId)
{
  if(isEcalBarrel(detId)) { //ecal
    if((iPhiBarrel(detId)%20)<=1) return true;
    else return false;
  }else if(isHcal(detId) && iEtaAbsHcal(detId)<=17){ //17 is in the endcap of hcal but covers barrel region of ecal
    if((iPhiHcal(detId)%4)<=1) return true;
    else return false;
  }
  return false;
}


bool DetIdTools::isNextToBarrelEtaGap(int detId,int maxDistToGap)
{
  if(isEcalBarrel(detId)) { //ecal
    int iEtaAbs = iEtaAbsBarrel(detId);
    const int nrCellsNextToGap = 8;
    const int cellsNextToGap[nrCellsNextToGap]={1,25,26,45,46,65,66,85};
    for(int cellNr=0;cellNr<nrCellsNextToGap;cellNr++){
      if(abs(iEtaAbs-cellsNextToGap[cellNr])<=maxDistToGap) return true;
    }
  }else if(isHcal(detId)) { //hcal
    int iEtaAbs = iEtaAbsHcal(detId);
    const int nrCellsNextToGap = 8;
    const int cellsNextToGap[nrCellsNextToGap]={1,5,6,9,10,13,14,17};
    for(int cellNr=0;cellNr<nrCellsNextToGap;cellNr++){
      if(abs(iEtaAbs-cellsNextToGap[cellNr])<=maxDistToGap) return true;
    }
  }
  return false;
}



int DetIdTools::nrOfNearestGap(int detId)
{
  const int errorCode=-10;
  if(isEcalBarrel(detId)) { //ecal
    int iEtaAbs = iEtaAbsBarrel(detId);
    int gapNr=0;
    if(iEtaAbs>=1 && iEtaAbs<=13) gapNr=0; //13 = special case as its the only one equidistant from two cracks, will go for the centre crack (as its so far away from gap it doesnt really matter)
    else if(iEtaAbs>=14 && iEtaAbs<=35) gapNr=1;
    else if(iEtaAbs>=36 && iEtaAbs<=55) gapNr=2;
    else if(iEtaAbs>=56 && iEtaAbs<=75) gapNr=3;
    else if(iEtaAbs>=76 && iEtaAbs<=85) gapNr=4;
    else return errorCode;
   

    return iEtaBarrel(detId) > 0 ? gapNr : gapNr*-1;
  }else if(isHcal(detId)){
    int iEtaAbs = iEtaAbsHcal(detId);
    int gapNr=0;
    if(iEtaAbs>=1 && iEtaAbs<=3) gapNr=0; //3 special case as equidistance, go for centre gap
    else if(iEtaAbs>=4 && iEtaAbs<=7) gapNr=1;
    else if(iEtaAbs>=8 && iEtaAbs<=11) gapNr=2;
    else if(iEtaAbs>=12 && iEtaAbs<=15) gapNr=3;
    else if(iEtaAbs>=16 && iEtaAbs<=17) gapNr=4;
    else return errorCode;
    
    return iEtaHcal(detId) > 0 ? gapNr : gapNr*-1;
  }
  return errorCode;
}


bool DetIdTools::isValidPhase0HcalId(int iEta,int iPhi,int depth)
{  
  return isValidPhase0HcalBarrelId(iEta,iPhi,depth) || isValidPhase0HcalEndcapId(iEta,iPhi,depth);
}

bool DetIdTools::isValidPhase1HcalId(int iEta,int iPhi,int depth)
{
  return isValidPhase1HcalBarrelId(iEta,iPhi,depth) || isValidPhase1HcalEndcapId(iEta,iPhi,depth);
}

bool DetIdTools::isValidPhase1HcalId(int subdetId,int iEta,int iPhi,int depth)
{
  if(subdetId==kBarrelCode) return isValidPhase1HcalBarrelId(iEta,iPhi,depth);
  else if(subdetId==kEndcapCode) return isValidPhase1HcalEndcapId(iEta,iPhi,depth);
  else return false;
}


bool DetIdTools::isValidCaloId(int iEta,int iPhi)
{
  if(abs(iEta)>=1 && abs(iEta)<=29 && iPhi>=1 && iPhi<=72){
    if(abs(iEta)>20){
      if(iPhi%2==1) return true;
      else return false;
    }else return true;
  }
  return false;
}

bool DetIdTools::isValidL1CaloId(int iEta,int iPhi)
{
  if(abs(iEta)>=1 && abs(iEta)<=28 && iPhi>=1 && iPhi<=72) return true;
  else return false;
}


//depth 1= all towers up to and including tower 17 + depth 1 18-29 + depth 2 27-29
//depth 2= depth 2 18-26, depth 3 27-29
//warning, emulating calo tower bug where depth 2 tower 27 is depth 2 and depth 3 is depth 1
int DetIdTools::getEffectiveHcalDepth(int detId)
{
  //if(DetIdTools::newFormatHcal(detId)){
    //LogErr<<" new format HCAL detId, new function not updated to this"<<std::endl;
    //}
  if(!isHcal(detId)){
    LogErr <<" : Warning detId "<<detId<<" is not hcal "<<std::endl;
    return 0;
  }else{
    int iEtaAbs = iEtaAbsHcal(detId);
    int depth = depthHcal(detId);
    if(iEtaAbs<=17 || 
       (iEtaAbs<=29 && depth==1) ||
       (iEtaAbs>=27 && iEtaAbs<=29 && depth==2)){   
      //(iEtaAbs>=28 && iEtaAbs<=29 && depth==2) || 
      // (iEtaAbs==27 && depth==3)){
      return 1;
    }else return 2;
  }
}  

//towers with one depth 1-15,17
//towers with two depths 18-26
//towers with three depths 16,27-29 (
int DetIdTools::getNrDepthsInHcalTower(int detId)
{

  if(!isHcal(detId)){
    LogErr <<" : Warning detId "<<detId<<" is not hcal,not implimented yet for ecal ids "<<std::endl;
    return 0;
  }else if(newFormatHcal(detId)){
    LogErr <<" : Warning detId "<<detId<<" is new format, not validated yet for this "<<std::endl;
    return 0;
  }else{
    int iEtaAbs = iEtaAbsHcal(detId);
    if(iEtaAbs<=15) return 1;
    else if(iEtaAbs==16) return 3;
    else if(iEtaAbs==17) return 1;
    else if(iEtaAbs<=26) return 2;
    else if(iEtaAbs<=29) return 3;
    else {
      LogErr <<" : Warning detId "<<detId<<" has invalid abs eta "<<iEtaAbs<<std::endl;
      return -1;
    }
  }//end valid hcal det id check
}

bool DetIdTools::isValidPhase0HcalBarrelId(int iEta,int iPhi,int depth)
{
  const int hcalBarrelIEtaMax = 16;
  const int hcalBarrelIPhiMin = 1;
  const int hcalBarrelIPhiMax = 72;
  
  if(abs(iEta)>=1 && abs(iEta)<=hcalBarrelIEtaMax){ //iEta good
    if(iPhi>=hcalBarrelIPhiMin && iPhi<=hcalBarrelIPhiMax){ //iPhi good
      if(depth==1 || (abs(iEta)>=15 && depth==2)){ //depth good
	return true;
      }
    }
  }
  return false;
}

bool DetIdTools::isValidPhase0HcalEndcapId(int iEta,int iPhi,int depth)
{
  const int iEtaMin = 16;
  const int iEtaMax = 29;
  const int iEtaPhiBoundary=21;
  const int iPhiMin = 1;
  const int iPhiMax = 72;// 


  //first check ieta and phi 
  if(abs(iEta)<iEtaMin || abs(iEta)>iEtaMax) return false;
  else if(iPhi<iPhiMin || iPhi>iPhiMax) return false;

  //final phi check, if its above the eta boundary, it goes like 1,3,5,..,69,71
  //so we check it isnt odd
  if(abs(iEta)>=iEtaPhiBoundary && iPhi%2==0) return false;

  //depth checks
  if(abs(iEta)==16){
    if(depth==3) return true;
    else return false;
  }else if(abs(iEta)==17){
    if(depth==1) return true;
    else return false;
  } else if(abs(iEta)<=26){
    if(depth>=1 && depth<=2) return true;
    else return false;
  }else if(abs(iEta)<=28){
    if(depth>=1 && depth<=3) return true;
    else return false;
  }else if(abs(iEta)<=29){
    if(depth>=1 && depth<=2) return true;
    else return false;
  }else return false;
}

bool DetIdTools::isValidPhase1HcalEndcapId(int iEta,int iPhi,int depth)
{
  const int iEtaMin = 16;
  const int iEtaMax = 29;
  const int iEtaPhiBoundary=21;
  const int iPhiMin = 1;
  const int iPhiMax = 72;// 


  //first check ieta and phi 
  if(abs(iEta)<iEtaMin || abs(iEta)>iEtaMax) return false;
  else if(iPhi<iPhiMin || iPhi>iPhiMax) return false;

  //final phi check, if its above the eta boundary, it goes like 1,3,5,..,69,71
  //so we check it isnt odd
  if(abs(iEta)>=iEtaPhiBoundary && iPhi%2==0) return false;

  //depth checks
  if(abs(iEta)==16){
    if(depth==4) return true;
    else return false;
  }else if(abs(iEta)==17){
    if(depth==2 || depth==3) return true;
    else return false;
  }else if(abs(iEta)==18){
    if(depth>=1 && depth<=5) return true;
    else return false;
  } else if(abs(iEta)<=25){
    if(depth>=1 && depth<=6) return true;
    else return false;
  }else if(abs(iEta)<=28){
    if(depth>=1 && depth<=7) return true;
    else return false;
  }else if(abs(iEta)<=29){
    if(depth>=1 && depth<=3) return true;
    else return false;
  }else return false;
}

bool DetIdTools::isValidEcalId(int detId)
{
  if(DetIdTools::isEcalBarrel(detId)) return isValidEcalBarrelId(detId);
  else if(DetIdTools::isEcalEndcap(detId)) return isValidEcalEndcapId(detId);
  else return false;
}

bool DetIdTools::isValidEcalBarrelId(int iEta, int iPhi) 
{
  bool valid = true;
  if (iEta < -kMaxIEtaBarrel || iEta == 0 || iEta > kMaxIEtaBarrel ||
      iPhi < kMinIPhiBarrel || iPhi > kMaxIPhiBarrel) {
    valid = false;
  }  
  return valid;

}
//stolen from EEDetId.cc in CMSSW
bool DetIdTools::isValidEcalEndcapId(int crystal_ix,int crystal_iy,int iz)
{
  bool valid = false;
  if (crystal_ix < kMinIXEndcap || crystal_ix > kMaxIXEndcap ||
      crystal_iy < kMinIYEndcap || crystal_iy > kMaxIYEndcap || abs(iz) != 1 ) 
    { 
      return valid ; 
    }
  if ( (crystal_ix >= 1 && crystal_ix <= 3 && (crystal_iy <= 40 || crystal_iy > 60) ) ||
       (crystal_ix >= 4 && crystal_ix <= 5 && (crystal_iy <= 35 || crystal_iy > 65) ) || 
       (crystal_ix >= 6 && crystal_ix <= 8 && (crystal_iy <= 25 || crystal_iy > 75) ) || 
       (crystal_ix >= 9 && crystal_ix <= 13 && (crystal_iy <= 20 || crystal_iy > 80) ) || 
       (crystal_ix >= 14 && crystal_ix <= 15 && (crystal_iy <= 15 || crystal_iy > 85) ) || 
       (crystal_ix >= 16 && crystal_ix <= 20 && (crystal_iy <= 13 || crystal_iy > 87) ) || 
       (crystal_ix >= 21 && crystal_ix <= 25 && (crystal_iy <= 8 || crystal_iy > 92) ) || 
       (crystal_ix >= 26 && crystal_ix <= 35 && (crystal_iy <= 5 || crystal_iy > 95) ) || 
       (crystal_ix >= 36 && crystal_ix <= 39 && (crystal_iy <= 3 || crystal_iy > 97) ) || 
       (crystal_ix >= 98 && crystal_ix <= 100 && (crystal_iy <= 40 || crystal_iy > 60) ) ||
       (crystal_ix >= 96 && crystal_ix <= 97 && (crystal_iy <= 35 || crystal_iy > 65) ) || 
       (crystal_ix >= 93 && crystal_ix <= 95 && (crystal_iy <= 25 || crystal_iy > 75) ) || 
       (crystal_ix >= 88 && crystal_ix <= 92 && (crystal_iy <= 20 || crystal_iy > 80) ) || 
       (crystal_ix >= 86 && crystal_ix <= 87 && (crystal_iy <= 15 || crystal_iy > 85) ) || 
       (crystal_ix >= 81 && crystal_ix <= 85 && (crystal_iy <= 13 || crystal_iy > 87) ) || 
       (crystal_ix >= 76 && crystal_ix <= 80 && (crystal_iy <= 8 || crystal_iy > 92) ) || 
       (crystal_ix >= 66 && crystal_ix <= 75 && (crystal_iy <= 5 || crystal_iy > 95) ) || 
       (crystal_ix >= 62 && crystal_ix <= 65 && (crystal_iy <= 3 || crystal_iy > 97) ) ||
       ( (crystal_ix == 40 || crystal_ix == 61) && ( (crystal_iy >= 46 && crystal_iy <= 55 ) || crystal_iy <= 3 || crystal_iy > 97 )) ||
       ( (crystal_ix == 41 || crystal_ix == 60) && crystal_iy >= 44 && crystal_iy <= 57 ) ||
       ( (crystal_ix == 42 || crystal_ix == 59) && crystal_iy >= 43 && crystal_iy <= 58 ) ||
       ( (crystal_ix == 43 || crystal_ix == 58) && crystal_iy >= 42 && crystal_iy <= 59 ) ||
       ( (crystal_ix == 44 || crystal_ix == 45 || crystal_ix == 57 || crystal_ix == 56) && crystal_iy >= 41 && crystal_iy <= 60 ) ||
       ( crystal_ix >= 46 && crystal_ix <= 55 && crystal_iy >= 40 && crystal_iy <= 61 ) 
       )
    { 
      return valid; 
    }
  valid = true;
  return valid;
  
}


int DetIdTools::endcapEtaRing(int detId)
{
  if(DetIdTools::isEcalEndcap(detId)){
    return endcapEtaRing(DetIdTools::iXEndcap(detId),DetIdTools::iYEndcap(detId));
  }else return -1;
}


//there are 40 rings in eta (i think)
int DetIdTools::endcapEtaRing(int ix,int iy)
{
  if(DetIdTools::isValidEcalEndcapId(ix,iy,1)){
    int iXNorm  = DetIdTools::normEndcapIXOrIY(ix);
    int iYNorm  = DetIdTools::normEndcapIXOrIY(iy);
    float iR = std::sqrt(iXNorm*iXNorm+iYNorm*iYNorm);

    int iRInt = static_cast<int>(iR);
    if(iR-iRInt>0.5) iRInt++; //rounding up

    return iRInt-10;
  }else return -1;
}

float DetIdTools::endcapEtaRingFloat(int ix,int iy)
{
  if(DetIdTools::isValidEcalEndcapId(ix,iy,1)){
    int iXNorm  = DetIdTools::normEndcapIXOrIY(ix);
    int iYNorm  = DetIdTools::normEndcapIXOrIY(iy);
    float iR = std::sqrt(iXNorm*iXNorm+iYNorm*iYNorm);
    
    return iR;//-10.;
  }else return -1;
}

float DetIdTools::endcapEtaRingFloat(int detId)
{
  if(DetIdTools::isEcalEndcap(detId)){
    return endcapEtaRingFloat(DetIdTools::iXEndcap(detId),DetIdTools::iYEndcap(detId));
  }else return -1.;
}

bool DetIdTools::isNextToRingBoundary(int detId)
{
  if(!isValidEcalEndcapId(detId)) return false;
     
  int ix=iXEndcap(detId);
  int iy=iYEndcap(detId);
  int zSide = zEndcap(detId);

  for (int i = -1; i <= 1; ++i) {
    for (int j = -1; j <= 1; ++j) {
      if ( ! isValidEcalEndcapId( ix + i, iy + j, zSide ) ) {
	return true;
      }
    }
  }
  return false;
}

int DetIdTools::normEndcapIXOrIY(int iXOrIY)
{
  int iXOrIYNorm = iXOrIY -50;
  //want to map 1=-50,50=-1,51=1 and 100 to 50 so sub off one if zero or neg
  if(iXOrIYNorm<=0) iXOrIYNorm--;
  return iXOrIYNorm;

}

const std::vector<int>& DetIdTools::getRecHitsOfTower(int towerId)
{
  if(isValidCaloId(towerId)){
    size_t hash = calHashCalo(towerId);
    if(hash<towerToRecHits_.size()) return towerToRecHits_[hash];
    else{
      LogErr <<" warning hash of "<<towerId<<" is "<<hash<<" which is bigger than "<<towerToRecHits_.size()<<std::endl;
    }
  }
  return dummyTowerHits_;
}


std::vector<int> DetIdTools::makeEBFastHashTable_()
{
  //not all entries will be filled, the 0x20000 is because the z, eta,phi max bit is 0x10000 so we are guaranteed to have enough space
  //note that later it might be better to work out the max value for more eff memory use
  const int nrEntries = 0x20000;
  std::vector<int> hashTable(nrEntries,-1);
  for(int index=0;index<nrEntries;index++){
    const int detId = index+ kEcalCode + kBarrelCode;
    if(isValidEcalBarrelId(detId)){ //its valid get the hash
      hashTable[index]=calHashEcalBarrel(detId);
    }
  }
  return hashTable;
}

std::vector<int> DetIdTools::makeEEFastHashTable_()
{
  //not all entries will be filled, the 0x8000 is because the z, eta,phi max bit is 0x4000 so we are guaranteed to have enough space
  //note that later it might be better to work out the max value for more eff memory use
  const int nrEntries = 0x8000;
  std::vector<int> hashTable(nrEntries,-1);
  for(int index=0;index<nrEntries;index++){
    const int detId = index+ kEcalCode + kEndcapCode;
    if(isValidEcalEndcapId(detId)){ //its valid get the hash
      hashTable[index]=calHashEcalEndcap(detId);
    }
  }
  return hashTable;
}

struct PairLessThan {

  bool operator()(int lhs,int rhs)const{return lhs<rhs;}
  bool operator()(int lhs,const std::pair<int,int>& rhs)const{return lhs<rhs.first;}
  bool operator()(const std::pair<int,int>&lhs,int rhs)const{return lhs.first<rhs;}
  bool operator()(const std::pair<int,int>&lhs,const std::pair<int,int>& rhs)const{return lhs.first<rhs.first;}
};

int DetIdTools::towerIdEndcap(int detId)
{
  if(eeDetIdToTowerId_.empty()) fillEEToTowerIdMap("input/CaloTowerEEGeometric.map");
  if(eeDetIdToTowerId_.empty() && std::getenv("CMSSW_BASE")){ //try for CMSSW setup
    fillEEToTowerIdMap(std::string(std::getenv("CMSSW_BASE"))+"/src/SHarper/SHNtupliser/data/CaloTowerEEGeometric.map");
  }
  if(eeDetIdToTowerId_.empty()){
    //    LogErr <<" Warning, you didnt load in a tower map, I didnt find the default, so I'm now going to put a dummy entry in which will stop all warnings and just fail silently "; //commented out to not give a confusing warning to general users
    eeDetIdToTowerId_.push_back(std::make_pair<int,int>(0,0)); 
  }
  typedef std::vector<std::pair<int,int> >::const_iterator ConstIt;
  std::pair<ConstIt,ConstIt> result = std::equal_range(eeDetIdToTowerId_.begin(),eeDetIdToTowerId_.end(),detId,PairLessThan());
						    
  if(result.second!=result.first){
    if(result.second-result.first>1) LogErr<<" Warrning multiple entries for det id "<<detId<<std::endl;
    return result.first->second;
  }
  return 0;
}
int DetIdTools::ecalToTowerId(int detId)
{
  if(isEcalBarrel(detId)) return makeCaloDetId(iEtaBarrelTower(detId),iPhiBarrelTower(detId));
  else if(isEcalEndcap(detId)) return towerIdEndcap(detId);
  else return 0;
}


int DetIdTools::newToOldFormatHcal_(int detId)
{
  return makeHcalDetIdOld(DetIdTools::iEtaHcal(detId),
			  DetIdTools::iPhiHcal(detId),
			  DetIdTools::depthHcal(detId));
}


std::unordered_map<int,int> DetIdTools::makeHcalFastHashTable_()
{
  std::unordered_map<int,int> hashTable;

  for(int side=1;side>=0;side--){   
    for(int iEtaAbs=kHcalIEtaAbsMin;iEtaAbs<=kHcalIEtaAbsMax;iEtaAbs++){
      for(int iPhi=kHcalIPhiMin;iPhi<=kHcalIPhiMax;iPhi++){
	for(int depth=kHcalDepthMin;depth<=kHcalDepthMax;depth++){
	  int iEta = iEtaAbs*(2*side-1);  
	  int detId =0;
	  if(isValidHcalBarrelId(iEta,iPhi,depth)){
	    detId = DetIdTools::makeHcalBarrelDetId(iEta,iPhi,depth);
	  }else if(isValidHcalEndcapId(iEta,iPhi,depth)){
	    detId = DetIdTools::makeHcalEndcapDetId(iEta,iPhi,depth);
	  }
	  if(detId){
	    auto res = hashTable.insert({detId,calHashHcal(detId)});
	    if(!res.second) LogErr<<" HCAL "<<iEta<<" "<<iPhi<<" "<<depth<<" already added "<<std::endl;
	  }
	}
      }
    }
  }
   
  return hashTable;
}
 

void DetIdTools::fillEEToTowerIdMap(const std::string& filename)
{
  
  eeDetIdToTowerId_.clear();
  if(!filename.empty()){
    std::ifstream file(filename.c_str());
    if(file.bad()){
      std::cout <<"file map not found"<<std::endl;
    } else if(file.eof()){
      std::cout <<"file map found but empty"<<std::endl;
    }
    for(std::string line;std::getline(file,line); ){
      std::stringstream ss(line);
      int detId,eta,phi;
      if(ss >> std::hex >> detId >> std::dec >> eta >> phi){
	int towerId = makeCaloDetId(eta,phi);
	 //	std::cout <<" detId "<<detId<<" "<<eta<<" "<<phi<<std::endl;
	eeDetIdToTowerId_.push_back(std::make_pair(detId,towerId));
      }   
    }
  }//empty file name check
  //  std::cout <<"end of the menu "<<std::endl;
  std::sort(eeDetIdToTowerId_.begin(),eeDetIdToTowerId_.end(),boost::bind( std::less<int>(),boost::bind(&std::pair<int,int>::first,_1),boost::bind(&std::pair<int,int>::first,_2)));

  //std::cout <<"here "<<std::endl;
  //for(size_t i=0;i<eeDetIdToTowerId_.size();i++){
    //std::cout <<" index "<<i<<" val "<<eeDetIdToTowerId_[i].first<<std::endl;
    
    //}
}
  
std::vector<std::vector<int> > DetIdTools::makeTowerToRecHitsHashTable_()
{
  std::vector<std::vector<int> > hashTable(kNrCaloTowers);
  
  auto filler=[&hashTable](int towerId,int hitId){
    size_t towerHash = calHashCalo(towerId);
    if(towerHash<hashTable.size()){
      hashTable[towerHash].push_back(hitId);
    }else {
      LogErr<<" error hash of "<<towerId<<" is "<<towerHash<<" which is larger than hashTable "<<hashTable.size()<<std::endl;
    }
  };

  for(int iEta=-1*kMaxIEtaBarrel;iEta<=kMaxIEtaBarrel;iEta++){
    for(int iPhi=kMinIPhiBarrel;iPhi<=kMaxIPhiBarrel;iPhi++){
      if(isValidEcalBarrelId(iEta,iPhi)){
	int id =makeEcalBarrelId(iEta,iPhi);
	int towerId = ecalToTowerId(id);
	if(towerId!=0) filler(towerId,id);
      }
    }
  }
 
  for(int iX=kMinIXEndcap;iX<=kMaxIXEndcap;iX++){
    for(int iY=kMinIYEndcap;iY<=kMaxIYEndcap;iY++){
      for(int iZ=-1;iZ<=1;iZ+=2){
	if(isValidEcalEndcapId(iX,iY,iZ)){
	  int id =makeEcalEndcapId(iX,iY,iZ);
	  int towerId = ecalToTowerId(id);
	  if(towerId!=0) filler(towerId,id);
	}
      }
    }
  }

  for(int iEta=-1*kHcalIEtaAbsMax;iEta<=kHcalIEtaAbsMax;iEta++){
    for(int iPhi=kHcalIPhiMin;iPhi<=kHcalIPhiMax;iPhi++){
      for(int depth=kHcalDepthMin;depth<=kHcalDepthMax;depth++){
	if(isValidHcalId(iEta,iPhi,depth)){
	  int id = makeHcalDetId(iEta,iPhi,depth);
	  int towerId= makeCaloDetId(iEta,iPhi);
	  if(towerId!=0 && isValidCaloId(iEta,iPhi)) filler(towerId,id);
	}
      }
    }
  }
  return hashTable;
}
