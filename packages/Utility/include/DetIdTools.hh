#ifndef DETIDTOOLS
#define DETIDTOOLS

//this class defines tools for manuiplating ints storing det id information
//in theory this will become a class one day

#include <iostream>
#include <vector>
#include <unordered_map>
#include <cstdlib>

class DetIdTools {

public:
  //a class that allows us to navigate the ecal
  class EcalNavigator {
  private:
    int startId_; //the initial position
    mutable int currEtaOrX_; //where we are currently
    mutable int currPhiOrY_;
    mutable int zSide_;
    bool isBarrel_; //if true, we are navigating barrel, else endcap
  public:
    EcalNavigator(int startId);
    ~EcalNavigator(){}
    
    //north/south is phi for barrel and y for endcap (north:increasing phi/Y
    //east/west is eta (eb) and x (ee), east is +pos

    //all return the detId or 0 if invalid
    int goToStart(){initPosCoords_();return startId_;}
    int goNorth(int nrSteps=1); //negative values are allowed, -2 north = 2 south
    int goSouth(int nrSteps=1){return goNorth(nrSteps*-1);}
    int goEast(int nrSteps=1);
    int goWest(int nrSteps=1){return goEast(nrSteps*-1);}
    int getIdAtPos(int nrStepsEast,int nrStepsNorth)const;
    int curPos()const; //where we currently are
    void look()const{std::cout <<"you are surrounded by crystals, in the distance, you see a shambling electron clawing its way out of the nearby chasm"<<std::endl;}

  private:
    void initPosCoords_();
    int getDetId_(int etaOrX,int phiOrY)const;

  };

  class CaloNavigator {
  private:
    int startId_; //the initial position
    int currEta_; //where we are currently
    int currPhi_;
    bool l1Navigator_; //if its a L1 navigator (L1 keeps 72 towers in phi upto iEta of 28)
  public:
    const static int kIEtaPhiSegChange=20;
    const static int kNrPhiSeg=72;

  public:
    CaloNavigator(int startId,bool l1Navigator=false);
    ~CaloNavigator(){}
    
    //north/south is phi for barrel and y for endcap (north:increasing phi/Y
    //east/west is eta (eb) and x (ee), east is +pos

    //all return the detId or 0 if invalid
    int moveToStart(){initPosCoords_();return startId_;}
    int moveInEta(int nrSteps=1);
    int moveInPhi(int nrSteps=1);
    int getIdAtPos(int nrStepsEta,int nrStepsPhi)const;
    int curPos()const{return l1Navigator_ ? DetIdTools::makeL1CaloDetId(currEta_,currPhi_) : DetIdTools::makeCaloDetId(currEta_,currPhi_);}
 
    static int offsetInEta(int startEta,int nrSteps);
    static int offsetInL1Phi(int startPhi,int nrSteps);

  private:
    void initPosCoords_();
    int changedPhiSeg(int newEta,int oldEta)const{ return ( (abs(newEta)<=kIEtaPhiSegChange) != (abs(oldEta)<=kIEtaPhiSegChange) ) && !l1Navigator_;} //l1 towers dont change in phi seg (atleast in HB, HE
    
  };

  
  

 private:
  DetIdTools();
  virtual ~DetIdTools(){}

private:
  static const int kIXMax_=100;
  static const int kICrFee_ = 7740;
  static const int nBegin_[kIXMax_];
  static const int nIntegral_[kIXMax_];

  //my hashing function is a little slow
  //this is a bandaid, a precomputed tables of all the hashes
  //it keys of the detId - detector+subdetector code
  //it wastes a little memory but gives fast results (1mb for hcal, 0.5mb for eb, 128k for ee, I think its well worth it)
  static const std::vector<int> ebFastHashTable_;
  static const std::vector<int> eeFastHashTable_;
  static const std::unordered_map<int,int> hcalFastHashTable_;

  static std::vector<std::pair<int,int> > eeDetIdToTowerId_;

  static const std::vector<std::vector<int> > towerToRecHits_;
  static const std::vector<int> dummyTowerHits_;

 public:
  static const int kDetOffset = 28;
  static const int kSubDetOffset = 25;

  static const int kDetMask = 0xF << kDetOffset;
  static const int kSubDetMask  = 0x7 << kSubDetOffset;
  static const int kBarrelCode = 1 << kSubDetOffset;
  static const int kEndcapCode = 2 << kSubDetOffset;
  static const int kTrackerCode = 1 << kDetOffset;
  static const int kEcalCode = 3 << kDetOffset;
  static const int kHcalCode = 4 << kDetOffset;
  static const int kCaloCode = 5 << kDetOffset;

  static const int kNrCrysBarrel = 61200;
  static const int kNrCrysEndcap = 7740*2;

  static const int kMaxIEtaBarrel = 85;
  static const int kMaxIPhiBarrel = 360;
  static const int kMinIEtaBarrel = 1;
  static const int kMinIPhiBarrel = 1;

  static const int kMinIXEndcap = 1;
  static const int kMaxIXEndcap = 100;
  static const int kMinIYEndcap = 1;
  static const int kMaxIYEndcap = 100;

  static const int kNrHcalCellsBarrel = 1296*2;
  static const int kNrHcalCellsEndcap = 3384*2;
  static const int kNrEcalCellsBarrel = 61200; 
  static const int kNrEcalCellsEndcap = 7740*2; 
  static const int kNrCaloTowers = 1800*2; //upto and including eta of 29

  static const int kHcalIPhiMin = 1;
  static const int kHcalIPhiMax = 72;
  static const int kHcalIEtaAbsMin = 1;
  static const int kHcalIEtaAbsMax = 29;
  static const int kHcalDepthMin = 1;
  static const int kHcalDepthMax = 7;

  static const int kHcalPhiMask1       = 0x7F;
  static const int kHcalPhiMask2       = 0x3FF;
  static const int kHcalEtaOffset1     = 7;
  static const int kHcalEtaOffset2     = 10;
  static const int kHcalEtaMask1       = 0x3F;
  static const int kHcalEtaMask2       = 0x1FF;
  static const int kHcalZsideMask1     = 0x2000;
  static const int kHcalZsideMask2     = 0x80000;
  static const int kHcalDepthOffset1   = 14;
  static const int kHcalDepthOffset2   = 20;
  static const int kHcalDepthMask1     = 0x1F;
  static const int kHcalDepthMask2     = 0xF;
  static const int kHcalDepthSet1      = 0x1C000;
  static const int kHcalDepthSet2      = 0xF00000;
  static const int kHcalIdFormat2      = 0x1000000;
  static const int kHcalIdMask         = 0xFE000000;



  static const int kNrL1CaloTowers = 64*72; //not strictly true, its actually 28*72 + 4*16 (HF towers have less phi segmentation) but easier for maths right now
  static const int kL1CaloIPhiMin =1;
  static const int kL1CaloIPhiMax =72;
  static const int kL1CaloIEtaAbsMin = 1;
  static const int kL1CaloIEtaAbsMax = 32;

  static const int kBPXLayerOffset     = 20;
  static const int kBPXLadderOffset    = 12;
  static const int kBPXModuleOffset    = 2;
  static const int kBPXLayerMask       = 0xF;
  static const int kBPXLadderMask      = 0xFF;
  static const int kBPXModuleMask      = 0xFF;
  
  static const int kFPXSideOffset      = 23;
  static const int kFPXDiskOffset      = 18;
  static const int kFPXBladeOffset     = 12;   
  static const int kFPXPanelOffset     = 10; 
  static const int kFPXModuleOffset    = 2;   
  static const int kFPXSideMask        = 0x3;
  static const int kFPXDiskMask        = 0xF;
  static const int kFPXBladeMask       = 0x3F;   
  static const int kFPXPanelMask       = 0x3; 
  static const int kFPXModuleMask      = 0xFF; 


  //detector id
  static bool isEcal(int detId){return (detId&kDetMask) == kEcalCode;}
  static bool isHcal(int detId){return (detId&kDetMask) == kHcalCode;}
  static bool isCalo(int detId){return (detId&kDetMask) == kCaloCode;}
  static bool isTracker(int detId){return (detId&kDetMask) == kTrackerCode;}
  static bool isBarrel(int detId){return (detId&kSubDetMask) == kBarrelCode;}
  static bool isEndcap(int detId){return (detId&kSubDetMask) == kEndcapCode;}
  static bool isEcalBarrel(int detId){return isEcal(detId) && isBarrel(detId);}
  static bool isEcalEndcap(int detId){return isEcal(detId) && isEndcap(detId);}
  static bool isHcalBarrel(int detId){return isHcal(detId) && isBarrel(detId);}
  static bool isHcalEndcap(int detId){return isHcal(detId) && isEndcap(detId);}
  static int detSubDet(int detId){return detId&(kDetMask | kSubDetMask);}

  //to make det ids
  static int makeEcalBarrelId(int iEta,int iPhi);
  static int makeEcalEndcapId(int ix,int iy,int iz);
  static int makeHcalDetIdOld(int subDetCode,int iEta,int iPhi,int depth);
  static int makeHcalDetIdOld(int iEta,int iPhi,int depth);
  static int makeHcalBarrelDetIdOld(int iEta,int iPhi,int depth){return makeHcalDetIdOld(kBarrelCode,iEta,iPhi,depth);}
  static int makeHcalEndcapDetIdOld(int iEta,int iPhi,int depth){return makeHcalDetIdOld(kEndcapCode,iEta,iPhi,depth);}
  static int makeHcalDetId(int subDetCode,int iEta,int iPhi,int depth);
  static int makeHcalDetId(int iEta,int iPhi,int depth);
  static int makeHcalBarrelDetId(int iEta,int iPhi,int depth){return makeHcalDetId(kBarrelCode,iEta,iPhi,depth);}
  static int makeHcalEndcapDetId(int iEta,int iPhi,int depth){return makeHcalDetId(kEndcapCode,iEta,iPhi,depth);}
  

  static int makeCaloDetId(int iEta,int iPhi);
  static int makeL1CaloDetId(int iEta,int iPhi);

  //check ids are valid
  static bool isValidEcalBarrelId(int iEta,int iPhi);
  static bool isValidEcalBarrelId(int detId){return isEcalBarrel(detId) && isValidEcalBarrelId(iEtaBarrel(detId),iPhiBarrel(detId));}
  static bool isValidEcalEndcapId(int crystal_ix,int crystal_iy,int iz);
  static bool isValidEcalEndcapId(int detId){return isEcalEndcap(detId) && isValidEcalEndcapId(iXEndcap(detId),iYEndcap(detId),zEndcap(detId));}
  static bool isValidEcalId(int detId);

  static bool isValidPhase0HcalId(int iEta,int iPhi,int depth); 
  static bool isValidPhase0HcalId(int detId){return isHcal(detId)&& isValidPhase0HcalId(iEtaHcal(detId),iPhiHcal(detId),depthHcal(detId));}
  static bool isValidPhase0HcalBarrelId(int iEta,int iPhi,int depth);
  static bool isValidPhase0HcalEndcapId(int iEta,int iPhi,int depth);
  
  static bool isValidPhase1HcalId(int iEta,int iPhi,int depth); 
  static bool isValidPhase1HcalId(int detId){return isHcal(detId)&& isValidPhase1HcalId(iEtaHcal(detId),iPhiHcal(detId),depthHcal(detId));}
  static bool isValidPhase1HcalId(int subdetId,int iEta,int iPhi,int depth);
  static bool isValidPhase1HcalBarrelId(int iEta,int iPhi,int depth){return isValidPhase0HcalBarrelId(iEta,iPhi,depth);} //phase0 and phase 1 are identical for HB
  static bool isValidPhase1HcalEndcapId(int iEta,int iPhi,int depth); 

  static bool isValidHcalId(int iEta,int iPhi,int depth){return isValidPhase1HcalId(iEta,iPhi,depth);}
  static bool isValidHcalId(int detId){return isValidPhase1HcalId(detId);}
  static bool isValidHcalId(int subdetId,int iEta,int iPhi,int depth){return isValidPhase1HcalId(subdetId,iEta,iPhi,depth);}
  static bool isValidHcalBarrelId(int iEta,int iPhi,int depth){return isValidPhase1HcalBarrelId(iEta,iPhi,depth);}
  static bool isValidHcalEndcapId(int iEta,int iPhi,int depth){return isValidPhase1HcalEndcapId(iEta,iPhi,depth);}


  static bool isValidCaloId(int detId){return isCalo(detId)&& isValidCaloId(iEtaCalo(detId),iPhiCalo(detId));}
  static bool isValidCaloId(int iEta,int iPhi);
  static bool isValidL1CaloId(int detId){return isCalo(detId)&& isValidL1CaloId(iEtaCalo(detId),iPhiCalo(detId));}
  static bool isValidL1CaloId(int iEta,int iPhi);

  //gap tools, maxDistToGap is number of crystals the gap is away (0=this crys is next to gap)
  static bool isNextToBarrelPhiGap(int detId);
  static bool isNextToBarrelEtaGap(int detId,int maxDistToGap=0);
  static int nrOfNearestGap(int detId);
 
  
  //ECAL barrel tools
  static int iEtaBarrel(int detId){return positiveZBarrel(detId) ? iEtaAbsBarrel(detId) : -1*iEtaAbsBarrel(detId) ;}
  static int iEtaAbsBarrel(int detId){return (detId>>9) & 0x7F ;}
  static int iEtaAbsBarrelTower(int detId){return (iEtaAbsBarrel(detId)-1)/5+1;}
  static int iEtaBarrelTower(int detId){return iEtaAbsBarrelTower(detId)*zSideBarrel(detId);}
  static int iPhiBarrel(int detId){return detId&0x1FF;}
  static int iPhiBarrelTower(int detId){int iPhi = (iPhiBarrel(detId)-1)/5+1;iPhi-=2;return iPhi<=0 ? iPhi+72 : iPhi;}
  static bool positiveZBarrel(int detId){return detId&0x10000;}
  static int zSideBarrel(int detId){return positiveZBarrel(detId) ? 1 : -1;}
  static std::ostream& printEBDetId(int detId,std::ostream& out);

  static int dIEtaBarrel(int iEta1,int iEta2);
  static int dAbsIEtaBarrel(int iEta1,int iEta2){return abs(dIEtaBarrel(iEta1,iEta2));}
  static int dIPhiBarrel(int iPhi1,int iPhi2);
  static int dAbsIPhiBarrel(int iPhi1,int iPhi2){return abs(dIEtaBarrel(iPhi1,iPhi2));}

  //ECAL endcap tools
  static int iYEndcap(int detId){return detId&0x7f;}
  static int iXEndcap(int detId){return (detId>>7)&0x7F;} 
  static int normEndcapIXOrIY(int iXorIY); //x/y distance from centre 
  static bool positiveZEndcap(int detId){return detId&0x4000;}
  static int zEndcap(int detId){return positiveZEndcap(detId) ? 1 : -1;}
  static int endcapEtaRing(int detId); //these tools calculate eta = sqrt(x*x+y*y)
  static int endcapEtaRing(int ix,int iy);
  static float endcapEtaRingFloat(int detId);
  static float endcapEtaRingFloat(int ix,int iy);
  static bool isNextToRingBoundary(int detId);
  
  static int towerIdEndcap(int detId);
  static int ecalToTowerId(int detId);
  //HCAL tools (unified for HB/HE)
  static bool newFormatHcal(int detId){return !oldFormatHcal(detId);}
  static bool oldFormatHcal(int detId){ return ((detId&kHcalIdFormat2)==0)?(true):(false); }
  static int convertToOldFormatHcal(int detId){ return oldFormatHcal(detId) ? detId : newToOldFormatHcal_(detId);}
  static int iPhiHcal(int detId){return newFormatHcal(detId) ? detId&kHcalPhiMask2 : detId&kHcalPhiMask1;}
  static int iEtaAbsHcal(int detId){return newFormatHcal(detId) ? (detId>>kHcalEtaOffset2)&kHcalEtaMask2 : (detId>>kHcalEtaOffset1)&kHcalEtaMask1;}
  static int iEtaHcal(int detId){return zSideHcal(detId)*iEtaAbsHcal(detId);}
  static int zSideHcal(int detId){int mask=newFormatHcal(detId) ? kHcalZsideMask2 : kHcalZsideMask1;return (detId&mask) ? (1) : (-1);}
  static int depthHcal(int detId){return newFormatHcal(detId) ? (detId>>kHcalDepthOffset2)&kHcalDepthMask2 : (detId>>kHcalDepthOffset1)&kHcalDepthMask1;}
  static int getEffectiveHcalDepth(int detId); //effective depth is for isolation and accounts for not all depths being equal
  static int getNrDepthsInHcalTower(int detId); 
  static void getMatchingIdsHcal(int etaAbs,int phi,int side,int depth,std::vector<int>& ids);
  static void printHcalDetId(int detId);

  //CaloTower tools
  static int iPhiCalo(int detId){return detId&0x7F;}
  static int iEtaAbsCalo(int detId){return (detId>>7)&0x3f;}
  static int zSideCalo(int detId){return (detId&0x2000) ? (1) : (-1);}
  static int iEtaCalo(int detId){return zSideCalo(detId)*iEtaAbsCalo(detId);}
  
  //preshower tools
  static int zSideES(int detId){return (detId&0x80000) ? 1 : -1;}
  static int planeES(int detId){return ((detId>>18)&0x1) + 1;}
  static int sixES(int detId){return (detId>>6)&0x3F;}
  static int siyES(int detId){return (detId>>12)&0x3F;}
  static int stripES(int detId){return detId&0x3F;}

  //pixel tools
  static int getVal(int detId,int offset,int mask){return (detId>>offset)&mask;}
  static int layerBPX(int detId){return getVal(detId,kBPXLayerOffset,kBPXLayerMask);}
  static int ladderBPX(int detId){return getVal(detId,kBPXLadderOffset,kBPXLadderMask);}
  static int moduleBPX(int detId){return getVal(detId,kBPXModuleOffset,kBPXModuleMask);}
  static int sideFPX(int detId){return getVal(detId,kFPXSideOffset,kFPXSideMask);} 
  static int diskFPX(int detId){return getVal(detId,kFPXDiskOffset,kFPXDiskMask);} 
  static int bladeFPX(int detId){return getVal(detId,kFPXBladeOffset,kFPXBladeMask);} 
  static int panelFPX(int detId){return getVal(detId,kFPXPanelOffset,kFPXPanelMask);} 
  static int moduleFPX(int detId){return getVal(detId,kFPXModuleOffset,kFPXModuleMask);} 

  
  //hashes for fast lookup

  //this functions convert an ECAL detId into an array index
  //note barrel is first then endcap
  static int calHashEcal(int detId){return isEcal(detId) ? isBarrel(detId) ? calHashEcalBarrel(detId) : calHashEcalEndcap(detId) + kNrEcalCellsBarrel : 0;}
  static int calHashEcalEndcap(int detId) {
    return iYEndcap(detId) - nBegin_[iXEndcap(detId)-1] + nIntegral_[iXEndcap(detId) -1 ] + (positiveZEndcap(detId) ? kICrFee_ : 0);
  }
  static int calHashEcalBarrel(int detId){
    int etaBandNr =  kMaxIEtaBarrel + (positiveZBarrel(detId) ? iEtaAbsBarrel(detId)-1 : -iEtaAbsBarrel(detId)); // 0 - 189 starting at lowest eta (~-1.5) to highest
    return etaBandNr* kMaxIPhiBarrel + iPhiBarrel(detId)-1;
  }
  static int calHashHcalLegacy(int detId);

  static int calHashHcal(int detId);
  
  static int calHashCalo(int detId);

  static int calHashL1Calo(int iEta,int iPhi);

  static int getHashHcal(int detId){
    auto res = hcalFastHashTable_.find(detId);
    if(res!=hcalFastHashTable_.end()) return res->second;
    else return -1;
  }
  static int getHashEcal(int detId){return isEcal(detId) ? isBarrel(detId) ? getHashEcalBarrel(detId) : getHashEcalEndcap(detId) + kNrEcalCellsBarrel : 0;}
  static int getHashEcalBarrel(int detId){return ebFastHashTable_[detId & ~(kDetMask | kSubDetMask)];}
  static int getHashEcalEndcap(int detId){return eeFastHashTable_[detId & ~(kDetMask | kSubDetMask)];}


  static void fillEEToTowerIdMap(const std::string& filename);
    
  static const std::vector<int>& getRecHitsOfTower(int towerId); //this is definately bugged for iEta 27/28/29 and possibly bugged for the endcap

private:
  static int newToOldFormatHcal_(int detId);
  //yes I'm aware these functions pass 1MB vectors by value, they
  //are only called once and its the only way I can see
  static std::vector<int> makeEBFastHashTable_();
  static std::vector<int> makeEEFastHashTable_();
  static std::unordered_map<int,int> makeHcalFastHashTable_();
  static std::vector<std::vector<int> > makeTowerToRecHitsHashTable_();

};

#endif
