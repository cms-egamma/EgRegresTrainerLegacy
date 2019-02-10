
#ifndef EGAMMAOBJECTS_GBREvent
#define EGAMMAOBJECTS_GBREvent

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// GBREvent                                                             //
//                                                                      //
// A fast minimal implementation of Gradient-Boosted Regression Trees   //
// which has been especially optimized for size on disk and in memory.  //                                                                  
//                                                                      //
// This is a helper class for GBRTrainer to store  needed information   //
// in memory and facilitate sorting according to target or input        //
// variables.                                                           //
//                                                                      //
//  Josh Bendavid - MIT                                                 //
//////////////////////////////////////////////////////////////////////////

#include <vector>
#include <string>
#include <stdio.h>
#include <algorithm>
#include <cmath>
  
  class GBREvent {

    public:

       GBREvent(int nvars);       
       ~GBREvent();
       
       float Var(int i) const { return fVars[i]; }
       unsigned short Quantile(int i) const { return fQuantiles[i]; }
       float Target() const   { return fTarget;  }
       float TransTarget() const { return fTransTarget;  }       
       float Weight() const   { return fWeight;  }
       float WeightedTransTarget() const { return fWeightedTransTarget; }
       float WeightedTransTarget2() const { return fWeightedTransTarget2; }
       
       void SetVar(int i, float x) { fVars[i] = x; }
       void SetQuantile(int i, int q) { fQuantiles[i] = q; }
       void SetTarget(float x)     { fTarget = x; }
       //cache computed qunatities needed for split-search
       void SetTransTarget(float x) { fTransTarget = x; fWeightedTransTarget = fWeight*fTransTarget; fWeightedTransTarget2 = fWeightedTransTarget*fTransTarget; }
       void SetWeight(float x)     { fWeight = x; }
       
       
    protected:
      float                    *fVars;
      int                      *fQuantiles;
      float                     fTarget;
      float                     fTransTarget;
      float                     fWeight;
      float                     fWeightedTransTarget;
      float                     fWeightedTransTarget2;
  };
  
  
  class GBRTargetCMP : public std::binary_function<GBREvent*, GBREvent*, bool> {
    public:
      GBRTargetCMP() {}
      bool operator() (const GBREvent *ev1, const GBREvent *ev2) const { return ev1->Target()<ev2->Target() ? true : false; }
  };
  
  class GBRAbsTargetCMP : public std::binary_function<GBREvent*, GBREvent*, bool> {
    public:
      GBRAbsTargetCMP() {}
      bool operator() (const GBREvent *ev1, const GBREvent *ev2) const { return std::abs(ev1->Target())<std::abs(ev2->Target()) ? true : false; }
  };  
  
  class GBRVarCMP : public std::binary_function<GBREvent*, GBREvent*, bool> {
    public:
      GBRVarCMP() {}
      GBRVarCMP(int idx) : fVarIdx(idx) {}      
      bool operator() (const GBREvent *ev1, const GBREvent *ev2) const { return ev1->Var(fVarIdx)<ev2->Var(fVarIdx) ? true : false; }
      
    protected:
      int fVarIdx;
  };  
  
#endif
