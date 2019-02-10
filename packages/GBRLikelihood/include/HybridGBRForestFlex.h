
#ifndef GBRLIKELIHOOD_HybridGBRForestFlex
#define GBRLIKELIHOOD_HybridGBRForestFlex

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// HybridGBRForestFlex                                                            //
//                                                                      //
// A fast minimal implementation of Gradient-Boosted Regression Trees   //
// which has been especially optimized for size on disk and in memory.  //                                                                  
//                                                                      //
// Designed to be built from TMVA-trained trees, but could also be      //
// generalized to otherwise-trained trees, classification,              //
//  or other boosting methods in the future                             //
//                                                                      //
//  Josh Bendavid - MIT                                                 //
//////////////////////////////////////////////////////////////////////////

#include <vector>
#include "HybridGBRTreeD.h"
#include <math.h>
#include <stdio.h>
#include "Rtypes.h"

  class HybridGBRForestFlex {
    public:
       typedef HybridGBRTreeD TreeT;

       HybridGBRForestFlex();   
       virtual ~HybridGBRForestFlex();
       
       int NTargets() const { return fTrees.size(); }
       
       double GetResponse(const float* vector) const;
       
       double InitialResponse() const { return fInitialResponse; }
       void SetInitialResponse(double response) { fInitialResponse = response; }
       
       std::vector<HybridGBRTreeD> &Trees() { return fTrees; }
       const std::vector<HybridGBRTreeD> &Trees() const { return fTrees; }
       
    protected:
      double fInitialResponse;
      std::vector<HybridGBRTreeD> fTrees;  
      
    private:

      ClassDef(HybridGBRForestFlex,1)       
      
  };

//_______________________________________________________________________
inline double HybridGBRForestFlex::GetResponse(const float* vector) const {
  double response = fInitialResponse;
  for (std::vector<HybridGBRTreeD>::const_iterator it=fTrees.begin(); it!=fTrees.end(); ++it) {
    int termidx = it->TerminalIndex(vector);
    response += it->GetResponse(termidx);
  }    
  return response;
}


#endif
