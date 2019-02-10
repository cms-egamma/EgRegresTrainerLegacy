
#ifndef GBRLIKELIHOOD_HybridGBRForestD
#define GBRLIKELIHOOD_HybridGBRForestD

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// HybridGBRForestD                                                            //
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

  class HybridGBRForestD {

    public:

       HybridGBRForestD() {}      
       HybridGBRForestD(int ntargets);
       virtual ~HybridGBRForestD();
       
       int NTargets() const { return fTrees.size(); }
       
       //void GetResponse(const float* vector) const;
       //double GetResponse(int idx) const { return fResponses[idx]; }
       
       double GetResponse(const float* vector, int idx) const;
       
       void SetInitialResponse(int idx, double response) { fInitialResponse[idx] = response; }
       
       std::vector<std::vector<HybridGBRTreeD> > &Trees() { return fTrees; }
       const std::vector<std::vector<HybridGBRTreeD> > &Trees() const { return fTrees; }
       
    protected:
      std::vector<double> fInitialResponse;
      //mutable std::vector<double> fResponses;
      std::vector<std::vector<HybridGBRTreeD> > fTrees;  
      
    private:

      ClassDef(HybridGBRForestD,2)       
      
  };

//_______________________________________________________________________
// inline void HybridGBRForestD::GetResponse(const float* vector) const {
//   for (unsigned int itgt=0; itgt<fResponses.size(); ++itgt) {
//     fResponses[itgt] = fInitialResponse[itgt];
//     for (std::vector<HybridGBRTreeD>::const_iterator it=fTrees[itgt].begin(); it!=fTrees[itgt].end(); ++it) {
//       int termidx = it->TerminalIndex(vector);
//       fResponses[itgt] += it->GetResponse(termidx);
//     }    
//   }
// }

//_______________________________________________________________________
inline double HybridGBRForestD::GetResponse(const float* vector, int idx) const {
  double response = fInitialResponse[idx];
  for (std::vector<HybridGBRTreeD>::const_iterator it=fTrees[idx].begin(); it!=fTrees[idx].end(); ++it) {
    int termidx = it->TerminalIndex(vector);
    response += it->GetResponse(termidx);
  }    
  return response;
}


#endif
