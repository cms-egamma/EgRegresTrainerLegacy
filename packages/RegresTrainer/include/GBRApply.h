
#ifndef EGAMMAOBJECTS_GBRApply
#define EGAMMAOBJECTS_GBRApply

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// GBRApply                                                           //
//                                                                      //
// A fast minimal implementation of Gradient-Boosted Regression Trees   //
// which has been especially optimized for size on disk and in memory.  //                                                                  
//                                                                      //
// This class implements a cpu-time optimized training routine          //
// designed especially for regression training with a very large,       //
//  number of input variables and large training samples                //
//                                                                      //
// Expensive loops are optimized to ensure vectorization by the         //
// compiler, and the pre-sorting, as well as  main split-finding        //
//  search is multithreaded using OpenMP at the level of the            //
// loop over input variables                                            //
//                                                                      //
// Split-finding algorithm is a reworked version of the fixed binning   //
// method used in TMVA.  Input variables are binned in terms of         //
// (weighted) quantiles of the training sample.  In this way the binning //
// is robust against outliers in the training sample, and much fewer    //
// bins are needed for the split search on each node.                   //
//                                                                      //
//  Josh Bendavid - MIT                                                 //
//////////////////////////////////////////////////////////////////////////

#include <vector>
#include <string>
#include <stdio.h>

  class GBRForest;
  class GBRForestD;
  class TTree;
  class ICorrection;
  
  class GBRApply {
    
    public:

       GBRApply();
       ~GBRApply();
       
       TTree *ApplyAsFriend(TTree *intree,
               const GBRForest *forestEB, const GBRForest *forestEE,
               const std::vector<std::string> &varsEB, const std::vector<std::string> &varsEE,
               const std::string& cutEB, const std::string& cutEE,
               std::string targetname="target",
               const ICorrection* correction=NULL) const;

       // TTree *ApplyAsFriendTransform(TTree *intree,
       //         const GBRForestD *forestEB, const GBRForestD *forestEE,
       //         const std::vector<std::string> &varsEB, const std::vector<std::string> &varsEE,
       //         const std::string& cutEB, const std::string& cutEE,
       //         std::string targetname,
       //         double low, double high) const;

       TTree *ApplyAsFriendTransform(TTree *intree,
               const GBRForestD *forestEB, const GBRForestD *forestEE,
               const std::vector<std::string> &varsEB, const std::vector<std::string> &varsEE,
               const std::string& cutEB, const std::string& cutEE,
               std::string targetname,
               double low, double high,
               const bool doEB
               ) const;
       
    protected:
      
  };
#endif
