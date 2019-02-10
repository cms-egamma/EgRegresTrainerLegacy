
#ifndef EGAMMAOBJECTS_GBRTrainer
#define EGAMMAOBJECTS_GBRTrainer

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// GBRTrainer                                                           //
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

#include <TRandom3.h>

  class GBRForest;
  class GBRTree;
  class TTree;
  class TCut;
  class GBREvent;
  
  class GBRTrainer {
    
    public:

       GBRTrainer();
       ~GBRTrainer();
       
       void AddInputVar(std::string var)    { fInputVars.push_back(var); }
       void SetTargetVar(std::string var)   { fTargetVar = var;          }
       //void SetTree(TTree *tree)            { fTree = tree;              }
       void AddTree(TTree *tree, double w=1.0) { fTrees.push_back(tree); fTreeWeights.push_back(w); }
       void SetTrainingCut(std::string cut) { fTrainingCut = cut;        }
       void SetMinEvents(int n)             { fMinEvents = n;            }
       void SetShrinkage(float x)           { fShrinkage = x;            }
       void SetMinCutSignificance(float x)  { fMinCutSignificance = x;   }
       void SetTransitionQuantile(float x)  { fTransitionQuantile = x;   }

       void SetRandomSeed(std::string seed) { fRandomSeedFormula = seed; }
       void SetEventWeight(std::string w)   { fEventWeightFormula = w;  }

       void PrintInputVars()
       {
           printf("Input variables:\n");
           for(unsigned int v=0; v<fInputVars.size(); v++)
            {
                printf("Var %i = %s\n", v, fInputVars[v].c_str());
            }
       }
      
       const GBRForest *TrainForest(int ntrees);
       
    protected:

       //float WeightedMedian(std::vector<GBREvent*> &evts);
      
      void TrainTree(const std::vector<GBREvent*> &evts, double sumwtotal, GBRTree &tree, int nvars, double transition);      
      void BuildLeaf(const std::vector<GBREvent*> &evts, double sumw, GBRTree &tree, double transition);
      
      std::vector<TTree*>       fTrees;
      std::vector<double>       fTreeWeights;
      std::string               fTrainingCut;
      std::vector<std::string>  fInputVars;  
      std::string               fTargetVar;
      int                       fMinEvents;
      float                     fShrinkage;
      //std::vector<std::vector<float> > fQuantileMaps;
      int                       fNQuantiles;
      unsigned int              fNBinsMax;
      float                     fTransitionQuantile;
      float                     fMinCutSignificance;

      std::string               fRandomSeedFormula;
      std::string               fEventWeightFormula;
      TRandom3                  fRandom;
      
      
      float *_sepgains;
      float *_sepgainsigs;
      float *_cutvals;  
      int *_nlefts;
      int *_nrights;
      float *_sumwlefts;
      float *_sumwrights;  
      float *_sumtgtlefts;
      float *_sumtgtrights;
      float *_leftvars;
      float *_rightvars;      
      float *_fullvars;
      int   *_bestbins;
      
      
      float **_ws;
      float **_ws2;
      int **_ns;
      float **_tgts;
      float **_tgt2s;
      float **_sumws;
      float **_sumws2;      
      int **_sumns;
      float **_sumtgts;
      float **_sumtgt2s;
      float **_varvals;
      float **_bsepgains;
      float **_bsepgainsigs;      
      
      int **_quants;
      int **_bins;
      
      float **fQuantileMaps;
      
  };
#endif
