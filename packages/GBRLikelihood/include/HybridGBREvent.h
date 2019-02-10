#ifndef EGAMMAOBJECTS_HybridGBREvent
#define EGAMMAOBJECTS_HybridGBREvent

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// HybridGBREvent                                                             //
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
#include "TVectorD.h"  
#include "TMatrixDSym.h"
  
  class HybridGBREvent {

    public:

       HybridGBREvent(int nvars, int ntargets = 1, int nparms=0); 
       HybridGBREvent(int nvars, int ntargets, const HybridGBREvent &other, int nparms=0) :
        fPdfVal(other.fPdfVal),
        fWeight(other.fWeight),
        fClass(other.fClass) {
	  
          fVars = new float[nvars];
          fQuantiles = new int[nvars];
          fTargets = new double[ntargets];
          fTransTargets = new float[ntargets];
	  fTransTargets2 = new float[ntargets];
          fDerivatives = new double[nparms];
	  fDerivatives2 = new double[nparms];
	  fCurrentNodes =  new unsigned int[ntargets];
	  
          for (int ivar=0; ivar<nvars; ++ivar) {
            fVars[ivar] = other.fVars[ivar];
            fQuantiles[ivar] = other.fQuantiles[ivar];
          }
          
          for (int itgt=0; itgt<ntargets; ++itgt) {
            fTargets[itgt] = other.fTargets[itgt];
            fTransTargets[itgt] = other.fTransTargets[itgt];
	    fTransTargets2[itgt] = other.fTransTargets2[itgt];	    
	    fCurrentNodes[itgt] = other.fCurrentNodes[itgt];
          }
          
          for (int iparm=0; iparm<nparms; ++iparm) {
	    fDerivatives[iparm] = other.fDerivatives[iparm];
	    fDerivatives2[iparm] = other.fDerivatives2[iparm];
	  }
          
        }
        
       ~HybridGBREvent();
       
       unsigned char Class() const { return fClass; }
       float Var(int i) const { return fVars[i]; }
       unsigned short Quantile(int i) const { return fQuantiles[i]; }   
       double Weight() const   { return fWeight;  }
       
       void SetClass(unsigned char i) { fClass = i; }
       void SetVar(int i, float x) { fVars[i] = x; }
       void SetQuantile(int i, int q) { fQuantiles[i] = q; }
       //cache computed qunatities needed for split-search
       void SetWeight(double x)     { fWeight = x; }
       
       double Target(int i) const { return fTargets[i]; }
       float TransTarget(int i) const { return fTransTargets[i]; }
       void SetTarget(int i, double x) { fTargets[i] = x; }
       void SetTransTarget(int i, float x) { fTransTargets[i] = x; }
       
       double PdfVal() const { return fPdfVal; }
       void SetPdfVal(double x) { fPdfVal = x; }

       double Derivative(int i) const { return fDerivatives[i];  }       
       void SetDerivative(int i, double x) { fDerivatives[i] = x; }     
       
       double Derivative2(int i) const { return fDerivatives2[i];  }       
       void SetDerivative2(int i, double x) { fDerivatives2[i] = x; }        

       double TransTarget2(int i) const { return fTransTargets2[i];  }       
       void SetTransTarget2(int i, double x) { fTransTargets2[i] = x; }        
       
       unsigned int CurrentNode(int i) const { return fCurrentNodes[i];  }       
       void SetCurrentNode(int i, unsigned int x) { fCurrentNodes[i] = x; }         
       
       
       
       
    protected:
      float                    *fVars;
      double                   *fTargets;
      float                    *fTransTargets;
      float                    *fTransTargets2;      
      int                      *fQuantiles;
      double                   *fDerivatives;
      double                   *fDerivatives2;
      unsigned int             *fCurrentNodes;
      double                    fPdfVal;
      double                    fWeight;
      unsigned char             fClass;
  };
  
  
/*  class GBRTargetCMP : public std::binary_function<HybridGBREvent*, HybridGBREvent*, bool> {
    public:
      GBRTargetCMP() {}
      bool operator() (const HybridGBREvent *ev1, const HybridGBREvent *ev2) const { return ev1->Target()<ev2->Target() ? true : false; }
  };
  
  class GBRAbsTargetCMP : public std::binary_function<HybridGBREvent*, HybridGBREvent*, bool> {
    public:
      GBRAbsTargetCMP() {}
      bool operator() (const HybridGBREvent *ev1, const HybridGBREvent *ev2) const { return fabs(ev1->Target())<fabs(ev2->Target()) ? true : false; }
  }; */ 
  
  class GBRVarCMP : public std::binary_function<HybridGBREvent*, HybridGBREvent*, bool> {
    public:
      GBRVarCMP() {}
      GBRVarCMP(int idx) : fVarIdx(idx) {}      
      bool operator() (const HybridGBREvent *ev1, const HybridGBREvent *ev2) const { return ev1->Var(fVarIdx)<ev2->Var(fVarIdx) ? true : false; }
      
    protected:
      int fVarIdx;
  };
  
  
#endif
