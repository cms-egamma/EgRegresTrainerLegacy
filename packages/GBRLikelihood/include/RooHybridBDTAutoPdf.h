/*****************************************************************************
 * Project: RooFit                                                           *
 * Package: RooFitModels                                                     *
 *    File: $Id: HGGRooPdfs.h,v 1.1 2012/02/10 15:10:48 gpetrucc Exp $
 * Authors:                                                                  *
 *   WV, Wouter Verkerke, UC Santa Barbara, verkerke@slac.stanford.edu       *
 *   DK, David Kirkby,    UC Irvine,         dkirkby@uci.edu                 *
 *                                                                           *
 * Copyright (c) 2000-2005, Regents of the University of California          *
 *                          and Stanford University. All rights reserved.    *
 *                                                                           *
 * Redistribution and use in source and binary forms,                        *
 * with or without modification, are permitted according to the terms        *
 * listed in LICENSE (http://roofit.sourceforge.net/license.txt)             *
 *****************************************************************************/
#ifndef ROO_HYBRIDBDTAUTOPDF
#define ROO_HYBRIDBDTAUTOPDF

#include "RooAbsPdf.h"
#include "RooRealProxy.h"
#include "RooListProxy.h"
#include "TMatrixDSym.h"
#include "GBRLikelihood/HybridGBREvent.h"




class RooRealVar;
class RooAbsReal;
class RooArgSet;
class HybridGBRForestD;
class HybridGBRTreeD;
class HybridGBRForestFlex;
class TNtuple;

class RooTreeConvert {
  
public:
  RooTreeConvert() {}
  virtual ~RooTreeConvert() {}
  
  static RooDataSet *CreateDataSet(std::string name, TTree *tree, std::vector<std::string> vars, std::string weight);
  static RooDataSet *CreateDataSet(std::string name, TTree *tree, RooArgList &vars, RooRealVar &weight, bool limitvals = true);
  
  
private:
  ClassDef(RooTreeConvert,1)  
  
};

class RooNormPdf : public RooAbsReal {
 
public:
  RooNormPdf() {}
  RooNormPdf(const char *name, const char *title, RooAbsPdf &pdf, const RooArgSet &forcednormset);
  RooNormPdf(const RooNormPdf &other, const char* name=0);
  virtual ~RooNormPdf() {}
  
  virtual TObject* clone(const char* newname) const { return new RooNormPdf(*this,newname); }
  
  
protected:
  Double_t evaluate() const;
  
  RooRealProxy _pdf;
  RooListProxy _forcednormset;
  
  
private:
  ClassDef(RooNormPdf,1)    
  
};

class RooRealConstraint : public RooAbsReal {
 
public:
  RooRealConstraint() : _low(0.), _high(1.) {}
  RooRealConstraint(const char *name, const char *title, RooAbsReal &real, double low, double high);
  RooRealConstraint(const RooRealConstraint &other, const char* name=0);
  virtual ~RooRealConstraint() {}
  
  virtual TObject* clone(const char* newname) const { return new RooRealConstraint(*this,newname); }
  
  
protected:
  Double_t evaluate() const;
  
  RooRealProxy _real;
  Double_t     _low;
  Double_t     _high;
  double        _scale;
  double        _offset;
  
  
private:
  ClassDef(RooRealConstraint,1)    
  
};


class RooPowerLaw : public RooAbsPdf {
  
public:
  RooPowerLaw() {}
  RooPowerLaw(const char *name, const char *title, RooAbsReal &x, RooAbsReal &p);
  RooPowerLaw(const RooPowerLaw& other, const char* name=0);
  virtual ~RooPowerLaw() {}
  
  virtual TObject* clone(const char* newname) const { return new RooPowerLaw(*this,newname); }
  
  //virtual Bool_t selfNormalized() const { return kTRUE; }
  
  //virtual ExtendMode extendMode() const { return CanBeExtended ; }  

//   virtual Double_t expectedEvents(const RooArgSet* nset) const { return _norm.arg().getVal(); }
//   virtual Double_t expectedEvents(const RooArgSet& nset) const { return _norm.arg().getVal(); }  
  
  virtual Int_t getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analVars, const char* rangeName = 0) const;
  virtual Double_t analyticalIntegral(Int_t code, const char* rangeName = 0) const;
  
protected:
  
  Double_t evaluate() const;
  
  
  RooRealProxy _x;
  RooRealProxy _p;  
 
  
private:
  ClassDef(RooPowerLaw,1)  
  
};

class RooCondAddPdf : public RooAbsPdf {
  
public:
  RooCondAddPdf() {}
  RooCondAddPdf(const char *name, const char *title, const RooArgList &pdfs, const RooArgList &coeffs);
  RooCondAddPdf(const RooCondAddPdf& other, const char* name=0);
  virtual ~RooCondAddPdf() {}
  
  virtual TObject* clone(const char* newname) const { return new RooCondAddPdf(*this,newname); }
  
  virtual Bool_t selfNormalized() const { return kTRUE; }
  
  //virtual ExtendMode extendMode() const { return CanBeExtended ; }  

//   virtual Double_t expectedEvents(const RooArgSet* nset) const { return _norm.arg().getVal(); }
//   virtual Double_t expectedEvents(const RooArgSet& nset) const { return _norm.arg().getVal(); }  
  
  virtual Int_t getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analVars, const char* rangeName = 0) const;
  virtual Double_t analyticalIntegral(Int_t code, const char* rangeName = 0) const;
  
  RooAbsReal *createCDF(const RooArgSet &iset, const RooArgSet &nset = RooArgSet());
  
protected:
  
  Double_t evaluate() const;
  
  
  RooListProxy _pdfs;
  RooListProxy _coeffs;  
  //RooRealProxy _norm;
  bool _selfnorm;
  
private:
  ClassDef(RooCondAddPdf,1)  
  
};

class RooCondRatioPdf : public RooAbsPdf {
  
public:
  RooCondRatioPdf() {}
  RooCondRatioPdf(const char *name, const char *title, RooAbsReal &ratio, RooAbsReal &pdfden);
  RooCondRatioPdf(const RooCondRatioPdf& other, const char* name=0);
  virtual ~RooCondRatioPdf() {}
  
  virtual TObject* clone(const char* newname) const { return new RooCondRatioPdf(*this,newname); }
  
  virtual Bool_t selfNormalized() const { return kTRUE; }
  
  //virtual ExtendMode extendMode() const { return CanBeExtended ; }  

//   virtual Double_t expectedEvents(const RooArgSet* nset) const { return _norm.arg().getVal(); }
//   virtual Double_t expectedEvents(const RooArgSet& nset) const { return _norm.arg().getVal(); }  
  
  //virtual Int_t getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analVars, const char* rangeName = 0) const;
  //virtual Double_t analyticalIntegral(Int_t code, const char* rangeName = 0) const;
  
protected:
  
  Double_t evaluate() const;
  
  RooRealProxy _ratio;
  RooRealProxy _pdfden;  
  
private:
  ClassDef(RooCondRatioPdf,1)  
  
};


class RooPdfAddReal : public RooAbsReal {
  
public:
  RooPdfAddReal() {}
  RooPdfAddReal(const char *name, const char *title, const RooArgList &pdfs, const RooArgList &coeffs);
  RooPdfAddReal(const RooPdfAddReal& other, const char* name=0);
  virtual ~RooPdfAddReal() {}
  
  virtual TObject* clone(const char* newname) const { return new RooPdfAddReal(*this,newname); }
  
  virtual Int_t getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analVars, const char* rangeName = 0) const;
  virtual Double_t analyticalIntegral(Int_t code, const char* rangeName = 0) const;
  
protected:
  
  Double_t evaluate() const;
  
  
  RooListProxy _pdfs;
  RooListProxy _coeffs;  
  //RooRealProxy _norm;
  
private:
  ClassDef(RooPdfAddReal,1)  
  
};

class RooGBRFunctionFlex : public RooAbsReal {
 
public:
  RooGBRFunctionFlex() : _forest(0) {}
  RooGBRFunctionFlex(const char *name, const char *title);
  RooGBRFunctionFlex(const RooGBRFunctionFlex& other, const char* name=0);
  virtual ~RooGBRFunctionFlex();
  
  virtual TObject* clone(const char* newname) const { return new RooGBRFunctionFlex(*this,newname); }
      
  HybridGBRForestFlex *Forest() { return _forest; }
  const HybridGBRForestFlex *Forest() const { return _forest; }
  
  void SetForest(HybridGBRForestFlex *forest);
  void ResetForest();
  
protected:
  virtual Double_t evaluate() const { return 0.; }
  
  HybridGBRForestFlex *_forest;
  
private:
  ClassDef(RooGBRFunctionFlex,1)

  
  
};

class RooGBRTargetFlex : public RooAbsReal {
  
public:
  RooGBRTargetFlex() {}
  RooGBRTargetFlex(const char *name, const char *title, RooGBRFunctionFlex &func, RooRealVar &var, const RooArgList &funcvars);
  RooGBRTargetFlex(const RooGBRTargetFlex& other, const char* name=0);
  
  virtual TObject* clone(const char* newname) const { return new RooGBRTargetFlex(*this,newname); }
  
  void SetUseFunc(bool b);
  void ClearFuncServers();
  
  RooRealVar *Var() { return (RooRealVar*)(&_var.arg()); }
  const RooArgList &FuncVars() const { return _funcvars; }  
  
  RooGBRFunctionFlex *Func() { return (RooGBRFunctionFlex*)(_func.absArg()); }
  HybridGBRForestFlex *Forest() { return Func()->Forest(); }
  
protected:
  virtual Double_t evaluate() const { return _usefunc ? EvalFunc() : _var.arg().getVal(); }
  //virtual Double_t evaluate() const { return _var.arg().getVal(); }
  
  double EvalFunc() const;
  
  RooArgProxy _func;
  int _itgt;
  RooRealProxy _var;
  bool _usefunc;
  
  RooListProxy _funcvars;  
  mutable std::vector<float> _eval;
  
private:
    ClassDef(RooGBRTargetFlex,1)
  
  
  
};


class RooGBRFunction : public RooAbsReal {
 
public:
  RooGBRFunction() : _forest(0) {}
  RooGBRFunction(const char *name, const char *title, const RooArgList &vars, int ntargets);
  RooGBRFunction(const RooGBRFunction& other, const char* name=0);
  virtual ~RooGBRFunction();
  
  virtual TObject* clone(const char* newname) const { return new RooGBRFunction(*this,newname); }
  
  
  double GetResponse(int itgt) const;
  
  const RooArgList &Vars() const { return _vars; }
  
  
  HybridGBRForestD *Forest() { return _forest; }
  void SetForest(HybridGBRForestD *forest);
  
protected:
  virtual Double_t evaluate() const { return 0.; }
  
  
  RooListProxy _vars;
  HybridGBRForestD *_forest;
  mutable std::vector<float> _eval;
  
private:
  ClassDef(RooGBRFunction,2)

  
  
};

class RooGBRTarget : public RooAbsReal {
  
public:
  RooGBRTarget() {}
  RooGBRTarget(const char *name, const char *title, RooGBRFunction &func, int itgt, RooRealVar &var);
  RooGBRTarget(const RooGBRTarget& other, const char* name=0);
  
  virtual TObject* clone(const char* newname) const { return new RooGBRTarget(*this,newname); }
  
  void SetUseFunc(bool b);
  
  RooRealVar *Var() { return (RooRealVar*)(&_var.arg()); }
  int Index() const { return _itgt; }
  
protected:
  virtual Double_t evaluate() const { return _usefunc ? static_cast<RooGBRFunction*>(_func.absArg())->GetResponse(_itgt) : _var.arg().getVal(); }
  //virtual Double_t evaluate() const { return _var.arg().getVal(); }
  
  RooArgProxy _func;
  int _itgt;
  RooRealProxy _var;
  bool _usefunc;

  
private:
    ClassDef(RooGBRTarget,1)
  
  
  
};

class RooHybridBDTAutoPdf : public TNamed {
public:
  
  //RooHybridBDTAutoPdf() {} ;
  RooHybridBDTAutoPdf(const char *name, const char *title, const RooArgList &tgtvars, RooAbsReal &n0, RooRealVar &r, const std::vector<RooAbsData*> &data, const std::vector<RooAbsReal*> &pdfs);

  ~RooHybridBDTAutoPdf();

  
  void AddInputVar(std::string var)    { fInputVars.push_back(var); }
  void SetTargetVar(std::string var)   { fTargetVar = var;          }
  //void SetTree(TTree *tree)            { fTree = tree;              }
  //void AddTree(TTree *tree, double w=1.0) { fTrees.push_back(tree); fTreeWeights.push_back(w); }
  void SetTrainingCut(std::string cut) { fTrainingCut = cut;        }
  void SetMinEvents(int n)             { fMinEvents = n;            }
  void SetShrinkage(double x)          { fShrinkage = x;            }
  void SetMinCutSignificance(double x); //  { fMinCutSignificance = x*x/2.0; }
  void SetMaxNSpurious(double x) { fMaxNSpurious = x; }
  void SetTransitionQuantile(float x)  { fTransitionQuantile = x;   }
  void SetMinWeights(const std::vector<double> &minweights) { fMinWeights = minweights; }
  void SetMinWeightTotal(double x) { fMinWeightTotal = x; }
  void SetMaxDepth(int depth) { fMaxDepth = depth; } 
  void SetMaxNodes(int max) { fMaxNodes = max; }
  void SetPrescaleInit(int n) { fPrescaleInit = n; }
  void SetDoInitialFit(bool b) { fDoInitialFit = b; }
 
  void TrainForest(int ntrees, bool reuseforest = false);  
  
  void fitWithMinosFast();
  void fitWithMinos();
  
  double DLdR() const { return fdLdR; }
  
  TNtuple *ResTree() { return fResTree; } 
  
  TGraph *DrvGraph() { return fDrvGraph; }
  TGraph *DrvGraphSmooth() { return fDrvGraphSmooth; }
  TGraph *GraphDelta() { return fGraphDelta; }
  
  double RMin() const { return fRMin; }
  double RHigh() const { return fRHigh; }
  double RLow() const { return fRLow; }
  
protected:

//   typedef alignas(32) int int;
//   typedef alignas(32) float float;
//   typedef alignas(32) float double;

//   typedef int __attribute__ ((aligned (32))) int;
//   typedef float __attribute__ ((aligned (32))) float;
//   typedef double __attribute__ ((aligned (32))) double;

//   typedef int int;
//   typedef float float;
//   typedef double double;  
  
  void BuildQuantiles(int nvars, double sumabsw);
  void UpdateTargets(int nvars, int selvar);
  void FillDerivatives();
  
  void TrainTree(const std::vector<HybridGBREvent*> &evts, double sumwtotal, HybridGBRTreeD &tree, double transition, int depth, std::vector<std::pair<float,float> > limits, int tgtidx);      
  void BuildLeaf(const std::vector<HybridGBREvent*> &evts, HybridGBRTreeD &tree, int tgtidx);
  void UpdateCurrentNodes(const std::vector<HybridGBREvent*> &evts, HybridGBRTreeD &tree, int tgtidx);
  

  //void FitResponses(const std::vector<HybridGBREvent*> &evts, double sumwtotal, HybridGBRTreeD &tree);
  void FitResponses(int selvar);  
  
  TMatrixD vmultT(const TVectorD &v, const TVectorD &vT) const;
  double vmult(const TVectorD &vT, const TVectorD &v) const;
  
  double findCrossing(RooRealVar &r, double level, double rStart, double rBound, double rerr);
  
  void RecomputeTargets();
  void GradientMinos();
  
  double EvalLossAvg();
  
  static double EvalLossNull(double dummy);
  
  double EvalLossRooFit();
  double EvalLoss(double lambda, const TVectorD &dL, int itree=-1);

  
  double Derivative1Fast(RooAbsReal *function, double currentval, RooRealVar *var, RooArgSet *nset=0, double step=1e-3);
  double Derivative2Fast(RooAbsReal *function, double currentval, RooRealVar *var, RooArgSet *nset=0, double step=1e-3);  
  double Derivative2Fast(RooAbsReal *function, RooRealVar *var1, RooRealVar *var2, RooArgSet *nset=0, double step1=1e-3, double step2=1e-3);  

  
  double Derivative1(RooAbsReal *function, RooRealVar *var, RooArgSet *nset=0, double step=1e-3);
  double Derivative2(RooAbsReal *function, RooRealVar *var, RooArgSet *nset=0, double step=1e-3);
  double Derivative2(RooAbsReal *function, RooRealVar *var1,RooRealVar *var2, RooArgSet *nset=0, double stepa=1e-3, double stepb=1e-3);
  
  
  RooArgList fCondVars;
  RooArgList fParmVars;
  RooArgList fTgtVars;
  RooArgList fExtVars;
  RooArgList fFullParms;
  RooArgList fFullFuncs;
  
  
  //RooGBRFunction *fFunc;

  TNtuple *fResTree;
  
  
  std::vector<RooAbsReal*> fPdfs;
  std::vector<RooAbsData*> fData;
  
  RooArgList fStaticTgts;
  RooArgList fStaticPdfs;
  RooArgList fFuncs;
  
  std::vector<RooArgList> fFuncTgts;
  std::vector<RooArgList> fTgtCondVars;
  std::vector<RooArgList> fFuncCondVars;
  
  std::vector<std::vector<RooRealVar*> > fCondVarsClones;
  std::vector<std::vector<RooRealVar*> > fParmVarsClones;
  std::vector<std::vector<RooRealVar*> > fStaticTgtsClones;
  std::vector<std::vector<RooAbsReal*> > fStaticPdfsClones;
  std::vector<std::vector<RooRealVar*> > fFullParmsClones;
  std::vector<std::vector<RooRealVar*> > fExtVarsClones;
  std::vector<RooArgSet> fParmSetClones;
  RooArgList fClones;
  
  int fNThreads;
  
  std::vector<RooAbsReal*> fLLRTargets;
  
  std::vector<std::vector<RooAbsReal*> > fDerivatives;
  std::vector<std::vector<std::vector<RooAbsReal*> > > f2Derivatives;
  
  double fLambdaVal;
  TMatrixDSym fHessian;
  //TVectorD fGradient;
  
  std::vector<TMatrixD> fHessians;
  std::vector<TVectorD> fGradients;

  RooAbsReal *fExternal;
  RooAbsReal *fN0;
  RooRealVar *fR;
  double fN0Obs;
  //double fNLLVal;
  RooRealVar *fLambda;  
  
  RooRealVar *fConstraintVal;
  RooRealVar *fConstraintCoeff;
  
  double fRMin;
  double fRHigh;
  double fRLow;
  
  TGraph *fDrvGraph;
  TGraph *fDrvGraphSmooth;
  TGraph *fGraphDelta;
  
  RooArgSet fGarbageCollection;
  
  double fProcNorm;
  
  std::vector<std::vector<int> > fOuterIndices;
  std::vector<std::set<std::pair<int,int> > > fIndices;
  
  
  //std::vector<TTree*>       fTrees;
  std::vector<double>       fTreeWeights;
  std::string               fTrainingCut;
  std::vector<std::string>  fInputVars;  
  std::string               fTargetVar;
  int                       fMinEvents;
  std::vector<double>       fMinWeights;
  double                    fMinWeightTotal;
  double                    fShrinkage;
  int                       fNTrees;
  const int                 fNQuantiles;
  const unsigned int        fNBinsMax;
  float                     fTransitionQuantile;
  double                    fMinCutSignificance;
  double                    fMinCutSignificanceMulti;
  double                    fMaxNSpurious;
  bool                      fDoInitialFit;
  
  double                    fSumWTimesNVars;
  
  
  int                       fMaxDepth;
  int                       fMaxNodes;

  std::vector<double>       fSigmaConsts;
  
  int                      fNTargets;
  //double                   fNSMCObs;
  double                   fNLLVal;
  double                   fdLdR;
  
  int                      fPrescaleInit;
  
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
  int *_bestbins; 
  
  
  double **_ws; 
  double **_ws2; 
  double ***_wscls; 
  int **_ns; 
  int **_nsd;   
  double **_tgts; 
  double **_tgt2s; 
  double **_sumws; 
  double **_sumws2;       
  double ***_sumwscls; 
  int **_sumns; 
  int **_sumnsd;   
  double **_sumtgts; 
  double **_sumtgt2s; 
  float **_varvals; 
  float **_bsepgains; 
  float **_bsepgainsigs;       
  
  int **_quants; 
  int **_binquants; 
  
  int *_clss; 
  double *_tgtvals;
  double *_tgt2vals;
  double *_weightvals;
  
  float **fQuantileMaps;     
  
  std::vector<int> sparserows;
  std::vector<int> sparsecols;
  std::vector<double> sparsedata;
  
  std::vector<double> fStepSizes;
  
//   float *a;
//   float *b;
//   float *c;
//   
//   static void vtest(const float **ra, const float **rb, float **rc);
   
//   float *_sepgains; //!
//   float *_sepgainsigs; //!
//   float *_cutvals;  //!
//   int *_nlefts; //!
//   int *_nrights; //!
//   float *_sumwlefts; //!
//   float *_sumwrights;   //!
//   float *_sumtgtlefts; //!
//   float *_sumtgtrights; //!
//   float *_leftvars; //!
//   float *_rightvars;    //!   
//   float *_fullvars; //!
//   int   *_bestbins; //!
//   
//   
//   float **_ws; //!
//   float **_ws2; //!
//   int **_ns; //!
//   int **_nsd; //!  
//   float **_tgts; //!
//   float **_tgt2s; //!
//   float **_sumws; //!
//   float **_sumws2;    //!   
//   int **_sumns; //!
//   int **_sumnsd; //!  
//   float **_sumtgts; //!
//   float **_sumtgt2s; //!
//   float **_varvals; //!
//   float **_bsepgains; //!
//   float **_bsepgainsigs;  //!     
//   
//   int **_quants; //!
//   int **_bins; //!
//   
//   float **fQuantileMaps;   //!
  
  std::vector<HybridGBREvent*> fEvts;
  
  
private:
  //ClassDef(RooHybridBDTAutoPdf,1) // Exponential PDF
};


#endif
