/*****************************************************************************
 * Project: RooFit                                                           *
 * Package: RooFitModels                                                     *
 * @(#)root/roofit:$Id: HGGRooPdfs.cc,v 1.1 2012/02/10 15:10:48 gpetrucc Exp $
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
 
//////////////////////////////////////////////////////////////////////////////
//
// BEGIN_HTML
// Power function p.d.f
// END_HTML
//
 
#include "RooFit.h"

#include "Riostream.h"
#include "Riostream.h"
#include <math.h>

#include "GBRLikelihood/RooHybridBDTAutoPdf.h"
#include "RooRealVar.h"
#include "RooAbsData.h"
#include "RooUniform.h"
#include "RooDataSet.h"
#include "RooDataHist.h"
#include "RooHistPdf.h"
#include "RooGaussian.h"
#include "RooExponential.h"
#include "RooConstVar.h"
#include "RooDerivative.h"
#include "TFitter.h"
#include "TCanvas.h"
#include "RooPlot.h"
#include "TRandom.h"
#include "TMath.h"
#include <Math/QuantFuncMathCore.h>
#include <Math/ProbFunc.h>
#include "GBRLikelihood/HybridGBRForestD.h"
#include "GBRLikelihood/HybridGBRForestFlex.h"
#include "TMatrixDSym.h"
#include "TVectorD.h"
#include "TDecompChol.h"
#include "TDecompBK.h"
#include "TDecompLU.h"
#include "TDecompSparse.h"
#include "TTree.h"
#include "TTreeFormula.h"
#include "RooCategory.h"
#include "RooSimultaneous.h"
#include "RooMinimizer.h"
#include "RooAddition.h"
#include "RooProduct.h"
#include "TNtuple.h"
#include "TDecompQRH.h"
#include "TDecompSVD.h"
#include "RooCFunction1Binding.h"
#include "TH1D.h"
#include "omp.h"
#include <malloc.h>
#include "GBRLikelihood/GBRArrayUtils.h"
#include "GBRLikelihood/GBRMath.h"
#include "TRandom3.h"
#include "RooFormulaVar.h"

ClassImp(RooTreeConvert)

RooDataSet *RooTreeConvert::CreateDataSet(std::string name, TTree *tree, std::vector<std::string> vars, std::string weight) {
 
  //initialize TTreeFormulas to read variables from TTree
  std::vector<TTreeFormula*> inputforms;
  RooArgList roovars;
  
  for (std::vector<std::string>::const_iterator it = vars.begin(); 
      it != vars.end(); ++it) {
    inputforms.push_back(new TTreeFormula(it->c_str(),it->c_str(),tree));
    RooRealVar *roovar = new RooRealVar(it->c_str(),"",0.);
    roovar->setConstant(false);
    roovars.add(*roovar);
  }  
  
  RooArgSet dvars(roovars);
  
  TTreeFormula cutform(weight.c_str(),weight.c_str(),tree);  
  RooRealVar *weightvar = new RooRealVar(TString::Format("%s_weight",name.c_str()),"",1.);
  roovars.add(*weightvar);
  
  RooDataSet *dset = new RooDataSet(name.c_str(),"",roovars,RooFit::WeightVar(*weightvar));
  
  std::vector<double> minvals(vars.size(),std::numeric_limits<double>::max());
  std::vector<double> maxvals(vars.size(),-std::numeric_limits<double>::max());
  
  for (Long64_t iev=0; iev<tree->GetEntries(); ++iev) {
    //    if (iev%100000==0) // printf("%i\n",int(iev));
    tree->LoadTree(iev);
    
    cutform.GetNdata();
    double weight = cutform.EvalInstance();
    
    if (weight==0.) continue; //skip events with 0 weight
    
    for (unsigned int ivar=0; ivar<vars.size(); ++ivar) {
      inputforms[ivar]->GetNdata();
      double val = inputforms[ivar]->EvalInstance();
      //      std::cout << "RCLSA " << inputforms[ivar]->GetExpFormula().Data() << " " << val << std::endl;
      if (val<minvals[ivar]) minvals[ivar] = val;
      if (val>maxvals[ivar]) maxvals[ivar] = val;
      static_cast<RooRealVar*>(roovars.at(ivar))->setVal(val);
    }
    dset->add(dvars,weight);

  }

  for (std::vector<TTreeFormula*>::const_iterator it = inputforms.begin(); 
      it != inputforms.end(); ++it) {
    delete *it;
  }  
  
//   for (unsigned int ivar=0; ivar<vars.size(); ++ivar) {
//     // printf("ivar = %i, min = %5f, max = %5f\n",ivar,minvals[ivar],maxvals[ivar]);
//     static_cast<RooRealVar*>(roovars.at(ivar))->setRange(minvals[ivar],maxvals[ivar]);
//   }
  
  
  return dset;
  
}

RooDataSet *RooTreeConvert::CreateDataSet(std::string name, TTree *tree, RooArgList &vars, RooRealVar &weight, bool limitvals) {
 
  //initialize TTreeFormulas to read variables from TTree
  std::vector<TTreeFormula*> inputforms;
  
  for (int ivar=0; ivar<vars.getSize(); ++ivar) {
    RooRealVar *var = static_cast<RooRealVar*>(vars.at(ivar));
    inputforms.push_back(new TTreeFormula(var->GetTitle(),var->GetTitle(),tree));
    var->setConstant(false);
  }  
  
  RooArgList roovars(vars);
  
  TTreeFormula cutform(weight.GetTitle(),weight.GetTitle(),tree);  
  weight.setConstant(false);

  roovars.add(weight);
  
  RooDataSet *dset = new RooDataSet(name.c_str(),"",roovars,RooFit::WeightVar(weight));
  
  const int nvars = vars.getSize();
  std::vector<RooRealVar*> varsv(vars.getSize());
  for (int ivar=0; ivar<nvars; ++ivar) {
    varsv[ivar] = static_cast<RooRealVar*>(vars.at(ivar));
  }
  
  RooArgSet varss(vars);
  
  
  int currenttree = -1;
  for (Long64_t iev=0; iev<tree->GetEntries(); ++iev) {
    tree->LoadTree(iev);
    int thistree = tree->GetTreeNumber();
    bool newtree = currenttree!=thistree;
    currenttree = thistree;    
    
    if (newtree) {
      cutform.Notify();
      for (int ivar=0; ivar<nvars; ++ivar) {
        inputforms[ivar]->Notify();      
      }
    }

    // RCLSA: I am forcing this to be 1... 
    double weight = cutform.EvalInstance(); 
    if (weight==0.) continue; //skip entries with 0 weight
    
    bool valid = true;
    for (int ivar=0; ivar<nvars; ++ivar) {
      inputforms[ivar]->GetNdata();
      RooRealVar *var = varsv[ivar];
      double val = inputforms[ivar]->EvalInstance();
      if (val<var->getMin() || val>var->getMax()) {
	valid = false;
      }
      var->setVal(val);
    }
    if (!valid && limitvals) continue;
    
    //dset->add(vars,weight);
    dset->addFast(varss,weight);

  }

  for (std::vector<TTreeFormula*>::const_iterator it = inputforms.begin(); 
      it != inputforms.end(); ++it) {
    delete *it;
  }  
  
  return dset;
  
}


ClassImp(RooNormPdf)

RooNormPdf::RooNormPdf(const char *name, const char *title, RooAbsPdf &pdf, const RooArgSet &forcednormset) :
  RooAbsReal(name,title),
  _pdf("pdf","",this,pdf),
  _forcednormset("forcednormset","",this)
{

  _forcednormset.add(forcednormset);
  
}
  
  
RooNormPdf::RooNormPdf(const RooNormPdf& other, const char* name) :
  RooAbsReal(other,name),
  _pdf("pdf",this,other._pdf),
  _forcednormset("pdf",this,other._forcednormset)
{
  
}

Double_t RooNormPdf::evaluate() const
{
 
  return static_cast<const RooAbsPdf*>(_pdf.absArg())->getNorm(_forcednormset);
  
}


ClassImp(RooRealConstraint)

RooRealConstraint::RooRealConstraint(const char *name, const char *title, RooAbsReal &real, double low, double high) :
  RooAbsReal(name,title),
  _real("real","",this,real),
  _low(low),
  _high(high),
  _scale(0.5*(_high-_low)),
  _offset(_low + 0.5*(_high-_low))
{

  RooGBRTarget *tgt = dynamic_cast<RooGBRTarget*>(&real);
  RooGBRTargetFlex *tgtflex = dynamic_cast<RooGBRTargetFlex*>(&real);
  RooRealVar *var = 0;
  
  
  if (tgt) var = tgt->Var();
  else if (tgtflex) var = tgtflex->Var();
  else var = dynamic_cast<RooRealVar*>(&real);
  
  
  if (var) {
    double oldval = var->getVal();
    double newval = asin(2.0*(oldval-_low)/(_high-_low)-1.0);
    //double newval = atanh( (oldval - _offset)/_scale );
    //double newval = -log(_scale/(oldval-_low) - 1.0);
    //double newval = tan( (oldval - _offset)/_scale );
    var->setVal(newval);
    
    // printf("oldval = %f, newval = %5f, evaluate = %5f\n",oldval,newval,evaluate());
  }
  
}
  
  
RooRealConstraint::RooRealConstraint(const RooRealConstraint& other, const char* name) :
  RooAbsReal(other,name),
  _real("real",this,other._real),
  _low(other._low),
  _high(other._high),
  _scale(other._scale),
  _offset(other._offset)
{
  
}

Double_t RooRealConstraint::evaluate() const
{
 
//   double hprd = (_high -_low);
//   double upbound = _low + hprd;
//   double lowbound = _low - hprd;
//   
//   double val = _real.arg().getVal();
//   
//   while (val>upbound) {
//     val -= 2.0*hprd;
//   }
//   while (val<lowbound) {
//     val += 2.0*hprd;
//   }
//   
//   return _low + std::abs(val-_low);
  
  //return (_low + std::abs((_real.arg().getVal() % (2.0*(_high-_low))) - _high));
  //return std::max(std::min(_real.arg().getVal(),_high),_low);
  //return std::max(std::min(_real.arg().getVal(),_high),_low);
  //return _low + 0.5*(_high-_low)*(sin(_real)+1.0);
  
  //return _offset + _scale*vdt::fast_sinf(_real);
  return _offset + _scale*vdt::fast_sin(_real);
  //return _offset + _scale*sin(_real);
  
  //return _offset + _scale*tanh(_real);
  //return _low + _scale/(1.0+exp(-_real));
  //return _offset + _scale*atan(_real);
  
}


ClassImp(RooPowerLaw)

RooPowerLaw::RooPowerLaw(const char *name, const char *title, RooAbsReal &x, RooAbsReal &p) :
  RooAbsPdf(name,title),
  _x("x","",this,x),
  _p("p","",this,p)
{

/*  double xmin = 100.;
  double xmax = 180.;
  double omp = 1e-9;
  double testint = ((gbrmath::fast_pow(xmax,omp)-gbrmath::fast_pow(xmin,omp))*vdt::fast_inv(omp));
  double testintd = ((pow(xmax,omp)-pow(xmin,omp))/omp);

  double realint1 = vdt::fast_log(xmax)-vdt::fast_log(xmin);
  double realint2 = vdt::fast_log(xmax)-vdt::fast_log(xmin) + 0.5*omp*(vdt::fast_log(xmax)*vdt::fast_log(xmax)-vdt::fast_log(xmin)*vdt::fast_log(xmin));
  
  // printf("testint = %10e, testintd = %10e, realint1 = %10e, realint2 = %10e\n",testint,testintd,realint1,realint2);
  return; */ 
  
}
  
  
RooPowerLaw::RooPowerLaw(const RooPowerLaw& other, const char* name) :
  RooAbsPdf(other,name),
  _x("x",this,other._x),
  _p("p",this,other._p)
{
  
}


Double_t RooPowerLaw::evaluate() const
{
 
  //const RooArgSet* nset = _normSet ;  
  
  return gbrmath::fast_pow(_x,_p);
    
}

Int_t RooPowerLaw::getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analVars, const char* rangeName) const
{
   
  if (matchArgs(allVars,analVars,_x)) return 1 ;
  //if (matchArgs(allVars,analVars,_p)) return 2 ;
  return 0;
  
}

Double_t RooPowerLaw::analyticalIntegral(Int_t code, const char* rangeName) const
{
  
  assert(code==1);
  
  double xmax = _x.max(rangeName);
  double xmin = _x.min(rangeName);
  
  double omp = 1.0 + _p;
  
  if (std::abs(omp)>1e-9) {
    return  ((gbrmath::fast_pow(xmax,omp)-gbrmath::fast_pow(xmin,omp))*vdt::fast_inv(omp));
  }
  else {
    double logxmax = vdt::fast_log(xmax);
    double logxmin = vdt::fast_log(xmin);
    //series expansion for normalization integral around omp=0, precise to O(omp^2)
    return ( logxmax - logxmin + 0.5*omp*(logxmax*logxmax - logxmin*logxmin) );
  }
}

ClassImp(RooCondAddPdf)

RooCondAddPdf::RooCondAddPdf(const char *name, const char *title, const RooArgList &pdfs, const RooArgList &coeffs) :
  RooAbsPdf(name,title),
  _pdfs("pdfs","",this),
  _coeffs("coeffs","",this),
  _selfnorm(coeffs.getSize()==pdfs.getSize() ? false : true)
{
  _pdfs.add(pdfs);
  _coeffs.add(coeffs);
}
  
  
RooCondAddPdf::RooCondAddPdf(const RooCondAddPdf& other, const char* name) :
  RooAbsPdf(other,name),
  _pdfs("pdfs",this,other._pdfs),
  _coeffs("coeffs",this,other._coeffs),
  _selfnorm(other._selfnorm)
{
  
}

Double_t RooCondAddPdf::evaluate() const
{
 
  //const RooArgSet* nset = _normSet ;  
  
  
  double sumcoeff = 0.;
  double val = 0.;
  //double finalcoeff = 1.0;
  for (int i=0; i<_coeffs.getSize(); ++i) {
    double coval = static_cast<RooAbsReal*>(_coeffs.at(i))->getVal();
    val += coval*static_cast<RooAbsPdf*>(_pdfs.at(i))->getValV(_normSet);
    sumcoeff += coval;
    //finalcoeff -= coval;
  }
  
  if (_selfnorm) {
    val += (1.0-sumcoeff)*static_cast<RooAbsPdf*>(_pdfs.at(_pdfs.getSize()-1))->getValV(_normSet);
  }
  else {
    val /= sumcoeff;
  }
  
  return val;
  
}

Int_t RooCondAddPdf::getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analVars, const char* rangeName) const
{
 
  int code = static_cast<RooAbsReal*>(_pdfs.at(0))->getAnalyticalIntegral(allVars,analVars,rangeName);
  for (int ipdf=1; ipdf<_pdfs.getSize(); ++ipdf) {
    int pdfcode = static_cast<RooAbsPdf*>(_pdfs.at(ipdf))->getAnalyticalIntegral(allVars,analVars,rangeName);
    if (pdfcode!=code) code = 0;
  }
  
  return code;
  
}

Double_t RooCondAddPdf::analyticalIntegral(Int_t code, const char* rangeName) const
{ 
  return 1.0;
}

RooAbsReal *RooCondAddPdf::createCDF(const RooArgSet &iset, const RooArgSet &nset) {
 
  RooArgList cdfcomps;
  
  RooAbsArg *lastcoeff = 0;
  RooAddition *sumcoeff = new RooAddition(TString::Format("%s_sumcoeff",GetName()),"",_coeffs);
  
  if (_selfnorm) {
    lastcoeff = new RooFormulaVar(TString::Format("%s_lastcoeff",GetName()),"","1.0 - @0",*sumcoeff);
  }
  else {
    lastcoeff = _coeffs.at(_pdfs.getSize()-1);
  }
  
  RooArgSet fullset;
  fullset.add(iset);
  fullset.add(nset);
  
  for (int ipdf=0; ipdf<_pdfs.getSize(); ++ipdf) {
    RooAbsReal *singlecdf = static_cast<RooAbsPdf*>(_pdfs.at(ipdf))->createRunningIntegral(iset,nset);
    //const RooAbsReal *singlenorm = static_cast<RooAbsPdf*>(_pdfs.at(ipdf))->getNormIntegral(fullset);
    RooAbsArg *coeff = ipdf==(_pdfs.getSize()-1) ? lastcoeff : _coeffs.at(ipdf);
    RooProduct *singlecdfnorm = new RooProduct(TString::Format("%s_singlecdfnorm_%i",GetName(),ipdf),"",RooArgList(*coeff,*singlecdf));
    //RooFormulaVar *singlecdfnorm = new RooFormulaVar(TString::Format("%s_singlecdfnorm_%i",GetName(),ipdf),"","@0*@1/@2",RooArgList(*coeff,*singlecdf,*singlenorm));
    cdfcomps.add(*singlecdfnorm);
  }
  
  RooAbsReal *cdf = new RooAddition(TString::Format("%s_cdf",GetName()),"",cdfcomps);
  
  RooAbsReal *cdfnorm = 0;
  if (_selfnorm) {
    cdfnorm = cdf;
  }
  else {
    RooFormulaVar *normfactor = new RooFormulaVar(TString::Format("%s_normfactor",GetName()),"","1.0/@0",*sumcoeff);
    cdfnorm = new RooProduct(TString::Format("%s_cdfnorm",GetName()),"",RooArgList(*normfactor,*cdf));
  }
    
  return cdfnorm;
  
}


ClassImp(RooCondRatioPdf)

RooCondRatioPdf::RooCondRatioPdf(const char *name, const char *title, RooAbsReal &ratio, RooAbsReal &pdfden) :
  RooAbsPdf(name,title),
  _ratio("ratio","",this,ratio),
  _pdfden("pdfden","",this,pdfden)
{

}
  
  
RooCondRatioPdf::RooCondRatioPdf(const RooCondRatioPdf& other, const char* name) :
  RooAbsPdf(other,name),
  _ratio("ratio",this,other._ratio),
  _pdfden("pdfden",this,other._pdfden)
{
  
}

Double_t RooCondRatioPdf::evaluate() const
{
 
  return _ratio.arg().getVal(_normSet)*_pdfden.arg().getVal(_normSet);
  
}

// Int_t RooCondRatioPdf::getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analVars, const char* rangeName) const
// {
//  
//   int code = _pdfden.arg().getAnalyticalIntegral(allVars,analVars,rangeName);
// 
//   return code;
//   
// }

// Double_t RooCondRatioPdf::analyticalIntegral(Int_t code, const char* rangeName) const
// {
//   
//   return 1.0;
//  
// //   double finalcoeff = 1.0;
// //   double integral=0.;
// //   for (int i=0; i<_coeffs.getSize(); ++i) {
// //     double coval = static_cast<RooAbsReal*>(_coeffs.at(i))->getVal();
// //     integral += coval*static_cast<RooAbsReal*>(_pdfs.at(i))->analyticalIntegral(code,rangeName);
// //     finalcoeff -= coval;
// //   }
// //   
// //   integral += finalcoeff*static_cast<RooAbsReal*>(_pdfs.at(_pdfs.getSize()-1))->analyticalIntegral(code,rangeName);
// //   
// //   return integral;
//   
// }


ClassImp(RooPdfAddReal)

RooPdfAddReal::RooPdfAddReal(const char *name, const char *title, const RooArgList &pdfs, const RooArgList &coeffs) :
  RooAbsReal(name,title),
  _pdfs("pdfs","",this),
  _coeffs("coeffs","",this)
{
  _pdfs.add(pdfs);
  _coeffs.add(coeffs);
}
  
  
RooPdfAddReal::RooPdfAddReal(const RooPdfAddReal& other, const char* name) :
  RooAbsReal(other,name),
  _pdfs("pdfs",this,other._pdfs),
  _coeffs("coeffs",this,other._coeffs)
{
  
}

Double_t RooPdfAddReal::evaluate() const
{
 
  const RooArgSet *nset = _pdfs.nset();  
  
  //if (nset) nset->Print("V");
  
  double val = 0.;
  for (int i=0; i<_coeffs.getSize(); ++i) {
    double coval = static_cast<RooAbsReal*>(_coeffs.at(i))->getVal();
    val += coval*static_cast<RooAbsPdf*>(_pdfs.at(i))->getValV(nset);
  }
    
  //// printf("val = %5f\n",val);
    
  return val;
  
}

Int_t RooPdfAddReal::getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analVars, const char* rangeName) const
{
 
  int code = static_cast<RooAbsReal*>(_pdfs.at(0))->getAnalyticalIntegral(allVars,analVars,rangeName);
  for (int ipdf=1; ipdf<_pdfs.getSize(); ++ipdf) {
    int pdfcode = static_cast<RooAbsPdf*>(_pdfs.at(ipdf))->getAnalyticalIntegral(allVars,analVars,rangeName);
    if (pdfcode!=code) code = 0;
  }
  
  //// printf("RooPdfAddReal code = %i\n",code);
  //return 1;
  return code;
  
}

Double_t RooPdfAddReal::analyticalIntegral(Int_t code, const char* rangeName) const
{
  
  //return 1.0;
 
  double integral=0.;
  for (int i=0; i<_coeffs.getSize(); ++i) {
    double coval = static_cast<RooAbsReal*>(_coeffs.at(i))->getVal();
    integral += coval*static_cast<RooAbsPdf*>(_pdfs.at(i))->analyticalIntegral(code,rangeName);
  }
    
  return integral;
  
}


ClassImp(RooGBRFunctionFlex)

//_____________________________________________________________________________
RooGBRFunctionFlex::RooGBRFunctionFlex(const char *name, const char *title) :
  RooAbsReal(name,title),
  _forest(new HybridGBRForestFlex())
{
  
}

//_____________________________________________________________________________
RooGBRFunctionFlex::RooGBRFunctionFlex(const RooGBRFunctionFlex& other, const char* name) :
  RooAbsReal(other,name), 
  _forest(new HybridGBRForestFlex(*other._forest))
{
  
}
  
//_____________________________________________________________________________  
RooGBRFunctionFlex::~RooGBRFunctionFlex()
{
  if (_forest) delete _forest;
}
  
//_____________________________________________________________________________  
void RooGBRFunctionFlex::SetForest(HybridGBRForestFlex *forest) {
 
  if (_forest) delete _forest;
  _forest = forest;
  setValueDirty();
  
}

//_____________________________________________________________________________  
void RooGBRFunctionFlex::ResetForest() {
 
  SetForest(new HybridGBRForestFlex());
  
}

ClassImp(RooGBRTargetFlex)

//_____________________________________________________________________________
RooGBRTargetFlex::RooGBRTargetFlex(const char *name, const char *title, RooGBRFunctionFlex &func, RooRealVar &var, const RooArgList &funcvars) :
  RooAbsReal(name,title),
  _func("func","",this,func,true,false),
  _var("var","",this,var),
  _usefunc(false),
  _funcvars("funcvars","",this),
  _eval(funcvars.getSize())
{
  _funcvars.add(funcvars);
  
  Var()->setConstant(_usefunc);
  
}

//_____________________________________________________________________________
RooGBRTargetFlex::RooGBRTargetFlex(const RooGBRTargetFlex& other, const char* name) :
  RooAbsReal(other,name), 
  _func("func",this,other._func),  
  _itgt(other._itgt),
  _var("var",this,other._var),
  _usefunc(other._usefunc),
  _funcvars("funcvars",this,other._funcvars),
  _eval(other._eval)
{
    
}

//_____________________________________________________________________________
void RooGBRTargetFlex::SetUseFunc(bool b)
{
    
  if (_usefunc==b) return;
  
  _usefunc = b;
  Var()->setConstant(_usefunc);

  setValueDirty();
  setShapeDirty();
  
}

//_____________________________________________________________________________
void RooGBRTargetFlex::ClearFuncServers() {
  
  //clear function and input variables from value servers,
  //not needed during training phase
  
  removeServer(*_func.absArg());
  for (int ivar=0; ivar<_funcvars.getSize(); ++ivar) {
    removeServer(*_funcvars.at(ivar));
  }  
  
}

//_____________________________________________________________________________
double RooGBRTargetFlex::EvalFunc() const
{
  
  RooFIter iter = _funcvars.fwdIterator();
  for (int ivar=0; ivar<_funcvars.getSize(); ++ivar) {
    _eval[ivar] = static_cast<RooAbsReal*>(iter.next())->getVal();
    //// printf("ivar = %i, var = %5f\n",ivar,_eval[ivar]);
  }  
  return static_cast<RooGBRFunctionFlex*>(_func.absArg())->Forest()->GetResponse(&_eval[0]);
  
}


ClassImp(RooGBRFunction)

//_____________________________________________________________________________
RooGBRFunction::RooGBRFunction(const char *name, const char *title, const RooArgList &vars, int ntargets) :
  RooAbsReal(name,title),
  _vars("vars","",this),
  _forest(new HybridGBRForestD(ntargets)),
  _eval(vars.getSize())
{
  _vars.add(vars);
}

//_____________________________________________________________________________
RooGBRFunction::RooGBRFunction(const RooGBRFunction& other, const char* name) :
  RooAbsReal(other,name), 
  _vars("vars",this,other._vars),  
  _forest(new HybridGBRForestD(*other._forest)),
  _eval(other._eval)
{
  
}
  
//_____________________________________________________________________________  
RooGBRFunction::~RooGBRFunction()
{
  if (_forest) delete _forest;
}
  
//_____________________________________________________________________________  
double RooGBRFunction::GetResponse(int itgt) const {
  
  //// printf("RooGBRFunction::GetResponse(%i)\n",itgt);
//   if (isValueDirtyAndClear()) {
//   //if (1) {
//     //// printf("recomputing bdt response\n");
//     for (int ivar=0; ivar<_vars.getSize(); ++ivar) {
//       _eval[ivar] = static_cast<RooAbsReal*>(_vars.at(ivar))->getVal();
//       //// printf("ivar = %i, var = %5f\n",ivar,_eval[ivar]);
//     }
//     _forest->GetResponse(&_eval[0]);
//   }
//   
//   //// printf("response %i = %5f\n",itgt,_forest->GetResponse(itgt));
//   return _forest->GetResponse(itgt);
  
  
  for (int ivar=0; ivar<_vars.getSize(); ++ivar) {
    _eval[ivar] = static_cast<RooAbsReal*>(_vars.at(ivar))->getVal();
    //// printf("ivar = %i, var = %5f\n",ivar,_eval[ivar]);
  }  
  return _forest->GetResponse(&_eval[0],itgt);

  
}

//_____________________________________________________________________________  
void RooGBRFunction::SetForest(HybridGBRForestD *forest) {
 
  if (_forest) delete _forest;
  _forest = forest;
  setValueDirty();
  
}

ClassImp(RooGBRTarget)

//_____________________________________________________________________________
RooGBRTarget::RooGBRTarget(const char *name, const char *title, RooGBRFunction &func, int itgt, RooRealVar &var) :
  RooAbsReal(name,title),
  _func("func","",this,func,true,false),
  _itgt(itgt),
  _var("var","",this,var),
  _usefunc(false)
{
  Var()->setConstant(_usefunc);
  //unRegisterProxy(_func);
  
}

//_____________________________________________________________________________
RooGBRTarget::RooGBRTarget(const RooGBRTarget& other, const char* name) :
  RooAbsReal(other,name), 
  _func("func",this,other._func),  
  _itgt(other._itgt),
  _var("var",this,other._var),
  _usefunc(other._usefunc)  
{
  
//   if (_usefunc) {
//     unRegisterProxy(_var);
//   }
//   else {
//     unRegisterProxy(_func);
//   }
//   
//   setValueDirty();
//   setShapeDirty();
  
}

//_____________________________________________________________________________
void RooGBRTarget::SetUseFunc(bool b)
{
  
  if (_usefunc==b) return;
  
  _usefunc = b;
  Var()->setConstant(_usefunc);
  
//   if (_usefunc) {
//     registerProxy(_func);
//     unRegisterProxy(_var);
//   }
//   else {
//     unRegisterProxy(_func);
//     registerProxy(_var);
//   }  
  
  setValueDirty();
  setShapeDirty();
  
}

//ClassImp(RooHybridBDTAutoPdf) 

RooHybridBDTAutoPdf *gHybridBDTAutoPointer;


//_____________________________________________________________________________
RooHybridBDTAutoPdf::RooHybridBDTAutoPdf(const char *name, const char *title, const RooArgList &tgtvars, RooAbsReal &n0, RooRealVar &r, const std::vector<RooAbsData*> &data, const std::vector<RooAbsReal*> &pdfs) :
  TNamed(name,title),
  fTgtVars(tgtvars),
  //fCondVars(fFunc->Vars()),
  fResTree(0),
  fPdfs(pdfs),
  fData(data),
   fNThreads(std::max(1,omp_get_max_threads())),
  // fNThreads(1),  
  fLambdaVal(1.0), 
  fExternal(&n0),  
  fR(&r),  
  fN0Obs(data.front()->sumEntries()),  
  fLambda(0),  
  fDrvGraph(0),
  fDrvGraphSmooth(0),
  fGraphDelta(0),
  //fLambda(new RooRealVar("lambda","",0.)),
  fMinEvents(-99),  
  //fMinWeights(std::vector<double>(data.size(),1000.)),
  fMinWeightTotal(-99.),
  fShrinkage(0.5),
  fNTrees(20),
  fNQuantiles(std::numeric_limits<unsigned short>::max()+1),
  //fNQuantiles(128),
  fNBinsMax(128),
  //fNBinsMax(fNQuantiles),
  fTransitionQuantile(0.7),
  fMinCutSignificance(-99.),
  fMinCutSignificanceMulti(-99.),
  fMaxNSpurious(-99.),
  fDoInitialFit(true),
  fSumWTimesNVars(0.),
  fMaxDepth(-1),
  fMaxNodes(-1),
  fNTargets(tgtvars.getSize()),  
  fPrescaleInit(-1),
  _sepgains(0),
  _ws(0)  
{
  
  omp_set_num_threads(fNThreads);
  
  //// printf("fN0 = %5f\n",fN0->getVal());
  
  //create constraint term for profile likelihood gradient scan and combine with external likelihood
  fConstraintVal = new RooRealVar(TString::Format("%s_constraintval",GetName()),"",fR->getVal());
  fConstraintCoeff = new RooRealVar(TString::Format("%s_constraintcoeff",GetName()),"",0.);
  RooFormulaVar *constraint = new RooFormulaVar(TString::Format("%s_constraint",GetName()),"","@0*pow(@1-@2,2)",RooArgList(*fConstraintCoeff,*fR,*fConstraintVal));
  
  fN0 = new RooAddition(TString::Format("%s_fullexternal",GetName()),"",RooArgList(*fExternal,*constraint));
  
  fGarbageCollection.addOwned(*fConstraintVal);
  fGarbageCollection.addOwned(*fConstraintCoeff);
  fGarbageCollection.addOwned(*constraint);
  fGarbageCollection.addOwned(*fN0);
  
  fTgtCondVars.resize(fTgtVars.getSize(),RooArgList());
  
  //Fill RooRealVars underlying RooGBRTargetFlex objects, as well as conditional variables and underlying RooGBRFunctionFlex's
  //Fill also corresponding maps between targets, functions, and variables
  for (int itgt=0; itgt<fTgtVars.getSize(); ++itgt) {
    RooGBRTargetFlex *target = static_cast<RooGBRTargetFlex*>(fTgtVars.at(itgt));
    target->SetUseFunc(false);
    
    RooRealVar *var = target->Var();
    fStaticTgts.add(*var);
    
    if (!fFuncs.contains(*target->Func())) {
      fFuncs.add(*target->Func());
      fFuncTgts.push_back(RooArgList());
      fFuncCondVars.push_back(RooArgList());
    }
    
    int ifunc = fFuncs.index(target->Func());
    
    fFuncTgts[ifunc].add(*target);
    
    
    const RooArgList &funcvars = target->FuncVars();
    for (int ifuncvar=0; ifuncvar<funcvars.getSize(); ++ifuncvar) {
      RooAbsArg *funcvar = funcvars.at(ifuncvar);
      if ( !fCondVars.contains(*funcvar) ) {
        fCondVars.add(*funcvar);
      }
      if (!fTgtCondVars[itgt].contains(*funcvar)) {
        fTgtCondVars[itgt].add(*funcvar);
      }
      if (!fFuncCondVars[ifunc].contains(*funcvar)) {
        fFuncCondVars[ifunc].add(*funcvar);
      }      
    }
  }
  
  for (unsigned int ipdf=0; ipdf<fPdfs.size(); ++ipdf) {
    fStaticPdfs.add(*fPdfs[ipdf]);
  }    
  
  //fCondVars.Print("V");
  
  // printf("filling observables\n");
  //fill observables for paramtric pdfs  
  RooArgSet *allvarss = fPdfs.front()->getObservables(*data.front());
  {
  RooArgList allvars(*allvarss);
  for (int ivar=0; ivar<allvars.getSize(); ++ivar) {
    if (!fParmVars.contains(*allvars.at(ivar)) && !fCondVars.contains(*allvars.at(ivar))) fParmVars.add(*allvars.at(ivar));
  }
  }
  delete allvarss;
  
  

  

  
//   for (int ipdf=0; ipdf<fStaticPdfs.getSize(); ++ipdf) {  
//     fStaticPdfs.at(ipdf)->recursiveRedirectServers(fStaticTgts);
//   }
  
  //fStaticTgts.Print("V");
  
  //targets corresponding to log likelihood ratios
  fLLRTargets.push_back(0);
  for (int ipdf=1; ipdf<fStaticPdfs.getSize(); ++ipdf) {  
    fLLRTargets.push_back(static_cast<RooAbsReal*>(fStaticTgts.at(ipdf-1)));
  }

  // printf("first loop, count events\n");

  //fStaticPdfs.at(0)->Print("V");
  //fStaticPdfs.at(0)->getParameters(*data.front())->Print("V");
  
  
  int nev = 0;
  for (unsigned int idata=0; idata<data.size(); ++idata) {
    nev += data[idata]->numEntries();
    // printf("idata = %i, sumentries = %5f, numentries = %i\n",idata,data[idata]->sumEntries(),data[idata]->numEntries());
  }
  
  int nvars = fCondVars.getSize() + fParmVars.getSize();
  
  
  //set up computation of first and second derivatives
//   fDerivatives.resize(fFullFuncs.getSize(),std::vector<RooAbsReal*>(nparms));
//   f2Derivatives.resize(fFullFuncs.getSize(),std::vector<std::vector<RooAbsReal*> >(nparms,std::vector<RooAbsReal*>(nparms)));
//   for (int ipdf=0; ipdf<fFullFuncs.getSize(); ++ipdf) {
//     for (int iparm=0; iparm<nparms; ++iparm) {
//       if (fFullFuncs.at(ipdf)->overlaps(*fFullParms.at(iparm))) {
// 	fDerivatives[ipdf][iparm] = static_cast<RooAbsReal*>(fFullFuncs.at(ipdf))->derivative(*static_cast<RooRealVar*>(fFullParms.at(iparm)),fParmVars,1);
//       }
//       else
// 	fDerivatives[ipdf][iparm] = new RooConstVar("constzero","",0.);
//     }
//     for (int iparm=0; iparm<nparms; ++iparm) {
//       for (int jparm=0; jparm<nparms; ++jparm) {
// 	if (fDerivatives[ipdf][iparm]->overlaps(*fFullParms.at(jparm))) {
// 	  f2Derivatives[ipdf][iparm][jparm] = fDerivatives[ipdf][iparm]->derivative(*static_cast<RooRealVar*>(fFullParms.at(jparm)),fParmVars,1);
// 	}
// 	else
// 	  f2Derivatives[ipdf][iparm][jparm] = new RooConstVar("constzero","",0.);
//       }
//     }    
//   }  
  


   
  
//first loop here  
  
  // printf("nev = %i, nvar = %i\n",int(nev),nvars);

  int ncls = fData.size();
  
//  fEvalVector.resize(fCondVars.getSize());
  
 
   
  //initialize arrays (ensure 32 byte alignment for avx vector instructions)

  _sepgains = (float*)memalign(32, nvars*sizeof(float));
  _sepgainsigs = (float*)memalign(32, nvars*sizeof(float));
  _cutvals = (float*)memalign(32, nvars*sizeof(float));
  _nlefts = (int*)memalign(32, nvars*sizeof(int));
  _nrights = (int*)memalign(32, nvars*sizeof(int));
  _sumwlefts = (float*)memalign(32, nvars*sizeof(float));
  _sumwrights = (float*)memalign(32, nvars*sizeof(float));
  _sumtgtlefts = (float*)memalign(32, nvars*sizeof(float));
  _sumtgtrights = (float*)memalign(32, nvars*sizeof(float));
  _leftvars = (float*)memalign(32, nvars*sizeof(float));
  _rightvars = (float*)memalign(32, nvars*sizeof(float));  
  _fullvars = (float*)memalign(32, nvars*sizeof(float));  
  _bestbins = (int*)memalign(32, nvars*sizeof(int));

 
  
  
  
  _ws = new double*[nvars];
  _ws2 = new double*[nvars];
  _wscls = new double**[nvars];
  _ns = new int*[nvars];
  _nsd = new int*[nvars];  
  _tgts = new double*[nvars];  
  _tgt2s = new double*[nvars];  
  _sumws = new double*[nvars];
  _sumws2 = new double*[nvars];
  _sumwscls = new double**[nvars];
  _sumns = new int*[nvars];
  _sumtgts = new double*[nvars];  
  _sumtgt2s = new double*[nvars];
  _varvals = new float*[nvars];    
  _bsepgains = new float*[nvars];
  _bsepgainsigs = new float*[nvars];
  
  _binquants  = new int*[nvars];  
  _quants  = new int*[nvars];  
  
  _clss  = (int*)memalign(32, nev*sizeof(int));
  _tgtvals  = (double*)memalign(32, nev*sizeof(double));
  _tgt2vals  = (double*)memalign(32, nev*sizeof(double));
  _weightvals  = (double*)memalign(32, nev*sizeof(double));
  
  fQuantileMaps = new float*[nvars];  
  
  
  for (int ivar=0; ivar<nvars; ++ivar) {
    _ws[ivar] = (double*)memalign(32, fNBinsMax*sizeof(double));
    _ws2[ivar] = (double*)memalign(32, fNBinsMax*sizeof(double));
    _ns[ivar] = (int*)memalign(32, fNBinsMax*sizeof(int));
    _tgts[ivar] = (double*)memalign(32, fNBinsMax*sizeof(double));
    _tgt2s[ivar] = (double*)memalign(32, fNBinsMax*sizeof(double));  
    _sumws[ivar] = (double*)memalign(32, fNBinsMax*sizeof(double));
    _sumws2[ivar] = (double*)memalign(32, fNBinsMax*sizeof(double));
    _sumns[ivar] = (int*)memalign(32, fNBinsMax*sizeof(int));
    _sumtgts[ivar] = (double*)memalign(32, fNBinsMax*sizeof(double));
    _sumtgt2s[ivar] = (double*)memalign(32, fNBinsMax*sizeof(double));
    _varvals[ivar] = (float*)memalign(32, fNBinsMax*sizeof(float));
    _bsepgains[ivar] = (float*)memalign(32, fNBinsMax*sizeof(float));
    _bsepgainsigs[ivar] = (float*)memalign(32, fNBinsMax*sizeof(float));
    
    _wscls[ivar] = new double*[ncls];
    _sumwscls[ivar] = new double*[ncls];    
    
    _binquants[ivar] = (int*)memalign(32, fNBinsMax*sizeof(int));
    _quants[ivar] = (int*)memalign(32, nev*sizeof(int));
    
    fQuantileMaps[ivar] = (float*)memalign(32, fNQuantiles*sizeof(float));
    

    for (int icls=0; icls<ncls; ++icls) {
      _wscls[ivar][icls] = (double*)memalign(32, fNBinsMax*sizeof(double));
      _sumwscls[ivar][icls] = (double*)memalign(32, fNBinsMax*sizeof(double));
    }
    
  }
  
/*  _sepgains = new float[nvars];
  _sepgainsigs = new float[nvars];
  _cutvals = new float[nvars];
  _nlefts = new int[nvars];
  _nrights = new int[nvars];
  _sumwlefts = new float[nvars];
  _sumwrights = new float[nvars];
  _sumtgtlefts = new float[nvars];
  _sumtgtrights = new float[nvars];
  _leftvars = new float[nvars];
  _rightvars = new float[nvars];
  _fullvars = new float[nvars];
  _bestbins = new int[nvars];

 
  
  
  
  _ws = new double*[nvars];
  _ws2 = new double*[nvars];
  _wscls = new double**[nvars];
  _ns = new int*[nvars];
  _nsd = new int*[nvars];
  _tgts = new double*[nvars];
  _tgt2s = new double*[nvars];
  _sumws = new double*[nvars];
  _sumws2 = new double*[nvars];
  _sumwscls = new double**[nvars];
  _sumns = new int*[nvars];
  _sumtgts = new double*[nvars];
  _sumtgt2s = new double*[nvars];
  _varvals = new float*[nvars];
  _bsepgains = new float*[nvars];
  _bsepgainsigs = new float*[nvars];
  
  _quants = new int*[nvars];
  _bins = new int*[nvars];
  _clss = new int[nvars];
  _tgtvals  = new double[nev];
  _tgt2vals  = new double[nev];
  _weightvals  = new double[nev];
  
  fQuantileMaps = new float*[nvars];
  
  
  for (int ivar=0; ivar<nvars; ++ivar) {
    _ws[ivar] = new double[fNBinsMax];
    _ws2[ivar] = new double[fNBinsMax];
    _ns[ivar] = new int[fNBinsMax];
    _tgts[ivar] = new double[fNBinsMax];
    _tgt2s[ivar] = new double[fNBinsMax];
    _sumws[ivar] = new double[fNBinsMax];
    _sumws2[ivar] = new double[fNBinsMax];
    _sumns[ivar] = new int[fNBinsMax];
    _sumtgts[ivar] = new double[fNBinsMax];
    _sumtgt2s[ivar] = new double[fNBinsMax];
    _varvals[ivar] = new float[fNBinsMax];
    _bsepgains[ivar] = new float[fNBinsMax];
    _bsepgainsigs[ivar] = new float[fNBinsMax];
    
    _wscls[ivar] = new double*[ncls];
    _sumwscls[ivar] = new double*[ncls];
    
    _quants[ivar] = new int[nev];
    _bins[ivar] = new int[fNBinsMax];
    
    
    fQuantileMaps[ivar] = new float[fNQuantiles];
    
    for (int icls=0; icls<ncls; ++icls) {
      _wscls[ivar][icls] = new double[fNBinsMax];
      _sumwscls[ivar][icls] = new double[fNBinsMax];
    }    
    
  }*/   
  
  
  
      


  //TrainForest(fNTrees);
//  TrainForest(fNTrees);
  
  
  fLambdaVal = 1.;
//   int nparms = fExtVars.getSize() + fNTargets*nterm; 
//   fHessian = TMatrixDSym(nparms);   
//   for (int iel=0; iel<nparms; ++iel) {
//     fHessian(iel,iel) = 1.0;
//   }
  

  
//   fExtVars.removeAll();
//   fFullParms.removeAll();
//   fFullFuncs.removeAll();
//   fOuterIndices.clear();
//   fIndices.clear();
  
  //fill list of global parameters
  std::vector<RooAbsReal*> sources;
  for (unsigned int ipdf=0; ipdf<fPdfs.size(); ++ipdf) {
    sources.push_back(fPdfs.at(ipdf));
  }
  sources.push_back(fN0);
  
  
  //// printf("filling parameters\n");
  for (unsigned int isrc=0; isrc<sources.size(); ++isrc) {
    RooArgSet *allparmss = sources[isrc]->getParameters(*fData.front());
    {
    RooArgList allparms(*allparmss);
    for (int ivar=0; ivar<allparms.getSize(); ++ivar) {
      if (!fExtVars.contains(*allparms.at(ivar)) && !fStaticTgts.contains(*allparms.at(ivar)) && !allparms.at(ivar)->getAttribute("Constant")) fExtVars.add(*allparms.at(ivar));
    }
    }
    delete allparmss;
  }
  
  if (dynamic_cast<RooRealVar*>(fN0) && !fExtVars.contains(*fN0) && !fStaticTgts.contains(*fN0) && !static_cast<RooRealVar*>(fN0)->getAttribute("Constant")) fExtVars.add(*fN0);
  
  //// printf("ExtVars:\n");
  //fExtVars.Print("V");  
  
  fFullParms.add(fExtVars);
  fFullParms.add(fStaticTgts);
//  int nparms = fFullParms.getSize();
  
  fFullFuncs.add(fStaticPdfs);
  fFullFuncs.add(*fN0);
  
  //fFullFuncs.Print("V");
  
  fOuterIndices.resize(fFullFuncs.getSize());
  fIndices.resize(fFullFuncs.getSize());
  
  for (int ipdf=0; ipdf<fFullFuncs.getSize(); ++ipdf) {
    for (int iparm=0; iparm<fFullParms.getSize(); ++iparm) {
      if (fFullFuncs.at(ipdf)->overlaps(*fFullParms.at(iparm))) {
	fOuterIndices[ipdf].push_back(iparm);
      }
    }
    for (unsigned int iidx=0; iidx<fOuterIndices[ipdf].size(); ++iidx) {
      for (unsigned int jidx=iidx; jidx<fOuterIndices[ipdf].size(); ++jidx) {
	fIndices[ipdf].insert(std::pair<int,int>(fOuterIndices[ipdf][iidx],fOuterIndices[ipdf][jidx]));
      }
    }
  }
  
  fCondVarsClones.resize(fNThreads);
  fParmVarsClones.resize(fNThreads);
  fStaticTgtsClones.resize(fNThreads);
  fStaticPdfsClones.resize(fNThreads);
  fFullParmsClones.resize(fNThreads);
  fExtVarsClones.resize(fNThreads);
  fParmSetClones.reserve(fNThreads);
  
  RooAddition fullpdfsum("fullpdfsum","",fFullFuncs);
  
  //fullpdfsum.getComponents()->Print("V");
  
  for (int ithread=0; ithread<fNThreads; ++ithread) {
    // printf("ithread = %i\n",ithread);
    RooAbsArg *clone = fullpdfsum.cloneTree();
    fClones.addOwned(*clone);
    
    RooArgSet *clonecomps = clone->getComponents();
    RooArgSet *clonevars = clone->getVariables();
    clonecomps->add(*clonevars);
    //clonecomps->Print("V");
    for (int ivar=0; ivar<fCondVars.getSize(); ++ivar) {
      fCondVarsClones[ithread].push_back(static_cast<RooRealVar*>(clonecomps->find(fCondVars.at(ivar)->GetName())));
    }
    for (int ivar=0; ivar<fParmVars.getSize(); ++ivar) {
      fParmVarsClones[ithread].push_back(static_cast<RooRealVar*>(clonecomps->find(fParmVars.at(ivar)->GetName())));
    }
    for (int ivar=0; ivar<fStaticTgts.getSize(); ++ivar) {
      fStaticTgtsClones[ithread].push_back(static_cast<RooRealVar*>(clonecomps->find(fStaticTgts.at(ivar)->GetName())));
    }
    for (int ivar=0; ivar<fStaticPdfs.getSize(); ++ivar) {
      fStaticPdfsClones[ithread].push_back(static_cast<RooAbsReal*>(clonecomps->find(fStaticPdfs.at(ivar)->GetName())));
    }    
    for (int ivar=0; ivar<fFullParms.getSize(); ++ivar) {
      fFullParmsClones[ithread].push_back(static_cast<RooRealVar*>(clonecomps->find(fFullParms.at(ivar)->GetName())));
    }   
    for (int ivar=0; ivar<fExtVars.getSize(); ++ivar) {
      fExtVarsClones[ithread].push_back(static_cast<RooRealVar*>(clonecomps->find(fExtVars.at(ivar)->GetName())));
    }
    //remove extraneous value servers for target objects
    for (int ivar=0; ivar<fTgtVars.getSize(); ++ivar) {
      RooGBRTargetFlex *target = static_cast<RooGBRTargetFlex*>(clonecomps->find(fTgtVars.at(ivar)->GetName()));
      target->ClearFuncServers();
    }
    
    delete clonecomps;
    delete clonevars;
    
    RooArgSet parmsetclone;
    for (unsigned int ivar=0; ivar<fParmVarsClones[ithread].size(); ++ivar) {
      parmsetclone.add(*fParmVarsClones[ithread][ivar]);
    }
    fParmSetClones.push_back(parmsetclone);
  }
  
//   fCondVarsClones[0].Print("V");
//   fParmVarsClones[0].Print("V");
//   fStaticTgtsClones[0].Print("V");
//   fStaticPdfsClones[0].Print("V");
//   fFullParmsClones[0].Print("V");
  
  fEvts.reserve(nev);
  
  double sumw = 0.;
  double sumabsw = 0.;
  
  // printf("second loop, fill events in memory\n");
  //loop over trees to fill arrays and event vector
   
  
  //second loop here
  for (unsigned int idata=0; idata<data.size(); ++idata) {
    std::vector<RooRealVar*> dcondvars(fCondVars.getSize());
    std::vector<RooRealVar*> dparmvars(fParmVars.getSize());

    const RooArgSet *dset = data[idata]->get();
    
    for (int ivar=0; ivar<fCondVars.getSize(); ++ivar) {
      dcondvars[ivar] = static_cast<RooRealVar*>(dset->find(fCondVars.at(ivar)->GetName()));
    }
      
    for (int ivar=0; ivar<fParmVars.getSize(); ++ivar) {
      dparmvars[ivar] = static_cast<RooRealVar*>(dset->find(fParmVars.at(ivar)->GetName()));
    }
      
    for (int iev=0; iev<data[idata]->numEntries(); ++iev) {
      data[idata]->get(iev);
      fEvts.push_back(new HybridGBREvent(nvars,fNTargets,fFullParms.getSize()));
      HybridGBREvent *evt = fEvts.back();
      evt->SetWeight(data[idata]->weight());
      evt->SetClass(idata);

      sumw += evt->Weight();
      sumabsw += std::abs(evt->Weight());
      
      for (unsigned int ivar=0; ivar<dcondvars.size(); ++ivar) {
	evt->SetVar(ivar,dcondvars[ivar]->getVal());
      }
      for (unsigned int ivar=0; ivar<dparmvars.size(); ++ivar) {
	evt->SetVar(dcondvars.size() + ivar, dparmvars[ivar]->getVal());
      }    
    }
  }
  
  fSumWTimesNVars = sumw*fCondVars.getSize();
  
  // printf("filled data\n");
    
  //int nevr = datau->numEntries();
  
  
  //map of input variable quantiles to values
  //fQuantileMaps.resize(nvars, std::vector<float>(fNQuantiles));
  
  //BuildQuantiles(nvars, sumw);
  BuildQuantiles(nvars, sumabsw);  
  
  
  
}

//_____________________________________________________________________________
RooHybridBDTAutoPdf::~RooHybridBDTAutoPdf() {
  
  int nvars = fCondVars.getSize() + fParmVars.getSize();
  int ncls = fData.size();  
  
  for (int ivar=0; ivar<nvars; ++ivar) {
    
    for (int icls=0; icls<ncls; ++icls) {
      free(_wscls[ivar][icls]);
      free(_sumwscls[ivar][icls]);
    }    
    
    free(_ws[ivar]);
    free(_ws2[ivar]);
    free(_ns[ivar]);
    free(_tgts[ivar]);
    free(_tgt2s[ivar]);
    free(_sumws[ivar]);
    free(_sumws2[ivar]);
    free(_sumns[ivar]);
    free(_sumtgts[ivar]);
    free(_sumtgt2s[ivar]);
    free(_varvals[ivar]);
    free(_bsepgains[ivar]);
    free(_bsepgainsigs[ivar]);
    
    delete[] _wscls[ivar];
    delete[] _sumwscls[ivar];
    
    free(_binquants[ivar]);
    free(_quants[ivar]);
    
    free(fQuantileMaps[ivar]);
  }    
  
  free(_sepgains);
  free(_sepgainsigs);
  free(_cutvals);
  free(_nlefts);
  free(_nrights);
  free(_sumwlefts);
  free(_sumwrights);
  free(_sumtgtlefts);
  free(_sumtgtrights);
  free(_leftvars);
  free(_rightvars);
  free(_fullvars);
  free(_bestbins);
  
  delete[] _ws;
  delete[] _ws2;
  delete[] _wscls;
  delete[] _ns;
  delete[] _nsd;
  delete[] _tgts;
  delete[] _tgt2s;
  delete[] _sumws;
  delete[] _sumws2;
  delete[] _sumwscls;
  delete[] _sumns;
  delete[] _sumtgts;
  delete[] _sumtgt2s;
  delete[] _varvals;
  delete[] _bsepgains;
  delete[] _bsepgainsigs;
  
  delete[] _binquants;
  delete[] _quants;
  
  free(_clss);
  free(_tgtvals);
  free(_tgt2vals);
  free(_weightvals);
  
  delete[] fQuantileMaps;


  for (unsigned int iev=0; iev<fEvts.size(); ++iev) {
    delete fEvts[iev];
    fEvts[iev] = 0;
  }
  
  

  
}
  


//_____________________________________________________________________________
void RooHybridBDTAutoPdf::SetMinCutSignificance(double x) {
  
  //fMinCutSignificance = TMath::ChisquareQuantile(TMath::Erf(x/sqrt(2)),1)/2.0;
  fMinCutSignificance = x*x/2.0;
  fMinCutSignificanceMulti = TMath::ChisquareQuantile(TMath::Erf(x/sqrt(2)),fNTargets)/2.0;
  
  
}

void RooHybridBDTAutoPdf::UpdateTargets(int nvars, int selvar = -1) {
 
  //// printf("UpdateTargets\n");
  
  
//  int tgtidx = itree%fNTargets;
  
  #pragma omp parallel for
  for (unsigned int iev=0; iev<fEvts.size(); ++iev) {
    
    int ithread =  omp_get_thread_num();
    
    for (unsigned int ivar=0; ivar<fCondVarsClones[ithread].size(); ++ivar) {
      fCondVarsClones[ithread][ivar]->setVal(fEvts.at(iev)->Var(ivar));
    }
    for (unsigned int ivar=0; ivar<fParmVarsClones[ithread].size(); ++ivar) {
      fParmVarsClones[ithread][ivar]->setVal(fEvts.at(iev)->Var(fCondVarsClones[ithread].size() + ivar));
    }
    
    for (unsigned int itgt=0; itgt<fStaticTgtsClones[ithread].size(); ++itgt) {
      fStaticTgtsClones[ithread][itgt]->setVal(fEvts.at(iev)->Target(itgt));
    }
    
    int evcls = fEvts.at(iev)->Class();
    double pdfval = fEvts.at(iev)->PdfVal();
    
    double invpdf = vdt::fast_inv(pdfval);
    double invpdfsq = invpdf*invpdf;
    double weight = fEvts.at(iev)->Weight();
      
    for (int itgt=0; itgt<fNTargets; ++itgt) {
      fEvts.at(iev)->SetTransTarget(itgt,0.);
      fEvts.at(iev)->SetTransTarget2(itgt,0.);
    }      

    for (unsigned int ivar=0; ivar<fFullParmsClones[ithread].size(); ++ivar) {
      fEvts.at(iev)->SetDerivative(ivar,0.);
      fEvts.at(iev)->SetDerivative2(ivar,0.);
    }          
    
    for (unsigned int iidx=0; iidx<fOuterIndices[evcls].size(); ++iidx) {
      
      int ivar = fOuterIndices[evcls][iidx];
      if (selvar>=0 && ivar!=selvar) continue;
      
      int itgt = ivar - fExtVars.getSize();
      
      RooRealVar *var = fFullParmsClones[ithread][ivar];
      double startval = var->getVal();
      //double step = 1e-3*var->getError();
      double step = 1e-3*var->getError();
      
      RooAbsReal *func = fStaticPdfsClones[ithread][evcls];
      
      var->setVal(startval + step);
      double upval = func->getValV(&fParmSetClones[ithread]);
      
      var->setVal(startval - step);
      double downval = func->getValV(&fParmSetClones[ithread]);
      
      var->setVal(startval);
      
      double drvval = (upval-downval)*vdt::fast_inv(2.0*step);
      

      fEvts.at(iev)->SetDerivative(ivar,drvval);
      
      
      double drv2val = (upval + downval - 2.0*pdfval)*vdt::fast_inv(step*step);
      fEvts.at(iev)->SetDerivative2(ivar,drv2val);
      
      if (itgt<0) continue;   
      
      fEvts.at(iev)->SetTransTarget(itgt,-weight*drvval*invpdf);
      fEvts.at(iev)->SetTransTarget2(itgt,-weight*drv2val*invpdf + weight*drvval*drvval*invpdfsq);

      
    }
    
  }  
  
  
}

void RooHybridBDTAutoPdf::BuildQuantiles(int nvars, double sumabsw) {
 
  //parallelize building of quantiles for each input variable
  //(sorting of event pointer vector is cpu-intensive)
  #pragma omp parallel for
  for (int ivar=0; ivar<nvars; ++ivar) {
    // printf("sorting var %i\n",ivar);
        
    
    std::map<int,float,std::greater<float> > tmpmap;
    std::vector<HybridGBREvent*> evtsvarsort(fEvts.begin(),fEvts.end());
    
    std::sort(evtsvarsort.begin(),evtsvarsort.end(),GBRVarCMP(ivar));
    
    double sumwq = 0;
    for (unsigned int iev=0; iev<evtsvarsort.size(); ++iev) {
      sumwq += std::abs(evtsvarsort[iev]->Weight());
      int quant = int((sumwq/sumabsw)*(fNQuantiles-1));
      float val = evtsvarsort[iev]->Var(ivar);
          
      //ensure that events with numerically identical values receive the same quantile
      if (iev>0 && val==evtsvarsort[iev-1]->Var(ivar)) quant = evtsvarsort[iev-1]->Quantile(ivar);
    
      evtsvarsort[iev]->SetQuantile(ivar,quant);
    
      tmpmap[quant] = val;
          
    }
    
    
    for (int i=0; i<fNQuantiles; ++i) {
      std::map<int,float,std::greater<float> >::const_iterator mit = tmpmap.lower_bound(i);
      
      float val;
      if (mit!=tmpmap.end()) val = mit->second;
      else val = -std::numeric_limits<float>::max();
      
      fQuantileMaps[ivar][i] = val;      
      
    }
    
    
    
  }    
  
}


void RooHybridBDTAutoPdf::TrainForest(int ntrees, bool reuseforest) {

  for (int ithread=0; ithread<fNThreads; ++ithread) {
    static_cast<RooAbsReal*>(fClones.at(ithread))->getValV(&fParmSetClones[ithread]);
  }
  
  for (int ivar=0; ivar<fTgtVars.getSize(); ++ivar) {
    static_cast<RooGBRTarget*>(fTgtVars.at(ivar))->SetUseFunc(false);  
  }
  
  for (unsigned int iev=0; iev<fEvts.size(); ++iev) {
    for (int itgt=0; itgt<fNTargets; ++itgt) {
      fEvts.at(iev)->SetCurrentNode(itgt,0);
    }
  }
  
  if (fResTree) delete fResTree;
  
  fResTree = new TNtuple("restree","","iter:r:nd:nllval:dldr");  
  
  
  
  int nvars = fCondVars.getSize() + fParmVars.getSize();
  int nvarstrain = fCondVars.getSize();
  
  double sumw = 0.;
  double sumwd = 0.;
  double sumwu = 0.;
  //  int nev = fEvts.size();
  for (std::vector<HybridGBREvent*>::iterator it=fEvts.begin(); it!=fEvts.end(); ++it) {
    sumw += (*it)->Weight();
    if ( (*it)->Class()==0 ) sumwd += (*it)->Weight();
    else if ( (*it)->Class()==1 ) sumwu += (*it)->Weight();
  }
  
  RooArgSet cvarset(fCondVars);
 
  if (!reuseforest) {
    for (int ifunc=0; ifunc<fFuncs.getSize(); ++ifunc) {
      static_cast<RooGBRFunctionFlex*>(fFuncs.at(ifunc))->ResetForest();
    }      
  }
 
  // printf("setting targets\n");
 
  //// printf("setting cloned vals\n");
  for (int ithread=0; ithread<fNThreads; ++ithread) {
    //// printf("ithread = %i\n",ithread);
    for (int ivar=0; ivar<fFullParms.getSize(); ++ivar) {
      //// printf("ivar = %i: %s\n",ivar, fFullParms.at(ivar)->GetName());
      RooRealVar *cloneparm = static_cast<RooRealVar*>(fFullParmsClones[ithread].at(ivar));
      cloneparm->setVal(static_cast<RooRealVar*>(fFullParms.at(ivar))->getVal());
      cloneparm->setError(static_cast<RooRealVar*>(fFullParms.at(ivar))->getError());
      // std::cout << "Printing " << (fFullParms.at(ivar))->GetTitle() << (fFullParms.at(ivar))->getVal() << std::endl;
    }
  }      
  
  if (!reuseforest) {
    {  
      gHybridBDTAutoPointer = this;
      RooAddition fullparmsum("fullparmsum","",fFullParms);
      RooCFunction1Binding<double,double> nllfunc("nllfunc","", &EvalLossNull,fullparmsum);
      
      nllfunc.Print();
      
      if (fDoInitialFit) {
        RooMinimizer *minim = new RooMinimizer(nllfunc);
        minim->setErrorLevel(0.5);
        minim->setStrategy(0);
        minim->minimize("Minuit2","minimize");  
        delete minim;  
      }
      
      //if (fConstraintCoeff->getVal()>0.) fR->setError(1e3*fR->getError());
    }
    
    
    for (int itgt = 0; itgt<fStaticTgts.getSize(); ++itgt) {
      double initF = static_cast<RooRealVar*>(fStaticTgts.at(itgt))->getVal();
      HybridGBRForestFlex *forest = static_cast<RooGBRTargetFlex*>(fTgtVars.at(itgt))->Forest();
      forest->SetInitialResponse(initF);
      for (std::vector<HybridGBREvent*>::iterator it=fEvts.begin(); it!=fEvts.end(); ++it) {
	(*it)->SetTarget(itgt,initF);
      }    
    }  
    
    for (int ithread=0; ithread<fNThreads; ++ithread) {
      for (int ivar=0; ivar<fFullParms.getSize(); ++ivar) {
	RooRealVar *cloneparm = static_cast<RooRealVar*>(fFullParmsClones[ithread].at(ivar));
	cloneparm->setVal(static_cast<RooRealVar*>(fFullParms.at(ivar))->getVal());
	cloneparm->setError(static_cast<RooRealVar*>(fFullParms.at(ivar))->getError());
      }
    }
  }
  
  
  fStepSizes.resize(fFullParms.getSize());  
  for (int iparm=0; iparm<fFullParms.getSize(); ++iparm) {
    fStepSizes[iparm] = 1e-3*static_cast<RooRealVar*>(fFullParms.at(iparm))->getError();
  }
  
  
  // printf("nev = %i, sumw = %5f\n",int(nev), sumw);
  
  std::vector<std::pair<float,float> > limits;
  for (int ivar=0; ivar<nvarstrain; ++ivar) {
    //limits.push_back(std::pair<float,float>(-std::numeric_limits<float>::max(),std::numeric_limits<float>::max()));
    limits.push_back(std::pair<float,float>(-10.,10.));
  }
  
  RooArgSet parmset(fParmVars);

  
  TVectorD dparnull(fFullParms.getSize());
  fNLLVal = EvalLoss(0.,dparnull);

  // printf("Initial fNLLVal = %5f\n",fNLLVal);
  
  // printf("fullparms:\n");
  fFullParms.Print("V");
  
  TRandom3 rndtgtsel(fEvts.size());
  
  //loop over requested number of trees
  int nunittrees = 0;
  int nsmalltrees = 0;
  std::vector<double> nllvals;
  std::vector<double> dldrvals;
  for (int itree=0; itree<ntrees; ++itree) {
    // printf("tree %i\n",itree);

    bool sharedfuncs = false;
    
    for (unsigned int iev=0; iev<fEvts.size(); ++iev) {
      for (int itgt=0; itgt<fNTargets; ++itgt) {
        fEvts.at(iev)->SetCurrentNode(itgt,0);
      }
    }    

    int maxtreesize = 0;
    //    int maxtreeidx = 0;
    
    std::vector<int> treesizes(fFuncs.getSize());

    for (int ifunc=0; ifunc<fFuncs.getSize(); ++ifunc) {
      HybridGBRForestFlex *forest = static_cast<RooGBRFunctionFlex*>(fFuncs.at(ifunc))->Forest();
      forest->Trees().push_back(HybridGBRTreeD()); 
    }
        
    UpdateTargets(nvars,-1);       
    
    for (int ifunc=0; ifunc<fFuncs.getSize(); ++ifunc) {
      HybridGBRForestFlex *forest = static_cast<RooGBRFunctionFlex*>(fFuncs.at(ifunc))->Forest();
      HybridGBRTreeD &tree = forest->Trees().back();     
      
      //modulate through targets associated to this function
      int nfunctgts = fFuncTgts[ifunc].getSize();
      int seltgt = itree%nfunctgts;
      //int seltgt = rndtgtsel.Integer(nfunctgts);
      int tgtidx = fTgtVars.index(fFuncTgts[ifunc].at(seltgt));
      TrainTree(fEvts,sumw,tree,0.,0,limits,tgtidx);
      
      int treesize = tree.Responses().size();
      treesizes[ifunc] = treesize;           
      
      //Set Terminal nodes for unselected targets belonging to this function
      for (int itgt=0; itgt<nfunctgts; ++itgt) {
        if (itgt!=seltgt) {
          int updtgt = fTgtVars.index(fFuncTgts[ifunc].at(itgt));
          UpdateCurrentNodes(fEvts,tree,updtgt);
          sharedfuncs = true;
        }
      }      
      

    }
    
    FitResponses(-1);
     
    for (int ifunc=0; ifunc<fFuncs.getSize(); ++ifunc) {
      int treesize = treesizes[ifunc];
      if (treesize>maxtreesize) {
	maxtreesize = treesize;
	//        maxtreeidx = ifunc;
      }      
    }
    
    
    // printf("maxtreesize = %i, maxtgtidx = %i\n",maxtreesize,maxtreeidx);

    if (maxtreesize==1) {
      ++nunittrees;
    }
    else {
      nunittrees = 0.;
    }
    
    if (maxtreesize>1 && maxtreesize<16) {
      ++nsmalltrees;
    }
    else {
      nsmalltrees = 0;
    }
    
    //FitResponses(forest);

    
    nllvals.push_back(fNLLVal);
    
    double dldrval = -2.0*fConstraintCoeff->getVal()*(fR->getVal()-fConstraintVal->getVal());
    dldrvals.push_back(dldrval);
    
    fResTree->Fill(itree,fR->getVal(),fN0->getVal(),fNLLVal,fdLdR);
   
    int oldnllidx = nllvals.size() - 5 - 1;
    //if (oldnllidx>=0 && (fNLLVal - nllvals[oldnllidx])>(-2e-3) && std::abs(dldrval-dldrvals[oldnllidx])<2e-1 && nunittrees>10) {
    if (oldnllidx>=0 && (fNLLVal - nllvals[oldnllidx])>(-2e-3) && nunittrees>10) {      
      break;
    }
    
    if (nunittrees>100) {
      // printf("Max number of unit trees %i exceeded, breaking\n",nunittrees);
      break;
    }

    int voldnllidx = nllvals.size() - 20 - 1;
    //if (oldnllidx>=0 && (fNLLVal - nllvals[oldnllidx])>(-2e-3) && std::abs(dldrval-dldrvals[oldnllidx])<2e-1 && nunittrees>10) {
    if (voldnllidx>=0 && (fNLLVal - nllvals[voldnllidx])>(-2e-3) && sharedfuncs) {      
      break;
    }    
    
    if (0) {
      if (sharedfuncs && oldnllidx>=0 && (fNLLVal - nllvals[oldnllidx])>(-1.0) && nsmalltrees>10) {
        fMinWeightTotal = 2.0*sumw;
      }
    }
    
    if (0) {
      oldnllidx = nllvals.size() - 3.0/fShrinkage - 1;
      //if (oldnllidx>=0 && (fNLLVal - nllvals[oldnllidx])>(-2e-3) && std::abs(dldrval-dldrvals[oldnllidx])<2e-1) {
      if (oldnllidx>=0 && (fNLLVal - nllvals[oldnllidx])>(-2e-3)) {
	// printf("breaking\n");
	break;
      }
    }

  }
   
  for (int ifunc=0; ifunc<fFuncs.getSize(); ++ifunc) {
    fFuncs.at(ifunc)->setValueDirty();
  }        
  
  for (int ivar=0; ivar<fTgtVars.getSize(); ++ivar) {
    static_cast<RooGBRTarget*>(fTgtVars.at(ivar))->SetUseFunc(true);  
  }  
  
  
}
 
    
// void RooHybridBDTAutoPdf::vtest(const float **__restrict__ ra, const float **__restrict__ rb, float **__restrict__ rc) {
//     
//   for (int iat=0; iat<128; ++iat) {
// //     const float *__restrict__ sra = (const float*)__builtin_assume_aligned(ra[iat],32);
// //     const float *__restrict__ srb = (const float*)__builtin_assume_aligned(rb[iat],32);
// //     float *__restrict__ src = (float*)__builtin_assume_aligned(rc[iat],32);    
//     //rc[iat] = ra[iat] + rb[iat];
//     for (int jat=0; jat<128; ++jat) {
//       //src[jat] = sra[jat] + srb[jat];
//       rc[iat][jat] = ra[iat][jat] + rb[iat][jat];
//     }
//   }
//        
// }  
    
  
//_______________________________________________________________________
void RooHybridBDTAutoPdf::TrainTree(const std::vector<HybridGBREvent*> &evts, double sumwtotal, HybridGBRTreeD &tree, double transition, int depth, const std::vector<std::pair<float,float> > limits, int tgtidx) {
  
  const RooGBRTargetFlex *target = static_cast<RooGBRTargetFlex*>(fTgtVars.at(tgtidx));
  
  int nvars = target->FuncVars().getSize();  

  int thisidx = tree.CutIndices().size();    
  
  //number of events input to node
  const int nev = evts.size();
  const int ncls = fStaticPdfs.getSize();
  
  //index of best cut variable
  int bestvar = 0;
   


  //build map to global variable index
  std::vector<int> varidxs(nvars);
  for (int ivar=0; ivar<nvars; ++ivar) {
    varidxs[ivar] = fCondVars.index(target->FuncVars().at(ivar));
  }
  
  //// printf("first parallel loop\n");
  //// printf("nev = %i\n",nev);
  //fill temporary array of quantiles (to allow auto-vectorization of later loops)
  #pragma omp parallel for
  for (int iev = 0; iev<nev; ++iev) {
    _clss[iev] = evts[iev]->Class();
    _tgtvals[iev] = evts[iev]->TransTarget(tgtidx);
    _tgt2vals[iev] = evts[iev]->TransTarget2(tgtidx);
    _weightvals[iev] = evts[iev]->Weight();
    for (int ivar=0; ivar<nvars; ++ivar) {
      _quants[ivar][iev] = evts[iev]->Quantile(varidxs[ivar]);   
    }
  }  
    
  //// printf("second parallel loop\n");
  //trivial open-mp based multithreading of loop over input variables
  //The loop is thread safe since each iteration writes into its own
  //elements of the 2-d arrays
  #pragma omp parallel for schedule(dynamic,1)
  for (int ivar=0; ivar<nvars; ++ivar) {
             
    
    int minquant;
    int maxquant;
    
    //find max and min quantiles in the input events
    //(this loop should be vectorized by gcc with reasonable optimization options)
    GBRArrayUtils::MinMaxQuants(minquant, maxquant, _quants[ivar], nev);

    //calculate offset and scaling (powers of 2) to reduce the total number of quantiles
    //to the fNBinsMax for the search for the best split value
    int offset = minquant;
    unsigned int bincount = maxquant-minquant+1;
    unsigned int pscale = 0;
    while (bincount>fNBinsMax) {
      ++pscale;
      //bincount >>= 1;
      bincount = ((maxquant-offset)>>pscale) + 1;
    }    

    const unsigned int nbins = ((maxquant-offset)>>pscale) + 1;
    assert(nbins<=fNBinsMax);
        
    //zero arrays where necessary and 
    GBRArrayUtils::InitArrays(_ns[ivar],_tgts[ivar],_tgt2s[ivar],_bsepgains[ivar],nbins);
    
    //// printf("touch wscls\n");
    //the inner loop here should also vectorize
    for (int icls=0; icls<ncls; ++icls) {
      GBRArrayUtils::ZeroArray(_wscls[ivar][icls],nbins);
    }    
    //// printf("done wscls\n");
    
    
    //compute map between bin numbers
    //and variable cut values
    
    //// printf("quant manipulation\n");
    //this loop should auto-vectorize
    GBRArrayUtils::FillBinQuants(_binquants[ivar], offset, pscale,fNQuantiles, nbins);
    
    //this loop won't auto-vectorize because it's another gather operation (maybe in avx2 and gcc>4.7)
    for (unsigned int ibin=0; ibin<nbins; ++ibin) { 
      int quant = _binquants[ivar][ibin];
      _varvals[ivar][ibin] = fQuantileMaps[varidxs[ivar]][quant];
    }
    //// printf("done quant manipulation\n");
    
    //// printf("compute bins\n");
   

      
    //// printf("filling histogram-style arrays\n");
     
    //compute summed quantities differential in each bin
    //(filling 'histograms')
    //This loop is one of the most expensive in the algorithm for large training samples
    //This loop can unfortunately not be vectorized because the memory addressed 
    //are computed within the loop iteration
    //JOSH: Once the data-dependency checks are appropriately bypassed, this should actually vectorize, but only
    //for avx2 targets which support vectorized gather instructions
//    int nevd = 0;
    for (int iev=0;iev<nev;++iev) {
      //if (!evts[iev]->IsPrimary()) continue;
      
      int icls = _clss[iev];
      int ibin = (_quants[ivar][iev]-offset)>>pscale;

      ++_ns[ivar][ibin];
      
      _wscls[ivar][icls][ibin] += _weightvals[iev];            
      //// printf("done incrementing wscls\n");
      
      _tgts[ivar][ibin] += _tgtvals[iev];
      _tgt2s[ivar][ibin] += _tgt2vals[iev];
      
    }
    
      
    //// printf("starting split search\n");
 
    //// printf("compute array integrals\n");
    //convert differential arrays to cumulative arrays by summing over
    //each element
    //loop cannot be vectorized because this is an iterative calculation
    _sumws[ivar][0] = _ws[ivar][0];
    _sumws2[ivar][0] = _ws2[ivar][0];
    for (int icls=0; icls<ncls; ++icls) {
      _sumwscls[ivar][icls][0] = _wscls[ivar][icls][0];
    }
    _sumns[ivar][0] = _ns[ivar][0];
    _sumtgts[ivar][0] = _tgts[ivar][0];
    _sumtgt2s[ivar][0] = _tgt2s[ivar][0];    
    
    for (unsigned int ibin=1; ibin<nbins; ++ibin) {      
      _sumns[ivar][ibin] = _sumns[ivar][ibin-1] + _ns[ivar][ibin];
      _sumtgts[ivar][ibin] = _sumtgts[ivar][ibin-1] + _tgts[ivar][ibin];  
      _sumtgt2s[ivar][ibin] = _sumtgt2s[ivar][ibin-1] + _tgt2s[ivar][ibin];  
    }
    
    //// printf("cls array integrals\n");
    for (int icls=0; icls<ncls; ++icls) {
      for (unsigned int ibin=1; ibin<nbins; ++ibin) {  
	_sumwscls[ivar][icls][ibin] = _sumwscls[ivar][icls][ibin-1] + _wscls[ivar][icls][ibin];
      }
    }    
   // // printf("done cls array integrals\n");
    

    const double sumtgt = _sumtgts[ivar][nbins-1];
    
    const double sumtgt2 = _sumtgt2s[ivar][nbins-1];      

    
    //// printf("short loop\n");
    float maxsepgain = 0.;
    float cutval = 0.;
    int nleft= 0;
    int nright = 0;
    int bestbin=0;
    
    const double fulldiff = std::min(0.,-0.5*sumtgt*sumtgt*vdt::fast_inv(sumtgt2));
    
    //// printf("start heavy loop\n");
    //loop over all bins and compute improvement in weighted variance of target for each split
    //This loop is relatively expensive and should auto-vectorize in the appropriate compiler/settings
    GBRArrayUtils::FillSepGains(_sumtgts[ivar], _sumtgt2s[ivar], _bsepgains[ivar], fulldiff, sumtgt, sumtgt2, nbins);
     
    //// printf("start final loop\n");
    //loop over computed variance improvements and select best split, respecting also minimum number of events per node
    //This loop cannot auto-vectorize, at least in gcc 4.6x due to the mixed type conditional, but it's relatively fast
    //in any case
    for (unsigned int ibin=0; ibin<nbins; ++ibin) {   

      if ( _bsepgains[ivar][ibin]>maxsepgain && std::isnormal(_bsepgains[ivar][ibin])) {
	
	bool passminweights = true;
	double minweights = std::numeric_limits<double>::max();
	double totalweightleft = 0;
	double totalweightright = 0;
	for (int icls=0; icls<ncls; ++icls) {
	  double minweightcls = std::min(_sumwscls[ivar][icls][ibin], _sumwscls[ivar][icls][nbins-1] - _sumwscls[ivar][icls][ibin]);
	  if (fMinWeights.size() && minweightcls < fMinWeights[icls]) {
	    passminweights = false;
	  }
	  if (minweightcls < minweights) {
	    minweights = minweightcls;
	  }
	  totalweightleft += _sumwscls[ivar][icls][ibin];
	  totalweightright += _sumwscls[ivar][icls][nbins-1] - _sumwscls[ivar][icls][ibin];
	}
	
	if (fMinWeightTotal>=0. && (totalweightleft<fMinWeightTotal || totalweightright<fMinWeightTotal) ) {
	  passminweights = false;
	}
	
	
	bool passminsig = true;
	
	passminsig &= fMinCutSignificance<0. || (_bsepgains[ivar][ibin] > fMinCutSignificance);
	
        if (fMaxNSpurious >= 0.) {
          double prob = 1.0 - TMath::Erf(sqrt(_bsepgains[ivar][ibin]));
          double nspurious = fSumWTimesNVars*prob;	
          passminsig &= nspurious<fMaxNSpurious;	
        }
	
	if (passminweights && passminsig) {
	  maxsepgain = _bsepgains[ivar][ibin];
	  bestbin = ibin;
	}
      }
      
    }
     
    cutval = _varvals[ivar][bestbin];
    nleft = _sumns[ivar][bestbin];
    nright = nev - nleft;
    
    _sepgains[ivar] = maxsepgain;
    _sepgainsigs[ivar] = _bsepgainsigs[ivar][bestbin];
    _cutvals[ivar] = cutval;
    _nlefts[ivar] = nleft;
    _nrights[ivar] = nright;
    _bestbins[ivar] = bestbin;
        
     //// printf("done var %i\n",ivar);
  }
  

  
  float globalsepgain = 0.;
  for (int ivar=0; ivar<nvars; ++ivar) {
    if (_sepgains[ivar]>globalsepgain) {
      globalsepgain = _sepgains[ivar];
      bestvar = ivar;
    }
  }    
    
  //if no appropriate split found, make this node terminal
  if (globalsepgain<=0.) {
    //no valid split found, making this node a leaf node
    //// printf("thisidx = %i, globalsepgain = %5f, no valid split\n",thisidx, globalsepgain);
    tree.CutIndices().push_back(0);
    tree.CutVals().push_back(0);
    tree.LeftIndices().push_back(0);   
    tree.RightIndices().push_back(0);    
    
    tree.RightIndices()[thisidx] = -tree.Responses().size();
    tree.LeftIndices()[thisidx] = -tree.Responses().size();
    
    BuildLeaf(evts,tree,tgtidx);
    return;
  }
  
  //fill vectors of event pointers for left and right nodes below this one
  std::vector<HybridGBREvent*> leftevts;
  std::vector<HybridGBREvent*> rightevts;
  
  leftevts.reserve(nev);
  rightevts.reserve(nev);
  
  int nleft = 0;
  int nright = 0;
  double sumwleft = 0.;
  double sumwright = 0.;
  
  for (std::vector<HybridGBREvent*>::const_iterator it = evts.begin(); it!=evts.end(); ++it) {
    if ((*it)->Var(varidxs[bestvar])>_cutvals[bestvar]) {
      ++nright;
      sumwright += (*it)->Weight();
      rightevts.push_back(*it);
    }
    else {
      ++nleft;
      sumwleft += (*it)->Weight();
      leftevts.push_back(*it);
    }
  }
  
//   // printf("bestvar = %i, nleft = %i, nright = %i, nlefts[bestvar] = %i, nrights[bestvar] = %i, cutvals[bestvar] = %5f, nmismatch = %i\n",bestvar,nleft,nright,_nlefts[bestvar],_nrights[bestvar],_cutvals[bestvar],nmismatch);

  //nleft and nright are just the number of events in a given bin, nleft < cut value, nright > cut value
  //
  if(_nlefts[bestvar]!=nleft || _nrights[bestvar]!=nright){
    printf("when applying the cut value %f for this bin on variable number %i, the count obtained previously and the count just recalculated do not match\n",_cutvals[bestvar],bestvar);
    printf("this can happen if variable number %i is NaN, check for this (you look need to look a the config to map the variable number to variable name)\n",bestvar);
  }
         
  assert(_nlefts[bestvar]==nleft); 
  assert(_nrights[bestvar]==nright);
  
  //// printf("nleft = %i, nright = %i\n",nleft,nright);
  
  
  double bestcutval = _cutvals[bestvar];  
  
  //fill intermediate node
  tree.CutIndices().push_back(bestvar);
  tree.CutVals().push_back(_cutvals[bestvar]);
  tree.LeftIndices().push_back(0);   
  tree.RightIndices().push_back(0);  
  
  //check if left node is terminal
  //bool termleft = nleft<=(2*fMinEvents) || depth==fMaxDepth;
  bool termleft = sumwleft<=(2*fMinEvents) || (fMaxDepth>=0 && depth==fMaxDepth) || (fMaxNodes>=0 && int(tree.Responses().size())>=fMaxNodes) ;
  if (termleft) tree.LeftIndices()[thisidx] = -tree.Responses().size();
  else tree.LeftIndices()[thisidx] = tree.CutIndices().size();
  
  //// printf("this idx = %i, termleft = %i, nleft = %i, fMinEvents = %i\n",thisidx,  termleft,nleft,fMinEvents);  
  
  //build left node as appropriate
  std::vector<std::pair<float,float> > limitsleft(limits);
  limitsleft[bestvar].second = bestcutval;  
  //// printf("bestvar = %i, limlow = %5f, limhigh = %5f, limleftlow = %5f, limlefthigh = %5f\n",bestvar,limits[bestvar].first,limits[bestvar].second,limitsleft[bestvar].first,limitsleft[bestvar].second);
  if (termleft) {  
    BuildLeaf(leftevts,tree,tgtidx);
  }
  else {  
    TrainTree(leftevts,sumwleft,tree,transition,depth+1,limitsleft,tgtidx);  
  }
  
  //check if right node is terminal
  //bool termright = nright<=(2*fMinEvents) || depth==fMaxDepth;
  bool termright = sumwright<=(2*fMinEvents) || (fMaxDepth>=0 && depth==fMaxDepth) || (fMaxNodes>=0 && int(tree.Responses().size())>=fMaxNodes);
  if (termright) tree.RightIndices()[thisidx] = -tree.Responses().size();
  else tree.RightIndices()[thisidx] = tree.CutIndices().size();
    
  //// printf("this idx = %i, termright = %i, nright = %i, fMinEvents = %i\n",thisidx,  termright,nright,fMinEvents);    
  
  //build right node as appropriate
  std::vector<std::pair<float,float> > limitsright(limits);
  limitsright[bestvar].first = bestcutval;  
  //// printf("bestvar = %i, limlow = %5f, limhigh = %5f, limrightlow = %5f, limrighthigh = %5f\n",bestvar, limits[bestvar].first,limits[bestvar].second,limitsright[bestvar].first,limitsright[bestvar].second);  
  if (termright) {  
    BuildLeaf(rightevts,tree,tgtidx);
  }
  else {      
    TrainTree(rightevts,sumwright,tree,transition,depth+1,limitsright, tgtidx);  
  }
  
}

  

  


//_______________________________________________________________________
void RooHybridBDTAutoPdf::BuildLeaf(const std::vector<HybridGBREvent*> &evts, HybridGBRTreeD &tree, int tgtidx) {

  //// printf("building leaf\n");
  
  int thisidx = -tree.Responses().size();
    
  tree.Responses().push_back(0.);
  
  for (std::vector<HybridGBREvent*>::const_iterator it = evts.begin(); it!=evts.end(); ++it) {
    (*it)->SetCurrentNode(tgtidx, -thisidx);
  }
  
}

//_______________________________________________________________________
void RooHybridBDTAutoPdf::UpdateCurrentNodes(const std::vector<HybridGBREvent*> &evts, HybridGBRTreeD &tree, int tgtidx) {
  
  //update current node directly from GBRTree (needed for targets which have been skipped during tree growth)
  
  const RooGBRTargetFlex *target = static_cast<RooGBRTargetFlex*>(fTgtVars.at(tgtidx));
  int nvars = target->FuncVars().getSize();
  
  //build map to global variable index
  std::vector<int> varidxs(nvars);
  for (int ivar=0; ivar<nvars; ++ivar) {
    varidxs[ivar] = fCondVars.index(target->FuncVars().at(ivar));
  }  
  
  std::vector<std::vector<float> > evals(fNThreads,std::vector<float>(nvars,0.));

  int nev = evts.size();  
  
  #pragma omp parallel for
  for (int iev=0; iev<nev; ++iev) {    
    int ithread =  omp_get_thread_num();

    for (int ivar=0; ivar<nvars; ++ivar) {
      evals[ithread][ivar] = evts[iev]->Var(varidxs[ivar]);
    }
    
    int termidx = tree.TerminalIndex(&evals[ithread][0]);
    
    assert(termidx<int(tree.Responses().size()));
    
    evts[iev]->SetCurrentNode(tgtidx, termidx);
  }
  
}



double RooHybridBDTAutoPdf::EvalLossNull(double dummy) {
 
  return gHybridBDTAutoPointer->EvalLossRooFit();
  
}

double RooHybridBDTAutoPdf::EvalLossRooFit() {
 
  RooAbsReal::setEvalErrorLoggingMode(RooAbsReal::CountErrors);
  RooAbsReal::clearEvalErrorLog() ;  

  for (int ithread=0; ithread<fNThreads; ++ithread) {
    for (int ivar=0; ivar<fFullParms.getSize(); ++ivar) {
      RooRealVar *cloneparm = static_cast<RooRealVar*>(fFullParmsClones[ithread].at(ivar));
      cloneparm->setVal(static_cast<RooRealVar*>(fFullParms.at(ivar))->getVal());
      cloneparm->setError(static_cast<RooRealVar*>(fFullParms.at(ivar))->getError());
    }
  } 

  double nllval = 0.;
  //RooArgSet parmset(fParmVars);
  
  std::vector<double> nllvals(fNThreads);
	
  #pragma omp parallel for
  for (unsigned int ievt=0; ievt<fEvts.size(); ++ievt) {

    //if (ievt%20!=0) continue;
    //if (ievt%100!=0) continue;
    
    if (fPrescaleInit>0 && ievt%fPrescaleInit!=0) continue;
    
    int ithread =  omp_get_thread_num();
    //int ithread =  0;

    //int termidx = fEvts[ievt]->CurrentNode();
    //if (seltermidx>=0 && termidx!=seltermidx) continue;
        
    for (unsigned int ivar=0; ivar<fCondVarsClones[ithread].size(); ++ivar) {
      fCondVarsClones[ithread][ivar]->setVal(fEvts[ievt]->Var(ivar));
    }
    for (unsigned int ivar=0; ivar<fParmVarsClones[ithread].size(); ++ivar) {
      fParmVarsClones[ithread][ivar]->setVal(fEvts[ievt]->Var(fCondVarsClones[ithread].size() + ivar));
    }    

    
    int evcls = fEvts[ievt]->Class(); 
    double weight = fEvts[ievt]->Weight();

    //nll value for minos constraint
    //if (evcls==0) fNLLVal += -log(static_cast<RooAbsReal*>(fStaticPdfs.at(evcls))->getValV(&parmset));
    
    double pdfval = fStaticPdfsClones[ithread][evcls]->getValV(&fParmSetClones[ithread]);
    
    //nllval += -weight*log(pdfval);
    //nllvals[ithread] += -weight*vdt::fast_logf(pdfval);
    nllvals[ithread] += -weight*log(pdfval);
    
    if (RooAbsReal::numEvalErrors()>0  || pdfval<0.) {
      nllvals[ithread] += std::numeric_limits<float>::max();
    }
    
/*    if (RooAbsReal::numEvalErrors()>0 || RooAbsPdf::evalError() || pdfval<0.) {
    //if (RooAbsReal::numEvalErrors()>0 || RooAbsPdf::evalError()) {  
      //// printf("numEvalErrors = %i\n",RooAbsReal::numEvalErrors());
      RooAbsReal::clearEvalErrorLog();
      RooAbsPdf::clearEvalError();
      // printf("eval error in EvalLoss\n");
      RooAbsReal::setEvalErrorLoggingMode(RooAbsReal::PrintErrors);
      return std::numeric_limits<double>::max();
    }   */ 
    
    
  }
  
  for (int ithread=0; ithread<fNThreads; ++ithread) {
    nllval += nllvals[ithread];
  }
  
  int infunc = fFullFuncs.getSize()-1;
  nllval += static_cast<RooAbsReal*>(fFullFuncs.at(infunc))->getVal();

  RooAbsReal::setEvalErrorLoggingMode(RooAbsReal::PrintErrors);
  
  return nllval;  
  
}

double RooHybridBDTAutoPdf::EvalLoss(double lambda, const TVectorD &dL, int itree) {
 
 // int tgtidx = itree%fNTargets;
  
  RooAbsReal::setEvalErrorLoggingMode(RooAbsReal::CountErrors);
  RooAbsReal::clearEvalErrorLog() ;  
  
//   int nterm = tree.Responses()[0].size();
// 
  int nextvars = fExtVars.getSize();
// 
  
  int nparms = nextvars;
  std::vector<int> localidxs(fFuncs.getSize());
  for (int ifunc=0; ifunc<fFuncs.getSize(); ++ifunc) {
    localidxs[ifunc] = nparms;
    HybridGBRForestFlex *forest = static_cast<RooGBRFunctionFlex*>(fFuncs.at(ifunc))->Forest();
        
    if (forest->Trees().size() && forest->Trees().back().Responses().size()) {
      nparms += forest->Trees().back().Responses().size();
    }
    else {
      nparms += 1;    
    }
  }   
   
//   int nparms = fExtVars.getSize() + fNTargets*nterm; 
  
  //// printf("nextvars = %i, nparms = %i\n",nextvars,nparms);
  
 // double dr = 0.;
  
  int idxglobal = 0;
  //int idxn0 = idxglobal + fExtVars.index(fN0);
  
  //// printf("idxn0 = %i, fN0 = %5f\n", idxn0,fN0->getVal());

  double nllval = 0.;
  
  
  //// printf("initval = %5f\n",static_cast<RooRealVar*>(fExtVars.at(1))->getVal());
  
  //global variables
//   RooArgList extbak = fExtVars;
//   
  std::vector<double> extvals(fExtVars.getSize());
  
  for (int ivar=0; ivar<fExtVars.getSize(); ++ivar) {
    int iel = idxglobal + ivar;
    extvals[ivar] = static_cast<RooRealVar*>(fExtVars.at(ivar))->getVal();
    static_cast<RooRealVar*>(fExtVars.at(ivar))->setVal(static_cast<RooRealVar*>(fExtVars.at(ivar))->getVal() + lambda*dL(iel));   
  }
  
  for (int ithread=0; ithread<fNThreads; ++ithread) {
    for (int ivar=0; ivar<fExtVars.getSize(); ++ivar) {
      RooRealVar *cloneparm = static_cast<RooRealVar*>(fExtVarsClones[ithread].at(ivar));
      cloneparm->setVal(static_cast<RooRealVar*>(fExtVars.at(ivar))->getVal());
      cloneparm->setError(static_cast<RooRealVar*>(fExtVars.at(ivar))->getError());
    }
  }      
  
//   // printf("intermediateval = %5f\n",static_cast<RooRealVar*>(fExtVars.at(1))->getVal());
//   
//   // printf("begin loop\n");

//   int pointidx = -1;
//   if (itree>=0) {
//     pointidx = itree%fFullParms.getSize();
//   }


  //map of function indices
  std::vector<int> funcidxs(fStaticTgts.getSize());
  for (int itgt=0; itgt<fStaticTgts.getSize(); ++itgt) {
    funcidxs[itgt] = fFuncs.index(static_cast<RooGBRTargetFlex*>(fTgtVars.at(itgt))->Func());
  }

  std::vector<double> nllvals(fNThreads);

  #pragma omp parallel for
  for (unsigned int ievt=0; ievt<fEvts.size(); ++ievt) {

    int ithread =  omp_get_thread_num();
    //int ithread =  0;
    
    
    //int termidx = fEvts[ievt]->CurrentNode();
    //if (seltermidx>=0 && termidx!=seltermidx) continue;
    
    //int idxlocal = nextvars + fNTargets*termidx;
    
    for (unsigned int ivar=0; ivar<fCondVarsClones[ithread].size(); ++ivar) {
      fCondVarsClones[ithread][ivar]->setVal(fEvts[ievt]->Var(ivar));
    }
    for (unsigned int ivar=0; ivar<fParmVarsClones[ithread].size(); ++ivar) {
      fParmVarsClones[ithread][ivar]->setVal(fEvts[ievt]->Var(fCondVarsClones[ithread].size() + ivar));
    }    

    for (unsigned int itgt=0; itgt<fStaticTgtsClones[ithread].size(); ++itgt) {
      int ifunc = funcidxs[itgt];
      int iel = localidxs[ifunc] + fEvts[ievt]->CurrentNode(itgt);
      fStaticTgtsClones[ithread][itgt]->setVal(fEvts[ievt]->Target(itgt) + lambda*dL[iel]);
    }
    
    int evcls = fEvts[ievt]->Class(); 
    double weight = fEvts[ievt]->Weight();

    //nll value for minos constraint
    //if (evcls==0) fNLLVal += -log(static_cast<RooAbsReal*>(fStaticPdfs.at(evcls))->getValV(&parmset));
    
    double pdfval = fStaticPdfsClones[ithread][evcls]->getValV(&fParmSetClones[ithread]);
    fEvts[ievt]->SetPdfVal(pdfval);
    
    //if (!std::isnormal(pdfval)) fStaticPdfsClones[ithread].at(evcls)->Print("V");
    
    //nllval += -weight*log(pdfval);
    //nllvals[ithread] += -weight*vdt::fast_logf(pdfval);
    nllvals[ithread] += -weight*log(pdfval); 
    
    
//     if (pointidx>=0) {
//       fEvts[ievt]->ValVector()(pointidx) = pdfval;
//       for (int iparm=0; iparm<fFullParms.getSize(); ++iparm) {
// 	fEvts[ievt]->ParmMatrix()(iparm,pointidx) = static_cast<RooRealVar*>(fFullParms.at(iparm))->getVal();
//       }
//     }
    
//     if (RooAbsReal::numEvalErrors()>0) {
//       RooAbsReal::clearEvalErrorLog();
//       return std::numeric_limits<double>::max();
//     }
    
//     double testmass = fEvts[ievt]->Var(fCondVars.getSize());
//     if (testmass>178.546 && testmass<178.548)
//       // printf("evcls = %i, pdfval = %5e, logmode = %i\n", evcls,pdfval,int(RooAbsReal::evalErrorLoggingMode()));
    
    if (RooAbsReal::numEvalErrors()>0  || pdfval<0.) {
      nllvals[ithread] += std::numeric_limits<float>::max();
    }
    
    
/*    if (RooAbsReal::numEvalErrors()>0 || RooAbsPdf::evalError() || pdfval<0.) {
    //if (RooAbsReal::numEvalErrors()>0 || RooAbsPdf::evalError()) {  
      //// printf("numEvalErrors = %i\n",RooAbsReal::numEvalErrors());
      RooAbsReal::clearEvalErrorLog();
      RooAbsPdf::clearEvalError();
      // printf("eval error in EvalLoss\n");
      for (int ivar=0; ivar<fExtVars.getSize(); ++ivar) {
        static_cast<RooRealVar*>(fExtVars.at(ivar))->setVal(extvals[ivar]);   
      }        
      for (int ithread=0; ithread<fNThreads; ++ithread) {
	for (int ivar=0; ivar<fExtVars.getSize(); ++ivar) {
	  RooRealVar *cloneparm = static_cast<RooRealVar*>(fExtVarsClones[ithread].at(ivar));
	  cloneparm->setVal(static_cast<RooRealVar*>(fExtVars.at(ivar))->getVal());
	  cloneparm->setError(static_cast<RooRealVar*>(fExtVars.at(ivar))->getError());
	}
      }          
      RooAbsReal::setEvalErrorLoggingMode(RooAbsReal::PrintErrors);
      return std::numeric_limits<double>::max();
    }*/    
    
    //if (std::isnan(nllval) || std::isinf(nllval) || nllval==0.) return std::numeric_limits<double>::max();
    
    //normalization terms
//     double normden = 1.0;
//     for (int ipdf=1; ipdf<fStaticPdfs.getSize(); ++ipdf) {
//       int itgt = ipdf - 1;
//       int iel = idxlocal + itgt;
//       double tgtval = fEvts[ievt]->Target(itgt) + lambda*dL[iel];
//       double expF = exp(tgtval);
//       normden += expF;
//       
//       if (ipdf==evcls) {
//         nllval += -weight*tgtval;
//       }
//     }
//     nllval += weight*log(normden);
    
    
/*    if (evcls==0) {
      for (int ipdf=1; ipdf<fStaticPdfs.getSize(); ++ipdf) {
	int itgt = ipdf - 1;
	int iel = idxlocal + itgt;
	double nexpF = exp(-fEvts[ievt]->Target(itgt) - lambda*dL[iel]);
	//nllval += weight*nexpF/fN0Obs;
	nllval += weight*nexpF;
      }
    }
    else {
      int itgt = evcls - 1;
      int iel = idxlocal + itgt;
      nllval += weight*(fEvts[ievt]->Target(itgt) + lambda*dL[iel]);      
    }*/    
    
  }
  
  for (int ithread=0; ithread<fNThreads; ++ithread) {
    nllval += nllvals[ithread];
  }
  
  //global terms
  //if (seltermidx<0) {
  int infunc = fFullFuncs.getSize()-1;
  nllval += static_cast<RooAbsReal*>(fFullFuncs.at(infunc))->getVal();
  //nllval += static_cast<RooAbsReal*>(fFullFuncs.at(infunc))->getVal() - fN0Obs*log(static_cast<RooAbsReal*>(fFullFuncs.at(infunc))->getVal());  
  //}
  
  //fExtVars = extbak;
  
  for (int ivar=0; ivar<fExtVars.getSize(); ++ivar) {
    static_cast<RooRealVar*>(fExtVars.at(ivar))->setVal(extvals[ivar]);  
  }  
  for (int ithread=0; ithread<fNThreads; ++ithread) {
    for (int ivar=0; ivar<fExtVars.getSize(); ++ivar) {
      RooRealVar *cloneparm = static_cast<RooRealVar*>(fExtVarsClones[ithread].at(ivar));
      cloneparm->setVal(static_cast<RooRealVar*>(fExtVars.at(ivar))->getVal());
      cloneparm->setError(static_cast<RooRealVar*>(fExtVars.at(ivar))->getError());
    }
  }      

//   for (unsigned int itgt=0; itgt<fStaticTgtsClones[0].size(); ++itgt) {
//     // printf("itgt = %i, statictgt = %5e +- %5e\n, const = %i", itgt, fStaticTgtsClones[0][itgt]->getVal(), fStaticTgtsClones[0][itgt]->getError(),int(fStaticTgtsClones[0][itgt]->getAttribute("Constant")));
//   }  
  
  //// printf("finalval = %5f\n",static_cast<RooRealVar*>(fExtVars.at(1))->getVal());

  RooAbsReal::setEvalErrorLoggingMode(RooAbsReal::PrintErrors);
  
  return nllval;
      
      
}


double RooHybridBDTAutoPdf::EvalLossAvg() {
 
  //double nllval = 0.;
  
  RooArgSet parmset(fParmVars);
	
  
  
  TH1D *hmass = new TH1D("hmass","",400,0.,100.);

    
  
  //std::vector<double> clsvals(fData.size());  
  
  for (std::vector<HybridGBREvent*>::const_iterator it = fEvts.begin(); it!=fEvts.end(); ++it) {
    
    int evcls = (*it)->Class(); 
    if (evcls!=0) continue;
    
    double weight = (*it)->Weight();    
    double pdfval = (*it)->PdfVal();
    
    double expF = exp((*it)->Target(0));
    double mpdfval = (1.0-expF)*pdfval;
    
    double mval = (*it)->Var(fCondVars.getSize() + 0);
    
    hmass->Fill(mval,weight*mpdfval);
    
    //clsvals[evcls] += clsweights[evcls]*weight*static_cast<RooAbsReal*>(fStaticPdfs.at(evcls))->getValV(&parmset);
        
    //nllval += -weight*log(pdfval); 
         
  }  
  
  
  hmass->Scale(1.0/hmass->GetSumOfWeights());
  
  double nllval = 0.;
  for (std::vector<HybridGBREvent*>::const_iterator it = fEvts.begin(); it!=fEvts.end(); ++it) {
    int evcls = (*it)->Class(); 
    if (evcls!=0) continue;
    
    double weight = (*it)->Weight();     
    
    double mval = (*it)->Var(fCondVars.getSize() + 0);

    
    double pdfval = hmass->GetBinContent(hmass->GetXaxis()->FindFixBin(mval));
    nllval += -weight*log(pdfval);
  }
  
  //delete hmass;
  
  new TCanvas;
  hmass->Draw();
  
  return nllval;

  
  
//   for (std::vector<HybridGBREvent*>::const_iterator it = fEvts.begin(); it!=fEvts.end(); ++it) {
//     
// 
//     for (int ivar=0; ivar<fParmVars.getSize(); ++ivar) {
//       static_cast<RooRealVar*>(fParmVars.at(ivar))->setVal((*it)->Var(fCondVars.getSize() + ivar));
//     }
//         
//     int evcls = (*it)->Class(); 
//     double weight = (*it)->Weight();    
//     
//     double pdfval = 0.;
//     
//     for (std::vector<HybridGBREvent*>::const_iterator jt = fEvts.begin(); jt!=fEvts.end(); ++jt) { 
//       int jcls = (*jt)->Class();
//       
//       if (jcls!=evcls) continue;
//       
//       int jweight = (*jt)->Weight();
//       
//       for (int ivar=0; ivar<fCondVars.getSize(); ++ivar) {
// 	static_cast<RooRealVar*>(fCondVars.at(ivar))->setVal((*jt)->Var(ivar));
//       }    
//       
//       for (int itgt=0; itgt<fNTargets; ++itgt) {
// 	static_cast<RooRealVar*>(fStaticTgts.at(itgt))->setVal((*jt)->Target(itgt));
//       }
//       
//       pdfval += clsweights[evcls]*jweight*static_cast<RooAbsReal*>(fStaticPdfs.at(evcls))->getValV(&parmset);
//       
//     }
//    
//     nllval += -weight*log(pdfval); 
//          
//   }
  
//   int infunc = fFullFuncs.getSize()-1;
//   nllval += static_cast<RooAbsReal*>(fFullFuncs.at(infunc))->getVal();
//   
//   double sumweight = 0.;
//   double sumnll = 0.;
//   double sumnllsq = 0.;
//   for (std::vector<HybridGBREvent*>::const_iterator it = fEvts.begin(); it!=fEvts.end(); ++it) {
//     double weight = (*it)->Weight();
//     double pdfval = (*it)->PdfVal();
//     
//     double modnll = nllval + weight*log(pdfval);
//     
//     sumweight += weight;
//     sumnll += modnll;
//     sumnllsq += modnll*modnll;
//     
// 
//   }
//   
//   
//   return sqrt(sumnllsq/sumweight - sumnll*sumnll/sumweight/sumweight);
      
      
}









TMatrixD RooHybridBDTAutoPdf::vmultT(const TVectorD &v, const TVectorD &vT) const
{
  int nrows = v.GetNrows();
  int ncols = vT.GetNrows();
  TMatrixD m(nrows,ncols);
  
  for (int irow=0; irow<nrows; ++irow) {
    for (int icol=0; icol<ncols; ++icol) {
      m(irow,icol) = v(irow)*vT(icol);
    }
  }
  
  return m;
  
}
  
double RooHybridBDTAutoPdf::vmult(const TVectorD &vT, const TVectorD &v) const
{
 
  assert(vT.GetNrows()==v.GetNrows());
  
  double val = 0.;
  for (int iel=0; iel<vT.GetNrows(); ++iel) {
    val += vT(iel)*v(iel);
  }
  
  return val;
  
}

void RooHybridBDTAutoPdf::FillDerivatives() {
 
  

//  int nextvars = fExtVars.getSize();
  
  int nparms = fFullParms.getSize();
  
  //// printf("nextvars = %i, nparms = %i\n",nextvars,nparms);
  
//  double dr = 0.;
  
//  int idxglobal = 0;
//  int idxlocal = nextvars;
  //int idxn0 = idxglobal + fExtVars.index(fN0);
  
  //// printf("idxn0 = %i, fN0 = %5f\n", idxn0,fN0->getVal());

  double nmc = 0.;
  double nmcalt = 0.;
  
    
  int msize = nparms;
  
  
  //double lambda = 1.0;
//    double nll = 0.;
  fNLLVal = 0.; 



  RooArgSet parmset(fParmVars);
	
  //// printf("begin loop\n");
  
  for (std::vector<HybridGBREvent*>::const_iterator it = fEvts.begin(); it!=fEvts.end(); ++it) {

    int ievt(it-fEvts.begin());
    
    TVectorD &dL = fGradients[ievt];
    TMatrixD &d2L = fHessians[ievt];
    
    //int termidx = (*it)->CurrentNode();
    //int idxlocal = nextvars + fNTargets*termidx;
    
    
    for (int ivar=0; ivar<fCondVars.getSize(); ++ivar) {
      static_cast<RooRealVar*>(fCondVars.find(*fCondVars.at(ivar)))->setVal((*it)->Var(ivar));
    }
    for (int ivar=0; ivar<fParmVars.getSize(); ++ivar) {
      static_cast<RooRealVar*>(fParmVars.find(*fParmVars.at(ivar)))->setVal((*it)->Var(fCondVars.getSize() + ivar));
    }    

    for (int itgt=0; itgt<fNTargets; ++itgt) {
      static_cast<RooRealVar*>(fStaticTgts.at(itgt))->setVal((*it)->Target(itgt));
    }
    
    int evcls = (*it)->Class(); 

    //nll value for minos constraint
    fNLLVal += -log(static_cast<RooAbsReal*>(fStaticPdfs.at(evcls))->getValV(&parmset));
    
    
    //derivatives for pdfs
    //for (int ivar=0; ivar<fFullParms.getSize(); ++ivar) {
    for (unsigned int iidx=0; iidx<fOuterIndices[evcls].size(); ++iidx) {
      
      int ivar = fOuterIndices[evcls][iidx];
      int idrv = ivar;
//      int itgt = ivar - fExtVars.getSize();
//      int iel;
//      if (itgt>=0) iel = idxlocal + itgt;
//      else iel = idxglobal + ivar;
      
      double drvi = Derivative1(static_cast<RooAbsReal*>(fStaticPdfs.at(evcls)),static_cast<RooRealVar*>(fFullParms.at(idrv)),&parmset,1e-3*static_cast<RooRealVar*>(fFullParms.at(idrv))->getError());
      //double drvialt = fDerivatives[evcls][idrv]->getVal();
      
      //// printf("drv%i = %5f, alt = %5f\n",idrv,drvi,drvialt);
      
      
      //dL[iel] += -fDerivatives[evcls][idrv]->getVal()/static_cast<RooAbsReal*>(fStaticPdfs.at(evcls))->getValV(&parmset);
      dL[ivar] = -drvi/static_cast<RooAbsReal*>(fStaticPdfs.at(evcls))->getValV(&parmset);
      
      
// 	if (ivar==0 && std::abs(fDerivatives[evcls][idrv]->getVal())>0.01) // printf("derivative for class %i, drv %i = %5f\n",evcls,idrv,fDerivatives[evcls][idrv]->getVal());
// 	if (static_cast<RooAbsReal*>(fStaticPdfs.at(evcls))->getValV(&parmset)<1e-4) // printf("pdf %i val = %5f\n",evcls,static_cast<RooAbsReal*>(fStaticPdfs.at(evcls))->getValV(&parmset));
      //if (ivar==0) // printf("dL[0] = %5f, fN0 = %5f\n", dL[0],fN0->getVal());
      //for (int jvar=0; jvar<fFullParms.getSize(); ++jvar) {
	
      for (unsigned int jidx=iidx; jidx<fOuterIndices[evcls].size(); ++jidx) {
	
	int jvar = fOuterIndices[evcls][jidx];
	int jdrv = jvar;
//	int jtgt = jvar - fExtVars.getSize();
	//int jel;
	//if (jtgt>=0) jel = idxlocal + jtgt;
	//else jel = idxglobal + jvar;
	
	double drvj; 
	double drv2ij;
	if (idrv==jdrv) {
	  drvj = drvi;
	  drv2ij = Derivative2(static_cast<RooAbsReal*>(fStaticPdfs.at(evcls)),static_cast<RooRealVar*>(fFullParms.at(idrv)),&parmset,1e-3*static_cast<RooRealVar*>(fFullParms.at(idrv))->getError());
	}
	else {
	  drvj = Derivative1(static_cast<RooAbsReal*>(fStaticPdfs.at(evcls)),static_cast<RooRealVar*>(fFullParms.at(jdrv)),&parmset,1e-3*static_cast<RooRealVar*>(fFullParms.at(jdrv))->getError());
	  drv2ij = Derivative2(static_cast<RooAbsReal*>(fStaticPdfs.at(evcls)),static_cast<RooRealVar*>(fFullParms.at(idrv)),static_cast<RooRealVar*>(fFullParms.at(jdrv)),&parmset,1e-3*static_cast<RooRealVar*>(fFullParms.at(idrv))->getError(),1e-3*static_cast<RooRealVar*>(fFullParms.at(jdrv))->getError());
	}
	
	d2L[ivar][jvar] = -drv2ij/static_cast<RooAbsReal*>(fStaticPdfs.at(evcls))->getValV(&parmset) + drvi*drvj/static_cast<RooAbsReal*>(fStaticPdfs.at(evcls))->getValV(&parmset)/static_cast<RooAbsReal*>(fStaticPdfs.at(evcls))->getValV(&parmset);	  
	
	//d2L[iel][jel] += -f2Derivatives[evcls][idrv][jdrv]->getVal()/static_cast<RooAbsReal*>(fStaticPdfs.at(evcls))->getValV(&parmset) + fDerivatives[evcls][idrv]->getVal()*fDerivatives[evcls][jdrv]->getVal()/static_cast<RooAbsReal*>(fStaticPdfs.at(evcls))->getValV(&parmset)/static_cast<RooAbsReal*>(fStaticPdfs.at(evcls))->getValV(&parmset);

	  
	  
      }
    }
    
    for (int i=0; i<msize; ++i) {
      for (int j=0; j<i; ++j) {
	//assert(d2L(i,j)==0.);
	d2L(i,j) = d2L(j,i);
      }
    }  
    
    //normalization terms
    if (evcls==0) {
      for (int ipdf=1; ipdf<fStaticPdfs.getSize(); ++ipdf) {
	int itgt = ipdf - 1;
//	int iel = idxlocal + itgt;
	double nexpF = exp(-(*it)->Target(itgt));
	fNLLVal += nexpF/fN0Obs;
	nmc += nexpF/fN0Obs;
	
// 	dL[iel] += -nexpF/fN0Obs;
// 	d2L[iel][iel] += nexpF/fN0Obs;
      }
    }
    else {
      int itgt = evcls - 1;
//      int iel = idxlocal + itgt;
      fNLLVal += (*it)->Target(itgt);
      double expF = exp((*it)->Target(itgt));
      nmcalt += expF;
      
      //dL[iel] += 1.0;
    }

  }
  
  int infunc = fFullFuncs.getSize()-1;
  fNLLVal += static_cast<RooAbsReal*>(fFullFuncs.at(infunc))->getVal() - fN0Obs*log(static_cast<RooAbsReal*>(fFullFuncs.at(infunc))->getVal());  
  
  
}


void RooHybridBDTAutoPdf::FitResponses(int selvar = -1) {

  // printf("fit responses, selvar = %i\n",selvar);

  //int seltgt = selvar - fExtVars.getSize();
  
  int nextvars = fExtVars.getSize();
    
  //// printf("build indexes\n");
  int nparms = nextvars;
  std::vector<int> localidxs(fFuncs.getSize());
  for (int ifunc=0; ifunc<fFuncs.getSize(); ++ifunc) {
    localidxs[ifunc] = nparms;
    HybridGBRForestFlex *forest = static_cast<RooGBRFunctionFlex*>(fFuncs.at(ifunc))->Forest();
        
    if (forest->Trees().size() && forest->Trees().back().Responses().size()) {
      nparms += forest->Trees().back().Responses().size();
    }
    else {
      nparms += 1;    
    }
  }  
  //// printf("done build indexes\n");
  

  //  double nmc = 0.;
  //double nmcalt = 0.;
  int idxglobal = 0;
  
  
      
  int msize = nparms;
  
  
  bool usematrix = (msize<8000);
  
  TMatrixDSym d2L;
  if (usematrix) {
    // printf("allocating matrix of size %i\n",msize);
    d2L.ResizeTo(msize,msize);
  }
  //std::map<std::pair<int,int>,double> d2Lmap;
  TVectorD dL(msize);
  //TVectorD d2Lv(msize);

  RooArgSet parmset(fParmVars);
	
  //// printf("begin loop\n");
  
  //std::vector<TMatrixDSym> d2Ls(fNThreads,TMatrixDSym(msize));
  //std::vector<std::map<std::pair<int,int>,double> > d2Lmaps(fNThreads);
  std::vector<TVectorD> dLs(fNThreads,TVectorD(msize));
 //std::vector<TVectorD> d2Lvs(fNThreads,TVectorD(msize));
  
  
  //map of function indices
  std::vector<int> funcidxs(fStaticTgts.getSize());
  for (int itgt=0; itgt<fStaticTgts.getSize(); ++itgt) {
    funcidxs[itgt] = fFuncs.index(static_cast<RooGBRTargetFlex*>(fTgtVars.at(itgt))->Func());
  }  
  
  std::vector<double> nmcs(fNThreads);
  
  //// printf("FitResponses start loop\n");
  
  #pragma omp parallel for
  for (unsigned int iev=0; iev<fEvts.size(); ++iev) {

    //int ithread = omp_get_thread_num();
    //int ithread = 0;
      
//    int termidx = fEvts[iev]->CurrentNode();
//    int idxlocal = nextvars + fNTargets*termidx;
    
    int evcls = fEvts[iev]->Class(); 
    double weight = fEvts[iev]->Weight();
    
    //double pdfval = static_cast<RooAbsReal*>(fStaticPdfs.at(evcls))->getValV(&parmset);
    double pdfval = fEvts[iev]->PdfVal();
    
    //double invpdf = 1.0/pdfval;
    double invpdf = vdt::fast_inv(pdfval);
    //double invpdf = vdt::fast_invf(pdfval);
    double invpdfsq = invpdf*invpdf;
    
    //double fval = fEvts[iev]->Target(0);
    //double flim = 0. + 0.5*(1.0-0.)*(sin(fval)+1.0);
    //double flim = 0. + 0.5*(1.0-0.)*(atan(fval)+1.0)*2.0/TMath::Pi();
    
    //double flim = TMath::Pi()/2.0 + atan(fval)/TMath::Pi();
    //double flim = 0.5 + atan(fval)/TMath::Pi();
    
    //double flim = 0. + 0.5*(1.0-0.)*(tanh(fval)+1.0);
    //double flim = 0.5 + atan(fval)/TMath::Pi();    
    //nmcs[ithread] += 1.0 - fval;
    //nmcs[ithread] += weight*1.0/(1.0+exp(fval));
    //nmcs[ithread] += weight*(1.0-flim);
    //nmcs[ithread] += weight*(1.0-fval);
      
    bool useseconddrv = fOuterIndices[evcls].size()==1;
    
    for (unsigned int iidx=0; iidx<fOuterIndices[evcls].size(); ++iidx) {
	
      int ivar = fOuterIndices[evcls][iidx];
      if (selvar>=0 && ivar!=selvar) continue;
      
      int idrv = ivar;
      int itgt = ivar - fExtVars.getSize(); 
      int ifunc = funcidxs[itgt];
      int iel;
      if (itgt>=0) iel = localidxs[ifunc] + fEvts[iev]->CurrentNode(itgt);
      else iel = idxglobal + ivar;

      
      double drvi = fEvts[iev]->Derivative(idrv);
      
      //dLs[ithread][iel] += -weight*drvi*invpdf;
      
      double dval = -weight*drvi*invpdf;
      #pragma omp atomic
      dL(iel) += dval;
      
      //double drv2i = fEvts[iev]->Derivative2(idrv);        
      
     // d2Lvs[ithread][iel] += -weight*drv2i*invpdf + weight*drvi*drvi*invpdfsq;
      //d2Lvs[ithread][iel] += weight*drvi*drvi*invpdfsq;
      
//       if (!std::isnormal(d2Lvs[ithread][iel]) && d2Lvs[ithread][iel]!=0.) {
//         // printf("d2LV = %5f, weight = %5e,drvi = %5e, drv2i = %5e, invpdf = %5e, invpdfsq = %5e\n",d2Lvs[ithread][iel], weight,drv2i,drvi,invpdf,invpdfsq);
//       }
      
      if (!usematrix) continue;
      
//      double drv2i = fEvts[iev]->Derivative2(idrv);  
      
      //double drv2i = fEvts[iev]->Derivative2(idrv);
      //d2Lmaps[ithread][std::pair<int,int>(iel,iel)] += -weight*drv2i*invpdf;
      
      //d2Ls[ithread][iel][iel] += -weight*drv2i*invpdf + weight*drvi*drvi*invpdfsq;
      //d2Ls[ithread][iel][iel] += -weight*drv2i*invpdf;
      

//       double drv2i = fEvts[iev]->Derivative2(idrv);
//       
//       double d2approxi = weight*drvi*drvi*invpdfsq;
//       double d2i = d2approxi - weight*drv2i*invpdf;
//       
//       if (std::abs(d2i)<(std::abs(1e-3*d2approxi))) continue;
      
      
      for (unsigned int jidx=iidx; jidx<fOuterIndices[evcls].size(); ++jidx) {
      
	
	int jvar = fOuterIndices[evcls][jidx];
	int jdrv = jvar;
	int jtgt = jvar - fExtVars.getSize();
        int jfunc = funcidxs[jtgt];     
	int jel;      
	if (jtgt>=0) jel = localidxs[jfunc] + fEvts[iev]->CurrentNode(jtgt);
	else jel = idxglobal + jvar;    
	
	
	
	double drvj = fEvts[iev]->Derivative(jdrv);
        //double drv2ij = fEvts[iev]->ParmMatrix()(idrv,jdrv);
        
        //double updval = d2L(iel,jel) + weight*drvi*drvj*invpdfsq;
        
        //if (fEvts[iev]->Derivative2(idrv)==0. || fEvts[iev]->Derivative2(jdrv)==0.) continue;
        
/*        double drv2j = fEvts[iev]->Derivative2(idrv);
        
        
        double d2approxj = weight*drvj*drvj*invpdfsq;
        double d2j = d2approxj - weight*drv2j*invpdf;
        
        if (std::abs(d2j)<(std::abs(1e-3*d2approxj))) continue;  */    
        
        double d2val = weight*drvi*drvj*invpdfsq;
        if (useseconddrv) {
          double drv2j = fEvts[iev]->Derivative2(idrv);
          d2val += -weight*drv2j*invpdf;
        }
        
        
        
//         if (idrv==jdrv) {
//           double drv2ij = fEvts[iev]->Derivative2(jdrv);
//           d2val += -weight*drv2ij*invpdf;
//         }
//         else {
//           continue;
//         }
        
        int ielf = iel;
        int jelf = jel;
        if (ielf>jelf) {
          ielf = jelf;
          jelf = ielf;
        }
        
        
        #pragma omp atomic
        d2L(ielf,jelf) += d2val;
        //d2L(iel,jel) += weight*drvi*drvj*invpdfsq;
        
        //d2Lmaps[ithread][std::pair<int,int>(iel,jel)] += weight*drvi*drvj*invpdfsq;
	//d2Ls[ithread][iel][jel] += weight*drvi*drvj*invpdfsq;
        //d2Ls[ithread][iel][jel] +=  -weight*drv2ij*invpdf + weight*drvi*drvj*invpdfsq;

      
      }
	    
    }

  }
  
//   for (int ithread=0; ithread<fNThreads; ++ithread) {
//     dL += dLs[ithread];
//     //d2Lv += d2Lvs[ithread];
//     //d2L += d2Ls[ithread];
//     nmc += nmcs[ithread];
//     
// //     for (int i=0; i<msize; ++i) {
// //       // printf("ithread = %i, i = %i, dLs[ithread] = %5e, d2Lvs[ithread] = %5e, d2Lv = %5e\n",ithread,i,dLs[ithread][i],d2Lvs[ithread][i],d2Lv[i]);
// //     }
//     
//     //if (!usematrix) continue;
// 
//     
// //     for (std::map<std::pair<int, int>, double>::const_iterator it=d2Lmaps[ithread].begin(); it!=d2Lmaps[ithread].end(); ++it) {
// //       d2Lmap[it->first] += it->second;
// //     }
//   }

  
  int infunc = fFullFuncs.getSize()-1;
   
  for (unsigned int iidx=0; iidx<fOuterIndices[infunc].size(); ++iidx) {
	int ivar = fOuterIndices[infunc][iidx];
        if (selvar>=0 && ivar!=selvar) continue;
        
	int idrv = ivar;      
	int iel = idxglobal + ivar;

	double drvi = Derivative1(static_cast<RooAbsReal*>(fFullFuncs.at(infunc)),static_cast<RooRealVar*>(fFullParms.at(idrv)),&parmset,1e-3*static_cast<RooRealVar*>(fFullParms.at(idrv))->getError());
	
	dL[iel] += drvi;

	//double drv2i = Derivative2(static_cast<RooAbsReal*>(fFullFuncs.at(infunc)),static_cast<RooRealVar*>(fFullParms.at(idrv)),&parmset,1e-3*static_cast<RooRealVar*>(fFullParms.at(idrv))->getError());
	
	//d2Lv[iel] += drv2i;
        
        if (!usematrix) continue;        
	


	
    for (unsigned int jidx=iidx; jidx<fOuterIndices[infunc].size(); ++jidx) {
	
	int jvar = fOuterIndices[infunc][jidx];
	int jdrv = jvar; 
	int jel = idxglobal + jvar;
	
	//double drvj;
	double drv2ij;
	if (idrv==jdrv) {
	  //drvj = drvi;
	  drv2ij = Derivative2(static_cast<RooAbsReal*>(fFullFuncs.at(infunc)),static_cast<RooRealVar*>(fFullParms.at(idrv)),&parmset,1e-3*static_cast<RooRealVar*>(fFullParms.at(idrv))->getError());
	}
	else {
	  //drvj = Derivative1(static_cast<RooAbsReal*>(fFullFuncs.at(infunc)),static_cast<RooRealVar*>(fFullParms.at(jdrv)),&parmset,1e-3*static_cast<RooRealVar*>(fFullParms.at(jdrv))->getError());
	  drv2ij = Derivative2(static_cast<RooAbsReal*>(fFullFuncs.at(infunc)),static_cast<RooRealVar*>(fFullParms.at(idrv)),static_cast<RooRealVar*>(fFullParms.at(jdrv)),&parmset,1e-3*static_cast<RooRealVar*>(fFullParms.at(idrv))->getError(),1e-3*static_cast<RooRealVar*>(fFullParms.at(jdrv))->getError());
	}
	
	d2L[iel][jel] += drv2ij;
	//d2Lmap[std::pair<int,int>(iel,jel)] += drv2ij;
	
    }
	
  }
  
 // // printf("FitRespones end loop\n");
  
//  double sumweight = 0.;
  
  //symmetrize matrix
  if (usematrix) {
    #pragma omp parallel for
    for (int i=0; i<msize; ++i) {
      for (int j=0; j<i; ++j) {
	assert(d2L(i,j)==0.);
	d2L(i,j) = d2L(j,i);
      }
    }
  }

  bool solved = false;
  //bool gradient = false;
  
  TVectorD dpar(msize);
  TVectorD dparm(msize);    
  TVectorD dparg = -1.0*dL;
  
//solve
  if (usematrix) {
    TVectorD dLc(dL);
    // printf("start matrix decomposition\n");  
    //TDecompLU dbk(d2L);
    TDecompBK dbk(d2L);
    //TDecompChol dbk(d2L);
    solved = dbk.Solve(dLc);
    // printf("done matrix decomposition\n");  
    dparm = -1.0*dLc;
  }
  
//   double dllrm = 0.;
//   double dllrg = 0.;
  
  double step = 0;
  double stepm = 0;
  double stepg = 0;
  
  //solved = false;
  
  if (solved) {
    for (int i=0; i<msize; ++i) {
      //protect against inf/nan elements and set them to zero
      if (!std::isnormal(dparm(i))) dparm(i) = 0.;
      
      //protect against elements which run against the gradient and invert them
      if ( (dL(i)*dparm(i))>0. ) dparm(i) = -dparm(i);
      //if ( (dL(i)*dparm(i))>0. ) dparm(i) = 0.;
    }
    
    stepm = 1.0;
    
    step = fShrinkage*stepm;
    dpar = dparm;    
  }
  else {
    dparg = -1.0*dL;
        
    double deltaL = 0;
    for (int i=0; i<msize; ++i) {
      deltaL += dL[i]*dparg[i];
    }   
   
    stepg = 0.;
    double drvstep = 0.1/deltaL;
   
    // printf("fallback to gradient descent: deltaL = %5e, initial stepg = %5e, drvstep = %5e\n",deltaL, stepg,drvstep);
    
    double nllval = fNLLVal;
    
    while (stepg==0.) {
      double upnllval = EvalLoss(drvstep,dparg);
      if (!std::isnormal(upnllval)) {
        drvstep /= 2.0;
        continue;
      }
        
      double downnllval = EvalLoss(-drvstep,dparg);
      if (!std::isnormal(downnllval)) {
        drvstep /= 2.0;
        continue;
      }    
      
      double nlldrv2 = (upnllval + downnllval - 2.0*nllval)/drvstep/drvstep;      

      stepg = -deltaL/nlldrv2;

      if (nlldrv2<0.) {
        stepg = -stepg;
      }      
      
      dpar = dparg;
      step = fShrinkage*stepg;       
    
    }
    
  }
    
  double nllval = fNLLVal;
  int stepiter = 0;
  do {

//     if (stepiter==1) {
//       // printf("fallback to gradient descent\n");
//       dpar = -1.0*dL;
//       step = 1e-8;      
//     }
    
    nllval = EvalLoss(step,dpar);
    // printf("step = %5f, nllval = %5f, fNLLVal = %5f\n",step,nllval,fNLLVal);
    //if (std::isnormal(nllval)) break;

    
    if ( std::isnormal(nllval) && (nllval-fNLLVal)<1e-3 ) {
      break;
    }
    else {
      step /= 2.0;
    }
    
//     if ( (nllval-fNLLVal)<1e-3 ) {
//       break; 
//     }
//     else {
//       step /= 2.0;
//     }
    ++stepiter;
  }
  while (stepiter<50);

//   if (!((nllval-fNLLVal)<1e-3)) {
//     step = 0.;
//     nllval = EvalLoss(forest, step,dpar);
//   }
  
  fNLLVal = nllval;
  
  
//   for (int itgt=0; itgt<fNTargets; ++itgt) {
//     double currentval = (*it)->Target(itgt);
//     double newval = currentval + step*dpar(ivar);
//     double minval = static_cast<RooRealVar*>(fStaticTgts.at(itgt))->getMin();
//     double maxval = static_cast<RooRealVar*>(fStaticTgts.at(itgt))->getMax();
//     if (newval<minval) {
//       dpar(ivar) = (minval-currentval)/fShrinkage;
//     }
//     if (newval>maxval) {
//       dpar(ivar) = (maxval-currentval)/fShrinkage;
//     }
//   }
  
  
  for (int ivar=0; ivar<fExtVars.getSize(); ++ivar) {
    if (selvar>=0 && ivar!=selvar) continue;
    static_cast<RooRealVar*>(fExtVars.at(ivar))->setVal(static_cast<RooRealVar*>(fExtVars.at(ivar))->getVal() + step*dpar(ivar));   
  }
  
  for (int ithread=0; ithread<fNThreads; ++ithread) {
    for (int ivar=0; ivar<fExtVars.getSize(); ++ivar) {
      RooRealVar *cloneparm = static_cast<RooRealVar*>(fExtVarsClones[ithread].at(ivar));
      cloneparm->setVal(static_cast<RooRealVar*>(fExtVars.at(ivar))->getVal());
      cloneparm->setError(static_cast<RooRealVar*>(fExtVars.at(ivar))->getError());
    }
  }    
    
  for (int ifunc=0; ifunc<fFuncs.getSize(); ++ifunc) {
    HybridGBRForestFlex *forest = static_cast<RooGBRFunctionFlex*>(fFuncs.at(ifunc))->Forest();
    int nterm = forest->Trees().back().Responses().size();
    for (int iterm=0; iterm<nterm; ++iterm) {
      forest->Trees().back().Responses()[iterm] = forest->Trees().back().Responses()[iterm] + step*dpar(localidxs[ifunc] + iterm);
    }
  }
  
  
  for (std::vector<HybridGBREvent*>::const_iterator it = fEvts.begin(); it!=fEvts.end(); ++it) {
    for (int itgt=0; itgt<fNTargets; ++itgt) {
      int ifunc = funcidxs[itgt];   
      int termidx = (*it)->CurrentNode(itgt);
      (*it)->SetTarget(itgt,(*it)->Target(itgt)+step*dpar(localidxs[ifunc] + termidx));
    }
  }  

  double dldr = -Derivative1(fN0,fR,&parmset,1e-3*fR->getError());
  fdLdR = dldr;
  
  //  double etermval = fN0->getVal();
  

  
  // printf("r = %5f +- %5f, dldr = %5f, dL[0] = %5e,, n0 = %5f, ndobs = %5f, fNLLVal = %5f, bnll = %5f, nmc = %5f, nmcalt = %5f\n",fR->getVal(),fR->getError(),dldr,dL[0],fN0->getVal(),fN0Obs,fNLLVal,fNLLVal-etermval,nmc,nmcalt);
    
  
  
}

void RooHybridBDTAutoPdf::RecomputeTargets() {
 
  std::vector<std::vector<float> > evalvectors(fNThreads, std::vector<float>(fCondVars.getSize()));
  
  #pragma omp parallel for
  for (unsigned int ievt=0; ievt<fEvts.size(); ++ievt) {
    int ithread =  omp_get_thread_num();
    
    for (int ivar=0; ivar<fCondVars.getSize(); ++ivar) {
      evalvectors[ithread][ivar] = fEvts.at(ievt)->Var(ivar);
    }
    
    for (int itgt=0; itgt<fStaticTgts.getSize(); ++itgt) {
      RooGBRTargetFlex *target = static_cast<RooGBRTargetFlex*>(fTgtVars.at(itgt));
      fEvts[ievt]->SetTarget(itgt, target->Forest()->GetResponse(&evalvectors[ithread][0]));
      //fEvts[ievt]->SetTarget(itgt, fFunc->Forest()->GetResponse(&evalvectors[ithread][0],treetgt));
    }
    
  }
   
  
}

void RooHybridBDTAutoPdf::GradientMinos() {
 
//     SetMinCutSignificance(1.0);
//     TrainForest(1e6,true);
//     return; 

#if 0
    //save initial state
    HybridGBRForestD *origforest = new HybridGBRForestD(*fFunc->Forest());
    std::vector<double> extvals(fExtVars.getSize());
    for (int ivar=0; ivar<fExtVars.getSize(); ++ivar) {
      extvals[ivar] = static_cast<RooRealVar*>(fExtVars.at(ivar))->getVal();
    }

  
    double rval = fR->getVal(); 
    double rerr = fR->getError();
        
    //double errscale = 1e-2;
    double errscale = 1e-3;
    
    int nstep = 3;
    double stepsize = 1.2*rerr*(1.0+errscale*errscale)/(double)nstep;
    
    std::vector<std::pair<double,double> > drvvals;
    drvvals.push_back(std::pair<double,double>(rval,0.));
    
    double coeffval = 0.5/pow(errscale*rerr,2);
    double errval = errscale*rerr;
    
    fConstraintCoeff->setVal(coeffval);  
    fR->setError(errval);     
    
    //upper uncertainty
    double rstep = rval;
    //fConstraintVal->setVal(rstep);
    //fConstraintCoeff->setVal(0.);  
    //TrainForest(0,false);
    //fConstraintCoeff->setVal(coeffval);  
    //fR->setError(errval);
    //TrainForest(1e6,true);
/*    TrainForest(e6,false);
    double drvval = -2.0*fConstraintCoeff->getVal()*(fR->getVal()-fConstraintVal->getVal());
    drvvals.push_back(std::pair<double,double>(fR->getVal(),drvval));*/    
    
    for (int istep=0; istep<nstep; ++istep) {
      rstep += stepsize;
      fConstraintVal->setVal(rstep);
      //fR->setVal(rstep);
//       fConstraintCoeff->setVal(0.);  
//       TrainForest(0,false);
//       fConstraintCoeff->setVal(coeffval);  
//       fR->setError(errval); 
//       TrainForest(1e6,true);
      TrainForest(1e6,false);
      double drvval = -2.0*fConstraintCoeff->getVal()*(fR->getVal()-fConstraintVal->getVal());
      drvvals.push_back(std::pair<double,double>(fR->getVal(),drvval));
    }   
    
//     //restore initial state
//     fFunc->SetForest(new HybridGBRForestD(*origforest));
//     for (int ivar=0; ivar<fExtVars.getSize(); ++ivar) {
//       static_cast<RooRealVar*>(fExtVars.at(ivar))->setVal(extvals[ivar]);
//     }
//     RecomputeTargets();
    
 
    
    //lower uncertainty
    rstep = rval;
    for (int istep=0; istep<nstep; ++istep) {
      rstep -= stepsize;
      fConstraintVal->setVal(rstep);
      //fR->setVal(rstep);
//       fConstraintCoeff->setVal(0.);  
//       TrainForest(0,false);
//       fConstraintCoeff->setVal(coeffval);  
//       fR->setError(errval);
//       TrainForest(1e6,true);
      TrainForest(1e6,false);
      double drvval = -2.0*fConstraintCoeff->getVal()*(fR->getVal()-fConstraintVal->getVal());
      drvvals.push_back(std::pair<double,double>(fR->getVal(),drvval));
    }    
    
    std::sort(drvvals.begin(),drvvals.end());
    
    if (fDrvGraph) delete fDrvGraph;
    fDrvGraph = new TGraph(drvvals.size());		
    for (unsigned int istep=0; istep<drvvals.size(); ++istep) {
      fDrvGraph->SetPoint(istep,drvvals[istep].first,2.0*drvvals[istep].second);
    }
    
    int nstepint = 1001;
    double stepsizeint = (drvvals.back().first-drvvals.front().first)/(double)(nstepint-1);
    rstep = drvvals.front().first;
    if (fDrvGraphSmooth) delete fDrvGraphSmooth;
    fDrvGraphSmooth = new TGraph(nstepint);
   
    TGraph *intgraph = new TGraph(nstepint);
    double intval = 0.; 
    
    for (int istep = 0; istep<nstepint; ++istep) {

      double drvval = fDrvGraph->Eval(rstep,0,"S");
      intval += drvval*stepsizeint;
      
      double rstepint = rstep + 0.5*stepsizeint;
      
      fDrvGraphSmooth->SetPoint(istep,rstep,drvval);
      intgraph->SetPoint(istep,rstepint,intval);
      
      rstep += stepsizeint;
    }
    
    double minval = std::numeric_limits<double>::max();
    double minr = 0.;
    int minstep = 0;
    for (int istep=0; istep<nstepint; ++istep) {
      double val;
      intgraph->GetPoint(istep,rstep,val);
      if (val<minval) {
	minval = val;
	minr = rstep;
	minstep = istep;
      }
    }
    
    double rlow = 0.;
    double rhigh = 0.;
    
    for (int istep=minstep; istep>=0; --istep) {
      double val;
      intgraph->GetPoint(istep,rstep,val);
      if ( (val-minval)>1.0 ) {
	rlow = rstep;
	break;
      }
    }

    for (int istep=minstep; istep<nstepint; ++istep) {
      double val;
      intgraph->GetPoint(istep,rstep,val);
      if ( (val-minval)>1.0 ) {
	rhigh = rstep;
	break;
      }
    }
    
    if (fGraphDelta) delete fGraphDelta;
    fGraphDelta = new TGraph(nstepint);			
    for (int istep=0; istep<nstepint; ++istep) {
      double val;
      intgraph->GetPoint(istep,rstep,val);
      fGraphDelta->SetPoint(istep,rstep,(val-minval));
    }
    delete intgraph;
    
    fRMin = minr;
    fRHigh = rhigh;
    fRLow = rlow;
    
    // printf("rfit = %5f, rerr = %5f, minr = %5f, rlow = %5f, rhigh = %5f, sigmu = %5f\n",rval,rerr,minr,rlow,rhigh,0.5*(rhigh-rlow));
    // printf("nomfit: r = %5f +%5f -%5f\n",rval,rhigh-rval,rval-rlow);
    // printf("minr  : r = %5f +%5f -%5f\n",minr,rhigh-minr,minr-rlow);      
    
   FILE * fp;
   fp = fopen ("confint.txt", "w");
   // fprintf(fp,"rfit = %5f, rerr = %5f, minr = %5f, rlow = %5f, rhigh = %5f, sigmu = %5f,\n",rval,rerr,minr,rlow,rhigh,0.5*(rhigh-rlow));
   // fprintf(fp,"nomfit: r = %5f +%5f -%5f\n",rval,rhigh-rval,rval-rlow);
   // fprintf(fp,"minr  : r = %5f +%5f -%5f\n",minr,rhigh-minr,minr-rlow);    
   fclose(fp);     
    
    //restore initial state
    fFunc->SetForest(new HybridGBRForestD(*origforest));
    delete origforest;
    for (int ivar=0; ivar<fExtVars.getSize(); ++ivar) {
      static_cast<RooRealVar*>(fExtVars.at(ivar))->setVal(extvals[ivar]);
    }
    fR->setError(rerr);
    
#endif
    
    return;
  
}

void RooHybridBDTAutoPdf::fitWithMinosFast() {
  TrainForest(fNTrees);
  //TrainForest(fNTrees);
  
  //return;

  
  
  //fLambda->setVal(1.0);
  fLambda->setVal(1.);
//  double thres = fNLLVal + 0.5;
  //fNLLVal = 0.;
  //double thres = -245e3;
//  TrainForest(10*fNTrees,true,thres,fR->getVal(),fR->getError());
  //TrainForest(1,true,fNLLVal+10e3,fR->getVal(),fR->getError());
  fLambda->setVal(0.);
}


void RooHybridBDTAutoPdf::fitWithMinos() {
  
    RooRealVar &r = *fR;
    bool isconstant = r.getAttribute("Constant");
    
    
    TrainForest(fNTrees);
    double rerr = r.getError();
  
    //RooAbsReal *nll = this->createNLL(data,RooFit::ConditionalObservables(fCondVars),RooFit::Optimize(false),RooFit::Extended());  

    double r0 = r.getVal();
    double rMax = r.getMax();
    double rMin = r.getMin();
    
    r.setConstant();
    
//    int ndim = 1;
    //double delta68 = 0.5*ROOT::Math::chisquared_quantile_c(1-0.68,ndim);
    double delta68 = 0.5;
    double nll0 = fNLLVal;
    double threshold68 = nll0 + delta68; 
    
    double hi68 = findCrossing(r, threshold68, r0,   rMax, rerr);
    
    fNLLVal = nll0;
    //double hi95 = do95_ ? findCrossing(minim2, *nll, r, threshold95, std::isnan(hi68) ? r0 : hi68, std::max(rMax, std::isnan(hi68*2-r0) ? r0 : hi68*2-r0)) : r0;
    // low error 
    double lo68 = findCrossing(r, threshold68, r0,   rMin, rerr); 
    //double lo95 = do95_ ? findCrossing(minim2, *nll, r, threshold95, std::isnan(lo68) ? r0 : lo68, rMin) : r0;

    fNLLVal = nll0;
    
    r.setVal(r0);
    r.setError(rerr);
    r.setAsymError(!std::isnan(lo68) ? lo68 - r0 : 0, !std::isnan(hi68) ? hi68 - r0 : 0);    
    
    r.setConstant(isconstant);
    
    
    
    // printf("r = %5f +%5f -%5f\n",r.getVal(),std::abs(r.getAsymErrorHi()),std::abs(r.getAsymErrorLo()));
    
    //delete nll;
}
 
double RooHybridBDTAutoPdf::findCrossing(RooRealVar &r, double level, double rStart, double rBound, double rerr) {

    r.setVal(rStart);
    //double rInc = 0.1*(rBound - rStart);
    double rInc = rerr*(rBound-rStart)/std::abs(rBound-rStart);
    
    //    int verbose = 9;
    
    double here = fNLLVal;
    do {
        rStart += rInc;
        if (rInc*(rStart - rBound) > 0) { // went beyond bounds
            rStart -= rInc;
            rInc    = 0.5*(rBound-rStart);
        }
        r.setVal(rStart);
        //nll.clearEvalErrorLog(); nll.getVal();
        TrainForest(fNTrees);
        double there = here;
        here = fNLLVal;
	//        if (verbose > 0) { // printf("%f    %+.5f  %+.5f    %f\n", rStart, level-here, level-there, rInc); fflush(stdout); }
        if ( fabs(here - level) < 4*1e-4 ) {
            // set to the right point with interpolation
            r.setVal(rStart + (level-here)*(level-there)/(here-there));
            return r.getVal();
        } else if (here > level) {
            // I'm above the level that I wanted, this means I stepped too long
            // First I take back all the step
            rStart -= rInc; 
            // Then I try to figure out a better step
            if (1) {
                if (fabs(there - level) > 0.05) { // If started from far away, I still have to step carefully
                    double rIncFactor = std::max(0.2, std::min(0.7, 0.75*(level-there)/(here-there)));
                    //// printf("\t\t\t\t\tCase A1: level-there = %f,  here-there = %f,   rInc(Old) = %f,  rInFactor = %f,  rInc(New) = %f\n", level-there, here-there, rInc, rIncFactor, rInc*rIncFactor);
                    rInc *= rIncFactor;
                } else { // close enough to jump straight to target
                    double rIncFactor = std::max(0.05, std::min(0.95, 0.95*(level-there)/(here-there)));
                    //// printf("\t\t\t\t\tCase A2: level-there = %f,  here-there = %f,   rInc(Old) = %f,  rInFactor = %f,  rInc(New) = %f\n", level-there, here-there, rInc, rIncFactor, rInc*rIncFactor);
                    rInc *= rIncFactor;
                }
            } else {
                rInc *= 0.3;
            }
            //if (allpars.get() == 0) allpars.reset(nll.getParameters((const RooArgSet *)0));
            //RooArgSet oldparams(checkpoint->floatParsFinal());
            //*allpars = oldparams;
        } else if ((here-there)*(level-there) < 0 && // went wrong
                   fabs(here-there) > 0.1) {         // by more than roundoff
//            if (allpars.get() == 0) allpars.reset(nll.getParameters((const RooArgSet *)0));
            //RooArgSet oldparams(checkpoint->floatParsFinal());
            //*allpars = oldparams;
            rStart -= rInc; rInc *= 0.5;
        } else {
            // I did this step, and I'm not there yet
            if (1) {
                if (fabs(here - level) > 0.05) { // we still have to step carefully
                    if ((here-there)*(level-there) > 0) { // if we went in the right direction
                        // then optimize step size
                        double rIncFactor = std::max(0.2, std::min(2.0, 0.75*(level-there)/(here-there)));
                        //// printf("\t\t\t\t\tCase B1: level-there = %f,  here-there = %f,   rInc(Old) = %f,  rInFactor = %f,  rInc(New) = %f\n", level-there, here-there, rInc, rIncFactor, rInc*rIncFactor);
                        rInc *= rIncFactor;
                    } //else // printf("\t\t\t\t\tCase B3: level-there = %f,  here-there = %f,   rInc(Old) = %f\n", level-there, here-there, rInc);
                } else { // close enough to jump straight to target
                    double rIncFactor = std::max(0.05, std::min(4.0, 0.95*(level-there)/(here-there)));
                    //// printf("\t\t\t\t\tCase B2: level-there = %f,  here-there = %f,   rInc(Old) = %f,  rInFactor = %f,  rInc(New) = %f\n", level-there, here-there, rInc, rIncFactor, rInc*rIncFactor);
                    rInc *= rIncFactor;
                }
            } else {
                //nothing?
            }
            //checkpoint.reset(minim.save());
        }
    } while (fabs(rInc) > 1e-4*0.1*std::max(1.0,rBound-rStart));
    return r.getVal();
//     if (fabs(here - level) > 0.01) {
//         std::cout << "Error: closed range without finding crossing." << std::endl;
//         return NAN;
//     } else {
//         return r.getVal();
//     } 
  
}



double RooHybridBDTAutoPdf::Derivative1Fast(RooAbsReal *function, double currentval, RooRealVar *var, RooArgSet *nset, double step) {
 
  //// printf("var = %s, step = %5e\n",var->GetName(),step);
  
  double startval = var->getVal();
  
  var->setVal(startval+step);
  double valup = function->getValV(nset);

  var->setVal(startval-step);
  double valdown = function->getValV(nset);  
  
  double drv = (valup-valdown)*vdt::fast_inv(2.0*step);
  
  var->setVal(startval);
  
  return drv;
  
}

double RooHybridBDTAutoPdf::Derivative2Fast(RooAbsReal *function, double currentval, RooRealVar *var, RooArgSet *nset, double step) {
 
  double startval = var->getVal();
  
  double valnom = currentval;
  
  var->setVal(startval+step);
  double valup1 = function->getValV(nset);
  
  var->setVal(startval-step);
  double valdown1 = function->getValV(nset);  

  double drv1 = (valup1+valdown1-2.0*valnom)*vdt::fast_inv(step*step);
  
  var->setVal(startval);
  
  return drv1;
  
}

double RooHybridBDTAutoPdf::Derivative2Fast(RooAbsReal *function, RooRealVar *var1, RooRealVar *var2, RooArgSet *nset, double step1, double step2) {

  double startval1 = var1->getVal();
  double startval2 = var2->getVal();
    
  //double valnom = function->getValV(nset);
  
  double stepa1 = step1;
  double stepb1 = step2;
  
  var1->setVal(startval1+stepa1);
  var2->setVal(startval2+stepb1);
  double valupup1 = function->getValV(nset);

  var1->setVal(startval1+stepa1);
  var2->setVal(startval2-stepb1);
  double valupdown1 = function->getValV(nset);  
  
  var1->setVal(startval1-stepa1);
  var2->setVal(startval2+stepb1);
  double valdownup1 = function->getValV(nset);    
  
  var1->setVal(startval1-stepa1);
  var2->setVal(startval2-stepb1);
  double valdowndown1 = function->getValV(nset);    
  
  double drv1 = (valupup1+valdowndown1-valupdown1-valdownup1)*vdt::fast_inv(4.0*stepa1*stepb1);
  
  var1->setVal(startval1);
  var2->setVal(startval2);
  
  return drv1;
  
}



double RooHybridBDTAutoPdf::Derivative1(RooAbsReal *function, RooRealVar *var, RooArgSet *nset, double step) {
 
  //// printf("var = %s, step = %5e\n",var->GetName(),step);
  
  double startval = var->getVal();
  
  var->setVal(startval+step);
  double valup1 = function->getValV(nset);
  
  var->setVal(startval-step);
  double valdown1 = function->getValV(nset);  

  double drv1 = (valup1-valdown1)/(2.0*step);
  
  double step2 = 0.5*step;
  
  var->setVal(startval+step2);
  double valup2 = function->getValV(nset);
  
  var->setVal(startval-step2);
  double valdown2 = function->getValV(nset);    
  
  double drv2 = (valup2-valdown2)/(2.0*step2);
  
  double drv = (4.0*drv2 - drv1)/3.0;
  
  var->setVal(startval);
  
  return drv;
  
}

double RooHybridBDTAutoPdf::Derivative2(RooAbsReal *function, RooRealVar *var, RooArgSet *nset, double step) {
 
  double startval = var->getVal();
  
  double valnom = function->getValV(nset);
  
  var->setVal(startval+step);
  double valup1 = function->getValV(nset);
  
  var->setVal(startval-step);
  double valdown1 = function->getValV(nset);  

  double drv1 = (valup1+valdown1-2.0*valnom)/(step*step);
  
  double step2 = 0.5*step;  
  
  var->setVal(startval+step2);
  double valup2 = function->getValV(nset);
  
  var->setVal(startval-step2);
  double valdown2 = function->getValV(nset);    
  
  double drv2 = (valup2+valdown2-2.0*valnom)/(step2*step2);
  
  double drv = (4.0*drv2 - drv1)/3.0;
  
  var->setVal(startval);
  
  return drv;
  
}

double RooHybridBDTAutoPdf::Derivative2(RooAbsReal *function, RooRealVar *var1, RooRealVar *var2, RooArgSet *nset, double stepa, double stepb) {
 
  //double step = std::min(stepa,stepb);
  
  double startval1 = var1->getVal();
  double startval2 = var2->getVal();
    
  //double valnom = function->getValV(nset);
  
  double stepa1 = stepa;
  double stepb1 = stepb;
  
  var1->setVal(startval1+stepa1);
  var2->setVal(startval2+stepb1);
  double valupup1 = function->getValV(nset);

  var1->setVal(startval1+stepa1);
  var2->setVal(startval2-stepb1);
  double valupdown1 = function->getValV(nset);  
  
  var1->setVal(startval1-stepa1);
  var2->setVal(startval2+stepb1);
  double valdownup1 = function->getValV(nset);    
  
  var1->setVal(startval1-stepa1);
  var2->setVal(startval2-stepb1);
  double valdowndown1 = function->getValV(nset);    
  
  double drv1 = (valupup1+valdowndown1-valupdown1-valdownup1)/(4.0*stepa1*stepb1);
  
  
      
  
  double stepa2 = 0.5*stepa;
  double stepb2 = 0.5*stepb;  
  
  var1->setVal(startval1+stepa2);
  var2->setVal(startval2+stepb2);
  double valupup2 = function->getValV(nset);

  var1->setVal(startval1+stepa2);
  var2->setVal(startval2-stepb2);
  double valupdown2 = function->getValV(nset);  
  
  var1->setVal(startval1-stepa2);
  var2->setVal(startval2+stepb2);
  double valdownup2 = function->getValV(nset);    
  
  var1->setVal(startval1-stepa2);
  var2->setVal(startval2-stepb2);
  double valdowndown2 = function->getValV(nset);    
  
  double drv2 = (valupup2+valdowndown2-valupdown2-valdownup2)/(4.0*stepa2*stepb2);
  
  double drv = (4.0*drv2 - drv1)/3.0;
  
  var1->setVal(startval1);
  var2->setVal(startval2);
  
  return drv;
  
}
