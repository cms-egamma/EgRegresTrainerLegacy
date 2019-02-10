#include "GBRLikelihood/GBRArrayUtils.h" 
#include "GBRLikelihood/GBRMath.h"
#include <limits>
#include <algorithm>    
    
void GBRArrayUtils::InitArrays(int *__restrict__ ns, double *__restrict__ tgts, double *__restrict__ tgt2s, float *__restrict__ bsepgains, const int nbins) {
 
#if  __GNUC__>4 || (__GNUC__==4 && __GNUC_MINOR__>=7)
  ns = (int*)__builtin_assume_aligned(ns,32);
  tgts = (double*)__builtin_assume_aligned(tgts,32);
  tgt2s = (double*)__builtin_assume_aligned(tgt2s,32);
  bsepgains = (float*)__builtin_assume_aligned(bsepgains,32);
#endif  
  
  for (int ibin=0; ibin<nbins; ++ibin) {
    ns[ibin] = 0;
    tgts[ibin] = 0.;
    tgt2s[ibin] = 0.;     
    
    bsepgains[ibin] = -std::numeric_limits<float>::max();
  }
   
}
    
void GBRArrayUtils::ZeroArray(double *__restrict__ wscls, const int nbins) {
  
#if  __GNUC__>4 || (__GNUC__==4 && __GNUC_MINOR__>=7)
  wscls = (double*)__builtin_assume_aligned(wscls,32);
#endif
  
  for (int ibin=0; ibin<nbins; ++ibin) {
    wscls[ibin] = 0.;
  }
}
  
void GBRArrayUtils::MinMaxQuants(int &__restrict__ minquant, int &__restrict__ maxquant, const int *__restrict__ quants, const int nev) {
  
#if  __GNUC__>4 || (__GNUC__==4 && __GNUC_MINOR__>=7)  
  quants = (const int*)__builtin_assume_aligned(quants,32);
#endif  
  
  minquant = std::numeric_limits<int>::max();
  maxquant = 0;
   
  for (int iev = 0; iev<nev; ++iev) {
    if (quants[iev]<minquant) minquant = quants[iev];
    if (quants[iev]>maxquant) maxquant = quants[iev];
  }      
  
} 
 
void GBRArrayUtils::FillBinQuants(int *__restrict__ binquants, const unsigned int offset, const unsigned int pscale, const unsigned int nquantiles, const unsigned int nbins) {

#if  __GNUC__>4 || (__GNUC__==4 && __GNUC_MINOR__>=7)    
  binquants = (int*)__builtin_assume_aligned(binquants,32);
#endif
    
  for (unsigned int ibin=0; ibin<nbins; ++ibin) { 
    //int scaledbin
    //int quant = ((ibin+1)<<pscale) + offset - 1;
    unsigned int quant = ((ibin+1)<<pscale) + offset - 1;
    //unsigned short quant = (ibin<<pscale) + offset - 1;
    //int quant = ((ibin+1)<<pscale) + offset - 1;
    binquants[ibin] = std::min(quant, nquantiles-1);
    //binquants[ibin] = quant < nquantiles ? quant : nquantiles-1;
  }  
  
}
 
void GBRArrayUtils::FillSepGains(const double *__restrict__ sumtgts, const double *__restrict__ sumtgt2s, float *__restrict__ bsepgains, const double fulldiff, const double sumtgt, const double sumtgt2, const int nbins) {
  
#if __GNUC__>4 || (__GNUC__==4 && __GNUC_MINOR__>=7)  
  sumtgts = (const double*)__builtin_assume_aligned(sumtgts,32);
  sumtgt2s = (const double*)__builtin_assume_aligned(sumtgt2s,32);
  bsepgains = (float*)__builtin_assume_aligned(bsepgains,32);
#endif  
  
  for (int ibin=0; ibin<nbins; ++ibin) {     
        
    double leftdiff = std::min(0.,-0.5*sumtgts[ibin]*sumtgts[ibin]*vdt::fast_inv(sumtgt2s[ibin]));
    //double leftdiff = std::min(0.,-0.5*sumtgts[ibin]*sumtgts[ibin]/sumtgt2s[ibin]);
    //double leftdiff = -0.5*sumtgts[ibin]*sumtgts[ibin]/sumtgt2s[ibin];

    double righttgtsum = sumtgt - sumtgts[ibin];
    double righttgt2sum = sumtgt2 - sumtgt2s[ibin];
    
    double rightdiff = std::min(0.,-0.5*righttgtsum*righttgtsum*vdt::fast_inv(righttgt2sum));
    //double rightdiff = std::min(0.,-0.5*righttgtsum*righttgtsum/righttgt2sum);
    //double rightdiff = -0.5*righttgtsum*righttgtsum/righttgt2sum;

	  
    //weighted improvement in variance from this split     
    //bsepgains[ibin] = std::max(0.,fulldiff - leftdiff - rightdiff);
    bsepgains[ibin] = fulldiff - leftdiff - rightdiff;
    
    //float valid =  sumtgt2s[ibin]==0.;// || righttgt2sum==0. || leftdiff>0. || rightdiff>0.);

    

    
  }  
  
  
}
