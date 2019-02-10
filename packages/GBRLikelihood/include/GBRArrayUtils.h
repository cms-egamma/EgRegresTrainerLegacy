#ifndef GBRARRAYUTILS
#define GBRARRAYUTILS

class GBRArrayUtils {
  
  friend class RooHybridBDTAutoPdf;
  
protected:
 static void InitArrays(int *ns, double *tgts, double *tgt2s, float *bsepgains, const int nbins);
 static void ZeroArray(double *wscls, const int nbins);
 static void MinMaxQuants(int &minquant, int &maxquant, const int *quants, const int nev);
 static void FillBinQuants(int *binquants, const unsigned int offset, const unsigned int pscale, const unsigned int nquantiles, const unsigned int nbins);
 static void FillSepGains(const double *sumtgts, const double *sumtgt2s, float *bsepgains, const double fulldiff, const double sumtgt, const double sumtgt2, const int nbins);
 
};
 
#endif
