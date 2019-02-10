#include "GBRLikelihood/HybridGBREvent.h"

//_______________________________________________________________________
HybridGBREvent::HybridGBREvent(int nvars, int ntargets, int nparms) : 
  fVars(new float[nvars]),
  fTargets(new double[ntargets]),
  fTransTargets(new float[ntargets]),
  fTransTargets2(new float[ntargets]),
  fQuantiles(new int[nvars]),
  fDerivatives(new double[nparms]),
  fDerivatives2(new double[nparms]),
  fCurrentNodes(new unsigned int[ntargets]), 
  fPdfVal(0.),
  fWeight(1.0),
  fClass(0)
{

  for (int itgt=0; itgt<ntargets; ++itgt) {
    fTargets[itgt] = 0.;
    fTransTargets[itgt] = 0.;
    fTransTargets2[itgt] = 0.;
    fCurrentNodes[itgt] = 0;
  }  
  
  for (int ivar=0; ivar<nparms; ++ivar) {
    fDerivatives[ivar] = 0.;
    fDerivatives2[ivar] = 0.;
  }   
  
}

//_______________________________________________________________________
HybridGBREvent::~HybridGBREvent() 
{
  if (fVars) delete[] fVars;
  if (fTargets) delete[] fTargets;
  if (fTransTargets) delete[] fTransTargets;
  if (fTransTargets2) delete[] fTransTargets2;
  if (fQuantiles) delete [] fQuantiles;
  if (fDerivatives) delete [] fDerivatives;
  if (fDerivatives2) delete [] fDerivatives2;
  if (fCurrentNodes) delete [] fCurrentNodes;
  
}
