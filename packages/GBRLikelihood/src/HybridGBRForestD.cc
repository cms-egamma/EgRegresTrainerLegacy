#include "GBRLikelihood/HybridGBRForestD.h"

ClassImp(HybridGBRForestD) 


//_______________________________________________________________________
HybridGBRForestD::HybridGBRForestD(int ntargets)
{
  fInitialResponse.resize(ntargets,0.);
//  fResponses.resize(ntargets);
  fTrees.resize(ntargets);
}

//_______________________________________________________________________
HybridGBRForestD::~HybridGBRForestD() 
{
}
