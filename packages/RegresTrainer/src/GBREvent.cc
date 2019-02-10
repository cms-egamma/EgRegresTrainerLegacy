#include "RegresTrainer/GBREvent.h"

//_______________________________________________________________________
GBREvent::GBREvent(int nvars) : 
  fVars(new float[nvars]),
  fQuantiles(new int[nvars]),
  fTarget(0.0),
  fTransTarget(0.0),
  fWeight(1.0),
  fWeightedTransTarget(0.),
  fWeightedTransTarget2(0.)
{

}

//_______________________________________________________________________
GBREvent::~GBREvent() 
{
  delete[] fVars;
  delete [] fQuantiles;
  
}
