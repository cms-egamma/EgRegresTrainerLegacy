
#include <assert.h>
#include "RegresTrainer/SmearingCorrection.h"
#include <iostream>
#include <math.h>



using namespace std;


/*****************************************************************/
float SmearingCorrection::call(const std::vector<float>& inputs) const
/*****************************************************************/
{
    bool isEB = (bool)inputs[0];
    float eta  = inputs[1];
    float r9   = inputs[2];
    
    float dsigMC = 0.;

    // Moriond smearing - old regression
    //if (isEB && fabs(eta)<1 && r9<0.94) dsigMC = 0.0109;
    //if (isEB && fabs(eta)<1 && r9>=0.94) dsigMC = 0.0099;
    //if (isEB && fabs(eta)>=1 && r9<0.94) dsigMC = 0.0182;
    //if (isEB && fabs(eta)>=1 && r9>=0.94) dsigMC = 0.0200;
    //if (!isEB && fabs(eta)<2 && r9<0.94) dsigMC = 0.0282;
    //if (!isEB && fabs(eta)<2 && r9>=0.94) dsigMC = 0.0309;
    //if (!isEB && fabs(eta)>=2 && r9<0.94) dsigMC = 0.0386;
    //if (!isEB && fabs(eta)>=2 && r9>=0.94) dsigMC = 0.0359;

    // Legacy smearing - new regression
    if (isEB && fabs(eta)<1 && r9<0.94) dsigMC    = 0.0094;
    if (isEB && fabs(eta)<1 && r9>=0.94) dsigMC   = 0.0092;
    if (isEB && fabs(eta)>=1 && r9<0.94) dsigMC   = 0.0182;
    if (isEB && fabs(eta)>=1 && r9>=0.94) dsigMC  = 0.0139;
    if (!isEB && fabs(eta)<2 && r9<0.94) dsigMC   = 0.0220;
    if (!isEB && fabs(eta)<2 && r9>=0.94) dsigMC  = 0.0229;
    if (!isEB && fabs(eta)>=2 && r9<0.94) dsigMC  = 0.0290;
    if (!isEB && fabs(eta)>=2 && r9>=0.94) dsigMC = 0.0234;

    float corr = m_random->Gaus(1., dsigMC);
    //cout << "Smearing correction (isEB=" << isEB << ",eta=" << eta << ",r9=" << r9 << ") = " << corr << ", dsig = " << dsigMC << "\n";
    return corr;
}
