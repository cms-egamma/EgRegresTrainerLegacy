
#include <assert.h>
#include "RegresTrainer/TrackMomentumCorrection.h"
#include <iostream>
#include <math.h>



using namespace std;


/*****************************************************************/
float TrackMomentumCorrection::call(const std::vector<float>& inputs) const
/*****************************************************************/
{
    bool isEB           = (bool)inputs[0];
    int elClass         = (int)inputs[1];
    float trackMomentum = inputs[2];

    // tracker momentum corr corrections (Mykhailo Dalchenko)
    double corr = 1.;
    if (isEB)
    {
        if (elClass==0) corr = 1./(0.00104*sqrt(trackMomentum)+1);
        if (elClass==1) corr = 1./(0.0017*sqrt(trackMomentum)+0.9986);
        if (elClass==3) corr = 1./(1.004 - 0.00021*trackMomentum);
        if (elClass==4) corr = 0.995;
    } 
    else
    {
        if (elClass==3) corr = 1./(1.01432-0.00201872*trackMomentum+0.0000142621*trackMomentum*trackMomentum);
        if (elClass==4) corr = 1./(0.996859-0.000345347*trackMomentum);
    }
    if (corr<0.) corr = 1.;  // CC added protection


    //cout << "Track momentum correction (isEB=" << isEB << ",class=" << elClass << ",p=" << trackMomentum << ") = " << corr << "\n";
    return corr;
}
