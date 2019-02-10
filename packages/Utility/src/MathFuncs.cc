#include "Utility/MathFuncs.hh"

#include "Math/ProbFuncMathCore.h"

const double MathFuncs::kPi = 2*asin(1.);
std::vector<double> MathFuncs::_logFacLookUp;
const double MathFuncs::kSqrt2Pi=sqrt(2*kPi);
TF1 MathFuncs::_probFunc("probFunc",MathFuncs::smearedProb,0,1000,4);
TF1 MathFuncs::_gausFunc("gausFunc",MathFuncs::gausFunc,-1000,1000,2);

TF1 MathFuncs::smearedDoublePeakUnBinnedLikelihoodFunc_("smearedDoublePeakUnBinnedLikelihoodFunc",
							MathFuncs::smearedDoublePeakUnBinnedLikelihood,
							0,1000,7);
TF1 MathFuncs::doublePeakUnBinnedLikeFunc_("doublePeakUnBinnedLikeFunc",
				       MathFuncs::calSmearedDoublePeakUnBinnedLikelihood,
				       0,1000,6);

// TF1 MathFuncs::smearedUnBinnedLikelihoodFunc_("smearedDoublePeakUnBinnedLikelihoodFunc",
// 					      MathFuncs::smearedUnBinnedLikelihood,
// 					      0,1000,7);
// TF1 MathFuncs::unBinnedLikeFunc_("doublePeakUnBinnedLikeFunc",
// 				 MathFuncs::calSmearedUnBinnedLikelihood,
// 				 0,1000,6);


double MathFuncs::detEtaFromEvnt(double evntEta,double z0)
{
 
  double thetaEvt = MathFuncs::etaToTheta(evntEta);
  double z = 129.4 / tan(thetaEvt); //129.4 is the average barrel radius
  double zTot = z+z0;

  if(fabs(zTot)<269){ //269 is an emperically derived number which means that < its likely in th barrel
    return zTot !=0 ? MathFuncs::thetaToEta(atan(129.4/zTot)) : 0.; //otherwise endcap time
  }
  double endcapZ = 319.2; //average z position of endcap
  if(evntEta<0) endcapZ*=-1;
  double rxy = tan(thetaEvt) * (endcapZ-z0);
  return MathFuncs::thetaToEta(atan(rxy/endcapZ));

}

double MathFuncs::evtEtaFromDet(double detEta,double z0)
{
//   double zcentroid = 31.44;
//   if(evntEta<0) zcentroid*=-1;
//   double evntCotTheta = sinh(evntEta);
//   double detCotTheta = zcentroid/(zcentroid-z0) * evntCotTheta;
//   return asinh(detCotTheta);

  double thetaDet = etaToTheta(detEta);
  double z = 129. / tan(thetaDet);

 
  double zTot = z-z0;
   if(detEta<0) z*=-1;
    if(detEta<0) zTot*=-1;
  return thetaToEta(atan(z*tan(thetaDet)/zTot));



}

double MathFuncs::detEtaToEvt(TVector3 caloPos,double z0)
{
  caloPos.SetZ(caloPos.Z()-z0);
  return caloPos.Eta();
}

double MathFuncs::evtEtaToDet(TVector3 caloPos,double z0)
{
  caloPos.SetZ(caloPos.Z()+z0);
  return caloPos.Eta();
}


double MathFuncs::etaToTheta(double eta)
{
  //  if(eta<0) return -2*atan(exp(eta));
  //  else return 2*atan(exp(-1*eta));
  return 2*atan(exp(-1*eta));
  //else return 2*atan(exp(-1*eta));

}

double MathFuncs::thetaToEta(double theta)
{
  //first bounds check theta to get into -pi/2 - pi/2 range
  while( fabs(theta) > MathFuncs::kPi/2.){
    if(theta>0) theta-=MathFuncs::kPi;
    else theta+=MathFuncs::kPi;
  }
  //now check sign
  if(theta<0) return log(tan(fabs(theta/2.)));
  else return -1.*log(tan(theta/2.));
}

double MathFuncs::cosThetaStar(float ele1Et,float ele1Eta,float ele1Phi,
			       float ele2Et,float ele2Eta,float ele2Phi,int ele1Charge)
{
  TLorentzVector ele1P4;
  ele1P4.SetPtEtaPhiM(ele1Et,ele1Eta,ele1Phi,0);
  TLorentzVector ele2P4;
  ele2P4.SetPtEtaPhiM(ele2Et,ele2Eta,ele2Phi,0);
 
  if(ele1Charge==1) return cosThetaStar(ele1P4,ele2P4);
  else return cosThetaStar(ele2P4,ele1P4);

}

double MathFuncs::cosThetaStar(const TLorentzVector& ele1P4,const TLorentzVector& ele2P4)
{
 
  TLorentzVector diEleP4 = ele1P4+ele2P4;
  
  
  double p1Plus =  (ele1P4.E()+ele1P4.Pz());
  double p1Minus = (ele1P4.E()-ele1P4.Pz());

  double p2Plus =  (ele2P4.E()+ele2P4.Pz());
  double p2Minus = (ele2P4.E()-ele2P4.Pz());

  return (p1Plus*p2Minus-p1Minus*p2Plus)/sqrt(diEleP4.Mag2()*(diEleP4.Mag2()+diEleP4.Perp2()));

}

double MathFuncs::phi(double px,double py)
{
  if(px!=0) return atan2(py,px); // pass as y , x
  else if(py!=0) return py > 0 ? kPi/2 : -kPi/2;
  else return 0.;
}

//return in range from 0 to 2*pi
double MathFuncs::degreeToRad(double degree)
{
  double rad = degree/180 * kPi;
  return normAngleRange(rad,0,2*kPi);

}

//casts the angle into the range specified adding the range difference to do it
double MathFuncs::normAngleRange(double angle,double minRange,double maxRange)
{
  double normedAngle = angle; 
  double range = maxRange-minRange;
  while (normedAngle<minRange) normedAngle+=range;
  while (normedAngle>maxRange) normedAngle-=range;
  return normedAngle;
}

//range is either 2pi or 360 but anything else could work
double MathFuncs::angleDiff(double angle1,double angle2,double range)
{
  double normedAngle1 = normAngleRange(angle1,0,range);
  double normedAngle2 = normAngleRange(angle2,0,range);
  double diff = normedAngle1 -normedAngle2;
  if(diff>range/2.) diff = range-diff;
  if(diff<-range/2.) diff = -range-diff;


  return diff;
}

bool MathFuncs::isAngleDiffInRadRange(double angleDiff,double minAngle,double maxAngle)
{
  double normedMinAngle = normAngleRange(minAngle,0,2*kPi);
  double normedMaxAngle = normAngleRange(maxAngle,0,2*kPi);
  if(fabs(angleDiff)>normedMinAngle && fabs(angleDiff)<normedMaxAngle) return true;
  else return false;

}

// double MathFuncs::angleDegMinDiff(double angle1,double angle2)
// {
  
//   double diff = angle1-angle2;
//   while (diff<-180) diff+=180;
//   while (diff>180) diff-=180;
//   return diff;
// }

int MathFuncs::getCoordQuadrant(double xCoord,double yCoord)
{
  int quadrantNr=0;
  if(xCoord>0){
    if(yCoord>0) quadrantNr=1; //x=+ve, y=+ve 0 < phi < pi/2
    else quadrantNr = 4; //x=+ve, y=-ve 3/2 < phi < 2*pi
  }else{
    if(yCoord>0) quadrantNr=2; //x=-ve, y=+ve pi/2 <phi< pi
    else quadrantNr=3; //x=-ve, y=-ve pi < phi < 3/2 pi
  }
  return quadrantNr;
} 

//returns phi in range from 0-2pi from the x and y coordinates
double MathFuncs::getAngleFromCoord(double xCoord,double yCoord)
{
  
  int quadrantNr=getCoordQuadrant(xCoord,yCoord);
  
  double angle = atan(yCoord/xCoord);
  // std::cout <<" angle "<<angle<<" x "<<xCoord<<" y "<<yCoord<<" quadrant "<<quadrantNr;
  //now need to convert angle from range of -pi/2 to pi/2 into 0 to 2*pi
  switch(quadrantNr){
  case 1:
    break; //quadrant 1 is already fine
  case 2:
    angle += MathFuncs::kPi;
    break;
  case 3:
    angle += MathFuncs::kPi;
    break;
  case 4:
    angle += 2*MathFuncs::kPi;
    break;
  default:
    angle = 999.;
  }
  // std::cout <<" corr angle "<<angle<<std::endl;

  return angle;
}

float MathFuncs::calEffFromPassFail(float nrPass,float nrFail)
{
  float nrTot = nrPass+nrFail;
  if(nrTot!=0) return nrPass/nrTot;
  else return 0.;
}

float MathFuncs::calEffErrFromPassFail(float nrPass,float nrPassErr,float nrFail,float nrFailErr)
{
  float nrTot = nrPass+nrFail;
  if(nrTot!=0){
    return std::sqrt(nrPass*nrPass*nrFailErr*nrFailErr + nrFail*nrFail*nrPassErr*nrPassErr )/(nrTot*nrTot);
  }else return 0.000000001;
}

//Using the polar form of the Box-Muller transformation to generate gausian values
//Intensionally discarding one of my randomly generated gausian values, inefficent 
//but as as the number of generated events is small can get away with this. 
double MathFuncs::randGaus()
{  
  //  double x1, x2, w, y1, y2;
  double x1, x2, w, y1;
  do {
    x1 = 2.0 * (double)rand()/((double)RAND_MAX+1.) - 1.0;
    x2 = 2.0 * (double)rand()/((double)RAND_MAX+1.) - 1.0;
    w = x1 * x1 + x2 * x2;
  } while ( w >= 1.0 );
  
  w = sqrt( (-2.0 * log( w ) ) / w );
  y1 = x1 * w;
  //y2 = x2 * w;
  
  return y1;
}


double MathFuncs::logFactorial(int n)
{
  if(n>((int)_logFacLookUp.size())-1) _fillLogFacLookup(n);
  return _logFacLookUp[n];
}

void MathFuncs::_fillLogFacLookup(int n)
{
  if(_logFacLookUp.size()==0) _logFacLookUp.push_back(0);
  for(unsigned i=_logFacLookUp.size();((int)i)<=n;i++){   
    _logFacLookUp.push_back(_logFacLookUp[i-1] + log((double)i));
  }
}


double MathFuncs::poisson(int n,double expect)
{ 
  if(expect<=0) return 1.;
  if(n<0) return 0.;
 
  
  double nLogFac = logFactorial(n);
  double logProb = -expect + n*log(expect) - nLogFac;
  return exp(logProb);

  
}

std::pair<float,float> MathFuncs::calEffAndErr(float nrPass,float nrPassW2,float nrTot,float nrTotW2)
{
  float eff = nrPass/nrTot;
  float nrFail = nrTot-nrPass;
  float nrFailW2 = nrTotW2-nrPassW2;
  
  float err = sqrt(nrPass*nrPass*nrFailW2+nrFail*nrFail*nrPassW2)/nrTot/nrTot;
  return std::pair<float,float>(eff,err);
}

std::pair<float,float> MathFuncs::calEffAndErr(float nrPass,float nrTot)
{
  float eff = nrPass/nrTot;
  float nrFail = nrTot-nrPass;
  
  float err = sqrt(nrPass*nrPass*nrFail+nrFail*nrFail*nrPass)/nrTot/nrTot;
  return std::pair<float,float>(eff,err);
}

double MathFuncs::logPoisson(int n,double expect)
{ 
  if(expect<=0) return 1.;
  if(n<0) return 0.;
 
  
  double nLogFac = logFactorial(n);
  double logProb = -expect + n*log(expect) - nLogFac;
  return logProb;
  
}

double MathFuncs::intPoissonProb(int minN,int maxN,double expect)
{
  if(expect<=0) return expect;
  double totalProb = 0.;
  for(int i=minN;i<maxN;i++){
    totalProb += poisson(i,expect);
  }
  return totalProb;
}

//para[0] = start value
//para[1] = end (data) value
//para[2] = bkg ground expect
//para[3] = bkg error
double MathFuncs::smearedProb(double *x,double *para)
{
  int paraInt[2];
  paraInt[0] = static_cast<int>(para[0]);
  paraInt[1] = static_cast<int>(para[1]);
  double prob = MathFuncs::intPoissonProb(paraInt[0],paraInt[1],*x)*MathFuncs::gaus(*x,para[2],para[3]);
  // std::cout <<"expect "<<*x<<" mean expect "<<para[2]<<" prob "<<prob<<std::endl;
 return prob;
}


//para[0] nrObs peak 1
//para[1] nrObs peak 2
//para[2] sig peak 1
//para[3] sig peak 2
//para[4] bkg peak 1
//para[5] bkg peak 2
//para[6] bkg peak 1 error
//assumes errors are correlated and an increase in bkg of x% will affect both peaks
double MathFuncs::smearedDoublePeakUnBinnedLikelihood(double* x,double *para)
{
  

   double probPeak1  = MathFuncs::poisson(static_cast<int>(para[0]),*x+para[2]);
  //  double probPeak1 =1.;
   //  std::cout <<" nr obs "<<para[1]<<" bkg "<<*x*para[5]/para[4]+para[3]<<std::endl;
   double probPeak2  = MathFuncs::poisson(static_cast<int>(para[1]),*x*para[5]/para[4]+para[3]);
   //  std::cout <<" nr obs "<<para[1]<<" bkg "<<*x*para[5]/para[4]+para[3]<<" prob "<<probPeak2<<std::endl;
   //double probPeak2 =  1.;
  double prob = probPeak1*probPeak2*MathFuncs::gaus(*x,para[4],para[6])*1E12;
  return prob;
}

//para[0] nrObs peak 1
//para[1] sig peak 1
//para[2] bkg peak 1
//para[3] bkg peak 1 error
//assumes errors are correlated and an increase in bkg of x% will affect both peaks
// double MathFuncs::smearedUnBinnedLikelihood(double* x,double *para)
// {
  

//    double probPeak1  = MathFuncs::poisson(static_cast<int>(para[0]),*x+para[2]);

//    //  std::cout <<" nr obs "<<para[1]<<" bkg "<<*x*para[5]/para[4]+para[3]<<" prob "<<probPeak2<<std::endl;
//    //double probPeak2 =  1.;
//   double prob = probPeak1*MathFuncs::gaus(*x,para[2],para[2])*1E12;
//   return prob;
// }

//x = n_sig
//para[0] = obs p1
//para[1] = obs p2
//para[2] = n_peak 2 / n_peak 1
//para[3] = bkg peak 1
//para[4] = bkg peak 2
//para[5] = bkg peak 1  error 
double MathFuncs::calSmearedDoublePeakUnBinnedLikelihood(double* x,double *para)
{
  float minBkg = para[3]-3*para[5];
  float maxBkg = para[3]+3*para[5];
  
  double paraToFunc[7];
  paraToFunc[0] = para[0];
  paraToFunc[1] = para[1];
  paraToFunc[2] = *x;
  paraToFunc[3] = *x*para[2];
  paraToFunc[4] = para[3];
  paraToFunc[5] = para[4];
  paraToFunc[6] = para[5];
  smearedDoublePeakUnBinnedLikelihoodFunc_.SetParameters(paraToFunc);
  return smearedDoublePeakUnBinnedLikelihoodFunc_.Integral(minBkg,maxBkg);
}

double MathFuncs::cal95CLDoublePeakUnBinnedLikelihood(double nrObs1,double nrObs2,double nrBkg1,double nrBkg2,double ratioPeak2ToPeak1,double bkgErr)
{
  double para[6];
  para[0] = nrObs1;
  para[1] = nrObs2;
  para[2] = ratioPeak2ToPeak1;
  para[3] = nrBkg1;
  para[4] = nrBkg2;
  para[5] = bkgErr;

  doublePeakUnBinnedLikeFunc_.SetParameters(para);
  float totIntegral = doublePeakUnBinnedLikeFunc_.Integral(0,nrBkg1*10);
  
  double conFidLvl = 0.95;
  double limit = 0;
  double stepSize = 1;
  double minStepSize = 0.001;
  bool beenBelow95Lvl = false; //as we dont chose 0 as our starting value, there is a change we might start above 95% level
                               //this varible checks that there has been atleast one result below the 95% level, ie we're converging on the limit
  do{
    
    double intResult = doublePeakUnBinnedLikeFunc_.Integral(0,limit);
    if(intResult>conFidLvl*totIntegral){
      limit -=stepSize;
      if(limit<0){
	limit=0;
       beenBelow95Lvl=true; //startCross has got to be below 95%....
      }
      if(beenBelow95Lvl) stepSize/=10.;
    }else{
      if(!beenBelow95Lvl) beenBelow95Lvl = true;
      limit +=stepSize;
      //std::cout <<"limit "<<limit<<" intResult "<<intResult<<" tot in "<<totIntegral<<std::endl;
    }
  }while(stepSize>=minStepSize && limit<100);
  if(limit >=100) std::cout <<"Integral didnt converge"<<std::endl;
  
  return limit;
}


//  //has precomputed values from 0 to 1000 in 10000 steps
//  //lets use that memory :)
//  double MathFuncs::smearedProbApprox(double *x,double *para)
//  {
//    int paraInt[2];
//    paraInt[0] = static_cast<int>(para[0]);
//    paraInt[1] = static_cast<int>(para[1]);
  
//    if(paraInt[0]==0 && x<100 && paraInt[2]<200){ //can approximate
//      double xLarge = x*10000;
//      int index = static_cast<int>(xLarge);
//      if(xLarge-index>0.5) index++; //rounding up if needs be
//      return intProbLookUpTable_[index][paraInt[1]]*MathFuncs::gaus(*x,para[2],para[3]); 
//    }else { //cant approximate
//      return MathFuncs::intPoissonProb(paraInt[0],paraInt[1],*x)*MathFuncs::gaus(*x,para[2],para[3]);
//    }
//  }
//   //  std::cout <<"expect "<<*x<<" mean expect "<<para[2]<<" prob "<<prob<<std::endl;
//  return prob;
// }

//with errors now
double MathFuncs::intPoissonProb(int minN,int maxN,double expect,double expectErr)
{
  if(expectErr==0) return intPoissonProb(minN,maxN,expect);

  double para[4];
  para[0] = minN;
  para[1] = maxN;
  para[2] = expect;
  para[3] = expectErr;
  //para[3] =0;
  //need to work out the integration limits as roots too dumb to do it
  double minInt = expect-3*expectErr;
  double maxInt = expect+3*expectErr;
  if(minInt<0) minInt=0; //constain to being non-negative
  //std::cout << "gaus int "<<intGaus(expect,expectErr,minInt,maxInt)<<std::endl;
  // std::cout <<"start prob calc, expect "<<expect<<std::endl;
  _probFunc.SetParameters(para);
  double prob =_probFunc.Integral(minInt,maxInt,1E-14)/0.997300203936739793; //intGaus(expect,expectErr,minInt,maxInt);
 // std::cout <<"end prob calc "<<std::endl;
 return prob;
}

double MathFuncs::getMedian(std::list<double> &list)
{
  double median=0.;
  std::list<double>::iterator listIt = list.begin();
   if(list.size() % 2 !=0) {
    int middleNr = (list.size()+1)/2;
    for(int i=0;i<middleNr;i++) listIt++;
    median = *listIt;
  }else{ //even number, take median as halfway between the two middle values
    int middleNr = (list.size()+1)/2;
    for(int i=0;i<middleNr;i++) listIt++;
    median= *listIt;
    listIt++;
    median+= *listIt;
    median/=2.;
  }
   return median;
}

double MathFuncs::getMedian(const std::vector<double> &vec)
{
  double median=0.;
  
  //odd number, definate median
  if(vec.size() % 2 !=0) {
    int middleNr = (vec.size()+1)/2;
    median = vec[middleNr];
  }else{ //even number, take median as halfway between the two middle values
    int middleNr = (vec.size()+1)/2;
    median= vec[middleNr];
    if(middleNr+1 <(int) vec.size()) median+= vec[middleNr+1];
    median/=2.;
  }
  return median;
}

double MathFuncs::getWeightedListValue(std::list<double> &list,double listNr)
{
  // if(listNr>list.size() || listNr<0) return list.back();
  if(listNr>=list.size()) return list.back();
  else if(listNr<=0) return list.front();

  double tempListNr; //just a tempory value to hold the int part of list nr before it is cast into an interger
  double fracNr = std::modf(listNr,&tempListNr); //decomposing the listNr into frac and int parts
  int intNr = (int) tempListNr;
  
  std::list<double>::iterator listIt = list.begin();
  for(int i=0;i<intNr;i++) listIt++; //incrementing the list pointer to the element we want
  double lowerValue = *listIt;
  listIt++;
  double upperValue = *listIt;

  return lowerValue*(1-fracNr) + upperValue*fracNr;
}

//I may be offset by one from where I want to be, need to think about it (12/09/12)
double MathFuncs::getWeightedListValue(const std::vector<double> &vec,double index)
{
  //if(index>=vec.size() || index<0) return -999;
  if(index>=vec.size()-1) return vec.back();
  else if(index<=0) return vec.front();

  double tempIndex; //just a tempory value to hold the int part of the index before it is cast into an interger
  double fracNr = std::modf(index,&tempIndex); //decomposing the index into frac and int parts
  size_t intNr = static_cast<size_t>(tempIndex);
  
  if(intNr+1>=vec.size() || fracNr==0) return vec[intNr];
  else{
    double lowerValue = vec[intNr];
    double upperValue = vec[intNr+1];

    return lowerValue*(1-fracNr) + upperValue*fracNr;
  }
}


double MathFuncs::round(double numToRound,float nrFigures,bool isDecPlaces)
{
  if(isDecPlaces){
    int powerOf10 = MathFuncs::power(10,nrFigures);
    if(numToRound>0) return std::floor( ( numToRound * powerOf10 ) +0.5 ) / powerOf10;
    else return std::ceil( ( numToRound * powerOf10 ) -0.5 ) / powerOf10; 
  }else{
    std::cout <<"not implimented yet"<<std::endl;
    return 0.;
  }
} 

int MathFuncs::roundToInt(double numToRound)
{
  if(numToRound>0) return std::floor(numToRound+0.5);
  else return std::ceil(numToRound-0.5);
}


double MathFuncs::randNr(double min,double max)
{
  return (double)rand()/(double) RAND_MAX * (max-min) + min;
}

//calcualates the factorial of an integer number
double MathFuncs::factorial(int n)
{
  double result=1.;
  for(int i=0;i<n;i++){
    result *= n-i;
  }
  return result;
}


double MathFuncs::cosThetaStarCS(const TLorentzVector& eleP4,const TLorentzVector& posP4)
{
  const TLorentzVector zP4(eleP4+posP4);
  double numerator = (eleP4.E()+eleP4.Pz())*(posP4.E()-posP4.Pz()) - (eleP4.E()-eleP4.Pz())*(posP4.E()+posP4.Pz());
  double denominator = zP4.Mag()*sqrt(zP4.Mag2() + zP4.Perp2());
  double cosThetaStar = denominator!=0 ? numerator/denominator : -999;
  if(zP4.Pz()!=0) cosThetaStar*=fabs(zP4.Pz())/zP4.Pz();
  return cosThetaStar;
}

double MathFuncs::cosThetaStarCS(const TLorentzVector& ele1P4,int ele1Charge,const TLorentzVector& ele2P4,int ele2Charge)
{
  if(ele1Charge*ele2Charge<0){ //os, easy
    return ele1Charge==-1 ? cosThetaStarCS(ele1P4,ele2P4) : cosThetaStarCS(ele2P4,ele1P4);
  }else{ //not so easy, pick barrel as best, then highest energy
    int ele1Region = fabs(ele1P4.Eta())<1.5 ? 0 : 1;
    int ele2Region = fabs(ele2P4.Eta())<1.5 ? 0 : 1;
    
    if(ele1Region+ele2Region!=1){ //so EB-EB, and EE-EE 
      if(ele1P4.Et()>ele2P4.Et()){ //ele1 is highest, trust that charge
	return ele1Charge==-1 ? cosThetaStarCS(ele1P4,ele2P4) : cosThetaStarCS(ele2P4,ele1P4);
      }else return ele2Charge==-1 ? cosThetaStarCS(ele2P4,ele1P4) : cosThetaStarCS(ele1P4,ele2P4);
    }else if(fabs(ele1P4.Eta())<1.5) {//one is barrel, one endcap, is it ele1
      return ele1Charge==-1 ? cosThetaStarCS(ele1P4,ele2P4) : cosThetaStarCS(ele2P4,ele1P4);
    }else return ele2Charge==-1 ? cosThetaStarCS(ele2P4,ele1P4) : cosThetaStarCS(ele1P4,ele2P4);//no must be ele2
  }
}

// double MathFuncs::calDeltaR2(double eta1,double phi1, double eta2,double phi2)
// {
//   double dEta = eta1-eta2;
//   double dPhi = deltaPhi(phi1,phi2);
//   return dEta*dEta + dPhi*dPhi;
// }

// double MathFuncs::deltaPhi(double phi1,double phi2)
// {
//   double dPhi = phi1-phi2;
//   if(dPhi>kPi) dPhi -= 2*kPi;
//   if(dPhi<=-kPi) dPhi += 2*kPi;

// }


double MathFuncs::intGaus(double mean,double sigma)
{
  return intGaus(mean,sigma,mean-3*sigma,mean+3*sigma);
}

double MathFuncs::intGaus(double mean,double sigma,double min,double max)
{
  double para[2];
  para[0] = mean;
  para[1] = sigma;
  _gausFunc.SetParameters(para);
  return _gausFunc.Integral(min,max);
}

double MathFuncs::power(double x,int n)
{
  double result = 1.;
  for(int i=0;i<n;i++){
    result *=x;
  }
  return result;
}

int MathFuncs::power(int x,int n)
{
  int result =1;
  for(int i=0;i<n;i++){
    result *=x;
  }
  return result;
}

std::pair<int,float> MathFuncs::estQuantile(const std::vector<int>& dataPoints,float quantile)
{
  if(dataPoints.empty()) return std::pair<int,float>(0,0.);
  size_t index = dataPoints.size()*quantile+0.5;
  if(index>0) index--; //-1 as we are offset by one
  float uncert = martizJarrettUncert(dataPoints,index);
 
 
  return std::pair<int,float>(dataPoints[index],uncert);
}

float MathFuncs::martizJarrettUncert(const std::vector<int>& dataPoints,size_t m)
{
  size_t n=dataPoints.size();
  int a=m-1;
  int b=n-m;
  
  float c1=0;
  float c2=0;

  for(size_t i=0;i<n;i++){
    c1+=(ROOT::Math::beta_cdf(static_cast<float>(i+1)/n,a,b) - ROOT::Math::beta_cdf(static_cast<float>(i)/n,a,b))*dataPoints[i];
    c2+=(ROOT::Math::beta_cdf(static_cast<float>(i+1)/n,a,b) - ROOT::Math::beta_cdf(static_cast<float>(i)/n,a,b))*dataPoints[i]*dataPoints[i];  
  }
  
  return sqrt(c2-c1*c1);
}
