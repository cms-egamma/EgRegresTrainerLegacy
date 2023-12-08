#include "RegresTrainer/GBRTrainer.h"
#include "RegresTrainer/GBREvent.h"
#include "CondFormats/GBRForest/interface/GBRForest.h"
#include "TTree.h"
#include "TTreeFormula.h"
#include <assert.h>
#include <malloc.h>
#include <limits>
#include <omp.h>
#include <algorithm>
#include <iostream>

//_______________________________________________________________________
GBRTrainer::GBRTrainer() : 
  fMinEvents(2000),
  fShrinkage(0.1),
  fNQuantiles(std::numeric_limits<unsigned short>::max()+1),
  fNBinsMax(128),
  fTransitionQuantile(0.7),
  fMinCutSignificance(-99.0),
  fRandomSeedFormula("1"),
  fEventWeightFormula("1"),
  _sepgains(0),
  _ws(0)
{

}

//_______________________________________________________________________
GBRTrainer::~GBRTrainer() 
{

  //clear arrays
  if (_sepgains) {
    delete[] _sepgains;
    delete[] _cutvals;
    delete[] _nlefts;
    delete[] _nrights;
    delete[] _sumwlefts;
    delete[] _sumwrights;
    delete[] _sumtgtlefts;
    delete[] _sumtgtrights;
    delete[] _leftvars;
    delete[] _rightvars;
    delete[] _fullvars;
    delete[] _bestbins;  
  }
  
  if (_ws) {

    for (unsigned int ivar=0; ivar<fInputVars.size(); ++ivar) {
      delete[] _ws[ivar];
      delete[] _ws2[ivar];
      delete[] _ns[ivar];
      delete[] _tgts[ivar];
      delete[] _tgt2s[ivar];
      
      delete[] _sumws[ivar];
      delete[] _sumws2[ivar];
      delete[] _sumns[ivar];
      delete[] _sumtgts[ivar];
      delete[] _sumtgt2s[ivar];
      delete[] _varvals[ivar];   
      delete[] _bsepgains[ivar];
      delete[] _bsepgainsigs[ivar];
      
      
      delete[] _quants[ivar];
      delete[] _bins[ivar];
      
      delete[] fQuantileMaps[ivar];
    }
    
    delete[] _ws;
    delete[] _ws2;
    delete[] _ns;
    delete[] _tgts;
    delete[] _tgt2s;
    
    delete[] _sumws;
    delete[] _sumws2;
    delete[] _sumns;
    delete[] _sumtgts;
    delete[] _sumtgt2s;
    delete[] _varvals;  
    delete[] _bsepgains;
    
    delete[] _quants;
    delete[] _bins;
    
    delete[] fQuantileMaps;
  }
  
}

//_______________________________________________________________________
const GBRForest *GBRTrainer::TrainForest(int ntrees)
{
  
  const int nvars = fInputVars.size();
    
  Long64_t nev = 0;  
  
  printf("first loop, count events\n");
  //loop over trees to count training events with non-zero weight;
  for (unsigned int itree=0; itree<fTrees.size(); ++itree) {
    TTree *tree = fTrees[itree];
    //printf("Target = %s\n",fTargetVar.c_str());
    //printf("Cuts   = %s\n",fTrainingCut.c_str());
    //printf("Seed   = %s\n",fRandomSeedFormula.c_str());
    //printf("Weight = %s\n",fEventWeightFormula.c_str());
    TTreeFormula targetform(fTargetVar.c_str(),fTargetVar.c_str(),tree);
    TTreeFormula cutform(fTrainingCut.c_str(),fTrainingCut.c_str(),tree);
    TTreeFormula seedform(fRandomSeedFormula.c_str(), fRandomSeedFormula.c_str(), tree);
    TTreeFormula weightform(fEventWeightFormula.c_str(), fEventWeightFormula.c_str(), tree);
    for (Long64_t iev=0; iev<tree->GetEntries(); ++iev) {
      if (iev%100000==0) printf("%i\n",int(iev));
      tree->LoadTree(iev);
      if ((fTreeWeights[itree]*cutform.EvalInstance())!=0.) {
        double rdm = 1.;
        double fraction = 1.;
        if(fEventWeightFormula!="")
        {
            int seed = (unsigned int)(seedform.EvalInstance());
            fRandom.SetSeed(seed);
            rdm = fRandom.Uniform();
            fraction = weightform.EvalInstance();
        }
        if(rdm<=fraction)
            ++nev;
      }
    }
  }
  
  printf("nev = %i, nvar = %i\n",int(nev),nvars);
  PrintInputVars();

  //initialize arrays

  _sepgains = new float[nvars];
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

  _ws = new float*[nvars];
  _ws2 = new float*[nvars];
  _ns = new int*[nvars];
  _tgts = new float*[nvars];
  _tgt2s = new float*[nvars];  
  _sumws = new float*[nvars];
  _sumws2 = new float*[nvars];
  _sumns = new int*[nvars];
  _sumtgts = new float*[nvars];
  _sumtgt2s = new float*[nvars];
  _varvals = new float*[nvars];    
  _bsepgains = new float*[nvars];
  _bsepgainsigs = new float*[nvars];
  
  _quants = new int*[nvars];
  _bins = new int*[nvars];
  
  fQuantileMaps = new float*[nvars];
  
  for (int ivar=0; ivar<nvars; ++ivar) {
    _ws[ivar] = new float[fNBinsMax];
    _ws2[ivar] = new float[fNBinsMax];    
    _ns[ivar] = new int[fNBinsMax];
    _tgts[ivar] = new float[fNBinsMax];
    _tgt2s[ivar] = new float[fNBinsMax];  
    _sumws[ivar] = new float[fNBinsMax];
    _sumws2[ivar] = new float[fNBinsMax];
    _sumns[ivar] = new int[fNBinsMax];
    _sumtgts[ivar] = new float[fNBinsMax];
    _sumtgt2s[ivar] = new float[fNBinsMax];
    _varvals[ivar] = new float[fNBinsMax];  
    _bsepgains[ivar] = new float[fNBinsMax];      
    _bsepgainsigs[ivar] = new float[fNBinsMax];      
    
    _quants[ivar] = new int[nev];
    _bins[ivar] = new int[nev];
    
    fQuantileMaps[ivar] = new float[fNQuantiles];
  }
      
  std::vector<GBREvent*> evts;
  evts.reserve(nev);
  
  double sumw = 0.;
  
  printf("second loop, fill events in memory\n");
  //loop over trees to fill arrays and event vector
  
  int nNaNwarning = 0;
  for (unsigned int itree=0; itree<fTrees.size(); ++itree) {
    TTree *tree = fTrees[itree];
    
    //initialize TTreeFormulas to read variables from TTree
    std::vector<TTreeFormula*> inputforms;
    for (std::vector<std::string>::const_iterator it = fInputVars.begin(); 
        it != fInputVars.end(); ++it) {
      inputforms.push_back(new TTreeFormula(it->c_str(),it->c_str(),tree));
    }
    
    TTreeFormula targetform(fTargetVar.c_str(),fTargetVar.c_str(),tree);
    TTreeFormula cutform(fTrainingCut.c_str(),fTrainingCut.c_str(),tree);  
    TTreeFormula seedform(fRandomSeedFormula.c_str(), fRandomSeedFormula.c_str(), tree);
    TTreeFormula weightform(fEventWeightFormula.c_str(), fEventWeightFormula.c_str(), tree);

    for (Long64_t iev=0; iev<tree->GetEntries(); ++iev) {
      if (iev%100000==0) printf("%i\n",int(iev));
      tree->LoadTree(iev);
      
      float weight = fTreeWeights[itree]*cutform.EvalInstance();
      
      if (weight==0.) continue; //skip events with 0 weight
      double rdm = 1.;
      double fraction = 1.;
      if(fEventWeightFormula!="")
      {
          int seed = (unsigned int)(seedform.EvalInstance());
          fRandom.SetSeed(seed);
          rdm = fRandom.Uniform();
          fraction = weightform.EvalInstance();
      }
      if(rdm>fraction) continue;
  
      sumw += weight;
      
      evts.push_back(new GBREvent(nvars));
      GBREvent *evt = evts.back();
      evt->SetWeight(weight);
      evt->SetTarget(targetform.EvalInstance());
      
      //if(targetform.EvalInstance()<=0.)
      //    printf("target = %5f\n",targetform.EvalInstance());
      
      for (int i=0; i<nvars; ++i) {
          float var = inputforms[i]->EvalInstance();
          if(var!=var)// protection against nan
          {
              var = 0.;
              if(nNaNwarning<=10)
                  std::cerr << "WARNING: found input variable " << i << " = NaN. Set it to 0. This may be inappropriate. Please fix this.\n";
              if(nNaNwarning==10)
                  std::cerr << "Last Warning...\n";
              nNaNwarning ++;
          }
          evt->SetVar(i,var);
      }

    }

    for (std::vector<TTreeFormula*>::const_iterator it = inputforms.begin(); 
        it != inputforms.end(); ++it) {
      delete *it;
    }

  }

  
  //map of input variable quantiles to values
  //fQuantileMaps.resize(nvars, std::vector<float>(fNQuantiles));
  
  //parallelize building of quantiles for each input variable
  //(sorting of event pointer vector is cpu-intensive)
#pragma omp parallel for
  for (int ivar=0; ivar<nvars; ++ivar) {
    //printf("sorting var %i\n",ivar);
        
    std::map<int,float,std::greater<float> > tmpmap;    
    std::vector<GBREvent*> evtsvarsort(evts.begin(),evts.end());
    
    std::sort(evtsvarsort.begin(),evtsvarsort.end(),GBRVarCMP(ivar));
    
    double sumwq = 0;
    for (unsigned int iev=0; iev<evtsvarsort.size(); ++iev) {
      sumwq += evtsvarsort[iev]->Weight();
      int quant = int((sumwq/sumw)*(fNQuantiles-1));
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
    
  //sort events by target and compute median
  std::sort(evts.begin(),evts.end(),GBRTargetCMP());
  double medsumw = 0;
  float median = 0.;

  std::vector<double> targets;
  targets.reserve(evts.size());
  std::vector<GBREvent*>::const_iterator medit=evts.begin();
  while(medsumw<(0.5*sumw) && medit!=evts.end()) {
    medsumw += (*medit)->Weight();
    median = (*medit)->Target();
    targets.push_back((*medit)->Target());
    ++medit;
  }
  
  //set initial response and recompute targets
  GBRForest *forest = new GBRForest;
  forest->SetInitialResponse(median);
  

  for (std::vector<GBREvent*>::iterator it=evts.begin(); it!=evts.end(); ++it) {
    (*it)->SetTarget((*it)->Target()-median);
  }  
  
  //sort by absolute value of the recomputed target and computed transformed target
  //according to huber loss function derivative (cutoff of outliers)
  std::sort(evts.begin(),evts.end(),GBRAbsTargetCMP());
  double transumw = 0.;
  float transition = 0.;
  std::vector<GBREvent*>::const_iterator transit=evts.begin();
  while(transumw<(fTransitionQuantile*sumw) && transit!=evts.end()) {
    transumw += (*transit)->Weight();
    transition = std::abs((*transit)->Target());  
    ++transit;
  } 
  
  for (std::vector<GBREvent*>::iterator it=evts.begin(); it!=evts.end(); ++it) {
    float tgt = (*it)->Target();
    if (std::abs(tgt)<transition) (*it)->SetTransTarget(tgt);
    else if (tgt<0.) (*it)->SetTransTarget(-transition);
    else (*it)->SetTransTarget(transition);
  }    

  
  printf("nev = %i, sumw = %5f, median = %5f, transition = %5f\n",int(nev), sumw, median,transition);
  //printf("Initial response = %5f\n",forest->GetInitialResponse());
  
  //loop over requested number of trees
  for (int itree=0; itree<ntrees; ++itree) {
    printf("tree %i\n",itree);

    //sort events by recomputed target, which is expected/required for correct computation
    //of median for each terminal mode
    std::sort(evts.begin(),evts.end(),GBRTargetCMP());
      
    forest->Trees().push_back(GBRTree());
    GBRTree &tree = forest->Trees().back();

    //train a single tree recursively from the root node
    TrainTree(evts,sumw,tree,nvars,transition);
    
    //stop training if root node is already terminal
    if (itree==(ntrees-1) || tree.LeftIndices().size()==0 || tree.LeftIndices().front()==tree.RightIndices().front()) break;
    
    //check for border case where training converges to a two node tree with a statistically significant mean offset but the same median
    //in which case significance cutoff never converges to the single node tree case we checked above
    //if (tree.Responses().size()==2) {
    //  int numidenticaltrees = 1;
    //  int jtree = itree-1;
    //  int bestvar = tree.CutIndices().front();
    //  float cutval = tree.CutVals().front();
    //  while (jtree>=0 && forest->Trees().at(jtree).Responses().size()==2 && forest->Trees().at(jtree).CutIndices().front()==bestvar && forest->Trees().at(jtree).CutVals().front()==cutval) {
	//++numidenticaltrees;
	//--jtree;
    //  }
    //  printf("numidenticaltrees = %i\n",numidenticaltrees);
    //  if ( numidenticaltrees > (1.0/fShrinkage) ) break;
    //}
    //else if (tree.Responses().size()==3) {
    //    int numidenticaltrees = 1;
    //    int jtree = itree-1;
    //    int bestvar1 = tree.CutIndices()[0];
    //    float cutval1 = tree.CutVals()[0];
    //    int bestvar2 = tree.CutIndices()[1];
    //    float cutval2 = tree.CutVals()[1];
    //    while (jtree>=0 && forest->Trees().at(jtree).Responses().size()==3 && 
    //            forest->Trees().at(jtree).CutIndices()[0]==bestvar1 && forest->Trees().at(jtree).CutVals()[0]==cutval1 &&
    //            forest->Trees().at(jtree).CutIndices()[1]==bestvar2 && forest->Trees().at(jtree).CutVals()[1]==cutval2
    //            ) {
    //        ++numidenticaltrees;
    //        --jtree;
    //    }
    //  printf("numidenticaltrees = %i\n",numidenticaltrees);
    //  if ( numidenticaltrees > (1.0/fShrinkage) ) break;
    //}

    int jtree = itree - 1;
    int numidenticaltrees = 1;
    while (jtree>=0 && forest->Trees().at(jtree).CutIndices().size()==tree.CutIndices().size())
    {
        bool diff = false;
        for(unsigned int cutid=0; cutid<tree.CutIndices().size(); cutid++)
        {
            if(forest->Trees().at(jtree).CutIndices()[cutid]!=tree.CutIndices()[cutid] &&
                    forest->Trees().at(jtree).CutVals()[cutid]!=tree.CutVals()[cutid]
              )
            {
                diff = true;
                break;
            }
        }
        if(!diff)
        {
            ++numidenticaltrees;
            --jtree;
        }
        else
            break;
    }
    if(numidenticaltrees>1)
        printf("numidenticaltrees = %i\n",numidenticaltrees);
    if ( numidenticaltrees > (1.0/fShrinkage) ) break;


    //recompute transition point and transformed target
    std::sort(evts.begin(),evts.end(),GBRAbsTargetCMP());
    transumw = 0.;
    transit=evts.begin();
    while(transumw<(fTransitionQuantile*sumw) && transit!=evts.end()) {
      transumw += (*transit)->Weight();
      transition = std::abs((*transit)->Target());  
      ++transit;
    } 
    
    for (std::vector<GBREvent*>::iterator it=evts.begin(); it!=evts.end(); ++it) {
      double tgt = (*it)->Target();
      if (std::abs(tgt)<transition) (*it)->SetTransTarget(tgt);
      else if (tgt<0.) (*it)->SetTransTarget(-transition);
      else (*it)->SetTransTarget(transition);
    }       
    
  }

  // clean all events when training is done
    for (std::vector<GBREvent*>::iterator it=evts.begin(); it!=evts.end(); ++it) {
        delete *it;
    }
    evts.clear();
  
  //return fully trained GBRForest
  return forest;
  
}

//_______________________________________________________________________
void GBRTrainer::TrainTree(const std::vector<GBREvent*> &evts, double sumwtotal, GBRTree &tree, int nvars, double transition) {
  //index of current intermediate node
  int thisidx = tree.CutIndices().size();    
  
  //number of events input to node
  int nev = evts.size();
  
  //index of best cut variable
  int bestvar = 0;

  //float *__restrict *__restrict__ sumws = _sumws;  
  
  //trivial open-mp based multithreading of loop over input variables
  //The loop is thread safe since each iteration writes into its own
  //elements of the 2-d arrays
#pragma omp parallel for
  for (int ivar=0; ivar<nvars; ++ivar) {
    
     //int thread_id = omp_get_thread_num(); 

    //fill temporary array of quantiles (to allow auto-vectorization of later loops)
    for (int iev = 0; iev<nev; ++iev) {
      _quants[ivar][iev] = evts[iev]->Quantile(ivar);
    }
    
    int minquant = std::numeric_limits<int>::max();
    int maxquant = 0;
    
    //find max and min quantiles in the input events
    //(this loop should be vectorized by gcc with reasonable optimization options)
    for (int iev = 0; iev<nev; ++iev) {
      if (_quants[ivar][iev]<minquant) minquant = _quants[ivar][iev];
      if (_quants[ivar][iev]>maxquant) maxquant = _quants[ivar][iev];
    }    
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

    //zero arrays where necessary and compute map between bin numbers
    //and variable cut values
    //This loop should auto-vectorize in appropriate compiler/settings
    for (unsigned int ibin=0; ibin<nbins; ++ibin) {
      _ws[ivar][ibin] = 0.;
      _ws2[ivar][ibin] = 0.;
      _ns[ivar][ibin] = 0;
      _tgts[ivar][ibin] = 0.;
      _tgt2s[ivar][ibin] = 0.;
      
      int quant = ((1+ibin) << pscale) + offset - 1;
      if (quant>=fNQuantiles) quant = fNQuantiles-1;
      
      _varvals[ivar][ibin] = fQuantileMaps[ivar][quant];

    }
    //compute reduced bin value for each event using bit-shift operations
    //This loop should auto-vectorize in appropriate compiler/settings
    for (int iev=0;iev<nev;++iev) {
      _bins[ivar][iev] = (_quants[ivar][iev]-offset)>>pscale;
    }

    //compute summed quantities differential in each bin
    //(filling 'histograms')
    //This loop is one of the most expensive in the algorithm for large training samples
    //This loop can unfortunately not be vectorized because the memory addressed 
    //are computed within the loop iteration
    //JOSH: Is this already fundamentally making vectorization impossible because the addresses to be incremented are
    //scattered, or is it just that the compiler can't resolve the dependencies?  If the latter, can we force gcc to vectorize
    //this loop)
    
    for (int iev=0;iev<nev;++iev) {
      int ibin = _bins[ivar][iev];
      
      _ws[ivar][ibin] += evts[iev]->Weight();
      _ws2[ivar][ibin] += evts[iev]->Weight()*evts[iev]->Weight();
      ++_ns[ivar][ibin];
      _tgts[ivar][ibin] += evts[iev]->WeightedTransTarget();
      _tgt2s[ivar][ibin] += evts[iev]->WeightedTransTarget2();

    } 
    //convert differential arrays to cumulative arrays by summing over
    //each element
    //loop cannot be vectorized because this is an iterative calculation
    _sumws[ivar][0] = _ws[ivar][0];
    _sumws2[ivar][0] = _ws2[ivar][0];
    _sumns[ivar][0] = _ns[ivar][0];
    _sumtgts[ivar][0] = _tgts[ivar][0];
    _sumtgt2s[ivar][0] = _tgt2s[ivar][0];    
    
    for (unsigned int ibin=1; ibin<nbins; ++ibin) {      
      _sumws[ivar][ibin] = _sumws[ivar][ibin-1] + _ws[ivar][ibin];
      _sumws2[ivar][ibin] = _sumws2[ivar][ibin-1] + _ws2[ivar][ibin];
      _sumns[ivar][ibin] = _sumns[ivar][ibin-1] + _ns[ivar][ibin];
      _sumtgts[ivar][ibin] = _sumtgts[ivar][ibin-1] + _tgts[ivar][ibin];
      _sumtgt2s[ivar][ibin] = _sumtgt2s[ivar][ibin-1] + _tgt2s[ivar][ibin];  
    }
    //int n = sumns[ivar][nbins-1];
    float sumw = _sumws[ivar][nbins-1];
    float sumw2 = _sumws2[ivar][nbins-1];
    float sumtgt = _sumtgts[ivar][nbins-1];
    float sumtgt2 = _sumtgt2s[ivar][nbins-1];      
    
    //weighted variance of target in full dataset
    float fullvariance = sumtgt2 - sumtgt*sumtgt/sumw;
    //    float fullvariancevar = fullvariance*fullvariance/sumw2/sumw2;
    
    _fullvars[ivar] = fullvariance;
    
   // printf("fullrms = %5f, sumtgt2 = %5f, sumtgt = %5f, sumw = %5f\n",fullrms,sumtgt2,sumtgt,sumw);
    
    //printf("short loop\n");
    float maxsepgain = -std::numeric_limits<float>::max();
    float cutval = 0.;
    int nleft= 0;
    int nright = 0;
    float sumwleft=0.;
    float sumwright=0.;
    int bestbin=0;
    //loop over all bins and compute improvement in weighted variance of target for each split
    //This loop is relatively expensive and should auto-vectorize in the appropriate compiler/settings
    for (unsigned int ibin=0; ibin<nbins; ++ibin) {      
      float leftvariance = _sumtgt2s[ivar][ibin] - _sumtgts[ivar][ibin]*_sumtgts[ivar][ibin]/_sumws[ivar][ibin];
      //float leftvariancevar = leftvariance*leftvariance/sumws2[ivar][ibin]/sumws2[ivar][ibin];
       
      float rightsumw = sumw - _sumws[ivar][ibin];
      float rightsumw2 = sumw2 - _sumws2[ivar][ibin];
      float righttgtsum = sumtgt - _sumtgts[ivar][ibin];
      float righttgt2sum = sumtgt2 - _sumtgt2s[ivar][ibin];
      float rightvariance = righttgt2sum - righttgtsum*righttgtsum/rightsumw;
      //float rightvariancevar = rightvariance*rightvariance/rightsumw2/rightsumw2;
      
      
      //weighted improvement in variance from this split
      _bsepgains[ivar][ibin] = fullvariance - rightvariance - leftvariance;
      //bsepgainsigs[ivar][ibin] = bsepgains[ivar][ibin]/sqrt(leftvariancevar+rightvariancevar+fullvariancevar);
      //_bsepgainsigs[ivar][ibin] = sqrt((_sumtgts[ivar][ibin]/_sumws[ivar][ibin] - righttgtsum/rightsumw)*(_sumtgts[ivar][ibin]/_sumws[ivar][ibin] - righttgtsum/rightsumw)/(leftvariance/_sumws[ivar][ibin]/_sumws2[ivar][ibin] + rightvariance/rightsumw/rightsumw2));
      _bsepgainsigs[ivar][ibin] = sqrt((_sumtgts[ivar][ibin]/_sumws[ivar][ibin] - righttgtsum/rightsumw)*(_sumtgts[ivar][ibin]/_sumws[ivar][ibin] - righttgtsum/rightsumw)/(leftvariance/_sumws[ivar][ibin]/_sumws[ivar][ibin]/_sumws[ivar][ibin]*_sumws2[ivar][ibin] + rightvariance/rightsumw/rightsumw/rightsumw*rightsumw2));
    }
    //loop over computed variance improvements and select best split, respecting also minimum number of events per node
    //This loop cannot auto-vectorize, at least in gcc 4.6x due to the mixed type conditional, but it's relatively fast
    //in any case
    for (unsigned int ibin=0; ibin<nbins; ++ibin) {   
        //printf("sumnleft = %i, sumnright = %i, sepgain = %5f, sepgainsig = %5f\n",_sumns[ivar][ibin], (nev-_sumns[ivar][ibin]),_bsepgains[ivar][ibin],_bsepgainsigs[ivar][ibin]);
      //if (sumns[ivar][ibin]>=fMinEvents && (nev-sumns[ivar][ibin])>=fMinEvents && bsepgains[ivar][ibin]>maxsepgain) {
	if (_sumns[ivar][ibin]>=fMinEvents && (nev-_sumns[ivar][ibin])>=fMinEvents && _bsepgains[ivar][ibin]>maxsepgain && _bsepgainsigs[ivar][ibin]>fMinCutSignificance) {
	maxsepgain = _bsepgains[ivar][ibin];
        bestbin = ibin;
      }
    }
    cutval = _varvals[ivar][bestbin];
    nleft = _sumns[ivar][bestbin];
    nright = nev - nleft;
    sumwleft = _sumws[ivar][bestbin];
    sumwright = sumw - sumwleft;        
    
    _sepgains[ivar] = maxsepgain;
    _sepgainsigs[ivar] = _bsepgainsigs[ivar][bestbin];
    _cutvals[ivar] = cutval;
    _nlefts[ivar] = nleft;
    _nrights[ivar] = nright;
    _sumwlefts[ivar] = sumwleft;
    _sumwrights[ivar] = sumwright;
    _sumtgtlefts[ivar] = _sumtgts[ivar][bestbin];
    _sumtgtrights[ivar] = sumtgt - _sumtgts[ivar][bestbin];
    _leftvars[ivar] = _sumtgt2s[ivar][bestbin] - _sumtgts[ivar][bestbin]*_sumtgts[ivar][bestbin]/_sumws[ivar][bestbin];
    _rightvars[ivar] = (sumtgt2-_sumtgt2s[ivar][bestbin]) - (sumtgt-_sumtgts[ivar][bestbin])*(sumtgt-_sumtgts[ivar][bestbin])/(sumw-_sumws[ivar][bestbin]);
    _bestbins[ivar] = bestbin;
  }
  

  
  float globalsepgain = -std::numeric_limits<float>::max();
  for (int ivar=0; ivar<nvars; ++ivar) {
    if (_sepgains[ivar]>globalsepgain) {
      globalsepgain = _sepgains[ivar];
      bestvar = ivar;
    }
  }    

  //printf("Best var = %i\n", bestvar);
  
  //if no appropriate split found, make this node terminal
  if (globalsepgain<=0.) {
    //no valid split found, making this node a leaf node
    //printf("thisidx = %i, globalsepgain = %5f, no valid split\n",thisidx, globalsepgain);
    tree.CutIndices().push_back(0);
    tree.CutVals().push_back(0);
    //tree.CutSepGains().push_back(0);
    tree.LeftIndices().push_back(0);   
    tree.RightIndices().push_back(0);    
    
    tree.RightIndices()[thisidx] = -tree.Responses().size();
    tree.LeftIndices()[thisidx] = -tree.Responses().size();
    
    BuildLeaf(evts,sumwtotal,tree,transition);
    return;
  }
  
  //fill vectors of event pointers for left and right nodes below this one
  std::vector<GBREvent*> leftevts;
  std::vector<GBREvent*> rightevts;
  
  leftevts.reserve(nev);
  rightevts.reserve(nev);
  
  int nleft = 0;
  int nright = 0;
  double sumwleft = 0.;
  double sumwright = 0.;
  
  for (std::vector<GBREvent*>::const_iterator it = evts.begin(); it!=evts.end(); ++it) {
    if ((*it)->Var(bestvar)>_cutvals[bestvar]) {
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
 
  //  float fullres = sqrt(_fullvars[bestvar]/sumwtotal);
  // float leftres = sqrt(_leftvars[bestvar]/sumwleft);
  //float rightres = sqrt(_rightvars[bestvar]/sumwright);
 
  // float fullmean = (_sumtgtlefts[bestvar] + _sumtgtrights[bestvar])/sumwtotal;
  //float leftmean = _sumtgtlefts[bestvar]/sumwleft;
  //float rightmean = _sumtgtrights[bestvar]/sumwright;
  
  
  //printf("thisidx = %i, bestvar = %i, cutval = %5f, n = %i, nleft = %i, nright = %i, fullres = %5f, leftres = %5f, rightres = %5f, fullmean = %5f, leftmean = %5f, rightmrean = %5f, leftsepgain = %5f, sepgainsig = %5f\n",thisidx,bestvar,_cutvals[bestvar],nev,_nlefts[bestvar],_nrights[bestvar],fullres,leftres,rightres,fullmean, leftmean, rightmean, _sepgains[bestvar],_sepgainsigs[bestvar]);
  

  assert(_nlefts[bestvar]==nleft);
  assert(_nrights[bestvar]==nright);
  
  //fill intermediate node
  tree.CutIndices().push_back(bestvar);
  tree.CutVals().push_back(_cutvals[bestvar]);
  //tree.CutSepGains().push_back(_sepgains[bestvar]);
  tree.LeftIndices().push_back(0);   
  tree.RightIndices().push_back(0);  
  
  //check if left node is terminal
  bool termleft = nleft<=(2*fMinEvents);
  if (termleft) tree.LeftIndices()[thisidx] = -tree.Responses().size();
  else tree.LeftIndices()[thisidx] = tree.CutIndices().size();
  
  //printf("this idx = %i, termleft = %i, nleft = %i, fMinEvents = %i\n",thisidx,  termleft,nleft,fMinEvents);  
  
  //build left node as appropriate
  if (termleft) {  
    BuildLeaf(leftevts,sumwleft,tree,transition);
  }
  else {  
      //printf("Training left node: sumwleft=%5f, nvars=%i, transition=%5f\n", sumwleft, nvars, transition);
    TrainTree(leftevts,sumwleft,tree,nvars,transition);  
  }
  
  //check if right node is terminal
  bool termright = nright<=(2*fMinEvents);
  if (termright) tree.RightIndices()[thisidx] = -tree.Responses().size();
  else tree.RightIndices()[thisidx] = tree.CutIndices().size();
    
  //printf("this idx = %i, termright = %i, nright = %i, fMinEvents = %i\n",thisidx,  termright,nright,fMinEvents);    
  
  //build right node as appropriate
  if (termright) {  
    BuildLeaf(rightevts,sumwright,tree,transition);
  }
  else {  
      //printf("Training right node: sumwleft=%5f, nvars=%i, transition=%5f\n", sumwleft, nvars, transition);
    TrainTree(rightevts,sumwright,tree,nvars,transition);  
  }
  
}

  
  


//_______________________________________________________________________
void GBRTrainer::BuildLeaf(const std::vector<GBREvent*> &evts, double sumw, GBRTree &tree, double transition) {

  //printf("building leaf\n");
  
  //  int thisidx = -tree.Responses().size();
  //printf("thisidx = %i\n",thisidx);
  
  float medsumw = 0;
  float median = 0.;  
  for (std::vector<GBREvent*>::const_iterator it = evts.begin(); it!=evts.end(); ++it) {    
    if (medsumw<(0.5*sumw)) {
      median = (*it)->Target();
      medsumw += (*it)->Weight();
    }
    else break;
  }

  
  float shift = 0.;
  const float invsumw = 1.0/sumw;
  for (std::vector<GBREvent*>::const_iterator it = evts.begin(); it!=evts.end(); ++it) {
    float weight = (*it)->Weight();
    float diff = (*it)->Target() - median;
    
    if (std::abs(diff) > transition) {
      if (diff<0.) diff = -transition;
      else diff = transition;
    }
    
    shift += weight*invsumw*diff; 
  
  }
  
  float response = fShrinkage*(median+shift);
  //float response = fShrinkage*median;


  tree.Responses().push_back(response);
  
  for (std::vector<GBREvent*>::const_iterator it = evts.begin(); it!=evts.end(); ++it) {
    (*it)->SetTarget((*it)->Target()-response);
  }
  
  //printf("thisidx = %i, n = %i, response = %5f, median=%5f, shift=%5f \n", thisidx, int(evts.size()) ,response, median, shift);
  
}


