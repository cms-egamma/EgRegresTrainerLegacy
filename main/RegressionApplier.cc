#include "CondFormats/EgammaObjects/interface/GBRForest.h"
#include "CondFormats/EgammaObjects/interface/GBRForestD.h"

#include "Utility/HistFuncs.hh"
#include "Utility/CmdLineInt.hh"

#include "TFile.h"
#include "TTree.h"
#include "TROOT.h"
#include "RooWorkspace.h"

#include <vdt/vdtMath.h>

#include <cmath>
#include <thread>
#include <future>
  

struct TreeData{
public:
  struct EvtData{
    int runnr,lumiSec,eventnr;
    static std::string contents(){return "runnr/I:lumiSec:eventnr";}
  };
  EvtData evt;
  float mean;
  float sigma;
  float invTar;
    
  void createBranches(TTree* tree){
    tree->Branch("evt",&evt,evt.contents().c_str());
    tree->Branch("mean",&mean,"mean/F");
    tree->Branch("sigma",&sigma,"sigma/F");
    tree->Branch("invTar",&invTar,"invTar/F");
  }
  
};

std::pair<double,double> getRes(const std::vector<float>& regData,const std::pair<const GBRForestD*,const GBRForestD*>& gbrForest);

class RegFunc {
public:
  std::pair<double,double> res;

  void operator()(const std::vector<float>& regData,
		const std::pair<const GBRForestD*,const GBRForestD*>& gbrForest){
    res = getRes(regData,gbrForest);
  }
};


void fillTree(const std::vector<std::vector<float> >& regDataEB,
	      const std::vector<std::vector<float> >& regDataEE,
	      const std::vector<std::vector<float> >& evtData,
	      const std::pair<const GBRForestD*,const GBRForestD*>& gbrForestEB,
	      const std::pair<const GBRForestD*,const GBRForestD*>& gbrForestEE,
	      TreeData& outTreeData,
	      TTree* outTree,
	      unsigned int nrThreads);



int main(int argc, char** argv)
{
  char inFilename[1024];
  char outFilename[1024];
  char gbrFilenameEB[1024];
  char gbrFilenameEE[1024];
  char treeName[1024];
  int nrThreads;
  bool writeFullTree;
  CmdLineInt cmdLineInt(argv[0]);
  cmdLineInt.addNonOption(inFilename,true," ","input files");
  cmdLineInt.addNonOption(outFilename,true," ","output filename");
  cmdLineInt.addOption("gbrForestFileEB",gbrFilenameEB,"test.root","gbrForestFile for barrel");
  cmdLineInt.addOption("gbrForestFileEE",gbrFilenameEE,"test.root","gbrForestFile for endcap");
  cmdLineInt.addOption("nrThreads",&nrThreads,1,"number of threads for reading tree");
  cmdLineInt.addOption("treeName",treeName,"egRegTree"," name of the tree");
  cmdLineInt.addOption("writeFullTree",&writeFullTree,false," writes the full tree to file");
  if(!cmdLineInt.processCmdLine(argc,argv)) return 0; //exit if we havnt managed to get required parameters
  //this appears to do very little...
  if(nrThreads>1){
    ROOT::EnableImplicitMT(nrThreads);
  }

  TTree* inTree = HistFuncs::makeChain(treeName,inFilename);
  TFile* outFile = new TFile(outFilename,"RECREATE");
  TTree* outTree = new TTree((std::string(treeName)+"Friend").c_str(),"");
    
  TreeData outTreeData;
  outTreeData.createBranches(outTree);

  TFile* gbrFileEB = TFile::Open(gbrFilenameEB,"READ");
  TFile* gbrFileEE = TFile::Open(gbrFilenameEE,"READ");
  auto readVarList = [](TFile* file,const std::string& listname){
    auto vars = reinterpret_cast<std::vector<std::string>* >(file->Get(listname.c_str()));
    if(!vars){
      std::cout <<"could not var list "<<listname<<" "<<vars<<" "<<file->Get(listname.c_str())<<std::endl;
      return std::string("");
    }
    std::string varsStr;
    for(const auto& var : *vars){
      if(!varsStr.empty()) varsStr+=":";
      varsStr+=var;
    }
    return varsStr;
  };
  const std::string varsEB = readVarList(gbrFileEB,"varlistEB");
  const std::string varsEE = readVarList(gbrFileEE,"varlistEE");
  const std::string* targetEB = reinterpret_cast<std::string*>(gbrFileEB->Get("targetEB"));
  const std::string* targetEE = reinterpret_cast<std::string*>(gbrFileEE->Get("targetEE"));
 
  GBRForestD* gbrMeanEB = reinterpret_cast<GBRForestD*>(gbrFileEB->Get("EBCorrection"));
  GBRForestD* gbrSigmaEB = reinterpret_cast<GBRForestD*>(gbrFileEB->Get("EBUncertainty"));
  GBRForestD* gbrMeanEE = reinterpret_cast<GBRForestD*>(gbrFileEE->Get("EECorrection"));
  GBRForestD* gbrSigmaEE = reinterpret_cast<GBRForestD*>(gbrFileEE->Get("EEUncertainty"));
  if(!gbrMeanEB || !gbrSigmaEB){
    std::cout <<" error, couldnt get GBRForests "<<gbrMeanEB<<" "<<gbrSigmaEB<<" from "<<gbrFilenameEB<<" for barrel"<<std::endl;
  }
  if(!gbrMeanEE || !gbrSigmaEE){
    std::cout <<" error, couldnt get GBRForests "<<gbrMeanEE<<" "<<gbrSigmaEE<<" from "<<gbrFilenameEE<<" for endcap"<<std::endl;
  }
  if(varsEB.empty() || varsEE.empty()){
    std::cout <<"vars not found, exiting"<<std::endl;
    outFile->Write();
    return 0;
  }

  const auto regDataEBAll = HistFuncs::readTree(inTree,varsEB+":"+*targetEB,"");
  const auto regDataEEAll = HistFuncs::readTree(inTree,varsEE+":"+*targetEE,"");
  const auto evtData = HistFuncs::readTree(inTree,"runnr:eventnr:lumiSec:sc.isEB","");
  fillTree(regDataEBAll,regDataEEAll,evtData,{gbrMeanEB,gbrSigmaEB},{gbrMeanEE,gbrMeanEE},
	   outTreeData,outTree,nrThreads);


  outFile->Write();

  if(writeFullTree){
    outFile->cd();
    TTree* fullTree = inTree->CloneTree();
    TBranch* fullMeanBranch=fullTree->Branch("regMean",&outTreeData.mean);
    TBranch* fullSigmaBranch=fullTree->Branch("regSigma",&outTreeData.sigma);
    TBranch* fullInvTarBranch=fullTree->Branch("regInvTar",&outTreeData.invTar);
  
    int nrEntries = outTree->GetEntries();
    std::cout <<"entryies "<<nrEntries<<std::endl;
    for(int entryNr=0;entryNr<nrEntries;entryNr++){
      outTree->GetEntry(entryNr);
      fullMeanBranch->Fill();
      fullSigmaBranch->Fill();
      fullInvTarBranch->Fill();
    }
  }
  outFile->Write();
  return 0;
}

void fillTree(const std::vector<std::vector<float> >& regDataEB,
	      const std::vector<std::vector<float> >& regDataEE,
	      const std::vector<std::vector<float> >& evtData,
	      const std::pair<const GBRForestD*,const GBRForestD*>& gbrForestEB,
	      const std::pair<const GBRForestD*,const GBRForestD*>& gbrForestEE,
	      TreeData& outTreeData,
	      TTree* outTree,
	      unsigned int nrThreads)
{
  AnaFuncs::Timer timer;
  using RegResult = std::future<std::pair<double,double> >;
  std::vector<RegResult> threads(nrThreads);
  
  std::vector<RegFunc> regFuncs(nrThreads);
  std::vector<unsigned int> threadEntryNrs;
  for(unsigned int i=0;i<nrThreads;i++) threadEntryNrs.push_back(i);
  
  auto initThread = [&threads,&regDataEB,&regDataEE,
		     &gbrForestEB,gbrForestEE,
		     &regFuncs,&evtData,&threadEntryNrs](unsigned int threadNr){
    const auto entryNr = threadEntryNrs[threadNr];
    if(entryNr<regDataEB.size()){
      auto& regData = evtData[entryNr][3] ? regDataEB[entryNr] : regDataEE[entryNr];
      auto& gbrForest = evtData[entryNr][3] ? gbrForestEB : gbrForestEE;
      threads[threadNr] = std::async(getRes,regData,gbrForest);
    }
  };
  
  for(size_t threadNr=0;threadNr<threads.size();threadNr++) initThread(threadNr);
  

  bool anyThreadActive = true;
  while(anyThreadActive){
    anyThreadActive = false;
    for(unsigned int threadNr=0;threadNr<nrThreads;threadNr++){
      if(threads[threadNr].valid()){
	anyThreadActive = true;
	
	auto regResult = threads[threadNr].get();
	
	auto entryNr = threadEntryNrs[threadNr];
	if(entryNr%100000==0) std::cout <<"entry "<<entryNr<< "/" <<regDataEB.size()<<" time "<<timer<<std::endl;
   
	const auto& evtDataEntry = evtData[entryNr];
	outTreeData.evt.runnr=evtDataEntry[0];
	outTreeData.evt.eventnr=evtDataEntry[1];
	outTreeData.evt.lumiSec=evtDataEntry[2];
	outTreeData.mean = regResult.first;
	outTreeData.sigma = regResult.second;
	outTreeData.invTar = evtDataEntry[3] ? 1./regDataEB[entryNr].back() : 1./regDataEE[entryNr].back();
	threadEntryNrs[threadNr]+=nrThreads;
	initThread(threadNr);

	outTree->Fill();
	
      }
    }
  }
}

std::pair<double,double> getRes(const std::vector<float>& regData,const std::pair<const GBRForestD*,const GBRForestD*>& gbrForest)
{
  
  //(These should be stored inside the conditions object in the future as well)
  constexpr double meanlimlow  = 0.2;
  constexpr double meanlimhigh = 2.0;
  constexpr double meanoffset  = meanlimlow + 0.5*(meanlimhigh-meanlimlow);
  constexpr double meanscale   = 0.5*(meanlimhigh-meanlimlow);
  
  constexpr double sigmalimlow  = 0.0002;
  constexpr double sigmalimhigh = 0.5;
  constexpr double sigmaoffset  = sigmalimlow + 0.5*(sigmalimhigh-sigmalimlow);
  constexpr double sigmascale   = 0.5*(sigmalimhigh-sigmalimlow);  
  
  
  //these are the actual BDT responses
  double rawmean = gbrForest.first->GetResponse(regData.data());
  double rawsigma = gbrForest.second->GetResponse(regData.data());
  
  //apply transformation to limited output range (matching the training)
  double mean = meanoffset + meanscale*vdt::fast_sin(rawmean);
  double sigma = sigmaoffset + sigmascale*vdt::fast_sin(rawsigma);
  return {mean,sigma};
}

