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
  void setBranchAddresses(TTree* tree){
    tree->SetBranchAddress("evt",&evt);
    tree->SetBranchAddress("mean",&mean);
    tree->SetBranchAddress("sigma",&sigma);
    tree->SetBranchAddress("invTar",&invTar);
  }
};

class GBRForestCont{
private:
  std::pair<const GBRForestD*,const GBRForestD*> eb_;
  std::pair<const GBRForestD*,const GBRForestD*> ee_;
  std::pair<const GBRForestD*,const GBRForestD*> ebHighEt_;
  std::pair<const GBRForestD*,const GBRForestD*> eeHighEt_;
  double highEtThres_;

public:
  GBRForestCont(const std::string& ebFilename,const std::string& eeFilename,
		const std::string& ebHighEtFilename,const std::string& eeHighEtFilename,
		double highEtThres=std::numeric_limits<double>::max());
  GBRForestCont(const GBRForestCont&)=delete;
  GBRForestCont& operator=(const GBRForestCont&)=delete;
  ~GBRForestCont();
  
  const std::pair<const GBRForestD*,const GBRForestD*>& operator()(double et,bool isEB)const{
    if(et<highEtThres_) return isEB ? eb_ : ee_;
    else return isEB ? ebHighEt_ : eeHighEt_;
  }
private:
  static void loadForests(std::pair<const GBRForestD*,const GBRForestD*>& forests,const std::string& filename,const std::string& tag);
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
	      const GBRForestCont& gbrForests,
	      TreeData& outTreeData,
	      TTree* outTree,
	      unsigned int nrThreads);



int main(int argc, char** argv)
{
  char inFilename[1024];
  char outFilename[1024];
  char gbrFilenameEB[1024];
  char gbrFilenameEE[1024];
  char gbrFilenameEBHighEt[1024];
  char gbrFilenameEEHighEt[1024];
  char etBinVar[1024];
  char treeName[1024];
  char regOutTagChar[1024];
  int nrThreads;
  double highEtThres;
  bool writeFullTree;
  CmdLineInt cmdLineInt(argv[0]);
  cmdLineInt.addNonOption(inFilename,true," ","input files");
  cmdLineInt.addNonOption(outFilename,true," ","output filename");
  cmdLineInt.addOption("gbrForestFileEB",gbrFilenameEB,"test.root","gbrForestFile for barrel");
  cmdLineInt.addOption("gbrForestFileEE",gbrFilenameEE,"test.root","gbrForestFile for endcap");
  cmdLineInt.addOption("gbrForestFileEBHighEt",gbrFilenameEBHighEt,"","gbrForestFile for barrel high et, highEtThres must be set for this to be read");
  cmdLineInt.addOption("gbrForestFileEEHighEt",gbrFilenameEEHighEt,"","gbrForestFile for endcap high et, highEtThres must be set for this to be read");
  cmdLineInt.addOption("highEtThres",&highEtThres,std::numeric_limits<double>::max(),"threshold at which to apply the high Et forests");
  cmdLineInt.addOption("etBinVar",etBinVar,"(sc.rawEnergy+sc.rawESEnergy)*sin(2*atan(exp(-1*sc.scEta)))","et variable to bin vs");
  cmdLineInt.addOption("nrThreads",&nrThreads,1,"number of threads for reading tree");
  cmdLineInt.addOption("treeName",treeName,"egRegTree"," name of the tree");
  cmdLineInt.addOption("regOutTag",regOutTagChar,"","tag of the output regression branches , eg \"reg{regOutTagChar}Mean\" if writing full tree");
  cmdLineInt.addOption("writeFullTree",&writeFullTree,false," writes the full tree to file");
  if(!cmdLineInt.processCmdLine(argc,argv)) return 0; //exit if we havnt managed to get required parameters
 
  const std::string regOutTag(regOutTagChar);

  TTree* inTree = HistFuncs::makeChain(treeName,inFilename);
  TFile* outFile = new TFile(outFilename,"RECREATE");
  TTree* outTree = new TTree((std::string(treeName)+"Friend").c_str(),"");
    
  TreeData outTreeData;
  outTreeData.createBranches(outTree);

  GBRForestCont gbrForests(gbrFilenameEB,gbrFilenameEE,
			   gbrFilenameEBHighEt,gbrFilenameEEHighEt,
			   highEtThres);

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
  
  //now we get the target variable from the files
  //it rather assumes the high pt ones use the same variables
  //which currently is a good assumption
  TFile* gbrFileEB = TFile::Open(gbrFilenameEB,"READ");
  TFile* gbrFileEE = TFile::Open(gbrFilenameEE,"READ");

  const std::string varsEB = readVarList(gbrFileEB,"varlistEB");
  const std::string varsEE = readVarList(gbrFileEE,"varlistEE");
  const std::string* targetEB = reinterpret_cast<std::string*>(gbrFileEB->Get("targetEB"));
  const std::string* targetEE = reinterpret_cast<std::string*>(gbrFileEE->Get("targetEE"));
 
  if(varsEB.empty() || varsEE.empty()){
    std::cout <<"vars not found, exiting"<<std::endl;
    outFile->Write();
    return 0;
  }

  const auto regDataEBAll = HistFuncs::readTree(inTree,varsEB+":"+*targetEB,"");
  const auto regDataEEAll = HistFuncs::readTree(inTree,varsEE+":"+*targetEE,"");
  const auto evtData = HistFuncs::readTree(inTree,"runnr:eventnr:lumiSec:sc.isEB:"+std::string(etBinVar),"");
  fillTree(regDataEBAll,regDataEEAll,evtData,gbrForests,outTreeData,outTree,nrThreads);


  outFile->Write();
  delete outFile;
  delete gbrFileEB;
  delete gbrFileEE;

  if(writeFullTree){
    outFile = new TFile(outFilename,"UPDATE");
    outTree = static_cast<TTree*>(outFile->Get((std::string(treeName)+"Friend").c_str()));
    outTreeData.setBranchAddresses(outTree);
    TTree* fullTree = inTree->CloneTree();
    TBranch* fullMeanBranch=fullTree->Branch(("reg"+regOutTag+"Mean").c_str(),&outTreeData.mean);
    TBranch* fullSigmaBranch=fullTree->Branch(("reg"+regOutTag+"Sigma").c_str(),&outTreeData.sigma);
    TBranch* fullInvTarBranch=fullTree->Branch(("reg"+regOutTag+"InvTar").c_str(),&outTreeData.invTar);
  
    int nrEntries = outTree->GetEntries();
    for(int entryNr=0;entryNr<nrEntries;entryNr++){
      outTree->GetEntry(entryNr);
      fullMeanBranch->Fill();
      fullSigmaBranch->Fill();
      fullInvTarBranch->Fill();
    }
    outTree->Write(0,TObject::kOverwrite);
    outFile->Write();
  }

  return 0;
}

void fillTree(const std::vector<std::vector<float> >& regDataEB,
	      const std::vector<std::vector<float> >& regDataEE,
	      const std::vector<std::vector<float> >& evtData,
	      const GBRForestCont& gbrForests,
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
		     &gbrForests,
		     &regFuncs,&evtData,&threadEntryNrs](unsigned int threadNr){
    const auto entryNr = threadEntryNrs[threadNr];
    if(entryNr<regDataEB.size()){
      auto& regData = evtData[entryNr][3] ? regDataEB[entryNr] : regDataEE[entryNr];
      auto& gbrForest = gbrForests(evtData[entryNr][4],evtData[entryNr][3]);
      threads[threadNr] = std::async(getRes,regData,gbrForest);
      return true;
    }else return false;
  };

  unsigned int nrActiveThreads=0;
  for(size_t threadNr=0;threadNr<threads.size();threadNr++){
    if(initThread(threadNr)) nrActiveThreads = threadNr+1; //works as in increasing threadnr
  }

  //we want to read them in order so we wait for each thread to be done
  while(nrActiveThreads!=0){
    unsigned int newNrActiveThreads = 0;
    threads[0].wait();
 
    for(unsigned int threadNr=0;threadNr<nrActiveThreads;threadNr++){
      threads[threadNr].wait();	
      auto regResult = threads[threadNr].get();
	
      auto entryNr = threadEntryNrs[threadNr];
      if(entryNr%100000==0) std::cout <<"entry "<<entryNr<< "/" <<regDataEB.size()<<" time "<<timer<<std::endl;
      
      // if(!evtData[entryNr][3] && regResult.second>0.2 && regResult.second<0.3){
      //  	std::cout <<"dumping reg data "<<evtData[entryNr][1]<<" isEB "<<evtData[entryNr][3]<<std::endl;
      //  	for(size_t i=0;i<regDataEE[entryNr].size();i++){
      //  	  std::cout <<"   "<<i<<" "<<regDataEE[entryNr][i]<<std::endl;
      // 	}
      // }
      const auto& evtDataEntry = evtData[entryNr];
      outTreeData.evt.runnr=evtDataEntry[0];
      outTreeData.evt.eventnr=evtDataEntry[1];
      outTreeData.evt.lumiSec=evtDataEntry[2];
      outTreeData.mean = regResult.first;
      outTreeData.sigma = regResult.second;
      outTreeData.invTar = evtDataEntry[3] ? 1./regDataEB[entryNr].back() : 1./regDataEE[entryNr].back();
      outTree->Fill();	
      threadEntryNrs[threadNr]+=nrThreads;
      if(initThread(threadNr)) newNrActiveThreads = threadNr+1;
    }
    nrActiveThreads = newNrActiveThreads;
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

GBRForestCont::GBRForestCont(const std::string& ebFilename,
			     const std::string& eeFilename,
			     const std::string& ebHighEtFilename,
			     const std::string& eeHighEtFilename,
			     double highEtThres):
  eb_(nullptr,nullptr), ee_(nullptr,nullptr),
  ebHighEt_(nullptr,nullptr), eeHighEt_(nullptr,nullptr),
  highEtThres_(highEtThres)
{
  loadForests(eb_,ebFilename,"EB");
  loadForests(ee_,eeFilename,"EE");
  if(highEtThres_!=std::numeric_limits<double>::max()){
    loadForests(ebHighEt_,ebHighEtFilename,"EB");
    loadForests(eeHighEt_,eeHighEtFilename,"EE");
  }
}

void GBRForestCont::loadForests(std::pair<const GBRForestD*,const GBRForestD*>& forests,const std::string& filename,const std::string& tag)
{ 

  TFile* file = TFile::Open(filename.c_str(),"READ");
  forests.first = reinterpret_cast<GBRForestD*>(file->Get((tag+"Correction").c_str()));
  forests.second = reinterpret_cast<GBRForestD*>(file->Get((tag+"Uncertainty").c_str()));
  if(!forests.first  || !forests.second){
    std::cout <<" error, couldnt get GBRForests "<<forests.first<<" "<<forests.second<<" from "<<filename<<" for "<<tag<<std::endl;
  }
  delete file;
}

GBRForestCont::~GBRForestCont()
{
  auto del = [](std::pair<const GBRForestD*,const GBRForestD*>& obj){
    delete obj.first;
    delete obj.second;
  };
  del(eb_);
  del(ee_);
  del(ebHighEt_);
  del(eeHighEt_);
}
