#include "Utility/HistFuncs.hh"
#include "TFile.h"
#include "TTree.h"
#include <iostream>


void addRegCorrToTree(const std::string& treeName,const std::string& inFilename,const std::string& corrFilename,const std::string& outFilename)
{
  
  TTree* inTree = HistFuncs::makeChain(treeName,inFilename);
  TTree* corrTree = HistFuncs::makeChain(treeName,corrFilename);
  
  TFile* outFile = TFile::Open(outFilename.c_str(),"RECREATE");
  TTree* outTree = inTree->CloneTree();

  float mean,sigma,invTar;
  corrTree->SetBranchAddress("mean", &mean);
  corrTree->SetBranchAddress("sigma", &sigma);
  corrTree->SetBranchAddress("invTar", &invTar);
  
  TBranch* outMeanBranch=outTree->Branch("regMean",&mean);
  TBranch* outSigmaBranch=outTree->Branch("regSigma",&sigma);
  TBranch* outInvTarBranch=outTree->Branch("regInvTar",&invTar);
  
  int nrEntries = inTree->GetEntries();
  for(int entryNr=0;entryNr<nrEntries;entryNr++){
    corrTree->GetEntry(entryNr);
    outMeanBranch->Fill();
    outSigmaBranch->Fill();
    outInvTarBranch->Fill();
  }
  
  outFile->Write();
  delete outFile;
}
