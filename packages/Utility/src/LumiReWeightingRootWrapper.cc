#include "Utility/LumiReWeightingRootWrapper.hh"

#include "Utility/Lumi3DReWeighting.hh"

LumiReWeightingRootWrapper::LumiReWeightingRootWrapper(const std::string& generatedFile,
						       const std::string& dataFile,
						       const std::string& genHistName,
						       const std::string& dataHistName):
  reWeighter_(new edm::Lumi3DReWeighting(generatedFile,dataFile,genHistName,dataHistName))
{

}

LumiReWeightingRootWrapper::~LumiReWeightingRootWrapper()
{
  delete reWeighter_;
}

void LumiReWeightingRootWrapper::weight3D_init(float scale)
{
  reWeighter_->weight3D_init(scale);
}

double LumiReWeightingRootWrapper::weight3D(int nrPUNeg,int nrPU,int nrPUPos)
{
  return reWeighter_->weight3D(nrPUNeg,nrPU,nrPUPos);
}

