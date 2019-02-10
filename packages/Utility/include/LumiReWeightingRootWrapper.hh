#ifndef LUMIREWEIGHTINGROOTWRAPPER
#define LUMIREWEIGHTINGROOTWRAPPER

//hides the big nasty boost from the poor little cint, poor thing it was destraught when boost told it what proper c++ coding and memory management looks like

#include <string>

namespace edm {
class Lumi3DReWeighting;
}

class LumiReWeightingRootWrapper {

private:
  edm::Lumi3DReWeighting* reWeighter_;

  //I fear to attempt to copy Lumi3DReWeighting, in theory it should copy fine but...
private:
  LumiReWeightingRootWrapper(const LumiReWeightingRootWrapper&){}
  LumiReWeightingRootWrapper& operator=(const LumiReWeightingRootWrapper& ){return *this;}
    
public:
  LumiReWeightingRootWrapper( const std::string& generatedFile,
				const std::string& dataFile,
				const std::string& genHistName,
				const std::string& dataHistName);
  ~LumiReWeightingRootWrapper();
  
 

  
  void weight3D_init(float scale);
  double weight3D(int nrPUNeg,int nrPU,int nrPUPos);
  
};

#endif
