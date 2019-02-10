#ifndef UTILITY_ROOTBOOST_H
#define UTILITY_ROOTBOOST_H

//root and boost sometimes dont get along, this is for root scripts
//nothing else should use this, just use boost directly

#include "TObject.h"
#include<string>
class RootBoost {
 
private:
  RootBoost(){}
  virtual ~RootBoost(){}
public:
  static void replace_all(std::string& input,const std::string& search, 
			  const std::string & format);
  ClassDef(RootBoost,1) 
};



#endif
