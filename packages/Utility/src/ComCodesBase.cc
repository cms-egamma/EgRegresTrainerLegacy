#include "Utility/ComCodesBase.hh"


void ComCodesBase::_setCode(const char* descript,int code)
{
  _codeDefs[descript] = code;
}


//multiple descriptions seperated by ":" 
int ComCodesBase::getCode(const char* descript)const
{ 
  //first copy the character string to a local array so we can manipulate it
  char localDescript[256];
  strcpy(localDescript,descript);

  int code = 0x0000; 
  char* codeKey = strtok(localDescript,":");
  std::map<std::string,int> ::const_iterator mapIt;
  while(codeKey!=NULL){
    mapIt = _codeDefs.find(codeKey);
    if(mapIt!=_codeDefs.end()) code |= mapIt->second;
    else std::cout<<"ComCodesBase::getCode : Error, Key "<<codeKey<<" not found"<<std::endl;
    codeKey = strtok(NULL,":"); //getting new substring
    
  }
  return code;
}

void ComCodesBase::getCodeName(int code,std::string& id)const
{
  id.clear();
  
  std::map<std::string,int> ::const_iterator mapIt;
  for(mapIt=_codeDefs.begin();mapIt!=_codeDefs.end();++mapIt){
    if((code&mapIt->second)==mapIt->second){
      if(!id.empty()) id+=":"; //seperating entries by a ':' 
      id+=mapIt->first;
    }
  }
}
