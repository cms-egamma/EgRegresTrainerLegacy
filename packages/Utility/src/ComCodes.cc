#include "Utility/ComCodes.hh"

#include <algorithm>

void ComCodes::setCode(const std::string& descript,int code)
{
  bool found=false;
  for(size_t i=0;i<_codeDefs.size() && !found;i++){
    if(_codeDefs[i].first.compare(descript)==0) found=true;
  }
  if(!found) _codeDefs.push_back(std::pair<std::string,int>(descript,code));
 
  //_codeDefs[descript] = code;
}


// //multiple descriptions seperated by ":" 
// int ComCodes::getCode(const char* descript)const
// { 
//   //first copy the character string to a local array so we can manipulate it
//   char localDescript[256];
//   strcpy(localDescript,descript);

//   int code = 0x0000; 
//   char* codeKey = strtok(localDescript,":");
//   std::map<std::string,int> ::const_iterator mapIt;
//   while(codeKey!=NULL){
//     mapIt = _codeDefs.find(codeKey);
//     if(mapIt!=_codeDefs.end()) code |= mapIt->second;
//     else std::cout<<"ComCodes::getCode : Error, Key "<<codeKey<<" not found"<<std::endl;
//     codeKey = strtok(NULL,":"); //getting new substring
    
//   }
//   return code;
// }

int ComCodes::getCode(const std::string& descript)const
{ 
  //first copy the character string to a local array so we can manipulate it
  char localDescript[256];
  strcpy(localDescript,descript.c_str());
  
  int code = 0x0000; 
  char* codeKey = strtok(localDescript,":");
  //  std::map<std::string,int> ::const_iterator mapIt;
  while(codeKey!=NULL){
    bool found=false;

    for(size_t i=0;i<_codeDefs.size() && !found;i++){
      if(_codeDefs[i].first.compare(codeKey)==0){
 	found=true;
 	code |= _codeDefs[i].second;

       }
    }
   
    if(!found)  std::cout<<"ComCodes::getCode : Error, Key "<<codeKey<<" not found"<<std::endl;
    codeKey = strtok(NULL,":"); //getting new substring
    
  }
  return code;
}

bool ComCodes::keyComp(const std::pair<std::string,int>& lhs,const std::pair<std::string,int>& rhs)
{
  return lhs.first < rhs.first;
}

void ComCodes::getCodeName(int code,std::string& id)const
{
  id.clear();
  for(size_t i=0;i<_codeDefs.size();i++){ 
    if((code&_codeDefs[i].second)==_codeDefs[i].second){
      if(!id.empty()) id+=":";//seperating entries by a ':'
      id+=_codeDefs[i].first;
    }
    //   std::cout <<std::hex<<"code "<<code<<" cut "<<_codeDefs[i].second<<" name "<<_codeDefs[i].first<<" id "<<id.c_str()<<std::endl;
  }
 //  std::map<std::string,int> ::const_iterator mapIt;
//   for(mapIt=_codeDefs.begin();mapIt!=_codeDefs.end();++mapIt){
//     if((code&mapIt->second)==mapIt->second){
//       if(!id.empty()) id+=":"; //seperating entries by a ':' 
//       id+=mapIt->first;
//     }
//   }
}
