#ifndef COMCODESBASE
#define COMCODESBASE

#include <cstring>
#include <map>
#include <string>
#include <iostream>

//class is the base comunication code class which all the other bitwise commication
//classes inherit from
class ComCodesBase { 

private:
  std::map<std::string,int> _codeDefs;

protected:
  ComCodesBase(){} //this class shouldnt be instanced,only those inheriting from it
  ~ComCodesBase(){} 
  
  void _setCode(const char *descript,int code);

public:
  int getCode(const char *descript)const;
  void getCodeName(int code,std::string& id)const;
};

#endif
