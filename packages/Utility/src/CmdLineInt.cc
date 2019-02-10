#include "Utility/CmdLineInt.hh"
#include <cstdlib>
CmdLineInt::CmdLineInt(const char* cmdName):
  _nrRequired(0),_cmdName(cmdName)
{

}

void CmdLineInt::addOption(const char* option,double *arg,double defaultArg,const char* desc)
{
  OptionArgPair optArgPair(option,arg,defaultArg,desc);
  _optionArgList.push_back(optArgPair);
}

void CmdLineInt::addOption(const char* option,int *arg,int defaultArg,const char* desc)
{
  OptionArgPair optArgPair(option,arg,defaultArg,desc);
  _optionArgList.push_back(optArgPair);
}

void CmdLineInt::addOption(const char* option,char *arg,const char* defaultArg,const char* desc)
{
  OptionArgPair optArgPair(option,arg,defaultArg,desc);
  _optionArgList.push_back(optArgPair);
}

void CmdLineInt::addOption(const char* option,bool *arg,bool defaultArg,const char* desc)
{
  OptionArgPair optArgPair(option,arg,defaultArg,desc);
  _optionArgList.push_back(optArgPair);
}

void CmdLineInt::addOption(const char* option,std::string *arg,const char* defaultArg,const char* desc)
{
  OptionArgPair optArgPair(option,arg,defaultArg,desc);
  _optionArgList.push_back(optArgPair);
}

void CmdLineInt::addNonOption(double *arg,bool required,double defaultArg,const char* desc)
{
  OptionArgPair optArgPair("null",arg,defaultArg,desc);
  _nonOptionArgList.push_back(optArgPair);
   if(required){
    _nrRequired++;
    if(_nrRequired!=(int)_nonOptionArgList.size()){
      std::cout<<"CmdLineInt::addNonOption : Error Required arguement added after non requred arguement"<<std::endl;
    }
   }
}


void CmdLineInt::addNonOption(int *arg,bool required,int defaultArg,const char* desc)
{
  OptionArgPair optArgPair("null",arg,defaultArg,desc);
  _nonOptionArgList.push_back(optArgPair);
   if(required){
    _nrRequired++;
    if(_nrRequired!=(int)_nonOptionArgList.size()){
      std::cout<<"CmdLineInt::addNonOption : Error Required arguement added after non requred arguement"<<std::endl;
    }
   }
}

void CmdLineInt::addNonOption(char *arg,bool required,const char* defaultArg,const char* desc)
{
  OptionArgPair optArgPair("null",arg,defaultArg,desc);
  _nonOptionArgList.push_back(optArgPair);
   if(required){
    _nrRequired++;
    if(_nrRequired!=(int)_nonOptionArgList.size()){
      std::cout<<"CmdLineInt::addNonOption : Error Required arguement added after non requred arguement"<<std::endl;
    }
   }
}

void CmdLineInt::addNonOption(bool *arg,bool required,bool defaultArg,const char* desc)
{
  OptionArgPair optArgPair("null",arg,defaultArg,desc);
  _nonOptionArgList.push_back(optArgPair);
   if(required){
    _nrRequired++;
    if(_nrRequired!=(int)_nonOptionArgList.size()){
      std::cout<<"CmdLineInt::addNonOption : Error Required arguement added after non requred arguement"<<std::endl;
    }
   }
}

void CmdLineInt::addNonOption(std::string *arg,bool required,const char* defaultArg,const char* desc)
{
  OptionArgPair optArgPair("null",arg,defaultArg,desc);
  _nonOptionArgList.push_back(optArgPair);
   if(required){
    _nrRequired++;
    if(_nrRequired!=(int)_nonOptionArgList.size()){
      std::cout<<"CmdLineInt::addNonOption : Error Required arguement added after non requred arguement"<<std::endl;
    }
   }
}

//I humbly applogise for my youngerselfs use of the variables i and j
//theres a small bug which 
bool CmdLineInt::processCmdLine(int argc,char* argv[])
{
  int nrNonOptionArg=0;
  for(int i=1;i<argc;i++){ 
    if(argv[i][0]=='-'){ //option arguement
      std::string option;
      //check if single letter or string option
      if(argv[i][1]!='-'){
	if(strlen(&argv[i][1])!=1){ //if single letter, check that its really a single letter and somebody didnt do -option instead of --option
	  std::cout <<"error: functionality to string multiple options together does not currently exist (yet), sorry"<<std::endl;
	  std::cout <<"you may have meant to type --"<<&argv[i][1]<<" but you actually typed -"<<&argv[i][1]<<std::endl;
	  return false;
	}
	option = argv[i][1];
      }
      else option = &argv[i][2]; 
      
      //check if the help option is selected and end if it is
      if(option=="h" || option=="help"){
	_printHelp();
	return false;
      }

      i++; //advancing to the next index which is the arguement
      bool found=false;
      for(int j=0;j<(int)_optionArgList.size() &&!found;j++){
	//finish if find the correct option
	if(_optionArgList[j].processOption(option.c_str(),argv[i])) found=true;
      }
      if(!found){
	std::cout <<"error: option \""<<option<<"\" not found "<<std::endl;
	_printHelp();
	return false;
      }
    }else{ //non optional arguement
      _nonOptionArgList[nrNonOptionArg].processOption("null",argv[i]);
      nrNonOptionArg++;
    }
  }

  if(nrNonOptionArg<_nrRequired){
    std::cout <<"insufficient parameters"<<std::endl;
    _printHelp();
    return false;
  }
  return true;
}


void CmdLineInt::_printHelp()
{
  std::cout <<"Useage : ";
  std::cout <<_cmdName.c_str();
  for(int i=0;i<(int)_nonOptionArgList.size();i++){
   if(i==_nrRequired) std::cout <<"  [optional]";
   std::cout <<"  "<<_nonOptionArgList[i].descript();
    
  }
  

  std::cout <<std::endl<<"Options : "<<std::endl;
  for(int i=0;i<(int)_optionArgList.size();i++){
    if(_optionArgList[i].singleCharOpt()) std::cout <<"  -";
    else std::cout <<" --";
    std::cout <<_optionArgList[i].option()<<" "<<_optionArgList[i].descript()
	      <<" (default value = "<<_optionArgList[i].defaultValue()<<")"<<std::endl;
  }
}

CmdLineInt::OptionArgPair::OptionArgPair(const char* option,double *arg,double defaultArg,const char* desc):
  _type("double"),_arg((void*)arg),_option(option),_desc(desc)
{
  *arg = defaultArg;
  std::ostringstream tempDefault;
  tempDefault <<defaultArg;
  _defaultValue = tempDefault.str();
}

CmdLineInt::OptionArgPair::OptionArgPair(const char* option,int *arg,int defaultArg,const char* desc):
  _type("int"),_arg((void*)arg),_option(option),_desc(desc)
{
  *arg = defaultArg;
  std::ostringstream tempDefault;
  tempDefault <<defaultArg;
  _defaultValue = tempDefault.str();
}

CmdLineInt::OptionArgPair::OptionArgPair(const char* option,char *arg,const char* defaultArg,const char* desc):
  _type("charString"),_arg((void*)arg),_option(option),_desc(desc)
{
  strcpy(arg,defaultArg);
  _defaultValue = defaultArg;
}

CmdLineInt::OptionArgPair::OptionArgPair(const char* option,bool *arg,bool defaultArg,const char* desc):
  _type("bool"),_arg((void*)arg),_option(option),_desc(desc)
{
  *arg = defaultArg;
  std::ostringstream tempDefault;
  tempDefault <<defaultArg;
  _defaultValue = tempDefault.str();
}

CmdLineInt::OptionArgPair::OptionArgPair(const char* option,std::string *arg,const char* defaultArg,const char* desc):
  _type("string"),_arg((void*)arg),_option(option),_desc(desc)
{
  *arg = defaultArg;
  std::ostringstream tempDefault;
  tempDefault <<defaultArg;
  _defaultValue = tempDefault.str();
}

CmdLineInt::OptionArgPair::OptionArgPair(const OptionArgPair& rhs):
  _type(rhs._type),_arg(rhs._arg),_option(rhs._option),_desc(rhs._desc),_defaultValue(rhs._defaultValue)
{
  
}

bool CmdLineInt::OptionArgPair::processOption(const char* option,const char* arg)
{
  if(_option==option){
    if(_type=="double"){
      double *argDouble = (double*) _arg;
      *argDouble = atof(arg);
    }else if(_type=="int"){
      int *argInt = (int*) _arg;
      *argInt = atoi(arg);
    }else if(_type=="charString"){
      char *argChar = (char*) _arg;
      strcpy(argChar,arg);
    }else if(_type=="bool"){
      bool *argBool = (bool*) _arg;
      *argBool = atoi(arg);
    }else if(_type=="string"){
      std::string* argString= (std::string*) _arg;
      *argString = arg;
    }else{
      std::cout <<"CmdLineInt::OptionArgPair::processOption : Error : type "<<_type<<" not recognised"<<std::endl;
    }
    return true;
  }else return false;
}
