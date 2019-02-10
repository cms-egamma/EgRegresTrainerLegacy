#ifndef CMDLINEINT
#define CMDLINEINT

//Class which handles the command line interface to the program

#include <string>
#include <cstring>
#include <iostream>
#include <vector>
#include <sstream>


//templates, I've heard of them... :)

class CmdLineInt {
private:

  class OptionArgPair{
  private:
    std::string _type;
    void *_arg; //little naughty....
    std::string _option;
    std::string _desc;
    std::string _defaultValue; //a string of the default option

  public:
    OptionArgPair(const char* option,double *arg,double defaultArg,const char* desc=NULL);
    OptionArgPair(const char* option,int *arg,int defaultArg,const char* desc=NULL);
    OptionArgPair(const char* option,char* arg,const char* defaultArg,const char* desc=NULL);
    OptionArgPair(const char* option,bool* arg,bool defaultArg,const char* desc=NULL);
    OptionArgPair(const char* option,std::string* arg,const char* defaultArg,const char* desc=NULL);
    OptionArgPair(const OptionArgPair& rhs);

    ~OptionArgPair(){} //dont own anything
    
    bool processOption(const char* option,const char* arg);

    const char* descript()const{return _desc.c_str();}
    const char* option()const{return _option.c_str();}
    const char* defaultValue()const{return _defaultValue.c_str();}
    bool singleCharOpt()const{return _option.length()==1;}

  };

private:
  std::vector<OptionArgPair> _optionArgList;
  std::vector<OptionArgPair> _nonOptionArgList;
  int _nrRequired;
  std::string _cmdName;

public:
  CmdLineInt(const char* cmdName=" ");
  ~CmdLineInt(){}

  void addOption(const char* option,double *arg,double defaultArg,const char* desc=NULL);
  void addOption(const char* option,int *arg,int defaultArg,const char* desc=NULL);
  void addOption(const char* option,char* arg,const char* defaultArg,const char* desc=NULL);
  void addOption(const char* option,bool* arg,bool defaultArg,const char* desc=NULL);
  void addOption(const char* option,std::string* arg,const char* defaultArg,const char* desc=NULL);

  void addNonOption(double *arg,bool required=true,double defaultArg=0.,const char* desc=NULL);
  void addNonOption(int *arg,bool required=true,int defaultArg=0,const char* desc=NULL);
  void addNonOption(char *arg,bool required=true,const char* defaultArg="",const char* desc=NULL);
  void addNonOption(bool *arg,bool required=true,bool defaultArg=false,const char* desc=NULL);
  void addNonOption(std::string *arg,bool required=true,const char* defaultArg="",const char* desc=NULL);

  bool processCmdLine(int argc,char* argv[]);
  
private:
  void _printHelp();

};

#endif
  



