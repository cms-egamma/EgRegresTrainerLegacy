
#include <vector>
#include <string>
  
class GoodLumiChecker{
  
public:
  struct GoodLumiData {
    int runnr;
    std::vector<std::pair<int,int> > allowedLumis;
    bool operator<(const GoodLumiData& rhs)const{return runnr<rhs.runnr;}
    bool operator<(const int rhsRunnr)const{return runnr<rhsRunnr;}
    bool operator==(const int rhsRunnr)const{return runnr==rhsRunnr;}
  };
  
  class GoodLumiDataComp {
    public:
    bool operator()(const GoodLumiData& lhs,const GoodLumiData& rhs)const{return lhs.runnr<rhs.runnr;}
    bool operator()(int lhs,const GoodLumiData& rhs)const{return lhs<rhs.runnr;}
    bool operator()(int lhs,int rhs)const{return lhs<rhs;}
    bool operator()(const GoodLumiData& lhs,int rhs)const{return lhs.runnr<rhs;}
  };
private:
  std::vector<GoodLumiData> goodLumis_; //sorted by runnr, each run number can only enter once
  mutable int indexOfLastRun_; //temporary cache hence mutable
public:
  GoodLumiChecker(const std::string& goodLumiList="");
  ~GoodLumiChecker(){}
  bool isGoodLumiSec(int runnr,int lumiSec)const;
  void setGoodLumiList(const std::string& goodLumiList);
  void clear(){goodLumis_.clear();indexOfLastRun_=-1;}
private:
  void addLumiRange_(int runnr,int lowerLumi,int upperLumi);
};
