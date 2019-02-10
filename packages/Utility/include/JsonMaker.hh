#ifndef JSONMAKER
#define JSONMAKER

//anticipated that there will be a mix of read/writes, thus suggesting a map/set is apppropriate

#include <map>
#include <set>
#include <iostream>
#include <vector>

class JsonMaker {
private:
  std::map<int,std::set<int> > data_; //first int is runnr, second is set of lumiSections
public:
  void fill(int runnr,int lumiSec);
    
  
  friend std::ostream &operator <<(std::ostream& output,const JsonMaker &jsonMaker);
  std::ostream &print(std::ostream& output)const;
private:
  static void fillLumiSecRangeVec_(std::vector<std::pair<int,int> >& lumiSecRanges,const std::set<int>& inputLumiSecs);
  void printRunLumi()const;
  void printRunLumiRanges()const;
};


#endif

