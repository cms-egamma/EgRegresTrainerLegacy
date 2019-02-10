#include "Utility/JsonMaker.hh"



void JsonMaker::fill(int runnr,int lumiSec)
{
  std::map<int,std::set<int> >::iterator mapIt;
  mapIt = data_.insert(std::make_pair(runnr,std::set<int>())).first;
  mapIt->second.insert(lumiSec);
  
}

std::ostream &operator <<(std::ostream& output,const JsonMaker &jsonMaker){return jsonMaker.print(output);}

std::ostream& JsonMaker::print(std::ostream& output)const
{
  typedef std::map<int,std::set<int> >::const_iterator MapIt;

  
  output <<"{";

  for(MapIt mapIt = data_.begin();mapIt!=data_.end();++mapIt){
    if(mapIt!=data_.begin()) output <<", "; //seperating the entries, so not needed for first one
    

    std::vector<std::pair<int,int> > lumiSecRanges;
    fillLumiSecRangeVec_(lumiSecRanges,mapIt->second);
    

    output <<"\""<<mapIt->first<<"\": [";
    for(size_t rangeNr=0;rangeNr<lumiSecRanges.size();rangeNr++){
      output  <<"["<<lumiSecRanges[rangeNr].first<<", "<<lumiSecRanges[rangeNr].second<<"]";
      if(rangeNr+1!=lumiSecRanges.size()) output <<", ";
    }
    output <<"]";
  }
  
  output <<" }";

  return output;
}

void JsonMaker::printRunLumiRanges()const
{
  typedef std::map<int,std::set<int> >::const_iterator MapIt;

  for(MapIt mapIt = data_.begin();mapIt!=data_.end();++mapIt){ 

    std::vector<std::pair<int,int> > lumiSecRanges;
    fillLumiSecRangeVec_(lumiSecRanges,mapIt->second);
    
    for(size_t rangeNr=0;rangeNr<lumiSecRanges.size();rangeNr++){
      std::cout  <<mapIt->first<<" "<<lumiSecRanges[rangeNr].first<<" "<<lumiSecRanges[rangeNr].second<<std::endl;
     
    }
   
  }

  
}

void JsonMaker::printRunLumi()const
{
  typedef std::map<int,std::set<int> >::const_iterator MapIt;

  
  //  output <<"{";

  for(MapIt mapIt = data_.begin();mapIt!=data_.end();++mapIt){
    //  if(mapIt!=data_.begin()) output <<", "; //seperating the entries, so not needed for first one
    

    //std::vector<std::pair<int,int> > lumiSecRanges;
    //fillLumiSecRangeVec_(lumiSecRanges,mapIt->second);
   
    std::set<int>::const_iterator setIt;
    for(setIt=mapIt->second.begin();setIt!=mapIt->second.end();++setIt) {
      std::cout <<mapIt->first<<" "<<*setIt<<std::endl;
    }
  }
  

}

   

void JsonMaker::fillLumiSecRangeVec_(std::vector<std::pair<int,int> >& lumiSecRanges,const std::set<int>& inputLumiSecs)
{
  std::set<int>::const_iterator setIt=inputLumiSecs.begin();
  
  int lastLumiSec =*setIt;
  int startOfRange= *setIt;
  for(++setIt;setIt!=inputLumiSecs.end();++setIt){
    if(*setIt-lastLumiSec>1) { //ie this is not the next sequential lumi sec
      lumiSecRanges.push_back(std::make_pair(startOfRange,lastLumiSec));
      startOfRange = *setIt;
      lastLumiSec = *setIt;
    }else{
      lastLumiSec = *setIt;
    }

  }
  lumiSecRanges.push_back(std::make_pair(startOfRange,lastLumiSec));

}
