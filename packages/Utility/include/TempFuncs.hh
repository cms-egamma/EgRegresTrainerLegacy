#ifndef TEMPFUNCS
#define TEMPFUNCS

//I appreciate this is an unforunate name. This is a collection of TEMPLATE functions not TEMPORARY functions which I find usefull
#include "Utility/AnaFuncs.hh"

#include "TClonesArray.h"

namespace TempFuncs {
  template<class T> std::string str(const T& val);
  template<class T> void copyTClonesArray(TClonesArray& lhs,const TClonesArray& rhs);
  template<class T> bool makeVec(const std::string& string,std::vector<T>& vec);
  template<class T> std::vector<T> convertStringToVec(const std::string& string){std::vector<T> vec;makeVec(string,vec); return vec;}
  template<class T,class Y,class CompType> struct PairComp : public std::binary_function<std::pair<T,Y>,std::pair<T,Y>,bool> {
    bool operator()(const std::pair<T,Y>& lhs,const std::pair<T,Y>& rhs,const CompType& compObj = CompType()){return compObj(lhs.first,rhs.first);} 
    bool operator()(const std::pair<T,Y>& lhs,const T& rhs,const CompType& compObj = CompType()){return compObj(lhs.first,rhs);}
    bool operator()(const T& lhs,const  std::pair<T,Y>& rhs,const CompType& compObj = CompType()){return compObj(lhs,rhs.first);}
  
  };

  template <class T1,class T2,typename Comp=std::less<T1> > struct PairSortBy1st : public std::binary_function<std::pair<T1,T2>,std::pair<T1,T2>,bool> { 
    Comp comp;
    bool operator()(const std::pair<T1,T2>& lhs,const std::pair<T1,T2>&rhs){return comp(lhs.first,rhs.first);}
    bool operator()(const T1& lhs,const std::pair<T1,T2>&rhs){return comp(lhs,rhs.first);}
    bool operator()(const std::pair<T1,T2>& lhs,const T1 &rhs){return comp(lhs.first,rhs);}
    bool operator()(const T1& lhs,const T1 &rhs){return comp(lhs,rhs);}
   

  };



  template <class T1,class T2,typename Comp=std::less<T2> > struct PairSortBy2nd : public std::binary_function<std::pair<T1,T2>,std::pair<T1,T2>,bool>  {
    Comp comp;
    bool operator()(const std::pair<T1,T2>& lhs,const std::pair<T1,T2>&rhs){return comp(lhs.second,rhs.second);}
    bool operator()(const T2& lhs,const std::pair<T1,T2>&rhs){return comp(lhs,rhs.second);}
    bool operator()(const std::pair<T1,T2>& lhs,const T2 &rhs){return comp(lhs.second,rhs);}
    bool operator()(const T2& lhs,const T2 &rhs){return comp(lhs,rhs);}
   
  };

   template <class T1,class T2> struct PairEq1st : public std::unary_function<std::pair<T1,T2>,bool> { 
     
     T1 val;
     PairEq1st(const T1& iVal):val(iVal){}
     bool operator()(const std::pair<T1,T2>& pair){return pair.first==val;}
     bool operator==(const std::pair<T1,T2>& pair){return pair.first==val;}
 
  };

   template <class T1,class T2> struct PairEq2nd : public std::unary_function<std::pair<T1,T2>,bool> { 
     T2 val;
     PairEq2nd(const T2& iVal):val(iVal){}
     bool operator()(const std::pair<T1,T2>& pair){return pair.second==val;} 
     bool operator==(const std::pair<T1,T2>& pair){return pair.second==val;}
 
  };

  template <typename TCont,typename TData,typename TVal> const TData* findSingleSorted(const TCont& cont,const TVal& val);


}

template<class T> void TempFuncs::copyTClonesArray(TClonesArray& lhs,const TClonesArray& rhs)
{
  lhs.Delete();
  for(int i=0;i<rhs.GetLast()+1;i++){
    T* obj = (T*) rhs[i];
    new(lhs[i]) T(*obj);
  }
}

template<class T> bool TempFuncs::makeVec(const std::string& string,std::vector<T>& vec)
{
  std::vector<std::string> splitValues;
  AnaFuncs::splitStrings(string.c_str(),splitValues,":");
  vec.clear();
  vec.reserve(splitValues.size());
  for(size_t valNr=0;valNr<splitValues.size();valNr++){
    std::istringstream iss(splitValues[valNr]);
    T value;
    if((iss >> value).fail()) {
      std::cout <<" TempFuncs::Error failed to convert "<<splitValues[valNr]<<std::endl;
      return false;
    }
    vec.push_back(value);
  }
  return true;

} 

template<class T> std::string TempFuncs::str(const T& val)
{
  std::ostringstream ss;
  ss << val;
  return ss.str();
}
  

template <typename TCont,typename TData,typename TVal> const TData* TempFuncs::findSingleSorted(const TCont& container,const TVal& val)
{
  // typename TCont::const_iterator;
  std::pair<typename TCont::const_iterator,typename TCont::const_iterator> result = std::equal_range(container.begin(),container.end(),val);
  int nrFound = result.second-result.first;
  if(nrFound==1) return &*result.first;
  else if(nrFound>=1) {
    std::cout <<"warning "<<nrFound<<" instead of 1 "<<std::endl;
    return &*result.first;
  }
  return 0;
}

#endif
