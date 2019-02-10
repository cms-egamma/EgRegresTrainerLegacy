#include "Utility/RootBoost.hh"

#ifndef __CINT__

#include <boost/algorithm/string/replace.hpp>
void RootBoost::replace_all(std::string& input,const std::string& search, 
			     const std::string & format)
{
  boost::replace_all(input,search,format);
}
#include <exception>
namespace boost {
  void throw_exception(std::exception const&){}
  
}


#endif

ClassImp(RootBoost)
