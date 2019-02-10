#ifndef EVENTLISTCOMP
#define EVENTLISTCOMP

//#include "TObject.h"

#include <iostream>
#include <vector>
#include <fstream>

class EventListComp {

public:
  struct EventListCompData{
    int runnr;
    int lumiSec;
    int eventnr;
    double userVar;
  
    bool operator<(const EventListCompData& rhs)const;
  };

private:
  std::vector<EventListComp::EventListCompData> eventListA1_; //I dont like having varibles differ by a letter hence A1 instead of 1
  std::vector<EventListComp::EventListCompData> eventListB2_;

public:
  EventListComp(){}
  virtual ~EventListComp(){}
 
  void compareLists(const char* filenameA1,const char* filenameB2);

  static void readList_(const char* filename,std::vector<EventListComp::EventListCompData>& eventList);

  // ClassDef(EventListComp,1)

};


std::ostream &operator <<(std::ostream& output,const EventListComp::EventListCompData &data);

#endif
