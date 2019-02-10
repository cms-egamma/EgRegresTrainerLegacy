#ifndef EVENTID_H
#define EVENTID_H

//sort by dataset, then run , then lumi then eventnr
struct EventId {
  
  int runnr;
  int lumiSec;
  int eventnr;
  int datasetCode;
  EventId():runnr(0),lumiSec(0),eventnr(0),datasetCode(-1){}
  EventId(int iRunnr,int iLumiSec,int iEventnr,int iDataset=-1):runnr(iRunnr),lumiSec(iLumiSec),eventnr(iEventnr),datasetCode(iDataset){}

  static std::string contents(){return "runnr/I:lumiSec:eventnr:datasetCode";}
  
  bool operator<(const EventId& rhs)const{
    if(datasetCode<rhs.datasetCode) return true;
    else if(datasetCode>rhs.datasetCode) return false;
    else { //same dataset
      if(runnr<rhs.runnr) return true;
      else if(runnr>rhs.runnr) return false;
      else{ //same run
	if(lumiSec<rhs.lumiSec) return true;
	else if(lumiSec>rhs.lumiSec) return false;
	else if(eventnr<rhs.eventnr) return true;
	else return false;
      }
    }
  }

  bool operator==(const EventId& rhs)const{
    return datasetCode==rhs.datasetCode && runnr==rhs.runnr && lumiSec==rhs.lumiSec && eventnr==rhs.eventnr;
  }
  std::ostream& print(std::ostream& output)const{output<<" event: "<<runnr<<" "<<lumiSec<<" "<<eventnr<<" "<<datasetCode;return output;}
};

std::ostream& operator <<(std::ostream& output,const EventId& obj){return obj.print(output);}
#endif
