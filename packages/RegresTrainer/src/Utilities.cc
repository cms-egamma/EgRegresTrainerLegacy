/**
 *  @file  Utilities.cxx
 *  @brief  
 *
 *
 *  @author  Jean-Baptiste Sauvan <sauvan@llr.in2p3.fr>
 *
 *  @date    03/27/2010
 *
 *  @internal
 *     Created :  03/27/2010
 * Last update :  03/27/2010 01:00:52 PM
 *          by :  JB Sauvan
 *
 * =====================================================================================
 */


#include "RegresTrainer/Utilities.h"


//--- STL

using namespace std;




///*****************************************************************/
//void tokenize(const string& str,
//        vector<string>& tokens,
//        const string& delimiters)
///*****************************************************************/
//{
//    // Skip delimiters at beginning.
//    string::size_type lastPos = str.find_first_not_of(delimiters, 0);
//    // Find first "non-delimiter".
//    string::size_type pos     = str.find_first_of(delimiters, lastPos);
//
//    while (string::npos != pos || string::npos != lastPos)
//    {
//        // Found a token, add it to the vector.
//        tokens.push_back(str.substr(lastPos, pos - lastPos));
//        // Skip delimiters.  Note the "not_of"
//        lastPos = str.find_first_not_of(delimiters, pos);
//        // Find next "non-delimiter"
//        pos = str.find_first_of(delimiters, lastPos);
//    }
//}

/*****************************************************************/
void tokenize(const string& str,
        vector<string>& tokens,
        const string& delimiter)
/*****************************************************************/
{
    string::size_type length = delimiter.size();
    string::size_type lastPos = 0;
    string::size_type pos     = str.find(delimiter, 0);


    while (string::npos != pos)
    {
        // Found a token, add it to the vector.
        if(str.substr(lastPos, pos - lastPos).size()>0)
            tokens.push_back(str.substr(lastPos, pos - lastPos));
        lastPos = pos + length;
        // Find next "non-delimiter"
        pos = str.find(delimiter, lastPos);
    }
    if(str.substr(lastPos).size()>0)
        tokens.push_back(str.substr(lastPos));
}



/*****************************************************************/
string intToString(int n)
/*****************************************************************/
{
    ostringstream oss;
    oss  <<  n;
    return oss.str();
}



/*****************************************************************/
void findAndReplace(string& sInput, string sFind, string sReplace )
/*****************************************************************/
{
    size_t itPos = 0; 
    size_t itFindLen = sFind.length();
    size_t itReplaceLen = sReplace.length();

    if( itFindLen == 0 )
        return;

    while( (itPos = sInput.find( sFind, itPos )) != std::string::npos )
    {
        sInput.replace( itPos, itFindLen, sReplace );
        itPos += itReplaceLen;
    }

}


/*****************************************************************/
void strip(std::string& sInput)
/*****************************************************************/
{
	//-- removing blanks at the beginning and at the end
	while(*sInput.begin()==' ' || *sInput.begin()=='\t') sInput.erase(sInput.begin());
	while(*(sInput.end()-1)==' ' || *(sInput.end()-1)=='\t') sInput.erase(sInput.end()-1);
}




