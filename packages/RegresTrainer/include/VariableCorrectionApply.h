

#ifndef VARIABLECORRECTIONAPPLY_H
#define VARIABLECORRECTIONAPPLY_H

#include <vector>
#include <string>

class TTree;
class ICorrection;

class VariableCorrectionApply
{

    public:

        VariableCorrectionApply();
        ~VariableCorrectionApply();

        TTree *ApplyAsFriend(TTree *intree,
                const std::string& varToCorrect,
                std::string targetname,
                const ICorrection* correction) const;

};

#endif
