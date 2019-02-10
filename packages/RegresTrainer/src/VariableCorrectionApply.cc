#include "RegresTrainer/VariableCorrectionApply.h"
#include "RegresTrainer/ICorrection.h"

#include "TTree.h"
#include "TTreeFormula.h"

#include <assert.h>
#include <iostream>

//_______________________________________________________________________
VariableCorrectionApply::VariableCorrectionApply()
{

}

//_______________________________________________________________________
VariableCorrectionApply::~VariableCorrectionApply() 
{

}



//_______________________________________________________________________
TTree *VariableCorrectionApply::ApplyAsFriend(TTree *intree, 
        const std::string& varToCorrect,
        std::string targetname,
        const ICorrection* correction) const
{

    //initialize TTreeFormula to read variable to be corrected
    TTreeFormula* varToCorrectForm = new TTreeFormula("varToCorrect",varToCorrect.c_str(),intree);

    // initialize TTreeFormulas to read input variables for correction
    std::vector<TTreeFormula*> inputformsCorr;
    if(correction)
    {
        const std::vector<std::string>& inputsCorr = correction->inputNames();
        std::vector<std::string>::const_iterator iti = inputsCorr.begin();
        std::vector<std::string>::const_iterator itiE = inputsCorr.end();
        for(;iti!=itiE;++iti)
            inputformsCorr.push_back(new TTreeFormula(iti->c_str(),iti->c_str(),intree));
    }

    Float_t target = 0.;
    std::vector<float> valsCorr(inputformsCorr.size());

    //initialize new friend tree
    TTree *friendtree = new TTree;
    friendtree->Branch(targetname.c_str(),&target,TString::Format("%s/F",targetname.c_str()));

    for (Long64_t iev=0; iev<intree->GetEntries(); ++iev) {
        target = 0.;
        if (iev%100000==0) printf("%i\n",int(iev));
        intree->LoadTree(iev);

        target = varToCorrectForm->EvalInstance();

        if(correction)
        {
            for (size_t i=0; i<inputformsCorr.size(); ++i)
            {
                valsCorr[i] = inputformsCorr[i]->EvalInstance();
                std::cout << inputformsCorr[i]->GetName() << "=" << valsCorr[i] << "\n";
            }
            std::cout << "Target before corr = " << target << "\n";
            target *= correction->call(valsCorr);
            std::cout << "Target after corr = " << target << "\n";
        }

        friendtree->Fill();

    }


    //clear TTreeFormulas
    delete varToCorrectForm;

    for (std::vector<TTreeFormula*>::const_iterator it = inputformsCorr.begin(); 
            it != inputformsCorr.end(); ++it) {
        delete *it;
    }


    intree->AddFriend(friendtree);

    //the branch addresses are set to local variables in this function
    //these local variables go out of scope after this function finishes
    //so we need to reset the branch addresses before returning
    friendtree->ResetBranchAddresses();
    return friendtree;

}


