#include "RegresTrainer/GBRApply.h"
#include "RegresTrainer/GBREvent.h"
#include "CondFormats/GBRForest/interface/GBRForest.h"
#include "CondFormats/GBRForest/interface/GBRForestD.h"
#include "RegresTrainer/ICorrection.h"


#include "TTree.h"
#include "TTreeFormula.h"

#include <assert.h>
#include <malloc.h>
#include <iostream>

//_______________________________________________________________________
GBRApply::GBRApply()
{

}

//_______________________________________________________________________
GBRApply::~GBRApply() 
{

}

//_______________________________________________________________________
//TTree *GBRApply::ApplyAsFriend(TTree *intree, const GBRForest *forest, const std::vector<std::string> &vars, std::string targetname) const
//{
//  
//  int nvars = vars.size();
//  
//  //initialize TTreeFormulas to read variables from TTree
//  std::vector<TTreeFormula*> inputforms;
//  for (std::vector<std::string>::const_iterator it = vars.begin(); 
//      it != vars.end(); ++it) {
//    inputforms.push_back(new TTreeFormula(it->c_str(),it->c_str(),intree));
//  }
//  
//  Float_t target = 0.;
//  Float_t *vals = new Float_t[nvars];
//  
//  //initialize new friend tree
//  TTree *friendtree = new TTree;
//  friendtree->Branch(targetname.c_str(),&target,TString::Format("%s/F",targetname.c_str()));
//  
//  for (Long64_t iev=0; iev<intree->GetEntries(); ++iev) {
//    //if (iev%100000==0) printf("%i\n",int(iev));
//    intree->LoadTree(iev);
//        
//    for (int i=0; i<nvars; ++i) {
//      vals[i] = inputforms[i]->EvalInstance();
//    }
//    
//    target = forest->GetResponse(vals);
//    
//    friendtree->Fill();
//
//  }
//  
//  
//  //clear TTreeFormulas
//  for (std::vector<TTreeFormula*>::const_iterator it = inputforms.begin(); 
//        it != inputforms.end(); ++it) {
//      delete *it;
//  }
//  
//  delete[] vals;
//  
//  intree->AddFriend(friendtree);
//
//  //the branch addresses are set to local variables in this function
//  //these local variables go out of scope after this function finishes
//  //so we need to reset the branch addresses before returning
//  friendtree->ResetBranchAddresses();
//  return friendtree;
//  
//}

//_______________________________________________________________________
TTree *GBRApply::ApplyAsFriend(TTree *intree, 
        const GBRForest *forestEB, const GBRForest *forestEE,
        const std::vector<std::string> &varsEB, const std::vector<std::string> &varsEE,
        const std::string& cutEB, const std::string& cutEE,
        std::string targetname,
        const ICorrection* correction) const
{

    int nvarsEB = varsEB.size();
    int nvarsEE = varsEE.size();

    // in EB or EE?
    TTreeFormula formIsEB(cutEB.c_str(), cutEB.c_str(), intree);
    TTreeFormula formIsEE(cutEE.c_str(), cutEE.c_str(), intree);

    //initialize TTreeFormulas to read variables from TTree
    std::vector<TTreeFormula*> inputformsEB;
    std::vector<TTreeFormula*> inputformsEE;
    for (std::vector<std::string>::const_iterator it = varsEB.begin(); 
            it != varsEB.end(); ++it) {
        inputformsEB.push_back(new TTreeFormula(it->c_str(),it->c_str(),intree));
    }
    for (std::vector<std::string>::const_iterator it = varsEE.begin(); 
            it != varsEE.end(); ++it) {
        inputformsEE.push_back(new TTreeFormula(it->c_str(),it->c_str(),intree));
    }

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
    Float_t *valsEB = new Float_t[nvarsEB];
    Float_t *valsEE = new Float_t[nvarsEE];
    std::vector<float> valsCorr(inputformsCorr.size());

    //initialize new friend tree
    TTree *friendtree = new TTree;
    friendtree->Branch(targetname.c_str(),&target,TString::Format("%s/F",targetname.c_str()));

    for (Long64_t iev=0; iev<intree->GetEntries(); ++iev) {
        target = 0.;
        if (iev%100000==0) printf("%i\n",int(iev));
        intree->LoadTree(iev);
        bool isEB = formIsEB.EvalInstance();
        bool isEE = formIsEE.EvalInstance();
        if(!isEB && !isEE)
            std::cout << "ERROR: GBRApply: not isEB nor isEE in event \n" << iev ;
        else if(isEB && isEE)
            std::cout << "ERROR: GBRApply: isEB and isEE in event \n" << iev ;
        
        if(isEB)
        {
            for (int i=0; i<nvarsEB; ++i) {
                valsEB[i] = inputformsEB[i]->EvalInstance();
            }
            target = forestEB->GetResponse(valsEB);
        }
        else if(isEE)
        {
            for (int i=0; i<nvarsEE; ++i) {
                valsEE[i] = inputformsEE[i]->EvalInstance();
            }
            target = forestEE->GetResponse(valsEE);
        }
        if(correction)
        {
            for (size_t i=0; i<inputformsCorr.size(); ++i)
            {
                valsCorr[i] = inputformsCorr[i]->EvalInstance();
                //std::cout << inputformsCorr[i]->GetName() << "=" << valsCorr[i] << "\n";
            }
            //std::cout << "Target before corr = " << target << "\n";
            target *= correction->call(valsCorr);
            //std::cout << "Target after corr = " << target << "\n";
        }

        friendtree->Fill();

    }


    //clear TTreeFormulas
    for (std::vector<TTreeFormula*>::const_iterator it = inputformsEB.begin(); 
            it != inputformsEB.end(); ++it) {
        delete *it;
    }
    for (std::vector<TTreeFormula*>::const_iterator it = inputformsEE.begin(); 
            it != inputformsEE.end(); ++it) {
        delete *it;
    }
    for (std::vector<TTreeFormula*>::const_iterator it = inputformsCorr.begin(); 
            it != inputformsCorr.end(); ++it) {
        delete *it;
    }

    delete[] valsEB;
    delete[] valsEE;

    intree->AddFriend(friendtree);

    //the branch addresses are set to local variables in this function
    //these local variables go out of scope after this function finishes
    //so we need to reset the branch addresses before returning
    friendtree->ResetBranchAddresses();
    return friendtree;

}

// //_______________________________________________________________________
// TTree *GBRApply::ApplyAsFriendTransform(TTree *intree, 
//         const GBRForestD *forestEB, const GBRForestD *forestEE,
//         const std::vector<std::string> &varsEB, const std::vector<std::string> &varsEE,
//         const std::string& cutEB, const std::string& cutEE,
//         std::string targetname,
//         double low, double high) const
// {
//     // transform parameters
//     double scale = 0.5*(high-low);
//     double offset = low + 0.5*(high-low);


//     int nvarsEB = varsEB.size();
//     int nvarsEE = varsEE.size();

//     // in EB or EE?
//     std::cout << "Making TTreeFormula for EB" << std::endl;
//     std::cout << "    cutEB = " << cutEB.c_str() << std::endl;
//     TTreeFormula formIsEB(cutEB.c_str(), cutEB.c_str(), intree);

//     std::cout << "Making TTreeFormula for EE" << std::endl;
//     std::cout << "    cutEE = " << cutEE.c_str() << std::endl;
//     TTreeFormula formIsEE(cutEE.c_str(), cutEE.c_str(), intree);

//     //initialize TTreeFormulas to read variables from TTree
//     std::vector<TTreeFormula*> inputformsEB;
//     std::vector<TTreeFormula*> inputformsEE;
//     for (std::vector<std::string>::const_iterator it = varsEB.begin(); 
//             it != varsEB.end(); ++it) {
//         inputformsEB.push_back(new TTreeFormula(it->c_str(),it->c_str(),intree));
//     }
//     for (std::vector<std::string>::const_iterator it = varsEE.begin(); 
//             it != varsEE.end(); ++it) {
//         inputformsEE.push_back(new TTreeFormula(it->c_str(),it->c_str(),intree));
//     }


//     Float_t target = 0.;
//     Float_t *valsEB = new Float_t[nvarsEB];
//     Float_t *valsEE = new Float_t[nvarsEE];

//     //initialize new friend tree
//     TTree *friendtree = new TTree;
//     friendtree->Branch(targetname.c_str(),&target,TString::Format("%s/F",targetname.c_str()));

//     for (Long64_t iev=0; iev<intree->GetEntries(); ++iev) {
//         target = 0.;
//         if (iev%100000==0) printf("%i\n",int(iev));
//         intree->LoadTree(iev);
//         bool isEB = formIsEB.EvalInstance();
//         bool isEE = formIsEE.EvalInstance();
//         if(!isEB && !isEE)
//             std::cout << "ERROR: GBRApply: not isEB nor isEE\n";
//         else if(isEB && isEE)
//             std::cout << "ERROR: GBRApply: isEB and isEE\n";
        
//         if(isEB)
//         {
//             for (int i=0; i<nvarsEB; ++i) {
//                 valsEB[i] = inputformsEB[i]->EvalInstance();
//             }
//             target = forestEB->GetResponse(valsEB);
//         }
//         else if(isEE)
//         {
//             for (int i=0; i<nvarsEE; ++i) {
//                 valsEE[i] = inputformsEE[i]->EvalInstance();
//             }
//             target = forestEE->GetResponse(valsEE);
//         }
//         target = offset + scale*sin(target);


//         friendtree->Fill();

//     }


//     //clear TTreeFormulas
//     for (std::vector<TTreeFormula*>::const_iterator it = inputformsEB.begin(); 
//             it != inputformsEB.end(); ++it) {
//         delete *it;
//     }
//     for (std::vector<TTreeFormula*>::const_iterator it = inputformsEE.begin(); 
//             it != inputformsEE.end(); ++it) {
//         delete *it;
//     }

//     delete[] valsEB;
//     delete[] valsEE;

//     intree->AddFriend(friendtree);

//     //the branch addresses are set to local variables in this function
//     //these local variables go out of scope after this function finishes
//     //so we need to reset the branch addresses before returning
//     friendtree->ResetBranchAddresses();
//     return friendtree;

// }



//_______________________________________________________________________
TTree *GBRApply::ApplyAsFriendTransform(TTree *intree, 
        const GBRForestD *forestEB, const GBRForestD *forestEE,
        const std::vector<std::string> &varsEB, const std::vector<std::string> &varsEE,
        const std::string& cutEB, const std::string& cutEE,
        std::string targetname,
        double low, double high,
        const bool doEB) const
{
    // transform parameters
    double scale = 0.5*(high-low);
    double offset = low + 0.5*(high-low);


    Float_t target = 0.;

    //initialize new friend tree
    TTree *friendtree = new TTree;
    friendtree->Branch(targetname.c_str(),&target,TString::Format("%s/F",targetname.c_str()));


    // ######################################
    // Process EB events
    // ######################################

    if(doEB){
        int nvarsEB = varsEB.size();

        // in EB or EE?
        std::cout << "Making TTreeFormula for EB" << std::endl;
        std::cout << "    cutEB = " << cutEB.c_str() << std::endl;
        TTreeFormula formIsEB(cutEB.c_str(), cutEB.c_str(), intree);

        //initialize TTreeFormulas to read variables from TTree
        std::vector<TTreeFormula*> inputformsEB;

        for (std::vector<std::string>::const_iterator it = varsEB.begin(); 
                it != varsEB.end(); ++it) {
            inputformsEB.push_back(new TTreeFormula(it->c_str(),it->c_str(),intree));
        }

        Float_t *valsEB = new Float_t[nvarsEB];

        for (Long64_t iev=0; iev<intree->GetEntries(); ++iev) {

            target = 0.;
            if (iev%100000==0) printf("%i\n",int(iev));

            intree->LoadTree(iev);

            bool isEB = formIsEB.EvalInstance();
            // bool isEE = formIsEE.EvalInstance();

            // if(!isEB && !isEE)
            //     std::cout << "ERROR: GBRApply: not isEB nor isEE\n";
            // else if(isEB && isEE)
            //     std::cout << "ERROR: GBRApply: isEB and isEE\n";
            
            if(isEB)
            {
                for (int i=0; i<nvarsEB; ++i) {
                    valsEB[i] = inputformsEB[i]->EvalInstance();
                }
                target = forestEB->GetResponse(valsEB);
            }
            // else if(isEE)
            // {
            //     for (int i=0; i<nvarsEE; ++i) {
            //         valsEE[i] = inputformsEE[i]->EvalInstance();
            //     }
            //     target = forestEE->GetResponse(valsEE);
            // }
            else {
                continue;
            }

            target = offset + scale*sin(target);

            friendtree->Fill();
        }


        //clear TTreeFormulas
        for (std::vector<TTreeFormula*>::const_iterator it = inputformsEB.begin(); 
                it != inputformsEB.end(); ++it) {
            delete *it;
        }

        delete[] valsEB;

    }

    // ######################################
    // Process EE events
    // ######################################

    if(!doEB){

        int nvarsEE = varsEE.size();

        std::cout << "Making TTreeFormula for EE" << std::endl;
        std::cout << "    cutEE = " << cutEE.c_str() << std::endl;
        TTreeFormula formIsEE(cutEE.c_str(), cutEE.c_str(), intree);

        std::vector<TTreeFormula*> inputformsEE;
        for (std::vector<std::string>::const_iterator it = varsEE.begin(); 
                it != varsEE.end(); ++it) {
            inputformsEE.push_back(new TTreeFormula(it->c_str(),it->c_str(),intree));
        }

        Float_t *valsEE = new Float_t[nvarsEE];

        for (Long64_t iev=0; iev<intree->GetEntries(); ++iev) {

            target = 0.;
            if (iev%100000==0) printf("%i\n",int(iev));

            intree->LoadTree(iev);

            // bool isEB = formIsEB.EvalInstance();
            bool isEE = formIsEE.EvalInstance();

            // if(!isEB && !isEE)
            //     std::cout << "ERROR: GBRApply: not isEB nor isEE\n";
            // else if(isEB && isEE)
            //     std::cout << "ERROR: GBRApply: isEB and isEE\n";
            
            // if(isEB)
            // {
            //     for (int i=0; i<nvarsEB; ++i) {
            //         valsEB[i] = inputformsEB[i]->EvalInstance();
            //     }
            //     target = forestEB->GetResponse(valsEB);
            // }
            if(isEE)
            {
                for (int i=0; i<nvarsEE; ++i) {
                    valsEE[i] = inputformsEE[i]->EvalInstance();
                }
                target = forestEE->GetResponse(valsEE);
            }
            else {
                continue;
            }

            target = offset + scale*sin(target);

            friendtree->Fill();
        }

        for (std::vector<TTreeFormula*>::const_iterator it = inputformsEE.begin(); 
                it != inputformsEE.end(); ++it) {
            delete *it;
        }

        delete[] valsEE;

    }

    // ######################################
    // Wrap up
    // ######################################

    intree->AddFriend(friendtree);

    //the branch addresses are set to local variables in this function
    //these local variables go out of scope after this function finishes
    //so we need to reset the branch addresses before returning
    friendtree->ResetBranchAddresses();
    return friendtree;

}
