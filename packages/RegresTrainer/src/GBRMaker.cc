/**
 *  @file  GBRMaker.cpp
 *  @brief  
 *
 *
 *  @author  Jean-Baptiste Sauvan <sauvan@llr.in2p3.fr>
 *
 *  @date    11/28/2012
 *
 *  @internal
 *     Created :  11/28/2012
 * Last update :  11/28/2012 10:15:12 PM
 *          by :  JB Sauvan
 *
 * =====================================================================================
 */



#include "RegresTrainer/GBRMaker.h"
#include "RegresTrainer/Utilities.h"
#include "RegresTrainer/GBRTrainer.h"
#include "CondFormats/GBRForest/interface/GBRForest.h"
#include "RegresTrainer/GBRApply.h"
#include "RegresTrainer/VariableCorrectionApply.h"
#include "RegresTrainer/SmearingCorrection.h"
#include "RegresTrainer/ErrorCorrection.h"
#include "RegresTrainer/TrackMomentumCorrection.h"

#include <TFile.h>
#include <TChain.h>
#include <TCut.h>
#include <TKey.h>

#include <iostream>
#include <sstream>
#include <time.h>


using namespace std;



/*****************************************************************/
GBRMaker::GBRMaker():
    m_name("regression"),
    m_fileOut(NULL),
    m_doErrors(true),
    m_doCombine(true),
    m_trainerEB(NULL),
    m_trainerEBVar(NULL),
    m_trainerEE(NULL),
    m_trainerEEVar(NULL),
    m_trainerComb(NULL),
    m_ntrees(10000)
/*****************************************************************/
{
}

/*****************************************************************/
GBRMaker::~GBRMaker()
/*****************************************************************/
{

}


/*****************************************************************/
void GBRMaker::close()
/*****************************************************************/
{
    vector<TFile*>::iterator it = m_filesIn.begin();
    vector<TFile*>::iterator itE = m_filesIn.end();
    for(;it!=itE;++it)
    {
        if(*it)
        {
            (*it)->Close();
            delete *it;
         }
    }
    m_filesIn.clear();
    if(m_fileOut) 
    {
        delete m_fileOut;
        m_fileOut = NULL;
    }
    //for(unsigned int i=0;i<m_trees.size();i++)
    //{
    //    //TTree* tree = m_trees[i];
    //    //if(tree) 
    //    //{
    //    //    tree->Delete();
    //    //    tree = NULL;
    //    //}
    //}
    m_trees.clear();
    if(m_trainerEB) 
    {
        delete m_trainerEB;
        m_trainerEB = NULL;
    }
    if(m_trainerEBVar) 
    {
        delete m_trainerEBVar;
        m_trainerEBVar = NULL;
    }
    if(m_trainerEE) 
    {
        delete m_trainerEE;
        m_trainerEE = NULL;
    }
    if(m_trainerEEVar) 
    {
        delete m_trainerEEVar;
        m_trainerEEVar = NULL;
    }
    if(m_trainerComb) 
    {
        delete m_trainerComb;
        m_trainerComb = NULL;
    }

}


/*****************************************************************/
bool GBRMaker::init(const string& name,
                           const string& fileNames,
                           const string& treeName,
                           const string& outputDirectory,
                           bool doErrors,
                           bool doCombine)
/*****************************************************************/
{
    cout << "INFO: init GBR regression " << name << "\n";
    m_name = name;
    m_doErrors = doErrors;
    m_doCombine = doCombine;
    if(m_doCombine && !m_doErrors) 
    {
        cout << "ERROR: combination requested, but without error estimation\n";
        return false;
    }
    // create GBR trainer
    m_trainerEB = new GBRTrainer();
    m_trainerEBVar = new GBRTrainer();
    m_trainerEE = new GBRTrainer();
    m_trainerEEVar = new GBRTrainer();
    m_trainerComb = new GBRTrainer();

    // fill GBRtrainer with TTrees
    vector<string> vectorFileNames;
    tokenize(fileNames, vectorFileNames, ":");
    vector<string>::iterator it = vectorFileNames.begin();
    vector<string>::iterator itE = vectorFileNames.end();
    for(;it!=itE;++it)
    {
        TFile* fileIn = TFile::Open(it->c_str());
        if(!fileIn || !fileIn->IsOpen())
        {
            cout << "FATAL: GBRMaker::init(): Cannot open input file " << *it << "\n";
            return false;
        }
        m_filesIn.push_back(fileIn);
        TTree* tree = (TTree*)fileIn->Get(treeName.c_str());
        if(!tree)
        {
            cout << "FATAL: GBRMaker::init(): Cannot find regression tree " << treeName << " in " << *it << "\n";
            return false;
        }
        m_trainerEB->AddTree(tree);
        m_trainerEBVar->AddTree(tree);
        m_trainerEE->AddTree(tree);
        m_trainerEEVar->AddTree(tree);
        m_trainerComb->AddTree(tree);
        m_trees.push_back(tree);
    }
    
    // open output file that will contain the GBRForest
    stringstream outFileName;
    outFileName  <<  outputDirectory  <<  "/"  <<  name  <<  "_results.root";
    m_fileOut = TFile::Open(outFileName.str().c_str(), "RECREATE");
    if(!m_fileOut || !m_fileOut->IsOpen())
    {
        cout << "FATAL: GBRMaker::init(): Cannot open output file " << outFileName.str() << "\n";
        return false;
    }
    m_fileOut->cd();


    return true;
}



/*****************************************************************/
void GBRMaker::addVariableEB(const string& name)
/*****************************************************************/
{
    if(!m_trainerEB || !m_trainerEBVar)
    {
        cout << "ERROR: GBRMaker::addVariable(): Cannot add variable: trainer doesn't exist\n";
        return;
    }
    m_variablesEB.push_back(name);
    m_trainerEB->AddInputVar(name);
    m_trainerEBVar->AddInputVar(name);
}

/*****************************************************************/
void GBRMaker::addVariableEE(const string& name)
/*****************************************************************/
{
    if(!m_trainerEE || !m_trainerEEVar)
    {
        cout << "ERROR: GBRMaker::addVariable(): Cannot add variable: trainer doesn't exist\n";
        return;
    }
    m_variablesEE.push_back(name);
    m_trainerEE->AddInputVar(name);
    m_trainerEEVar->AddInputVar(name);
}

/*****************************************************************/
void GBRMaker::addVariableComb(const string& name)
/*****************************************************************/
{
    if(!m_trainerComb)
    {
        cout << "ERROR: GBRMaker::addVariable(): Cannot add variable: trainer doesn't exist\n";
        return;
    }
    m_variablesComb.push_back(name);
    m_trainerComb->AddInputVar(name);
}


/*****************************************************************/
void GBRMaker::addTarget(const string& target, const string& targetError, const string& targetComb)
/*****************************************************************/
{
    if(!m_trainerEB || !m_trainerEBVar ||
            !m_trainerEE || !m_trainerEEVar ||
            !m_trainerComb)
    {
        cout << "ERROR: GBRMaker::addTarget(): Cannot add target: trainer doesn't exist\n";
        return;
    }
    //stringstream targetVar;
    //targetVar  <<  "1.4826*abs(BDTresponse - "  <<  name  <<  ")";
    m_trainerEB->SetTargetVar(target);
    m_trainerEBVar->SetTargetVar(targetError);
    m_trainerEE->SetTargetVar(target);
    m_trainerEEVar->SetTargetVar(targetError);
    m_trainerComb->SetTargetVar(targetComb);
}


/*****************************************************************/
void GBRMaker::prepareTraining(const string& cutBase, const string& cutVar, const string& cutComb, const string& cutEB, const string& cutEE, const string& options)
/*****************************************************************/
{
    cout << "INFO: Prepare training for " << m_name << "\n";
    m_cutEB = cutEB;
    m_cutEE = cutEE;
    TCut cutCentral(cutBase.c_str());
    TCut cutVariation(cutVar.c_str());
    TCut cutCombination(cutComb.c_str());
    TCut cutBarrel(cutEB.c_str());
    TCut cutEndcap(cutEE.c_str());
    cout << "INFO: Cuts for EB central value training = '" << string(cutCentral && cutBarrel) << "'\n";
    cout << "INFO: Cuts for EB uncertainty training    = '" << string(cutVariation && cutBarrel) << "'\n";
    cout << "INFO: Cuts for EE central value training = '" << string(cutCentral && cutEndcap) << "'\n";
    cout << "INFO: Cuts for EE uncertainty training    = '" << string(cutVariation && cutEndcap) << "'\n";
    cout << "INFO: Cuts for combination training    = '" << string(cutCombination) << "'\n";
    // set cut for training events
    m_trainerEB->SetTrainingCut(string(cutCentral && cutBarrel)); 
    m_trainerEBVar->SetTrainingCut(string(cutVariation && cutBarrel));
    m_trainerEE->SetTrainingCut(string(cutCentral && cutEndcap)); 
    m_trainerEEVar->SetTrainingCut(string(cutVariation && cutEndcap));
    m_trainerComb->SetTrainingCut(string(cutCombination)); 
    // loop over options
    vector<string> optionTokens;
    tokenize(options, optionTokens, ":");
    vector<string>::iterator it = optionTokens.begin();
    vector<string>::iterator itE = optionTokens.end();
    for(;it!=itE;++it)
    {
        string token = *it;
        vector<string> tagAndValue;
        tokenize(token, tagAndValue, "=");
        if(tagAndValue.size()!=2)
        {
            cout << "ERROR: GBRMaker::prepareTraining(): option " << token << " cannot be processed. Should be of the form tag=value.\n";
            continue;
        }
        string tag = tagAndValue[0];
        string value = tagAndValue[1];
        if(tag=="MinEvents")
        {
            int minEvents = 0;
            fromString(minEvents, value);
            m_trainerEB->SetMinEvents(minEvents);
            m_trainerEBVar->SetMinEvents(minEvents);
            m_trainerEE->SetMinEvents(minEvents);
            m_trainerEEVar->SetMinEvents(minEvents);
            m_trainerComb->SetMinEvents(minEvents);
            cout << "INFO: MinEvents = " << minEvents << "\n";
        }
        else if(tag=="Shrinkage")
        {
            float shrink = 0.;
            fromString(shrink, value);
            m_trainerEB->SetShrinkage(shrink);
            m_trainerEBVar->SetShrinkage(shrink);
            m_trainerEE->SetShrinkage(shrink);
            m_trainerEEVar->SetShrinkage(shrink);
            m_trainerComb->SetShrinkage(shrink);
            cout << "INFO: Shrinkage = " << shrink << "\n";
        }
        else if(tag=="MinSignificance")
        {
            float sig = 0.;
            fromString(sig, value);
            m_trainerEB->SetMinCutSignificance(sig);
            m_trainerEBVar->SetMinCutSignificance(sig);
            m_trainerEE->SetMinCutSignificance(sig);
            m_trainerEEVar->SetMinCutSignificance(sig);
            m_trainerComb->SetMinCutSignificance(sig);
            cout << "INFO: MinSignificance = " << sig << "\n";
        }
        else if(tag=="TransitionQuantile")
        {
            float trans = 0.;
            fromString(trans, value);
            m_trainerEB->SetTransitionQuantile(trans);
            m_trainerEBVar->SetTransitionQuantile(trans);
            m_trainerEE->SetTransitionQuantile(trans);
            m_trainerEEVar->SetTransitionQuantile(trans);
            m_trainerComb->SetTransitionQuantile(trans);
            cout << "INFO: TransitionQuantile = " << trans << "\n";
        }
        else if(tag=="NTrees")
        {
            int ntrees= 0;
            fromString(ntrees, value);
            m_ntrees = ntrees;
            cout << "INFO: NTrees = " << ntrees << "\n";
        }
        else if(tag=="RandomSeed")
        {
            m_trainerEB->SetRandomSeed(value);
            m_trainerEBVar->SetRandomSeed(value);
            m_trainerEE->SetRandomSeed(value);
            m_trainerEEVar->SetRandomSeed(value);
            m_trainerComb->SetRandomSeed(value);
            cout << "INFO: RandomSeed = " << value << "\n";
        }
        else if(tag=="EventWeight")
        {
            m_trainerEB->SetEventWeight(value);
            m_trainerEBVar->SetEventWeight(value);
            m_trainerEE->SetEventWeight(value);
            m_trainerEEVar->SetEventWeight(value);
            m_trainerComb->SetEventWeight(value);
            cout << "INFO: EventWeight = " << value << "\n";
        }
        else
        {
            cout << "ERROR: GBRMaker::prepareTraining(): Unknown option " << tag << "\n";
            cout << "ERROR: Possibilities are: MinEvents, Shrinkage, MinSignificance, TransitionQuantile, NTrees\n";
        }
    }
}



/*****************************************************************/
void GBRMaker::run()
/*****************************************************************/
{
    cout << "INFO: run GBR regression " << m_name << "\n";
    // run training for EB and EE central values
    cout << "INFO: train BDT for EB central value estimation\n";
    const GBRForest* forestEB = m_trainerEB->TrainForest(m_ntrees);
    delete m_trainerEB;
    m_trainerEB = NULL;

    cout << "INFO: train BDT for EE central value estimation\n";
    const GBRForest* forestEE = m_trainerEE->TrainForest(m_ntrees);
    delete m_trainerEE;
    m_trainerEE = NULL;

    //const SmearingCorrection* energyCorrection = new SmearingCorrection("");
    const GBRForest* forestEBVar = NULL;
    const GBRForest* forestEEVar = NULL;
    if(m_doErrors)
    {
        // add BDT response to the tree
        cout << "INFO: filling BDT response in tree\n";
        GBRApply gbrApply;
        for(unsigned int t=0; t<m_trees.size(); t++)
            gbrApply.ApplyAsFriend(m_trees[t], forestEB, forestEE,
                    m_variablesEB, m_variablesEE,
                    m_cutEB, m_cutEE,
                    "BDTresponse");
        //cout << "INFO: filling corrected BDT response in tree\n";
        //for(unsigned int t=0; t<m_trees.size(); t++)
        //    gbrApply.ApplyAsFriend(m_trees[t], forestEB, forestEE,
        //            m_variablesEB, m_variablesEE,
        //            m_cutEB, m_cutEE,
        //            "BDTresponseCorr",
        //            energyCorrection);

        // run training for EB and EE errors
        cout << "INFO: train BDT for EB uncertainty estimation\n";
        forestEBVar = m_trainerEBVar->TrainForest(m_ntrees);
        delete m_trainerEBVar;
        m_trainerEBVar = NULL;

        cout << "INFO: train BDT for EE uncertainty estimation\n";
        forestEEVar = m_trainerEEVar->TrainForest(m_ntrees);
        delete m_trainerEEVar;
        m_trainerEEVar = NULL;

        //const ErrorCorrection* errorCorrection = new ErrorCorrection("");

        // add BDT error response to the tree
        cout << "INFO: filling BDT error in tree\n";
        for(unsigned int t=0; t<m_trees.size(); t++)
            gbrApply.ApplyAsFriend(m_trees[t], forestEBVar, forestEEVar,
                    m_variablesEB, m_variablesEE,
                    m_cutEB, m_cutEE,
                    "BDTerror");//,
        //errorCorrection);

    }
    // add corrected momentum to the tree
    //VariableCorrectionApply varCorrApply;
    //const TrackMomentumCorrection* trackMomentumCorrection = new TrackMomentumCorrection("");
    //cout << "INFO: filling corrected track momentum in tree\n";
    //for(unsigned int t=0; t<m_trees.size(); t++)
    //    varCorrApply.ApplyAsFriend(m_trees[t], "el_gsftrk_pAtVtx",
    //            "el_gsftrk_pAtVtxCorr",
    //            trackMomentumCorrection);

    // run training for combination
    const GBRForest* forestComb = NULL;
    if(m_doCombine)
    {
        cout << "INFO: train BDT for combination\n";
        forestComb = m_trainerComb->TrainForest(m_ntrees);
        delete m_trainerComb;
        m_trainerComb = NULL;
    }

    // write GBRForest in output file
    m_fileOut->cd();
    if(forestEB)
        m_fileOut->WriteObject(forestEB, "EBCorrection");
    else
        cout << "WARNING: GBRMaker::run(): NULL EB forest\n";
    if(m_doErrors)
    {
        if(forestEBVar)
            m_fileOut->WriteObject(forestEBVar, "EBUncertainty");
        else
            cout << "WARNING: GBRMaker::run(): NULL forestEBVar\n";
    }
    if(forestEE)
        m_fileOut->WriteObject(forestEE, "EECorrection");
    else
        cout << "WARNING: GBRMaker::run(): NULL EE forest\n";
    if(m_doErrors)
    {
        if(forestEEVar)
            m_fileOut->WriteObject(forestEEVar, "EEUncertainty");
        else
            cout << "WARNING: GBRMaker::run(): NULL forestEEVar\n";
    }
    if(m_doCombine)
    {
        if(forestComb)
            m_fileOut->WriteObject(forestComb, "CombinationWeight");
        else
            cout << "WARNING: GBRMaker::run(): NULL comb forest\n";
    }

    m_fileOut->WriteObject(&m_variablesEB, "varlistEB");
    m_fileOut->WriteObject(&m_variablesEE, "varlistEE");
    if(m_doCombine)
        m_fileOut->WriteObject(&m_variablesComb, "varlistComb");
    m_fileOut->Close();
}



