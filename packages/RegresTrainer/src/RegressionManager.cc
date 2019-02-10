/**
 *  @file  RegressionManager.cpp
 *  @brief  
 *
 *
 *  @author  Jean-Baptiste Sauvan <sauvan@lal.in2p3.fr>
 *
 *  @date    11/12/2012
 *
 *  @internal
 *     Created :  11/12/2012
 * Last update :  11/30/2012 03:40:03 PM
 *          by :  JB Sauvan
 *
 * =====================================================================================
 */


#include "RegresTrainer/Utilities.h"
#include "RegresTrainer/RegressionManager.h"
#include "RegresTrainer/TMVAMaker.h"
#include "RegresTrainer/GBRMaker.h"
#include "RegresTrainer/GBRForest.h"
#include "RegresTrainer/HybridGBRMaker.h"


#include <TSystem.h>
#include <TFile.h>

#include <iostream>
#include <sstream>
#include <algorithm>


using namespace std;

/*****************************************************************/
RegressionManager::RegressionManager()
/*****************************************************************/
{
}



/*****************************************************************/
RegressionManager::~RegressionManager()
/*****************************************************************/
{
}



/*****************************************************************/
bool RegressionManager::init(const string& parFileName)
/*****************************************************************/
{
    bool status = true;

    // load dictionary in order to be able to read/write GBRForest objects from/to ROOT file
    //gSystem->Load("obj/libDictionary_C.so");
    // retrieve parameters from parameter file
    status = m_reader.read(parFileName);

    return status;
}



/*****************************************************************/
bool RegressionManager::makeRegression()
/*****************************************************************/
{
    bool status = true;
    // retrieve type of trainer: TMVA or GBRTrain
    string trainerType = m_reader.trainer();
    // loop over all regressions
    vector<RegressionParameters>::iterator it = m_reader.regressionBegin();
    vector<RegressionParameters>::iterator itE = m_reader.regressionEnd();
    for(;it!=itE;++it)
    {
        cout << "INFO: Creating new " << trainerType << " regression\n";

        // measure start time
        time_t tStart, tEnd;
        time(&tStart);
        if(trainerType=="TMVA")
        {
            // initialize TMVA regressionMaker
            TMVAMaker* regMaker = new TMVAMaker();
            status = regMaker->init(it->name,
                    it->inputFileNames,
                    it->treeName,
                    m_reader.factoryOptions(),
                    m_reader.outputDirectory()
                    );
            if(!status)
                break;
            // feed regressionMaker with input variable names
            string variables = it->variablesEB;
            vector<string> tokens;
            tokenize(variables, tokens, ":");
            vector<string>::iterator itv  = tokens.begin();
            vector<string>::iterator itvE = tokens.end();
            for(;itv!=itvE;++itv)
                regMaker->addVariable(*itv);
            // give target, book regression method, options, etc.
            regMaker->addTarget(it->target);
            regMaker->prepareTrainingAndTest(it->cutBase, it->tmvaTrainingOptions);
            regMaker->bookMethod(it->method, it->options);

            // run regression
            regMaker->run();
            regMaker->close();
            delete regMaker;
        }
        else if(trainerType=="GBRTrain")
        {
            // initialize GBR regression
            GBRMaker* regMaker = new GBRMaker();
            status = regMaker->init(it->name,
                    it->inputFileNames,
                    it->treeName,
                    m_reader.outputDirectory(),
                    it->doErrors,
                    it->doCombine
                    );
            if(!status)
                break;
            // feed with variables
            string variablesEB = it->variablesEB;
            string variablesEE = it->variablesEE;
            string variablesComb = it->variablesComb;
            vector<string> tokensEB;
            vector<string> tokensEE;
            vector<string> tokensComb;
            tokenize(variablesEB, tokensEB, ":");
            tokenize(variablesEE, tokensEE, ":");
            tokenize(variablesComb, tokensComb, ":");
            vector<string>::iterator itv  = tokensEB.begin();
            vector<string>::iterator itvE = tokensEB.end();
            for(;itv!=itvE;++itv)
                regMaker->addVariableEB(*itv);
            itv  = tokensEE.begin();
            itvE = tokensEE.end();
            for(;itv!=itvE;++itv)
                regMaker->addVariableEE(*itv);
            itv  = tokensComb.begin();
            itvE = tokensComb.end();
            for(;itv!=itvE;++itv)
                regMaker->addVariableComb(*itv);
            tokensEB.clear();
            tokensEE.clear();
            tokensComb.clear();
            // specify target, cuts, fill options
            regMaker->addTarget(it->target, it->targetError, it->targetComb);
            regMaker->prepareTraining(it->cutBase, it->cutError, it->cutComb, it->cutEB, it->cutEE, it->options);

            // run regression
            regMaker->run();
            regMaker->close();
            delete regMaker;
        }
        else if(trainerType=="GBRLikelihoodTrain")
        {
            // initialize GBR regression
            HybridGBRMaker* regMaker = new HybridGBRMaker();
            status = regMaker->init(it->name,
                    it->inputFileNames,
                    it->treeName,
                    m_reader.outputDirectory(),
                    it->doCombine,
                    it->doEB
                    );
            if(!status)
                break;
            // feed with variables
            string variablesEB = it->variablesEB;
            string variablesEE = it->variablesEE;
            string variablesComb = it->variablesComb;
            vector<string> tokensEB;
            vector<string> tokensEE;
            vector<string> tokensComb;
            tokenize(variablesEB, tokensEB, ":");
            tokenize(variablesEE, tokensEE, ":");
            tokenize(variablesComb, tokensComb, ":");
            vector<string>::iterator itv  = tokensEB.begin();
            vector<string>::iterator itvE = tokensEB.end();
            for(;itv!=itvE;++itv)
                regMaker->addVariableEB(*itv);
            itv  = tokensEE.begin();
            itvE = tokensEE.end();
            for(;itv!=itvE;++itv)
                regMaker->addVariableEE(*itv);
            itv  = tokensComb.begin();
            itvE = tokensComb.end();
            for(;itv!=itvE;++itv)
                regMaker->addVariableComb(*itv);
            tokensEB.clear();
            tokensEE.clear();
            tokensComb.clear();
            // specify target, cuts, fill options
            regMaker->addTarget(it->target, it->targetComb);
            //regMaker->prepareTraining(it->cutBase, it->cutComb, it->cutEB, it->cutEE, it->options);

            // run regression
            regMaker->run(it->cutBase, it->cutComb, it->cutEB, it->cutEE, it->options);
            regMaker->close();
            delete regMaker;
        }
        else
        {
            cout << "WARNING: RegressionManager::makeRegression(): unknown trainer " << trainerType << "\n";
        }
        // measure end time
        time(&tEnd);
        double t = difftime(tEnd, tStart);
        int hours = int(t)/3600;
        int min = (int(t)/60)%60;

        cout << "INFO: RegressionManager::makeRegression(): Elapsed time = " << t << " s\n";
        cout << "                                                        = " << hours << " h " << min << " min\n"; 
    }
    return status;
}

