/**
 *  @file  ParReader.cpp
 *  @brief  
 *
 *
 *  @author  Jean-Baptiste Sauvan <sauvan@llr.in2p3.fr>
 *
 *  @date    11/11/2012
 *
 *  @internal
 *     Created :  11/11/2012
 * Last update :  11/11/2012 05:33:22 PM
 *          by :  JB Sauvan
 *
 * =====================================================================================
 */



#include "RegresTrainer/ParReader.h"
#include "RegresTrainer/Utilities.h"

#include <TEnv.h>

#include <iostream>
#include <sstream>

using namespace std;


/*****************************************************************/
ParReader::ParReader():
    m_trainer("TMVA"),
    m_outputDirectory("./"),
    m_factoryOptions("")
/*****************************************************************/
{
}


/*****************************************************************/
ParReader::~ParReader()
/*****************************************************************/
{
}


/*****************************************************************/
bool ParReader::read(const string& parFileName)
/*****************************************************************/
{
    TEnv params;
    int status = params.ReadFile(parFileName.c_str(),EEnvLevel(0));
    if(status!=0) 
    {
        cout << "FATAL: ParReader::read(): Cannot read file " << parFileName << "\n"; 
        return false;
    }

    m_trainer = params.GetValue("Trainer", "TMVA");
    m_factoryOptions = params.GetValue("TMVAFactoryOptions", "!V:!Silent:Color:DrawProgressBar");
    m_outputDirectory = params.GetValue("OutputDirectory", "./");
    int nRegressions = params.GetValue("NumberOfRegressions", 0);

    if(nRegressions==0)
        cout << "WARNING: ParReader::read(): 0 regressions requested\n";

    string baseName("Regression");
    for(int i=1;i<=nRegressions;i++)
    {

       stringstream keyInputFiles;
       stringstream keyTree;
       stringstream keyName;
       stringstream keyMethod;
       stringstream keyTmvaTrainingOptions;
       stringstream keyOptions;
       stringstream keyVariablesEB;
       stringstream keyVariablesEE;
       stringstream keyVariablesComb;
       stringstream keyTarget;
       stringstream keyTargetError;
       stringstream keyTargetComb;
       stringstream keyCutBase;
       stringstream keyCutError;
       stringstream keyCutComb;
       stringstream keyCutEB;
       stringstream keyCutEE;
       stringstream keyDoErrors;
       stringstream keyDoCombine;

       stringstream keyDoEB;
       stringstream meanMin;
       stringstream meanMax;
       stringstream fixMean;

       keyInputFiles           <<  baseName  <<  "."  <<  i  <<  ".InputFiles";
       keyTree                 <<  baseName  <<  "."  <<  i  <<  ".Tree";
       keyName                 <<  baseName  <<  "."  <<  i  <<  ".Name";
       keyMethod               <<  baseName  <<  "."  <<  i  <<  ".Method";
       keyTmvaTrainingOptions  <<  baseName  <<  "."  <<  i  <<  ".TMVATrainingOptions";
       keyOptions              <<  baseName  <<  "."  <<  i  <<  ".Options";
       keyVariablesEB          <<  baseName  <<  "."  <<  i  <<  ".VariablesEB";
       keyVariablesEE          <<  baseName  <<  "."  <<  i  <<  ".VariablesEE";
       keyVariablesComb        <<  baseName  <<  "."  <<  i  <<  ".VariablesComb";
       keyTarget               <<  baseName  <<  "."  <<  i  <<  ".Target";
       keyTargetError          <<  baseName  <<  "."  <<  i  <<  ".TargetError";
       keyTargetComb           <<  baseName  <<  "."  <<  i  <<  ".TargetComb";
       keyCutBase              <<  baseName  <<  "."  <<  i  <<  ".CutBase";
       keyCutError             <<  baseName  <<  "."  <<  i  <<  ".CutError";
       keyCutComb              <<  baseName  <<  "."  <<  i  <<  ".CutComb";
       keyCutEB                <<  baseName  <<  "."  <<  i  <<  ".CutEB";
       keyCutEE                <<  baseName  <<  "."  <<  i  <<  ".CutEE";
       keyDoErrors             <<  baseName  <<  "."  <<  i  <<  ".DoErrors";
       keyDoCombine            <<  baseName  <<  "."  <<  i  <<  ".DoCombine";

       keyDoEB                 <<  baseName  <<  "."  <<  i  <<  ".DoEB";
       meanMin                 <<  baseName  <<  "."  <<  i  <<  ".MeanMin";
       meanMax                 <<  baseName  <<  "."  <<  i  <<  ".MeanMax";
       fixMean                 <<  baseName  <<  "."  <<  i  <<  ".FixMean";


       RegressionParameters par;
       par.name                = string(params.GetValue(keyName.str().c_str(),               "regression"));
       par.inputFileNames      = string(params.GetValue(keyInputFiles.str().c_str(),         "ntuple.root"));
       par.treeName            = string(params.GetValue(keyTree.str().c_str(),               "myTree"));
       par.variablesEB         = string(params.GetValue(keyVariablesEB.str().c_str(),        ""));
       par.variablesEE         = string(params.GetValue(keyVariablesEE.str().c_str(),        ""));
       par.variablesComb       = string(params.GetValue(keyVariablesComb.str().c_str(),      ""));
       par.target              = string(params.GetValue(keyTarget.str().c_str(),             ""));
       par.targetError         = string(params.GetValue(keyTargetError.str().c_str(),        ""));
       par.targetComb          = string(params.GetValue(keyTargetComb.str().c_str(),         ""));
       par.method              = string(params.GetValue(keyMethod.str().c_str(),             "BDT"));
       par.tmvaTrainingOptions = string(params.GetValue(keyTmvaTrainingOptions.str().c_str(),"SplitMode=random:!V"));
       par.options             = string(params.GetValue(keyOptions.str().c_str(),            "!H:!V:BoostType=Grad"));
       par.cutBase             = string(params.GetValue(keyCutBase.str().c_str(),            ""));      
       par.cutError            = string(params.GetValue(keyCutError.str().c_str(),             ""));
       par.cutComb             = string(params.GetValue(keyCutComb.str().c_str(),            ""));
       par.cutEB               = string(params.GetValue(keyCutEB.str().c_str(),              ""));
       par.cutEE               = string(params.GetValue(keyCutEE.str().c_str(),              ""));
       par.doErrors            =        params.GetValue(keyDoErrors.str().c_str(),           false);
       par.doCombine           =        params.GetValue(keyDoCombine.str().c_str(),          false);

       par.doEB                =        params.GetValue(keyDoEB.str().c_str(),               true);
       
       par.meanMin             =        params.GetValue(meanMin.str().c_str(),               0.2);
       par.meanMax             =        params.GetValue(meanMax.str().c_str(),               2.0);
       par.fixMean             =        params.GetValue(fixMean.str().c_str(),               false);
       
       m_regParams.push_back(par);

    }

    return true;
}
