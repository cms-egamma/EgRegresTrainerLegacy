/**
 *  @file  TMVAMaker.cpp
 *  @brief  
 *
 *
 *  @author  Jean-Baptiste Sauvan <sauvan@llr.in2p3.fr>
 *
 *  @date    11/11/2012
 *
 *  @internal
 *     Created :  11/11/2012
 * Last update :  11/11/2012 11:19:05 AM
 *          by :  JB Sauvan
 *
 * =====================================================================================
 */

#include "RegresTrainer/TMVAMaker.h"
#include "RegresTrainer/Utilities.h"

#include <TFile.h>
#include <TChain.h>
#include <TCut.h>
#include <TKey.h>
#include <TMVA/Ranking.h>
#include <TMVA/MethodBDT.h>

#include <iostream>
#include <sstream>


using namespace std;
using namespace TMVA;



/*****************************************************************/
TMVAMaker::TMVAMaker():
    m_name("regression"),
    m_tree(NULL),
    m_factory(NULL),
    m_dataloader(NULL),
    //m_fileIn(NULL),
    m_fileOut(NULL)
/*****************************************************************/
{
}

/*****************************************************************/
TMVAMaker::~TMVAMaker()
/*****************************************************************/
{

}


/*****************************************************************/
void TMVAMaker::close()
/*****************************************************************/
{
    if(m_fileOut) 
    {
        m_fileOut->Close();
        delete m_fileOut;
        m_fileOut = NULL;
    }
    //if(m_fileIn) 
    //{
    //    m_fileIn->Close();
    //    delete m_fileIn;
    //    m_fileIn = NULL;
    //}
    if(m_tree) 
    {
        //m_tree->Delete();
        delete m_tree;
        m_tree = NULL;
    }
    if(m_factory) 
    {
        //m_factory->Delete();
        delete m_factory;
        m_factory = NULL;
    }
    if(m_dataloader) 
    {
        //m_dataloader->Delete();
        delete m_dataloader;
        m_dataloader = NULL;
    }
    


}


/*****************************************************************/
bool TMVAMaker::init(const string& name,
                           const string& fileNames,
                           const string& treeName,
                           const string& options,
                           const string& outputDirectory)
/*****************************************************************/
{
    cout << "INFO: init regression " << name << "\n";
    m_name = name;
    //vector<string> treeNameTokens;
    //tokenize(treeName, treeNameTokens, "/");
    m_tree = new TChain(treeName.c_str());
    vector<string> vectorFileNames;
    tokenize(fileNames, vectorFileNames, ":");
    vector<string>::iterator it = vectorFileNames.begin();
    vector<string>::iterator itE = vectorFileNames.end();
    for(;it!=itE;++it)
    {
        TFile* fileIn = TFile::Open(it->c_str());
        if(!fileIn || !fileIn->IsOpen())
        {
            cout << "FATAL: RegressionMake::init(): Cannot open input file " << *it << "\n";
            return false;
        }
        TObject* key = fileIn->Get(treeName.c_str());
        if(!key || strcmp(key->ClassName(), "TTree")!=0)
        {
            cout << "FATAL: RegressionMake::init(): Cannot find regression tree " << treeName << " in " << *it << "\n";
            return false;
        }
        fileIn->Close();
        m_tree->Add(it->c_str());
    }
    
    stringstream outFileName;
    outFileName  <<  outputDirectory  <<  "/"  <<  name  <<  "_results.root";
    m_fileOut = TFile::Open(outFileName.str().c_str(), "RECREATE");
    if(!m_fileOut || !m_fileOut->IsOpen())
    {
        cout << "FATAL: RegressionMake::init(): Cannot open output file " << outFileName.str() << "\n";
        return false;
    }
    m_fileOut->cd();

    m_factory = new TMVA::Factory(name.c_str(), m_fileOut,options.c_str());
    m_dataloader = new TMVA::DataLoader(name.c_str());
    m_dataloader->AddRegressionTree(m_tree);

    return true;
}



/*****************************************************************/
void TMVAMaker::addVariable(const string& name)
/*****************************************************************/
{
    if(!m_dataloader)
    {
        cout << "ERROR: TMVAMaker::addVariable(): Cannot add variable: factory doesn't exist\n";
        return;
    }
    m_dataloader->AddVariable(name.c_str());
}

/*****************************************************************/
void TMVAMaker::addSpectator(const string& name)
/*****************************************************************/
{
    if(!m_dataloader)
    {
        cout << "ERROR: TMVAMaker::addSpectator(): Cannot add spectator: factory doesn't exist\n";
        return;
    }
    m_dataloader->AddSpectator(name.c_str());
}

/*****************************************************************/
void TMVAMaker::addTarget(const string& name)
/*****************************************************************/
{
    if(!m_dataloader)
    {
        cout << "ERROR: TMVAMaker::addTarget(): Cannot add target: factory doesn't exist\n";
        return;
    }
    m_dataloader->AddTarget(name.c_str());
}


/*****************************************************************/
void TMVAMaker::prepareTrainingAndTest(const string& cut, const string& options)
/*****************************************************************/
{
    TCut tcut(cut.c_str());
    m_dataloader->PrepareTrainingAndTestTree(tcut, options.c_str() ); 
}


/*****************************************************************/
void TMVAMaker::bookMethod(const string& method, const string& options)
/*****************************************************************/
{
  m_factory->BookMethod(m_dataloader, method.c_str(), m_name.c_str(), options.c_str()); 
}


/*****************************************************************/
void TMVAMaker::run()
/*****************************************************************/
{
    cout << "INFO: run regression " << m_name << "\n";
    m_fileOut->cd();
    m_factory->TrainAllMethods();
    MethodBDT* method = dynamic_cast<MethodBDT*>(m_factory->GetMethod(m_name.c_str(), ""));
    const Ranking* ranking = method->CreateRanking();
    cout << "INFO: TMVAMaker::run(): Printing variable ranking\n";
    ranking->Print();
    m_factory->TestAllMethods();
    m_factory->EvaluateAllMethods();
}
