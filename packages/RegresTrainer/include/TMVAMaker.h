/**
 *  @file  TMVAMaker.h
 *  @brief  
 *
 *
 *  @author  Jean-Baptiste Sauvan <sauvan@llr.in2p3.fr>
 *
 *  @date    11/11/2012
 *
 *  @internal
 *     Created :  11/11/2012
 * Last update :  11/11/2012 10:52:53 AM
 *          by :  JB Sauvan
 *
 * =====================================================================================
 */



#ifndef TMVAMAKER_H
#define TMVAMAKER_H

// STD
#include <string>

// ROOT
#include <TMVA/Factory.h>
#include <TMVA/DataLoader.h>

class TChain;

class TMVAMaker
{
    public:
        TMVAMaker();
        ~TMVAMaker();

        std::string name(){return m_name;};

        bool init(const std::string& name,
                  const std::string& fileNames,
                  const std::string& treeName,
                  const std::string& options,
                  const std::string& outputDirectory
                  );
        void addVariable(const std::string& name);
        void addSpectator(const std::string& name);
        void addTarget(const std::string& name);
        void prepareTrainingAndTest(const std::string& cut, const std::string& options);
        void bookMethod(const std::string& method, const std::string& options);
        void run();
        void close();

    private:
        std::string m_name;
        TChain* m_tree;
        TMVA::Factory* m_factory;
	TMVA::DataLoader *m_dataloader;
        //TFile* m_fileIn;
        TFile* m_fileOut;

};


#endif
