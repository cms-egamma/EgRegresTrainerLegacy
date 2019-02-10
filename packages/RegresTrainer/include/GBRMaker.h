/**
 *  @file  GBRMaker.h
 *  @brief  
 *
 *
 *  @author  Jean-Baptiste Sauvan <sauvan@llr.in2p3.fr>
 *
 *  @date    11/28/2012
 *
 *  @internal
 *     Created :  11/28/2012
 * Last update :  11/28/2012 09:54:38 PM
 *          by :  JB Sauvan
 *
 * =====================================================================================
 */



#ifndef GBRMAKER_H
#define GBRMAKER_H

// STD
#include <string>
#include <vector>



class TTree;
class TFile;
class GBRTrainer;

class GBRMaker
{
    public:
        GBRMaker();
        ~GBRMaker();


        std::string name(){return m_name;};

        bool init(const std::string& name,
                  const std::string& fileNames,
                  const std::string& treeName,
                  const std::string& outputDirectory,
                  bool doErrors,
                  bool doCombine
                  );
        void addVariableEB(const std::string& name);
        void addVariableEE(const std::string& name);
        void addVariableComb(const std::string& name);
        void addTarget(const std::string& target, const std::string& targetError, const std::string& targetComb);
        void prepareTraining(const std::string& cutBase, const std::string& cutVar, const std::string& cutComb, const std::string& cutEB, const std::string& cutEE, const std::string& options);
        void run();
        void close();

    private:
        std::string m_name;
        std::vector<TFile*> m_filesIn;
        std::vector<TTree*> m_trees;
        TFile* m_fileOut;

        std::string m_cutEB;
        std::string m_cutEE;

        bool m_doErrors;
        bool m_doCombine;

        GBRTrainer* m_trainerEB;
        GBRTrainer* m_trainerEBVar;
        GBRTrainer* m_trainerEE;
        GBRTrainer* m_trainerEEVar;
        GBRTrainer* m_trainerComb;
        std::vector<std::string> m_variablesEB;
        std::vector<std::string> m_variablesEE;
        std::vector<std::string> m_variablesComb;
        int m_ntrees;

};


#endif
