/**
 *  @file  HybridGBRMaker.h
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



#ifndef HYBRIDGBRMAKER_H
#define HYBRIDGBRMAKER_H

// STD
#include <string>
#include <vector>



class TTree;
class TChain;
class TFile;
class RooHybridBDTAutoPdf;
class RooDoubleCBFast;
class GBRTrainer;
class GBRForestD;

class HybridGBRMaker
{
    public:
        HybridGBRMaker();
        ~HybridGBRMaker();


        std::string name(){return m_name;};

        bool init(const std::string& name,
                  const std::string& fileNames,
                  const std::string& treeName,
                  const std::string& outputDirectory,
                  bool doCombine,
                  bool doEB,
		  float meanMin,
		  float meanMax, 
                  bool fixedMean
                  );
        void addVariableEB(const std::string& name);
        void addVariableEE(const std::string& name);
        void addVariableComb(const std::string& name);
        void addTarget(const std::string& target, const std::string& targetComb);
        void prepareTraining(const std::string& cutBase, const std::string& cutComb, const std::string& cutEB, const std::string& cutEE, const std::string& options);
        void run(const std::string& cutBase, const std::string& cutComb, const std::string& cutEB, const std::string& cutEE, const std::string& options);
        void runEB(const std::string& cutBase, const std::string& cutEB, const std::string& options);
        void runEE(const std::string& cutBase, const std::string& cutEE, const std::string& options);
        void runComb(const std::string& cutComb, const std::string& options);
        void close();

    private:
        std::string m_name;
        std::vector<TFile*> m_filesIn;
        std::vector<TTree*> m_trees;
        TChain* m_tree;
        std::string m_fileOutName;

        std::string m_target;
        std::string m_targetComb;
        std::string m_cutEB;
        std::string m_cutEE;
        std::string m_weight;

        bool m_doCombine;
        bool m_doEB;

        GBRTrainer* m_trainerComb;
        GBRForestD* m_forestEBmean;
        GBRForestD* m_forestEEmean;
        GBRForestD* m_forestEBwidth;
        GBRForestD* m_forestEEwidth;
        std::vector<std::string> m_variablesEB;
        std::vector<std::string> m_variablesEE;
        std::vector<std::string> m_variablesComb;
        int m_ntrees;

        float m_meanMin;
        float m_meanMax;
        bool m_fixedMean;

};


#endif
