/**
 *  @file  ParReader.h
 *  @brief  
 *
 *
 *  @author  Jean-Baptiste Sauvan <sauvan@llr.in2p3.fr>
 *
 *  @date    11/11/2012
 *
 *  @internal
 *     Created :  11/11/2012
 * Last update :  11/11/2012 05:12:03 PM
 *          by :  JB Sauvan
 *
 * =====================================================================================
 */



#ifndef PARREADER_H
#define PARREADER_H

#include <vector>
#include <string>


struct RegressionParameters
{
    std::string name;
    std::string inputFileNames;
    std::string treeName;
    std::string variablesEB;
    std::string variablesEE;
    std::string variablesComb;
    std::string target;
    std::string targetError;
    std::string targetComb;
    std::string method;
    std::string tmvaTrainingOptions;
    std::string options;
    std::string cutBase;
    std::string cutError;
    std::string cutComb;
    std::string cutEB;
    std::string cutEE;
    bool        doErrors;
    bool        doCombine;

    bool        doEB;
    float       meanMin;
    float       meanMax;
    bool        fixMean;
};


class ParReader
{
    public:
        ParReader();
        ~ParReader();

        bool read(const std::string& parFileName);

        inline std::vector<RegressionParameters>::iterator regressionBegin()
        {
            return m_regParams.begin();
        }
        inline std::vector<RegressionParameters>::iterator regressionEnd()
        {
            return m_regParams.end();
        }

        inline std::string& trainer()
        {
            return m_trainer;
        }
        inline std::string& outputDirectory()
        {
            return m_outputDirectory;
        }
        inline std::string& factoryOptions()
        {
            return m_factoryOptions;
        }



        std::string m_trainer;
        std::vector<RegressionParameters> m_regParams;
        std::string m_outputDirectory;
        std::string m_factoryOptions;


    private:
        // std::string m_trainer;
        // std::vector<RegressionParameters> m_regParams;
        // std::string m_outputDirectory;
        // std::string m_factoryOptions;
        


};

#endif
