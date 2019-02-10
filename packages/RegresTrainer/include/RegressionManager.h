/**
 *  @file  RegressionManager.h
 *  @brief  
 *
 *
 *  @author  Jean-Baptiste Sauvan <sauvan@lal.in2p3.fr>
 *
 *  @date    11/12/2012
 *
 *  @internal
 *     Created :  11/12/2012
 * Last update :  11/12/2012 09:37:57 AM
 *          by :  JB Sauvan
 *
 * =====================================================================================
 */

#ifndef REGRESSIONMANAGER_H
#define REGRESSIONMANAGER_H

#include "ParReader.h"

#include <string>
#include <vector>

class RegressionManager
{
    public:
        RegressionManager();
        ~RegressionManager();

        bool init(const std::string& parFileName);
        bool makeRegression();

    private:
        ParReader m_reader;

        std::vector<std::string> m_variables;
        // REMOVE
        // std::vector<std::string> m_variablesDouble;
        // std::vector<std::string> m_variablesInt;
        // std::vector<std::string> m_variablesBool;

};

#endif
