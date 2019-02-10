/**
 *  @file  main.cpp
 *  @brief  
 *
 *
 *  @author  Jean-Baptiste Sauvan <sauvan@llr.in2p3.fr>
 *
 *  @date    11/12/2012
 *
 *  @internal
 *     Created :  11/12/2012
 * Last update :  11/12/2012 10:56:02 AM
 *          by :  JB Sauvan
 *
 * =====================================================================================
 */




#include <iostream>
#include <string>

#include "RegresTrainer/RegressionManager.h"


int main(int argc, char** argv)
{
    if(argc!=2)
    {
        std::cout << "Usage: regression.exe configurationFile\n";
        return 1;
    }

    std::string parameterFile(argv[1]);
    RegressionManager manager;
    bool status = true;
    status = manager.init(parameterFile);
    if(status)
    {
        manager.makeRegression();
    }

    if(!status)
        std::cout << "FATAL: A fatal error occured - QUIT -\n";
    else
        std::cout << "- Finish - All good -\n";


    return status;
}
