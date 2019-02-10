/**
 *  @file  ICorrection.h
 *  @brief  
 *
 *
 *  @author  Jean-Baptiste Sauvan <sauvan@llr.in2p3.fr>
 *
 *  @date    03/15/2013
 *
 *  @internal
 *     Created :  03/15/2013
 * Last update :  03/15/2013 11:41:06 AM
 *          by :  JB Sauvan
 *
 * =====================================================================================
 */

#ifndef ICORRECTION_H
#define ICORRECTION_H

#include <vector>
#include <string>


class ICorrection
{
    public:
        ICorrection(const std::string& parameters){};
        float operator()(const std::vector<float>& inputs) const
        {
            return call(inputs);
        }
        virtual float call(const std::vector<float>& inputs) const =0;         
        const std::vector<std::string>& inputNames() const
        {
            return m_inputNames;
        }
    protected:
        std::vector<std::string> m_inputNames;
};


#endif
