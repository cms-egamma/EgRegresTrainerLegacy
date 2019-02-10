/**
 *  @file  ErrorCorrection.h
 *  @brief  
 *
 *
 *  @author  Jean-Baptiste Sauvan <sauvan@llr.in2p3.fr>
 *
 *  @date    03/15/2013
 *
 *  @internal
 *     Created :  03/15/2013
 * Last update :  03/15/2013 11:44:43 AM
 *          by :  JB Sauvan
 *
 * =====================================================================================
 */

#ifndef ERRORCORRECTION_H
#define ERRORCORRECTION_H

#include "ICorrection.h"
#include <TGraphErrors.h>
#include <vector>


class ErrorCorrection: public ICorrection
{
   public:
      ErrorCorrection(const std::string& parameters):
          ICorrection(parameters)
      {
          m_inputNames.push_back("el_isEB");
          m_inputNames.push_back("el_scl_eta");
          m_inputNames.push_back("el_scl_r9");
      };
      ~ErrorCorrection(){};

      virtual float call(const std::vector<float>& inputs) const;

};

#endif
