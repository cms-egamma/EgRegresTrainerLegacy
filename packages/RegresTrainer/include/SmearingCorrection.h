/**
 *  @file  SmearingCorrection.h
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

#ifndef SMEARINGCORRECTION_H
#define SMEARINGCORRECTION_H

#include "ICorrection.h"
#include <TGraphErrors.h>
#include <vector>
#include <TRandom3.h>


class SmearingCorrection: public ICorrection
{
   public:
      SmearingCorrection(const std::string& parameters):
          ICorrection(parameters)
      {
          m_inputNames.push_back("el_isEB");
          m_inputNames.push_back("el_scl_eta");
          m_inputNames.push_back("el_scl_r9");
          m_random = new TRandom3();
          m_random->SetSeed(1234);
      };
      ~SmearingCorrection()
      {
          delete m_random;
      };

      virtual float call(const std::vector<float>& inputs) const;

    private:
      TRandom3* m_random;

};

#endif
