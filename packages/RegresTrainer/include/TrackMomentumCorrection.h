/**
 *  @file  TrackMomentumCorrection.h
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

#ifndef TRACKMOMENTUMCORRECTION_H
#define TRACKMOMENTUMCORRECTION_H

#include "ICorrection.h"
#include <vector>


class TrackMomentumCorrection: public ICorrection
{
   public:
      TrackMomentumCorrection(const std::string& parameters):
          ICorrection(parameters)
      {
          m_inputNames.push_back("el_isEB");
          m_inputNames.push_back("el_classification");
          m_inputNames.push_back("el_gsftrk_pAtVtx");
      };
      ~TrackMomentumCorrection(){};

      virtual float call(const std::vector<float>& inputs) const;

};

#endif
