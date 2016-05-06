/* **************************************************
 *  Authors:  Mustafa Mustafa (mmustafa@lbl.gov)
 *            Hao Qiu         (hqiu@lbl.gov)
 *            **Guannan Xie     (guannanxie@lbl.gov)
 *
 *            ** code maintainer
 *
 * **************************************************
 */


#ifndef StEventPlaneConstants_H
#define StEventPlaneConstants_H

#include "TString.h"

namespace EventPlaneConstants
{
  int const nTrig = 5;
  int const mTriggerId[nTrig] = {450050, 450060,
				 450005, 450015,
				 450025 };
   TString qVectorRunDir = "/global/homes/q/qiuh/myProject/D0v2/recenter3/qVectorRun";
   TString qVectorDayDir = "/global/homes/q/qiuh/myProject/D0v2/recenter3/qVectorDay";
   //Event Cuts
   float const vzMax = 6.0;
   float const deltaVzMax = 3.0;

   //Track Cuts
   int const nHitsFitMin = 15;

   //Track cuts for event plane
   float const etaMaxEventPlane = 1.0;
   float const ptMinEventPlane  = 0.15;
   float const ptMaxEventPlane  = 2.0;
   float const dcaMaxEventPlane = 3.0;
}
#endif
