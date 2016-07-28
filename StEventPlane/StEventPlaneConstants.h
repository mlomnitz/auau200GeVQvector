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

  //SL16d prduction
  TString qVectorRunDir ="/global/homes/m/mlomnitz/mlomnitz_projectdir/Run14_Reprod/recenter_v2_2/qVectorRun_v2_v3_2016-07-28";
   TString qVectorDayDir = "/global/homes/m/mlomnitz/mlomnitz_projectdir/Run14_Reprod/recenter_v2_2/qVectorDay_v2_v3_2016-07-28";
   //Event Cuts
   float const vzMax = 6.0;
   float const deltaVzMax = 3.0;

   //Track Cuts
   int const nHitsFitMin = 15;
   float const mNHitsFitRatioMin = 0.52;
   //Track cuts for event plane
   float const etaMaxEventPlane = 1.0;
   float const ptMinEventPlane  = 0.2;
   float const ptMaxEventPlane  = 2.0;
   float const dcaMaxEventPlane = 1.0;
}
#endif
