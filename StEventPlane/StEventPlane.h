/* **************************************************
 *  Authors:  Mustafa Mustafa (mmustafa@lbl.gov)
 *            Hao Qiu         (hqiu@lbl.gov)
 *            **Guannan Xie     (guannanxie@lbl.gov)
 *
 *            ** code maintainer
 *
 * **************************************************
 */


#ifndef STAR_StEventPlane
#define STAR_StEventPlane

#include "StMaker.h"
#include "TVector2.h"
#include "StThreeVectorF.hh"

class StPicoDst;
class StPicoEvent;
class StPicoDstMaker;
class StRefMultCorr;
class TH1I;
class TH1F;
class TH2F;
class TH3F;
class THn;
class TProfile;

const int maxNTracks = 20000;

class StEventPlane : public StMaker
{
public:
  StEventPlane(const char* name, StPicoDstMaker* picoMaker, StRefMultCorr* grefmultCorrUtil, int harmonic = 2);

   Int_t Init();
   virtual Int_t Make();
   Int_t Finish();
   void setFileOut(TFile* fileOut);

   int   getCentrality() const;
   float getEventPlane(int nTracksToExclude = 0, int* indexTracksToExclude = 0) const;
   float getEventPlane1() const;
   float getEventPlane2() const;
   float getEventPlaneEtaPlus() const;
   float getEventPlaneEtaMinus() const;
   float getResolutionRandom() const;
   float getResolutionEta() const;
   void calculateHadronVn() const;
   int eventPlaneStatus() const;
   int harmonic() const;
   TVector2 Q() const;
   TVector2 QEtaPlusGap005() const;
   TVector2 QEtaMinusGap005() const;
   TVector2 QEtaGap(int iEta, int nEtaGaps) const;
   TVector2 QEta(int iEta) const;
   TVector2 q(int iTrack) const;

   int   getRunId() const;
   bool getAcceptEvent() const;

private:
   void getEventInfo();
   void getRunInfo(int runNumber);
   int calculateEventPlane();

   StPicoDstMaker* mPicoDstMaker;
   StPicoDst*      mPicoDst;
   StPicoEvent*    mPicoEvent;
   StRefMultCorr* mgrefmultCorrUtil;

   TFile* mFileOut;

   //track level qa
   TH1I*      hNHitsFit;
   TH1F*      hDca;
   TH1F*      hEta;
   TH1F*      hPt;

   TH2F*      hOneOverBetaDiffKaonP;
   TH2F*      hOneOverBetaDiffPionP;

   //event plane and Q vector
   TH2F*      hPhiCentEtaPlusZPlus;
   TH2F*      hPhiCentEtaPlusZMinus;
   TH2F*      hPhiCentEtaMinusZPlus;
   TH2F*      hPhiCentEtaMinusZMinus;

   TH2F*      hEventPlaneCent;
   TH2F*      hEventPlane1Cent;
   TH2F*      hEventPlane2Cent;
   TH2F*      hEventPlaneEtaPlusCent;
   TH2F*      hEventPlaneEtaMinusCent;
   TH3F*      hQyQxCent;
   TH3F*      hQyQx1Cent;
   TH3F*      hQyQx2Cent;
   TH3F*      hQyQxEtaPlusCent;
   TH3F*      hQyQxEtaMinusCent;
   TProfile*  prfCosResolutionRandomCent;
   TProfile*  prfCosResolutionEtaCent;

   TH3F*      hHadronVnPtCent;
   TH3F*      hHadronHftVnPtCent;
   TH3F*      hHadronPrimaryVnPtCent;
   TH3F*      hHadronHftPrimaryVnPtCent;
   THn*       hHadronVnPtCentEtaGap;

   bool   mAcceptEvent;
   bool   mAcceptQvectorFile;
   bool   mAcceptQvectorFiletmp;

   int         mHarmonic;
   int         mCent;
   int         mRunNumber;
   float       mBField;
   StThreeVectorF mVertexPos;
   int mEventPlaneStatus;
   float       mEventPlane, mEventPlane1, mEventPlane2, mEventPlaneEtaPlus, mEventPlaneEtaMinus;
   float       mResolutionRandom, mResolutionEta;
   TVector2    mQ, mQ1, mQ2, mQEtaPlus, mQEtaMinus, mQEtaPlusGap005, mQEtaMinusGap005;
   TVector2 mQEta[20];
   TProfile*  prfQxCentEtaPlus;
   TProfile*  prfQyCentEtaPlus;
   TProfile*  prfQxCentEtaMinus;
   TProfile*  prfQyCentEtaMinus;

   float      qxTracks[maxNTracks];
   float      qyTracks[maxNTracks];

   ClassDef(StEventPlane, 0)
};
inline int   StEventPlane::getCentrality() const
{
   return mCent;
}
inline float StEventPlane::getEventPlane1() const
{
   return mEventPlane1;
}
inline float StEventPlane::getEventPlane2() const
{
   return mEventPlane2;
}
inline float StEventPlane::getEventPlaneEtaPlus() const
{
   return mEventPlaneEtaPlus;
}
inline float StEventPlane::getEventPlaneEtaMinus() const
{
   return mEventPlaneEtaMinus;
}
inline float StEventPlane::getResolutionRandom() const
{
   return mResolutionRandom;
}
inline float StEventPlane::getResolutionEta() const
{
   return mResolutionEta;
}
inline int   StEventPlane::getRunId() const
{
   return mRunNumber;
}
inline bool  StEventPlane::getAcceptEvent() const
{
   return mAcceptQvectorFile && mAcceptQvectorFiletmp;
}
inline int StEventPlane::eventPlaneStatus() const
{
   return mEventPlaneStatus;
}
inline int StEventPlane::harmonic() const
{
  return mHarmonic;
}
inline TVector2 StEventPlane::Q() const
{
   return mQ;
}
inline TVector2 StEventPlane::QEtaPlusGap005() const
{
  return mQEtaPlusGap005;
}
inline TVector2 StEventPlane::QEtaMinusGap005() const
{
  return mQEtaMinusGap005;
}
inline TVector2 StEventPlane::QEta(int iEta) const
{
  return mQEta[iEta];
}
inline TVector2 StEventPlane::q(int iTrack) const
{
   TVector2 q_(qxTracks[iTrack], qyTracks[iTrack]);
   return q_;
}
#endif
