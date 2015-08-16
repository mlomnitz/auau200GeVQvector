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
class TProfile;

const int maxNTracks = 20000;

class StEventPlane : public StMaker {
  public:
  StEventPlane(const char* name, StPicoDstMaker* picoMaker, StRefMultCorr* grefmultCorrUtil);
  
  Int_t Init();
  virtual Int_t Make();
    void setFileOut(TFile* fileOut);

  int   getCentrality() const;
  float getEventPlane(int nTracksToExclude=0, int* indexTracksToExclude=0) const;
  float getEventPlane1() const;
  float getEventPlane2() const;
  float getEventPlaneEtaPlus() const;
  float getEventPlaneEtaMinus() const;
  float getResolutionRandom() const;
  float getResolutionEta() const;
  void calculateHadronV2() const;
  int eventPlaneStatus() const;
  
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

  TH3F*      hHadronV2PtCent;
  TH3F*      hHadronHftV2PtCent;
  TH3F*      hHadronPrimaryV2PtCent;
  TH3F*      hHadronHftPrimaryV2PtCent;

  bool	mAcceptEvent;
  bool	mAcceptQvectorFile;
  bool	mAcceptQvectorFiletmp;

  int         mCent;
  int         mRunNumber;
  float       mBField;
  StThreeVectorF mVertexPos;
  int mEventPlaneStatus;
  float       mEventPlane, mEventPlane1, mEventPlane2, mEventPlaneEtaPlus, mEventPlaneEtaMinus;
  float       mResolutionRandom, mResolutionEta;
  TVector2    mQ, mQ1, mQ2, mQEtaPlus, mQEtaMinus;
  TProfile*  prfQxCentEtaPlus;
  TProfile*  prfQyCentEtaPlus;
  TProfile*  prfQxCentEtaMinus;
  TProfile*  prfQyCentEtaMinus;    
  
  float      qxTracks[maxNTracks];
  float      qyTracks[maxNTracks];
  
  ClassDef(StEventPlane, 0)
};
inline int   StEventPlane::getCentrality() const { return mCent; }
inline float StEventPlane::getEventPlane1() const { return mEventPlane1; }
inline float StEventPlane::getEventPlane2() const { return mEventPlane2; }
inline float StEventPlane::getEventPlaneEtaPlus() const  { return mEventPlaneEtaPlus; }
inline float StEventPlane::getEventPlaneEtaMinus() const { return mEventPlaneEtaMinus; }
inline float StEventPlane::getResolutionRandom() const { return mResolutionRandom; }
inline float StEventPlane::getResolutionEta() const { return mResolutionEta; }
inline int   StEventPlane::getRunId() const { return mRunNumber; }
inline bool  StEventPlane::getAcceptEvent() const { return mAcceptQvectorFile && mAcceptQvectorFiletmp; }
inline int StEventPlane::eventPlaneStatus() const { return mEventPlaneStatus; }
#endif
