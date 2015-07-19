#ifndef STAR_StEventPlane
#define STAR_StEventPlane

#include "TChain.h"
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
  virtual ~StEventPlane();
  
  virtual Int_t Init();
  virtual Int_t Make();
  virtual void  Clear(Option_t* opt="");
  virtual Int_t Finish();
  
  int getEntries() const;

  int   getCentrality() const;
  float getEventPlane() const;
  float getEventPlane1() const;
  float getEventPlane2() const;
  float getEventPlaneEtaPlus() const;
  float getEventPlaneEtaMinus() const;
  float getResolutionRandom() const;
  float getResolutionEta() const;

  int   getRunId() const;
  bool getAcceptEvent() const;

 private:
  void getEventInfo();
  int calculateEventPlane();

  StPicoDstMaker* mPicoDstMaker;
  StPicoDst*      mPicoDst;
  StPicoEvent*    mPicoEvent;
  StRefMultCorr* mgrefmultCorrUtil;
  
  bool	mAcceptEvent;
  bool	mAcceptQvectorFile;
  bool	mAcceptQvectorFiletmp;
  int         mCent;
  int         mRunnumber;
  float       mBField;
  StThreeVectorF mVertexPos;
  float       mEventPlane, mEventPlane1, mEventPlane2, mEventPlaneEtaPlus, mEventPlaneEtaMinus;
  float       mResolutionRandom, mResolutionEta;
  TVector2    mQ, mQ1, mQ2, mQEtaPlus, mQEtaMinus;
  TString     mQVectorDir;
  TProfile*  prfQxCentEtaPlus;
  TProfile*  prfQyCentEtaPlus;
  TProfile*  prfQxCentEtaMinus;
  TProfile*  prfQyCentEtaMinus;    
  
  float      qxTracks[maxNTracks];
  float      qyTracks[maxNTracks];
  
  ClassDef(StEventPlane, 0)
};
inline int   StEventPlane::getCentrality() const { return mCent; }
inline float StEventPlane::getEventPlane() const { return mEventPlane; }
inline float StEventPlane::getEventPlane1() const { return mEventPlane1; }
inline float StEventPlane::getEventPlane2() const { return mEventPlane2; }
inline float StEventPlane::getEventPlaneEtaPlus() const  { return mEventPlaneEtaPlus; }
inline float StEventPlane::getEventPlaneEtaMinus() const { return mEventPlaneEtaMinus; }
inline float StEventPlane::getResolutionRandom() const { return mResolutionRandom; }
inline float StEventPlane::getResolutionEta() const { return mResolutionEta; }
inline int   StEventPlane::getRunId() const { return mRunnumber; }
inline bool  StEventPlane::getAcceptEvent() const { return mAcceptQvectorFile && mAcceptQvectorFiletmp; }
#endif
