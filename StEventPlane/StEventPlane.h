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
  void getEventInfo();
  int calculateEventPlane();

  Int_t   getCentrality() const;
  Float_t getEventPlane() const;
  Float_t getEventPlane1() const;
  Float_t getEventPlane2() const;
  Float_t getEventPlaneEtaPlus() const;
  Float_t getEventPlaneEtaMinus() const;
  Float_t getResolutionRandom() const;
  Float_t getResolutionEta() const;

  Int_t   getRunId() const;
  Bool_t getAcceptEvent() const;

 private:
  StPicoDstMaker* mPicoDstMaker;
  StPicoDst*      mPicoDst;
  StPicoEvent*    mPicoEvent;
  StRefMultCorr* mgrefmultCorrUtil;
  
  bool	mAcceptEvent;
  bool	mAcceptQvectorFile;
  bool	mAcceptQvectorFiletmp;
  Int_t         mCent;
  Int_t         mRunnumber;
  Float_t       mBField;
  StThreeVectorF mVertexPos;
  Float_t       mEventPlane, mEventPlane1, mEventPlane2, mEventPlaneEtaPlus, mEventPlaneEtaMinus;
  Float_t       mResolutionRandom, mResolutionEta;
  TVector2    mQ, mQ1, mQ2, mQEtaPlus, mQEtaMinus;
  TString     mQVectorDir;
  TProfile*  prfQxCentEtaPlus;
  TProfile*  prfQyCentEtaPlus;
  TProfile*  prfQxCentEtaMinus;
  TProfile*  prfQyCentEtaMinus;    
  
  //Event cuts
  Float_t mVzMax;
  Float_t mRefMultMin;
  Float_t mDeltaVzMax;
  
  //Track cuts
  UShort_t mNHitsFitMin;
  
  //Track cuts for event plane
  Float_t mEtaMaxEventPlane;
  Float_t mPtMinEventPlane;
  Float_t mPtMaxEventPlane;
  Float_t mDcaMaxEventPlane;
  
  Float_t      qxTracks[maxNTracks];
  Float_t      qyTracks[maxNTracks];
  
  ClassDef(StEventPlane, 0)
};
inline Int_t   StEventPlane::getCentrality() const { return (Int_t) mCent; }
inline Float_t StEventPlane::getEventPlane() const { return (Float_t) mEventPlane; }
inline Float_t StEventPlane::getEventPlane1() const { return (Float_t) mEventPlane1; }
inline Float_t StEventPlane::getEventPlane2() const { return (Float_t) mEventPlane2; }
inline Float_t StEventPlane::getEventPlaneEtaPlus() const { return (Float_t) mEventPlaneEtaPlus; }
inline Float_t StEventPlane::getEventPlaneEtaMinus() const { return (Float_t) mEventPlaneEtaMinus; }
inline Float_t StEventPlane::getResolutionRandom() const { return (Float_t) mResolutionRandom; }
inline Float_t StEventPlane::getResolutionEta() const { return (Float_t) mResolutionEta; }
inline Int_t   StEventPlane::getRunId() const { return (Int_t) mRunnumber; }
inline Bool_t  StEventPlane::getAcceptEvent() const { return (Bool_t)  mAcceptQvectorFile && mAcceptQvectorFiletmp; }
#endif
