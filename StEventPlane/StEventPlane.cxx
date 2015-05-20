#include "StEventPlane.h"
#include "StPicoDstMaker/StPicoDst.h"
#include "StPicoDstMaker/StPicoEvent.h"
#include "StPicoDstMaker/StPicoTrack.h"
#include "StPicoDstMaker/StPicoDstMaker.h"
#include "StPicoDstMaker/StPicoBTofPidTraits.h"
#include "StPicoD0EventMaker/StPicoD0Event.h"
#include "StPicoD0EventMaker/StKaonPion.h"
#include "PhysicalConstants.h"
#include "StRoot/StRefMultCorr/StRefMultCorr.h"
#include "StRoot/StRefMultCorr/CentralityMaker.h"

#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TFile.h"
#include "TProfile.h"
#include "TClonesArray.h"
#include "TTree.h"

#include <iostream>
#include <fstream>
#include <string>
#include <cmath>

ClassImp(StEventPlane)

//-----------------------------------------------------------------------------
StEventPlane::StEventPlane(const char* name, StPicoDstMaker *picoMaker, StRefMultCorr* grefmultCorrUtil)
: StMaker(name), mPicoDstMaker(picoMaker), mPicoDst(0),  mgrefmultCorrUtil(grefmultCorrUtil),
  mEventPlane(0), mEventPlane1(0), mEventPlane2(0), mEventPlaneEtaPlus(0), mEventPlaneEtaMinus(0), mResolutionRandom(0), mResolutionEta(0)
{
}

//-----------------------------------------------------------------------------
StEventPlane::~StEventPlane()
{ /*  */ }

//-----------------------------------------------------------------------------
Int_t StEventPlane::Init()
{
  
  mAcceptEvent = false;
  mRunnumber = 0;
  mQVectorDir = "/global/homes/q/qiuh/myEliza17/D0v2/recenter/qVectorRun";

  //Event Cuts 
  mVzMax = 6.0;
  mDeltaVzMax = 3.0;
  
  //Track Cuts
  mNHitsFitMin = 15;
  
  //Track cuts for event plane
  mEtaMaxEventPlane = 1.0;
  mPtMinEventPlane = 0.15;
  mPtMaxEventPlane = 2.;
  mDcaMaxEventPlane = 3.0;
  

  
 // // event plane and Q vector
  float PI = TMath::Pi();

 // // D0 v2 histograms
  
  StRefMultCorr* mgrefmultCorrUtil = new StRefMultCorr("grefmult");

  return kStOK;
}

//-----------------------------------------------------------------------------
Int_t StEventPlane::Finish() {
//  mFileOut->Write();
  return kStOK;
}

//-----------------------------------------------------------------------------
void StEventPlane::Clear(Option_t *opt) {
}

//-----------------------------------------------------------------------------
Int_t StEventPlane::Make() {
  if(!mPicoDstMaker) {
    LOG_WARN << " No PicoDstMaker! Skip! " << endm;
    return kStWarn;
    }
  
  mPicoDst = mPicoDstMaker->picoDst();
  if(!mPicoDst) {
        LOG_WARN << " No PicoDst! Skip! " << endm;
        return kStWarn;
  }
  
  getEventInfo();//get event info
  
  if(mAcceptEvent){
    int eventPlaneStatus = calculateEventPlane();
    if(!eventPlaneStatus)
      {
//	calculateHadronV2();
//	calculateD0V2();
      }
    }
  
  return kStOK;
}

/*--------------------------------------------------------------------------------------------------------------------------------------------------------*/
void StEventPlane::getEventInfo()
{
  mAcceptEvent = false;
  
  if(!mPicoDst) return;
  
  //Load event
  mPicoEvent = (StPicoEvent*)mPicoDst->event();
  if(!mPicoEvent){
    cerr<<"Error opening picoDst Event, skip!"<<endl;
    return;
  }
  
  //Remove bad vertices
  mVertexPos = mPicoEvent->primaryVertex();

  if(!mgrefmultCorrUtil) {  
    LOG_WARN << " No mgrefmultCorrUtil! Skip! " << endl;
    return;
  } 
  mgrefmultCorrUtil->init(mPicoDst->event()->runId());
  mgrefmultCorrUtil->initEvent(mPicoDst->event()->grefMult(),mVertexPos.z(),mPicoDst->event()->ZDCx()) ;
  mCent  = mgrefmultCorrUtil->getCentralityBin9();

  
  mAcceptEvent = true;

  if(mRunnumber != mPicoEvent->runId())
  {
    mRunnumber = mPicoEvent->runId();
    char fileName[256];
    sprintf(fileName, "%s/%i.qVector.root", mQVectorDir.Data(), mRunnumber);
    cout<<"load qVector file: "<<fileName<<endl;
    TFile* fQVector = new TFile(fileName);
    prfQxCentEtaPlus = (TProfile*)fQVector->Get("prfQxCentEtaPlus");
    prfQyCentEtaPlus = (TProfile*)fQVector->Get("prfQyCentEtaPlus");
    prfQxCentEtaMinus = (TProfile*)fQVector->Get("prfQxCentEtaMinus");
    prfQyCentEtaMinus = (TProfile*)fQVector->Get("prfQyCentEtaMinus");
  }

  mBField = mPicoEvent->bField();
}

/*----------------------------------------------------------------------------------------------------------------------*/
int StEventPlane::calculateEventPlane()
{
  memset(qxTracks, 0, maxNTracks*sizeof(float));
  memset(qyTracks, 0, maxNTracks*sizeof(float));

  // pre-loop to count tracks for event plane, prepare for shuffle
  int nTracksForEventPlane = 0;
  for(int iTrack=0; iTrack<mPicoDst->numberOfTracks(); iTrack++)
  {
    StPicoTrack* picoTrack =(StPicoTrack*) mPicoDst->track(iTrack);
    if(!picoTrack){
      break;
    }

    if(picoTrack->nHitsFit() < mNHitsFitMin) continue;

    StPhysicalHelix* helix = &picoTrack->dcaGeometry().helix();
    float dca = helix->geometricSignedDistance(mVertexPos);
    if(TMath::Abs(dca) > mDcaMaxEventPlane) continue;

    float pathLengthToPrimaryVertex =helix->pathLength(mVertexPos.x(), mVertexPos.y());
    StThreeVectorF momentum = helix->momentumAt(pathLengthToPrimaryVertex, mBField*kilogauss);
    float pt = momentum.perp();
    float eta = momentum.pseudoRapidity();
    if(fabs(eta) > mEtaMaxEventPlane) continue;
    if(pt<mPtMinEventPlane || pt>mPtMaxEventPlane) continue;

    nTracksForEventPlane++;
  }

  int indexTrack[nTracksForEventPlane];
  int Scount = nTracksForEventPlane/2;
  for(int q=0;q<nTracksForEventPlane;q++) indexTrack[q] = q;
  random_shuffle(indexTrack,indexTrack+nTracksForEventPlane);
  int iTrackForEventPlane = 0;

  // track loop
  float Qx=0., Qy=0.;
  float Qx1=0., Qy1=0., Qx2=0., Qy2=0.;
  float QxEtaPlus=0., QyEtaPlus=0., QxEtaMinus=0., QyEtaMinus=0.;
  float vertexZ = mVertexPos.z();
  for(int iTrack=0; iTrack<mPicoDst->numberOfTracks(); iTrack++)
  {
    StPicoTrack* picoTrack =(StPicoTrack*) mPicoDst->track(iTrack);
    if(!picoTrack){
      break;
    }

//    hNHitsFit->Fill(picoTrack->nHitsFit());
    if(picoTrack->nHitsFit() < mNHitsFitMin) continue;

    StPhysicalHelix* helix = &picoTrack->dcaGeometry().helix();
    float dca = helix->geometricSignedDistance(mVertexPos);
//    hDca->Fill(dca);
    if(TMath::Abs(dca) > mDcaMaxEventPlane) continue;

    float pathLengthToPrimaryVertex =helix->pathLength(mVertexPos.x(), mVertexPos.y());
    StThreeVectorF momentum = helix->momentumAt(pathLengthToPrimaryVertex, mBField*kilogauss);
    float pt = momentum.perp();
    float eta = momentum.pseudoRapidity();
    float phi = momentum.phi();
//    hEta->Fill(eta);
//    hPt->Fill(pt);
    if(fabs(eta) > mEtaMaxEventPlane) continue;
    if(pt<mPtMinEventPlane || pt>mPtMaxEventPlane) continue;


    float qx = cos(2*phi)*pt;
    float qy = sin(2*phi)*pt;

    if(eta>0)
    {
      qx -= prfQxCentEtaPlus->GetBinContent(mCent+1);
      qy -= prfQyCentEtaPlus->GetBinContent(mCent+1);
    }
    else
    {
      qx -= prfQxCentEtaMinus->GetBinContent(mCent+1);
      qy -= prfQyCentEtaMinus->GetBinContent(mCent+1);
    }

    Qx += qx;
    Qy += qy;

    if(indexTrack[iTrackForEventPlane] >= Scount)
    {
      Qx1 += qx;
      Qy1 += qy;
    }
    else
    {
      Qx2 += qx;
      Qy2 += qy;
    }

    if(eta>0)
    {
      QxEtaPlus += qx;
      QyEtaPlus += qy;
    }
    else
    {
      QxEtaMinus += qx;
      QyEtaMinus += qy;
    }

    qxTracks[iTrack] = qx;
    qyTracks[iTrack] = qy;

    iTrackForEventPlane++;
  }//loop thru picoTracks

  assert(iTrackForEventPlane == nTracksForEventPlane);

  mQ.Set(Qx, Qy);
  mQ1.Set(Qx1, Qy1);
  mQ2.Set(Qx2, Qy2);
  mQEtaPlus.Set(QxEtaPlus, QyEtaPlus);
  mQEtaMinus.Set(QxEtaMinus, QyEtaMinus);

  if(mQ.Mod2()==0 || mQ1.Mod2()==0 || mQ2.Mod2()==0 || mQEtaPlus.Mod2()==0 || mQEtaMinus.Mod2()==0)
  {
    cout<<"0 Q: "<<mQ.Mod2()<<" "<<mQ1.Mod2()<<" "<<mQ2.Mod2()<<" "<<mQEtaPlus.Mod2()<<" "<<mQEtaMinus.Mod2()<<" "<<mCent<<" "<<nTracksForEventPlane<<" "<<mPicoEvent->refMult()<<" "<<mPicoEvent->grefMult()<<" "<<mPicoDst->numberOfTracks()<<" "<<mPicoEvent->eventId()<<endl;
    return 1;
  }

  mEventPlane = mQ.Phi()*0.5;
  mEventPlane1 = mQ1.Phi()*0.5;
  mEventPlane2 = mQ2.Phi()*0.5;
  mEventPlaneEtaPlus = mQEtaPlus.Phi()*0.5;
  mEventPlaneEtaMinus = mQEtaMinus.Phi()*0.5;
//  hQyQxCent->Fill(mCent, Qx, Qy);
//  hEventPlaneCent->Fill(mCent, mEventPlane);

//  hQyQx1Cent->Fill(mCent, Qx1, Qy1);
//  hEventPlane1Cent->Fill(mCent, mEventPlane1);
//
//  hQyQx2Cent->Fill(mCent, Qx2, Qy2);
//  hEventPlane2Cent->Fill(mCent, mEventPlane2);
//
//  hQyQxEtaPlusCent->Fill(mCent, QxEtaPlus, QyEtaPlus);
//  hEventPlaneEtaPlusCent->Fill(mCent, mEventPlaneEtaPlus);
//
//  hQyQxEtaMinusCent->Fill(mCent, QxEtaMinus, QyEtaMinus);
//  hEventPlaneEtaMinusCent->Fill(mCent, mEventPlaneEtaMinus);
//
//  prfCosResolutionRandomCent->Fill(mCent, cos(2.*(mEventPlane1-mEventPlane2)));
//  prfCosResolutionEtaCent->Fill(mCent, cos(2.*(mEventPlaneEtaPlus-mEventPlaneEtaMinus)));
  mResolutionRandom = cos(2.*(mEventPlane1-mEventPlane2));
  mResolutionEta = cos(2.*(mEventPlaneEtaPlus-mEventPlaneEtaMinus));

  return 0;
}

