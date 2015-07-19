#include <iostream>
#include <cmath>

#include "TFile.h"
#include "TProfile.h"

#include "StPicoDstMaker/StPicoDst.h"
#include "StPicoDstMaker/StPicoEvent.h"
#include "StPicoDstMaker/StPicoTrack.h"
#include "StPicoDstMaker/StPicoDstMaker.h"
#include "StPicoDstMaker/StPicoBTofPidTraits.h"
#include "PhysicalConstants.h"
#include "StRoot/StRefMultCorr/StRefMultCorr.h"
#include "StRoot/StRefMultCorr/CentralityMaker.h"

#include "StEventPlaneConstants.h"
#include "StEventPlane.h"

ClassImp(StEventPlane)

//-----------------------------------------------------------------------------
StEventPlane::StEventPlane(const char* name, StPicoDstMaker *picoMaker, StRefMultCorr* grefmultCorrUtil)
   : StMaker(name), mPicoDstMaker(picoMaker), mPicoDst(NULL),  mPicoEvent(NULL), mgrefmultCorrUtil(grefmultCorrUtil),
     mAcceptEvent(false), mAcceptQvectorFile(false), mAcceptQvectorFiletmp(true), mCent(-1), mRunNumber(0), mBField(-999.),mVertexPos(-999,-999,-999),
     mEventPlane(0), mEventPlane1(0), mEventPlane2(0), mEventPlaneEtaPlus(0), mEventPlaneEtaMinus(0), mResolutionRandom(0), mResolutionEta(0),
     mQ(-999,-999), mQ1(-999,-999), mQ2(-999,-999), mQEtaPlus(-999,-999), mQEtaMinus(-999,-999),
     prfQxCentEtaPlus(NULL), prfQyCentEtaPlus(NULL), prfQxCentEtaMinus(NULL), prfQyCentEtaMinus(NULL)
{
}

//-----------------------------------------------------------------------------
Int_t StEventPlane::Init()
{
   // StRefMultCorr* mgrefmultCorrUtil = new StRefMultCorr("grefmult");

   return kStOK;
}

//-----------------------------------------------------------------------------
Int_t StEventPlane::Make()
{
   if (!mPicoDstMaker)
   {
      LOG_ERROR << " No PicoDstMaker! Skip! " << endm;
      return kStErr;
   }

   mPicoDst = mPicoDstMaker->picoDst();
   if (!mPicoDst)
   {
      LOG_ERROR << " No PicoDst! Skip! " << endm;
      return kStErr;
   }

   mPicoEvent = (StPicoEvent*)mPicoDst->event();
   if (!mPicoEvent)
   {
      LOG_ERROR << "Error opening picoDst Event, skip!" << endm;
      return kStErr;
   }
   if (mRunNumber != mPicoEvent->runId()) getRunInfo(mPicoEvent->runId());
   else mAcceptQvectorFile = true;

   getEventInfo();//get event info

 //  if (mAcceptEvent)
   if (mAcceptQvectorFile && mAcceptQvectorFiletmp)
   {
      int eventPlaneStatus = calculateEventPlane();
      if (!eventPlaneStatus)
      {
// calculateHadronV2();
// calculateD0V2();
      }
   }

   return kStOK;
}

/*--------------------------------------------------------------------------------------------------------------------------------------------------------*/
void StEventPlane::getEventInfo()
{
   mAcceptEvent = false;

   //Remove bad vertices
   mVertexPos = mPicoEvent->primaryVertex();

   mgrefmultCorrUtil->init(mPicoDst->event()->runId());
   mgrefmultCorrUtil->initEvent(mPicoDst->event()->grefMult(), mVertexPos.z(), mPicoDst->event()->ZDCx()) ;
   mCent  = mgrefmultCorrUtil->getCentralityBin9();

   mAcceptEvent = true;

   mBField = mPicoEvent->bField();
}

void StEventPlane::getRunInfo(int const runNumber)
{
   mRunNumber = runNumber;
   char fileName[256];
   sprintf(fileName, "%s/%i.qVector.root", EventPlaneConstants::qVectorDir.Data(), mRunNumber);
   cout << "load qVector file: " << fileName << endl;
   TFile fQVector(fileName);

   fQVector.GetObject("prfQxCentEtaPlus", prfQxCentEtaPlus);
   if (!prfQxCentEtaPlus)
   {
     LOG_INFO << "StEventPlane::THistograms and TProiles NOT found! shoudl check the files From HaoQiu" << endm;
     mAcceptQvectorFile = false;
     mAcceptQvectorFiletmp = false;
     return;
   }
   else
   {
     mAcceptQvectorFile = true;
     mAcceptQvectorFiletmp = true;
   }

   prfQxCentEtaPlus =  (TProfile*)fQVector.Get("prfQxCentEtaPlus")->Clone("prfQxCentEtaPlus");
   prfQyCentEtaPlus =  (TProfile*)fQVector.Get("prfQyCentEtaPlus")->Clone("prfQyCentEtaPlus");
   prfQxCentEtaMinus = (TProfile*)fQVector.Get("prfQxCentEtaMinus")->Clone("prfQxCentEtaMinus");
   prfQyCentEtaMinus = (TProfile*)fQVector.Get("prfQyCentEtaMinus")->Clone("prfQyCentEtaMinus");

   fQVector.Close();
}

/*----------------------------------------------------------------------------------------------------------------------*/
int StEventPlane::calculateEventPlane()
{
   memset(qxTracks, 0, maxNTracks * sizeof(float));
   memset(qyTracks, 0, maxNTracks * sizeof(float));

   // pre-loop to count tracks for event plane, prepare for shuffle
   int nTracksForEventPlane = 0;
   for (int iTrack = 0; iTrack < mPicoDst->numberOfTracks(); ++iTrack)
   {
      StPicoTrack* picoTrack = (StPicoTrack*) mPicoDst->track(iTrack);
      if (!picoTrack)
      {
         break;
      }

      if (picoTrack->nHitsFit() < EventPlaneConstants::nHitsFitMin) continue;

      StPhysicalHelix* helix = &picoTrack->dcaGeometry().helix();
      float dca = helix->geometricSignedDistance(mVertexPos);
      if (TMath::Abs(dca) > EventPlaneConstants::dcaMaxEventPlane) continue;

      float pathLengthToPrimaryVertex = helix->pathLength(mVertexPos.x(), mVertexPos.y());
      StThreeVectorF momentum = helix->momentumAt(pathLengthToPrimaryVertex, mBField * kilogauss);
      float pt = momentum.perp();
      float eta = momentum.pseudoRapidity();
      if (fabs(eta) > EventPlaneConstants::etaMaxEventPlane) continue;
      if (pt < EventPlaneConstants::ptMinEventPlane || pt > EventPlaneConstants::ptMaxEventPlane) continue;

      nTracksForEventPlane++;
   }

   int indexTrack[nTracksForEventPlane];
   int Scount = nTracksForEventPlane / 2;
   for (int q = 0; q < nTracksForEventPlane; ++q) indexTrack[q] = q;
   random_shuffle(indexTrack, indexTrack + nTracksForEventPlane);
   int iTrackForEventPlane = 0;

   // track loop
   float Qx = 0., Qy = 0.;
   float Qx1 = 0., Qy1 = 0., Qx2 = 0., Qy2 = 0.;
   float QxEtaPlus = 0., QyEtaPlus = 0., QxEtaMinus = 0., QyEtaMinus = 0.;
   float vertexZ = mVertexPos.z();
   for (int iTrack = 0; iTrack < mPicoDst->numberOfTracks(); ++iTrack)
   {
      StPicoTrack* picoTrack = (StPicoTrack*) mPicoDst->track(iTrack);
      if (!picoTrack)
      {
         break;
      }

      if (picoTrack->nHitsFit() < EventPlaneConstants::nHitsFitMin) continue;

      StPhysicalHelix* helix = &picoTrack->dcaGeometry().helix();
      float dca = helix->geometricSignedDistance(mVertexPos);
      if (TMath::Abs(dca) > EventPlaneConstants::dcaMaxEventPlane) continue;

      float pathLengthToPrimaryVertex = helix->pathLength(mVertexPos.x(), mVertexPos.y());
      StThreeVectorF momentum = helix->momentumAt(pathLengthToPrimaryVertex, mBField * kilogauss);
      float pt = momentum.perp();
      float eta = momentum.pseudoRapidity();
      float phi = momentum.phi();
      if (fabs(eta) > EventPlaneConstants::etaMaxEventPlane) continue;
      if (pt < EventPlaneConstants::ptMinEventPlane || pt > EventPlaneConstants::ptMaxEventPlane) continue;


      float qx = cos(2 * phi) * pt;
      float qy = sin(2 * phi) * pt;

      if (eta > 0)
      {
         qx -= prfQxCentEtaPlus->GetBinContent(mCent + 1);
         qy -= prfQyCentEtaPlus->GetBinContent(mCent + 1);
      }
      else
      {
         qx -= prfQxCentEtaMinus->GetBinContent(mCent + 1);
         qy -= prfQyCentEtaMinus->GetBinContent(mCent + 1);
      }

      Qx += qx;
      Qy += qy;

      if (indexTrack[iTrackForEventPlane] >= Scount)
      {
         Qx1 += qx;
         Qy1 += qy;
      }
      else
      {
         Qx2 += qx;
         Qy2 += qy;
      }

      if (eta > 0)
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

   if (mQ.Mod2() == 0 || mQ1.Mod2() == 0 || mQ2.Mod2() == 0 || mQEtaPlus.Mod2() == 0 || mQEtaMinus.Mod2() == 0)
   {
      // cout << "0 Q: " << mQ.Mod2() << " " << mQ1.Mod2() << " " << mQ2.Mod2() << " " << mQEtaPlus.Mod2() << " " << mQEtaMinus.Mod2() << " " << mCent << " " << nTracksForEventPlane << " " << mPicoEvent->refMult() << " " << mPicoEvent->grefMult() << " " << mPicoDst->numberOfTracks() << " " << mPicoEvent->eventId() << " " << mPicoDst->event()->primaryVertex().z() << " " << mPicoDst->event()->vzVpd() << " " << mPicoDst->event()->ZDCx() / 1.0e3 << endl;
      return 1;
   }

   mEventPlane = mQ.Phi() * 0.5;
   mEventPlane1 = mQ1.Phi() * 0.5;
   mEventPlane2 = mQ2.Phi() * 0.5;
   mEventPlaneEtaPlus = mQEtaPlus.Phi() * 0.5;
   mEventPlaneEtaMinus = mQEtaMinus.Phi() * 0.5;
   mResolutionRandom = cos(2.*(mEventPlane1 - mEventPlane2));
   mResolutionEta = cos(2.*(mEventPlaneEtaPlus - mEventPlaneEtaMinus));

   return 0;
}

