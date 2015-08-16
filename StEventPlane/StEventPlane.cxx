#include <iostream>
#include <cmath>

#include "TFile.h"
#include "TH1I.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
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
  if(mFileOut)
    {
      mFileOut->cd();
      // track level QA 
      hNHitsFit = new TH1I("hNHitsFit", "hNHitsFit", 50, 0, 50);
      hDca = new TH1F("hDca", "hDca", 20000, -10, 10);
      hEta = new TH1F("hEta", "hEta", 30, -1.5, 1.5);
      hPt = new TH1F("hPt", "hPt", 200, 0, 10);

      // event plane and Q vector 
      float PI = TMath::Pi();

      hPhiCentEtaPlusZPlus = new TH2F("hPhiCentEtaPlusZPlus","hPhiCentEtaPlusZPlus",9,0,9,120,-PI,PI);
      hPhiCentEtaPlusZMinus = new TH2F("hPhiCentEtaPlusZMinus","hPhiCentEtaPlusZMinus",9,0,9,120,-PI,PI);
      hPhiCentEtaMinusZPlus = new TH2F("hPhiCentEtaMinusZPlus","hPhiCentEtaMinusZPlus",9,0,9,120,-PI,PI);
      hPhiCentEtaMinusZMinus = new TH2F("hPhiCentEtaMinusZMinus","hPhiCentEtaMinusZMinus",9,0,9,120,-PI,PI);

      hEventPlaneCent = new TH2F("hEventPlaneCent","hEventPlaneCent",9,0,9,60,0,PI);
      hEventPlane1Cent = new TH2F("hEventPlane1Cent","hEventPlane1Cent",9,0,9,60,0,PI);
      hEventPlane2Cent = new TH2F("hEventPlane2Cent","hEventPlane2Cent",9,0,9,60,0,PI);
      hEventPlaneEtaPlusCent = new TH2F("hEventPlaneEtaPlusCent","hEventPlaneEtaPlusCent",9,0,9,60,0,PI);
      hEventPlaneEtaMinusCent = new TH2F("hEventPlaneEtaMinusCent","hEventPlaneEtaMinusCent",9,0,9,60,0,PI);
      hQyQxCent = new TH3F("hQyQxCent", "hQyQxCent", 9,0,9,1000,-50,50,1000,-50,50);
      hQyQx1Cent = new TH3F("hQyQx1Cent", "hQyQx1Cent", 9,0,9,1000,-50,50,1000,-50,50);
      hQyQx2Cent = new TH3F("hQyQx2Cent", "hQyQx2Cent", 9,0,9,1000,-50,50,1000,-50,50);
      hQyQxEtaPlusCent = new TH3F("hQyQxEtaPlusCent", "hQyQxEtaPlusCent", 9,0,9,1000,-50,50,1000,-50,50);
      hQyQxEtaMinusCent = new TH3F("hQyQxEtaMinusCent", "hQyQxEtaMinusCent", 9,0,9,1000,-50,50,1000,-50,50);
      prfCosResolutionRandomCent = new TProfile("prfCosResolutionRandomCent", "prfCosResolutionRandomCent", 9, 0, 9);
      prfCosResolutionEtaCent = new TProfile("prfCosResolutionEtaCent", "prfCosResolutionEtaCent", 9, 0, 9);

      hHadronV2PtCent = new TH3F("hHadronV2PtCent", "hHadronV2PtCent", 100, 0., 5., 9, 0, 9, 200, -1., 1.);
      hHadronHftV2PtCent = new TH3F("hHadronHftV2PtCent", "hHadronHftV2PtCent", 100, 0., 5., 9, 0, 9, 200, -1., 1.);
      hHadronPrimaryV2PtCent = new TH3F("hHadronPrimaryV2PtCent", "hHadronPrimaryV2PtCent", 100, 0., 5., 9, 0, 9, 200, -1., 1.);
      hHadronHftPrimaryV2PtCent = new TH3F("hHadronHftPrimaryV2PtCent", "hHadronHftPrimaryV2PtCent", 100, 0., 5., 9, 0, 9, 200, -1., 1.);
     
    }
  return kStOk;
}
//-----------------------------------------------------------------------------
void StEventPlane::setFileOut(TFile* fileOut)
{
  mFileOut = fileOut;
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
       mEventPlaneStatus = calculateEventPlane();
       if (!mEventPlaneStatus && mAcceptEvent)
	 {
	   calculateHadronV2();
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

   bool isVPDMB5=kFALSE;
   for(int i=0;i<9;i++) { if(mPicoEvent->triggerWord() & (1<<i)) isVPDMB5=kTRUE ;}//Select MB trigger
   if (!(isVPDMB5)) {
     return;
   }

   if(TMath::Abs(mVertexPos.z()) > EventPlaneConstants::vzMax) return;
   if(TMath::Abs(mVertexPos.z() - mPicoEvent->vzVpd()) > EventPlaneConstants::deltaVzMax) return;
   if(mCent<0 || mCent>9) return;
   
   mAcceptEvent = true;
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

   prfQxCentEtaPlus->SetDirectory(0);
   prfQyCentEtaPlus->SetDirectory(0);
   prfQxCentEtaMinus->SetDirectory(0);
   prfQyCentEtaMinus->SetDirectory(0);

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
float StEventPlane::getEventPlane(int nTracksToExclude, int* indexTracksToExclude) const 
{
  TVector2 Qsub = mQ;
  for(int i=0; i<nTracksToExclude; i++)
    {
      TVector2 qTrack(qxTracks[indexTracksToExclude[i]], qyTracks[indexTracksToExclude[i]]);
      Qsub -= qTrack;
    }
  return Qsub.Phi() * 0.5;
}
void StEventPlane::calculateHadronV2() const
{
  for(int iTrack=0; iTrack<mPicoDst->numberOfTracks(); iTrack++)
    {
      StPicoTrack* picoTrack =(StPicoTrack*) mPicoDst->track(iTrack);
      if(!picoTrack){
        break;
      }

      if(picoTrack->nHitsFit() < EventPlaneConstants::nHitsFitMin) continue;

      StPhysicalHelix* helix = &picoTrack->dcaGeometry().helix();
      float dca = helix->geometricSignedDistance(mVertexPos);
      if(TMath::Abs(dca) > EventPlaneConstants::dcaMaxEventPlane) continue;

      float pathLengthToPrimaryVertex =helix->pathLength(mVertexPos.x(), mVertexPos.y());
      StThreeVectorF momentum = helix->momentumAt(pathLengthToPrimaryVertex, mBField*kilogauss);
      float pt = momentum.perp();
      float eta = momentum.pseudoRapidity();
      float phi = momentum.phi();
      if(fabs(eta) > EventPlaneConstants::etaMaxEventPlane) continue;

      float qx = qxTracks[iTrack];
      float qy = qyTracks[iTrack];
      TVector2 qTrack(qx, qy);
      TVector2 QSub = mQ - qTrack;
      float psi = QSub.Phi()/2;

      float weight = mgrefmultCorrUtil->getWeight();

      hHadronV2PtCent->Fill(pt, mCent, cos(2.*(phi-psi)), weight);
      
      if(picoTrack->isHFTTrack())
	hHadronHftV2PtCent->Fill(pt, mCent, cos(2.*(phi-psi)), weight);
      
      if(picoTrack->pMom().mag()>0)
	hHadronPrimaryV2PtCent->Fill(pt, mCent, cos(2.*(phi-psi)), weight);
      
      if(picoTrack->isHFTTrack() && picoTrack->pMom().mag()>0)
	hHadronHftPrimaryV2PtCent->Fill(pt, mCent, cos(2.*(phi-psi)), weight);

    }
}

