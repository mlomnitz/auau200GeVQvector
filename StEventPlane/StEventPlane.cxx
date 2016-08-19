#define only__primary
#include <iostream>
#include <cmath>

#include "TFile.h"
#include "TH1I.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
#include "THn.h"
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
StEventPlane::StEventPlane(const char* name, StPicoDstMaker *picoMaker, StRefMultCorr* grefmultCorrUtil, int harmonic)
   : StMaker(name), mPicoDstMaker(picoMaker), mPicoDst(NULL),  mPicoEvent(NULL), mgrefmultCorrUtil(grefmultCorrUtil),
  mAcceptEvent(false), mAcceptQvectorFile(false), mAcceptQvectorFiletmp(true), mHarmonic(harmonic), mCent(-1), mRunNumber(0), mBField(-999.), 
     mVertexPos(-999, -999, -999), mEventPlane(0), mEventPlane1(0), mEventPlane2(0), mEventPlaneEtaPlus(0), mEventPlaneEtaMinus(0), 
     mResolutionRandom(0), mResolutionEta(0), mQ(-999, -999), mQ1(-999, -999), mQ2(-999, -999), mQEtaPlus(-999, -999), mQEtaMinus(-999, -999),
     prfQxCentEtaPlus(NULL), prfQyCentEtaPlus(NULL), prfQxCentEtaMinus(NULL), prfQyCentEtaMinus(NULL)
{
  TFile* outFile = new TFile(Form("%s.eventplane.root",name),"RECREATE");
  setFileOut(outFile);
}
//-----------------------------------------------------------------------------
Int_t StEventPlane::Init()
{
  //cout<<"Lomnitz : pre enter"<<endl;
   if (mFileOut)
   {
     // cout<<"Lomnitz if passed"<<endl;
      mFileOut->cd();
      // cout<<"Lomnitz all passed"<<endl;
      // track level QA
      hNHitsFit = new TH1I("hNHitsFit", "hNHitsFit", 50, 0, 50);
      hDca = new TH1F("hDca", "hDca", 20000, -10, 10);
      hEta = new TH1F("hEta", "hEta", 30, -1.5, 1.5);
      hPt = new TH1F("hPt", "hPt", 200, 0, 10);

      // event plane and Q vector
      float PI = TMath::Pi();

      hPhiCentEtaPlusZPlus = new TH2F("hPhiCentEtaPlusZPlus", "hPhiCentEtaPlusZPlus", 9, 0, 9, 120, -PI, PI);
      hPhiCentEtaPlusZMinus = new TH2F("hPhiCentEtaPlusZMinus", "hPhiCentEtaPlusZMinus", 9, 0, 9, 120, -PI, PI);
      hPhiCentEtaMinusZPlus = new TH2F("hPhiCentEtaMinusZPlus", "hPhiCentEtaMinusZPlus", 9, 0, 9, 120, -PI, PI);
      hPhiCentEtaMinusZMinus = new TH2F("hPhiCentEtaMinusZMinus", "hPhiCentEtaMinusZMinus", 9, 0, 9, 120, -PI, PI);

      hEventPlaneCent = new TH2F("hEventPlaneCent", "hEventPlaneCent", 9, 0, 9, 60, 0, (2.0/mHarmonic)*PI);
      hEventPlane1Cent = new TH2F("hEventPlane1Cent", "hEventPlane1Cent", 9, 0, 9, 60, 0, (2.0/mHarmonic)*PI);
      hEventPlane2Cent = new TH2F("hEventPlane2Cent", "hEventPlane2Cent", 9, 0, 9, 60, 0, (2.0/mHarmonic)*PI);
      hEventPlaneEtaPlusCent = new TH2F("hEventPlaneEtaPlusCent", "hEventPlaneEtaPlusCent", 9, 0, 9, 60, 0, (2.0/mHarmonic)*PI);
      hEventPlaneEtaMinusCent = new TH2F("hEventPlaneEtaMinusCent", "hEventPlaneEtaMinusCent", 9, 0, 9, 60, 0, (2.0/mHarmonic)*PI);
      hQyQxCent = new TH3F("hQyQxCent", "hQyQxCent", 9, 0, 9, 1000, -50, 50, 1000, -50, 50);
      hQyQx1Cent = new TH3F("hQyQx1Cent", "hQyQx1Cent", 9, 0, 9, 1000, -50, 50, 1000, -50, 50);
      hQyQx2Cent = new TH3F("hQyQx2Cent", "hQyQx2Cent", 9, 0, 9, 1000, -50, 50, 1000, -50, 50);
      hQyQxEtaPlusCent = new TH3F("hQyQxEtaPlusCent", "hQyQxEtaPlusCent", 9, 0, 9, 1000, -50, 50, 1000, -50, 50);
      hQyQxEtaMinusCent = new TH3F("hQyQxEtaMinusCent", "hQyQxEtaMinusCent", 9, 0, 9, 1000, -50, 50, 1000, -50, 50);
      prfCosResolutionRandomCent = new TProfile("prfCosResolutionRandomCent", "prfCosResolutionRandomCent", 9, 0, 9);
      prfCosResolutionEtaCent = new TProfile("prfCosResolutionEtaCent", "prfCosResolutionEtaCent", 9, 0, 9);

      hHadronVnPtCent = new TH3F(Form("hHadronV%iPtCent",mHarmonic), Form("hHadronV%iPtCent",mHarmonic), 100, 0., 5., 9, 0, 9, 200, -1., 1.);
      hHadronHftVnPtCent = new TH3F(Form("hHadronHftV%iPtCent",mHarmonic), Form("hHadronHftV%iPtCent",mHarmonic), 100, 0., 5., 9, 0, 9, 200, -1., 1.);
      hHadronPrimaryVnPtCent = new TH3F(Form("hHadronPrimaryV%iPtCent",mHarmonic), Form("hHadronPrimaryV%iPtCent",mHarmonic), 100, 0., 5., 9, 0, 9, 200, -1., 1.);
      hHadronHftPrimaryVnPtCent = new TH3F(Form("hHadronHftPrimaryV%iPtCent",mHarmonic), Form("hHadronHftPrimaryV%iPtCent",mHarmonic), 100, 0., 5., 9, 0, 9, 200, -1., 1.);

      const int nDim = 4;
      int nBins[nDim] = {100, 9, 200, 8};//pt, cent, v2, etaGap 
      double xMin[nDim] = {0, 0, -1, 0};
      double xMax[nDim] = {5, 9, 1, 0.8};
      hHadronVnPtCentEtaGap = new THnF(Form("hHadronV%iPtCentEtaGap",mHarmonic), Form("hHadronV%iPtCentEtaGap",mHarmonic), nDim, nBins, xMin, xMax);
   }
   return kStOk;
}
//----------------------------------------------------------------------------- 
Int_t StEventPlane::Finish()
{
  cout<<"StEventPlane::Finish()"<<endl;
  mFileOut->cd();

  hPhiCentEtaPlusZPlus->Write();
  hPhiCentEtaPlusZMinus->Write();
  hPhiCentEtaMinusZPlus->Write();
  hPhiCentEtaMinusZMinus->Write();

  hEventPlaneCent->Write();
  hEventPlane1Cent->Write();
  hEventPlane2Cent->Write();
  hEventPlaneEtaPlusCent->Write();
  hEventPlaneEtaMinusCent->Write();
  hQyQxCent->Write();
  hQyQx1Cent->Write();
  hQyQx2Cent->Write();
  hQyQxEtaPlusCent->Write();
  hQyQxEtaMinusCent->Write();
  prfCosResolutionRandomCent->Write();
  prfCosResolutionEtaCent->Write();
  
  hHadronVnPtCent->Write();
  hHadronHftVnPtCent->Write();
  hHadronPrimaryVnPtCent->Write();
  hHadronHftPrimaryVnPtCent->Write();
  
  hHadronVnPtCentEtaGap->Write();

  //  mFileOut->Write();
  //  mFileOut->Close();

  return kStOK;
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

   if (mAcceptQvectorFile && mAcceptQvectorFiletmp)
   {
      mEventPlaneStatus = calculateEventPlane();
      if (!mEventPlaneStatus && mAcceptEvent)
      {
         calculateHadronVn();
      }
   }

   return kStOK;
}

/*--------------------------------------------------------------------------------------------------------------------------------------------------------*/
void StEventPlane::getEventInfo()
{

   //Remove bad vertices
   mVertexPos = mPicoEvent->primaryVertex();

   mgrefmultCorrUtil->init(mPicoDst->event()->runId());
   mgrefmultCorrUtil->initEvent(mPicoDst->event()->grefMult(), mVertexPos.z(), mPicoDst->event()->ZDCx()) ;
   mCent  = mgrefmultCorrUtil->getCentralityBin9();

   mAcceptEvent = false;

   mBField = mPicoEvent->bField();

   bool isVPDMB5 = kFALSE;
   for (int i = 0; i < EventPlaneConstants::nTrig; i++)
   {
     if ( mPicoEvent->isTrigger(EventPlaneConstants::mTriggerId[i]) ) isVPDMB5 = kTRUE ;  //Select MB trigger
   }
   if (!(isVPDMB5))
   {
      return;
   }

   if (TMath::Abs(mVertexPos.z()) > EventPlaneConstants::vzMax) return;
   if (TMath::Abs(mVertexPos.z() - mPicoEvent->vzVpd()) > EventPlaneConstants::deltaVzMax) return;
   if (mCent < 0 || mCent > 9) return;

   mAcceptEvent = true;
}

void StEventPlane::getRunInfo(int const runNumber)
{
   mRunNumber = runNumber;

   char fileName[256];
   sprintf(fileName, "%s/%i.qVector.root", EventPlaneConstants::qVectorRunDir.Data(), mRunNumber);
   TFile* fQVector = new TFile(fileName);
   if(fQVector->IsZombie())
     {
       int dayNumber = mRunNumber%1000000/1000;
       sprintf(fileName, "%s/%03d.qVector.root", EventPlaneConstants::qVectorDayDir.Data(), dayNumber);
       delete fQVector;
       fQVector = new TFile(fileName);
       if(fQVector->IsZombie())
	 {
	   cout<<"can not load run or day qVector file: "<<mRunNumber<<endl;
	   return;
	 }
     }
   cout << "load qVector file: " << fileName << endl;

   fQVector->GetObject(Form("prfQxCentEtaPlus_v%i",mHarmonic), prfQxCentEtaPlus);
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

   prfQxCentEtaPlus = (TProfile*)fQVector->Get(Form("prfQxCentEtaPlus_v%i",mHarmonic))->Clone(Form("prfQxCentEtaPlus_v%i",mHarmonic));
   prfQyCentEtaPlus = (TProfile*)fQVector->Get(Form("prfQyCentEtaPlus_v%i",mHarmonic))->Clone(Form("prfQyCentEtaPlus_v%i",mHarmonic));
   prfQxCentEtaMinus = (TProfile*)fQVector->Get(Form("prfQxCentEtaMinus_v%i",mHarmonic))->Clone(Form("prfQxCentEtaMinus_v%i",mHarmonic));
   prfQyCentEtaMinus = (TProfile*)fQVector->Get(Form("prfQyCentEtaMinus_v%i",mHarmonic))->Clone(Form("prfQyCentEtaMinus_v%i",mHarmonic));

   prfQxCentEtaPlus->SetDirectory(0);
   prfQyCentEtaPlus->SetDirectory(0);
   prfQxCentEtaMinus->SetDirectory(0);
   prfQyCentEtaMinus->SetDirectory(0);

   fQVector->Close();
   delete fQVector;
}

/*----------------------------------------------------------------------------------------------------------------------*/
int StEventPlane::calculateEventPlane()
{
   memset(qxTracks, 0, maxNTracks * sizeof(float));
   memset(qyTracks, 0, maxNTracks * sizeof(float));

   // pre-loop to count tracks for event plane, prepare for shuffle
   int nTracksForEventPlane = 0;
   for (unsigned short iTrack = 0; iTrack < mPicoDst->numberOfTracks(); ++iTrack)
   {
      StPicoTrack* picoTrack = (StPicoTrack*) mPicoDst->track(iTrack);
      if (!picoTrack)
      {
         break;
      }

      if (picoTrack->nHitsFit() < EventPlaneConstants::nHitsFitMin) continue;
      if(1.*picoTrack->nHitsFit()/picoTrack->nHitsMax() < EventPlaneConstants::mNHitsFitRatioMin) continue;

      StPhysicalHelix helix = picoTrack->dcaGeometry().helix();
      float dca = helix.geometricSignedDistance(mVertexPos);
      if (TMath::Abs(dca) > EventPlaneConstants::dcaMaxEventPlane) continue;

      float pathLengthToPrimaryVertex = helix.pathLength(mVertexPos.x(), mVertexPos.y());
      StThreeVectorF momentum = helix.momentumAt(pathLengthToPrimaryVertex, mBField * kilogauss);
      //SL16d 
      double eta, pt;
      if( picoTrack->pMom().perp() > 0 ){
	eta =picoTrack->pMom().pseudoRapidity();
	pt =picoTrack->pMom().perp();
      }
      else{
#ifdef only__primary
	continue;
#endif
	pt = momentum.perp();
	eta = momentum.pseudoRapidity();
      }
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
   float QxEtaPlusGap005 = 0., QyEtaPlusGap005 = 0., QxEtaMinusGap005 = 0., QyEtaMinusGap005 = 0.;
   float QxEta[20], QyEta[20];
   memset(QxEta, 0, 20*sizeof(float));
   memset(QyEta, 0, 20*sizeof(float));
   float vertexZ = mVertexPos.z();
   for (unsigned short iTrack = 0; iTrack < mPicoDst->numberOfTracks(); ++iTrack)
   {
      StPicoTrack* picoTrack = (StPicoTrack*) mPicoDst->track(iTrack);
      if (!picoTrack)
      {
         break;
      }

      if (mAcceptEvent) hNHitsFit->Fill(picoTrack->nHitsFit());
      if (picoTrack->nHitsFit() < EventPlaneConstants::nHitsFitMin) continue;
      if(1.*picoTrack->nHitsFit()/picoTrack->nHitsMax() < EventPlaneConstants::mNHitsFitRatioMin) continue;

      StPhysicalHelix helix = picoTrack->dcaGeometry().helix();
      float dca = helix.geometricSignedDistance(mVertexPos);
      if (mAcceptEvent) hDca->Fill(dca);
      if (TMath::Abs(dca) > EventPlaneConstants::dcaMaxEventPlane) continue;

      float pathLengthToPrimaryVertex = helix.pathLength(mVertexPos.x(), mVertexPos.y());
      StThreeVectorF momentum = helix.momentumAt(pathLengthToPrimaryVertex, mBField * kilogauss);
      //SL16d
      double eta, pt, phi;
      if( picoTrack->pMom().perp() > 0 ){
	eta =picoTrack->pMom().pseudoRapidity();
	pt =picoTrack->pMom().perp();
	phi =picoTrack->pMom().phi();
      }
      else{
#ifdef only__primary
	continue;
#endif
	pt = momentum.perp();
	eta = momentum.pseudoRapidity();
	phi = momentum.phi();
      }

      if (mAcceptEvent) hEta->Fill(eta);
      if (mAcceptEvent) hPt->Fill(pt);
      if (fabs(eta) > EventPlaneConstants::etaMaxEventPlane) continue;
      if (pt < EventPlaneConstants::ptMinEventPlane || pt > EventPlaneConstants::ptMaxEventPlane) continue;

      if (mAcceptEvent && eta > 0 && vertexZ > 0) hPhiCentEtaPlusZPlus->Fill(mCent, phi);
      if (mAcceptEvent && eta > 0 && vertexZ < 0) hPhiCentEtaPlusZMinus->Fill(mCent, phi);
      if (mAcceptEvent && eta < 0 && vertexZ > 0) hPhiCentEtaMinusZPlus->Fill(mCent, phi);
      if (mAcceptEvent && eta < 0 && vertexZ < 0) hPhiCentEtaMinusZMinus->Fill(mCent, phi);

      float qx = cos(mHarmonic * phi) * pt;
      float qy = sin(mHarmonic * phi) * pt;

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
	  if(eta > 0.05)
	    {
	      QxEtaPlusGap005 += qx;
	      QyEtaPlusGap005 += qy;
	    }
	}
      else
	{
	  QxEtaMinus += qx;
	  QyEtaMinus += qy;
	  if(eta < -0.05)
	    {
	      QxEtaMinusGap005 += qx;
	      QyEtaMinusGap005 += qy;
	    }
	}

      int iEta = int(eta*10+10);
      QxEta[iEta] += qx;
      QyEta[iEta] += qy;

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
   mQEtaPlusGap005.Set(QxEtaPlusGap005, QyEtaPlusGap005);
   mQEtaMinusGap005.Set(QxEtaMinusGap005,QyEtaMinusGap005);
   for(int i=0; i<20; i++)
     mQEta[i].Set(QxEta[i], QyEta[i]);

   if (mQ.Mod2() == 0 || mQ1.Mod2() == 0 || mQ2.Mod2() == 0 || mQEtaPlus.Mod2() == 0 || mQEtaMinus.Mod2() == 0)
   {
      // cout << "0 Q: " << mQ.Mod2() << " " << mQ1.Mod2() << " " << mQ2.Mod2() << " " << mQEtaPlus.Mod2() << " " << mQEtaMinus.Mod2() << " " << mCent << " " << nTracksForEventPlane << " " << mPicoEvent->refMult() << " " << mPicoEvent->grefMult() << " " << mPicoDst->numberOfTracks() << " " << mPicoEvent->eventId() << " " << mPicoDst->event()->primaryVertex().z() << " " << mPicoDst->event()->vzVpd() << " " << mPicoDst->event()->ZDCx() / 1.0e3 << endl;
      return 1;
   }

   mEventPlane = mQ.Phi()/mHarmonic;
   mEventPlane1 = mQ1.Phi()/mHarmonic;
   mEventPlane2 = mQ2.Phi()/mHarmonic;
   mEventPlaneEtaPlus = mQEtaPlus.Phi()/mHarmonic;
   mEventPlaneEtaMinus = mQEtaMinus.Phi()/mHarmonic;
   mResolutionRandom = cos(mHarmonic*(mEventPlane1 - mEventPlane2));
   mResolutionEta = cos(mHarmonic*(mEventPlaneEtaPlus - mEventPlaneEtaMinus));

   if (mAcceptEvent)
   {
      hQyQxCent->Fill(mCent, Qx, Qy);
      hEventPlaneCent->Fill(mCent, getEventPlane());

      hQyQx1Cent->Fill(mCent, Qx1, Qy1);
      hEventPlane1Cent->Fill(mCent, mEventPlane1);

      hQyQx2Cent->Fill(mCent, Qx2, Qy2);
      hEventPlane2Cent->Fill(mCent, mEventPlane2);

      hQyQxEtaPlusCent->Fill(mCent, QxEtaPlus, QyEtaPlus);
      hEventPlaneEtaPlusCent->Fill(mCent, mEventPlaneEtaPlus);

      hQyQxEtaMinusCent->Fill(mCent, QxEtaMinus, QyEtaMinus);
      hEventPlaneEtaMinusCent->Fill(mCent, mEventPlaneEtaMinus);

      prfCosResolutionRandomCent->Fill(mCent, mResolutionRandom);
      prfCosResolutionEtaCent->Fill(mCent, mResolutionEta);
   }
   return 0;
}

float StEventPlane::getEventPlane(int nTracksToExclude, int* indexTracksToExclude) const
{
   TVector2 Qsub = mQ;
   for (int i = 0; i < nTracksToExclude; i++)
   {
      TVector2 qTrack(qxTracks[indexTracksToExclude[i]], qyTracks[indexTracksToExclude[i]]);
      Qsub -= qTrack;
   }
   return Qsub.Phi()/mHarmonic;
}

TVector2 StEventPlane::QEtaGap(int iEta, int nEtaGaps) const
{
  TVector2 QEtaGap_(0, 0);
  int iEta_ = iEta;
  if(iEta_ < nEtaGaps-1) iEta_ = nEtaGaps-1;
  if(iEta_ > 20-nEtaGaps) iEta_= 20-nEtaGaps;
  for(int i=0; i<20; i++)
    {
      if(fabs(i-iEta_) >= nEtaGaps) QEtaGap_ += mQEta[i];
    }
  return QEtaGap_;
}

void StEventPlane::calculateHadronVn() const
{
   for (unsigned short iTrack = 0; iTrack < mPicoDst->numberOfTracks(); iTrack++)
   {
      StPicoTrack* picoTrack = (StPicoTrack*) mPicoDst->track(iTrack);
      if (!picoTrack)
      {
         break;
      }

      if (picoTrack->nHitsFit() < EventPlaneConstants::nHitsFitMin) continue;

      StPhysicalHelix helix = picoTrack->dcaGeometry().helix();
      float dca = helix.geometricSignedDistance(mVertexPos);
      if (TMath::Abs(dca) > EventPlaneConstants::dcaMaxEventPlane) continue;

      float pathLengthToPrimaryVertex = helix.pathLength(mVertexPos.x(), mVertexPos.y());
      StThreeVectorF momentum = helix.momentumAt(pathLengthToPrimaryVertex, mBField * kilogauss);
      float pt = momentum.perp();
      float eta = momentum.pseudoRapidity();
      float phi = momentum.phi();
      if (fabs(eta) > EventPlaneConstants::etaMaxEventPlane) continue;

      float qx = qxTracks[iTrack];
      float qy = qyTracks[iTrack];
      TVector2 qTrack(qx, qy);
      TVector2 QSub = mQ - qTrack;
      float psi = QSub.Phi()/mHarmonic;

      float weight = mgrefmultCorrUtil->getWeight();

      hHadronVnPtCent->Fill(pt, mCent, cos(mHarmonic*(phi - psi)), weight);

      if (picoTrack->isHFTTrack())
	hHadronHftVnPtCent->Fill(pt, mCent, cos(mHarmonic*(phi - psi)), weight);

      if (picoTrack->pMom().mag() > 0)
	hHadronPrimaryVnPtCent->Fill(pt, mCent, cos(mHarmonic*(phi - psi)), weight);

      if (picoTrack->isHFTTrack() && picoTrack->pMom().mag() > 0)
         hHadronHftPrimaryVnPtCent->Fill(pt, mCent, cos(mHarmonic*(phi - psi)), weight);

      int iEta = (int)(eta*10+10);
      for(int nEtaGaps=0; nEtaGaps<8; nEtaGaps++)
	{
	  TVector2 QSubEtaGap = QEtaGap(iEta, nEtaGaps);
	  int iEta_ = iEta;
	  if(iEta_ < nEtaGaps) iEta_ = nEtaGaps-1;
	  if(iEta_ > 20-nEtaGaps) iEta_= 20-nEtaGaps;
	  if(fabs(iEta-iEta_) >= nEtaGaps)
	    QSubEtaGap -= qTrack;
	  if(QSubEtaGap.Mod()==0) {cout<<"QSubEtaGap.Mod()==0  nEtaGaps: "<<nEtaGaps<<"  cent: "<<mCent<<endl; continue;}
	  float dPhiEtaGap = phi-QSubEtaGap.Phi()/mHarmonic;
	  double toFill[4] = {pt, mCent+0.5, cos(mHarmonic*dPhiEtaGap), 0.1*nEtaGaps+0.05};
	  hHadronVnPtCentEtaGap->Fill(toFill, weight);

	}
   }
}

