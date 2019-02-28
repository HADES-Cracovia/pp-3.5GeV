//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Tue Feb 26 11:21:19 2019 by ROOT version 5.34/34
// from TTree T/T.2
// found on file: /lustre/nyx/hades/user/iciepal/pNb/dst/FILES/lambda1520_100k_01_dst.root
//////////////////////////////////////////////////////////

#ifndef T_h
#define T_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.
#include "hrecevent.h"
#include <TObject.h>
#include "htrack.h"
#include "hpartialevent.h"
#include "hlinearcategory.h"
#include "hcategory.h"
#include "hgeantkine.h"
#include "hmatrixcategory.h"
#include "hlinkeddataobject.h"
#include "hgeantmdc.h"
#include "hgeantshower.h"
#include "hgeanttof.h"
#include "hpidhitdata.h"
#include "hpidtrackdata.h"
#include "hpidtrackcand.h"
#include "hpidgeanttrackset.h"
#include "hpidcandidate.h"
#include <TVector3.h>
#include <TLorentzVector.h>
#include "hpidparticle.h"
#include "heventheader.h"
#include "hgeomvector.h"
#include "heventheader.h"


// Fixed size dimensions of array or collections stored in the TTree if any.
   const Int_t kMaxEvent = 1;
   const Int_t kMaxEvent_fTracks = 1;
   const Int_t kMaxMdc = 1;
   const Int_t kMaxRich = 1;
   const Int_t kMaxShower = 1;
   const Int_t kMaxTof = 1;
   const Int_t kMaxTofino = 1;
   const Int_t kMaxSimul = 1;
   const Int_t kMaxHGeantKine = 1;
   const Int_t kMaxHGeantKine_fData = 98;
   const Int_t kMaxHGeantMdc = 1;
   const Int_t kMaxHGeantMdc_fData = 295;
   const Int_t kMaxHGeantShower = 1;
   const Int_t kMaxHGeantShower_fData = 43;
   const Int_t kMaxHGeantTof = 1;
   const Int_t kMaxHGeantTof_fData = 38;
   const Int_t kMaxTracks = 1;
   const Int_t kMaxPid = 1;
   const Int_t kMaxHPidTrackCandSim = 1;
   const Int_t kMaxHPidTrackCandSim_fData = 48;
   const Int_t kMaxHPidCandidate = 1;
   const Int_t kMaxHPidCandidate_fData = 11;
   const Int_t kMaxHPidParticleSim = 1;
   const Int_t kMaxHPidParticleSim_fData = 11;
   const Int_t kMaxEventHeader = 1;

class T {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types
 //HRecEvent       *Event_;
   UInt_t          Event_HEvent_fUniqueID;
   UInt_t          Event_HEvent_fBits;
   Int_t           Event_fRecLevel;
   Int_t           Event_fNTracks;
   Int_t           Event_fTracks_;
   UInt_t          Event_fTracks_fUniqueID[kMaxEvent_fTracks];   //[Event.fTracks_]
   UInt_t          Event_fTracks_fBits[kMaxEvent_fTracks];   //[Event.fTracks_]
   Float_t         Event_fTracks_fP[kMaxEvent_fTracks];   //[Event.fTracks_]
 //HPartialEvent   *Mdc_;
   UInt_t          Mdc_HEvent_fUniqueID;
   UInt_t          Mdc_HEvent_fBits;
   Int_t           Mdc_fRecLevel;
   Short_t         Mdc_fBaseCategory;
 //HPartialEvent   *Rich_;
   UInt_t          Rich_HEvent_fUniqueID;
   UInt_t          Rich_HEvent_fBits;
   Int_t           Rich_fRecLevel;
   Short_t         Rich_fBaseCategory;
 //HPartialEvent   *Shower_;
   UInt_t          Shower_HEvent_fUniqueID;
   UInt_t          Shower_HEvent_fBits;
   Int_t           Shower_fRecLevel;
   Short_t         Shower_fBaseCategory;
 //HPartialEvent   *Tof_;
   UInt_t          Tof_HEvent_fUniqueID;
   UInt_t          Tof_HEvent_fBits;
   Int_t           Tof_fRecLevel;
   Short_t         Tof_fBaseCategory;
 //HPartialEvent   *Tofino_;
   UInt_t          Tofino_HEvent_fUniqueID;
   UInt_t          Tofino_HEvent_fBits;
   Int_t           Tofino_fRecLevel;
   Short_t         Tofino_fBaseCategory;
 //HPartialEvent   *Simul_;
   UInt_t          Simul_HEvent_fUniqueID;
   UInt_t          Simul_HEvent_fBits;
   Int_t           Simul_fRecLevel;
   Short_t         Simul_fBaseCategory;
 //HLinearCategory *HGeantKine_;
   UInt_t          HGeantKine_HCategory_fUniqueID;
   UInt_t          HGeantKine_HCategory_fBits;
   Short_t         HGeantKine_HCategory_fCat;
   Int_t           HGeantKine_HCategory_fBranchingLevel;
   Int_t           HGeantKine_fData_;
   UInt_t          HGeantKine_fData_fUniqueID[kMaxHGeantKine_fData];   //[HGeantKine.fData_]
   UInt_t          HGeantKine_fData_fBits[kMaxHGeantKine_fData];   //[HGeantKine.fData_]
   Int_t           HGeantKine_fData_trackNumber[kMaxHGeantKine_fData];   //[HGeantKine.fData_]
   Int_t           HGeantKine_fData_parentTrack[kMaxHGeantKine_fData];   //[HGeantKine.fData_]
   Int_t           HGeantKine_fData_particleID[kMaxHGeantKine_fData];   //[HGeantKine.fData_]
   Int_t           HGeantKine_fData_mediumNumber[kMaxHGeantKine_fData];   //[HGeantKine.fData_]
   Int_t           HGeantKine_fData_creationMechanism[kMaxHGeantKine_fData];   //[HGeantKine.fData_]
   Float_t         HGeantKine_fData_xVertex[kMaxHGeantKine_fData];   //[HGeantKine.fData_]
   Float_t         HGeantKine_fData_yVertex[kMaxHGeantKine_fData];   //[HGeantKine.fData_]
   Float_t         HGeantKine_fData_zVertex[kMaxHGeantKine_fData];   //[HGeantKine.fData_]
   Float_t         HGeantKine_fData_xMom[kMaxHGeantKine_fData];   //[HGeantKine.fData_]
   Float_t         HGeantKine_fData_yMom[kMaxHGeantKine_fData];   //[HGeantKine.fData_]
   Float_t         HGeantKine_fData_zMom[kMaxHGeantKine_fData];   //[HGeantKine.fData_]
   Float_t         HGeantKine_fData_generatorInfo[kMaxHGeantKine_fData];   //[HGeantKine.fData_]
   Float_t         HGeantKine_fData_generatorInfo1[kMaxHGeantKine_fData];   //[HGeantKine.fData_]
   Float_t         HGeantKine_fData_generatorInfo2[kMaxHGeantKine_fData];   //[HGeantKine.fData_]
   Float_t         HGeantKine_fData_generatorWeight[kMaxHGeantKine_fData];   //[HGeantKine.fData_]
   Short_t         HGeantKine_fData_firstRichHit[kMaxHGeantKine_fData];   //[HGeantKine.fData_]
   Short_t         HGeantKine_fData_firstMdcHit[kMaxHGeantKine_fData];   //[HGeantKine.fData_]
   Short_t         HGeantKine_fData_firstTofHit[kMaxHGeantKine_fData];   //[HGeantKine.fData_]
   Short_t         HGeantKine_fData_firstRpcHit[kMaxHGeantKine_fData];   //[HGeantKine.fData_]
   Short_t         HGeantKine_fData_firstShowerHit[kMaxHGeantKine_fData];   //[HGeantKine.fData_]
   Short_t         HGeantKine_fData_firstWallHit[kMaxHGeantKine_fData];   //[HGeantKine.fData_]
   Bool_t          HGeantKine_fData_active[kMaxHGeantKine_fData];   //[HGeantKine.fData_]
   Bool_t          HGeantKine_fData_suppressed[kMaxHGeantKine_fData];   //[HGeantKine.fData_]
   Float_t         HGeantKine_fData_userVal[kMaxHGeantKine_fData];   //[HGeantKine.fData_]
   Int_t           HGeantKine_fNDataObjs;
   Bool_t          HGeantKine_hasDynamicObjects;
 //HMatrixCategory *HGeantMdc_;
   UInt_t          HGeantMdc_HCategory_fUniqueID;
   UInt_t          HGeantMdc_HCategory_fBits;
   Short_t         HGeantMdc_HCategory_fCat;
   Int_t           HGeantMdc_HCategory_fBranchingLevel;
   Int_t           HGeantMdc_fNDataObjs;
   Int_t           HGeantMdc_fData_;
   UInt_t          HGeantMdc_fData_fUniqueID[kMaxHGeantMdc_fData];   //[HGeantMdc.fData_]
   UInt_t          HGeantMdc_fData_fBits[kMaxHGeantMdc_fData];   //[HGeantMdc.fData_]
   Short_t         HGeantMdc_fData_nextHit[kMaxHGeantMdc_fData];   //[HGeantMdc.fData_]
   Int_t           HGeantMdc_fData_trackNumber[kMaxHGeantMdc_fData];   //[HGeantMdc.fData_]
   Float_t         HGeantMdc_fData_xHit[kMaxHGeantMdc_fData];   //[HGeantMdc.fData_]
   Float_t         HGeantMdc_fData_yHit[kMaxHGeantMdc_fData];   //[HGeantMdc.fData_]
   Float_t         HGeantMdc_fData_thetaHit[kMaxHGeantMdc_fData];   //[HGeantMdc.fData_]
   Float_t         HGeantMdc_fData_phiHit[kMaxHGeantMdc_fData];   //[HGeantMdc.fData_]
   Float_t         HGeantMdc_fData_tofHit[kMaxHGeantMdc_fData];   //[HGeantMdc.fData_]
   Float_t         HGeantMdc_fData_momHit[kMaxHGeantMdc_fData];   //[HGeantMdc.fData_]
   Char_t          HGeantMdc_fData_sector[kMaxHGeantMdc_fData];   //[HGeantMdc.fData_]
   Char_t          HGeantMdc_fData_module[kMaxHGeantMdc_fData];   //[HGeantMdc.fData_]
   Char_t          HGeantMdc_fData_layer[kMaxHGeantMdc_fData];   //[HGeantMdc.fData_]
 //HMatrixCategory *HGeantShower_;
   UInt_t          HGeantShower_HCategory_fUniqueID;
   UInt_t          HGeantShower_HCategory_fBits;
   Short_t         HGeantShower_HCategory_fCat;
   Int_t           HGeantShower_HCategory_fBranchingLevel;
   Int_t           HGeantShower_fNDataObjs;
   Int_t           HGeantShower_fData_;
   UInt_t          HGeantShower_fData_fUniqueID[kMaxHGeantShower_fData];   //[HGeantShower.fData_]
   UInt_t          HGeantShower_fData_fBits[kMaxHGeantShower_fData];   //[HGeantShower.fData_]
   Short_t         HGeantShower_fData_nextHit[kMaxHGeantShower_fData];   //[HGeantShower.fData_]
   Int_t           HGeantShower_fData_trackNumber[kMaxHGeantShower_fData];   //[HGeantShower.fData_]
   Float_t         HGeantShower_fData_eHit[kMaxHGeantShower_fData];   //[HGeantShower.fData_]
   Float_t         HGeantShower_fData_xHit[kMaxHGeantShower_fData];   //[HGeantShower.fData_]
   Float_t         HGeantShower_fData_yHit[kMaxHGeantShower_fData];   //[HGeantShower.fData_]
   Float_t         HGeantShower_fData_thetaHit[kMaxHGeantShower_fData];   //[HGeantShower.fData_]
   Float_t         HGeantShower_fData_phiHit[kMaxHGeantShower_fData];   //[HGeantShower.fData_]
   Float_t         HGeantShower_fData_betaHit[kMaxHGeantShower_fData];   //[HGeantShower.fData_]
   Char_t          HGeantShower_fData_sector[kMaxHGeantShower_fData];   //[HGeantShower.fData_]
   Char_t          HGeantShower_fData_module[kMaxHGeantShower_fData];   //[HGeantShower.fData_]
 //HMatrixCategory *HGeantTof_;
   UInt_t          HGeantTof_HCategory_fUniqueID;
   UInt_t          HGeantTof_HCategory_fBits;
   Short_t         HGeantTof_HCategory_fCat;
   Int_t           HGeantTof_HCategory_fBranchingLevel;
   Int_t           HGeantTof_fNDataObjs;
   Int_t           HGeantTof_fData_;
   UInt_t          HGeantTof_fData_fUniqueID[kMaxHGeantTof_fData];   //[HGeantTof.fData_]
   UInt_t          HGeantTof_fData_fBits[kMaxHGeantTof_fData];   //[HGeantTof.fData_]
   Short_t         HGeantTof_fData_nextHit[kMaxHGeantTof_fData];   //[HGeantTof.fData_]
   Int_t           HGeantTof_fData_trackNumber[kMaxHGeantTof_fData];   //[HGeantTof.fData_]
   Float_t         HGeantTof_fData_trackLength[kMaxHGeantTof_fData];   //[HGeantTof.fData_]
   Float_t         HGeantTof_fData_eHit[kMaxHGeantTof_fData];   //[HGeantTof.fData_]
   Float_t         HGeantTof_fData_xHit[kMaxHGeantTof_fData];   //[HGeantTof.fData_]
   Float_t         HGeantTof_fData_yHit[kMaxHGeantTof_fData];   //[HGeantTof.fData_]
   Float_t         HGeantTof_fData_tofHit[kMaxHGeantTof_fData];   //[HGeantTof.fData_]
   Float_t         HGeantTof_fData_momHit[kMaxHGeantTof_fData];   //[HGeantTof.fData_]
   Char_t          HGeantTof_fData_sector[kMaxHGeantTof_fData];   //[HGeantTof.fData_]
   Char_t          HGeantTof_fData_module[kMaxHGeantTof_fData];   //[HGeantTof.fData_]
   Char_t          HGeantTof_fData_cell[kMaxHGeantTof_fData];   //[HGeantTof.fData_]
 //HPartialEvent   *Tracks_;
   UInt_t          Tracks_HEvent_fUniqueID;
   UInt_t          Tracks_HEvent_fBits;
   Int_t           Tracks_fRecLevel;
   Short_t         Tracks_fBaseCategory;
 //HPartialEvent   *Pid_;
   UInt_t          Pid_HEvent_fUniqueID;
   UInt_t          Pid_HEvent_fBits;
   Int_t           Pid_fRecLevel;
   Short_t         Pid_fBaseCategory;
 //HLinearCategory *HPidTrackCandSim_;
   UInt_t          HPidTrackCandSim_HCategory_fUniqueID;
   UInt_t          HPidTrackCandSim_HCategory_fBits;
   Short_t         HPidTrackCandSim_HCategory_fCat;
   Int_t           HPidTrackCandSim_HCategory_fBranchingLevel;
   Int_t           HPidTrackCandSim_fData_;
   UInt_t          HPidTrackCandSim_fData_fUniqueID[kMaxHPidTrackCandSim_fData];   //[HPidTrackCandSim.fData_]
   UInt_t          HPidTrackCandSim_fData_fBits[kMaxHPidTrackCandSim_fData];   //[HPidTrackCandSim.fData_]
   UInt_t          HPidTrackCandSim_fData_itsHitData_fUniqueID[kMaxHPidTrackCandSim_fData];   //[HPidTrackCandSim.fData_]
   UInt_t          HPidTrackCandSim_fData_itsHitData_fBits[kMaxHPidTrackCandSim_fData];   //[HPidTrackCandSim.fData_]
   Int_t           HPidTrackCandSim_fData_itsHitData_nSector[kMaxHPidTrackCandSim_fData];   //[HPidTrackCandSim.fData_]
   Short_t         HPidTrackCandSim_fData_itsHitData_iSystem[kMaxHPidTrackCandSim_fData];   //[HPidTrackCandSim.fData_]
   Int_t           HPidTrackCandSim_fData_itsHitData_nRingPadNr[kMaxHPidTrackCandSim_fData];   //[HPidTrackCandSim.fData_]
   Float_t         HPidTrackCandSim_fData_itsHitData_fRingCentroid[kMaxHPidTrackCandSim_fData];   //[HPidTrackCandSim.fData_]
   Float_t         HPidTrackCandSim_fData_itsHitData_fRichTheta[kMaxHPidTrackCandSim_fData];   //[HPidTrackCandSim.fData_]
   Float_t         HPidTrackCandSim_fData_itsHitData_fRichPhi[kMaxHPidTrackCandSim_fData];   //[HPidTrackCandSim.fData_]
   Int_t           HPidTrackCandSim_fData_itsHitData_nRingPatMat[kMaxHPidTrackCandSim_fData];   //[HPidTrackCandSim.fData_]
   Int_t           HPidTrackCandSim_fData_itsHitData_nRingHouTra[kMaxHPidTrackCandSim_fData];   //[HPidTrackCandSim.fData_]
   Int_t           HPidTrackCandSim_fData_itsHitData_nRingAmplitude[kMaxHPidTrackCandSim_fData];   //[HPidTrackCandSim.fData_]
   Int_t           HPidTrackCandSim_fData_itsHitData_nRingLocalMax4[kMaxHPidTrackCandSim_fData];   //[HPidTrackCandSim.fData_]
   Float_t         HPidTrackCandSim_fData_itsHitData_fInnerMdcChiSquare[kMaxHPidTrackCandSim_fData];   //[HPidTrackCandSim.fData_]
   Float_t         HPidTrackCandSim_fData_itsHitData_fInnerMdcdEdx[kMaxHPidTrackCandSim_fData];   //[HPidTrackCandSim.fData_]
   Float_t         HPidTrackCandSim_fData_itsHitData_fInnerMdcdEdxSigma[kMaxHPidTrackCandSim_fData];   //[HPidTrackCandSim.fData_]
   Float_t         HPidTrackCandSim_fData_itsHitData_fMdcRCoord[kMaxHPidTrackCandSim_fData];   //[HPidTrackCandSim.fData_]
   Float_t         HPidTrackCandSim_fData_itsHitData_fMdcZCoord[kMaxHPidTrackCandSim_fData];   //[HPidTrackCandSim.fData_]
   Float_t         HPidTrackCandSim_fData_itsHitData_fMdcTheta[kMaxHPidTrackCandSim_fData];   //[HPidTrackCandSim.fData_]
   Float_t         HPidTrackCandSim_fData_itsHitData_fMdcPhi[kMaxHPidTrackCandSim_fData];   //[HPidTrackCandSim.fData_]
   Float_t         HPidTrackCandSim_fData_itsHitData_fOuterMdcChiSquare[kMaxHPidTrackCandSim_fData];   //[HPidTrackCandSim.fData_]
   Float_t         HPidTrackCandSim_fData_itsHitData_fOuterMdcdEdx[kMaxHPidTrackCandSim_fData];   //[HPidTrackCandSim.fData_]
   Float_t         HPidTrackCandSim_fData_itsHitData_fOuterMdcdEdxSigma[kMaxHPidTrackCandSim_fData];   //[HPidTrackCandSim.fData_]
   Float_t         HPidTrackCandSim_fData_itsHitData_fCombinedMdcdEdx[kMaxHPidTrackCandSim_fData];   //[HPidTrackCandSim.fData_]
   Float_t         HPidTrackCandSim_fData_itsHitData_fCombinedMdcdEdxSigma[kMaxHPidTrackCandSim_fData];   //[HPidTrackCandSim.fData_]
   Int_t           HPidTrackCandSim_fData_itsHitData_iIPURingQuality[kMaxHPidTrackCandSim_fData];   //[HPidTrackCandSim.fData_]
   Int_t           HPidTrackCandSim_fData_itsHitData_iIPUVetoQuality[kMaxHPidTrackCandSim_fData];   //[HPidTrackCandSim.fData_]
   Float_t         HPidTrackCandSim_fData_itsHitData_fShowerSum[kMaxHPidTrackCandSim_fData][3];   //[HPidTrackCandSim.fData_]
   Int_t           HPidTrackCandSim_fData_itsHitData_nShowerClS[kMaxHPidTrackCandSim_fData][3];   //[HPidTrackCandSim.fData_]
   Int_t           HPidTrackCandSim_fData_itsHitData_nShowerRow[kMaxHPidTrackCandSim_fData];   //[HPidTrackCandSim.fData_]
   Int_t           HPidTrackCandSim_fData_itsHitData_nShowerCol[kMaxHPidTrackCandSim_fData];   //[HPidTrackCandSim.fData_]
   Float_t         HPidTrackCandSim_fData_itsHitData_fShowerTimeOfFlight[kMaxHPidTrackCandSim_fData];   //[HPidTrackCandSim.fData_]
   Float_t         HPidTrackCandSim_fData_itsHitData_fMetaLocalX[kMaxHPidTrackCandSim_fData];   //[HPidTrackCandSim.fData_]
   Float_t         HPidTrackCandSim_fData_itsHitData_fMetaLocalY[kMaxHPidTrackCandSim_fData];   //[HPidTrackCandSim.fData_]
   Float_t         HPidTrackCandSim_fData_itsHitData_fTOFTimeOfFlight[kMaxHPidTrackCandSim_fData];   //[HPidTrackCandSim.fData_]
   Float_t         HPidTrackCandSim_fData_itsHitData_fTOFLeftAmplitude[kMaxHPidTrackCandSim_fData];   //[HPidTrackCandSim.fData_]
   Float_t         HPidTrackCandSim_fData_itsHitData_fTOFRightAmplitude[kMaxHPidTrackCandSim_fData];   //[HPidTrackCandSim.fData_]
   Float_t         HPidTrackCandSim_fData_itsHitData_fTofEloss[kMaxHPidTrackCandSim_fData];   //[HPidTrackCandSim.fData_]
   Int_t           HPidTrackCandSim_fData_itsHitData_iTofinoMult[kMaxHPidTrackCandSim_fData];   //[HPidTrackCandSim.fData_]
   Int_t           HPidTrackCandSim_fData_itsHitData_nTofClsSize[kMaxHPidTrackCandSim_fData];   //[HPidTrackCandSim.fData_]
   Int_t           HPidTrackCandSim_fData_itsHitData_nMetaCell[kMaxHPidTrackCandSim_fData];   //[HPidTrackCandSim.fData_]
   Int_t           HPidTrackCandSim_fData_itsHitData_nTofCell[kMaxHPidTrackCandSim_fData];   //[HPidTrackCandSim.fData_]
   Int_t           HPidTrackCandSim_fData_itsHitData_nTofModule[kMaxHPidTrackCandSim_fData];   //[HPidTrackCandSim.fData_]
   Short_t         HPidTrackCandSim_fData_itsHitData_iIndRICH[kMaxHPidTrackCandSim_fData];   //[HPidTrackCandSim.fData_]
   Short_t         HPidTrackCandSim_fData_itsHitData_iIndRICHIPU[kMaxHPidTrackCandSim_fData];   //[HPidTrackCandSim.fData_]
   Short_t         HPidTrackCandSim_fData_itsHitData_iIndInnerSeg[kMaxHPidTrackCandSim_fData];   //[HPidTrackCandSim.fData_]
   Short_t         HPidTrackCandSim_fData_itsHitData_iIndOuterSeg[kMaxHPidTrackCandSim_fData];   //[HPidTrackCandSim.fData_]
   Short_t         HPidTrackCandSim_fData_itsHitData_iIndTOF[kMaxHPidTrackCandSim_fData];   //[HPidTrackCandSim.fData_]
   Short_t         HPidTrackCandSim_fData_itsHitData_iIndShower[kMaxHPidTrackCandSim_fData];   //[HPidTrackCandSim.fData_]
   Short_t         HPidTrackCandSim_fData_itsHitData_iIndClusInf0[kMaxHPidTrackCandSim_fData];   //[HPidTrackCandSim.fData_]
   Short_t         HPidTrackCandSim_fData_itsHitData_iIndClusInf1[kMaxHPidTrackCandSim_fData];   //[HPidTrackCandSim.fData_]
   Short_t         HPidTrackCandSim_fData_itsHitData_iIndClusInf2[kMaxHPidTrackCandSim_fData];   //[HPidTrackCandSim.fData_]
   Short_t         HPidTrackCandSim_fData_itsHitData_iIndClusInf3[kMaxHPidTrackCandSim_fData];   //[HPidTrackCandSim.fData_]
   Short_t         HPidTrackCandSim_fData_itsHitData_iIndMatch[kMaxHPidTrackCandSim_fData];   //[HPidTrackCandSim.fData_]
   Float_t         HPidTrackCandSim_fData_itsHitData_fDistanceToVertex[kMaxHPidTrackCandSim_fData][10];   //[HPidTrackCandSim.fData_]
   Bool_t          HPidTrackCandSim_fData_itsHitData_hasRingCorrelation[kMaxHPidTrackCandSim_fData][10];   //[HPidTrackCandSim.fData_]
   Bool_t          HPidTrackCandSim_fData_itsHitData_hasMetaTrackCorrelation[kMaxHPidTrackCandSim_fData][10];   //[HPidTrackCandSim.fData_]
   UInt_t          HPidTrackCandSim_fData_itsTrackData_fUniqueID[kMaxHPidTrackCandSim_fData];   //[HPidTrackCandSim.fData_]
   UInt_t          HPidTrackCandSim_fData_itsTrackData_fBits[kMaxHPidTrackCandSim_fData];   //[HPidTrackCandSim.fData_]
   HSymMat6        HPidTrackCandSim_fData_itsTrackData_cov[10][kMaxHPidTrackCandSim_fData];
   Int_t           HPidTrackCandSim_fData_itsTrackData_nBestMomAlg[kMaxHPidTrackCandSim_fData];   //[HPidTrackCandSim.fData_]
   Short_t         HPidTrackCandSim_fData_itsTrackData_nRKTrackInd[kMaxHPidTrackCandSim_fData];   //[HPidTrackCandSim.fData_]
   Short_t         HPidTrackCandSim_fData_itsTrackData_nKickTrackInd[kMaxHPidTrackCandSim_fData];   //[HPidTrackCandSim.fData_]
   Short_t         HPidTrackCandSim_fData_itsTrackData_nKickTrack123Ind[kMaxHPidTrackCandSim_fData];   //[HPidTrackCandSim.fData_]
   Short_t         HPidTrackCandSim_fData_itsTrackData_nRefTrackInd[kMaxHPidTrackCandSim_fData];   //[HPidTrackCandSim.fData_]
   Short_t         HPidTrackCandSim_fData_itsTrackData_nSplineTrackInd[kMaxHPidTrackCandSim_fData];   //[HPidTrackCandSim.fData_]
   Float_t         HPidTrackCandSim_fData_itsTrackData_fMetaMatchingQuality[kMaxHPidTrackCandSim_fData];   //[HPidTrackCandSim.fData_]
   Float_t         HPidTrackCandSim_fData_itsTrackData_fRKRichMatchingQuality[kMaxHPidTrackCandSim_fData];   //[HPidTrackCandSim.fData_]
   Float_t         HPidTrackCandSim_fData_itsTrackData_dxRkMeta[kMaxHPidTrackCandSim_fData];   //[HPidTrackCandSim.fData_]
   Float_t         HPidTrackCandSim_fData_itsTrackData_dyRkMeta[kMaxHPidTrackCandSim_fData];   //[HPidTrackCandSim.fData_]
   Float_t         HPidTrackCandSim_fData_itsTrackData_dzRkMeta[kMaxHPidTrackCandSim_fData];   //[HPidTrackCandSim.fData_]
   Float_t         HPidTrackCandSim_fData_itsTrackData_dxMdcMeta[kMaxHPidTrackCandSim_fData];   //[HPidTrackCandSim.fData_]
   Float_t         HPidTrackCandSim_fData_itsTrackData_dyMdcMeta[kMaxHPidTrackCandSim_fData];   //[HPidTrackCandSim.fData_]
   Float_t         HPidTrackCandSim_fData_itsTrackData_dzMdcMeta[kMaxHPidTrackCandSim_fData];   //[HPidTrackCandSim.fData_]
   Float_t         HPidTrackCandSim_fData_itsTrackData_xMeta[kMaxHPidTrackCandSim_fData];   //[HPidTrackCandSim.fData_]
   Float_t         HPidTrackCandSim_fData_itsTrackData_yMeta[kMaxHPidTrackCandSim_fData];   //[HPidTrackCandSim.fData_]
   Float_t         HPidTrackCandSim_fData_itsTrackData_zMeta[kMaxHPidTrackCandSim_fData];   //[HPidTrackCandSim.fData_]
   Float_t         HPidTrackCandSim_fData_itsTrackData_errXMeta[kMaxHPidTrackCandSim_fData];   //[HPidTrackCandSim.fData_]
   Float_t         HPidTrackCandSim_fData_itsTrackData_errYMeta[kMaxHPidTrackCandSim_fData];   //[HPidTrackCandSim.fData_]
   Float_t         HPidTrackCandSim_fData_itsTrackData_errZMeta[kMaxHPidTrackCandSim_fData];   //[HPidTrackCandSim.fData_]
   Int_t           HPidTrackCandSim_fData_itsTrackData_nCloseTracklets[kMaxHPidTrackCandSim_fData];   //[HPidTrackCandSim.fData_]
   Float_t         HPidTrackCandSim_fData_itsTrackData_fPull[kMaxHPidTrackCandSim_fData];   //[HPidTrackCandSim.fData_]
   Float_t         HPidTrackCandSim_fData_itsTrackData_fSplineChiSquare[kMaxHPidTrackCandSim_fData];   //[HPidTrackCandSim.fData_]
   Float_t         HPidTrackCandSim_fData_itsTrackData_fRKChiSquare[kMaxHPidTrackCandSim_fData];   //[HPidTrackCandSim.fData_]
   Int_t           HPidTrackCandSim_fData_itsTrackData_iIndexClosestTracklet[kMaxHPidTrackCandSim_fData];   //[HPidTrackCandSim.fData_]
   TArrayI         HPidTrackCandSim_fData_itsTrackData_aTrackletClusInf0[kMaxHPidTrackCandSim_fData];
   TArrayI         HPidTrackCandSim_fData_itsTrackData_aTrackletClusInf1[kMaxHPidTrackCandSim_fData];
   TArrayF         HPidTrackCandSim_fData_itsTrackData_aTrackletDistances[kMaxHPidTrackCandSim_fData];
   Float_t         HPidTrackCandSim_fData_itsTrackData_qIOMatching[kMaxHPidTrackCandSim_fData][10];   //[HPidTrackCandSim.fData_]
   Short_t         HPidTrackCandSim_fData_itsTrackData_nPolarity[kMaxHPidTrackCandSim_fData][10];   //[HPidTrackCandSim.fData_]
   Float_t         HPidTrackCandSim_fData_itsTrackData_fMomenta[kMaxHPidTrackCandSim_fData][10];   //[HPidTrackCandSim.fData_]
   Float_t         HPidTrackCandSim_fData_itsTrackData_fMomError[kMaxHPidTrackCandSim_fData][10];   //[HPidTrackCandSim.fData_]
   Float_t         HPidTrackCandSim_fData_itsTrackData_fTrackR[kMaxHPidTrackCandSim_fData][10];   //[HPidTrackCandSim.fData_]
   Float_t         HPidTrackCandSim_fData_itsTrackData_fTrackZ[kMaxHPidTrackCandSim_fData][10];   //[HPidTrackCandSim.fData_]
   Float_t         HPidTrackCandSim_fData_itsTrackData_fRKPhi[kMaxHPidTrackCandSim_fData];   //[HPidTrackCandSim.fData_]
   Float_t         HPidTrackCandSim_fData_itsTrackData_fRKTheta[kMaxHPidTrackCandSim_fData];   //[HPidTrackCandSim.fData_]
   Float_t         HPidTrackCandSim_fData_itsTrackData_fCorrectedEloss[kMaxHPidTrackCandSim_fData][10];   //[HPidTrackCandSim.fData_]
   Float_t         HPidTrackCandSim_fData_itsTrackData_fCorrectedBeta[kMaxHPidTrackCandSim_fData][10];   //[HPidTrackCandSim.fData_]
   Float_t         HPidTrackCandSim_fData_itsTrackData_fPathLength[kMaxHPidTrackCandSim_fData][10];   //[HPidTrackCandSim.fData_]
   Float_t         HPidTrackCandSim_fData_itsTrackData_fMassSquared[kMaxHPidTrackCandSim_fData][10];   //[HPidTrackCandSim.fData_]
   Bool_t          HPidTrackCandSim_fData_itsTrackData_bIsAccepted[kMaxHPidTrackCandSim_fData][10];   //[HPidTrackCandSim.fData_]
   Short_t         HPidTrackCandSim_fData_itsTrackData_nTofRecFlag[kMaxHPidTrackCandSim_fData][10];   //[HPidTrackCandSim.fData_]
   Float_t         HPidTrackCandSim_fData_itsTrackData_fTofRecTof[kMaxHPidTrackCandSim_fData][10];   //[HPidTrackCandSim.fData_]
   Float_t         HPidTrackCandSim_fData_itsTrackData_fTofRecBeta[kMaxHPidTrackCandSim_fData][10];   //[HPidTrackCandSim.fData_]
   Float_t         HPidTrackCandSim_fData_itsTrackData_fTofRecMassSquared[kMaxHPidTrackCandSim_fData][10];   //[HPidTrackCandSim.fData_]
   Int_t           HPidTrackCandSim_fData_flags[kMaxHPidTrackCandSim_fData];   //[HPidTrackCandSim.fData_]
   UInt_t          HPidTrackCandSim_fData_itsGeantTrackSet_fUniqueID[kMaxHPidTrackCandSim_fData];   //[HPidTrackCandSim.fData_]
   UInt_t          HPidTrackCandSim_fData_itsGeantTrackSet_fBits[kMaxHPidTrackCandSim_fData];   //[HPidTrackCandSim.fData_]
   Bool_t          HPidTrackCandSim_fData_itsGeantTrackSet_isSorted[kMaxHPidTrackCandSim_fData];   //[HPidTrackCandSim.fData_]
   Short_t         HPidTrackCandSim_fData_itsGeantTrackSet_sNCorrTrackIds[kMaxHPidTrackCandSim_fData];   //[HPidTrackCandSim.fData_]
   UInt_t          HPidTrackCandSim_fData_itsGeantTrackSet_nRICHTracks[kMaxHPidTrackCandSim_fData];   //[HPidTrackCandSim.fData_]
   UInt_t          HPidTrackCandSim_fData_itsGeantTrackSet_nRICHIPUTracks[kMaxHPidTrackCandSim_fData];   //[HPidTrackCandSim.fData_]
   UInt_t          HPidTrackCandSim_fData_itsGeantTrackSet_nInnerMdcTracks[kMaxHPidTrackCandSim_fData];   //[HPidTrackCandSim.fData_]
   UInt_t          HPidTrackCandSim_fData_itsGeantTrackSet_nOuterMdcTracks[kMaxHPidTrackCandSim_fData];   //[HPidTrackCandSim.fData_]
   UInt_t          HPidTrackCandSim_fData_itsGeantTrackSet_nShowerTracks[kMaxHPidTrackCandSim_fData];   //[HPidTrackCandSim.fData_]
   UInt_t          HPidTrackCandSim_fData_itsGeantTrackSet_nTOFTracks[kMaxHPidTrackCandSim_fData];   //[HPidTrackCandSim.fData_]
   Bool_t          HPidTrackCandSim_fData_itsGeantTrackSet_bIsLepFromPrimary[kMaxHPidTrackCandSim_fData];   //[HPidTrackCandSim.fData_]
   TArrayI         HPidTrackCandSim_fData_itsGeantTrackSet_correlatedTrackIds[kMaxHPidTrackCandSim_fData];
   TArrayI         HPidTrackCandSim_fData_itsGeantTrackSet_correlationFlags[kMaxHPidTrackCandSim_fData];
   TArrayI         HPidTrackCandSim_fData_itsGeantTrackSet_ProcessIds[kMaxHPidTrackCandSim_fData];
   TArrayI         HPidTrackCandSim_fData_itsGeantTrackSet_ParentIds[kMaxHPidTrackCandSim_fData];
   TArrayI         HPidTrackCandSim_fData_itsGeantTrackSet_Parents[kMaxHPidTrackCandSim_fData];
   TArrayI         HPidTrackCandSim_fData_itsGeantTrackSet_GrandParentIds[kMaxHPidTrackCandSim_fData];
   TArrayI         HPidTrackCandSim_fData_itsGeantTrackSet_GrandParents[kMaxHPidTrackCandSim_fData];
   TArrayF         HPidTrackCandSim_fData_itsGeantTrackSet_GenInfo[kMaxHPidTrackCandSim_fData];
   TArrayF         HPidTrackCandSim_fData_itsGeantTrackSet_GenInfo1[kMaxHPidTrackCandSim_fData];
   TArrayF         HPidTrackCandSim_fData_itsGeantTrackSet_GenInfo2[kMaxHPidTrackCandSim_fData];
   TArrayF         HPidTrackCandSim_fData_itsGeantTrackSet_GenWeight[kMaxHPidTrackCandSim_fData];
   TArrayF         HPidTrackCandSim_fData_itsGeantTrackSet_VertexX[kMaxHPidTrackCandSim_fData];
   TArrayF         HPidTrackCandSim_fData_itsGeantTrackSet_VertexY[kMaxHPidTrackCandSim_fData];
   TArrayF         HPidTrackCandSim_fData_itsGeantTrackSet_VertexZ[kMaxHPidTrackCandSim_fData];
   TArrayI         HPidTrackCandSim_fData_itsGeantTrackSet_GeantPIDs[kMaxHPidTrackCandSim_fData];
   TArrayI         HPidTrackCandSim_fData_itsGeantTrackSet_MediumIds[kMaxHPidTrackCandSim_fData];
   TArrayF         HPidTrackCandSim_fData_itsGeantTrackSet_GeantMomX[kMaxHPidTrackCandSim_fData];
   TArrayF         HPidTrackCandSim_fData_itsGeantTrackSet_GeantMomY[kMaxHPidTrackCandSim_fData];
   TArrayF         HPidTrackCandSim_fData_itsGeantTrackSet_GeantMomZ[kMaxHPidTrackCandSim_fData];
   TArrayF         HPidTrackCandSim_fData_itsGeantTrackSet_ShowerWeights[kMaxHPidTrackCandSim_fData];
   TArrayF         HPidTrackCandSim_fData_itsGeantTrackSet_TOFWeights[kMaxHPidTrackCandSim_fData];
   TArrayF         HPidTrackCandSim_fData_itsGeantTrackSet_RICHWeights[kMaxHPidTrackCandSim_fData];
   TArrayF         HPidTrackCandSim_fData_itsGeantTrackSet_RICHIPUWeights[kMaxHPidTrackCandSim_fData];
   TArrayF         HPidTrackCandSim_fData_itsGeantTrackSet_InnerMDCWeights[kMaxHPidTrackCandSim_fData];
   TArrayF         HPidTrackCandSim_fData_itsGeantTrackSet_OuterMDCWeights[kMaxHPidTrackCandSim_fData];
   TArrayI         HPidTrackCandSim_fData_itsGeantTrackSet_hasHitInShower[kMaxHPidTrackCandSim_fData];
   TArrayI         HPidTrackCandSim_fData_itsGeantTrackSet_hasHitInTOF[kMaxHPidTrackCandSim_fData];
   Int_t           HPidTrackCandSim_fNDataObjs;
   Bool_t          HPidTrackCandSim_hasDynamicObjects;
 //HLinearCategory *HPidCandidate_;
   UInt_t          HPidCandidate_HCategory_fUniqueID;
   UInt_t          HPidCandidate_HCategory_fBits;
   Short_t         HPidCandidate_HCategory_fCat;
   Int_t           HPidCandidate_HCategory_fBranchingLevel;
   Int_t           HPidCandidate_fData_;
   UInt_t          HPidCandidate_fData_fUniqueID[kMaxHPidCandidate_fData];   //[HPidCandidate.fData_]
   UInt_t          HPidCandidate_fData_fBits[kMaxHPidCandidate_fData];   //[HPidCandidate.fData_]
   Short_t         HPidCandidate_fData_iTrackCandIndex[kMaxHPidCandidate_fData];   //[HPidCandidate.fData_]
   UInt_t          HPidCandidate_fData_NUM_ALGORITHMS[kMaxHPidCandidate_fData];   //[HPidCandidate.fData_]
   UInt_t          HPidCandidate_fData_NUM_PARTICLES[kMaxHPidCandidate_fData];   //[HPidCandidate.fData_]
   UInt_t          HPidCandidate_fData_NUM_VALUES[kMaxHPidCandidate_fData];   //[HPidCandidate.fData_]
   TArrayS         HPidCandidate_fData_aAlgorithms[kMaxHPidCandidate_fData];
   TArrayS         HPidCandidate_fData_aParticles[kMaxHPidCandidate_fData];
   TArrayF         HPidCandidate_fData_aValues[kMaxHPidCandidate_fData];
   Int_t           HPidCandidate_fData_nMomAlgIndex[kMaxHPidCandidate_fData];   //[HPidCandidate.fData_]
   Int_t           HPidCandidate_fData_flags[kMaxHPidCandidate_fData];   //[HPidCandidate.fData_]
   Int_t           HPidCandidate_fNDataObjs;
   Bool_t          HPidCandidate_hasDynamicObjects;
 //HLinearCategory *HPidParticleSim_;
   UInt_t          HPidParticleSim_HCategory_fUniqueID;
   UInt_t          HPidParticleSim_HCategory_fBits;
   Short_t         HPidParticleSim_HCategory_fCat;
   Int_t           HPidParticleSim_HCategory_fBranchingLevel;
   Int_t           HPidParticleSim_fData_;
   UInt_t          HPidParticleSim_fData_fUniqueID[kMaxHPidParticleSim_fData];   //[HPidParticleSim.fData_]
   UInt_t          HPidParticleSim_fData_fBits[kMaxHPidParticleSim_fData];   //[HPidParticleSim.fData_]
   UInt_t          HPidParticleSim_fData_fP_fUniqueID[kMaxHPidParticleSim_fData];   //[HPidParticleSim.fData_]
   UInt_t          HPidParticleSim_fData_fP_fBits[kMaxHPidParticleSim_fData];   //[HPidParticleSim.fData_]
   Double_t        HPidParticleSim_fData_fP_fX[kMaxHPidParticleSim_fData];   //[HPidParticleSim.fData_]
   Double_t        HPidParticleSim_fData_fP_fY[kMaxHPidParticleSim_fData];   //[HPidParticleSim.fData_]
   Double_t        HPidParticleSim_fData_fP_fZ[kMaxHPidParticleSim_fData];   //[HPidParticleSim.fData_]
   Double_t        HPidParticleSim_fData_fE[kMaxHPidParticleSim_fData];   //[HPidParticleSim.fData_]
   Bool_t          HPidParticleSim_fData_kUsesIdealMass[kMaxHPidParticleSim_fData];   //[HPidParticleSim.fData_]
   Int_t           HPidParticleSim_fData_nPossibleSpecies[kMaxHPidParticleSim_fData];   //[HPidParticleSim.fData_]
   Short_t         HPidParticleSim_fData_momAlgIndex[kMaxHPidParticleSim_fData];   //[HPidParticleSim.fData_]
   Short_t         HPidParticleSim_fData_nPidCandidateIndex[kMaxHPidParticleSim_fData];   //[HPidParticleSim.fData_]
   TArrayS         HPidParticleSim_fData_possibleSpecies[kMaxHPidParticleSim_fData];
   TArrayF         HPidParticleSim_fData_assignedWeights[kMaxHPidParticleSim_fData];
   Int_t           HPidParticleSim_fData_nAssignedPID[kMaxHPidParticleSim_fData];   //[HPidParticleSim.fData_]
   Float_t         HPidParticleSim_fData_fTestVal[kMaxHPidParticleSim_fData];   //[HPidParticleSim.fData_]
   Float_t         HPidParticleSim_fData_fWeight[kMaxHPidParticleSim_fData];   //[HPidParticleSim.fData_]
   Float_t         HPidParticleSim_fData_fMomRescal[kMaxHPidParticleSim_fData];   //[HPidParticleSim.fData_]
   UInt_t          HPidParticleSim_fData_itsHitData_fUniqueID[kMaxHPidParticleSim_fData];   //[HPidParticleSim.fData_]
   UInt_t          HPidParticleSim_fData_itsHitData_fBits[kMaxHPidParticleSim_fData];   //[HPidParticleSim.fData_]
   Int_t           HPidParticleSim_fData_itsHitData_nSector[kMaxHPidParticleSim_fData];   //[HPidParticleSim.fData_]
   Short_t         HPidParticleSim_fData_itsHitData_iSystem[kMaxHPidParticleSim_fData];   //[HPidParticleSim.fData_]
   Int_t           HPidParticleSim_fData_itsHitData_nRingPadNr[kMaxHPidParticleSim_fData];   //[HPidParticleSim.fData_]
   Float_t         HPidParticleSim_fData_itsHitData_fRingCentroid[kMaxHPidParticleSim_fData];   //[HPidParticleSim.fData_]
   Float_t         HPidParticleSim_fData_itsHitData_fRichTheta[kMaxHPidParticleSim_fData];   //[HPidParticleSim.fData_]
   Float_t         HPidParticleSim_fData_itsHitData_fRichPhi[kMaxHPidParticleSim_fData];   //[HPidParticleSim.fData_]
   Int_t           HPidParticleSim_fData_itsHitData_nRingPatMat[kMaxHPidParticleSim_fData];   //[HPidParticleSim.fData_]
   Int_t           HPidParticleSim_fData_itsHitData_nRingHouTra[kMaxHPidParticleSim_fData];   //[HPidParticleSim.fData_]
   Int_t           HPidParticleSim_fData_itsHitData_nRingAmplitude[kMaxHPidParticleSim_fData];   //[HPidParticleSim.fData_]
   Int_t           HPidParticleSim_fData_itsHitData_nRingLocalMax4[kMaxHPidParticleSim_fData];   //[HPidParticleSim.fData_]
   Float_t         HPidParticleSim_fData_itsHitData_fInnerMdcChiSquare[kMaxHPidParticleSim_fData];   //[HPidParticleSim.fData_]
   Float_t         HPidParticleSim_fData_itsHitData_fInnerMdcdEdx[kMaxHPidParticleSim_fData];   //[HPidParticleSim.fData_]
   Float_t         HPidParticleSim_fData_itsHitData_fInnerMdcdEdxSigma[kMaxHPidParticleSim_fData];   //[HPidParticleSim.fData_]
   Float_t         HPidParticleSim_fData_itsHitData_fMdcRCoord[kMaxHPidParticleSim_fData];   //[HPidParticleSim.fData_]
   Float_t         HPidParticleSim_fData_itsHitData_fMdcZCoord[kMaxHPidParticleSim_fData];   //[HPidParticleSim.fData_]
   Float_t         HPidParticleSim_fData_itsHitData_fMdcTheta[kMaxHPidParticleSim_fData];   //[HPidParticleSim.fData_]
   Float_t         HPidParticleSim_fData_itsHitData_fMdcPhi[kMaxHPidParticleSim_fData];   //[HPidParticleSim.fData_]
   Float_t         HPidParticleSim_fData_itsHitData_fOuterMdcChiSquare[kMaxHPidParticleSim_fData];   //[HPidParticleSim.fData_]
   Float_t         HPidParticleSim_fData_itsHitData_fOuterMdcdEdx[kMaxHPidParticleSim_fData];   //[HPidParticleSim.fData_]
   Float_t         HPidParticleSim_fData_itsHitData_fOuterMdcdEdxSigma[kMaxHPidParticleSim_fData];   //[HPidParticleSim.fData_]
   Float_t         HPidParticleSim_fData_itsHitData_fCombinedMdcdEdx[kMaxHPidParticleSim_fData];   //[HPidParticleSim.fData_]
   Float_t         HPidParticleSim_fData_itsHitData_fCombinedMdcdEdxSigma[kMaxHPidParticleSim_fData];   //[HPidParticleSim.fData_]
   Int_t           HPidParticleSim_fData_itsHitData_iIPURingQuality[kMaxHPidParticleSim_fData];   //[HPidParticleSim.fData_]
   Int_t           HPidParticleSim_fData_itsHitData_iIPUVetoQuality[kMaxHPidParticleSim_fData];   //[HPidParticleSim.fData_]
   Float_t         HPidParticleSim_fData_itsHitData_fShowerSum[kMaxHPidParticleSim_fData][3];   //[HPidParticleSim.fData_]
   Int_t           HPidParticleSim_fData_itsHitData_nShowerClS[kMaxHPidParticleSim_fData][3];   //[HPidParticleSim.fData_]
   Int_t           HPidParticleSim_fData_itsHitData_nShowerRow[kMaxHPidParticleSim_fData];   //[HPidParticleSim.fData_]
   Int_t           HPidParticleSim_fData_itsHitData_nShowerCol[kMaxHPidParticleSim_fData];   //[HPidParticleSim.fData_]
   Float_t         HPidParticleSim_fData_itsHitData_fShowerTimeOfFlight[kMaxHPidParticleSim_fData];   //[HPidParticleSim.fData_]
   Float_t         HPidParticleSim_fData_itsHitData_fMetaLocalX[kMaxHPidParticleSim_fData];   //[HPidParticleSim.fData_]
   Float_t         HPidParticleSim_fData_itsHitData_fMetaLocalY[kMaxHPidParticleSim_fData];   //[HPidParticleSim.fData_]
   Float_t         HPidParticleSim_fData_itsHitData_fTOFTimeOfFlight[kMaxHPidParticleSim_fData];   //[HPidParticleSim.fData_]
   Float_t         HPidParticleSim_fData_itsHitData_fTOFLeftAmplitude[kMaxHPidParticleSim_fData];   //[HPidParticleSim.fData_]
   Float_t         HPidParticleSim_fData_itsHitData_fTOFRightAmplitude[kMaxHPidParticleSim_fData];   //[HPidParticleSim.fData_]
   Float_t         HPidParticleSim_fData_itsHitData_fTofEloss[kMaxHPidParticleSim_fData];   //[HPidParticleSim.fData_]
   Int_t           HPidParticleSim_fData_itsHitData_iTofinoMult[kMaxHPidParticleSim_fData];   //[HPidParticleSim.fData_]
   Int_t           HPidParticleSim_fData_itsHitData_nTofClsSize[kMaxHPidParticleSim_fData];   //[HPidParticleSim.fData_]
   Int_t           HPidParticleSim_fData_itsHitData_nMetaCell[kMaxHPidParticleSim_fData];   //[HPidParticleSim.fData_]
   Int_t           HPidParticleSim_fData_itsHitData_nTofCell[kMaxHPidParticleSim_fData];   //[HPidParticleSim.fData_]
   Int_t           HPidParticleSim_fData_itsHitData_nTofModule[kMaxHPidParticleSim_fData];   //[HPidParticleSim.fData_]
   Short_t         HPidParticleSim_fData_itsHitData_iIndRICH[kMaxHPidParticleSim_fData];   //[HPidParticleSim.fData_]
   Short_t         HPidParticleSim_fData_itsHitData_iIndRICHIPU[kMaxHPidParticleSim_fData];   //[HPidParticleSim.fData_]
   Short_t         HPidParticleSim_fData_itsHitData_iIndInnerSeg[kMaxHPidParticleSim_fData];   //[HPidParticleSim.fData_]
   Short_t         HPidParticleSim_fData_itsHitData_iIndOuterSeg[kMaxHPidParticleSim_fData];   //[HPidParticleSim.fData_]
   Short_t         HPidParticleSim_fData_itsHitData_iIndTOF[kMaxHPidParticleSim_fData];   //[HPidParticleSim.fData_]
   Short_t         HPidParticleSim_fData_itsHitData_iIndShower[kMaxHPidParticleSim_fData];   //[HPidParticleSim.fData_]
   Short_t         HPidParticleSim_fData_itsHitData_iIndClusInf0[kMaxHPidParticleSim_fData];   //[HPidParticleSim.fData_]
   Short_t         HPidParticleSim_fData_itsHitData_iIndClusInf1[kMaxHPidParticleSim_fData];   //[HPidParticleSim.fData_]
   Short_t         HPidParticleSim_fData_itsHitData_iIndClusInf2[kMaxHPidParticleSim_fData];   //[HPidParticleSim.fData_]
   Short_t         HPidParticleSim_fData_itsHitData_iIndClusInf3[kMaxHPidParticleSim_fData];   //[HPidParticleSim.fData_]
   Short_t         HPidParticleSim_fData_itsHitData_iIndMatch[kMaxHPidParticleSim_fData];   //[HPidParticleSim.fData_]
   Float_t         HPidParticleSim_fData_itsHitData_fDistanceToVertex[kMaxHPidParticleSim_fData][10];   //[HPidParticleSim.fData_]
   Bool_t          HPidParticleSim_fData_itsHitData_hasRingCorrelation[kMaxHPidParticleSim_fData][10];   //[HPidParticleSim.fData_]
   Bool_t          HPidParticleSim_fData_itsHitData_hasMetaTrackCorrelation[kMaxHPidParticleSim_fData][10];   //[HPidParticleSim.fData_]
   UInt_t          HPidParticleSim_fData_itsTrackData_fUniqueID[kMaxHPidParticleSim_fData];   //[HPidParticleSim.fData_]
   UInt_t          HPidParticleSim_fData_itsTrackData_fBits[kMaxHPidParticleSim_fData];   //[HPidParticleSim.fData_]
   HSymMat6        HPidParticleSim_fData_itsTrackData_cov[10][kMaxHPidParticleSim_fData];
   Int_t           HPidParticleSim_fData_itsTrackData_nBestMomAlg[kMaxHPidParticleSim_fData];   //[HPidParticleSim.fData_]
   Short_t         HPidParticleSim_fData_itsTrackData_nRKTrackInd[kMaxHPidParticleSim_fData];   //[HPidParticleSim.fData_]
   Short_t         HPidParticleSim_fData_itsTrackData_nKickTrackInd[kMaxHPidParticleSim_fData];   //[HPidParticleSim.fData_]
   Short_t         HPidParticleSim_fData_itsTrackData_nKickTrack123Ind[kMaxHPidParticleSim_fData];   //[HPidParticleSim.fData_]
   Short_t         HPidParticleSim_fData_itsTrackData_nRefTrackInd[kMaxHPidParticleSim_fData];   //[HPidParticleSim.fData_]
   Short_t         HPidParticleSim_fData_itsTrackData_nSplineTrackInd[kMaxHPidParticleSim_fData];   //[HPidParticleSim.fData_]
   Float_t         HPidParticleSim_fData_itsTrackData_fMetaMatchingQuality[kMaxHPidParticleSim_fData];   //[HPidParticleSim.fData_]
   Float_t         HPidParticleSim_fData_itsTrackData_fRKRichMatchingQuality[kMaxHPidParticleSim_fData];   //[HPidParticleSim.fData_]
   Float_t         HPidParticleSim_fData_itsTrackData_dxRkMeta[kMaxHPidParticleSim_fData];   //[HPidParticleSim.fData_]
   Float_t         HPidParticleSim_fData_itsTrackData_dyRkMeta[kMaxHPidParticleSim_fData];   //[HPidParticleSim.fData_]
   Float_t         HPidParticleSim_fData_itsTrackData_dzRkMeta[kMaxHPidParticleSim_fData];   //[HPidParticleSim.fData_]
   Float_t         HPidParticleSim_fData_itsTrackData_dxMdcMeta[kMaxHPidParticleSim_fData];   //[HPidParticleSim.fData_]
   Float_t         HPidParticleSim_fData_itsTrackData_dyMdcMeta[kMaxHPidParticleSim_fData];   //[HPidParticleSim.fData_]
   Float_t         HPidParticleSim_fData_itsTrackData_dzMdcMeta[kMaxHPidParticleSim_fData];   //[HPidParticleSim.fData_]
   Float_t         HPidParticleSim_fData_itsTrackData_xMeta[kMaxHPidParticleSim_fData];   //[HPidParticleSim.fData_]
   Float_t         HPidParticleSim_fData_itsTrackData_yMeta[kMaxHPidParticleSim_fData];   //[HPidParticleSim.fData_]
   Float_t         HPidParticleSim_fData_itsTrackData_zMeta[kMaxHPidParticleSim_fData];   //[HPidParticleSim.fData_]
   Float_t         HPidParticleSim_fData_itsTrackData_errXMeta[kMaxHPidParticleSim_fData];   //[HPidParticleSim.fData_]
   Float_t         HPidParticleSim_fData_itsTrackData_errYMeta[kMaxHPidParticleSim_fData];   //[HPidParticleSim.fData_]
   Float_t         HPidParticleSim_fData_itsTrackData_errZMeta[kMaxHPidParticleSim_fData];   //[HPidParticleSim.fData_]
   Int_t           HPidParticleSim_fData_itsTrackData_nCloseTracklets[kMaxHPidParticleSim_fData];   //[HPidParticleSim.fData_]
   Float_t         HPidParticleSim_fData_itsTrackData_fPull[kMaxHPidParticleSim_fData];   //[HPidParticleSim.fData_]
   Float_t         HPidParticleSim_fData_itsTrackData_fSplineChiSquare[kMaxHPidParticleSim_fData];   //[HPidParticleSim.fData_]
   Float_t         HPidParticleSim_fData_itsTrackData_fRKChiSquare[kMaxHPidParticleSim_fData];   //[HPidParticleSim.fData_]
   Int_t           HPidParticleSim_fData_itsTrackData_iIndexClosestTracklet[kMaxHPidParticleSim_fData];   //[HPidParticleSim.fData_]
   TArrayI         HPidParticleSim_fData_itsTrackData_aTrackletClusInf0[kMaxHPidParticleSim_fData];
   TArrayI         HPidParticleSim_fData_itsTrackData_aTrackletClusInf1[kMaxHPidParticleSim_fData];
   TArrayF         HPidParticleSim_fData_itsTrackData_aTrackletDistances[kMaxHPidParticleSim_fData];
   Float_t         HPidParticleSim_fData_itsTrackData_qIOMatching[kMaxHPidParticleSim_fData][10];   //[HPidParticleSim.fData_]
   Short_t         HPidParticleSim_fData_itsTrackData_nPolarity[kMaxHPidParticleSim_fData][10];   //[HPidParticleSim.fData_]
   Float_t         HPidParticleSim_fData_itsTrackData_fMomenta[kMaxHPidParticleSim_fData][10];   //[HPidParticleSim.fData_]
   Float_t         HPidParticleSim_fData_itsTrackData_fMomError[kMaxHPidParticleSim_fData][10];   //[HPidParticleSim.fData_]
   Float_t         HPidParticleSim_fData_itsTrackData_fTrackR[kMaxHPidParticleSim_fData][10];   //[HPidParticleSim.fData_]
   Float_t         HPidParticleSim_fData_itsTrackData_fTrackZ[kMaxHPidParticleSim_fData][10];   //[HPidParticleSim.fData_]
   Float_t         HPidParticleSim_fData_itsTrackData_fRKPhi[kMaxHPidParticleSim_fData];   //[HPidParticleSim.fData_]
   Float_t         HPidParticleSim_fData_itsTrackData_fRKTheta[kMaxHPidParticleSim_fData];   //[HPidParticleSim.fData_]
   Float_t         HPidParticleSim_fData_itsTrackData_fCorrectedEloss[kMaxHPidParticleSim_fData][10];   //[HPidParticleSim.fData_]
   Float_t         HPidParticleSim_fData_itsTrackData_fCorrectedBeta[kMaxHPidParticleSim_fData][10];   //[HPidParticleSim.fData_]
   Float_t         HPidParticleSim_fData_itsTrackData_fPathLength[kMaxHPidParticleSim_fData][10];   //[HPidParticleSim.fData_]
   Float_t         HPidParticleSim_fData_itsTrackData_fMassSquared[kMaxHPidParticleSim_fData][10];   //[HPidParticleSim.fData_]
   Bool_t          HPidParticleSim_fData_itsTrackData_bIsAccepted[kMaxHPidParticleSim_fData][10];   //[HPidParticleSim.fData_]
   Short_t         HPidParticleSim_fData_itsTrackData_nTofRecFlag[kMaxHPidParticleSim_fData][10];   //[HPidParticleSim.fData_]
   Float_t         HPidParticleSim_fData_itsTrackData_fTofRecTof[kMaxHPidParticleSim_fData][10];   //[HPidParticleSim.fData_]
   Float_t         HPidParticleSim_fData_itsTrackData_fTofRecBeta[kMaxHPidParticleSim_fData][10];   //[HPidParticleSim.fData_]
   Float_t         HPidParticleSim_fData_itsTrackData_fTofRecMassSquared[kMaxHPidParticleSim_fData][10];   //[HPidParticleSim.fData_]
   Int_t           HPidParticleSim_fData_flags[kMaxHPidParticleSim_fData];   //[HPidParticleSim.fData_]
   UInt_t          HPidParticleSim_fData_itsGeantTrackSet_fUniqueID[kMaxHPidParticleSim_fData];   //[HPidParticleSim.fData_]
   UInt_t          HPidParticleSim_fData_itsGeantTrackSet_fBits[kMaxHPidParticleSim_fData];   //[HPidParticleSim.fData_]
   Bool_t          HPidParticleSim_fData_itsGeantTrackSet_isSorted[kMaxHPidParticleSim_fData];   //[HPidParticleSim.fData_]
   Short_t         HPidParticleSim_fData_itsGeantTrackSet_sNCorrTrackIds[kMaxHPidParticleSim_fData];   //[HPidParticleSim.fData_]
   UInt_t          HPidParticleSim_fData_itsGeantTrackSet_nRICHTracks[kMaxHPidParticleSim_fData];   //[HPidParticleSim.fData_]
   UInt_t          HPidParticleSim_fData_itsGeantTrackSet_nRICHIPUTracks[kMaxHPidParticleSim_fData];   //[HPidParticleSim.fData_]
   UInt_t          HPidParticleSim_fData_itsGeantTrackSet_nInnerMdcTracks[kMaxHPidParticleSim_fData];   //[HPidParticleSim.fData_]
   UInt_t          HPidParticleSim_fData_itsGeantTrackSet_nOuterMdcTracks[kMaxHPidParticleSim_fData];   //[HPidParticleSim.fData_]
   UInt_t          HPidParticleSim_fData_itsGeantTrackSet_nShowerTracks[kMaxHPidParticleSim_fData];   //[HPidParticleSim.fData_]
   UInt_t          HPidParticleSim_fData_itsGeantTrackSet_nTOFTracks[kMaxHPidParticleSim_fData];   //[HPidParticleSim.fData_]
   Bool_t          HPidParticleSim_fData_itsGeantTrackSet_bIsLepFromPrimary[kMaxHPidParticleSim_fData];   //[HPidParticleSim.fData_]
   TArrayI         HPidParticleSim_fData_itsGeantTrackSet_correlatedTrackIds[kMaxHPidParticleSim_fData];
   TArrayI         HPidParticleSim_fData_itsGeantTrackSet_correlationFlags[kMaxHPidParticleSim_fData];
   TArrayI         HPidParticleSim_fData_itsGeantTrackSet_ProcessIds[kMaxHPidParticleSim_fData];
   TArrayI         HPidParticleSim_fData_itsGeantTrackSet_ParentIds[kMaxHPidParticleSim_fData];
   TArrayI         HPidParticleSim_fData_itsGeantTrackSet_Parents[kMaxHPidParticleSim_fData];
   TArrayI         HPidParticleSim_fData_itsGeantTrackSet_GrandParentIds[kMaxHPidParticleSim_fData];
   TArrayI         HPidParticleSim_fData_itsGeantTrackSet_GrandParents[kMaxHPidParticleSim_fData];
   TArrayF         HPidParticleSim_fData_itsGeantTrackSet_GenInfo[kMaxHPidParticleSim_fData];
   TArrayF         HPidParticleSim_fData_itsGeantTrackSet_GenInfo1[kMaxHPidParticleSim_fData];
   TArrayF         HPidParticleSim_fData_itsGeantTrackSet_GenInfo2[kMaxHPidParticleSim_fData];
   TArrayF         HPidParticleSim_fData_itsGeantTrackSet_GenWeight[kMaxHPidParticleSim_fData];
   TArrayF         HPidParticleSim_fData_itsGeantTrackSet_VertexX[kMaxHPidParticleSim_fData];
   TArrayF         HPidParticleSim_fData_itsGeantTrackSet_VertexY[kMaxHPidParticleSim_fData];
   TArrayF         HPidParticleSim_fData_itsGeantTrackSet_VertexZ[kMaxHPidParticleSim_fData];
   TArrayI         HPidParticleSim_fData_itsGeantTrackSet_GeantPIDs[kMaxHPidParticleSim_fData];
   TArrayI         HPidParticleSim_fData_itsGeantTrackSet_MediumIds[kMaxHPidParticleSim_fData];
   TArrayF         HPidParticleSim_fData_itsGeantTrackSet_GeantMomX[kMaxHPidParticleSim_fData];
   TArrayF         HPidParticleSim_fData_itsGeantTrackSet_GeantMomY[kMaxHPidParticleSim_fData];
   TArrayF         HPidParticleSim_fData_itsGeantTrackSet_GeantMomZ[kMaxHPidParticleSim_fData];
   TArrayF         HPidParticleSim_fData_itsGeantTrackSet_ShowerWeights[kMaxHPidParticleSim_fData];
   TArrayF         HPidParticleSim_fData_itsGeantTrackSet_TOFWeights[kMaxHPidParticleSim_fData];
   TArrayF         HPidParticleSim_fData_itsGeantTrackSet_RICHWeights[kMaxHPidParticleSim_fData];
   TArrayF         HPidParticleSim_fData_itsGeantTrackSet_RICHIPUWeights[kMaxHPidParticleSim_fData];
   TArrayF         HPidParticleSim_fData_itsGeantTrackSet_InnerMDCWeights[kMaxHPidParticleSim_fData];
   TArrayF         HPidParticleSim_fData_itsGeantTrackSet_OuterMDCWeights[kMaxHPidParticleSim_fData];
   TArrayI         HPidParticleSim_fData_itsGeantTrackSet_hasHitInShower[kMaxHPidParticleSim_fData];
   TArrayI         HPidParticleSim_fData_itsGeantTrackSet_hasHitInTOF[kMaxHPidParticleSim_fData];
   Int_t           HPidParticleSim_fNDataObjs;
   Bool_t          HPidParticleSim_hasDynamicObjects;
 //HEventHeader    *EventHeader_;
   UInt_t          EventHeader_TObject_fUniqueID;
   UInt_t          EventHeader_TObject_fBits;
   UInt_t          EventHeader_fVertex_fUniqueID;
   UInt_t          EventHeader_fVertex_fBits;
   UInt_t          EventHeader_fVertex_pos_fUniqueID;
   UInt_t          EventHeader_fVertex_pos_fBits;
   Double_t        EventHeader_fVertex_pos_x;
   Double_t        EventHeader_fVertex_pos_y;
   Double_t        EventHeader_fVertex_pos_z;
   Short_t         EventHeader_fVertex_iterations;
   Float_t         EventHeader_fVertex_chi2;
   Float_t         EventHeader_fVertex_sumOfWeights;
   Int_t           EventHeader_timeInSpill;
   UInt_t          EventHeader_downscaling;
   UInt_t          EventHeader_downscalingFlag;
   UInt_t          EventHeader_fDate;
   UInt_t          EventHeader_fErrorBit;
   UInt_t          EventHeader_fEventDecoding;
   UInt_t          EventHeader_fEventPad;
   UInt_t          EventHeader_fEventRunNumber;
   UInt_t          EventHeader_fEventSeqNumber;
   UInt_t          EventHeader_fEventSize;
   UInt_t          EventHeader_fId;
   UInt_t          EventHeader_fTBit;
   UInt_t          EventHeader_fTime;
   UInt_t          EventHeader_fVersion;
   UInt_t          EventHeader_triggerDecision;
   UInt_t          EventHeader_triggerDecisionEmu;

   // List of branches
   TBranch        *b_Event_HEvent_fUniqueID;   //!
   TBranch        *b_Event_HEvent_fBits;   //!
   TBranch        *b_Event_fRecLevel;   //!
   TBranch        *b_Event_fNTracks;   //!
   TBranch        *b_Event_fTracks_;   //!
   TBranch        *b_Event_fTracks_fUniqueID;   //!
   TBranch        *b_Event_fTracks_fBits;   //!
   TBranch        *b_Event_fTracks_fP;   //!
   TBranch        *b_Mdc_HEvent_fUniqueID;   //!
   TBranch        *b_Mdc_HEvent_fBits;   //!
   TBranch        *b_Mdc_fRecLevel;   //!
   TBranch        *b_Mdc_fBaseCategory;   //!
   TBranch        *b_Rich_HEvent_fUniqueID;   //!
   TBranch        *b_Rich_HEvent_fBits;   //!
   TBranch        *b_Rich_fRecLevel;   //!
   TBranch        *b_Rich_fBaseCategory;   //!
   TBranch        *b_Shower_HEvent_fUniqueID;   //!
   TBranch        *b_Shower_HEvent_fBits;   //!
   TBranch        *b_Shower_fRecLevel;   //!
   TBranch        *b_Shower_fBaseCategory;   //!
   TBranch        *b_Tof_HEvent_fUniqueID;   //!
   TBranch        *b_Tof_HEvent_fBits;   //!
   TBranch        *b_Tof_fRecLevel;   //!
   TBranch        *b_Tof_fBaseCategory;   //!
   TBranch        *b_Tofino_HEvent_fUniqueID;   //!
   TBranch        *b_Tofino_HEvent_fBits;   //!
   TBranch        *b_Tofino_fRecLevel;   //!
   TBranch        *b_Tofino_fBaseCategory;   //!
   TBranch        *b_Simul_HEvent_fUniqueID;   //!
   TBranch        *b_Simul_HEvent_fBits;   //!
   TBranch        *b_Simul_fRecLevel;   //!
   TBranch        *b_Simul_fBaseCategory;   //!
   TBranch        *b_HGeantKine_HCategory_fUniqueID;   //!
   TBranch        *b_HGeantKine_HCategory_fBits;   //!
   TBranch        *b_HGeantKine_HCategory_fCat;   //!
   TBranch        *b_HGeantKine_HCategory_fBranchingLevel;   //!
   TBranch        *b_HGeantKine_fData_;   //!
   TBranch        *b_HGeantKine_fData_fUniqueID;   //!
   TBranch        *b_HGeantKine_fData_fBits;   //!
   TBranch        *b_HGeantKine_fData_trackNumber;   //!
   TBranch        *b_HGeantKine_fData_parentTrack;   //!
   TBranch        *b_HGeantKine_fData_particleID;   //!
   TBranch        *b_HGeantKine_fData_mediumNumber;   //!
   TBranch        *b_HGeantKine_fData_creationMechanism;   //!
   TBranch        *b_HGeantKine_fData_xVertex;   //!
   TBranch        *b_HGeantKine_fData_yVertex;   //!
   TBranch        *b_HGeantKine_fData_zVertex;   //!
   TBranch        *b_HGeantKine_fData_xMom;   //!
   TBranch        *b_HGeantKine_fData_yMom;   //!
   TBranch        *b_HGeantKine_fData_zMom;   //!
   TBranch        *b_HGeantKine_fData_generatorInfo;   //!
   TBranch        *b_HGeantKine_fData_generatorInfo1;   //!
   TBranch        *b_HGeantKine_fData_generatorInfo2;   //!
   TBranch        *b_HGeantKine_fData_generatorWeight;   //!
   TBranch        *b_HGeantKine_fData_firstRichHit;   //!
   TBranch        *b_HGeantKine_fData_firstMdcHit;   //!
   TBranch        *b_HGeantKine_fData_firstTofHit;   //!
   TBranch        *b_HGeantKine_fData_firstRpcHit;   //!
   TBranch        *b_HGeantKine_fData_firstShowerHit;   //!
   TBranch        *b_HGeantKine_fData_firstWallHit;   //!
   TBranch        *b_HGeantKine_fData_active;   //!
   TBranch        *b_HGeantKine_fData_suppressed;   //!
   TBranch        *b_HGeantKine_fData_userVal;   //!
   TBranch        *b_HGeantKine_fNDataObjs;   //!
   TBranch        *b_HGeantKine_hasDynamicObjects;   //!
   TBranch        *b_HGeantMdc_HCategory_fUniqueID;   //!
   TBranch        *b_HGeantMdc_HCategory_fBits;   //!
   TBranch        *b_HGeantMdc_HCategory_fCat;   //!
   TBranch        *b_HGeantMdc_HCategory_fBranchingLevel;   //!
   TBranch        *b_HGeantMdc_fNDataObjs;   //!
   TBranch        *b_HGeantMdc_fData_;   //!
   TBranch        *b_HGeantMdc_fData_fUniqueID;   //!
   TBranch        *b_HGeantMdc_fData_fBits;   //!
   TBranch        *b_HGeantMdc_fData_nextHit;   //!
   TBranch        *b_HGeantMdc_fData_trackNumber;   //!
   TBranch        *b_HGeantMdc_fData_xHit;   //!
   TBranch        *b_HGeantMdc_fData_yHit;   //!
   TBranch        *b_HGeantMdc_fData_thetaHit;   //!
   TBranch        *b_HGeantMdc_fData_phiHit;   //!
   TBranch        *b_HGeantMdc_fData_tofHit;   //!
   TBranch        *b_HGeantMdc_fData_momHit;   //!
   TBranch        *b_HGeantMdc_fData_sector;   //!
   TBranch        *b_HGeantMdc_fData_module;   //!
   TBranch        *b_HGeantMdc_fData_layer;   //!
   TBranch        *b_HGeantShower_HCategory_fUniqueID;   //!
   TBranch        *b_HGeantShower_HCategory_fBits;   //!
   TBranch        *b_HGeantShower_HCategory_fCat;   //!
   TBranch        *b_HGeantShower_HCategory_fBranchingLevel;   //!
   TBranch        *b_HGeantShower_fNDataObjs;   //!
   TBranch        *b_HGeantShower_fData_;   //!
   TBranch        *b_HGeantShower_fData_fUniqueID;   //!
   TBranch        *b_HGeantShower_fData_fBits;   //!
   TBranch        *b_HGeantShower_fData_nextHit;   //!
   TBranch        *b_HGeantShower_fData_trackNumber;   //!
   TBranch        *b_HGeantShower_fData_eHit;   //!
   TBranch        *b_HGeantShower_fData_xHit;   //!
   TBranch        *b_HGeantShower_fData_yHit;   //!
   TBranch        *b_HGeantShower_fData_thetaHit;   //!
   TBranch        *b_HGeantShower_fData_phiHit;   //!
   TBranch        *b_HGeantShower_fData_betaHit;   //!
   TBranch        *b_HGeantShower_fData_sector;   //!
   TBranch        *b_HGeantShower_fData_module;   //!
   TBranch        *b_HGeantTof_HCategory_fUniqueID;   //!
   TBranch        *b_HGeantTof_HCategory_fBits;   //!
   TBranch        *b_HGeantTof_HCategory_fCat;   //!
   TBranch        *b_HGeantTof_HCategory_fBranchingLevel;   //!
   TBranch        *b_HGeantTof_fNDataObjs;   //!
   TBranch        *b_HGeantTof_fData_;   //!
   TBranch        *b_HGeantTof_fData_fUniqueID;   //!
   TBranch        *b_HGeantTof_fData_fBits;   //!
   TBranch        *b_HGeantTof_fData_nextHit;   //!
   TBranch        *b_HGeantTof_fData_trackNumber;   //!
   TBranch        *b_HGeantTof_fData_trackLength;   //!
   TBranch        *b_HGeantTof_fData_eHit;   //!
   TBranch        *b_HGeantTof_fData_xHit;   //!
   TBranch        *b_HGeantTof_fData_yHit;   //!
   TBranch        *b_HGeantTof_fData_tofHit;   //!
   TBranch        *b_HGeantTof_fData_momHit;   //!
   TBranch        *b_HGeantTof_fData_sector;   //!
   TBranch        *b_HGeantTof_fData_module;   //!
   TBranch        *b_HGeantTof_fData_cell;   //!
   TBranch        *b_Tracks_HEvent_fUniqueID;   //!
   TBranch        *b_Tracks_HEvent_fBits;   //!
   TBranch        *b_Tracks_fRecLevel;   //!
   TBranch        *b_Tracks_fBaseCategory;   //!
   TBranch        *b_Pid_HEvent_fUniqueID;   //!
   TBranch        *b_Pid_HEvent_fBits;   //!
   TBranch        *b_Pid_fRecLevel;   //!
   TBranch        *b_Pid_fBaseCategory;   //!
   TBranch        *b_HPidTrackCandSim_HCategory_fUniqueID;   //!
   TBranch        *b_HPidTrackCandSim_HCategory_fBits;   //!
   TBranch        *b_HPidTrackCandSim_HCategory_fCat;   //!
   TBranch        *b_HPidTrackCandSim_HCategory_fBranchingLevel;   //!
   TBranch        *b_HPidTrackCandSim_fData_;   //!
   TBranch        *b_HPidTrackCandSim_fData_fUniqueID;   //!
   TBranch        *b_HPidTrackCandSim_fData_fBits;   //!
   TBranch        *b_HPidTrackCandSim_fData_itsHitData_fUniqueID;   //!
   TBranch        *b_HPidTrackCandSim_fData_itsHitData_fBits;   //!
   TBranch        *b_HPidTrackCandSim_fData_itsHitData_nSector;   //!
   TBranch        *b_HPidTrackCandSim_fData_itsHitData_iSystem;   //!
   TBranch        *b_HPidTrackCandSim_fData_itsHitData_nRingPadNr;   //!
   TBranch        *b_HPidTrackCandSim_fData_itsHitData_fRingCentroid;   //!
   TBranch        *b_HPidTrackCandSim_fData_itsHitData_fRichTheta;   //!
   TBranch        *b_HPidTrackCandSim_fData_itsHitData_fRichPhi;   //!
   TBranch        *b_HPidTrackCandSim_fData_itsHitData_nRingPatMat;   //!
   TBranch        *b_HPidTrackCandSim_fData_itsHitData_nRingHouTra;   //!
   TBranch        *b_HPidTrackCandSim_fData_itsHitData_nRingAmplitude;   //!
   TBranch        *b_HPidTrackCandSim_fData_itsHitData_nRingLocalMax4;   //!
   TBranch        *b_HPidTrackCandSim_fData_itsHitData_fInnerMdcChiSquare;   //!
   TBranch        *b_HPidTrackCandSim_fData_itsHitData_fInnerMdcdEdx;   //!
   TBranch        *b_HPidTrackCandSim_fData_itsHitData_fInnerMdcdEdxSigma;   //!
   TBranch        *b_HPidTrackCandSim_fData_itsHitData_fMdcRCoord;   //!
   TBranch        *b_HPidTrackCandSim_fData_itsHitData_fMdcZCoord;   //!
   TBranch        *b_HPidTrackCandSim_fData_itsHitData_fMdcTheta;   //!
   TBranch        *b_HPidTrackCandSim_fData_itsHitData_fMdcPhi;   //!
   TBranch        *b_HPidTrackCandSim_fData_itsHitData_fOuterMdcChiSquare;   //!
   TBranch        *b_HPidTrackCandSim_fData_itsHitData_fOuterMdcdEdx;   //!
   TBranch        *b_HPidTrackCandSim_fData_itsHitData_fOuterMdcdEdxSigma;   //!
   TBranch        *b_HPidTrackCandSim_fData_itsHitData_fCombinedMdcdEdx;   //!
   TBranch        *b_HPidTrackCandSim_fData_itsHitData_fCombinedMdcdEdxSigma;   //!
   TBranch        *b_HPidTrackCandSim_fData_itsHitData_iIPURingQuality;   //!
   TBranch        *b_HPidTrackCandSim_fData_itsHitData_iIPUVetoQuality;   //!
   TBranch        *b_HPidTrackCandSim_fData_itsHitData_fShowerSum;   //!
   TBranch        *b_HPidTrackCandSim_fData_itsHitData_nShowerClS;   //!
   TBranch        *b_HPidTrackCandSim_fData_itsHitData_nShowerRow;   //!
   TBranch        *b_HPidTrackCandSim_fData_itsHitData_nShowerCol;   //!
   TBranch        *b_HPidTrackCandSim_fData_itsHitData_fShowerTimeOfFlight;   //!
   TBranch        *b_HPidTrackCandSim_fData_itsHitData_fMetaLocalX;   //!
   TBranch        *b_HPidTrackCandSim_fData_itsHitData_fMetaLocalY;   //!
   TBranch        *b_HPidTrackCandSim_fData_itsHitData_fTOFTimeOfFlight;   //!
   TBranch        *b_HPidTrackCandSim_fData_itsHitData_fTOFLeftAmplitude;   //!
   TBranch        *b_HPidTrackCandSim_fData_itsHitData_fTOFRightAmplitude;   //!
   TBranch        *b_HPidTrackCandSim_fData_itsHitData_fTofEloss;   //!
   TBranch        *b_HPidTrackCandSim_fData_itsHitData_iTofinoMult;   //!
   TBranch        *b_HPidTrackCandSim_fData_itsHitData_nTofClsSize;   //!
   TBranch        *b_HPidTrackCandSim_fData_itsHitData_nMetaCell;   //!
   TBranch        *b_HPidTrackCandSim_fData_itsHitData_nTofCell;   //!
   TBranch        *b_HPidTrackCandSim_fData_itsHitData_nTofModule;   //!
   TBranch        *b_HPidTrackCandSim_fData_itsHitData_iIndRICH;   //!
   TBranch        *b_HPidTrackCandSim_fData_itsHitData_iIndRICHIPU;   //!
   TBranch        *b_HPidTrackCandSim_fData_itsHitData_iIndInnerSeg;   //!
   TBranch        *b_HPidTrackCandSim_fData_itsHitData_iIndOuterSeg;   //!
   TBranch        *b_HPidTrackCandSim_fData_itsHitData_iIndTOF;   //!
   TBranch        *b_HPidTrackCandSim_fData_itsHitData_iIndShower;   //!
   TBranch        *b_HPidTrackCandSim_fData_itsHitData_iIndClusInf0;   //!
   TBranch        *b_HPidTrackCandSim_fData_itsHitData_iIndClusInf1;   //!
   TBranch        *b_HPidTrackCandSim_fData_itsHitData_iIndClusInf2;   //!
   TBranch        *b_HPidTrackCandSim_fData_itsHitData_iIndClusInf3;   //!
   TBranch        *b_HPidTrackCandSim_fData_itsHitData_iIndMatch;   //!
   TBranch        *b_HPidTrackCandSim_fData_itsHitData_fDistanceToVertex;   //!
   TBranch        *b_HPidTrackCandSim_fData_itsHitData_hasRingCorrelation;   //!
   TBranch        *b_HPidTrackCandSim_fData_itsHitData_hasMetaTrackCorrelation;   //!
   TBranch        *b_HPidTrackCandSim_fData_itsTrackData_fUniqueID;   //!
   TBranch        *b_HPidTrackCandSim_fData_itsTrackData_fBits;   //!
   TBranch        *b_HPidTrackCandSim_fData_itsTrackData_cov;   //!
   TBranch        *b_HPidTrackCandSim_fData_itsTrackData_nBestMomAlg;   //!
   TBranch        *b_HPidTrackCandSim_fData_itsTrackData_nRKTrackInd;   //!
   TBranch        *b_HPidTrackCandSim_fData_itsTrackData_nKickTrackInd;   //!
   TBranch        *b_HPidTrackCandSim_fData_itsTrackData_nKickTrack123Ind;   //!
   TBranch        *b_HPidTrackCandSim_fData_itsTrackData_nRefTrackInd;   //!
   TBranch        *b_HPidTrackCandSim_fData_itsTrackData_nSplineTrackInd;   //!
   TBranch        *b_HPidTrackCandSim_fData_itsTrackData_fMetaMatchingQuality;   //!
   TBranch        *b_HPidTrackCandSim_fData_itsTrackData_fRKRichMatchingQuality;   //!
   TBranch        *b_HPidTrackCandSim_fData_itsTrackData_dxRkMeta;   //!
   TBranch        *b_HPidTrackCandSim_fData_itsTrackData_dyRkMeta;   //!
   TBranch        *b_HPidTrackCandSim_fData_itsTrackData_dzRkMeta;   //!
   TBranch        *b_HPidTrackCandSim_fData_itsTrackData_dxMdcMeta;   //!
   TBranch        *b_HPidTrackCandSim_fData_itsTrackData_dyMdcMeta;   //!
   TBranch        *b_HPidTrackCandSim_fData_itsTrackData_dzMdcMeta;   //!
   TBranch        *b_HPidTrackCandSim_fData_itsTrackData_xMeta;   //!
   TBranch        *b_HPidTrackCandSim_fData_itsTrackData_yMeta;   //!
   TBranch        *b_HPidTrackCandSim_fData_itsTrackData_zMeta;   //!
   TBranch        *b_HPidTrackCandSim_fData_itsTrackData_errXMeta;   //!
   TBranch        *b_HPidTrackCandSim_fData_itsTrackData_errYMeta;   //!
   TBranch        *b_HPidTrackCandSim_fData_itsTrackData_errZMeta;   //!
   TBranch        *b_HPidTrackCandSim_fData_itsTrackData_nCloseTracklets;   //!
   TBranch        *b_HPidTrackCandSim_fData_itsTrackData_fPull;   //!
   TBranch        *b_HPidTrackCandSim_fData_itsTrackData_fSplineChiSquare;   //!
   TBranch        *b_HPidTrackCandSim_fData_itsTrackData_fRKChiSquare;   //!
   TBranch        *b_HPidTrackCandSim_fData_itsTrackData_iIndexClosestTracklet;   //!
   TBranch        *b_HPidTrackCandSim_fData_itsTrackData_aTrackletClusInf0;   //!
   TBranch        *b_HPidTrackCandSim_fData_itsTrackData_aTrackletClusInf1;   //!
   TBranch        *b_HPidTrackCandSim_fData_itsTrackData_aTrackletDistances;   //!
   TBranch        *b_HPidTrackCandSim_fData_itsTrackData_qIOMatching;   //!
   TBranch        *b_HPidTrackCandSim_fData_itsTrackData_nPolarity;   //!
   TBranch        *b_HPidTrackCandSim_fData_itsTrackData_fMomenta;   //!
   TBranch        *b_HPidTrackCandSim_fData_itsTrackData_fMomError;   //!
   TBranch        *b_HPidTrackCandSim_fData_itsTrackData_fTrackR;   //!
   TBranch        *b_HPidTrackCandSim_fData_itsTrackData_fTrackZ;   //!
   TBranch        *b_HPidTrackCandSim_fData_itsTrackData_fRKPhi;   //!
   TBranch        *b_HPidTrackCandSim_fData_itsTrackData_fRKTheta;   //!
   TBranch        *b_HPidTrackCandSim_fData_itsTrackData_fCorrectedEloss;   //!
   TBranch        *b_HPidTrackCandSim_fData_itsTrackData_fCorrectedBeta;   //!
   TBranch        *b_HPidTrackCandSim_fData_itsTrackData_fPathLength;   //!
   TBranch        *b_HPidTrackCandSim_fData_itsTrackData_fMassSquared;   //!
   TBranch        *b_HPidTrackCandSim_fData_itsTrackData_bIsAccepted;   //!
   TBranch        *b_HPidTrackCandSim_fData_itsTrackData_nTofRecFlag;   //!
   TBranch        *b_HPidTrackCandSim_fData_itsTrackData_fTofRecTof;   //!
   TBranch        *b_HPidTrackCandSim_fData_itsTrackData_fTofRecBeta;   //!
   TBranch        *b_HPidTrackCandSim_fData_itsTrackData_fTofRecMassSquared;   //!
   TBranch        *b_HPidTrackCandSim_fData_flags;   //!
   TBranch        *b_HPidTrackCandSim_fData_itsGeantTrackSet_fUniqueID;   //!
   TBranch        *b_HPidTrackCandSim_fData_itsGeantTrackSet_fBits;   //!
   TBranch        *b_HPidTrackCandSim_fData_itsGeantTrackSet_isSorted;   //!
   TBranch        *b_HPidTrackCandSim_fData_itsGeantTrackSet_sNCorrTrackIds;   //!
   TBranch        *b_HPidTrackCandSim_fData_itsGeantTrackSet_nRICHTracks;   //!
   TBranch        *b_HPidTrackCandSim_fData_itsGeantTrackSet_nRICHIPUTracks;   //!
   TBranch        *b_HPidTrackCandSim_fData_itsGeantTrackSet_nInnerMdcTracks;   //!
   TBranch        *b_HPidTrackCandSim_fData_itsGeantTrackSet_nOuterMdcTracks;   //!
   TBranch        *b_HPidTrackCandSim_fData_itsGeantTrackSet_nShowerTracks;   //!
   TBranch        *b_HPidTrackCandSim_fData_itsGeantTrackSet_nTOFTracks;   //!
   TBranch        *b_HPidTrackCandSim_fData_itsGeantTrackSet_bIsLepFromPrimary;   //!
   TBranch        *b_HPidTrackCandSim_fData_itsGeantTrackSet_correlatedTrackIds;   //!
   TBranch        *b_HPidTrackCandSim_fData_itsGeantTrackSet_correlationFlags;   //!
   TBranch        *b_HPidTrackCandSim_fData_itsGeantTrackSet_ProcessIds;   //!
   TBranch        *b_HPidTrackCandSim_fData_itsGeantTrackSet_ParentIds;   //!
   TBranch        *b_HPidTrackCandSim_fData_itsGeantTrackSet_Parents;   //!
   TBranch        *b_HPidTrackCandSim_fData_itsGeantTrackSet_GrandParentIds;   //!
   TBranch        *b_HPidTrackCandSim_fData_itsGeantTrackSet_GrandParents;   //!
   TBranch        *b_HPidTrackCandSim_fData_itsGeantTrackSet_GenInfo;   //!
   TBranch        *b_HPidTrackCandSim_fData_itsGeantTrackSet_GenInfo1;   //!
   TBranch        *b_HPidTrackCandSim_fData_itsGeantTrackSet_GenInfo2;   //!
   TBranch        *b_HPidTrackCandSim_fData_itsGeantTrackSet_GenWeight;   //!
   TBranch        *b_HPidTrackCandSim_fData_itsGeantTrackSet_VertexX;   //!
   TBranch        *b_HPidTrackCandSim_fData_itsGeantTrackSet_VertexY;   //!
   TBranch        *b_HPidTrackCandSim_fData_itsGeantTrackSet_VertexZ;   //!
   TBranch        *b_HPidTrackCandSim_fData_itsGeantTrackSet_GeantPIDs;   //!
   TBranch        *b_HPidTrackCandSim_fData_itsGeantTrackSet_MediumIds;   //!
   TBranch        *b_HPidTrackCandSim_fData_itsGeantTrackSet_GeantMomX;   //!
   TBranch        *b_HPidTrackCandSim_fData_itsGeantTrackSet_GeantMomY;   //!
   TBranch        *b_HPidTrackCandSim_fData_itsGeantTrackSet_GeantMomZ;   //!
   TBranch        *b_HPidTrackCandSim_fData_itsGeantTrackSet_ShowerWeights;   //!
   TBranch        *b_HPidTrackCandSim_fData_itsGeantTrackSet_TOFWeights;   //!
   TBranch        *b_HPidTrackCandSim_fData_itsGeantTrackSet_RICHWeights;   //!
   TBranch        *b_HPidTrackCandSim_fData_itsGeantTrackSet_RICHIPUWeights;   //!
   TBranch        *b_HPidTrackCandSim_fData_itsGeantTrackSet_InnerMDCWeights;   //!
   TBranch        *b_HPidTrackCandSim_fData_itsGeantTrackSet_OuterMDCWeights;   //!
   TBranch        *b_HPidTrackCandSim_fData_itsGeantTrackSet_hasHitInShower;   //!
   TBranch        *b_HPidTrackCandSim_fData_itsGeantTrackSet_hasHitInTOF;   //!
   TBranch        *b_HPidTrackCandSim_fNDataObjs;   //!
   TBranch        *b_HPidTrackCandSim_hasDynamicObjects;   //!
   TBranch        *b_HPidCandidate_HCategory_fUniqueID;   //!
   TBranch        *b_HPidCandidate_HCategory_fBits;   //!
   TBranch        *b_HPidCandidate_HCategory_fCat;   //!
   TBranch        *b_HPidCandidate_HCategory_fBranchingLevel;   //!
   TBranch        *b_HPidCandidate_fData_;   //!
   TBranch        *b_HPidCandidate_fData_fUniqueID;   //!
   TBranch        *b_HPidCandidate_fData_fBits;   //!
   TBranch        *b_HPidCandidate_fData_iTrackCandIndex;   //!
   TBranch        *b_HPidCandidate_fData_NUM_ALGORITHMS;   //!
   TBranch        *b_HPidCandidate_fData_NUM_PARTICLES;   //!
   TBranch        *b_HPidCandidate_fData_NUM_VALUES;   //!
   TBranch        *b_HPidCandidate_fData_aAlgorithms;   //!
   TBranch        *b_HPidCandidate_fData_aParticles;   //!
   TBranch        *b_HPidCandidate_fData_aValues;   //!
   TBranch        *b_HPidCandidate_fData_nMomAlgIndex;   //!
   TBranch        *b_HPidCandidate_fData_flags;   //!
   TBranch        *b_HPidCandidate_fNDataObjs;   //!
   TBranch        *b_HPidCandidate_hasDynamicObjects;   //!
   TBranch        *b_HPidParticleSim_HCategory_fUniqueID;   //!
   TBranch        *b_HPidParticleSim_HCategory_fBits;   //!
   TBranch        *b_HPidParticleSim_HCategory_fCat;   //!
   TBranch        *b_HPidParticleSim_HCategory_fBranchingLevel;   //!
   TBranch        *b_HPidParticleSim_fData_;   //!
   TBranch        *b_HPidParticleSim_fData_fUniqueID;   //!
   TBranch        *b_HPidParticleSim_fData_fBits;   //!
   TBranch        *b_HPidParticleSim_fData_fP_fUniqueID;   //!
   TBranch        *b_HPidParticleSim_fData_fP_fBits;   //!
   TBranch        *b_HPidParticleSim_fData_fP_fX;   //!
   TBranch        *b_HPidParticleSim_fData_fP_fY;   //!
   TBranch        *b_HPidParticleSim_fData_fP_fZ;   //!
   TBranch        *b_HPidParticleSim_fData_fE;   //!
   TBranch        *b_HPidParticleSim_fData_kUsesIdealMass;   //!
   TBranch        *b_HPidParticleSim_fData_nPossibleSpecies;   //!
   TBranch        *b_HPidParticleSim_fData_momAlgIndex;   //!
   TBranch        *b_HPidParticleSim_fData_nPidCandidateIndex;   //!
   TBranch        *b_HPidParticleSim_fData_possibleSpecies;   //!
   TBranch        *b_HPidParticleSim_fData_assignedWeights;   //!
   TBranch        *b_HPidParticleSim_fData_nAssignedPID;   //!
   TBranch        *b_HPidParticleSim_fData_fTestVal;   //!
   TBranch        *b_HPidParticleSim_fData_fWeight;   //!
   TBranch        *b_HPidParticleSim_fData_fMomRescal;   //!
   TBranch        *b_HPidParticleSim_fData_itsHitData_fUniqueID;   //!
   TBranch        *b_HPidParticleSim_fData_itsHitData_fBits;   //!
   TBranch        *b_HPidParticleSim_fData_itsHitData_nSector;   //!
   TBranch        *b_HPidParticleSim_fData_itsHitData_iSystem;   //!
   TBranch        *b_HPidParticleSim_fData_itsHitData_nRingPadNr;   //!
   TBranch        *b_HPidParticleSim_fData_itsHitData_fRingCentroid;   //!
   TBranch        *b_HPidParticleSim_fData_itsHitData_fRichTheta;   //!
   TBranch        *b_HPidParticleSim_fData_itsHitData_fRichPhi;   //!
   TBranch        *b_HPidParticleSim_fData_itsHitData_nRingPatMat;   //!
   TBranch        *b_HPidParticleSim_fData_itsHitData_nRingHouTra;   //!
   TBranch        *b_HPidParticleSim_fData_itsHitData_nRingAmplitude;   //!
   TBranch        *b_HPidParticleSim_fData_itsHitData_nRingLocalMax4;   //!
   TBranch        *b_HPidParticleSim_fData_itsHitData_fInnerMdcChiSquare;   //!
   TBranch        *b_HPidParticleSim_fData_itsHitData_fInnerMdcdEdx;   //!
   TBranch        *b_HPidParticleSim_fData_itsHitData_fInnerMdcdEdxSigma;   //!
   TBranch        *b_HPidParticleSim_fData_itsHitData_fMdcRCoord;   //!
   TBranch        *b_HPidParticleSim_fData_itsHitData_fMdcZCoord;   //!
   TBranch        *b_HPidParticleSim_fData_itsHitData_fMdcTheta;   //!
   TBranch        *b_HPidParticleSim_fData_itsHitData_fMdcPhi;   //!
   TBranch        *b_HPidParticleSim_fData_itsHitData_fOuterMdcChiSquare;   //!
   TBranch        *b_HPidParticleSim_fData_itsHitData_fOuterMdcdEdx;   //!
   TBranch        *b_HPidParticleSim_fData_itsHitData_fOuterMdcdEdxSigma;   //!
   TBranch        *b_HPidParticleSim_fData_itsHitData_fCombinedMdcdEdx;   //!
   TBranch        *b_HPidParticleSim_fData_itsHitData_fCombinedMdcdEdxSigma;   //!
   TBranch        *b_HPidParticleSim_fData_itsHitData_iIPURingQuality;   //!
   TBranch        *b_HPidParticleSim_fData_itsHitData_iIPUVetoQuality;   //!
   TBranch        *b_HPidParticleSim_fData_itsHitData_fShowerSum;   //!
   TBranch        *b_HPidParticleSim_fData_itsHitData_nShowerClS;   //!
   TBranch        *b_HPidParticleSim_fData_itsHitData_nShowerRow;   //!
   TBranch        *b_HPidParticleSim_fData_itsHitData_nShowerCol;   //!
   TBranch        *b_HPidParticleSim_fData_itsHitData_fShowerTimeOfFlight;   //!
   TBranch        *b_HPidParticleSim_fData_itsHitData_fMetaLocalX;   //!
   TBranch        *b_HPidParticleSim_fData_itsHitData_fMetaLocalY;   //!
   TBranch        *b_HPidParticleSim_fData_itsHitData_fTOFTimeOfFlight;   //!
   TBranch        *b_HPidParticleSim_fData_itsHitData_fTOFLeftAmplitude;   //!
   TBranch        *b_HPidParticleSim_fData_itsHitData_fTOFRightAmplitude;   //!
   TBranch        *b_HPidParticleSim_fData_itsHitData_fTofEloss;   //!
   TBranch        *b_HPidParticleSim_fData_itsHitData_iTofinoMult;   //!
   TBranch        *b_HPidParticleSim_fData_itsHitData_nTofClsSize;   //!
   TBranch        *b_HPidParticleSim_fData_itsHitData_nMetaCell;   //!
   TBranch        *b_HPidParticleSim_fData_itsHitData_nTofCell;   //!
   TBranch        *b_HPidParticleSim_fData_itsHitData_nTofModule;   //!
   TBranch        *b_HPidParticleSim_fData_itsHitData_iIndRICH;   //!
   TBranch        *b_HPidParticleSim_fData_itsHitData_iIndRICHIPU;   //!
   TBranch        *b_HPidParticleSim_fData_itsHitData_iIndInnerSeg;   //!
   TBranch        *b_HPidParticleSim_fData_itsHitData_iIndOuterSeg;   //!
   TBranch        *b_HPidParticleSim_fData_itsHitData_iIndTOF;   //!
   TBranch        *b_HPidParticleSim_fData_itsHitData_iIndShower;   //!
   TBranch        *b_HPidParticleSim_fData_itsHitData_iIndClusInf0;   //!
   TBranch        *b_HPidParticleSim_fData_itsHitData_iIndClusInf1;   //!
   TBranch        *b_HPidParticleSim_fData_itsHitData_iIndClusInf2;   //!
   TBranch        *b_HPidParticleSim_fData_itsHitData_iIndClusInf3;   //!
   TBranch        *b_HPidParticleSim_fData_itsHitData_iIndMatch;   //!
   TBranch        *b_HPidParticleSim_fData_itsHitData_fDistanceToVertex;   //!
   TBranch        *b_HPidParticleSim_fData_itsHitData_hasRingCorrelation;   //!
   TBranch        *b_HPidParticleSim_fData_itsHitData_hasMetaTrackCorrelation;   //!
   TBranch        *b_HPidParticleSim_fData_itsTrackData_fUniqueID;   //!
   TBranch        *b_HPidParticleSim_fData_itsTrackData_fBits;   //!
   TBranch        *b_HPidParticleSim_fData_itsTrackData_cov;   //!
   TBranch        *b_HPidParticleSim_fData_itsTrackData_nBestMomAlg;   //!
   TBranch        *b_HPidParticleSim_fData_itsTrackData_nRKTrackInd;   //!
   TBranch        *b_HPidParticleSim_fData_itsTrackData_nKickTrackInd;   //!
   TBranch        *b_HPidParticleSim_fData_itsTrackData_nKickTrack123Ind;   //!
   TBranch        *b_HPidParticleSim_fData_itsTrackData_nRefTrackInd;   //!
   TBranch        *b_HPidParticleSim_fData_itsTrackData_nSplineTrackInd;   //!
   TBranch        *b_HPidParticleSim_fData_itsTrackData_fMetaMatchingQuality;   //!
   TBranch        *b_HPidParticleSim_fData_itsTrackData_fRKRichMatchingQuality;   //!
   TBranch        *b_HPidParticleSim_fData_itsTrackData_dxRkMeta;   //!
   TBranch        *b_HPidParticleSim_fData_itsTrackData_dyRkMeta;   //!
   TBranch        *b_HPidParticleSim_fData_itsTrackData_dzRkMeta;   //!
   TBranch        *b_HPidParticleSim_fData_itsTrackData_dxMdcMeta;   //!
   TBranch        *b_HPidParticleSim_fData_itsTrackData_dyMdcMeta;   //!
   TBranch        *b_HPidParticleSim_fData_itsTrackData_dzMdcMeta;   //!
   TBranch        *b_HPidParticleSim_fData_itsTrackData_xMeta;   //!
   TBranch        *b_HPidParticleSim_fData_itsTrackData_yMeta;   //!
   TBranch        *b_HPidParticleSim_fData_itsTrackData_zMeta;   //!
   TBranch        *b_HPidParticleSim_fData_itsTrackData_errXMeta;   //!
   TBranch        *b_HPidParticleSim_fData_itsTrackData_errYMeta;   //!
   TBranch        *b_HPidParticleSim_fData_itsTrackData_errZMeta;   //!
   TBranch        *b_HPidParticleSim_fData_itsTrackData_nCloseTracklets;   //!
   TBranch        *b_HPidParticleSim_fData_itsTrackData_fPull;   //!
   TBranch        *b_HPidParticleSim_fData_itsTrackData_fSplineChiSquare;   //!
   TBranch        *b_HPidParticleSim_fData_itsTrackData_fRKChiSquare;   //!
   TBranch        *b_HPidParticleSim_fData_itsTrackData_iIndexClosestTracklet;   //!
   TBranch        *b_HPidParticleSim_fData_itsTrackData_aTrackletClusInf0;   //!
   TBranch        *b_HPidParticleSim_fData_itsTrackData_aTrackletClusInf1;   //!
   TBranch        *b_HPidParticleSim_fData_itsTrackData_aTrackletDistances;   //!
   TBranch        *b_HPidParticleSim_fData_itsTrackData_qIOMatching;   //!
   TBranch        *b_HPidParticleSim_fData_itsTrackData_nPolarity;   //!
   TBranch        *b_HPidParticleSim_fData_itsTrackData_fMomenta;   //!
   TBranch        *b_HPidParticleSim_fData_itsTrackData_fMomError;   //!
   TBranch        *b_HPidParticleSim_fData_itsTrackData_fTrackR;   //!
   TBranch        *b_HPidParticleSim_fData_itsTrackData_fTrackZ;   //!
   TBranch        *b_HPidParticleSim_fData_itsTrackData_fRKPhi;   //!
   TBranch        *b_HPidParticleSim_fData_itsTrackData_fRKTheta;   //!
   TBranch        *b_HPidParticleSim_fData_itsTrackData_fCorrectedEloss;   //!
   TBranch        *b_HPidParticleSim_fData_itsTrackData_fCorrectedBeta;   //!
   TBranch        *b_HPidParticleSim_fData_itsTrackData_fPathLength;   //!
   TBranch        *b_HPidParticleSim_fData_itsTrackData_fMassSquared;   //!
   TBranch        *b_HPidParticleSim_fData_itsTrackData_bIsAccepted;   //!
   TBranch        *b_HPidParticleSim_fData_itsTrackData_nTofRecFlag;   //!
   TBranch        *b_HPidParticleSim_fData_itsTrackData_fTofRecTof;   //!
   TBranch        *b_HPidParticleSim_fData_itsTrackData_fTofRecBeta;   //!
   TBranch        *b_HPidParticleSim_fData_itsTrackData_fTofRecMassSquared;   //!
   TBranch        *b_HPidParticleSim_fData_flags;   //!
   TBranch        *b_HPidParticleSim_fData_itsGeantTrackSet_fUniqueID;   //!
   TBranch        *b_HPidParticleSim_fData_itsGeantTrackSet_fBits;   //!
   TBranch        *b_HPidParticleSim_fData_itsGeantTrackSet_isSorted;   //!
   TBranch        *b_HPidParticleSim_fData_itsGeantTrackSet_sNCorrTrackIds;   //!
   TBranch        *b_HPidParticleSim_fData_itsGeantTrackSet_nRICHTracks;   //!
   TBranch        *b_HPidParticleSim_fData_itsGeantTrackSet_nRICHIPUTracks;   //!
   TBranch        *b_HPidParticleSim_fData_itsGeantTrackSet_nInnerMdcTracks;   //!
   TBranch        *b_HPidParticleSim_fData_itsGeantTrackSet_nOuterMdcTracks;   //!
   TBranch        *b_HPidParticleSim_fData_itsGeantTrackSet_nShowerTracks;   //!
   TBranch        *b_HPidParticleSim_fData_itsGeantTrackSet_nTOFTracks;   //!
   TBranch        *b_HPidParticleSim_fData_itsGeantTrackSet_bIsLepFromPrimary;   //!
   TBranch        *b_HPidParticleSim_fData_itsGeantTrackSet_correlatedTrackIds;   //!
   TBranch        *b_HPidParticleSim_fData_itsGeantTrackSet_correlationFlags;   //!
   TBranch        *b_HPidParticleSim_fData_itsGeantTrackSet_ProcessIds;   //!
   TBranch        *b_HPidParticleSim_fData_itsGeantTrackSet_ParentIds;   //!
   TBranch        *b_HPidParticleSim_fData_itsGeantTrackSet_Parents;   //!
   TBranch        *b_HPidParticleSim_fData_itsGeantTrackSet_GrandParentIds;   //!
   TBranch        *b_HPidParticleSim_fData_itsGeantTrackSet_GrandParents;   //!
   TBranch        *b_HPidParticleSim_fData_itsGeantTrackSet_GenInfo;   //!
   TBranch        *b_HPidParticleSim_fData_itsGeantTrackSet_GenInfo1;   //!
   TBranch        *b_HPidParticleSim_fData_itsGeantTrackSet_GenInfo2;   //!
   TBranch        *b_HPidParticleSim_fData_itsGeantTrackSet_GenWeight;   //!
   TBranch        *b_HPidParticleSim_fData_itsGeantTrackSet_VertexX;   //!
   TBranch        *b_HPidParticleSim_fData_itsGeantTrackSet_VertexY;   //!
   TBranch        *b_HPidParticleSim_fData_itsGeantTrackSet_VertexZ;   //!
   TBranch        *b_HPidParticleSim_fData_itsGeantTrackSet_GeantPIDs;   //!
   TBranch        *b_HPidParticleSim_fData_itsGeantTrackSet_MediumIds;   //!
   TBranch        *b_HPidParticleSim_fData_itsGeantTrackSet_GeantMomX;   //!
   TBranch        *b_HPidParticleSim_fData_itsGeantTrackSet_GeantMomY;   //!
   TBranch        *b_HPidParticleSim_fData_itsGeantTrackSet_GeantMomZ;   //!
   TBranch        *b_HPidParticleSim_fData_itsGeantTrackSet_ShowerWeights;   //!
   TBranch        *b_HPidParticleSim_fData_itsGeantTrackSet_TOFWeights;   //!
   TBranch        *b_HPidParticleSim_fData_itsGeantTrackSet_RICHWeights;   //!
   TBranch        *b_HPidParticleSim_fData_itsGeantTrackSet_RICHIPUWeights;   //!
   TBranch        *b_HPidParticleSim_fData_itsGeantTrackSet_InnerMDCWeights;   //!
   TBranch        *b_HPidParticleSim_fData_itsGeantTrackSet_OuterMDCWeights;   //!
   TBranch        *b_HPidParticleSim_fData_itsGeantTrackSet_hasHitInShower;   //!
   TBranch        *b_HPidParticleSim_fData_itsGeantTrackSet_hasHitInTOF;   //!
   TBranch        *b_HPidParticleSim_fNDataObjs;   //!
   TBranch        *b_HPidParticleSim_hasDynamicObjects;   //!
   TBranch        *b_EventHeader_TObject_fUniqueID;   //!
   TBranch        *b_EventHeader_TObject_fBits;   //!
   TBranch        *b_EventHeader_fVertex_fUniqueID;   //!
   TBranch        *b_EventHeader_fVertex_fBits;   //!
   TBranch        *b_EventHeader_fVertex_pos_fUniqueID;   //!
   TBranch        *b_EventHeader_fVertex_pos_fBits;   //!
   TBranch        *b_EventHeader_fVertex_pos_x;   //!
   TBranch        *b_EventHeader_fVertex_pos_y;   //!
   TBranch        *b_EventHeader_fVertex_pos_z;   //!
   TBranch        *b_EventHeader_fVertex_iterations;   //!
   TBranch        *b_EventHeader_fVertex_chi2;   //!
   TBranch        *b_EventHeader_fVertex_sumOfWeights;   //!
   TBranch        *b_EventHeader_timeInSpill;   //!
   TBranch        *b_EventHeader_downscaling;   //!
   TBranch        *b_EventHeader_downscalingFlag;   //!
   TBranch        *b_EventHeader_fDate;   //!
   TBranch        *b_EventHeader_fErrorBit;   //!
   TBranch        *b_EventHeader_fEventDecoding;   //!
   TBranch        *b_EventHeader_fEventPad;   //!
   TBranch        *b_EventHeader_fEventRunNumber;   //!
   TBranch        *b_EventHeader_fEventSeqNumber;   //!
   TBranch        *b_EventHeader_fEventSize;   //!
   TBranch        *b_EventHeader_fId;   //!
   TBranch        *b_EventHeader_fTBit;   //!
   TBranch        *b_EventHeader_fTime;   //!
   TBranch        *b_EventHeader_fVersion;   //!
   TBranch        *b_EventHeader_triggerDecision;   //!
   TBranch        *b_EventHeader_triggerDecisionEmu;   //!

   T(TTree *tree=0);
   virtual ~T();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef T_cxx
T::T(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("/lustre/nyx/hades/user/iciepal/pNb/dst/FILES/lambda1520_100k_01_dst.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("/lustre/nyx/hades/user/iciepal/pNb/dst/FILES/lambda1520_100k_01_dst.root");
      }
      f->GetObject("T",tree);

   }
   Init(tree);
}

T::~T()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t T::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t T::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void T::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("Event.HEvent.fUniqueID", &Event_HEvent_fUniqueID, &b_Event_HEvent_fUniqueID);
   fChain->SetBranchAddress("Event.HEvent.fBits", &Event_HEvent_fBits, &b_Event_HEvent_fBits);
   fChain->SetBranchAddress("Event.fRecLevel", &Event_fRecLevel, &b_Event_fRecLevel);
   fChain->SetBranchAddress("Event.fNTracks", &Event_fNTracks, &b_Event_fNTracks);
   fChain->SetBranchAddress("Event.fTracks", &Event_fTracks_, &b_Event_fTracks_);
   fChain->SetBranchAddress("Event.fTracks.fUniqueID", &Event_fTracks_fUniqueID, &b_Event_fTracks_fUniqueID);
   fChain->SetBranchAddress("Event.fTracks.fBits", &Event_fTracks_fBits, &b_Event_fTracks_fBits);
   fChain->SetBranchAddress("Event.fTracks.fP", &Event_fTracks_fP, &b_Event_fTracks_fP);
   fChain->SetBranchAddress("Mdc.HEvent.fUniqueID", &Mdc_HEvent_fUniqueID, &b_Mdc_HEvent_fUniqueID);
   fChain->SetBranchAddress("Mdc.HEvent.fBits", &Mdc_HEvent_fBits, &b_Mdc_HEvent_fBits);
   fChain->SetBranchAddress("Mdc.fRecLevel", &Mdc_fRecLevel, &b_Mdc_fRecLevel);
   fChain->SetBranchAddress("Mdc.fBaseCategory", &Mdc_fBaseCategory, &b_Mdc_fBaseCategory);
   fChain->SetBranchAddress("Rich.HEvent.fUniqueID", &Rich_HEvent_fUniqueID, &b_Rich_HEvent_fUniqueID);
   fChain->SetBranchAddress("Rich.HEvent.fBits", &Rich_HEvent_fBits, &b_Rich_HEvent_fBits);
   fChain->SetBranchAddress("Rich.fRecLevel", &Rich_fRecLevel, &b_Rich_fRecLevel);
   fChain->SetBranchAddress("Rich.fBaseCategory", &Rich_fBaseCategory, &b_Rich_fBaseCategory);
   fChain->SetBranchAddress("Shower.HEvent.fUniqueID", &Shower_HEvent_fUniqueID, &b_Shower_HEvent_fUniqueID);
   fChain->SetBranchAddress("Shower.HEvent.fBits", &Shower_HEvent_fBits, &b_Shower_HEvent_fBits);
   fChain->SetBranchAddress("Shower.fRecLevel", &Shower_fRecLevel, &b_Shower_fRecLevel);
   fChain->SetBranchAddress("Shower.fBaseCategory", &Shower_fBaseCategory, &b_Shower_fBaseCategory);
   fChain->SetBranchAddress("Tof.HEvent.fUniqueID", &Tof_HEvent_fUniqueID, &b_Tof_HEvent_fUniqueID);
   fChain->SetBranchAddress("Tof.HEvent.fBits", &Tof_HEvent_fBits, &b_Tof_HEvent_fBits);
   fChain->SetBranchAddress("Tof.fRecLevel", &Tof_fRecLevel, &b_Tof_fRecLevel);
   fChain->SetBranchAddress("Tof.fBaseCategory", &Tof_fBaseCategory, &b_Tof_fBaseCategory);
   fChain->SetBranchAddress("Tofino.HEvent.fUniqueID", &Tofino_HEvent_fUniqueID, &b_Tofino_HEvent_fUniqueID);
   fChain->SetBranchAddress("Tofino.HEvent.fBits", &Tofino_HEvent_fBits, &b_Tofino_HEvent_fBits);
   fChain->SetBranchAddress("Tofino.fRecLevel", &Tofino_fRecLevel, &b_Tofino_fRecLevel);
   fChain->SetBranchAddress("Tofino.fBaseCategory", &Tofino_fBaseCategory, &b_Tofino_fBaseCategory);
   fChain->SetBranchAddress("Simul.HEvent.fUniqueID", &Simul_HEvent_fUniqueID, &b_Simul_HEvent_fUniqueID);
   fChain->SetBranchAddress("Simul.HEvent.fBits", &Simul_HEvent_fBits, &b_Simul_HEvent_fBits);
   fChain->SetBranchAddress("Simul.fRecLevel", &Simul_fRecLevel, &b_Simul_fRecLevel);
   fChain->SetBranchAddress("Simul.fBaseCategory", &Simul_fBaseCategory, &b_Simul_fBaseCategory);
   fChain->SetBranchAddress("HGeantKine.HCategory.fUniqueID", &HGeantKine_HCategory_fUniqueID, &b_HGeantKine_HCategory_fUniqueID);
   fChain->SetBranchAddress("HGeantKine.HCategory.fBits", &HGeantKine_HCategory_fBits, &b_HGeantKine_HCategory_fBits);
   fChain->SetBranchAddress("HGeantKine.HCategory.fCat", &HGeantKine_HCategory_fCat, &b_HGeantKine_HCategory_fCat);
   fChain->SetBranchAddress("HGeantKine.HCategory.fBranchingLevel", &HGeantKine_HCategory_fBranchingLevel, &b_HGeantKine_HCategory_fBranchingLevel);
   fChain->SetBranchAddress("HGeantKine.fData", &HGeantKine_fData_, &b_HGeantKine_fData_);
   fChain->SetBranchAddress("HGeantKine.fData.fUniqueID", HGeantKine_fData_fUniqueID, &b_HGeantKine_fData_fUniqueID);
   fChain->SetBranchAddress("HGeantKine.fData.fBits", HGeantKine_fData_fBits, &b_HGeantKine_fData_fBits);
   fChain->SetBranchAddress("HGeantKine.fData.trackNumber", HGeantKine_fData_trackNumber, &b_HGeantKine_fData_trackNumber);
   fChain->SetBranchAddress("HGeantKine.fData.parentTrack", HGeantKine_fData_parentTrack, &b_HGeantKine_fData_parentTrack);
   fChain->SetBranchAddress("HGeantKine.fData.particleID", HGeantKine_fData_particleID, &b_HGeantKine_fData_particleID);
   fChain->SetBranchAddress("HGeantKine.fData.mediumNumber", HGeantKine_fData_mediumNumber, &b_HGeantKine_fData_mediumNumber);
   fChain->SetBranchAddress("HGeantKine.fData.creationMechanism", HGeantKine_fData_creationMechanism, &b_HGeantKine_fData_creationMechanism);
   fChain->SetBranchAddress("HGeantKine.fData.xVertex", HGeantKine_fData_xVertex, &b_HGeantKine_fData_xVertex);
   fChain->SetBranchAddress("HGeantKine.fData.yVertex", HGeantKine_fData_yVertex, &b_HGeantKine_fData_yVertex);
   fChain->SetBranchAddress("HGeantKine.fData.zVertex", HGeantKine_fData_zVertex, &b_HGeantKine_fData_zVertex);
   fChain->SetBranchAddress("HGeantKine.fData.xMom", HGeantKine_fData_xMom, &b_HGeantKine_fData_xMom);
   fChain->SetBranchAddress("HGeantKine.fData.yMom", HGeantKine_fData_yMom, &b_HGeantKine_fData_yMom);
   fChain->SetBranchAddress("HGeantKine.fData.zMom", HGeantKine_fData_zMom, &b_HGeantKine_fData_zMom);
   fChain->SetBranchAddress("HGeantKine.fData.generatorInfo", HGeantKine_fData_generatorInfo, &b_HGeantKine_fData_generatorInfo);
   fChain->SetBranchAddress("HGeantKine.fData.generatorInfo1", HGeantKine_fData_generatorInfo1, &b_HGeantKine_fData_generatorInfo1);
   fChain->SetBranchAddress("HGeantKine.fData.generatorInfo2", HGeantKine_fData_generatorInfo2, &b_HGeantKine_fData_generatorInfo2);
   fChain->SetBranchAddress("HGeantKine.fData.generatorWeight", HGeantKine_fData_generatorWeight, &b_HGeantKine_fData_generatorWeight);
   fChain->SetBranchAddress("HGeantKine.fData.firstRichHit", HGeantKine_fData_firstRichHit, &b_HGeantKine_fData_firstRichHit);
   fChain->SetBranchAddress("HGeantKine.fData.firstMdcHit", HGeantKine_fData_firstMdcHit, &b_HGeantKine_fData_firstMdcHit);
   fChain->SetBranchAddress("HGeantKine.fData.firstTofHit", HGeantKine_fData_firstTofHit, &b_HGeantKine_fData_firstTofHit);
   fChain->SetBranchAddress("HGeantKine.fData.firstRpcHit", HGeantKine_fData_firstRpcHit, &b_HGeantKine_fData_firstRpcHit);
   fChain->SetBranchAddress("HGeantKine.fData.firstShowerHit", HGeantKine_fData_firstShowerHit, &b_HGeantKine_fData_firstShowerHit);
   fChain->SetBranchAddress("HGeantKine.fData.firstWallHit", HGeantKine_fData_firstWallHit, &b_HGeantKine_fData_firstWallHit);
   fChain->SetBranchAddress("HGeantKine.fData.active", HGeantKine_fData_active, &b_HGeantKine_fData_active);
   fChain->SetBranchAddress("HGeantKine.fData.suppressed", HGeantKine_fData_suppressed, &b_HGeantKine_fData_suppressed);
   fChain->SetBranchAddress("HGeantKine.fData.userVal", HGeantKine_fData_userVal, &b_HGeantKine_fData_userVal);
   fChain->SetBranchAddress("HGeantKine.fNDataObjs", &HGeantKine_fNDataObjs, &b_HGeantKine_fNDataObjs);
   fChain->SetBranchAddress("HGeantKine.hasDynamicObjects", &HGeantKine_hasDynamicObjects, &b_HGeantKine_hasDynamicObjects);
   fChain->SetBranchAddress("HGeantMdc.HCategory.fUniqueID", &HGeantMdc_HCategory_fUniqueID, &b_HGeantMdc_HCategory_fUniqueID);
   fChain->SetBranchAddress("HGeantMdc.HCategory.fBits", &HGeantMdc_HCategory_fBits, &b_HGeantMdc_HCategory_fBits);
   fChain->SetBranchAddress("HGeantMdc.HCategory.fCat", &HGeantMdc_HCategory_fCat, &b_HGeantMdc_HCategory_fCat);
   fChain->SetBranchAddress("HGeantMdc.HCategory.fBranchingLevel", &HGeantMdc_HCategory_fBranchingLevel, &b_HGeantMdc_HCategory_fBranchingLevel);
   fChain->SetBranchAddress("HGeantMdc.fNDataObjs", &HGeantMdc_fNDataObjs, &b_HGeantMdc_fNDataObjs);
   fChain->SetBranchAddress("HGeantMdc.fData", &HGeantMdc_fData_, &b_HGeantMdc_fData_);
   fChain->SetBranchAddress("HGeantMdc.fData.fUniqueID", HGeantMdc_fData_fUniqueID, &b_HGeantMdc_fData_fUniqueID);
   fChain->SetBranchAddress("HGeantMdc.fData.fBits", HGeantMdc_fData_fBits, &b_HGeantMdc_fData_fBits);
   fChain->SetBranchAddress("HGeantMdc.fData.nextHit", HGeantMdc_fData_nextHit, &b_HGeantMdc_fData_nextHit);
   fChain->SetBranchAddress("HGeantMdc.fData.trackNumber", HGeantMdc_fData_trackNumber, &b_HGeantMdc_fData_trackNumber);
   fChain->SetBranchAddress("HGeantMdc.fData.xHit", HGeantMdc_fData_xHit, &b_HGeantMdc_fData_xHit);
   fChain->SetBranchAddress("HGeantMdc.fData.yHit", HGeantMdc_fData_yHit, &b_HGeantMdc_fData_yHit);
   fChain->SetBranchAddress("HGeantMdc.fData.thetaHit", HGeantMdc_fData_thetaHit, &b_HGeantMdc_fData_thetaHit);
   fChain->SetBranchAddress("HGeantMdc.fData.phiHit", HGeantMdc_fData_phiHit, &b_HGeantMdc_fData_phiHit);
   fChain->SetBranchAddress("HGeantMdc.fData.tofHit", HGeantMdc_fData_tofHit, &b_HGeantMdc_fData_tofHit);
   fChain->SetBranchAddress("HGeantMdc.fData.momHit", HGeantMdc_fData_momHit, &b_HGeantMdc_fData_momHit);
   fChain->SetBranchAddress("HGeantMdc.fData.sector", HGeantMdc_fData_sector, &b_HGeantMdc_fData_sector);
   fChain->SetBranchAddress("HGeantMdc.fData.module", HGeantMdc_fData_module, &b_HGeantMdc_fData_module);
   fChain->SetBranchAddress("HGeantMdc.fData.layer", HGeantMdc_fData_layer, &b_HGeantMdc_fData_layer);
   fChain->SetBranchAddress("HGeantShower.HCategory.fUniqueID", &HGeantShower_HCategory_fUniqueID, &b_HGeantShower_HCategory_fUniqueID);
   fChain->SetBranchAddress("HGeantShower.HCategory.fBits", &HGeantShower_HCategory_fBits, &b_HGeantShower_HCategory_fBits);
   fChain->SetBranchAddress("HGeantShower.HCategory.fCat", &HGeantShower_HCategory_fCat, &b_HGeantShower_HCategory_fCat);
   fChain->SetBranchAddress("HGeantShower.HCategory.fBranchingLevel", &HGeantShower_HCategory_fBranchingLevel, &b_HGeantShower_HCategory_fBranchingLevel);
   fChain->SetBranchAddress("HGeantShower.fNDataObjs", &HGeantShower_fNDataObjs, &b_HGeantShower_fNDataObjs);
   fChain->SetBranchAddress("HGeantShower.fData", &HGeantShower_fData_, &b_HGeantShower_fData_);
   fChain->SetBranchAddress("HGeantShower.fData.fUniqueID", HGeantShower_fData_fUniqueID, &b_HGeantShower_fData_fUniqueID);
   fChain->SetBranchAddress("HGeantShower.fData.fBits", HGeantShower_fData_fBits, &b_HGeantShower_fData_fBits);
   fChain->SetBranchAddress("HGeantShower.fData.nextHit", HGeantShower_fData_nextHit, &b_HGeantShower_fData_nextHit);
   fChain->SetBranchAddress("HGeantShower.fData.trackNumber", HGeantShower_fData_trackNumber, &b_HGeantShower_fData_trackNumber);
   fChain->SetBranchAddress("HGeantShower.fData.eHit", HGeantShower_fData_eHit, &b_HGeantShower_fData_eHit);
   fChain->SetBranchAddress("HGeantShower.fData.xHit", HGeantShower_fData_xHit, &b_HGeantShower_fData_xHit);
   fChain->SetBranchAddress("HGeantShower.fData.yHit", HGeantShower_fData_yHit, &b_HGeantShower_fData_yHit);
   fChain->SetBranchAddress("HGeantShower.fData.thetaHit", HGeantShower_fData_thetaHit, &b_HGeantShower_fData_thetaHit);
   fChain->SetBranchAddress("HGeantShower.fData.phiHit", HGeantShower_fData_phiHit, &b_HGeantShower_fData_phiHit);
   fChain->SetBranchAddress("HGeantShower.fData.betaHit", HGeantShower_fData_betaHit, &b_HGeantShower_fData_betaHit);
   fChain->SetBranchAddress("HGeantShower.fData.sector", HGeantShower_fData_sector, &b_HGeantShower_fData_sector);
   fChain->SetBranchAddress("HGeantShower.fData.module", HGeantShower_fData_module, &b_HGeantShower_fData_module);
   fChain->SetBranchAddress("HGeantTof.HCategory.fUniqueID", &HGeantTof_HCategory_fUniqueID, &b_HGeantTof_HCategory_fUniqueID);
   fChain->SetBranchAddress("HGeantTof.HCategory.fBits", &HGeantTof_HCategory_fBits, &b_HGeantTof_HCategory_fBits);
   fChain->SetBranchAddress("HGeantTof.HCategory.fCat", &HGeantTof_HCategory_fCat, &b_HGeantTof_HCategory_fCat);
   fChain->SetBranchAddress("HGeantTof.HCategory.fBranchingLevel", &HGeantTof_HCategory_fBranchingLevel, &b_HGeantTof_HCategory_fBranchingLevel);
   fChain->SetBranchAddress("HGeantTof.fNDataObjs", &HGeantTof_fNDataObjs, &b_HGeantTof_fNDataObjs);
   fChain->SetBranchAddress("HGeantTof.fData", &HGeantTof_fData_, &b_HGeantTof_fData_);
   fChain->SetBranchAddress("HGeantTof.fData.fUniqueID", HGeantTof_fData_fUniqueID, &b_HGeantTof_fData_fUniqueID);
   fChain->SetBranchAddress("HGeantTof.fData.fBits", HGeantTof_fData_fBits, &b_HGeantTof_fData_fBits);
   fChain->SetBranchAddress("HGeantTof.fData.nextHit", HGeantTof_fData_nextHit, &b_HGeantTof_fData_nextHit);
   fChain->SetBranchAddress("HGeantTof.fData.trackNumber", HGeantTof_fData_trackNumber, &b_HGeantTof_fData_trackNumber);
   fChain->SetBranchAddress("HGeantTof.fData.trackLength", HGeantTof_fData_trackLength, &b_HGeantTof_fData_trackLength);
   fChain->SetBranchAddress("HGeantTof.fData.eHit", HGeantTof_fData_eHit, &b_HGeantTof_fData_eHit);
   fChain->SetBranchAddress("HGeantTof.fData.xHit", HGeantTof_fData_xHit, &b_HGeantTof_fData_xHit);
   fChain->SetBranchAddress("HGeantTof.fData.yHit", HGeantTof_fData_yHit, &b_HGeantTof_fData_yHit);
   fChain->SetBranchAddress("HGeantTof.fData.tofHit", HGeantTof_fData_tofHit, &b_HGeantTof_fData_tofHit);
   fChain->SetBranchAddress("HGeantTof.fData.momHit", HGeantTof_fData_momHit, &b_HGeantTof_fData_momHit);
   fChain->SetBranchAddress("HGeantTof.fData.sector", HGeantTof_fData_sector, &b_HGeantTof_fData_sector);
   fChain->SetBranchAddress("HGeantTof.fData.module", HGeantTof_fData_module, &b_HGeantTof_fData_module);
   fChain->SetBranchAddress("HGeantTof.fData.cell", HGeantTof_fData_cell, &b_HGeantTof_fData_cell);
   fChain->SetBranchAddress("Tracks.HEvent.fUniqueID", &Tracks_HEvent_fUniqueID, &b_Tracks_HEvent_fUniqueID);
   fChain->SetBranchAddress("Tracks.HEvent.fBits", &Tracks_HEvent_fBits, &b_Tracks_HEvent_fBits);
   fChain->SetBranchAddress("Tracks.fRecLevel", &Tracks_fRecLevel, &b_Tracks_fRecLevel);
   fChain->SetBranchAddress("Tracks.fBaseCategory", &Tracks_fBaseCategory, &b_Tracks_fBaseCategory);
   fChain->SetBranchAddress("Pid.HEvent.fUniqueID", &Pid_HEvent_fUniqueID, &b_Pid_HEvent_fUniqueID);
   fChain->SetBranchAddress("Pid.HEvent.fBits", &Pid_HEvent_fBits, &b_Pid_HEvent_fBits);
   fChain->SetBranchAddress("Pid.fRecLevel", &Pid_fRecLevel, &b_Pid_fRecLevel);
   fChain->SetBranchAddress("Pid.fBaseCategory", &Pid_fBaseCategory, &b_Pid_fBaseCategory);
   fChain->SetBranchAddress("HPidTrackCandSim.HCategory.fUniqueID", &HPidTrackCandSim_HCategory_fUniqueID, &b_HPidTrackCandSim_HCategory_fUniqueID);
   fChain->SetBranchAddress("HPidTrackCandSim.HCategory.fBits", &HPidTrackCandSim_HCategory_fBits, &b_HPidTrackCandSim_HCategory_fBits);
   fChain->SetBranchAddress("HPidTrackCandSim.HCategory.fCat", &HPidTrackCandSim_HCategory_fCat, &b_HPidTrackCandSim_HCategory_fCat);
   fChain->SetBranchAddress("HPidTrackCandSim.HCategory.fBranchingLevel", &HPidTrackCandSim_HCategory_fBranchingLevel, &b_HPidTrackCandSim_HCategory_fBranchingLevel);
   fChain->SetBranchAddress("HPidTrackCandSim.fData", &HPidTrackCandSim_fData_, &b_HPidTrackCandSim_fData_);
   fChain->SetBranchAddress("HPidTrackCandSim.fData.fUniqueID", HPidTrackCandSim_fData_fUniqueID, &b_HPidTrackCandSim_fData_fUniqueID);
   fChain->SetBranchAddress("HPidTrackCandSim.fData.fBits", HPidTrackCandSim_fData_fBits, &b_HPidTrackCandSim_fData_fBits);
   fChain->SetBranchAddress("HPidTrackCandSim.fData.itsHitData.fUniqueID", HPidTrackCandSim_fData_itsHitData_fUniqueID, &b_HPidTrackCandSim_fData_itsHitData_fUniqueID);
   fChain->SetBranchAddress("HPidTrackCandSim.fData.itsHitData.fBits", HPidTrackCandSim_fData_itsHitData_fBits, &b_HPidTrackCandSim_fData_itsHitData_fBits);
   fChain->SetBranchAddress("HPidTrackCandSim.fData.itsHitData.nSector", HPidTrackCandSim_fData_itsHitData_nSector, &b_HPidTrackCandSim_fData_itsHitData_nSector);
   fChain->SetBranchAddress("HPidTrackCandSim.fData.itsHitData.iSystem", HPidTrackCandSim_fData_itsHitData_iSystem, &b_HPidTrackCandSim_fData_itsHitData_iSystem);
   fChain->SetBranchAddress("HPidTrackCandSim.fData.itsHitData.nRingPadNr", HPidTrackCandSim_fData_itsHitData_nRingPadNr, &b_HPidTrackCandSim_fData_itsHitData_nRingPadNr);
   fChain->SetBranchAddress("HPidTrackCandSim.fData.itsHitData.fRingCentroid", HPidTrackCandSim_fData_itsHitData_fRingCentroid, &b_HPidTrackCandSim_fData_itsHitData_fRingCentroid);
   fChain->SetBranchAddress("HPidTrackCandSim.fData.itsHitData.fRichTheta", HPidTrackCandSim_fData_itsHitData_fRichTheta, &b_HPidTrackCandSim_fData_itsHitData_fRichTheta);
   fChain->SetBranchAddress("HPidTrackCandSim.fData.itsHitData.fRichPhi", HPidTrackCandSim_fData_itsHitData_fRichPhi, &b_HPidTrackCandSim_fData_itsHitData_fRichPhi);
   fChain->SetBranchAddress("HPidTrackCandSim.fData.itsHitData.nRingPatMat", HPidTrackCandSim_fData_itsHitData_nRingPatMat, &b_HPidTrackCandSim_fData_itsHitData_nRingPatMat);
   fChain->SetBranchAddress("HPidTrackCandSim.fData.itsHitData.nRingHouTra", HPidTrackCandSim_fData_itsHitData_nRingHouTra, &b_HPidTrackCandSim_fData_itsHitData_nRingHouTra);
   fChain->SetBranchAddress("HPidTrackCandSim.fData.itsHitData.nRingAmplitude", HPidTrackCandSim_fData_itsHitData_nRingAmplitude, &b_HPidTrackCandSim_fData_itsHitData_nRingAmplitude);
   fChain->SetBranchAddress("HPidTrackCandSim.fData.itsHitData.nRingLocalMax4", HPidTrackCandSim_fData_itsHitData_nRingLocalMax4, &b_HPidTrackCandSim_fData_itsHitData_nRingLocalMax4);
   fChain->SetBranchAddress("HPidTrackCandSim.fData.itsHitData.fInnerMdcChiSquare", HPidTrackCandSim_fData_itsHitData_fInnerMdcChiSquare, &b_HPidTrackCandSim_fData_itsHitData_fInnerMdcChiSquare);
   fChain->SetBranchAddress("HPidTrackCandSim.fData.itsHitData.fInnerMdcdEdx", HPidTrackCandSim_fData_itsHitData_fInnerMdcdEdx, &b_HPidTrackCandSim_fData_itsHitData_fInnerMdcdEdx);
   fChain->SetBranchAddress("HPidTrackCandSim.fData.itsHitData.fInnerMdcdEdxSigma", HPidTrackCandSim_fData_itsHitData_fInnerMdcdEdxSigma, &b_HPidTrackCandSim_fData_itsHitData_fInnerMdcdEdxSigma);
   fChain->SetBranchAddress("HPidTrackCandSim.fData.itsHitData.fMdcRCoord", HPidTrackCandSim_fData_itsHitData_fMdcRCoord, &b_HPidTrackCandSim_fData_itsHitData_fMdcRCoord);
   fChain->SetBranchAddress("HPidTrackCandSim.fData.itsHitData.fMdcZCoord", HPidTrackCandSim_fData_itsHitData_fMdcZCoord, &b_HPidTrackCandSim_fData_itsHitData_fMdcZCoord);
   fChain->SetBranchAddress("HPidTrackCandSim.fData.itsHitData.fMdcTheta", HPidTrackCandSim_fData_itsHitData_fMdcTheta, &b_HPidTrackCandSim_fData_itsHitData_fMdcTheta);
   fChain->SetBranchAddress("HPidTrackCandSim.fData.itsHitData.fMdcPhi", HPidTrackCandSim_fData_itsHitData_fMdcPhi, &b_HPidTrackCandSim_fData_itsHitData_fMdcPhi);
   fChain->SetBranchAddress("HPidTrackCandSim.fData.itsHitData.fOuterMdcChiSquare", HPidTrackCandSim_fData_itsHitData_fOuterMdcChiSquare, &b_HPidTrackCandSim_fData_itsHitData_fOuterMdcChiSquare);
   fChain->SetBranchAddress("HPidTrackCandSim.fData.itsHitData.fOuterMdcdEdx", HPidTrackCandSim_fData_itsHitData_fOuterMdcdEdx, &b_HPidTrackCandSim_fData_itsHitData_fOuterMdcdEdx);
   fChain->SetBranchAddress("HPidTrackCandSim.fData.itsHitData.fOuterMdcdEdxSigma", HPidTrackCandSim_fData_itsHitData_fOuterMdcdEdxSigma, &b_HPidTrackCandSim_fData_itsHitData_fOuterMdcdEdxSigma);
   fChain->SetBranchAddress("HPidTrackCandSim.fData.itsHitData.fCombinedMdcdEdx", HPidTrackCandSim_fData_itsHitData_fCombinedMdcdEdx, &b_HPidTrackCandSim_fData_itsHitData_fCombinedMdcdEdx);
   fChain->SetBranchAddress("HPidTrackCandSim.fData.itsHitData.fCombinedMdcdEdxSigma", HPidTrackCandSim_fData_itsHitData_fCombinedMdcdEdxSigma, &b_HPidTrackCandSim_fData_itsHitData_fCombinedMdcdEdxSigma);
   fChain->SetBranchAddress("HPidTrackCandSim.fData.itsHitData.iIPURingQuality", HPidTrackCandSim_fData_itsHitData_iIPURingQuality, &b_HPidTrackCandSim_fData_itsHitData_iIPURingQuality);
   fChain->SetBranchAddress("HPidTrackCandSim.fData.itsHitData.iIPUVetoQuality", HPidTrackCandSim_fData_itsHitData_iIPUVetoQuality, &b_HPidTrackCandSim_fData_itsHitData_iIPUVetoQuality);
   fChain->SetBranchAddress("HPidTrackCandSim.fData.itsHitData.fShowerSum[3]", HPidTrackCandSim_fData_itsHitData_fShowerSum, &b_HPidTrackCandSim_fData_itsHitData_fShowerSum);
   fChain->SetBranchAddress("HPidTrackCandSim.fData.itsHitData.nShowerClS[3]", HPidTrackCandSim_fData_itsHitData_nShowerClS, &b_HPidTrackCandSim_fData_itsHitData_nShowerClS);
   fChain->SetBranchAddress("HPidTrackCandSim.fData.itsHitData.nShowerRow", HPidTrackCandSim_fData_itsHitData_nShowerRow, &b_HPidTrackCandSim_fData_itsHitData_nShowerRow);
   fChain->SetBranchAddress("HPidTrackCandSim.fData.itsHitData.nShowerCol", HPidTrackCandSim_fData_itsHitData_nShowerCol, &b_HPidTrackCandSim_fData_itsHitData_nShowerCol);
   fChain->SetBranchAddress("HPidTrackCandSim.fData.itsHitData.fShowerTimeOfFlight", HPidTrackCandSim_fData_itsHitData_fShowerTimeOfFlight, &b_HPidTrackCandSim_fData_itsHitData_fShowerTimeOfFlight);
   fChain->SetBranchAddress("HPidTrackCandSim.fData.itsHitData.fMetaLocalX", HPidTrackCandSim_fData_itsHitData_fMetaLocalX, &b_HPidTrackCandSim_fData_itsHitData_fMetaLocalX);
   fChain->SetBranchAddress("HPidTrackCandSim.fData.itsHitData.fMetaLocalY", HPidTrackCandSim_fData_itsHitData_fMetaLocalY, &b_HPidTrackCandSim_fData_itsHitData_fMetaLocalY);
   fChain->SetBranchAddress("HPidTrackCandSim.fData.itsHitData.fTOFTimeOfFlight", HPidTrackCandSim_fData_itsHitData_fTOFTimeOfFlight, &b_HPidTrackCandSim_fData_itsHitData_fTOFTimeOfFlight);
   fChain->SetBranchAddress("HPidTrackCandSim.fData.itsHitData.fTOFLeftAmplitude", HPidTrackCandSim_fData_itsHitData_fTOFLeftAmplitude, &b_HPidTrackCandSim_fData_itsHitData_fTOFLeftAmplitude);
   fChain->SetBranchAddress("HPidTrackCandSim.fData.itsHitData.fTOFRightAmplitude", HPidTrackCandSim_fData_itsHitData_fTOFRightAmplitude, &b_HPidTrackCandSim_fData_itsHitData_fTOFRightAmplitude);
   fChain->SetBranchAddress("HPidTrackCandSim.fData.itsHitData.fTofEloss", HPidTrackCandSim_fData_itsHitData_fTofEloss, &b_HPidTrackCandSim_fData_itsHitData_fTofEloss);
   fChain->SetBranchAddress("HPidTrackCandSim.fData.itsHitData.iTofinoMult", HPidTrackCandSim_fData_itsHitData_iTofinoMult, &b_HPidTrackCandSim_fData_itsHitData_iTofinoMult);
   fChain->SetBranchAddress("HPidTrackCandSim.fData.itsHitData.nTofClsSize", HPidTrackCandSim_fData_itsHitData_nTofClsSize, &b_HPidTrackCandSim_fData_itsHitData_nTofClsSize);
   fChain->SetBranchAddress("HPidTrackCandSim.fData.itsHitData.nMetaCell", HPidTrackCandSim_fData_itsHitData_nMetaCell, &b_HPidTrackCandSim_fData_itsHitData_nMetaCell);
   fChain->SetBranchAddress("HPidTrackCandSim.fData.itsHitData.nTofCell", HPidTrackCandSim_fData_itsHitData_nTofCell, &b_HPidTrackCandSim_fData_itsHitData_nTofCell);
   fChain->SetBranchAddress("HPidTrackCandSim.fData.itsHitData.nTofModule", HPidTrackCandSim_fData_itsHitData_nTofModule, &b_HPidTrackCandSim_fData_itsHitData_nTofModule);
   fChain->SetBranchAddress("HPidTrackCandSim.fData.itsHitData.iIndRICH", HPidTrackCandSim_fData_itsHitData_iIndRICH, &b_HPidTrackCandSim_fData_itsHitData_iIndRICH);
   fChain->SetBranchAddress("HPidTrackCandSim.fData.itsHitData.iIndRICHIPU", HPidTrackCandSim_fData_itsHitData_iIndRICHIPU, &b_HPidTrackCandSim_fData_itsHitData_iIndRICHIPU);
   fChain->SetBranchAddress("HPidTrackCandSim.fData.itsHitData.iIndInnerSeg", HPidTrackCandSim_fData_itsHitData_iIndInnerSeg, &b_HPidTrackCandSim_fData_itsHitData_iIndInnerSeg);
   fChain->SetBranchAddress("HPidTrackCandSim.fData.itsHitData.iIndOuterSeg", HPidTrackCandSim_fData_itsHitData_iIndOuterSeg, &b_HPidTrackCandSim_fData_itsHitData_iIndOuterSeg);
   fChain->SetBranchAddress("HPidTrackCandSim.fData.itsHitData.iIndTOF", HPidTrackCandSim_fData_itsHitData_iIndTOF, &b_HPidTrackCandSim_fData_itsHitData_iIndTOF);
   fChain->SetBranchAddress("HPidTrackCandSim.fData.itsHitData.iIndShower", HPidTrackCandSim_fData_itsHitData_iIndShower, &b_HPidTrackCandSim_fData_itsHitData_iIndShower);
   fChain->SetBranchAddress("HPidTrackCandSim.fData.itsHitData.iIndClusInf0", HPidTrackCandSim_fData_itsHitData_iIndClusInf0, &b_HPidTrackCandSim_fData_itsHitData_iIndClusInf0);
   fChain->SetBranchAddress("HPidTrackCandSim.fData.itsHitData.iIndClusInf1", HPidTrackCandSim_fData_itsHitData_iIndClusInf1, &b_HPidTrackCandSim_fData_itsHitData_iIndClusInf1);
   fChain->SetBranchAddress("HPidTrackCandSim.fData.itsHitData.iIndClusInf2", HPidTrackCandSim_fData_itsHitData_iIndClusInf2, &b_HPidTrackCandSim_fData_itsHitData_iIndClusInf2);
   fChain->SetBranchAddress("HPidTrackCandSim.fData.itsHitData.iIndClusInf3", HPidTrackCandSim_fData_itsHitData_iIndClusInf3, &b_HPidTrackCandSim_fData_itsHitData_iIndClusInf3);
   fChain->SetBranchAddress("HPidTrackCandSim.fData.itsHitData.iIndMatch", HPidTrackCandSim_fData_itsHitData_iIndMatch, &b_HPidTrackCandSim_fData_itsHitData_iIndMatch);
   fChain->SetBranchAddress("HPidTrackCandSim.fData.itsHitData.fDistanceToVertex[10]", HPidTrackCandSim_fData_itsHitData_fDistanceToVertex, &b_HPidTrackCandSim_fData_itsHitData_fDistanceToVertex);
   fChain->SetBranchAddress("HPidTrackCandSim.fData.itsHitData.hasRingCorrelation[10]", HPidTrackCandSim_fData_itsHitData_hasRingCorrelation, &b_HPidTrackCandSim_fData_itsHitData_hasRingCorrelation);
   fChain->SetBranchAddress("HPidTrackCandSim.fData.itsHitData.hasMetaTrackCorrelation[10]", HPidTrackCandSim_fData_itsHitData_hasMetaTrackCorrelation, &b_HPidTrackCandSim_fData_itsHitData_hasMetaTrackCorrelation);
   fChain->SetBranchAddress("HPidTrackCandSim.fData.itsTrackData.fUniqueID", HPidTrackCandSim_fData_itsTrackData_fUniqueID, &b_HPidTrackCandSim_fData_itsTrackData_fUniqueID);
   fChain->SetBranchAddress("HPidTrackCandSim.fData.itsTrackData.fBits", HPidTrackCandSim_fData_itsTrackData_fBits, &b_HPidTrackCandSim_fData_itsTrackData_fBits);
   fChain->SetBranchAddress("HPidTrackCandSim.fData.itsTrackData.cov[10]", HPidTrackCandSim_fData_itsTrackData_cov, &b_HPidTrackCandSim_fData_itsTrackData_cov);
   fChain->SetBranchAddress("HPidTrackCandSim.fData.itsTrackData.nBestMomAlg", HPidTrackCandSim_fData_itsTrackData_nBestMomAlg, &b_HPidTrackCandSim_fData_itsTrackData_nBestMomAlg);
   fChain->SetBranchAddress("HPidTrackCandSim.fData.itsTrackData.nRKTrackInd", HPidTrackCandSim_fData_itsTrackData_nRKTrackInd, &b_HPidTrackCandSim_fData_itsTrackData_nRKTrackInd);
   fChain->SetBranchAddress("HPidTrackCandSim.fData.itsTrackData.nKickTrackInd", HPidTrackCandSim_fData_itsTrackData_nKickTrackInd, &b_HPidTrackCandSim_fData_itsTrackData_nKickTrackInd);
   fChain->SetBranchAddress("HPidTrackCandSim.fData.itsTrackData.nKickTrack123Ind", HPidTrackCandSim_fData_itsTrackData_nKickTrack123Ind, &b_HPidTrackCandSim_fData_itsTrackData_nKickTrack123Ind);
   fChain->SetBranchAddress("HPidTrackCandSim.fData.itsTrackData.nRefTrackInd", HPidTrackCandSim_fData_itsTrackData_nRefTrackInd, &b_HPidTrackCandSim_fData_itsTrackData_nRefTrackInd);
   fChain->SetBranchAddress("HPidTrackCandSim.fData.itsTrackData.nSplineTrackInd", HPidTrackCandSim_fData_itsTrackData_nSplineTrackInd, &b_HPidTrackCandSim_fData_itsTrackData_nSplineTrackInd);
   fChain->SetBranchAddress("HPidTrackCandSim.fData.itsTrackData.fMetaMatchingQuality", HPidTrackCandSim_fData_itsTrackData_fMetaMatchingQuality, &b_HPidTrackCandSim_fData_itsTrackData_fMetaMatchingQuality);
   fChain->SetBranchAddress("HPidTrackCandSim.fData.itsTrackData.fRKRichMatchingQuality", HPidTrackCandSim_fData_itsTrackData_fRKRichMatchingQuality, &b_HPidTrackCandSim_fData_itsTrackData_fRKRichMatchingQuality);
   fChain->SetBranchAddress("HPidTrackCandSim.fData.itsTrackData.dxRkMeta", HPidTrackCandSim_fData_itsTrackData_dxRkMeta, &b_HPidTrackCandSim_fData_itsTrackData_dxRkMeta);
   fChain->SetBranchAddress("HPidTrackCandSim.fData.itsTrackData.dyRkMeta", HPidTrackCandSim_fData_itsTrackData_dyRkMeta, &b_HPidTrackCandSim_fData_itsTrackData_dyRkMeta);
   fChain->SetBranchAddress("HPidTrackCandSim.fData.itsTrackData.dzRkMeta", HPidTrackCandSim_fData_itsTrackData_dzRkMeta, &b_HPidTrackCandSim_fData_itsTrackData_dzRkMeta);
   fChain->SetBranchAddress("HPidTrackCandSim.fData.itsTrackData.dxMdcMeta", HPidTrackCandSim_fData_itsTrackData_dxMdcMeta, &b_HPidTrackCandSim_fData_itsTrackData_dxMdcMeta);
   fChain->SetBranchAddress("HPidTrackCandSim.fData.itsTrackData.dyMdcMeta", HPidTrackCandSim_fData_itsTrackData_dyMdcMeta, &b_HPidTrackCandSim_fData_itsTrackData_dyMdcMeta);
   fChain->SetBranchAddress("HPidTrackCandSim.fData.itsTrackData.dzMdcMeta", HPidTrackCandSim_fData_itsTrackData_dzMdcMeta, &b_HPidTrackCandSim_fData_itsTrackData_dzMdcMeta);
   fChain->SetBranchAddress("HPidTrackCandSim.fData.itsTrackData.xMeta", HPidTrackCandSim_fData_itsTrackData_xMeta, &b_HPidTrackCandSim_fData_itsTrackData_xMeta);
   fChain->SetBranchAddress("HPidTrackCandSim.fData.itsTrackData.yMeta", HPidTrackCandSim_fData_itsTrackData_yMeta, &b_HPidTrackCandSim_fData_itsTrackData_yMeta);
   fChain->SetBranchAddress("HPidTrackCandSim.fData.itsTrackData.zMeta", HPidTrackCandSim_fData_itsTrackData_zMeta, &b_HPidTrackCandSim_fData_itsTrackData_zMeta);
   fChain->SetBranchAddress("HPidTrackCandSim.fData.itsTrackData.errXMeta", HPidTrackCandSim_fData_itsTrackData_errXMeta, &b_HPidTrackCandSim_fData_itsTrackData_errXMeta);
   fChain->SetBranchAddress("HPidTrackCandSim.fData.itsTrackData.errYMeta", HPidTrackCandSim_fData_itsTrackData_errYMeta, &b_HPidTrackCandSim_fData_itsTrackData_errYMeta);
   fChain->SetBranchAddress("HPidTrackCandSim.fData.itsTrackData.errZMeta", HPidTrackCandSim_fData_itsTrackData_errZMeta, &b_HPidTrackCandSim_fData_itsTrackData_errZMeta);
   fChain->SetBranchAddress("HPidTrackCandSim.fData.itsTrackData.nCloseTracklets", HPidTrackCandSim_fData_itsTrackData_nCloseTracklets, &b_HPidTrackCandSim_fData_itsTrackData_nCloseTracklets);
   fChain->SetBranchAddress("HPidTrackCandSim.fData.itsTrackData.fPull", HPidTrackCandSim_fData_itsTrackData_fPull, &b_HPidTrackCandSim_fData_itsTrackData_fPull);
   fChain->SetBranchAddress("HPidTrackCandSim.fData.itsTrackData.fSplineChiSquare", HPidTrackCandSim_fData_itsTrackData_fSplineChiSquare, &b_HPidTrackCandSim_fData_itsTrackData_fSplineChiSquare);
   fChain->SetBranchAddress("HPidTrackCandSim.fData.itsTrackData.fRKChiSquare", HPidTrackCandSim_fData_itsTrackData_fRKChiSquare, &b_HPidTrackCandSim_fData_itsTrackData_fRKChiSquare);
   fChain->SetBranchAddress("HPidTrackCandSim.fData.itsTrackData.iIndexClosestTracklet", HPidTrackCandSim_fData_itsTrackData_iIndexClosestTracklet, &b_HPidTrackCandSim_fData_itsTrackData_iIndexClosestTracklet);
   fChain->SetBranchAddress("HPidTrackCandSim.fData.itsTrackData.aTrackletClusInf0", HPidTrackCandSim_fData_itsTrackData_aTrackletClusInf0, &b_HPidTrackCandSim_fData_itsTrackData_aTrackletClusInf0);
   fChain->SetBranchAddress("HPidTrackCandSim.fData.itsTrackData.aTrackletClusInf1", HPidTrackCandSim_fData_itsTrackData_aTrackletClusInf1, &b_HPidTrackCandSim_fData_itsTrackData_aTrackletClusInf1);
   fChain->SetBranchAddress("HPidTrackCandSim.fData.itsTrackData.aTrackletDistances", HPidTrackCandSim_fData_itsTrackData_aTrackletDistances, &b_HPidTrackCandSim_fData_itsTrackData_aTrackletDistances);
   fChain->SetBranchAddress("HPidTrackCandSim.fData.itsTrackData.qIOMatching[10]", HPidTrackCandSim_fData_itsTrackData_qIOMatching, &b_HPidTrackCandSim_fData_itsTrackData_qIOMatching);
   fChain->SetBranchAddress("HPidTrackCandSim.fData.itsTrackData.nPolarity[10]", HPidTrackCandSim_fData_itsTrackData_nPolarity, &b_HPidTrackCandSim_fData_itsTrackData_nPolarity);
   fChain->SetBranchAddress("HPidTrackCandSim.fData.itsTrackData.fMomenta[10]", HPidTrackCandSim_fData_itsTrackData_fMomenta, &b_HPidTrackCandSim_fData_itsTrackData_fMomenta);
   fChain->SetBranchAddress("HPidTrackCandSim.fData.itsTrackData.fMomError[10]", HPidTrackCandSim_fData_itsTrackData_fMomError, &b_HPidTrackCandSim_fData_itsTrackData_fMomError);
   fChain->SetBranchAddress("HPidTrackCandSim.fData.itsTrackData.fTrackR[10]", HPidTrackCandSim_fData_itsTrackData_fTrackR, &b_HPidTrackCandSim_fData_itsTrackData_fTrackR);
   fChain->SetBranchAddress("HPidTrackCandSim.fData.itsTrackData.fTrackZ[10]", HPidTrackCandSim_fData_itsTrackData_fTrackZ, &b_HPidTrackCandSim_fData_itsTrackData_fTrackZ);
   fChain->SetBranchAddress("HPidTrackCandSim.fData.itsTrackData.fRKPhi", HPidTrackCandSim_fData_itsTrackData_fRKPhi, &b_HPidTrackCandSim_fData_itsTrackData_fRKPhi);
   fChain->SetBranchAddress("HPidTrackCandSim.fData.itsTrackData.fRKTheta", HPidTrackCandSim_fData_itsTrackData_fRKTheta, &b_HPidTrackCandSim_fData_itsTrackData_fRKTheta);
   fChain->SetBranchAddress("HPidTrackCandSim.fData.itsTrackData.fCorrectedEloss[10]", HPidTrackCandSim_fData_itsTrackData_fCorrectedEloss, &b_HPidTrackCandSim_fData_itsTrackData_fCorrectedEloss);
   fChain->SetBranchAddress("HPidTrackCandSim.fData.itsTrackData.fCorrectedBeta[10]", HPidTrackCandSim_fData_itsTrackData_fCorrectedBeta, &b_HPidTrackCandSim_fData_itsTrackData_fCorrectedBeta);
   fChain->SetBranchAddress("HPidTrackCandSim.fData.itsTrackData.fPathLength[10]", HPidTrackCandSim_fData_itsTrackData_fPathLength, &b_HPidTrackCandSim_fData_itsTrackData_fPathLength);
   fChain->SetBranchAddress("HPidTrackCandSim.fData.itsTrackData.fMassSquared[10]", HPidTrackCandSim_fData_itsTrackData_fMassSquared, &b_HPidTrackCandSim_fData_itsTrackData_fMassSquared);
   fChain->SetBranchAddress("HPidTrackCandSim.fData.itsTrackData.bIsAccepted[10]", HPidTrackCandSim_fData_itsTrackData_bIsAccepted, &b_HPidTrackCandSim_fData_itsTrackData_bIsAccepted);
   fChain->SetBranchAddress("HPidTrackCandSim.fData.itsTrackData.nTofRecFlag[10]", HPidTrackCandSim_fData_itsTrackData_nTofRecFlag, &b_HPidTrackCandSim_fData_itsTrackData_nTofRecFlag);
   fChain->SetBranchAddress("HPidTrackCandSim.fData.itsTrackData.fTofRecTof[10]", HPidTrackCandSim_fData_itsTrackData_fTofRecTof, &b_HPidTrackCandSim_fData_itsTrackData_fTofRecTof);
   fChain->SetBranchAddress("HPidTrackCandSim.fData.itsTrackData.fTofRecBeta[10]", HPidTrackCandSim_fData_itsTrackData_fTofRecBeta, &b_HPidTrackCandSim_fData_itsTrackData_fTofRecBeta);
   fChain->SetBranchAddress("HPidTrackCandSim.fData.itsTrackData.fTofRecMassSquared[10]", HPidTrackCandSim_fData_itsTrackData_fTofRecMassSquared, &b_HPidTrackCandSim_fData_itsTrackData_fTofRecMassSquared);
   fChain->SetBranchAddress("HPidTrackCandSim.fData.flags", HPidTrackCandSim_fData_flags, &b_HPidTrackCandSim_fData_flags);
   fChain->SetBranchAddress("HPidTrackCandSim.fData.itsGeantTrackSet.fUniqueID", HPidTrackCandSim_fData_itsGeantTrackSet_fUniqueID, &b_HPidTrackCandSim_fData_itsGeantTrackSet_fUniqueID);
   fChain->SetBranchAddress("HPidTrackCandSim.fData.itsGeantTrackSet.fBits", HPidTrackCandSim_fData_itsGeantTrackSet_fBits, &b_HPidTrackCandSim_fData_itsGeantTrackSet_fBits);
   fChain->SetBranchAddress("HPidTrackCandSim.fData.itsGeantTrackSet.isSorted", HPidTrackCandSim_fData_itsGeantTrackSet_isSorted, &b_HPidTrackCandSim_fData_itsGeantTrackSet_isSorted);
   fChain->SetBranchAddress("HPidTrackCandSim.fData.itsGeantTrackSet.sNCorrTrackIds", HPidTrackCandSim_fData_itsGeantTrackSet_sNCorrTrackIds, &b_HPidTrackCandSim_fData_itsGeantTrackSet_sNCorrTrackIds);
   fChain->SetBranchAddress("HPidTrackCandSim.fData.itsGeantTrackSet.nRICHTracks", HPidTrackCandSim_fData_itsGeantTrackSet_nRICHTracks, &b_HPidTrackCandSim_fData_itsGeantTrackSet_nRICHTracks);
   fChain->SetBranchAddress("HPidTrackCandSim.fData.itsGeantTrackSet.nRICHIPUTracks", HPidTrackCandSim_fData_itsGeantTrackSet_nRICHIPUTracks, &b_HPidTrackCandSim_fData_itsGeantTrackSet_nRICHIPUTracks);
   fChain->SetBranchAddress("HPidTrackCandSim.fData.itsGeantTrackSet.nInnerMdcTracks", HPidTrackCandSim_fData_itsGeantTrackSet_nInnerMdcTracks, &b_HPidTrackCandSim_fData_itsGeantTrackSet_nInnerMdcTracks);
   fChain->SetBranchAddress("HPidTrackCandSim.fData.itsGeantTrackSet.nOuterMdcTracks", HPidTrackCandSim_fData_itsGeantTrackSet_nOuterMdcTracks, &b_HPidTrackCandSim_fData_itsGeantTrackSet_nOuterMdcTracks);
   fChain->SetBranchAddress("HPidTrackCandSim.fData.itsGeantTrackSet.nShowerTracks", HPidTrackCandSim_fData_itsGeantTrackSet_nShowerTracks, &b_HPidTrackCandSim_fData_itsGeantTrackSet_nShowerTracks);
   fChain->SetBranchAddress("HPidTrackCandSim.fData.itsGeantTrackSet.nTOFTracks", HPidTrackCandSim_fData_itsGeantTrackSet_nTOFTracks, &b_HPidTrackCandSim_fData_itsGeantTrackSet_nTOFTracks);
   fChain->SetBranchAddress("HPidTrackCandSim.fData.itsGeantTrackSet.bIsLepFromPrimary", HPidTrackCandSim_fData_itsGeantTrackSet_bIsLepFromPrimary, &b_HPidTrackCandSim_fData_itsGeantTrackSet_bIsLepFromPrimary);
   fChain->SetBranchAddress("HPidTrackCandSim.fData.itsGeantTrackSet.correlatedTrackIds", HPidTrackCandSim_fData_itsGeantTrackSet_correlatedTrackIds, &b_HPidTrackCandSim_fData_itsGeantTrackSet_correlatedTrackIds);
   fChain->SetBranchAddress("HPidTrackCandSim.fData.itsGeantTrackSet.correlationFlags", HPidTrackCandSim_fData_itsGeantTrackSet_correlationFlags, &b_HPidTrackCandSim_fData_itsGeantTrackSet_correlationFlags);
   fChain->SetBranchAddress("HPidTrackCandSim.fData.itsGeantTrackSet.ProcessIds", HPidTrackCandSim_fData_itsGeantTrackSet_ProcessIds, &b_HPidTrackCandSim_fData_itsGeantTrackSet_ProcessIds);
   fChain->SetBranchAddress("HPidTrackCandSim.fData.itsGeantTrackSet.ParentIds", HPidTrackCandSim_fData_itsGeantTrackSet_ParentIds, &b_HPidTrackCandSim_fData_itsGeantTrackSet_ParentIds);
   fChain->SetBranchAddress("HPidTrackCandSim.fData.itsGeantTrackSet.Parents", HPidTrackCandSim_fData_itsGeantTrackSet_Parents, &b_HPidTrackCandSim_fData_itsGeantTrackSet_Parents);
   fChain->SetBranchAddress("HPidTrackCandSim.fData.itsGeantTrackSet.GrandParentIds", HPidTrackCandSim_fData_itsGeantTrackSet_GrandParentIds, &b_HPidTrackCandSim_fData_itsGeantTrackSet_GrandParentIds);
   fChain->SetBranchAddress("HPidTrackCandSim.fData.itsGeantTrackSet.GrandParents", HPidTrackCandSim_fData_itsGeantTrackSet_GrandParents, &b_HPidTrackCandSim_fData_itsGeantTrackSet_GrandParents);
   fChain->SetBranchAddress("HPidTrackCandSim.fData.itsGeantTrackSet.GenInfo", HPidTrackCandSim_fData_itsGeantTrackSet_GenInfo, &b_HPidTrackCandSim_fData_itsGeantTrackSet_GenInfo);
   fChain->SetBranchAddress("HPidTrackCandSim.fData.itsGeantTrackSet.GenInfo1", HPidTrackCandSim_fData_itsGeantTrackSet_GenInfo1, &b_HPidTrackCandSim_fData_itsGeantTrackSet_GenInfo1);
   fChain->SetBranchAddress("HPidTrackCandSim.fData.itsGeantTrackSet.GenInfo2", HPidTrackCandSim_fData_itsGeantTrackSet_GenInfo2, &b_HPidTrackCandSim_fData_itsGeantTrackSet_GenInfo2);
   fChain->SetBranchAddress("HPidTrackCandSim.fData.itsGeantTrackSet.GenWeight", HPidTrackCandSim_fData_itsGeantTrackSet_GenWeight, &b_HPidTrackCandSim_fData_itsGeantTrackSet_GenWeight);
   fChain->SetBranchAddress("HPidTrackCandSim.fData.itsGeantTrackSet.VertexX", HPidTrackCandSim_fData_itsGeantTrackSet_VertexX, &b_HPidTrackCandSim_fData_itsGeantTrackSet_VertexX);
   fChain->SetBranchAddress("HPidTrackCandSim.fData.itsGeantTrackSet.VertexY", HPidTrackCandSim_fData_itsGeantTrackSet_VertexY, &b_HPidTrackCandSim_fData_itsGeantTrackSet_VertexY);
   fChain->SetBranchAddress("HPidTrackCandSim.fData.itsGeantTrackSet.VertexZ", HPidTrackCandSim_fData_itsGeantTrackSet_VertexZ, &b_HPidTrackCandSim_fData_itsGeantTrackSet_VertexZ);
   fChain->SetBranchAddress("HPidTrackCandSim.fData.itsGeantTrackSet.GeantPIDs", HPidTrackCandSim_fData_itsGeantTrackSet_GeantPIDs, &b_HPidTrackCandSim_fData_itsGeantTrackSet_GeantPIDs);
   fChain->SetBranchAddress("HPidTrackCandSim.fData.itsGeantTrackSet.MediumIds", HPidTrackCandSim_fData_itsGeantTrackSet_MediumIds, &b_HPidTrackCandSim_fData_itsGeantTrackSet_MediumIds);
   fChain->SetBranchAddress("HPidTrackCandSim.fData.itsGeantTrackSet.GeantMomX", HPidTrackCandSim_fData_itsGeantTrackSet_GeantMomX, &b_HPidTrackCandSim_fData_itsGeantTrackSet_GeantMomX);
   fChain->SetBranchAddress("HPidTrackCandSim.fData.itsGeantTrackSet.GeantMomY", HPidTrackCandSim_fData_itsGeantTrackSet_GeantMomY, &b_HPidTrackCandSim_fData_itsGeantTrackSet_GeantMomY);
   fChain->SetBranchAddress("HPidTrackCandSim.fData.itsGeantTrackSet.GeantMomZ", HPidTrackCandSim_fData_itsGeantTrackSet_GeantMomZ, &b_HPidTrackCandSim_fData_itsGeantTrackSet_GeantMomZ);
   fChain->SetBranchAddress("HPidTrackCandSim.fData.itsGeantTrackSet.ShowerWeights", HPidTrackCandSim_fData_itsGeantTrackSet_ShowerWeights, &b_HPidTrackCandSim_fData_itsGeantTrackSet_ShowerWeights);
   fChain->SetBranchAddress("HPidTrackCandSim.fData.itsGeantTrackSet.TOFWeights", HPidTrackCandSim_fData_itsGeantTrackSet_TOFWeights, &b_HPidTrackCandSim_fData_itsGeantTrackSet_TOFWeights);
   fChain->SetBranchAddress("HPidTrackCandSim.fData.itsGeantTrackSet.RICHWeights", HPidTrackCandSim_fData_itsGeantTrackSet_RICHWeights, &b_HPidTrackCandSim_fData_itsGeantTrackSet_RICHWeights);
   fChain->SetBranchAddress("HPidTrackCandSim.fData.itsGeantTrackSet.RICHIPUWeights", HPidTrackCandSim_fData_itsGeantTrackSet_RICHIPUWeights, &b_HPidTrackCandSim_fData_itsGeantTrackSet_RICHIPUWeights);
   fChain->SetBranchAddress("HPidTrackCandSim.fData.itsGeantTrackSet.InnerMDCWeights", HPidTrackCandSim_fData_itsGeantTrackSet_InnerMDCWeights, &b_HPidTrackCandSim_fData_itsGeantTrackSet_InnerMDCWeights);
   fChain->SetBranchAddress("HPidTrackCandSim.fData.itsGeantTrackSet.OuterMDCWeights", HPidTrackCandSim_fData_itsGeantTrackSet_OuterMDCWeights, &b_HPidTrackCandSim_fData_itsGeantTrackSet_OuterMDCWeights);
   fChain->SetBranchAddress("HPidTrackCandSim.fData.itsGeantTrackSet.hasHitInShower", HPidTrackCandSim_fData_itsGeantTrackSet_hasHitInShower, &b_HPidTrackCandSim_fData_itsGeantTrackSet_hasHitInShower);
   fChain->SetBranchAddress("HPidTrackCandSim.fData.itsGeantTrackSet.hasHitInTOF", HPidTrackCandSim_fData_itsGeantTrackSet_hasHitInTOF, &b_HPidTrackCandSim_fData_itsGeantTrackSet_hasHitInTOF);
   fChain->SetBranchAddress("HPidTrackCandSim.fNDataObjs", &HPidTrackCandSim_fNDataObjs, &b_HPidTrackCandSim_fNDataObjs);
   fChain->SetBranchAddress("HPidTrackCandSim.hasDynamicObjects", &HPidTrackCandSim_hasDynamicObjects, &b_HPidTrackCandSim_hasDynamicObjects);
   fChain->SetBranchAddress("HPidCandidate.HCategory.fUniqueID", &HPidCandidate_HCategory_fUniqueID, &b_HPidCandidate_HCategory_fUniqueID);
   fChain->SetBranchAddress("HPidCandidate.HCategory.fBits", &HPidCandidate_HCategory_fBits, &b_HPidCandidate_HCategory_fBits);
   fChain->SetBranchAddress("HPidCandidate.HCategory.fCat", &HPidCandidate_HCategory_fCat, &b_HPidCandidate_HCategory_fCat);
   fChain->SetBranchAddress("HPidCandidate.HCategory.fBranchingLevel", &HPidCandidate_HCategory_fBranchingLevel, &b_HPidCandidate_HCategory_fBranchingLevel);
   fChain->SetBranchAddress("HPidCandidate.fData", &HPidCandidate_fData_, &b_HPidCandidate_fData_);
   fChain->SetBranchAddress("HPidCandidate.fData.fUniqueID", HPidCandidate_fData_fUniqueID, &b_HPidCandidate_fData_fUniqueID);
   fChain->SetBranchAddress("HPidCandidate.fData.fBits", HPidCandidate_fData_fBits, &b_HPidCandidate_fData_fBits);
   fChain->SetBranchAddress("HPidCandidate.fData.iTrackCandIndex", HPidCandidate_fData_iTrackCandIndex, &b_HPidCandidate_fData_iTrackCandIndex);
   fChain->SetBranchAddress("HPidCandidate.fData.NUM_ALGORITHMS", HPidCandidate_fData_NUM_ALGORITHMS, &b_HPidCandidate_fData_NUM_ALGORITHMS);
   fChain->SetBranchAddress("HPidCandidate.fData.NUM_PARTICLES", HPidCandidate_fData_NUM_PARTICLES, &b_HPidCandidate_fData_NUM_PARTICLES);
   fChain->SetBranchAddress("HPidCandidate.fData.NUM_VALUES", HPidCandidate_fData_NUM_VALUES, &b_HPidCandidate_fData_NUM_VALUES);
   fChain->SetBranchAddress("HPidCandidate.fData.aAlgorithms", HPidCandidate_fData_aAlgorithms, &b_HPidCandidate_fData_aAlgorithms);
   fChain->SetBranchAddress("HPidCandidate.fData.aParticles", HPidCandidate_fData_aParticles, &b_HPidCandidate_fData_aParticles);
   fChain->SetBranchAddress("HPidCandidate.fData.aValues", HPidCandidate_fData_aValues, &b_HPidCandidate_fData_aValues);
   fChain->SetBranchAddress("HPidCandidate.fData.nMomAlgIndex", HPidCandidate_fData_nMomAlgIndex, &b_HPidCandidate_fData_nMomAlgIndex);
   fChain->SetBranchAddress("HPidCandidate.fData.flags", HPidCandidate_fData_flags, &b_HPidCandidate_fData_flags);
   fChain->SetBranchAddress("HPidCandidate.fNDataObjs", &HPidCandidate_fNDataObjs, &b_HPidCandidate_fNDataObjs);
   fChain->SetBranchAddress("HPidCandidate.hasDynamicObjects", &HPidCandidate_hasDynamicObjects, &b_HPidCandidate_hasDynamicObjects);
   fChain->SetBranchAddress("HPidParticleSim.HCategory.fUniqueID", &HPidParticleSim_HCategory_fUniqueID, &b_HPidParticleSim_HCategory_fUniqueID);
   fChain->SetBranchAddress("HPidParticleSim.HCategory.fBits", &HPidParticleSim_HCategory_fBits, &b_HPidParticleSim_HCategory_fBits);
   fChain->SetBranchAddress("HPidParticleSim.HCategory.fCat", &HPidParticleSim_HCategory_fCat, &b_HPidParticleSim_HCategory_fCat);
   fChain->SetBranchAddress("HPidParticleSim.HCategory.fBranchingLevel", &HPidParticleSim_HCategory_fBranchingLevel, &b_HPidParticleSim_HCategory_fBranchingLevel);
   fChain->SetBranchAddress("HPidParticleSim.fData", &HPidParticleSim_fData_, &b_HPidParticleSim_fData_);
   fChain->SetBranchAddress("HPidParticleSim.fData.fUniqueID", HPidParticleSim_fData_fUniqueID, &b_HPidParticleSim_fData_fUniqueID);
   fChain->SetBranchAddress("HPidParticleSim.fData.fBits", HPidParticleSim_fData_fBits, &b_HPidParticleSim_fData_fBits);
   fChain->SetBranchAddress("HPidParticleSim.fData.fP.fUniqueID", HPidParticleSim_fData_fP_fUniqueID, &b_HPidParticleSim_fData_fP_fUniqueID);
   fChain->SetBranchAddress("HPidParticleSim.fData.fP.fBits", HPidParticleSim_fData_fP_fBits, &b_HPidParticleSim_fData_fP_fBits);
   fChain->SetBranchAddress("HPidParticleSim.fData.fP.fX", HPidParticleSim_fData_fP_fX, &b_HPidParticleSim_fData_fP_fX);
   fChain->SetBranchAddress("HPidParticleSim.fData.fP.fY", HPidParticleSim_fData_fP_fY, &b_HPidParticleSim_fData_fP_fY);
   fChain->SetBranchAddress("HPidParticleSim.fData.fP.fZ", HPidParticleSim_fData_fP_fZ, &b_HPidParticleSim_fData_fP_fZ);
   fChain->SetBranchAddress("HPidParticleSim.fData.fE", HPidParticleSim_fData_fE, &b_HPidParticleSim_fData_fE);
   fChain->SetBranchAddress("HPidParticleSim.fData.kUsesIdealMass", HPidParticleSim_fData_kUsesIdealMass, &b_HPidParticleSim_fData_kUsesIdealMass);
   fChain->SetBranchAddress("HPidParticleSim.fData.nPossibleSpecies", HPidParticleSim_fData_nPossibleSpecies, &b_HPidParticleSim_fData_nPossibleSpecies);
   fChain->SetBranchAddress("HPidParticleSim.fData.momAlgIndex", HPidParticleSim_fData_momAlgIndex, &b_HPidParticleSim_fData_momAlgIndex);
   fChain->SetBranchAddress("HPidParticleSim.fData.nPidCandidateIndex", HPidParticleSim_fData_nPidCandidateIndex, &b_HPidParticleSim_fData_nPidCandidateIndex);
   fChain->SetBranchAddress("HPidParticleSim.fData.possibleSpecies", HPidParticleSim_fData_possibleSpecies, &b_HPidParticleSim_fData_possibleSpecies);
   fChain->SetBranchAddress("HPidParticleSim.fData.assignedWeights", HPidParticleSim_fData_assignedWeights, &b_HPidParticleSim_fData_assignedWeights);
   fChain->SetBranchAddress("HPidParticleSim.fData.nAssignedPID", HPidParticleSim_fData_nAssignedPID, &b_HPidParticleSim_fData_nAssignedPID);
   fChain->SetBranchAddress("HPidParticleSim.fData.fTestVal", HPidParticleSim_fData_fTestVal, &b_HPidParticleSim_fData_fTestVal);
   fChain->SetBranchAddress("HPidParticleSim.fData.fWeight", HPidParticleSim_fData_fWeight, &b_HPidParticleSim_fData_fWeight);
   fChain->SetBranchAddress("HPidParticleSim.fData.fMomRescal", HPidParticleSim_fData_fMomRescal, &b_HPidParticleSim_fData_fMomRescal);
   fChain->SetBranchAddress("HPidParticleSim.fData.itsHitData.fUniqueID", HPidParticleSim_fData_itsHitData_fUniqueID, &b_HPidParticleSim_fData_itsHitData_fUniqueID);
   fChain->SetBranchAddress("HPidParticleSim.fData.itsHitData.fBits", HPidParticleSim_fData_itsHitData_fBits, &b_HPidParticleSim_fData_itsHitData_fBits);
   fChain->SetBranchAddress("HPidParticleSim.fData.itsHitData.nSector", HPidParticleSim_fData_itsHitData_nSector, &b_HPidParticleSim_fData_itsHitData_nSector);
   fChain->SetBranchAddress("HPidParticleSim.fData.itsHitData.iSystem", HPidParticleSim_fData_itsHitData_iSystem, &b_HPidParticleSim_fData_itsHitData_iSystem);
   fChain->SetBranchAddress("HPidParticleSim.fData.itsHitData.nRingPadNr", HPidParticleSim_fData_itsHitData_nRingPadNr, &b_HPidParticleSim_fData_itsHitData_nRingPadNr);
   fChain->SetBranchAddress("HPidParticleSim.fData.itsHitData.fRingCentroid", HPidParticleSim_fData_itsHitData_fRingCentroid, &b_HPidParticleSim_fData_itsHitData_fRingCentroid);
   fChain->SetBranchAddress("HPidParticleSim.fData.itsHitData.fRichTheta", HPidParticleSim_fData_itsHitData_fRichTheta, &b_HPidParticleSim_fData_itsHitData_fRichTheta);
   fChain->SetBranchAddress("HPidParticleSim.fData.itsHitData.fRichPhi", HPidParticleSim_fData_itsHitData_fRichPhi, &b_HPidParticleSim_fData_itsHitData_fRichPhi);
   fChain->SetBranchAddress("HPidParticleSim.fData.itsHitData.nRingPatMat", HPidParticleSim_fData_itsHitData_nRingPatMat, &b_HPidParticleSim_fData_itsHitData_nRingPatMat);
   fChain->SetBranchAddress("HPidParticleSim.fData.itsHitData.nRingHouTra", HPidParticleSim_fData_itsHitData_nRingHouTra, &b_HPidParticleSim_fData_itsHitData_nRingHouTra);
   fChain->SetBranchAddress("HPidParticleSim.fData.itsHitData.nRingAmplitude", HPidParticleSim_fData_itsHitData_nRingAmplitude, &b_HPidParticleSim_fData_itsHitData_nRingAmplitude);
   fChain->SetBranchAddress("HPidParticleSim.fData.itsHitData.nRingLocalMax4", HPidParticleSim_fData_itsHitData_nRingLocalMax4, &b_HPidParticleSim_fData_itsHitData_nRingLocalMax4);
   fChain->SetBranchAddress("HPidParticleSim.fData.itsHitData.fInnerMdcChiSquare", HPidParticleSim_fData_itsHitData_fInnerMdcChiSquare, &b_HPidParticleSim_fData_itsHitData_fInnerMdcChiSquare);
   fChain->SetBranchAddress("HPidParticleSim.fData.itsHitData.fInnerMdcdEdx", HPidParticleSim_fData_itsHitData_fInnerMdcdEdx, &b_HPidParticleSim_fData_itsHitData_fInnerMdcdEdx);
   fChain->SetBranchAddress("HPidParticleSim.fData.itsHitData.fInnerMdcdEdxSigma", HPidParticleSim_fData_itsHitData_fInnerMdcdEdxSigma, &b_HPidParticleSim_fData_itsHitData_fInnerMdcdEdxSigma);
   fChain->SetBranchAddress("HPidParticleSim.fData.itsHitData.fMdcRCoord", HPidParticleSim_fData_itsHitData_fMdcRCoord, &b_HPidParticleSim_fData_itsHitData_fMdcRCoord);
   fChain->SetBranchAddress("HPidParticleSim.fData.itsHitData.fMdcZCoord", HPidParticleSim_fData_itsHitData_fMdcZCoord, &b_HPidParticleSim_fData_itsHitData_fMdcZCoord);
   fChain->SetBranchAddress("HPidParticleSim.fData.itsHitData.fMdcTheta", HPidParticleSim_fData_itsHitData_fMdcTheta, &b_HPidParticleSim_fData_itsHitData_fMdcTheta);
   fChain->SetBranchAddress("HPidParticleSim.fData.itsHitData.fMdcPhi", HPidParticleSim_fData_itsHitData_fMdcPhi, &b_HPidParticleSim_fData_itsHitData_fMdcPhi);
   fChain->SetBranchAddress("HPidParticleSim.fData.itsHitData.fOuterMdcChiSquare", HPidParticleSim_fData_itsHitData_fOuterMdcChiSquare, &b_HPidParticleSim_fData_itsHitData_fOuterMdcChiSquare);
   fChain->SetBranchAddress("HPidParticleSim.fData.itsHitData.fOuterMdcdEdx", HPidParticleSim_fData_itsHitData_fOuterMdcdEdx, &b_HPidParticleSim_fData_itsHitData_fOuterMdcdEdx);
   fChain->SetBranchAddress("HPidParticleSim.fData.itsHitData.fOuterMdcdEdxSigma", HPidParticleSim_fData_itsHitData_fOuterMdcdEdxSigma, &b_HPidParticleSim_fData_itsHitData_fOuterMdcdEdxSigma);
   fChain->SetBranchAddress("HPidParticleSim.fData.itsHitData.fCombinedMdcdEdx", HPidParticleSim_fData_itsHitData_fCombinedMdcdEdx, &b_HPidParticleSim_fData_itsHitData_fCombinedMdcdEdx);
   fChain->SetBranchAddress("HPidParticleSim.fData.itsHitData.fCombinedMdcdEdxSigma", HPidParticleSim_fData_itsHitData_fCombinedMdcdEdxSigma, &b_HPidParticleSim_fData_itsHitData_fCombinedMdcdEdxSigma);
   fChain->SetBranchAddress("HPidParticleSim.fData.itsHitData.iIPURingQuality", HPidParticleSim_fData_itsHitData_iIPURingQuality, &b_HPidParticleSim_fData_itsHitData_iIPURingQuality);
   fChain->SetBranchAddress("HPidParticleSim.fData.itsHitData.iIPUVetoQuality", HPidParticleSim_fData_itsHitData_iIPUVetoQuality, &b_HPidParticleSim_fData_itsHitData_iIPUVetoQuality);
   fChain->SetBranchAddress("HPidParticleSim.fData.itsHitData.fShowerSum[3]", HPidParticleSim_fData_itsHitData_fShowerSum, &b_HPidParticleSim_fData_itsHitData_fShowerSum);
   fChain->SetBranchAddress("HPidParticleSim.fData.itsHitData.nShowerClS[3]", HPidParticleSim_fData_itsHitData_nShowerClS, &b_HPidParticleSim_fData_itsHitData_nShowerClS);
   fChain->SetBranchAddress("HPidParticleSim.fData.itsHitData.nShowerRow", HPidParticleSim_fData_itsHitData_nShowerRow, &b_HPidParticleSim_fData_itsHitData_nShowerRow);
   fChain->SetBranchAddress("HPidParticleSim.fData.itsHitData.nShowerCol", HPidParticleSim_fData_itsHitData_nShowerCol, &b_HPidParticleSim_fData_itsHitData_nShowerCol);
   fChain->SetBranchAddress("HPidParticleSim.fData.itsHitData.fShowerTimeOfFlight", HPidParticleSim_fData_itsHitData_fShowerTimeOfFlight, &b_HPidParticleSim_fData_itsHitData_fShowerTimeOfFlight);
   fChain->SetBranchAddress("HPidParticleSim.fData.itsHitData.fMetaLocalX", HPidParticleSim_fData_itsHitData_fMetaLocalX, &b_HPidParticleSim_fData_itsHitData_fMetaLocalX);
   fChain->SetBranchAddress("HPidParticleSim.fData.itsHitData.fMetaLocalY", HPidParticleSim_fData_itsHitData_fMetaLocalY, &b_HPidParticleSim_fData_itsHitData_fMetaLocalY);
   fChain->SetBranchAddress("HPidParticleSim.fData.itsHitData.fTOFTimeOfFlight", HPidParticleSim_fData_itsHitData_fTOFTimeOfFlight, &b_HPidParticleSim_fData_itsHitData_fTOFTimeOfFlight);
   fChain->SetBranchAddress("HPidParticleSim.fData.itsHitData.fTOFLeftAmplitude", HPidParticleSim_fData_itsHitData_fTOFLeftAmplitude, &b_HPidParticleSim_fData_itsHitData_fTOFLeftAmplitude);
   fChain->SetBranchAddress("HPidParticleSim.fData.itsHitData.fTOFRightAmplitude", HPidParticleSim_fData_itsHitData_fTOFRightAmplitude, &b_HPidParticleSim_fData_itsHitData_fTOFRightAmplitude);
   fChain->SetBranchAddress("HPidParticleSim.fData.itsHitData.fTofEloss", HPidParticleSim_fData_itsHitData_fTofEloss, &b_HPidParticleSim_fData_itsHitData_fTofEloss);
   fChain->SetBranchAddress("HPidParticleSim.fData.itsHitData.iTofinoMult", HPidParticleSim_fData_itsHitData_iTofinoMult, &b_HPidParticleSim_fData_itsHitData_iTofinoMult);
   fChain->SetBranchAddress("HPidParticleSim.fData.itsHitData.nTofClsSize", HPidParticleSim_fData_itsHitData_nTofClsSize, &b_HPidParticleSim_fData_itsHitData_nTofClsSize);
   fChain->SetBranchAddress("HPidParticleSim.fData.itsHitData.nMetaCell", HPidParticleSim_fData_itsHitData_nMetaCell, &b_HPidParticleSim_fData_itsHitData_nMetaCell);
   fChain->SetBranchAddress("HPidParticleSim.fData.itsHitData.nTofCell", HPidParticleSim_fData_itsHitData_nTofCell, &b_HPidParticleSim_fData_itsHitData_nTofCell);
   fChain->SetBranchAddress("HPidParticleSim.fData.itsHitData.nTofModule", HPidParticleSim_fData_itsHitData_nTofModule, &b_HPidParticleSim_fData_itsHitData_nTofModule);
   fChain->SetBranchAddress("HPidParticleSim.fData.itsHitData.iIndRICH", HPidParticleSim_fData_itsHitData_iIndRICH, &b_HPidParticleSim_fData_itsHitData_iIndRICH);
   fChain->SetBranchAddress("HPidParticleSim.fData.itsHitData.iIndRICHIPU", HPidParticleSim_fData_itsHitData_iIndRICHIPU, &b_HPidParticleSim_fData_itsHitData_iIndRICHIPU);
   fChain->SetBranchAddress("HPidParticleSim.fData.itsHitData.iIndInnerSeg", HPidParticleSim_fData_itsHitData_iIndInnerSeg, &b_HPidParticleSim_fData_itsHitData_iIndInnerSeg);
   fChain->SetBranchAddress("HPidParticleSim.fData.itsHitData.iIndOuterSeg", HPidParticleSim_fData_itsHitData_iIndOuterSeg, &b_HPidParticleSim_fData_itsHitData_iIndOuterSeg);
   fChain->SetBranchAddress("HPidParticleSim.fData.itsHitData.iIndTOF", HPidParticleSim_fData_itsHitData_iIndTOF, &b_HPidParticleSim_fData_itsHitData_iIndTOF);
   fChain->SetBranchAddress("HPidParticleSim.fData.itsHitData.iIndShower", HPidParticleSim_fData_itsHitData_iIndShower, &b_HPidParticleSim_fData_itsHitData_iIndShower);
   fChain->SetBranchAddress("HPidParticleSim.fData.itsHitData.iIndClusInf0", HPidParticleSim_fData_itsHitData_iIndClusInf0, &b_HPidParticleSim_fData_itsHitData_iIndClusInf0);
   fChain->SetBranchAddress("HPidParticleSim.fData.itsHitData.iIndClusInf1", HPidParticleSim_fData_itsHitData_iIndClusInf1, &b_HPidParticleSim_fData_itsHitData_iIndClusInf1);
   fChain->SetBranchAddress("HPidParticleSim.fData.itsHitData.iIndClusInf2", HPidParticleSim_fData_itsHitData_iIndClusInf2, &b_HPidParticleSim_fData_itsHitData_iIndClusInf2);
   fChain->SetBranchAddress("HPidParticleSim.fData.itsHitData.iIndClusInf3", HPidParticleSim_fData_itsHitData_iIndClusInf3, &b_HPidParticleSim_fData_itsHitData_iIndClusInf3);
   fChain->SetBranchAddress("HPidParticleSim.fData.itsHitData.iIndMatch", HPidParticleSim_fData_itsHitData_iIndMatch, &b_HPidParticleSim_fData_itsHitData_iIndMatch);
   fChain->SetBranchAddress("HPidParticleSim.fData.itsHitData.fDistanceToVertex[10]", HPidParticleSim_fData_itsHitData_fDistanceToVertex, &b_HPidParticleSim_fData_itsHitData_fDistanceToVertex);
   fChain->SetBranchAddress("HPidParticleSim.fData.itsHitData.hasRingCorrelation[10]", HPidParticleSim_fData_itsHitData_hasRingCorrelation, &b_HPidParticleSim_fData_itsHitData_hasRingCorrelation);
   fChain->SetBranchAddress("HPidParticleSim.fData.itsHitData.hasMetaTrackCorrelation[10]", HPidParticleSim_fData_itsHitData_hasMetaTrackCorrelation, &b_HPidParticleSim_fData_itsHitData_hasMetaTrackCorrelation);
   fChain->SetBranchAddress("HPidParticleSim.fData.itsTrackData.fUniqueID", HPidParticleSim_fData_itsTrackData_fUniqueID, &b_HPidParticleSim_fData_itsTrackData_fUniqueID);
   fChain->SetBranchAddress("HPidParticleSim.fData.itsTrackData.fBits", HPidParticleSim_fData_itsTrackData_fBits, &b_HPidParticleSim_fData_itsTrackData_fBits);
   fChain->SetBranchAddress("HPidParticleSim.fData.itsTrackData.cov[10]", HPidParticleSim_fData_itsTrackData_cov, &b_HPidParticleSim_fData_itsTrackData_cov);
   fChain->SetBranchAddress("HPidParticleSim.fData.itsTrackData.nBestMomAlg", HPidParticleSim_fData_itsTrackData_nBestMomAlg, &b_HPidParticleSim_fData_itsTrackData_nBestMomAlg);
   fChain->SetBranchAddress("HPidParticleSim.fData.itsTrackData.nRKTrackInd", HPidParticleSim_fData_itsTrackData_nRKTrackInd, &b_HPidParticleSim_fData_itsTrackData_nRKTrackInd);
   fChain->SetBranchAddress("HPidParticleSim.fData.itsTrackData.nKickTrackInd", HPidParticleSim_fData_itsTrackData_nKickTrackInd, &b_HPidParticleSim_fData_itsTrackData_nKickTrackInd);
   fChain->SetBranchAddress("HPidParticleSim.fData.itsTrackData.nKickTrack123Ind", HPidParticleSim_fData_itsTrackData_nKickTrack123Ind, &b_HPidParticleSim_fData_itsTrackData_nKickTrack123Ind);
   fChain->SetBranchAddress("HPidParticleSim.fData.itsTrackData.nRefTrackInd", HPidParticleSim_fData_itsTrackData_nRefTrackInd, &b_HPidParticleSim_fData_itsTrackData_nRefTrackInd);
   fChain->SetBranchAddress("HPidParticleSim.fData.itsTrackData.nSplineTrackInd", HPidParticleSim_fData_itsTrackData_nSplineTrackInd, &b_HPidParticleSim_fData_itsTrackData_nSplineTrackInd);
   fChain->SetBranchAddress("HPidParticleSim.fData.itsTrackData.fMetaMatchingQuality", HPidParticleSim_fData_itsTrackData_fMetaMatchingQuality, &b_HPidParticleSim_fData_itsTrackData_fMetaMatchingQuality);
   fChain->SetBranchAddress("HPidParticleSim.fData.itsTrackData.fRKRichMatchingQuality", HPidParticleSim_fData_itsTrackData_fRKRichMatchingQuality, &b_HPidParticleSim_fData_itsTrackData_fRKRichMatchingQuality);
   fChain->SetBranchAddress("HPidParticleSim.fData.itsTrackData.dxRkMeta", HPidParticleSim_fData_itsTrackData_dxRkMeta, &b_HPidParticleSim_fData_itsTrackData_dxRkMeta);
   fChain->SetBranchAddress("HPidParticleSim.fData.itsTrackData.dyRkMeta", HPidParticleSim_fData_itsTrackData_dyRkMeta, &b_HPidParticleSim_fData_itsTrackData_dyRkMeta);
   fChain->SetBranchAddress("HPidParticleSim.fData.itsTrackData.dzRkMeta", HPidParticleSim_fData_itsTrackData_dzRkMeta, &b_HPidParticleSim_fData_itsTrackData_dzRkMeta);
   fChain->SetBranchAddress("HPidParticleSim.fData.itsTrackData.dxMdcMeta", HPidParticleSim_fData_itsTrackData_dxMdcMeta, &b_HPidParticleSim_fData_itsTrackData_dxMdcMeta);
   fChain->SetBranchAddress("HPidParticleSim.fData.itsTrackData.dyMdcMeta", HPidParticleSim_fData_itsTrackData_dyMdcMeta, &b_HPidParticleSim_fData_itsTrackData_dyMdcMeta);
   fChain->SetBranchAddress("HPidParticleSim.fData.itsTrackData.dzMdcMeta", HPidParticleSim_fData_itsTrackData_dzMdcMeta, &b_HPidParticleSim_fData_itsTrackData_dzMdcMeta);
   fChain->SetBranchAddress("HPidParticleSim.fData.itsTrackData.xMeta", HPidParticleSim_fData_itsTrackData_xMeta, &b_HPidParticleSim_fData_itsTrackData_xMeta);
   fChain->SetBranchAddress("HPidParticleSim.fData.itsTrackData.yMeta", HPidParticleSim_fData_itsTrackData_yMeta, &b_HPidParticleSim_fData_itsTrackData_yMeta);
   fChain->SetBranchAddress("HPidParticleSim.fData.itsTrackData.zMeta", HPidParticleSim_fData_itsTrackData_zMeta, &b_HPidParticleSim_fData_itsTrackData_zMeta);
   fChain->SetBranchAddress("HPidParticleSim.fData.itsTrackData.errXMeta", HPidParticleSim_fData_itsTrackData_errXMeta, &b_HPidParticleSim_fData_itsTrackData_errXMeta);
   fChain->SetBranchAddress("HPidParticleSim.fData.itsTrackData.errYMeta", HPidParticleSim_fData_itsTrackData_errYMeta, &b_HPidParticleSim_fData_itsTrackData_errYMeta);
   fChain->SetBranchAddress("HPidParticleSim.fData.itsTrackData.errZMeta", HPidParticleSim_fData_itsTrackData_errZMeta, &b_HPidParticleSim_fData_itsTrackData_errZMeta);
   fChain->SetBranchAddress("HPidParticleSim.fData.itsTrackData.nCloseTracklets", HPidParticleSim_fData_itsTrackData_nCloseTracklets, &b_HPidParticleSim_fData_itsTrackData_nCloseTracklets);
   fChain->SetBranchAddress("HPidParticleSim.fData.itsTrackData.fPull", HPidParticleSim_fData_itsTrackData_fPull, &b_HPidParticleSim_fData_itsTrackData_fPull);
   fChain->SetBranchAddress("HPidParticleSim.fData.itsTrackData.fSplineChiSquare", HPidParticleSim_fData_itsTrackData_fSplineChiSquare, &b_HPidParticleSim_fData_itsTrackData_fSplineChiSquare);
   fChain->SetBranchAddress("HPidParticleSim.fData.itsTrackData.fRKChiSquare", HPidParticleSim_fData_itsTrackData_fRKChiSquare, &b_HPidParticleSim_fData_itsTrackData_fRKChiSquare);
   fChain->SetBranchAddress("HPidParticleSim.fData.itsTrackData.iIndexClosestTracklet", HPidParticleSim_fData_itsTrackData_iIndexClosestTracklet, &b_HPidParticleSim_fData_itsTrackData_iIndexClosestTracklet);
   fChain->SetBranchAddress("HPidParticleSim.fData.itsTrackData.aTrackletClusInf0", HPidParticleSim_fData_itsTrackData_aTrackletClusInf0, &b_HPidParticleSim_fData_itsTrackData_aTrackletClusInf0);
   fChain->SetBranchAddress("HPidParticleSim.fData.itsTrackData.aTrackletClusInf1", HPidParticleSim_fData_itsTrackData_aTrackletClusInf1, &b_HPidParticleSim_fData_itsTrackData_aTrackletClusInf1);
   fChain->SetBranchAddress("HPidParticleSim.fData.itsTrackData.aTrackletDistances", HPidParticleSim_fData_itsTrackData_aTrackletDistances, &b_HPidParticleSim_fData_itsTrackData_aTrackletDistances);
   fChain->SetBranchAddress("HPidParticleSim.fData.itsTrackData.qIOMatching[10]", HPidParticleSim_fData_itsTrackData_qIOMatching, &b_HPidParticleSim_fData_itsTrackData_qIOMatching);
   fChain->SetBranchAddress("HPidParticleSim.fData.itsTrackData.nPolarity[10]", HPidParticleSim_fData_itsTrackData_nPolarity, &b_HPidParticleSim_fData_itsTrackData_nPolarity);
   fChain->SetBranchAddress("HPidParticleSim.fData.itsTrackData.fMomenta[10]", HPidParticleSim_fData_itsTrackData_fMomenta, &b_HPidParticleSim_fData_itsTrackData_fMomenta);
   fChain->SetBranchAddress("HPidParticleSim.fData.itsTrackData.fMomError[10]", HPidParticleSim_fData_itsTrackData_fMomError, &b_HPidParticleSim_fData_itsTrackData_fMomError);
   fChain->SetBranchAddress("HPidParticleSim.fData.itsTrackData.fTrackR[10]", HPidParticleSim_fData_itsTrackData_fTrackR, &b_HPidParticleSim_fData_itsTrackData_fTrackR);
   fChain->SetBranchAddress("HPidParticleSim.fData.itsTrackData.fTrackZ[10]", HPidParticleSim_fData_itsTrackData_fTrackZ, &b_HPidParticleSim_fData_itsTrackData_fTrackZ);
   fChain->SetBranchAddress("HPidParticleSim.fData.itsTrackData.fRKPhi", HPidParticleSim_fData_itsTrackData_fRKPhi, &b_HPidParticleSim_fData_itsTrackData_fRKPhi);
   fChain->SetBranchAddress("HPidParticleSim.fData.itsTrackData.fRKTheta", HPidParticleSim_fData_itsTrackData_fRKTheta, &b_HPidParticleSim_fData_itsTrackData_fRKTheta);
   fChain->SetBranchAddress("HPidParticleSim.fData.itsTrackData.fCorrectedEloss[10]", HPidParticleSim_fData_itsTrackData_fCorrectedEloss, &b_HPidParticleSim_fData_itsTrackData_fCorrectedEloss);
   fChain->SetBranchAddress("HPidParticleSim.fData.itsTrackData.fCorrectedBeta[10]", HPidParticleSim_fData_itsTrackData_fCorrectedBeta, &b_HPidParticleSim_fData_itsTrackData_fCorrectedBeta);
   fChain->SetBranchAddress("HPidParticleSim.fData.itsTrackData.fPathLength[10]", HPidParticleSim_fData_itsTrackData_fPathLength, &b_HPidParticleSim_fData_itsTrackData_fPathLength);
   fChain->SetBranchAddress("HPidParticleSim.fData.itsTrackData.fMassSquared[10]", HPidParticleSim_fData_itsTrackData_fMassSquared, &b_HPidParticleSim_fData_itsTrackData_fMassSquared);
   fChain->SetBranchAddress("HPidParticleSim.fData.itsTrackData.bIsAccepted[10]", HPidParticleSim_fData_itsTrackData_bIsAccepted, &b_HPidParticleSim_fData_itsTrackData_bIsAccepted);
   fChain->SetBranchAddress("HPidParticleSim.fData.itsTrackData.nTofRecFlag[10]", HPidParticleSim_fData_itsTrackData_nTofRecFlag, &b_HPidParticleSim_fData_itsTrackData_nTofRecFlag);
   fChain->SetBranchAddress("HPidParticleSim.fData.itsTrackData.fTofRecTof[10]", HPidParticleSim_fData_itsTrackData_fTofRecTof, &b_HPidParticleSim_fData_itsTrackData_fTofRecTof);
   fChain->SetBranchAddress("HPidParticleSim.fData.itsTrackData.fTofRecBeta[10]", HPidParticleSim_fData_itsTrackData_fTofRecBeta, &b_HPidParticleSim_fData_itsTrackData_fTofRecBeta);
   fChain->SetBranchAddress("HPidParticleSim.fData.itsTrackData.fTofRecMassSquared[10]", HPidParticleSim_fData_itsTrackData_fTofRecMassSquared, &b_HPidParticleSim_fData_itsTrackData_fTofRecMassSquared);
   fChain->SetBranchAddress("HPidParticleSim.fData.flags", HPidParticleSim_fData_flags, &b_HPidParticleSim_fData_flags);
   fChain->SetBranchAddress("HPidParticleSim.fData.itsGeantTrackSet.fUniqueID", HPidParticleSim_fData_itsGeantTrackSet_fUniqueID, &b_HPidParticleSim_fData_itsGeantTrackSet_fUniqueID);
   fChain->SetBranchAddress("HPidParticleSim.fData.itsGeantTrackSet.fBits", HPidParticleSim_fData_itsGeantTrackSet_fBits, &b_HPidParticleSim_fData_itsGeantTrackSet_fBits);
   fChain->SetBranchAddress("HPidParticleSim.fData.itsGeantTrackSet.isSorted", HPidParticleSim_fData_itsGeantTrackSet_isSorted, &b_HPidParticleSim_fData_itsGeantTrackSet_isSorted);
   fChain->SetBranchAddress("HPidParticleSim.fData.itsGeantTrackSet.sNCorrTrackIds", HPidParticleSim_fData_itsGeantTrackSet_sNCorrTrackIds, &b_HPidParticleSim_fData_itsGeantTrackSet_sNCorrTrackIds);
   fChain->SetBranchAddress("HPidParticleSim.fData.itsGeantTrackSet.nRICHTracks", HPidParticleSim_fData_itsGeantTrackSet_nRICHTracks, &b_HPidParticleSim_fData_itsGeantTrackSet_nRICHTracks);
   fChain->SetBranchAddress("HPidParticleSim.fData.itsGeantTrackSet.nRICHIPUTracks", HPidParticleSim_fData_itsGeantTrackSet_nRICHIPUTracks, &b_HPidParticleSim_fData_itsGeantTrackSet_nRICHIPUTracks);
   fChain->SetBranchAddress("HPidParticleSim.fData.itsGeantTrackSet.nInnerMdcTracks", HPidParticleSim_fData_itsGeantTrackSet_nInnerMdcTracks, &b_HPidParticleSim_fData_itsGeantTrackSet_nInnerMdcTracks);
   fChain->SetBranchAddress("HPidParticleSim.fData.itsGeantTrackSet.nOuterMdcTracks", HPidParticleSim_fData_itsGeantTrackSet_nOuterMdcTracks, &b_HPidParticleSim_fData_itsGeantTrackSet_nOuterMdcTracks);
   fChain->SetBranchAddress("HPidParticleSim.fData.itsGeantTrackSet.nShowerTracks", HPidParticleSim_fData_itsGeantTrackSet_nShowerTracks, &b_HPidParticleSim_fData_itsGeantTrackSet_nShowerTracks);
   fChain->SetBranchAddress("HPidParticleSim.fData.itsGeantTrackSet.nTOFTracks", HPidParticleSim_fData_itsGeantTrackSet_nTOFTracks, &b_HPidParticleSim_fData_itsGeantTrackSet_nTOFTracks);
   fChain->SetBranchAddress("HPidParticleSim.fData.itsGeantTrackSet.bIsLepFromPrimary", HPidParticleSim_fData_itsGeantTrackSet_bIsLepFromPrimary, &b_HPidParticleSim_fData_itsGeantTrackSet_bIsLepFromPrimary);
   fChain->SetBranchAddress("HPidParticleSim.fData.itsGeantTrackSet.correlatedTrackIds", HPidParticleSim_fData_itsGeantTrackSet_correlatedTrackIds, &b_HPidParticleSim_fData_itsGeantTrackSet_correlatedTrackIds);
   fChain->SetBranchAddress("HPidParticleSim.fData.itsGeantTrackSet.correlationFlags", HPidParticleSim_fData_itsGeantTrackSet_correlationFlags, &b_HPidParticleSim_fData_itsGeantTrackSet_correlationFlags);
   fChain->SetBranchAddress("HPidParticleSim.fData.itsGeantTrackSet.ProcessIds", HPidParticleSim_fData_itsGeantTrackSet_ProcessIds, &b_HPidParticleSim_fData_itsGeantTrackSet_ProcessIds);
   fChain->SetBranchAddress("HPidParticleSim.fData.itsGeantTrackSet.ParentIds", HPidParticleSim_fData_itsGeantTrackSet_ParentIds, &b_HPidParticleSim_fData_itsGeantTrackSet_ParentIds);
   fChain->SetBranchAddress("HPidParticleSim.fData.itsGeantTrackSet.Parents", HPidParticleSim_fData_itsGeantTrackSet_Parents, &b_HPidParticleSim_fData_itsGeantTrackSet_Parents);
   fChain->SetBranchAddress("HPidParticleSim.fData.itsGeantTrackSet.GrandParentIds", HPidParticleSim_fData_itsGeantTrackSet_GrandParentIds, &b_HPidParticleSim_fData_itsGeantTrackSet_GrandParentIds);
   fChain->SetBranchAddress("HPidParticleSim.fData.itsGeantTrackSet.GrandParents", HPidParticleSim_fData_itsGeantTrackSet_GrandParents, &b_HPidParticleSim_fData_itsGeantTrackSet_GrandParents);
   fChain->SetBranchAddress("HPidParticleSim.fData.itsGeantTrackSet.GenInfo", HPidParticleSim_fData_itsGeantTrackSet_GenInfo, &b_HPidParticleSim_fData_itsGeantTrackSet_GenInfo);
   fChain->SetBranchAddress("HPidParticleSim.fData.itsGeantTrackSet.GenInfo1", HPidParticleSim_fData_itsGeantTrackSet_GenInfo1, &b_HPidParticleSim_fData_itsGeantTrackSet_GenInfo1);
   fChain->SetBranchAddress("HPidParticleSim.fData.itsGeantTrackSet.GenInfo2", HPidParticleSim_fData_itsGeantTrackSet_GenInfo2, &b_HPidParticleSim_fData_itsGeantTrackSet_GenInfo2);
   fChain->SetBranchAddress("HPidParticleSim.fData.itsGeantTrackSet.GenWeight", HPidParticleSim_fData_itsGeantTrackSet_GenWeight, &b_HPidParticleSim_fData_itsGeantTrackSet_GenWeight);
   fChain->SetBranchAddress("HPidParticleSim.fData.itsGeantTrackSet.VertexX", HPidParticleSim_fData_itsGeantTrackSet_VertexX, &b_HPidParticleSim_fData_itsGeantTrackSet_VertexX);
   fChain->SetBranchAddress("HPidParticleSim.fData.itsGeantTrackSet.VertexY", HPidParticleSim_fData_itsGeantTrackSet_VertexY, &b_HPidParticleSim_fData_itsGeantTrackSet_VertexY);
   fChain->SetBranchAddress("HPidParticleSim.fData.itsGeantTrackSet.VertexZ", HPidParticleSim_fData_itsGeantTrackSet_VertexZ, &b_HPidParticleSim_fData_itsGeantTrackSet_VertexZ);
   fChain->SetBranchAddress("HPidParticleSim.fData.itsGeantTrackSet.GeantPIDs", HPidParticleSim_fData_itsGeantTrackSet_GeantPIDs, &b_HPidParticleSim_fData_itsGeantTrackSet_GeantPIDs);
   fChain->SetBranchAddress("HPidParticleSim.fData.itsGeantTrackSet.MediumIds", HPidParticleSim_fData_itsGeantTrackSet_MediumIds, &b_HPidParticleSim_fData_itsGeantTrackSet_MediumIds);
   fChain->SetBranchAddress("HPidParticleSim.fData.itsGeantTrackSet.GeantMomX", HPidParticleSim_fData_itsGeantTrackSet_GeantMomX, &b_HPidParticleSim_fData_itsGeantTrackSet_GeantMomX);
   fChain->SetBranchAddress("HPidParticleSim.fData.itsGeantTrackSet.GeantMomY", HPidParticleSim_fData_itsGeantTrackSet_GeantMomY, &b_HPidParticleSim_fData_itsGeantTrackSet_GeantMomY);
   fChain->SetBranchAddress("HPidParticleSim.fData.itsGeantTrackSet.GeantMomZ", HPidParticleSim_fData_itsGeantTrackSet_GeantMomZ, &b_HPidParticleSim_fData_itsGeantTrackSet_GeantMomZ);
   fChain->SetBranchAddress("HPidParticleSim.fData.itsGeantTrackSet.ShowerWeights", HPidParticleSim_fData_itsGeantTrackSet_ShowerWeights, &b_HPidParticleSim_fData_itsGeantTrackSet_ShowerWeights);
   fChain->SetBranchAddress("HPidParticleSim.fData.itsGeantTrackSet.TOFWeights", HPidParticleSim_fData_itsGeantTrackSet_TOFWeights, &b_HPidParticleSim_fData_itsGeantTrackSet_TOFWeights);
   fChain->SetBranchAddress("HPidParticleSim.fData.itsGeantTrackSet.RICHWeights", HPidParticleSim_fData_itsGeantTrackSet_RICHWeights, &b_HPidParticleSim_fData_itsGeantTrackSet_RICHWeights);
   fChain->SetBranchAddress("HPidParticleSim.fData.itsGeantTrackSet.RICHIPUWeights", HPidParticleSim_fData_itsGeantTrackSet_RICHIPUWeights, &b_HPidParticleSim_fData_itsGeantTrackSet_RICHIPUWeights);
   fChain->SetBranchAddress("HPidParticleSim.fData.itsGeantTrackSet.InnerMDCWeights", HPidParticleSim_fData_itsGeantTrackSet_InnerMDCWeights, &b_HPidParticleSim_fData_itsGeantTrackSet_InnerMDCWeights);
   fChain->SetBranchAddress("HPidParticleSim.fData.itsGeantTrackSet.OuterMDCWeights", HPidParticleSim_fData_itsGeantTrackSet_OuterMDCWeights, &b_HPidParticleSim_fData_itsGeantTrackSet_OuterMDCWeights);
   fChain->SetBranchAddress("HPidParticleSim.fData.itsGeantTrackSet.hasHitInShower", HPidParticleSim_fData_itsGeantTrackSet_hasHitInShower, &b_HPidParticleSim_fData_itsGeantTrackSet_hasHitInShower);
   fChain->SetBranchAddress("HPidParticleSim.fData.itsGeantTrackSet.hasHitInTOF", HPidParticleSim_fData_itsGeantTrackSet_hasHitInTOF, &b_HPidParticleSim_fData_itsGeantTrackSet_hasHitInTOF);
   fChain->SetBranchAddress("HPidParticleSim.fNDataObjs", &HPidParticleSim_fNDataObjs, &b_HPidParticleSim_fNDataObjs);
   fChain->SetBranchAddress("HPidParticleSim.hasDynamicObjects", &HPidParticleSim_hasDynamicObjects, &b_HPidParticleSim_hasDynamicObjects);
   fChain->SetBranchAddress("EventHeader.TObject.fUniqueID", &EventHeader_TObject_fUniqueID, &b_EventHeader_TObject_fUniqueID);
   fChain->SetBranchAddress("EventHeader.TObject.fBits", &EventHeader_TObject_fBits, &b_EventHeader_TObject_fBits);
   fChain->SetBranchAddress("EventHeader.fVertex.fUniqueID", &EventHeader_fVertex_fUniqueID, &b_EventHeader_fVertex_fUniqueID);
   fChain->SetBranchAddress("EventHeader.fVertex.fBits", &EventHeader_fVertex_fBits, &b_EventHeader_fVertex_fBits);
   fChain->SetBranchAddress("EventHeader.fVertex.pos.fUniqueID", &EventHeader_fVertex_pos_fUniqueID, &b_EventHeader_fVertex_pos_fUniqueID);
   fChain->SetBranchAddress("EventHeader.fVertex.pos.fBits", &EventHeader_fVertex_pos_fBits, &b_EventHeader_fVertex_pos_fBits);
   fChain->SetBranchAddress("EventHeader.fVertex.pos.x", &EventHeader_fVertex_pos_x, &b_EventHeader_fVertex_pos_x);
   fChain->SetBranchAddress("EventHeader.fVertex.pos.y", &EventHeader_fVertex_pos_y, &b_EventHeader_fVertex_pos_y);
   fChain->SetBranchAddress("EventHeader.fVertex.pos.z", &EventHeader_fVertex_pos_z, &b_EventHeader_fVertex_pos_z);
   fChain->SetBranchAddress("EventHeader.fVertex.iterations", &EventHeader_fVertex_iterations, &b_EventHeader_fVertex_iterations);
   fChain->SetBranchAddress("EventHeader.fVertex.chi2", &EventHeader_fVertex_chi2, &b_EventHeader_fVertex_chi2);
   fChain->SetBranchAddress("EventHeader.fVertex.sumOfWeights", &EventHeader_fVertex_sumOfWeights, &b_EventHeader_fVertex_sumOfWeights);
   fChain->SetBranchAddress("EventHeader.timeInSpill", &EventHeader_timeInSpill, &b_EventHeader_timeInSpill);
   fChain->SetBranchAddress("EventHeader.downscaling", &EventHeader_downscaling, &b_EventHeader_downscaling);
   fChain->SetBranchAddress("EventHeader.downscalingFlag", &EventHeader_downscalingFlag, &b_EventHeader_downscalingFlag);
   fChain->SetBranchAddress("EventHeader.fDate", &EventHeader_fDate, &b_EventHeader_fDate);
   fChain->SetBranchAddress("EventHeader.fErrorBit", &EventHeader_fErrorBit, &b_EventHeader_fErrorBit);
   fChain->SetBranchAddress("EventHeader.fEventDecoding", &EventHeader_fEventDecoding, &b_EventHeader_fEventDecoding);
   fChain->SetBranchAddress("EventHeader.fEventPad", &EventHeader_fEventPad, &b_EventHeader_fEventPad);
   fChain->SetBranchAddress("EventHeader.fEventRunNumber", &EventHeader_fEventRunNumber, &b_EventHeader_fEventRunNumber);
   fChain->SetBranchAddress("EventHeader.fEventSeqNumber", &EventHeader_fEventSeqNumber, &b_EventHeader_fEventSeqNumber);
   fChain->SetBranchAddress("EventHeader.fEventSize", &EventHeader_fEventSize, &b_EventHeader_fEventSize);
   fChain->SetBranchAddress("EventHeader.fId", &EventHeader_fId, &b_EventHeader_fId);
   fChain->SetBranchAddress("EventHeader.fTBit", &EventHeader_fTBit, &b_EventHeader_fTBit);
   fChain->SetBranchAddress("EventHeader.fTime", &EventHeader_fTime, &b_EventHeader_fTime);
   fChain->SetBranchAddress("EventHeader.fVersion", &EventHeader_fVersion, &b_EventHeader_fVersion);
   fChain->SetBranchAddress("EventHeader.triggerDecision", &EventHeader_triggerDecision, &b_EventHeader_triggerDecision);
   fChain->SetBranchAddress("EventHeader.triggerDecisionEmu", &EventHeader_triggerDecisionEmu, &b_EventHeader_triggerDecisionEmu);
   Notify();
}

Bool_t T::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void T::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t T::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
  bool from_lambda=(HGeantKine_fData_parentTrack[0]==18);
  bool is_pim=(HGeantKine_fData_particleID[0]==9);
  bool is_p=(HGeantKine_fData_particleID[0]==9);

  for(int i=0; i<kMaxHGeantKine_fData;i++)
    for(int j=i; j<kMaxHGeantKine_fData;j++)
  cout<<"from lambda "<<from_lambda<<" ID: "<< HGeantKine_fData_particleID[0] << endl;
  if(from_lambda && (is_p || is_pim))
    return 1;
  else
    return -1;
}
#endif // #ifdef T_cxx
