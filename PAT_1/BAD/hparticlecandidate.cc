#include "hparticlecandidate.h"
#include "hpidtrackcandsim.h"


ClassImp(HParticleCandidate)

    extern int simflag;

HParticleCandidate::HParticleCandidate(HPidTrackCand* ptr) 
{
   if (ptr)
   {
      HitData   = ptr->getHitData();
      TrackData = ptr->getTrackData();

      set("id", 0. );
      set("sector", HitData->getSector() );
      set("system", HitData->getSystem() );
      set("p", TrackData->getMomenta(4) );
      set("theta", TrackData->getRKTheta() );
      set("phi", TrackData->getRKPhi() );
      set("q", TrackData->getPolarity(4) );
      set("theta_rich", HitData->getRichTheta() );
      set("phi_rich", HitData->getRichPhi() );
      set("r", TrackData->getTrackR(4) );
      set("z", TrackData->getTrackZ(4) );
      set("rkchi2", TrackData->getRKChiSquare() );
      set("mdcchi2", HitData->getInnerMdcChiSquare() );
      set("oa_lept", TrackData->getAngleWithClosestLeptonCandidate()*(int(TrackData->closestLeptonCandidateIsFitted())-0.5)*2 );
      set("oa_hadr", TrackData->getAngleWithClosestHadronCandidate()*(int(TrackData->closestHadronCandidateIsFitted())-0.5)*2 );
      set("isring", HitData->getRingCorrelation(4) );
      set("beta", TrackData->getCorrectedBeta(4) );
      set("dedx_in", HitData->getInnerMdcdEdx() );
      set("dedx_in_sigma", HitData->getInnerMdcdEdxSigma() );
      set("dedx_out", HitData->getOuterMdcdEdx() );
      set("dedx_out_sigma", HitData->getOuterMdcdEdxSigma() );
      set("dedx_mdc", HitData->getCombinedMdcdEdxSigma() );
      set("dedx_mdc_sigma", HitData->getCombinedMdcdEdxSigma() );
      set("dedx_tof", TrackData->getCorrectedEloss(4) );
      set("shw_sum0", HitData->getShowerSum(0) );
      set("shw_sum1", HitData->getShowerSum(1) );
      set("shw_sum2", HitData->getShowerSum(2) );
      set("rich_amp", HitData->getRingAmplitude() );
      set("rich_padnum", HitData->getRingPadNr() );
      set("rich_centr", HitData->getRingCentroid() );
      set("rich_patmat", HitData->getRingPatMat() );
      set("rich_houtra", HitData->getRingHouTra() );
      set("tofino_mult", HitData->getTofinoMult() );

	  if (HitData->getSystem() == 0)
	  {
	     set("resolution", 450); 
             set("tof_exp", HitData->getShowerTimeOfFlight() );
	  }
	  else
	  {
	     set("resolution", 150); 
             set("tof_exp", HitData->getTOFTimeOfFlight() );
	  }

/*
      set("kIsDoubleHitRICH", ptr->isFlagBit(HPidTrackCand::kIsDoubleHitRICH) );
      set("kIsDoubleHitInnerMDC", ptr->isFlagBit(HPidTrackCand::kIsDoubleHitInnerMDC) );
      set("kIsDoubleHitOuterMDC", ptr->isFlagBit(HPidTrackCand::kIsDoubleHitOuterMDC) );
      set("kIsDoubleHitMETA", ptr->isFlagBit(HPidTrackCand::kIsDoubleHitMETA) );
      set("kIsBestHitRICH", ptr->isFlagBit(HPidTrackCand::kIsBestHitRICH) );
      set("kIsBestHitInnerMDC", ptr->isFlagBit(HPidTrackCand::kIsBestHitInnerMDC) );
      set("kIsBestHitOuterMDC", ptr->isFlagBit(HPidTrackCand::kIsBestHitOuterMDC) );
      set("kIsBestHitMETA", ptr->isFlagBit(HPidTrackCand::kIsBestHitMETA) );
      set("kIsBestRKMETA", ptr->isFlagBit(HPidTrackCand::kIsBestRKMETA) );
      set("kIsBestRKRICH", ptr->isFlagBit(HPidTrackCand::kIsBestRKRICH) );
      set("kIsBestRK", ptr->isFlagBit(HPidTrackCand::kIsBestRK) );
      set("kIsBestSPLINE", ptr->isFlagBit(HPidTrackCand::kIsBestSPLINE) );
      set("kIsBestKICK", ptr->isFlagBit(HPidTrackCand::kIsBestKICK) );
      set("kIsBestRKRKMETA", ptr->isFlagBit(HPidTrackCand::kIsBestRKRKMETA) );
      set("kIsBestRKRKRICH", ptr->isFlagBit(HPidTrackCand::kIsBestRKRKRICH) );
      set("kIsBestRKMETAQA", ptr->isFlagBit(HPidTrackCand::kIsBestRKMETAQA) );
      set("kIsAcceptedHitRICH", ptr->isFlagBit(HPidTrackCand::kIsAcceptedHitRICH) );
      set("kIsAcceptedHitInnerMDC", ptr->isFlagBit(HPidTrackCand::kIsAcceptedHitInnerMDC) );
      set("kIsAcceptedHitOuterMDC", ptr->isFlagBit(HPidTrackCand::kIsAcceptedHitOuterMDC) );
      set("kIsAcceptedHitMETA", ptr->isFlagBit(HPidTrackCand::kIsAcceptedHitMETA) );
      set("kIsAcceptedRKMETA", ptr->isFlagBit(HPidTrackCand::kIsAcceptedRKMETA) );
      set("kIsAcceptedRKRICH", ptr->isFlagBit(HPidTrackCand::kIsAcceptedRKRICH) );
      set("kIsAcceptedRK", ptr->isFlagBit(HPidTrackCand::kIsAcceptedRK) );
      set("kIsAcceptedSPLINE", ptr->isFlagBit(HPidTrackCand::kIsAcceptedSPLINE) );
      set("kIsAcceptedKICK", ptr->isFlagBit(HPidTrackCand::kIsAcceptedKICK) );
      set("kIsAcceptedRKRKMETA", ptr->isFlagBit(HPidTrackCand::kIsAcceptedRKRKMETA) );
      set("kIsAcceptedRKRKRICH", ptr->isFlagBit(HPidTrackCand::kIsAcceptedRKRKRICH) );
      set("kIsAcceptedRKMETAQA", ptr->isFlagBit(HPidTrackCand::kIsAcceptedRKMETAQA) );
*/
      set("kIsLepton", ptr->isFlagBit(HPidTrackCand::kIsLepton) );
      set("kIsUsed", ptr->isFlagBit(HPidTrackCand::kIsUsed) );
//      set("kIsRejected", ptr->isFlagBit(HPidTrackCand::kIsRejected) );
//     set("isCorrelated", 0 );

      // *********************** below part of data from simulation ******************
      HPidTrackCandSim *PidCandSim = dynamic_cast<HPidTrackCandSim*>( ptr );
      if ( PidCandSim != 0 )
      {
		 HPidGeantTrackSet* ptrSim = PidCandSim->getGeantTrackSet();
		 set("sim_iscommon", ptrSim->getMostCommonCorrelation() );
		 set("sim_id", ptrSim->getGeantPID() );
		 set("sim_corrflag", ptrSim->getCorrelationFlag() );
		 set("sim_parentid", ptrSim->getGeantParentID() );
		 set("sim_geninfo", ptrSim->getGenInfo() );
		 set("sim_geninfo1", ptrSim->getGenInfo1() );
		 set("sim_geninfo2", ptrSim->getGenInfo2() );
		 set("sim_genweight", ptrSim->getGenWeight() );
		 set("sim_vertexx", ptrSim->getGeantVertexX() );
		 set("sim_vertexy", ptrSim->getGeantVertexY() );
		 set("sim_vertexz", ptrSim->getGeantVertexZ() );
		 set("sim_mediumid", ptrSim->getGeantMediumID() );
		 set("sim_processid", ptrSim->getGeantProcessID() );
		 set("sim_px", ptrSim->getGeantMomX() );
		 set("sim_py", ptrSim->getGeantMomY() );
		 set("sim_pz", ptrSim->getGeantMomZ() );
		 set("sim_p", ptrSim->getTotalMomentum() );
		 set("sim_primaryflag", ptrSim->getPrimaryFlag() );
      }

   }

}
