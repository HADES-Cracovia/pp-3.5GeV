#include "hades.h"
#include "hevent.h"
#include "heventheader.h"
#include "hparticlepool.h"
#include "hwallhit.h"
#include <algorithm>
#include <iostream>

using namespace std;


// ****************************************************************************
ClassImp(HParticlePool)



//-----------------------------------------------------------------------------
HParticlePool::HParticlePool(HOutputFile *ptr) : HPool(ptr), HParticleDataPool() 
{
   reset();
}


//-----------------------------------------------------------------------------
HParticlePool::~HParticlePool() 
{ 
   reset(); 
}


//-----------------------------------------------------------------------------
void HParticlePool::reset()
{
   partNum[ePositron] = 0;
   partNum[eElectron] = 0;
   partNum[ePiPlus] = 0;
   partNum[ePiMinus] = 0;
   partNum[eProton] = 0;
   partNum[eLeptonPos] = 0;
   partNum[eLeptonNeg] = 0;
   partNum[eHadronPos] = 0;
   partNum[eHadronNeg] = 0;

   std::for_each( partCand.begin(), partCand.end(), DelParticleCand() ); 
   partCand.clear();
}

//-----------------------------------------------------------------------------
void HParticlePool::loop(HIterator* dataIt)
{
   HPidTrackCand *PidCand    = 0;

   Bool_t isPositive = kFALSE;
   Bool_t isRing = kFALSE;

   dataIt->Reset();
   while (((PidCand = (HPidTrackCand *) dataIt->Next()) != 0))
   {
      isPositive = ( PidCand->getTrackData()->getPolarity(4) > 0 ) ? kTRUE : kFALSE;
      isRing = PidCand->getHitData()->getRingCorrelation(4);

      EParticle myEId = eUnknown;
      if (isPositive && isRing) myEId = eLeptonPos;
      else if (isPositive && !isRing) myEId = eHadronPos;
      else if (!isPositive && isRing) myEId = eLeptonNeg;
      else if (!isPositive && !isRing) myEId = eHadronNeg;

      if (getNum(myEId)<30)
      {
         // track must be fully reconstructed otherwise rejected
         if (PidCand->isFlagBit(HPidTrackCand::kIsUsed) == 1)
           addPartCand(myEId, PidCand);
      }
      else
      {
         cerr << "Particle " << convertId( myEId ) << ": more than 30 tracks - event rejected!" << endl;
         break;
      }
   }

}


//-----------------------------------------------------------------------------
void HParticlePool::fill(HIterator* dataIt) 
{
    HParticleCandidate *pc = 0;
    EParticle ePcTab[] = { eHadronPos, eHadronNeg, eLeptonPos, eLeptonNeg };
	unsigned pcTabSize = sizeof(ePcTab)/sizeof(EParticle);
    int count = 0;
    int totalcount = 0;
    int totalmult = 0;

	eventData.update();

    // header information
    HEventHeader *evHeader = gHades->getCurrentEvent()->getHeader();
		 
    eventData.set("trigdownscale", evHeader->getDownscaling() );
    eventData.set("trigdownscaleflag", evHeader->getDownscalingFlag() );
    eventData.set("trigdec", evHeader->getTriggerDecision() );
    eventData.set("trigbit", evHeader->getTBit() );
    eventData.set("event", evHeader->getEventSeqNumber() );
    eventData.set("runnumber", evHeader->getEventRunNumber() );

    for (unsigned i=0; i<pcTabSize; ++i)
	{
	   eventData.set( convertId( ePcTab[i] ) + "_mult" , getNum( ePcTab[i] ) );
	   totalmult += getNum( ePcTab[i] );
	}
	eventData.set("totalmult", totalmult);

	HVertex& vertex = gHades->getCurrentEvent()->getHeader()->getVertex();

	eventData.set("eVert_x", vertex.getX());
	eventData.set("eVert_y", vertex.getY());
	eventData.set("eVert_z", vertex.getZ());
	eventData.set("eVert_chi2", vertex.getChi2());

	// ----------- forward wall data ----------------------------------
	if ( dataIt != 0)
	{
	   // ------------- data reset --------------------------
             eventData.set("fw_mult", -100. );
             eventData.set("fw_time", -100. );
             eventData.set("fw_charge", -100. );
             eventData.set("fw_cell", -100. );
             eventData.set("fw_theta", -100. );
             eventData.set("fw_phi", -100. );
             eventData.set("fw_distance", -100. );
             eventData.set("fw_x_lab", -100. );
             eventData.set("fw_y_lab", -100. );
             eventData.set("fw_z_lab", -100. );
             eventData.set("fw_beta", -100. );
             eventData.set("fw_p", -100. );
       // ---------------------------------------------------
       HWallHit *WallHitCand    = 0;
	   int counter = 0;
	   int i_hit = 0;
	   float theta_smallest = 100.;
	   double fw_beta, fw_gamma, fw_mom;
	   double fw_beta_OK = -100., fw_mom_OK = -100., fw_mult = -100.;
	   
	   dataIt->Reset();
	   while (((WallHitCand = (HWallHit *) dataIt->Next()) != 0))
	   {
           ++counter;
		   if (WallHitCand->getTheta() <= theta_smallest)
		   {
		      // ---------------- momentum calculation -------------------
	          fw_beta = ((WallHitCand->getDistance()*cos(WallHitCand->getTheta()*TMath::DegToRad()))/(3*WallHitCand->getTime()*1e2));
			  fw_gamma = sqrt(1 - (fw_beta*fw_beta));
			  fw_mom = (fw_beta*0.93827231)/fw_gamma;
              if (fw_mom > 1.6 && fw_mom < 2.6) 
			  {
		         theta_smallest = WallHitCand->getTheta();
			     i_hit = counter;
				 fw_beta_OK = fw_beta;
				 fw_mom_OK = fw_mom;
			  }
		   }
       }
	   fw_mult = counter;
	   counter = 0;
	   dataIt->Reset();
	   while (((WallHitCand = (HWallHit *) dataIt->Next()) != 0))
	   {
          ++counter;
		  if (i_hit == counter)
		  {
		     // ----------------- store FWall hit to ntuple ------------------------- 
			 eventData.set("fw_mult", fw_mult );
			 eventData.set("fw_time", WallHitCand->getTime() );
			 eventData.set("fw_charge", WallHitCand->getCharge() );
			 eventData.set("fw_cell", WallHitCand->getCell() );
			 eventData.set("fw_theta", WallHitCand->getTheta() );
			 eventData.set("fw_phi", WallHitCand->getPhi() );
			 eventData.set("fw_distance", WallHitCand->getDistance() );
			 float x_lab, y_lab, z_lab;
			 WallHitCand->getXYZLab( x_lab, y_lab, z_lab );
			 eventData.set("fw_x_lab", x_lab );
			 eventData.set("fw_y_lab", y_lab );
			 eventData.set("fw_z_lab", z_lab );
			 eventData.set("fw_beta", fw_beta_OK );
			 eventData.set("fw_p", fw_mom_OK );
		  }
	   }
	}
	// ----------------------------------------------------------------
             
    if (ptrFile != 0) 
	{
	   totalcount = 0;
       for (unsigned i=0; i<pcTabSize; ++i)
       for (count = 0; count<getNum( ePcTab[i] ); ++count)
       {  
          // internal loop over all particle (set) patterns
          for (objIt = objectList.begin(); objIt != objectList.end(); ++objIt)
          {
       	     if ( (**objIt)[ ePcTab[i] ] > 0 )
       	     {
                (*objIt)->set("cid", count+1 );
                (*objIt)->set("mult", getNum(ePcTab[i]) );
   
                pc = getParticle( ePcTab[i], count );
                for (NtuplePairIter it = pc->begin(); it != pc->end(); ++it)
                {
                   (*objIt)->set( it->first, it->second );
                }
                
                for (NtuplePairIter it = eventData.begin(); it != eventData.end(); ++it)
                {
                   (*objIt)->set( it->first, it->second );
                }
                (*objIt)->set( "totalcid", ++totalcount );
                
                (*objIt)->fill();
       	     }
          }
          // end of internal loop
       }
	} // if ptrFile
}

