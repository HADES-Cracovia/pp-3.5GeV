#include "htrackplayer.h"
#include "walldef.h"
#include "hgeantkine.h"
#include "heditor.h"
#include "hparticledef.h"

ClassImp(HTrackPlayer)

extern int wallflag;
extern int btflag;
//-----------------------------------------------------------------
  HTrackPlayer::HTrackPlayer(HParticlePool& ref) : 
      HReconstructor( const_cast<Text_t*>("trackplayer"), const_cast<Text_t*>("trackplayer") ), 
	  partPool(ref)
  {
  }

//-----------------------------------------------------------------
  HTrackPlayer::~HTrackPlayer() 
  {
  }



//-----------------------------------------------------------------
Bool_t HTrackPlayer::init()
{
	m_pContItPart = 0;         // Iterator
	m_pContCatPart = 0;        // Category
	m_pFWallCat = 0;
	m_pFWallCatIter = 0;
	m_pGeantKineCat = 0;
	m_pGeantKineCatIter = 0;
	m_pBTCat=0;
	m_pBTCatIter=0;

	if ((m_pContCatPart =
		gHades->getCurrentEvent()->getCategory(catPidTrackCand)) == 0) {
		ErrorMsg(ERROR, "HTrackPlayer::init", 1, "Cannot get catPidTrackCand cat");
		return kFALSE;
	}
	m_pContItPart = (HIterator *) m_pContCatPart->MakeIterator();

        if ( wallflag == 1 )
        {
	   if ((m_pFWallCat =
		   gHades->getCurrentEvent()->getCategory(catWallHit)) == 0) {
		   ErrorMsg(ERROR, "HTrackPlayer::init", 1, "Cannot get catWallHit cat");
		   return kFALSE;
	   }
	   m_pFWallCatIter = (HIterator *) m_pFWallCat->MakeIterator();
        }

        if ( simflag > 0 )
        {
	   if ((m_pGeantKineCat =
		   gHades->getCurrentEvent()->getCategory(catGeantKine)) == 0) {
		   ErrorMsg(ERROR, "HTrackPlayer::init", 1, "Cannot get catGeantKine cat");
		   return kFALSE;
	   }
	   m_pGeantKineCatIter = (HIterator *) m_pGeantKineCat->MakeIterator();
        }

	if ( btflag == 1 )
        {
	   if ((m_pBTCat =
		   gHades->getCurrentEvent()->getCategory(catParticleBtRing)) == 0) {
		   ErrorMsg(ERROR, "HTrackPlayer::init", 1, "Cannot get catParticleBtRing cat");
		   return kFALSE;
	   }
	   m_pBTCatIter = (HIterator *) m_pBTCat->MakeIterator();
        }

	return kTRUE;
}

//-----------------------------------------------------------------
Int_t HTrackPlayer::execute()
{
   partPool.reset();
   partPool.loop( m_pContItPart );
#ifdef DEBUG
   partPool.dump();
#endif
   itCut = vCut.begin();
   while( itCut != vCut.end() )
   {
      (*itCut)->select(*this);
      ++itCut;
   }   
   
   partPool.fill( m_pFWallCatIter, m_pGeantKineCatIter, m_pBTCatIter );

   return 0;
}



