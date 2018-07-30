#ifndef HTRACKPLAYER_H
#define HTRACKPLAYER_H

#include <vector>

#include "hreconstructor.h"
#include "hiterator.h"
#include "hevent.h"
#include "hcategory.h"
#include "hades.h"

#include "hntuple.h"
#include "hparticlepool.h"
#include "hcut.h"


class HTrackPlayer : public HReconstructor
{

public:

  HTrackPlayer(HParticlePool& ref);
  ~HTrackPlayer();

  Bool_t          init();
  Bool_t          reinit() { return kTRUE; }
  Int_t           execute();
  Bool_t          finalize() { return kTRUE; }

  HParticlePool* getPool() { return &partPool; }
  void add(HCut &cut) { vCut.push_back( &cut ); }


 private:


  HCategory      *m_pContCatPart;       //!
  HIterator      *m_pContItPart;        //!
  HCategory      *m_pFWallCat;       //!
  HIterator      *m_pFWallCatIter;        //!
  HCategory      *m_pGeantKineCat;    //!
  HIterator      *m_pGeantKineCatIter;   //!


   HParticlePool &partPool;  

   std::vector<HCut*> vCut;
   std::vector<HCut*>::iterator itCut;

  ClassDef(HTrackPlayer, 0)   // Fills HParticlePool with all particle candidates
};

#endif
