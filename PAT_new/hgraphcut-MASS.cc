//-----------------------------------------------------------------------------
#include "hgraphcut.h"
#include "hhypplayer.h"
#include "hcommondef.h"
#include <TError.h>
#include <TMath.h>
#include <string>


// ****************************************************************************
ClassImp(HGraphCut)

using namespace std;
using namespace CommonDefinitions;


//-----------------------------------------------------------------------------
Bool_t HGraphCut::select(HReconstructor& rec)
{
   HHypPlayer& dataRec = dynamic_cast<HHypPlayer&>(rec);
   HPidPool *pPool = dataRec.getPool();
   //cout << "PID POOL : " << pPool->size() << endl;
   HHypCandidate *pCand = 0;
   if ( getName() == string("all") )
   {
      int i = 0;
      while (  ( pCand = pPool->getHyp( i ) ) != 0 )
      {
         graphCut( pCand );
         ++i;
      }
   }
   else
   {
      int i = 0;
      while ( ( pCand = pPool->getHyp( getName(), i ) ) != 0 )
      {
         graphCut( pCand );
         ++i;
      }
   }

return kTRUE;
}


//-----------------------------------------------------------------------------
void HGraphCut::graphCut(HHypCandidate *pHyp)
{
   for (int i=0; i<pHyp->getSize(); ++i)
   {
      if ( graphCut( pHyp->getPart(i) ) == false )
	  {
	     pHyp->setActive( false );
	  }
   }
}


//-----------------------------------------------------------------------------
bool HGraphCut::graphCut(HParticle *pPart)
{
   id = pPart->getId();
   mom = pPart->get("p");
   beta = pPart->get("beta_new");
   double mass = mom / ( beta * TMath::Sqrt( 1. / ( 1. - (beta*beta) ) ) );

   switch ( id )
   {
       case 2: return true; // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! always ok
               if (ep_cut) return ep_cut->IsInside( beta, mom );
	           break;
       case 3: return true; // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! always ok
               if (em_cut) return em_cut->IsInside( beta, mom );
	           break;
       case 8: //return true;
               if (pip_cut) return pip_cut->IsInside( mom, mass );
	           break;
       case 9: return true;
               if (pim_cut) return pim_cut->IsInside( mom, mass ); // same for pi-
       //case 9: if (pip_cut) return pip_cut->IsInside( beta, mom );
	           break;
       case 14: //return true; 
                if (p_cut) return p_cut->IsInside( mom, mass );
	            break;
       case 45: return true;
                if (d_cut) return d_cut->IsInside( beta, mom );
	            break;
       default: return true;
   }

return true;
}



