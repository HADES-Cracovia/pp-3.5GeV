//-----------------------------------------------------------------------------
#include "htrackcut.h"
#include "hhypdata.h"
#include "htrackplayer.h"
#include "hcommondef.h"
#include <TError.h>
#include <TMath.h>


// ****************************************************************************
ClassImp(HTrackCut)

using namespace CommonDefinitions;
using namespace AnalysisParameters;


//-----------------------------------------------------------------------------
Bool_t HTrackCut::select(HReconstructor& rec)
{
   HTrackPlayer& dataRec = dynamic_cast<HTrackPlayer&>(rec);
   HParticlePool *pPool = dataRec.getPool();
   HParticleCandidate *pCand = 0;
   std::map<HParticleCandidate*, EParticle> replaceCand;
   std::map<HParticleCandidate*, EParticle>::iterator replaceIt;
   replaceCand.clear();

   EParticle oldEId;
   EParticle newEId;
   Bool_t isPositive = kFALSE;
   Bool_t isCorrelated = kFALSE;
   int leptonId = 0;

   for(int i=0; i<pPool->getSize(); ++i)
   {
      pCand = pPool->getParticle(i);
      if ( pCand->isLepton() )
      {
         double delta_theta = pCand->get("theta") - pCand->get("theta_rich");
         double delta_phi = pCand->get("phi") - pCand->get("phi_rich");
         double sin_theta = TMath::Sin( TMath::DegToRad() * pCand->get("theta") );
         isPositive = ( pCand->get("q") > 0 ) ? kTRUE : kFALSE;
	 leptonId = (isPositive==kTRUE) ? 2 : 3;
         oldEId = pCand->getEId();
         isCorrelated = isInside(delta_theta, sin_theta*delta_phi, leptonId, 
                                 (short)pCand->get("system"), (short)pCand->get("sector"), pCand->get("p") );
         // if not correlated we have to remove lepton and add hadron
         if (isCorrelated == false)
         {
            if (isPositive) newEId = eHadronPos;
            else if (!isPositive) newEId = eHadronNeg;
            replaceCand.insert( pair<HParticleCandidate*, EParticle>(pCand, newEId) );
         }
      }
   }

   replaceIt = replaceCand.begin();
   while ( replaceIt != replaceCand.end() )
   {
      pPool->removeParticle(replaceIt->first);
      ++replaceIt;
   }
   replaceIt = replaceCand.begin();
   while ( replaceIt != replaceCand.end() )
   {
      pPool->addPartCand(replaceIt->second, replaceIt->first);
      ++replaceIt;
   }


return kTRUE;
}


//-----------------------------------------------------------------------------
Double_t HTrackCut::getValue(double p1, double p2, double p3, double p4, double p5, double x)
{
   return (p1/x)+p2*x*x*x+p3*x*x+p4*x+p5;
}

//-----------------------------------------------------------------------------
Double_t HTrackCut::getXtheta(Short_t pid, Short_t system, Short_t sector, Float_t mom)
{
    if (system == -1) return -1;
    Float_t p = (mom<50.) ? 50. : mom;
    if (pid==2)
    {
        p = (mom>300.) ? 300. : mom;

            double retVal =  getValue(positrontheta[sector][system][0], positrontheta[sector][system][1],
                            positrontheta[sector][system][2], positrontheta[sector][system][3],
                            positrontheta[sector][system][4], p);
        if(retVal< 1)   cout << retVal << ", mom: " << mom<<", sec: " << sector<< ", syst: " << system<<endl;
        return retVal;
    } else
    if (pid==3)
    {
        p = (mom>400.) ? 400. : mom;
            double retVal =  getValue(electrontheta[sector][system][0], electrontheta[sector][system][1],
                            electrontheta[sector][system][2], electrontheta[sector][system][3],
                            electrontheta[sector][system][4], p);
        if(retVal<1)    cout << retVal << ", mom: " << mom<<", sec: " << sector<< ", syst: " << system<<endl;
        return retVal;

    }

return 0.; // dummy
}


//-----------------------------------------------------------------------------
Double_t HTrackCut::getYphi(Short_t pid, Short_t system, Short_t sector, Float_t mom)
{
    if (system == -1) return -1;
    Float_t p = (mom<50.) ? 50. : mom;
    if (pid==2)
    {
        p = (mom>300.) ? 300. : mom;
            double retVal =  getValue(positronphi[sector][system][0], positronphi[sector][system][1],
                            positronphi[sector][system][2], positronphi[sector][system][3],
                            positronphi[sector][system][4], p);
        if(retVal<1)    cout << retVal << ", mom: " << mom<<", sec: " << sector<< ", syst: " << system<<endl;
        return retVal;
    } else
    if (pid==3)
    {
        p = (mom>400.) ? 400. : mom;
            double retVal =  getValue(electronphi[sector][system][0], electronphi[sector][system][1],
                            electronphi[sector][system][2], electronphi[sector][system][3],
                            electronphi[sector][system][4], p);
        if(retVal<1)    cout << retVal << ", mom: " << mom<<", sec: " << sector<< ", syst: " << system<<endl;

        return retVal;
    }

return 0.; // dummy
}



//-----------------------------------------------------------------------------
Bool_t HTrackCut::isInside(Double_t x1, Double_t y1, Short_t pid, Short_t system, Short_t sector, Float_t mom)
{
 // first parameter: theta, second: sin(theta)*phi
   Double_t a = getXtheta(pid, system, sector, mom);
   Double_t b = getYphi(pid, system, sector, mom);
   if (a==-1 || b==-1) return kFALSE;
   Double_t ell_value = b*TMath::Sqrt(1-x1*x1/(a*a));
   if (a<1 ||b<1) cout << "ell_val: " << ell_value << ", a: " << a << ", b: " << b<<  ", mom: " << mom<<", sec: " << sector<< ", syst: " << system<<endl;
   if (y1 <= ell_value && y1 >= -1.*ell_value)
       return kTRUE;
return kFALSE;
}

