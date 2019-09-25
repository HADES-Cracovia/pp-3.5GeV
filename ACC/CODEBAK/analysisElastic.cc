#include "analysisElastic.h"
#include "hntuple.h"
#include <TH1F.h>
#include <TComplex.h>
#include <TLorentzVector.h>



using namespace std;


   const double D2R = 1.74532925199432955e-02;
   const double R2D = 57.2957795130823229;

   double openingangle(TLorentzVector a, TLorentzVector b)
   {
         return TMath::ACos( (a.Px()*b.Px() + a.Py()*b.Py() +  a.Pz()*b.Pz() ) / ( a.Vect().Mag() * b.Vect().Mag() ) );
   }

   void normalize(TH1* hist)
   {
      for (Int_t j=1; j<hist->GetNbinsX()+1; ++j)
      {
         hist->SetBinContent( j, hist->GetBinContent(j) / hist->GetBinWidth(j) );
//         hist->SetBinError( j, TMath::Sqrt( hist->GetBinContent(j) ) );
         hist->SetBinError( j, hist->GetBinError(j) / hist->GetBinWidth(j) );
      }
   }


Float_t transformPhi(Float_t Phi)
{
   Float_t dPhi;

   if( (dPhi = TMath::RadToDeg() * Phi) < 0.0 )
      dPhi += 360.0;

   return dPhi;
}

Bool_t isGeantTrackInAcceptance(HGeantKine *pG, HIterator *pitMdcGeant, HIterator *pitTofGeant, HIterator *pitShowerGeant)
{
     Int_t lTrack, lId;
     pG->getParticle(lTrack,lId);

     Int_t nStatMDC1 = 0, nStatMDC2 = 0, nStatMDC3 = 0,  nStatMDC4 = 0;
     Int_t nStatTof = 0, nStatShower = 0;
     HGeantMdc *pMdc;
     pitMdcGeant->Reset();
     while((pMdc = (HGeantMdc*) pitMdcGeant->Next()) != NULL){

             if (pMdc->getTrack() == lTrack)
             {
                switch (((Int_t)pMdc->getModule()))
                {
                    case 0: nStatMDC1 = 1;
                            break;
                    case 1: nStatMDC2 = 1;
                            break;
                    case 2: nStatMDC3 = 1;
                            break;
                    case 3: nStatMDC4 = 1;
                            break;
                    default: cerr << "WRONG MDC module number!" << endl;
                }
             }
     }
     HGeantTof *pTof;
     pitTofGeant->Reset();
     while((pTof = (HGeantTof*) pitTofGeant->Next()) != NULL){

            if (pTof->getTrack() == lTrack)
            {
               nStatTof = 1;
            }
     }
     HGeantShower *pShower;
     pitShowerGeant->Reset();
     while((pShower = (HGeantShower*) pitShowerGeant->Next()) != NULL){

            if (pShower->getTrack() == lTrack)
            {
               nStatShower = 1;
            }

     }
   if (nStatMDC1 && nStatMDC2 && nStatMDC3 &&  nStatMDC4 &&  (nStatTof || nStatShower)){
       return kTRUE;
   }

   // below is a kind of debug in case of wrong data
   if (pG->getNMdcHits(0) && pG->getNMdcHits(1) && pG->getNMdcHits(2) && (pG->getNMdcHits(0)>7 || pG->getNMdcHits(1)>7 || pG->getNMdcHits(2)>7))
	{
	cout << " momentum = " << pG->getTotalMomentum();
	cout << " nStatMDC1: " << pG->getNMdcHits(0);
	cout << " nStatMDC2: "<< pG->getNMdcHits(1);
	cout << " nStatMDC3: "<< pG->getNMdcHits(2);
	cout << " nStatMDC4: "<< pG->getNMdcHits(3) <<endl;
	}

   return kFALSE;
}


// ------------- function called in the main ---------------------
Int_t analysis(TString infile = "")
{ 

 Hades* myHades = new Hades();
 gHades -> setQuietMode(2);
 
 // ****  INPUT ****
   TString inputDir  = "/hera/hades/user/przygoda/PION/GEANT/FILES/ELASTIC/";
   //TString inputDir  = "/hera/hades/user/przygoda/apr06/pluto_simFiz/PROGRAM/PION/GEANT/652/OUTPUT/";
   //TString inputDir  = "/hera/hades/user/przygoda/apr06/pluto_simFiz/PROGRAM/PION/GEANT/685/OUTPUT/";
   //TString inputDir  = "/hera/hades/user/przygoda/apr06/pluto_simFiz/PROGRAM/PION/GEANT/740/OUTPUT/";
   //TString inputDir  = "/hera/hades/user/przygoda/apr06/pluto_simFiz/PROGRAM/PION/GEANT/790/OUTPUT/";

 // **** OUTPUT ****
   TString outputDir = "/hera/hades/user/przygoda/PION/ACC/ELASTIC/";
   //TString outputDir = "/hera/hades/user/przygoda/PION/ACC/652/OUTPUT/";
   //TString outputDir = "/hera/hades/user/przygoda/PION/ACC/685/OUTPUT/";
   //TString outputDir = "/hera/hades/user/przygoda/PION/ACC/740/OUTPUT/";
   //TString outputDir = "/hera/hades/user/przygoda/PION/ACC/790/OUTPUT/";
 
 TString infile1 = infile+"1.root";
 TString outfile1 = infile+"_acc.root";
 
 HRootSource *source = new HRootSource;
 source -> setDirectory((char*)inputDir.Data());
 source -> addFile((char*)infile1.Data());
 gHades -> setDataSource(source);
 if(!gHades->init())
   {
     cout << "gHades->init() ERROR " << endl;
     exit(1);
   }
 
 /// Getting the categories
 HCategory* pGeantCat = (HCategory*)gHades->getCurrentEvent()->getCategory(catGeantKine);
 HIterator* pitGeant = 0;

 HCategory* pMdcGeantCat = gHades->getCurrentEvent()->getCategory(catMdcGeantRaw);
 HIterator* pitMdcGeant = (HIterator *)pMdcGeantCat->MakeIterator();
 pitGeant = (HIterator *)pGeantCat->MakeIterator();
 HCategory* pTofGeantCat = gHades->getCurrentEvent()->getCategory(catTofGeantRaw);
 HIterator* pitTofGeant = (HIterator *)pTofGeantCat->MakeIterator();
 HCategory* pShowerGeantCat = gHades->getCurrentEvent()->getCategory(catShowerGeantRaw);
 HIterator* pitShowerGeant = (HIterator *)pShowerGeantCat->MakeIterator();

 if(pGeantCat)
   {
   pitGeant = (HIterator *)pGeantCat->MakeIterator();
   } 
 else 
   {
     cout<<"category pGeantCat does not exist"<<endl;
      exit(1);
   }
 
 TFile ff(outputDir+outfile1,"recreate");

 HNtuple *nt = new HNtuple("ppel_acc","ppel_acc");
 nt->setFile( &ff );

 TH1F *PiM_theta_cm   = new TH1F("PiM_theta_cm" ,"PiM theta cm",53,11.5,170.5); PiM_theta_cm->Sumw2();

   //PION BEAM SCANNER
   // double pion_momentum = 656.0;
   // double pion_momentum = 690.0; 
   // double pion_momentum = 748.0;
   // double pion_momentum = 800.0;
   // -----------------------------
   // double pion_momentum = 652.0;
   // double pion_momentum = 685.0; 
   // double pion_momentum = 740.0;
   // double pion_momentum = 790.0;
 //---------------------------------------
 // double pion_momentum = 654.1; // PT
 // double pion_momentum = 683.5; // PT
 // double pion_momentum = 738.9; // PT
 // double pion_momentum = 791.1; // PT
 //---------------------------------------
  double pion_momentum = 685.8; // A 
 // double pion_momentum = 685.2; // B 
 // double pion_momentum = 688.1; // C 
 // double pion_momentum = 691.4; // D 


   double pion_energy = sqrt( pion_momentum*pion_momentum + 139.56995*139.56995 );

   TLorentzVector* proj = new TLorentzVector(0,0, pion_momentum, pion_energy); // PION BEAM momentum as above
   TLorentzVector* targ = new TLorentzVector(0,0,0,938.27231); // PROTON
/*******************************************************************************************************/
   TLorentzVector* beam = new TLorentzVector(0,0,0,0);
   *beam = *proj + *targ;

   TLorentzVector*  p = new TLorentzVector(0,0,0,0);
   TLorentzVector*  p_LAB = new TLorentzVector(0,0,0,0);
   TLorentzVector*  p_CM = new TLorentzVector(0,0,0,0);
   TLorentzVector*  pim = new TLorentzVector(0,0,0,0);
   TLorentzVector*  pim_LAB = new TLorentzVector(0,0,0,0);
   TLorentzVector*  pim_CM = new TLorentzVector(0,0,0,0);
   TLorentzVector*  ppim = new TLorentzVector(0,0,0,0);
   TLorentzVector*  ppim_CM = new TLorentzVector(0,0,0,0);
   TLorentzVector*  ppim_miss = new TLorentzVector(0,0,0,0);
   TLorentzVector*  ppim_miss_CM = new TLorentzVector(0,0,0,0);


 HGeantKine  *pKine = 0;

 TVector3    vec1, vec2, vec3, vec4;

 Float_t xMom=0.,yMom=0.,zMom=0.;
 Float_t xMom1=0.,yMom1=0.,zMom1=0.;
 Float_t xMom2=0.,yMom2=0.,zMom2=0.;
 Int_t aTrackLepton, aIDLepton, aParLepton, aMedLepton, aMechLepton;
 Int_t aTrackLepton1, aIDLepton1, aParLepton1, aMedLepton1, aMechLepton1, aGenInfo1_1 = 0;
 Int_t aTrackLepton2, aIDLepton2, aParLepton2, aMedLepton2, aMechLepton2, aGenInfo1_2 = 0;
 Float_t genInfo=0., genWeight=0.;
 Float_t genInfo1=0., genWeight1=0.;
 Float_t genInfo2=0., genWeight2=0.;
 bool protonAcc = false;
 bool pionAcc = false;
 
 Int_t  nentries = 100000;
 for (Int_t i = 0; i<nentries; i++)
   {//event loop
     if(!gHades->eventLoop(1)) 
     {
        if(i<100000) cout<<"TOO SMALL FILE "<<i<<endl;
        break;
     }
    
     if(i%1000==0) 
     cout<<" ---Current event: "<<i<<endl;
     //counter = pGeantCat->getEntries();
     //cout<<"main counter ---------"<<counter<<endl;
     //if(counter > numPart || counter < 2) continue;
     //if(counter > numPart) continue;
     

     // **** 
     // ****   genInfo
     // pi0
     // ****   1436  D+
     // ****   1438  N1440
     // ****
     //

     pipAcc = false;
     pimAcc = false;

     pitGeant->Reset();
     while((pKine = (HGeantKine*)pitGeant->Next()) != 0)
	 {

         //cout << " --------------------- " << endl;
	       pKine->getParticle( aTrackLepton, aIDLepton );
	       pKine->getCreator ( aParLepton, aMedLepton, aMechLepton );
	       pKine->getGenerator( genInfo, genInfo1, genInfo2);
	       pKine->getGenerator( genInfo, genWeight );
	       pKine->getMomentum( xMom, yMom, zMom );
           //if ( aParLepton==0 )
           //cout << " parent = " << aParLepton << " id = " << aIDLepton << " geninf = " << genInfo << endl;
       /************************************************************ PID ***************************/
	   if(aParLepton == 0 && aIDLepton == 8) 
       {
          aIDLepton1 = aIDLepton;
          xMom1 = xMom;
          yMom1 = yMom;
          zMom1 = zMom;
          aGenInfo1_1 = genInfo1;
          genWeight1 = genWeight;
	      if ( isGeantTrackInAcceptance(pKine, pitMdcGeant, pitTofGeant, pitShowerGeant) ) pipAcc = true;
          //cout << " mech = " << aMechLepton << " id = " << aIDLepton << " geninf = " << genInfo << " licznik " << licznik << " g1 " << aGenInfo1_1 << endl;
       }
	   if(aParLepton == 0 && aIDLepton == 9) // pim elastic
       {
          aIDLepton2 = aIDLepton;
          xMom2 = xMom;
          yMom2 = yMom;
          zMom2 = zMom;
          aGenInfo1_2 = genInfo1;
          genWeight2 = genWeight;
	      if ( isGeantTrackInAcceptance(pKine, pitMdcGeant, pitTofGeant, pitShowerGeant) ) pimAcc = true;
          //cout << " parent = " << aParLepton << " id = " << aIDLepton << " geninf = " << genInfo << " licznik " << licznik << " g2 " << aGenInfo1_2 << endl;
       }

    } // GeantKine loop

    vec1.SetXYZ( xMom1, yMom1, zMom1 );
    vec2.SetXYZ( xMom2, yMom2, zMom2 );
    pip_LAB->SetVectM(vec2, 139.56995); // pip
    pim_LAB->SetVectM(vec2, 139.56995); // pim

    double pip_p = pip_LAB->P();
    double pip_theta = pip_LAB->Theta()*TMath::RadToDeg();
    double pip_phi = transformPhi(pip_LAB->Phi());
    double pim_p = pim_LAB->P();
    double pim_theta = pim_LAB->Theta()*TMath::RadToDeg();
    double pim_phi = transformPhi(pim_LAB->Phi());

    *pip = *pip_LAB;
    *pim = *pim_LAB;

    *pip_CM = *pip;
    *pim_CM = *pim;

    *pippim = *pip + *pim;
    *pippim_CM = *pip + *pim;
    *pippim_miss = *beam - (*pip + *pim);
    *pippim_miss_CM = *beam - (*pip + *pim);

    pip_CM->Boost( -(*beam).BoostVector() );
    pim_CM->Boost( -(*beam).BoostVector() );
    pippim_CM->Boost( -(*beam).BoostVector() );
    pippim_miss_CM->Boost( -(*beam).BoostVector() );

    double oa = R2D * openingangle(*pip, *pim);

    double tantan = TMath::Tan( D2R * pip_theta ) * TMath::Tan( D2R * pim_theta );
    double dphi = TMath::Abs( pip_phi - pim_phi );

    (*nt)["pip_p"] = pip_p;
    (*nt)["pip_theta"] = pip_theta;
    (*nt)["pip_cm_theta"] = pip->Theta()*TMath::RadToDeg();
    (*nt)["pip_cm_costh"] = TMath::Cos( pip->Theta() );
    (*nt)["pip_phi"] = pip_phi;
    (*nt)["pip_acc"] = (int)pipAcc;

    (*nt)["pim_p"] = pim_p;
    (*nt)["pim_theta"] = pim_theta;
    (*nt)["pim_cm_theta"] = pim_CM->Theta()*TMath::RadToDeg();
    (*nt)["pim_cm_costh"] = TMath::Cos( pim_CM->Theta() );
    (*nt)["pim_phi"] = pim_phi;
    (*nt)["pim_acc"] = (int)pionAcc;

    (*nt)["tantan"] = tantan;
    (*nt)["dphi"] = dphi;
    (*nt)["oa"] = oa;

    (*nt)["pippim_miss_mass"] = pippim_miss->M();


    (*nt)["w"] = genWeight;

    nt->fill();

    PiM_theta_cm->Fill( R2D * pim_CM->Theta(), genWeight );


   } // event loop
     
 ff.cd();
 nt->Write();
 PiM_theta_cm->Write();
 ff.Close();

 delete myHades;
 
 return 1;
}

//------------------------------------- prog ----------------------------------------
#ifndef __CINT__
int main(int argc, char **argv)
{
  analysis(TString(argv[1]));
  return 0;
}
#endif

