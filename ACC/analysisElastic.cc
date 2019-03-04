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


// ------------- function called in the main ---------------------
Int_t analysis(TString infile = "")
{ 

  Hades* myHades = new Hades();
  gHades -> setQuietMode(2);
 
  // ****  INPUT ****
  TString inputDir  = "";

  // **** OUTPUT ****
  TString outputDir = "/lustre/nyx/hades/user/knowakow/PP/ACC/FILES/ver1/";
 
  TString inputFile = infile;
  
  inputDir=inputFile;
  inputDir.Resize(inputFile.Last('/')+1);
    
  inputFile=inputFile(inputFile.Last('/')+1,inputFile.Length()-inputFile.Last('/')-1);
  TString outfile1 = inputFile+"_acc.root";
  
  HRootSource *source = new HRootSource;
  source -> setDirectory((char*)inputDir.Data());
  source -> addFile((char*)inputFile.Data());
  gHades -> setDataSource(source);

  if(!gHades->init())
    {
      cout << "gHades->init() ERROR " << endl;
      exit(1);
    }
 
  /// Getting the categories
  HCategory* pGeantCat = (HCategory*)gHades->getCurrentEvent()->getCategory(catGeantKine);
  HIterator* pitGeant = (HIterator *)pGeantCat->MakeIterator();

  if(!pGeantCat)
    {
      cout<<"category pGeantCat does not exist"<<endl;
      exit(1);
    }
 
  TFile ff(outputDir+outfile1,"recreate");

  HNtuple *nt = new HNtuple("pippim","pippim");
  nt->setFile( &ff );

  TH1F *PiM_theta_cm   = new TH1F("PiM_theta_cm" ,"PiM theta cm",53,11.5,170.5); PiM_theta_cm->Sumw2();
  TH1F *beam_mom   = new TH1F("beam_mom" ,"beam_mom theta cm",1000,500,1000); beam_mom->Sumw2();

  //PION BEAM SCANNER
  //---------------------------------------
  double pion_momentum = 654.1; // PT
  // double pion_momentum = 683.5; // PT
  // double pion_momentum = 738.9; // PT
  // double pion_momentum = 791.1; // PT
  //---------------------------------------
  // double pion_momentum = 685.8; // A 
  // double pion_momentum = 685.2; // B 
  // double pion_momentum = 688.1; // C 
  // double pion_momentum = 691.4; // D 


  double pion_energy = sqrt( pion_momentum*pion_momentum + 139.56995*139.56995 );

  TLorentzVector* proj = new TLorentzVector(0,0, pion_momentum, pion_energy); // PION BEAM momentum as above
  TLorentzVector* targ = new TLorentzVector(0,0,0,938.27231); // PROTON
  /*******************************************************************************************************/
  TLorentzVector* beam = new TLorentzVector(0,0,0,0);
  TLorentzVector* beam_PT = new TLorentzVector(0,0,0,0);
  *beam = *proj + *targ;

  TLorentzVector*  n = new TLorentzVector(0,0,0,0);
  TLorentzVector*  n_LAB = new TLorentzVector(0,0,0,0);
  TLorentzVector*  n_CM = new TLorentzVector(0,0,0,0);
  TLorentzVector*  pim = new TLorentzVector(0,0,0,0);
  TLorentzVector*  pim_LAB = new TLorentzVector(0,0,0,0);
  TLorentzVector*  pim_CM = new TLorentzVector(0,0,0,0);
  TLorentzVector*  pip = new TLorentzVector(0,0,0,0);
  TLorentzVector*  pip_LAB = new TLorentzVector(0,0,0,0);
  TLorentzVector*  pip_CM = new TLorentzVector(0,0,0,0);
  TLorentzVector*  pippim_miss = new TLorentzVector(0,0,0,0);


  HGeantKine  *pKine = 0;

  TVector3    vec1, vec2, vec3, vec4;

  Float_t xMom=0.,yMom=0.,zMom=0.;
  Float_t xMom1=0.,yMom1=0.,zMom1=0.;
  Float_t xMom2=0.,yMom2=0.,zMom2=0.;
  Float_t xMom3=0.,yMom3=0.,zMom3=0.;
  Int_t aTrackLepton, aIDLepton, aParLepton, aMedLepton, aMechLepton;
  Int_t aTrackLepton1, aIDLepton1, aParLepton1, aMedLepton1, aMechLepton1, aGenInfo1_1 = 0;
  Int_t aTrackLepton2, aIDLepton2, aParLepton2, aMedLepton2, aMechLepton2, aGenInfo1_2 = 0;
  Int_t aTrackLepton3, aIDLepton3, aParLepton3, aMedLepton3, aMechLepton3, aGenInfo1_3 = 0;
  Float_t genInfo=0., genWeight=0.;
  Float_t genInfo1=0., genWeight1=0.;
  Float_t genInfo2=0., genWeight2=0.;
  Float_t genWeight3=0.;
 
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
	  if(aParLepton == 0 && aIDLepton == 8) // pip
	    {
	      aIDLepton1 = aIDLepton;
	      xMom1 = xMom;
	      yMom1 = yMom;
	      zMom1 = zMom;
	      aGenInfo1_1 = genInfo1;
	      genWeight1 = genWeight;
	      //cout << " mech = " << aMechLepton << " id = " << aIDLepton << " geninf = " << genInfo << " licznik " << licznik << " g1 " << aGenInfo1_1 << endl;
	    }
	  if(aParLepton == 0 && aIDLepton == 9) // pim 
	    {
	      aIDLepton2 = aIDLepton;
	      xMom2 = xMom;
	      yMom2 = yMom;
	      zMom2 = zMom;
	      aGenInfo1_2 = genInfo1;
	      genWeight2 = genWeight;
	      //cout << " parent = " << aParLepton << " id = " << aIDLepton << " geninf = " << genInfo << " licznik " << licznik << " g2 " << aGenInfo1_2 << endl;
	    }
	  if(aParLepton == 0 && aIDLepton == 13) // n
	    {
	      aIDLepton3 = aIDLepton;
	      xMom3 = xMom;
	      yMom3 = yMom;
	      zMom3 = zMom;
	      aGenInfo1_3 = genInfo1;
	      genWeight3 = genWeight;
	      //cout << " parent = " << aParLepton << " id = " << aIDLepton << " geninf = " << genInfo << " licznik " << licznik << " g2 " << aGenInfo1_2 << endl;
	    }


	} // GeantKine loop

      vec1.SetXYZ( xMom1, yMom1, zMom1 );
      vec2.SetXYZ( xMom2, yMom2, zMom2 );
      vec3.SetXYZ( xMom3, yMom3, zMom3 );
      pip_LAB->SetVectM(vec1, 139.56995); // pip
      pim_LAB->SetVectM(vec2, 139.56995); // pim
      n_LAB->SetVectM(vec3, 939.56563); // n

      double pip_p = pip_LAB->P();
      double pip_theta = pip_LAB->Theta()*TMath::RadToDeg();
      double pip_phi = transformPhi(pip_LAB->Phi());
      double pim_p = pim_LAB->P();
      double pim_theta = pim_LAB->Theta()*TMath::RadToDeg();
      double pim_phi = transformPhi(pim_LAB->Phi());
      double n_p = pim_LAB->P();
      double n_theta = pim_LAB->Theta()*TMath::RadToDeg();
      double n_phi = transformPhi(pim_LAB->Phi());

      *pip = *pip_LAB;
      *pim = *pim_LAB;
      *n = *n_LAB;

      *beam_PT = *pip_LAB + *pim_LAB + *n_LAB;

      beam_mom->Fill( beam_PT->P() );

      *pip_CM = *pip;
      *pim_CM = *pim;
      *n_CM = *n;

      *pippim_miss = *beam_PT - (*pip + *pim);

      pip_CM->Boost( -(*beam_PT).BoostVector() );
      pim_CM->Boost( -(*beam_PT).BoostVector() );
      n_CM->Boost( -(*beam_PT).BoostVector() );

      double oa = R2D * openingangle(*pip, *pim);

      double tantan = TMath::Tan( D2R * pip_theta ) * TMath::Tan( D2R * pim_theta );
      double dphi = TMath::Abs( pip_phi - pim_phi );

      (*nt)["pip_p"] = pip_p;
      (*nt)["pip_theta"] = pip_theta;
      (*nt)["pip_cm_theta"] = pip_CM->Theta()*TMath::RadToDeg();
      (*nt)["pip_cm_costh"] = TMath::Cos( pip_CM->Theta() );
      (*nt)["pip_phi"] = pip_phi;

      (*nt)["pim_p"] = pim_p;
      (*nt)["pim_theta"] = pim_theta;
      (*nt)["pim_cm_theta"] = pim_CM->Theta()*TMath::RadToDeg();
      (*nt)["pim_cm_costh"] = TMath::Cos( pim_CM->Theta() );
      (*nt)["pim_phi"] = pim_phi;

      (*nt)["n_p"] = n_p;
      (*nt)["n_theta"] = n_theta;
      (*nt)["n_cm_theta"] = n_CM->Theta()*TMath::RadToDeg();
      (*nt)["n_cm_costh"] = TMath::Cos( n_CM->Theta() );
      (*nt)["n_phi"] = n_phi;

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
  beam_mom->Write();
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

