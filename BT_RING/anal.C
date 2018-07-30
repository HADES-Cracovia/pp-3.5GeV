#define anal_cxx
#include "anal.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

void anal::Loop()
{
//   In a ROOT session, you can do:
//      Root > .L anal.C
//      Root > anal t
//      Root > t.GetEntry(12); // Fill t data members with entry number 12
//      Root > t.Show();       // Show values of entry 12
//      Root > t.Show(16);     // Read and show values of entry 16
//      Root > t.Loop();       // Loop on all entries
//

//     This is the loop skeleton where:
//    jentry is the global entry number in the chain
//    ientry is the entry number in the current Tree
//  Note that the argument to GetEntry must be:
//    jentry for TChain::GetEntry
//    ientry for TTree::GetEntry and TBranch::GetEntry
//
//       To read only selected branches, Insert statements like:
// METHOD1:
//    fChain->SetBranchStatus("*",0);  // disable all branches
//    fChain->SetBranchStatus("branchname",1);  // activate branchname
// METHOD2: replace line
//    fChain->GetEntry(jentry);       //read all branches
//by  b_branchname->GetEntry(ientry); //read only this branch

// declaration
  double openingangle(const TLorentzVector a, const TLorentzVector b)
  {
    return TMath::ACos( (a.Px()*b.Px() + a.Py()*b.Py() +  a.Pz()*b.Pz() ) / ( a.Vect().Mag() * b.Vect().Mag() ) );
  } 
  
  Float_t mass2_el=0.0;
  Float_t mass2_ep=0.0;
  Float_t mass2_p=0.0;
  Float_t mass_el=0.511;
  Float_t mass_p=939;
  Float_t deg2rad=180./3.14;
  //
  TLorentzVector el;
  TLorentzVector ep;
  TLorentzVector p;
  TLorentzVector beamT;
  TLorentzVector dilep;
  TLorentzVector p_miss;

  Int_t icount;
  //histograms
  TH1F *con=new TH1F("con","conditions",10,0,10);
  TH1F *inv_mass=new TH1F("inv_mass","invariant mass e+e-",50,0,1000);
  TH1F *pee_miss=new TH1F("pee_miss","pee missing mass",60,400,1600);
  TH1F *mom=new TH1F("mom","lepton momentum",200,-1200,1200);
  TH2F* p_mom_beta=new TH2F("p_mom_beta","",400,0,4000,100,0,1.2);
  TH2F* l_mom_beta=new TH2F("l_mom_beta","",400,-2000,2000,100,0,1.2);

  beamT.SetPxPyPzE(0,0,4537,5378); //in MeV E_kin=3500 , mass=939  
    
  if (fChain == 0) return;
  Long64_t nentries = fChain->GetEntriesFast();
  Long64_t nbytes = 0, nb = 0;
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;
    //if (Cut(ientry) < 0) continue;
    icount=0;

    el.SetXYZM(em_p*sin(em_theta/deg2rad)*cos(em_phi/deg2rad),em_p*sin(em_theta/deg2rad)*sin(em_phi/deg2rad),em_p*cos(em_theta/deg2rad),mass_el);
    ep.SetXYZM(ep_p*sin(ep_theta/deg2rad)*cos(ep_phi/deg2rad),ep_p*sin(ep_theta/deg2rad)*sin(ep_phi/deg2rad),ep_p*cos(ep_theta/deg2rad),mass_el);
    p.SetXYZM(p_p_corr_p*sin(p_theta/deg2rad)*cos(p_phi/deg2rad),p_p_corr_p*sin(p_theta/deg2rad)*sin(p_phi/deg2rad),p_p_corr_p*cos(p_theta/deg2rad),mass_p);

    double oa=180/3.14*openingangle(el,ep);
    //cout<<oa<<endl;    
    // ring finder/backtracking
    if (!(em_isring>0 || (em_isBT>=0 && em_btPadsRing>1))) continue;
    if (!(ep_isring>0 || (ep_isBT>=0 && ep_btPadsRing>1))) continue;
    //if (!(em_isring>0)) continue;
    //if (!(ep_isring>0)) continue;
    con->Fill(icount); 
    if (!isBest>0) continue;
    icount++;
    con->Fill(icount);
    mass2_el=em_p*em_p*(1/em_beta_new/em_beta_new -1);
    mass2_ep=ep_p*ep_p*(1/ep_beta_new/ep_beta_new -1);
    // electron mass 
    if(ep_beta_new>.85 && em_beta_new>.85 && ep_p<2000 && em_p<2000 ) {
      icount++;
      con->Fill(icount);
      // proton mass
      mass2_p=p_p*p_p*(1/p_beta_new/p_beta_new -1); 
      if(mass2_p>450000.)
	{
	  icount++;
	  con->Fill(icount);
	  p_mom_beta->Fill(p_p,p_beta_new);
	  
	  l_mom_beta->Fill(-em_p,em_beta_new);
	  l_mom_beta->Fill(ep_p,ep_beta_new);
	  if(/*em_oa_lept>4.0 && ep_oa_lept>4.0*/oa>4)
	    {
	      icount++;
	      con->Fill(icount);
	      mom->Fill(em_p*(-1.));
	      mom->Fill(ep_p);
	      
	      dilep=el+ep;
	      inv_mass->Fill(dilep.M());

	      //above pi0
	      if(dilep.M()>140.)
		{
		  p_miss=beamT-dilep-p;
		  pee_miss->Fill(p_miss.M());
		} 
	    }
	} 
    } 
  } 
}
