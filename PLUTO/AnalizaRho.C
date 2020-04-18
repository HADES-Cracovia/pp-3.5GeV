#define AnalizaRho_cxx
#include "AnalizaRho.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

void AnalizaRho::Loop()
{
//   In a ROOT session, you can do:
//      Root > .L AnalizaRho.C
//      Root > AnalizaRho t
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

  if (fChain == 0) return;
  TFile *MyFile = new TFile("PicturesRho.root","RECREATE");
  TH1F *hparticlePID=new TH1F("hparticlePID","PID for all particles;pid",50,0,50);
  TH1F *hparticlePID_L1520=new TH1F("hparticlePID_L1520","PID for particles from L1520;pid",50,0,50);
  TH1F *hparticlePID_reaction=new TH1F("hparticlePID_reaction","PID for particles from primary vertex;pid",50,0,50);
  TH1F *hMInvLPipPim=new TH1F("hMInvLPipPim","Invariant mass #Lambda #pi^{+} #pi^{-};M^{inv}_{#pi^{+} #pi^{-} #Lambda}[GeV]",500,1,2);
  TH1F *hMInvPipPim=new TH1F("hMInvPipPim","Invariant mass #pi^{+} #pi^{-};M^{inv}_{#pi^{+} #pi^{-}}[GeV]",500,0,0.5);
    
  Long64_t nentries = fChain->GetEntries();

  Long64_t nbytes = 0, nb = 0;
  for (Long64_t jentry=0; jentry<nentries;jentry++)
    {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      if(jentry % 1000==0)
	{
	  cout<<(100.0*jentry)/(1.0*nentries)<<"%"<<endl;
	  //cout<<"nentries= "<<nentries<<endl;
	  //cout<<"jentry= "<<jentry<<endl;
	  //cout<<"Particles_: "<<Particles_<<endl;
	  //cout<<"Npart: "<< Npart<<endl;
	}
      TLorentzVector Lambda;
      TLorentzVector Pip;
      TLorentzVector Pim;
      TLorentzVector PipPim;
      TLorentzVector LPipPim;
      for(int n=0;n<Particles_;n++)//loop over particles in event
	{
	  hparticlePID->Fill(Particles_pid[n]);
	  //cout<<Particles_pid[n]<<endl;
	  if(Particles_parentId[n]==72)//fromL(1520)
	    {
	      hparticlePID_L1520->Fill(Particles_pid[n]);
	      if(Particles_pid[n]==8)//pi+
		Pip.SetPxPyPzE(Particles_fP_fX[n],Particles_fP_fY[n],Particles_fP_fX[n],Particles_fE[n]);
	      if(Particles_pid[n]==9)//pi+
		Pim.SetPxPyPzE(Particles_fP_fX[n],Particles_fP_fY[n],Particles_fP_fX[n],Particles_fE[n]);
	      if(Particles_pid[n]==18)//pi+
		Lambda.SetPxPyPzE(Particles_fP_fX[n],Particles_fP_fY[n],Particles_fP_fX[n],Particles_fE[n]);
	    }
	  if(Particles_parentId[n]==14014)//primary vertex
	    {
	      hparticlePID_reaction->Fill(Particles_pid[n]);
	    }
	}
      PipPim=Pip+Pim;
      LPipPim=Pip+Pim+Lambda;
      hMInvLPipPim->Fill(LPipPim.M());
      hMInvPipPim->Fill(PipPim.M());
    }
  MyFile->Write();
  MyFile->Close();
}
