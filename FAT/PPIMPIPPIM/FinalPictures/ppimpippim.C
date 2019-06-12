#define ppimpippim_cxx
#include "ppimpippim.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

void ppimpippim::Loop()
{
  //   In a ROOT session, you can do:
  //      Root > .L ppimpippim.C
  //      Root > ppimpippim t
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
  TFile *MyFile = new TFile("pictures.root","RECREATE");
  if ( MyFile->IsOpen() ) printf("File opened successfully\n");
  else
    printf("Cann't open a file\n");
    
  const int npt=10;
  const int nw=10;

  int pt_points=800;
  int w_points=500;
  double ptmax=1600;
  double wmax=2;

  char L_ptname[20];
  char L_wname[20];
  char L_pttitle[20];
  char L_ptname[20];
  char K0_ptname[20];
  char K0_wname[20];
  char K0_pttitle[20];
  char K0_ptname[20];
    
  TH1F *h1Lambda_pt[npt];
  TH1F *h1Lambda_w[nw];
  TH1F *h1K0_pt[npt];
  TH1F *h1K0_w[nw];

  TH2F *h2Lambda_wpt=new TH2F("h2Lambda_wpt","Rapidity vs. p_{t} for #Lambda",w_points,0,wmax,pt_points,0,ptmax);
  TH2F *h2K0_wpt=new TH2F("h2K0_wpt","Rapidity vs. p_{t} for K^{0}",w_points,0,wmax,pt_points,0,ptmax);

  for(int i=0;i<npt;i++)
    {
      sprintf(L_pt_name,"Lambda_pt_%d_%d",ptmax*(double)i/npt,ptmax*(double)(i+1)/npt);
      sprintf(L_pt_title,"P^{T}_{#Lambda(1116)} (%d, %d)",ptmax*(double)i/npt,ptmax*(double)(i+1)/npt);

      sprintf(K0_pt_name,"K0_pt_%d_%d",ptmax*(double)i/npt,ptmax*(double)(i+1)/npt);
      sprintf(K0_pt_title,"P^{T}_{K^{0}} (%d, %d)",ptmax*(double)i/npt,ptmax*(double)(i+1)/npt);
      
      h1Lambda_pt[i]=new TH1F(L_pt_name,L_pt_title,pt_points,0,ptmax);
      h1K0_pt[i]=new TH1F(K0_pt_name,K0_pt_title,pt_points,0,ptmax);
    }

  for(int i=0;i<npt;i++)
    {
      sprintf(L_w_name,"Lambda_w_%d_%d",wmax*(double)i/nw,wmax*(double)(i+1)/nw);
      sprintf(L_w_title,"P^{T}_{#Lambda(1116)} (%d, %d)",wmax*(double)i/nw,wmax*(double)(i+1)/nw);

      sprintf(K0_w_name,"K0_w_%d_%d",wmax*(double)i/nw,wmax*(double)(i+1)/nw);
      sprintf(K0_w_title,"P^{T}_{K^{0}} (%d, %d)",wmax*(double)i/nw,wmax*(double)(i+1)/nw);
      
      h1Lambda_w[i]=new TH1F(L_w_name,L_w_title,w_points,0,wmax);
      h1K0_w[i]=new TH1F(K0_w_name,K0_w_title,w_points,0,wmax);
    }

  //Main loop
  if (fChain == 0) return;
  Long64_t nentries = fChain->GetEntriesFast();
  Long64_t nbytes = 0, nb = 0;

  for (Long64_t jentry=0; jentry<nentries;jentry++)
    {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      if (Cut(ientry) < 0
	  || isBest_new!=1
	  || mlp_output<0.4
	  )
	continue;

      if(jentry%100==0)
	cout<<(double)jentry/nentries<<"% "<<endl;

      h2Lambda_wpt->Fill(lambda_w,lambda_pt);
      h2K0_wpt->Fill(k0_w,k0_pt);
    }
  //End of main loop
  
  //Save histograms
  h2K0_wpt->Write();
  h2Lambda_wpt->Write();

  for(int i=0;i<pt_points;i++)
    {
      h1K0_pt[i]->Write();
      h1Lambda_pt[i]->Write();
    }
  for(int i=0;i<w_points;i++)
    {
      h1K0_w[i]->Write();
      h1Lambda_w[i]->Write();
    }
    
  MyFile->Close();
}
