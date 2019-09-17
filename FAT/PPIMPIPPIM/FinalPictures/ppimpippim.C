#define ppimpippim_cxx
#include "ppimpippim.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <iostream>
#include <TKey.h>
#include <TMath.h>
#include <TObject.h>
#include <TH1F.h>
#include <TF1.h>

void normalize(TH1* hist)
{
  for (Int_t j=1; j<hist->GetNbinsX()+1; ++j)
    {
      double scale=1.0/(3.13 * TMath::Power(10,8)) *1000; /*to get micro barn*/
      hist->SetBinContent(j, hist->GetBinContent(j) / hist->GetBinWidth(j) *scale);
      //hist->SetBinError( j, TMath::Sqrt( hist->GetBinContent(j) ) );
      hist->SetBinError( j, hist->GetBinError(j) / hist->GetBinWidth(j) * scale );
      hist->GetYaxis()->SetTitle("#frac{dN}{dE} [#frac{#mu b}{MeV}]");
    }
  
}


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
  if ( MyFile->IsOpen() )
    printf("File opened successfully\n");
  else
    printf("Cann't open a file\n");

  TDirectory *MyDirectory=new TDirectory("finalHistograms","Final Histograms destination");
  
  const int npt=7;
  const int nw=7;

  int pt_points=500;
  int w_points=500;
  double ptmax=1000;
  double wmax=2;

  char L_pt_name[50];
  char L_w_name[50];
  char L_pt_title[50];
  char L_w_title[50];

  char K0_pt_name[50];
  char K0_w_name[50];
  char K0_pt_title[50];
  char K0_w_title[50];
    
  TH1F* h1Lambda_pt[npt];
  TH1F* h1Lambda_pt_k0cut[nw];
  TH1F* h1Lambda_w[nw];
  TH1F* h1Lambda_w_k0cut[nw];
  TH1F* h1K0_pt[npt];
  TH1F* h1K0_pt_Lcut[nw];
  TH1F* h1K0_w[nw];
  TH1F* h1K0_w_Lcut[nw];

  TH2F *h2Lambda_wpt=new TH2F("h2Lambda_wpt","Rapidity vs. p_{t} for #Lambda",w_points,0,wmax,pt_points,0,ptmax);
  TH2F *h2K0_wpt=new TH2F("h2K0_wpt","Rapidity vs. p_{t} for K^{0}",w_points,0,wmax,pt_points,0,ptmax);
  TH2F *h2_m_inv=new TH2F("h2_m_inv","M^{inv}_{#pi^{-} #pi^{+}} vs. M^{inv}_{#pi^{-} p}",300,1000,1600,300,200,800); 
  
  TH1F *h1Lambda_pt_all=new TH1F("h1lambda_pt_all","p_{t} for #Lambda",pt_points,0,ptmax);
  TH1F *h1K0_pt_all=new TH1F("h1K0_pt_all","p_{t} for K^{0}",pt_points,0,ptmax);

  TH1F *h1Lambda_w_all=new TH1F("h1lambda_w_all","w for #Lambda",w_points,0,wmax);
  TH1F *h1K0_w_all=new TH1F("h1K0_w_all","w for K^{0}",w_points,0,wmax);

  TH1F *h1Lambda_m_all=new TH1F("h1Lambda_m_all","M^{inv}_{p #pi^{-}}",1000,1000,2000);
  TH1F *h1k0_m_all=new TH1F("h1k0_m_all","M^{inv}_{#pi^{+} #pi^{-}}",1000,200,1200);
  TH1F *h1Lambda_m_all_cut=new TH1F("h1Lambda_m_all_cut","M^{inv}_{p #pi^{-}}",1000,1000,2000);
  TH1F *h1k0_m_all_cut=new TH1F("h1k0_m_all_cut","M^{inv}_{#pi^{+} #pi^{-}}",1000,200,1200);

  TF1* fVoigt_bg_L1115= new TF1("fVoigt_bg_L1115","[0]*TMath::Voigt(x-[1],[2],[3])+pol5(4)",1090.00,1156.67);
  TF1* fVoigt_L1115= new TF1("fVoigt_L1115","[0]*TMath::Voigt(x-[1],[2],[3])",1090.00,1156.67);
  TF1* fbg_L1115= new TF1("fbg_L1115","pol5(0)",1090.00,1156.67);

  TF1* fVoigt_bg_K0= new TF1("fVoigt_bg_K0","[0]*TMath::Voigt(x-[1],[2],[3])+pol5(4)",250,700);
  TF1* fVoigt_K0= new TF1("fVoigt_K0","[0]*TMath::Voigt(x-[1],[2],[3])",250,700);
  TF1* fbg_K0= new TF1("fbg_K0","pol5(0)",250,700);

  cout<<endl<<"Init pt histograms"<<endl;
  for(int i=0;i<npt;i++)
    {
      sprintf(L_pt_name,"Lambda_pt_%.0f_%.0f",ptmax*(double)i/npt,ptmax*(double)(i+1)/npt);
      sprintf(L_pt_title,"P^{T}_{#Lambda(1116)} (%.2f, %.2f)",ptmax*(double)i/npt,ptmax*(double)(i+1)/npt);

      sprintf(K0_pt_name,"K0_pt_%.0f_%.0f",ptmax*(double)i/npt,ptmax*(double)(i+1)/npt);
      sprintf(K0_pt_title,"P^{T}_{K^{0}} (%.2f, %.2f)",ptmax*(double)i/npt,ptmax*(double)(i+1)/npt);

      cout<<"pt_name: "<<L_pt_name<<" "<<K0_pt_name<<endl;
      cout<<"pt_title: "<<L_pt_title<<" "<<K0_pt_title<<endl;
      
      h1Lambda_pt[i]=new TH1F(L_pt_name,L_pt_title,pt_points,1000,2000);
      h1K0_pt[i]=new TH1F(K0_pt_name,K0_pt_title,pt_points,0,1000);

      sprintf(L_pt_name,"Lambda_pt_%.2f_%.2f_k0cut",ptmax*(double)i/npt,ptmax*(double)(i+1)/npt);
      sprintf(L_pt_title,"pt_{#Lambda(1116)}(%.2f, %.2f) K^{0} cut;M^{inv}_{p #pim^{-}}[MeV]",ptmax*(double)i/npt,ptmax*(double)(i+1)/npt);
      
      sprintf(K0_pt_name,"K0_pt_%.2f_%.2f_Lcut",ptmax*(double)i/npt,ptmax*(double)(i+1)/npt);
      sprintf(K0_pt_title,"pt_{K^{0}} (%.2f, %.2f) L(1116) cut; M^{inv}_{p #pim^{-}}[MeV]",ptmax*(double)i/npt,ptmax*(double)(i+1)/npt);

      h1Lambda_pt_k0cut[i]=new TH1F(L_pt_name,L_pt_title,pt_points,1000,2000);
      h1K0_pt_Lcut[i]=new TH1F(K0_pt_name,K0_pt_title,pt_points,0,1000);
    }
  /*
    for(int i=0;i<npt;i++)
    {
    cout<<"i: "<<i<<endl;
    h1Lambda_pt[i]->Print();
    h1K0_pt[i]->Print();
    }
  */
  cout<<endl<<"Init w histograms"<<endl;
  for(int i=0;i<nw;i++)
    {
      sprintf(L_w_name,"Lambda_w_%.2f_%.2f",wmax*(double)i/nw,wmax*(double)(i+1)/nw);
      sprintf(L_w_title,"w_{#Lambda(1116)} (%.2f, %.2f)",wmax*(double)i/nw,wmax*(double)(i+1)/nw);
      
      sprintf(K0_w_name,"K0_w_%.2f_%.2f",wmax*(double)i/nw,wmax*(double)(i+1)/nw);
      sprintf(K0_w_title,"w_{K^{0}} (%.2f, %.2f)",wmax*(double)i/nw,wmax*(double)(i+1)/nw);

      cout<<"w_name: "<<L_w_name<<" "<<K0_w_name<<endl;
      cout<<"w_title: "<<L_w_title<<" "<<K0_w_title<<endl;
      
      h1Lambda_w[i]=new TH1F(L_w_name,L_w_title,w_points,1000,2000);
      h1K0_w[i]=new TH1F(K0_w_name,K0_w_title,w_points,0,1000);
      
      sprintf(L_w_name,"Lambda_w_%.2f_%.2f_k0cut",wmax*(double)i/nw,wmax*(double)(i+1)/nw);
      sprintf(L_w_title,"w_{#Lambda(1116)}(%.2f, %.2f) K^{0} cut; M^{inv}_{p #pim^{-}}[MeV]",wmax*(double)i/nw,wmax*(double)(i+1)/nw);
      
      sprintf(K0_w_name,"K0_w_%.2f_%.2f_Lcut",wmax*(double)i/nw,wmax*(double)(i+1)/nw);
      sprintf(K0_w_title,"w_{K^{0}} (%.2f, %.2f) L(1116) cut; M^{inv}_{p #pim^{-}}[MeV]",wmax*(double)i/nw,wmax*(double)(i+1)/nw);

      h1Lambda_w_k0cut[i]=new TH1F(L_w_name,L_w_title,w_points,1000,2000);
      h1K0_w_Lcut[i]=new TH1F(K0_w_name,K0_w_title,w_points,0,1000);
          
    }
  /*
    for(int ii=0;ii<nw;ii++)
    {
    cout<<"ii: "<<ii<<endl;
    h1K0_w[ii]->Print();
    h1Lambda_w[ii]->Print();
    //h1K0_w_Lcut[ii]->Print();
    //h1Lambda_w_k0cut[ii]->Print();
    }
  */

  //Main loop
  cout<<"start main loop"<<endl;
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
	  || mlp_output<0.58
	  || miss_mass_kp<1077
	  )
	continue;

      if(jentry%10000==0)
	std::cout<<(double)jentry/nentries*100<<"% "<<endl;
      
      h2Lambda_wpt->Fill(lambda_w,lambda_pt);
      h2K0_wpt->Fill(k0_w,k0_pt);

      h1Lambda_pt_all->Fill(lambda_pt);
      h1K0_pt_all->Fill(k0_pt);
      
      h1Lambda_w_all->Fill(lambda_w);
      h1K0_w_all->Fill(k0_w);

      h1k0_m_all->Fill(m_inv_pip_pim);
      h1Lambda_m_all->Fill(m_inv_p_pim);
      
      h2_m_inv->Fill(m_inv_p_pim,m_inv_pip_pim);

      if(m_inv_pip_pim<500 && m_inv_pip_pim>480)
	h1Lambda_m_all_cut->Fill(m_inv_p_pim);
      if(m_inv_p_pim<1120 && m_inv_p_pim>1110)
	h1k0_m_all_cut->Fill(m_inv_pip_pim);
	
      //fill histograms for pt and w slides
      for(int i=0; i<npt; i++)
	{
	  //cout<<"m_inv_p_pim: "<<m_inv_p_pim<<endl;
	  //cout<<"npt:"<<npt<<" i:"<<i<<endl;
	  double pt_min=ptmax*((double)i)/npt;
	  double pt_max=ptmax*((double)(i+1))/npt;
	  //h1Lambda_pt[i]->Print();
	  //h1K0_pt[i]->Print();
	  
	  if(lambda_pt<pt_max && lambda_pt>pt_min)
	    {
	      h1Lambda_pt[i]->Fill(m_inv_p_pim);
	      if(m_inv_pip_pim<510 && m_inv_pip_pim>490)
		h1Lambda_pt_k0cut[i]->Fill(m_inv_p_pim);
	    }
	  //cout<<"aab"<<endl;
	  if(k0_pt<pt_max && k0_pt>pt_min)
	    {
	      h1K0_pt[i]->Fill(m_inv_pip_pim);
	      if(m_inv_p_pim<1120 && m_inv_p_pim>1110)
		h1K0_pt_Lcut[i]->Fill(m_inv_pip_pim);
	    }
	}
      for(int i=0;i<nw;i++)
	{
	  
	  double w_min=(double)i/nw*wmax;
	  double w_max=(double)(i+1)/nw*wmax;
	  if(lambda_w<w_max && lambda_w>w_min)
	    {
	      h1Lambda_w[i]->Fill(m_inv_p_pim);
	      if(m_inv_pip_pim<500 && m_inv_pip_pim>480)
		h1Lambda_w_k0cut[i]->Fill(m_inv_p_pim);
	    }
	  if(k0_w<w_max && k0_w>w_min)
	    {
	      h1K0_w[i]->Fill(m_inv_pip_pim);
	      if(m_inv_p_pim<1120 && m_inv_p_pim>1110)
		h1K0_w_Lcut[i]->Fill(m_inv_pip_pim);
	    }
	}
    }
  //End of main loop
  cout<<"End of main loop"<<endl;

  //Save histograms
    
  h2K0_wpt->Write();
  h2Lambda_wpt->Write();

  h2Lambda_wpt->Write();
  h2K0_wpt->Write();

  h2_m_inv->Write();
  
  h1Lambda_pt_all->Write();
  h1K0_pt_all->Write();

  h1Lambda_w_all->Write();
  h1K0_w_all->Write();

  h1Lambda_m_all->Write();
  h1k0_m_all->Write();

  h1Lambda_m_all_cut->Write();
  h1k0_m_all_cut->Write();

  for(int i=0;i<npt;i++)
    {
      h1K0_pt[i]->Write();
      h1Lambda_pt[i]->Write();
      h1K0_pt_Lcut[i]->Write();
      h1Lambda_pt_k0cut[i]->Write();
    }
  
  for(int i=0;i<nw;i++)
    {
      h1K0_w[i]->Write();
      h1Lambda_w[i]->Write();
      h1K0_w_Lcut[i]->Write();
      h1Lambda_w_k0cut[i]->Write();
    }

  //Fit histograms

  //Lambda 1115
  fVoigt_bg_L1115->SetParameters(0.00951603,1114.13,1.44053,2.99181,-0.122525,2.18901e-5,7.87049e-8,7.05739e-11,1.59077e-14,-7.04671e-17);
  fVoigt_bg_L1115->SetParLimits(3,0,1);
  fVoigt_bg_L1115->SetParLimits(2,0,2);
  fVoigt_bg_L1115->SetParLimits(1,1112,1117);
  fVoigt_bg_L1115->SetRange(1106,1126);
  h1Lambda_m_all_cut->Fit(fVoigt_bg_L1115,"R");
  h1Lambda_m_all_cut->Fit(fVoigt_bg_L1115,"R");
  fVoigt_bg_L1115->SetRange(1099,1134);
  h1Lambda_m_all_cut->Fit(fVoigt_bg_L1115,"R");
  fVoigt_bg_L1115->SetRange(1094,1153);
  h1Lambda_m_all_cut->Fit(fVoigt_bg_L1115,"R");
  fVoigt_bg_L1115->SetRange(1085,1160);
  h1Lambda_m_all_cut->Fit(fVoigt_bg_L1115,"R");
  h1Lambda_m_all_cut->Fit(fVoigt_bg_L1115,"R");

  fbg_L1115->SetParameters(fVoigt_bg_L1115->GetParameter(4),fVoigt_bg_L1115->GetParameter(5),fVoigt_bg_L1115->GetParameter(6),fVoigt_bg_L1115->GetParameter(7),fVoigt_bg_L1115->GetParameter(8),fVoigt_bg_L1115->GetParameter(9));
  fVoigt_L1115->SetParameters(fVoigt_bg_L1115->GetParameter(0),fVoigt_bg_L1115->GetParameter(1),fVoigt_bg_L1115->GetParameter(2),fVoigt_bg_L1115->GetParameter(3));
  fVoigt_L1115->Draw("same");
  fVoigt_L1115->SetLineColor(kGreen);
  fbg_L1115->Draw("same");
  fbg_L1115->SetLineColor(kBlue);
  h1Lambda_m_all_cut->Draw();
  
  //K0
  fVoigt_bg_K0->SetParameters(0.010109,491.385,0.031725,11.091,0.0043912,-1.30034e-6,-1.26503e-8,-2.57051e-11,-1.37871e-14,1.23134e-16);
  fVoigt_bg_K0->SetParLimits(3,0,10);
  fVoigt_bg_K0->SetParLimits(2,0,10);
  fVoigt_bg_K0->SetParLimits(1,490,500);
  fVoigt_bg_K0->SetRange(485,512);
  h1k0_m_all_cut->Fit(fVoigt_bg_K0,"R");
  h1k0_m_all_cut->Fit(fVoigt_bg_K0,"R");
  fVoigt_bg_K0->SetRange(449,535);
  h1k0_m_all_cut->Fit(fVoigt_bg_K0,"R");
  fVoigt_bg_K0->SetRange(436,554);
  h1k0_m_all_cut->Fit(fVoigt_bg_K0,"R");
  fVoigt_bg_K0->SetRange(412,597);
  h1k0_m_all_cut->Fit(fVoigt_bg_K0,"R");
  h1k0_m_all_cut->Fit(fVoigt_bg_K0,"R");
  
  fbg_K0->SetParameters(fVoigt_bg_K0->GetParameter(4),fVoigt_bg_K0->GetParameter(5),fVoigt_bg_K0->GetParameter(6),fVoigt_bg_K0->GetParameter(7),fVoigt_bg_K0->GetParameter(8),fVoigt_bg_K0->GetParameter(9));
  fVoigt_K0->SetParameters(fVoigt_bg_K0->GetParameter(0),fVoigt_bg_K0->GetParameter(1),fVoigt_bg_K0->GetParameter(2),fVoigt_bg_K0->GetParameter(3));
  fVoigt_K0->Draw("same");
  fVoigt_K0->SetLineColor(kGreen);
  fbg_K0->Draw("same");
  fbg_K0->SetLineColor(kBlue);
  h1k0_m_all_cut->Draw();

  fVoigt_bg_L1115->Write();
  fVoigt_L1115->Write();
  fbg_L1115->Write();

  fVoigt_bg_K0->Write();
  fVoigt_K0->Write();
  fbg_K0->Write();

 //scale histograms 
  TIter next(MyFile->GetListOfKeys());
  TKey *key;
 
  
  while ((key = (TKey*)next()))
    {
      TClass *cl = gROOT->GetClass(key->GetClassName());
      if (!cl->InheritsFrom("TH1") || cl->InheritsFrom("TH2"))
	continue;
      TH1 *h = (TH1*)key->ReadObj();
      cout<<"histogram name:" <<h->GetName()<<endl;
      normalize(h);
      h->Write(0,TObject::kWriteDelete);
    }
  
  
  MyFile->Close();
}
