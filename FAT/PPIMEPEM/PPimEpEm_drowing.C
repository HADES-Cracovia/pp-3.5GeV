#define PPimEpEm_drowing_cxx
#include "PPimEpEm_drowing.h"
//#include "cut.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

void PPimEpEm_drowing::Loop()
{
  //   In a ROOT session, you can do:
  //      Root > .L PPimEpEm_drowing.C
  //      Root > PPimEpEm_drowing t
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
#if FORSIMUL==1
  TFile* out2=new TFile("pictures_comparison_simulation_ep_em.root","recreate");
#else
  TFile* out2=new TFile("pictures_Lambda_ep_em.root","recreate");
#endif
  cout<<"output file"<<endl;
  out2->Print();
  
#if FORSIMUL==1  
  Long64_t nentries = fChain->GetEntries();
  double SBmin=4;
  double SBmax=30;
  double SBasymetry=1.;
  int i_epem_min=0;
  int i_epem_max=600;
  int i_epem_n=30;
  int i_ppim_min=1000;
  int i_ppim_max=2000;
  int i_ppim_n=400;
  int i_ppimepem_min=1100;
  int i_ppimepem_max=2000;
  int i_ppimepem_n=40;
#else
  Long64_t nentries = fChain->GetEntries();
  double SBmin=5;
  double SBmax=40;
  double SBasymetry=1;
  int i_epem_min=0;
  int i_epem_max=600;
  int i_epem_n=15;
  int i_ppim_min=1000;
  int i_ppim_max=2000;
  int i_ppim_n=300;
  int i_ppimepem_min=1100;
  int i_ppimepem_max=2000;
  int i_ppimepem_n=40;
#endif

  TH1F* hMinv_p_pim=new TH1F("hMinv_p_pim","Mass without any cut;M^{inv}_{p #pi^{-}}[MeV];counts",i_ppim_n,i_ppim_min,i_ppim_max);
  TH1F* hMinv_p_pim_cut=new TH1F("hMinv_p_pim_cut","Mass after distance and Z-coordinate cut;M^{inv}_{p #pi^{-}}[MeV];counts",i_ppim_n,i_ppim_min,i_ppim_max);
  TH1F* hMinv_p_pim_cut_scale=new TH1F("hMinv_p_pim_cut_scale","Mass after distance and Z-coordinate cut;M^{inv}_{p #pi^{-}}[MeV];#frac{d #sigma}{d m} #left[#frac{#mu b}{MeV}#right]",i_ppim_n,i_ppim_min,i_ppim_max);
  TH1F* hMinv_ep_em_cut=new TH1F("hMinv_ep_em_cut","Di-lepton mass for #Lambda^{0};M^{inv}_{e^{+}e^{-}}[MeV];counts",i_epem_n,i_epem_min,i_epem_max);
  TH1F* hMinv_ep_em_cut_scale=new TH1F("hMinv_ep_em_cut_scale","Di-lepton mass for #Lambda^{0};M^{inv}_{e^{+}e^{-}}[MeV];#frac{d #sigma}{d m} #left[#frac{#mu b}{MeV}#right]",i_epem_n,i_epem_min,i_epem_max);
  TH1F* hMinv_ep_em_SB=new TH1F("hMinv_ep_em_SB","Di-lepton mass for #Lambda^{0};M^{inv}_{e^{+}e^{-}}[MeV];counts",i_epem_n,i_epem_min,i_epem_max);
  TH1F* hMinv_ep_em_SB_scale=new TH1F("hMinv_ep_em_SB_scale","Di-lepton mass for #Lambda^{0};M^{inv}_{e^{+}e^{-}}[MeV];#frac{d #sigma}{d m} #left[#frac{#mu b}{MeV}#right]",i_epem_n,i_epem_min,i_epem_max);
  TH1F* hMinv_ep_em_clean_scale=new TH1F("hMinv_ep_em_clean_scale","Di-lepton mass for #Lambda^{0};M^{inv}_{e^{+}e^{-}}[MeV];#frac{d #sigma}{d m} #left[#frac{#mu b}{MeV}#right]",i_epem_n,i_epem_min,i_epem_max);
  TH1F* hMinv_ep_em_clean_scale_extrapolated=new TH1F("hMinv_ep_em_clean_scale_extrapolated","Di-lepton mass for #Lambda^{0};M^{inv}_{e^{+}e^{-}}[MeV];#frac{d #sigma}{d m} #left[#frac{#mu b}{MeV}#right]",i_epem_n,i_epem_min,i_epem_max);
  TH1F* hMinv_ep_em_clean=new TH1F("hMinv_ep_em_clean","Di-lepton mass for #Lambda^{0};M^{inv}_{e^{+}e^{-}}[MeV];counts",i_epem_n,i_epem_min,i_epem_max);
  TH1F* hMinv_p_pim_ep_em_cut=new TH1F("hMinv_p_pim_ep_em_cut","Di-lepton + #Lambda^{0};M^{inv}_{p #pi^{-} e^{+}e^{-}}[MeV];counts",i_ppimepem_n,i_ppimepem_min,i_ppimepem_max);
  TH1F* hMinv_p_pim_ep_em_cut_scale=new TH1F("hMinv_p_pim_ep_em_cut_scale","Di-lepton + #Lambda^{0};M^{inv}_{p #pi^{-} e^{+}e^{-}}[MeV];#frac{d #sigma}{d m} #left[#frac{#mu b}{MeV}#right]",i_ppimepem_n,i_ppimepem_min,i_ppimepem_max);
  TH1F* hMinv_p_pim_ep_em_clean_scale=new TH1F("hMinv_p_pim_ep_em_clean_scale","Di-lepton + #Lambda^{0};M^{inv}_{p #pi^{-} e^{+}e^{-}}[MeV];#frac{d #sigma}{d m} #left[#frac{#mu b}{MeV}#right]",i_ppimepem_n,i_ppimepem_min,i_ppimepem_max);
  TH1F* hMinv_p_pim_ep_em_SB=new TH1F("hMinv_p_pim_ep_em_SB","Di-lepton + #Lambda^{0};M^{inv}_{p #pi^{-} e^{+}e^{-}}[MeV];counts",i_ppimepem_n,i_ppimepem_min,i_ppimepem_max);
  TH1F* hMinv_p_pim_ep_em_SB_scale=new TH1F("hMinv_p_pim_ep_em_SB_scale","Di-lepton + #Lambda^{0};M^{inv}_{p #pi^{-} e^{+}e^{-}}[MeV];#frac{d #sigma}{d m} #left[#frac{#mu b}{MeV}#right]",i_ppimepem_n,i_ppimepem_min,i_ppimepem_max);
TH1F* hMinv_p_pim_ep_em_clean=new TH1F("hMinv_p_pim_ep_em_clean","Di-lepton + #Lambda^{0};M^{inv}_{p #pi^{-} e^{+}e^{-}}[MeV];counts",i_ppimepem_n,i_ppimepem_min,i_ppimepem_max);
  TF1* fsum=new TF1("fsum","gaus(0)+pol2(3)",1100,1130);
  TF1* fsig=new TF1("fsig","gaus(0)",1084,1276);
  TF1* fbg=new TF1("fbg","pol2(0)",1084,1276);
  fsum->SetParameters(23.96,1115.15,3.03,2342.22,-3.83,0.0016);
  fsum->SetNpx(300);
  fsig->SetNpx(300);
  fbg->SetNpx(300);
  //lepton identif
  TH1F* hEp_beta=new TH1F("hEp_beta","#beta for e^{+};#beta_{e^{+}}",500,0,1.5);
  TH1F* hEm_beta=new TH1F("hEm_beta","#beta for e^{-};#beta_{e^{+}}",500,0,1.5);
  TH1F* hP_beta=new TH1F("hP_beta","#beta for p;#beta_{p}",500,0,1.5);
  TH1F* hPim_beta=new TH1F("hPim_beta","#beta for #pi^{-};#beta_{#pi^{-}}",500,0,1.5);

  TH2F* hEp_DeDx_p=new TH2F("hEp_DeDx_p","#frac{de}{dx} for e^{+};p[MeV];#frac{de}{dx}",100,0,1000,100,0,30);
  TH2F* hEm_DeDx_p=new TH2F("hEm_DeDx_p","#frac{de}{dx} for e^{-};p[MeV];#frac{de}{dx}",100,0,1000,100,0,30);
  TH2F* hP_DeDx_p=new TH2F("hP_DeDx_p","#frac{de}{dx} for p;p[MeV];#frac{de}{dx}",100,0,2000,100,0,30);
  TH2F* hPim_DeDx_p=new TH2F("hPim_DeDx_p","#frac{de}{dx} for #pi^{-};p[MeV];#frac{de}{dx}",100,0,1500,100,0,30);

  TH2F* hEp_Beta_p=new TH2F("hEp_Beta_p","#beta for e^{+};p[MeV];#beta",200,0,1000,200,0,2);
  TH2F* hEm_Beta_p=new TH2F("hEm_Beta_p","#beta for e^{-};p[MeV];#beta",200,0,1000,200,0,2);
  TH2F* hP_Beta_p=new TH2F("hP_Beta_p","#beta for p;p[MeV];#beta",200,0,2000,200,0,2);
  TH2F* hPim_Beta_p=new TH2F("hPim_Beta_p","#beta for #pi^{-};p[MeV];#beta",200,0,1500,200,0,2);
   
  TH2F* hEp_dPhi_p=new TH2F("hEp_dTheta_p","#Delta #Theta vs. p for Ep;#Delta #Theta [#circle];p[MeV]",100,-40,40,100,0,1000);

  //double sim_adjust=1./0.312;
#if FORSIMUL==1
  //double sim_adjust=3.5;//sigma pp@3.5/sigma pp@4.5 for L(1520)
  double sim_adjust=1.87;//sigma pp@3.5/sigma pp@4.5 for L(1116)
  
#else
  double sim_adjust=1;
#endif
  double lum=3.13*TMath::Power(10,8);
  double scale=1/lum *1000*sim_adjust /0.64;//to get #mu b from mb *L1116 BR
   
  Long64_t nbytes = 0, nb = 0;
  for (Long64_t jentry=0; jentry<nentries;jentry++)
    {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;

      hEp_beta->Fill(ep_beta);
      hEm_beta->Fill(em_beta);
      hP_beta->Fill(p_beta);
      hPim_beta->Fill(pim_beta);

      hEp_DeDx_p->Fill(ep_p,ep_dedx);
      hEm_DeDx_p->Fill(em_p,em_dedx);
      hP_DeDx_p->Fill(p_p,p_dedx);
      hPim_DeDx_p->Fill(pim_p,pim_dedx);

      hEp_Beta_p->Fill(ep_p,ep_beta);
      hEm_Beta_p->Fill(em_p,em_beta);
      hP_Beta_p->Fill(p_p,p_beta);
      hPim_Beta_p->Fill(pim_p,pim_beta);

      if(isBest==1 && trigdec==1)
	hMinv_p_pim->Fill(m_inv_p_pim);
		       
      if (!Cut(ientry)) continue;

      hMinv_p_pim_cut->Fill(m_inv_p_pim);
      hMinv_p_pim_cut_scale->Fill(m_inv_p_pim);

      if((m_inv_p_pim < 1116+SBmin) && (m_inv_p_pim>1116-SBmin))
	{
	  hMinv_ep_em_cut->Fill(m_inv_ep_em);
	  hMinv_ep_em_cut_scale->Fill(m_inv_ep_em);
	  if(m_inv_ep_em>140)
	    {
	      hMinv_p_pim_ep_em_cut->Fill(m_inv_p_pim_ep_em);
	      hMinv_p_pim_ep_em_cut_scale->Fill(m_inv_p_pim_ep_em);
	    }
	}
      
      if(((m_inv_p_pim > 1116+SBmin) && (m_inv_p_pim < 1116+SBmax*SBasymetry))
	 ||((m_inv_p_pim < 1116-SBmin) && (m_inv_p_pim > 1116-SBmax)))
	{
	  hMinv_ep_em_SB->Fill(m_inv_ep_em);
	  hMinv_ep_em_SB_scale->Fill(m_inv_ep_em);
	  if(m_inv_ep_em>140)
	    {
	      hMinv_p_pim_ep_em_SB->Fill(m_inv_p_pim_ep_em);
	      hMinv_p_pim_ep_em_SB_scale->Fill(m_inv_p_pim_ep_em);
	    }
	}
    }

  //p pim spectrum
  hMinv_p_pim_cut->Sumw2();
  hMinv_p_pim_cut_scale->Sumw2();
  hMinv_p_pim_cut_scale->Scale(scale/hMinv_p_pim_cut_scale->GetBinWidth(3));
  hMinv_p_pim_cut->Fit(fsum,"R");
  hMinv_p_pim_cut->Fit(fsum,"R");
  fsum->SetRange(1080,1160);
  hMinv_p_pim_cut->Fit(fsum,"R");
  hMinv_p_pim_cut->Fit(fsum,"R");
  fsum->SetRange(1080,1200);
  hMinv_p_pim_cut->Fit(fsum,"R");
  fsig->SetParameters(fsum->GetParameter(0),fsum->GetParameter(1),fsum->GetParameter(2));
  fbg->SetParameters(fsum->GetParameter(3),fsum->GetParameter(4),fsum->GetParameter(5));

  
  //ep em specrtum   
  hMinv_ep_em_cut->Sumw2();
  hMinv_ep_em_cut_scale->Sumw2();
  hMinv_ep_em_cut_scale->Scale(scale/hMinv_ep_em_cut_scale->GetBinWidth(3));

  //p pim ep em spectrum
  hMinv_p_pim_ep_em_cut->Sumw2();
  hMinv_p_pim_ep_em_cut_scale->Sumw2();
  hMinv_p_pim_ep_em_cut_scale->Scale(scale/hMinv_p_pim_ep_em_cut_scale->GetBinWidth(3));
   
  //side-band spectrum
  double scaleSB=fbg->Integral(1116-SBmin,1116+SBmin)/(fbg->Integral(1116-SBmax,1116-SBmin)+fbg->Integral(1116+SBmin,1116+SBmax*SBasymetry));
  hMinv_p_pim_ep_em_SB->Sumw2();
  hMinv_p_pim_ep_em_SB->Scale(scaleSB);
  hMinv_p_pim_ep_em_SB_scale->Sumw2();
  hMinv_p_pim_ep_em_SB_scale->Scale(scale/hMinv_p_pim_ep_em_SB_scale->GetBinWidth(3)*scaleSB);

  hMinv_ep_em_SB->Sumw2();
  hMinv_ep_em_SB->Scale(scaleSB);
  hMinv_ep_em_SB_scale->Sumw2();
  hMinv_ep_em_SB_scale->Scale(scale/hMinv_ep_em_SB_scale->GetBinWidth(3)*scaleSB);
   
  TCanvas* cCounts=new TCanvas("cCounts");
  hMinv_p_pim_cut->GetXaxis()->SetRangeUser(1000,1400);
  hMinv_p_pim_cut->Draw();
  fsum->Draw("same");
  fsig->Draw("same");
  fsig->SetLineColor(kGreen);
  fbg->Draw("same");
  fbg->SetLineColor(kBlue);

  TLatex *printFormula3 = new TLatex();                  
  char text12[10000];      
  sprintf(text12, "I_{exp}=%.1f",fsig->Integral(1080,1160)/hMinv_p_pim_cut->GetBinWidth(3));       
  printFormula3->SetNDC();                
  printFormula3->SetTextFont(32);        
  printFormula3->SetTextColor(1);          
  printFormula3->SetTextSize(0.05);         
  printFormula3->SetTextAlign(13);                                            
  printFormula3->DrawLatex(0.4,0.8, text12);                                                                                     

  TLine* l1=new TLine(1116-SBmax,0,1116-SBmax,fsum(1116-SBmax));
  TLine* l2=new TLine(1116-SBmin,0,1116-SBmin,fsum(1116-SBmin));
  TLine* l3=new TLine(1116+SBmin,0,1116+SBmin,fsum(1116+SBmin));
  TLine* l4=new TLine(1116+SBmax*SBasymetry,0,1116+SBmax*SBasymetry,fsum(1116+SBmax));

  l1->Draw("same");
  l2->Draw("same");
  l3->Draw("same");
  l4->Draw("same");
   
  TCanvas* cCS=new TCanvas("cCS");
  hMinv_p_pim_cut_scale->Draw();

  TCanvas* cBeta=new TCanvas("cBeta");
  cBeta->Divide(2);
  cBeta->cd(1);
  hEp_beta->Draw();
  hEm_beta->SetLineColor(kRed);
  hEm_beta->Draw("same");
  cBeta->cd(2);
  hP_beta->Draw();
  hPim_beta->SetLineColor(kRed);
  hPim_beta->Draw("same");

  TCanvas* cDeDx=new TCanvas("cDeDx");
  cDeDx->Divide(2,2);
  cDeDx->cd(1);
  hEp_DeDx_p->Draw("colz");
  cDeDx->cd(2);
  hEm_DeDx_p->Draw("colz");
  cDeDx->cd(3);
  hP_DeDx_p->Draw("colz");
  cDeDx->cd(4);
  hPim_DeDx_p->Draw("colz");

  TCanvas* cBetaP=new TCanvas("cBetaP");
  cBetaP->Divide(2,2);
  cBetaP->cd(1);
  hEp_Beta_p->Draw("colz");
  cBetaP->cd(2);
  hEm_Beta_p->Draw("colz");
  cBetaP->cd(3);
  hP_Beta_p->Draw("colz");
  cBetaP->cd(4);
  hPim_Beta_p->Draw("colz");

  TCanvas* cEpEm=new TCanvas("cEpEm");
  //cEpEm->SetLogy();
  cEpEm->Divide(2);
  cEpEm->cd(1);
  gPad->SetLogy();
  hMinv_ep_em_cut->Draw();
  hMinv_ep_em_SB->Draw("same");
  hMinv_ep_em_SB->SetLineColor(kRed);
  cEpEm->cd(2);
  gPad->SetLogy();
  hMinv_ep_em_cut_scale->Draw();
  hMinv_ep_em_SB_scale->Draw("same");
  hMinv_ep_em_SB_scale->SetLineColor(kRed);

  TCanvas* cPPimEpEm=new TCanvas("cPPimEpEm");
  //cEpEm->SetLogy();
  cPPimEpEm->Divide(2);
  cPPimEpEm->cd(1);
  //gPad->SetLogy();
  hMinv_p_pim_ep_em_cut->Draw();
  hMinv_p_pim_ep_em_SB->Draw("same");
  hMinv_p_pim_ep_em_SB->SetLineColor(kRed);
  cPPimEpEm->cd(2);
  //gPad->SetLogy();
  hMinv_p_pim_ep_em_cut_scale->Draw();
  hMinv_p_pim_ep_em_SB_scale->Draw("same");
  hMinv_p_pim_ep_em_SB_scale->SetLineColor(kRed);

 TCanvas* cPPimEpEmClean=new TCanvas("cPpimEpEmClean");
  cPPimEpEmClean->Divide(2);
  cPPimEpEmClean->cd(1);
  hMinv_p_pim_ep_em_clean=(TH1F*)hMinv_p_pim_ep_em_cut->Clone("hMinv_p_pim_ep_em_clean");
  hMinv_p_pim_ep_em_clean->Add(hMinv_p_pim_ep_em_SB,-1);
  hMinv_p_pim_ep_em_clean->Draw();
  cPPimEpEmClean->cd(2);
  hMinv_p_pim_ep_em_clean_scale=(TH1F*)hMinv_p_pim_ep_em_cut_scale->Clone("hMinv_p_pim_ep_em_clean_scale");
  hMinv_p_pim_ep_em_clean_scale->Add(hMinv_p_pim_ep_em_SB_scale,-1);
  hMinv_p_pim_ep_em_clean_scale->Draw();
  TCanvas* cEpEmClean=new TCanvas("cEpEmClean");
  cEpEmClean->Divide(2);
  cEpEmClean->cd(1);
  hMinv_ep_em_clean=(TH1F*)hMinv_ep_em_cut->Clone("hMinv_ep_em_clean");
  hMinv_ep_em_clean->Add(hMinv_ep_em_SB,-1);
  hMinv_ep_em_clean->Draw();
  cEpEmClean->cd(2);
  hMinv_ep_em_clean_scale=(TH1F*)hMinv_ep_em_cut_scale->Clone("hMinv_ep_em_clean_scale");
  hMinv_ep_em_clean_scale->Add(hMinv_ep_em_SB_scale,-1);
  hMinv_ep_em_clean_scale->Draw();

  TCanvas* cExtrapolated=new TCanvas("cExtrapolated");
  hMinv_ep_em_clean_scale_extrapolated=(TH1F*)hMinv_ep_em_clean_scale->Clone("hMinv_ep_em_clean_scale_extrapolated");
  double scale_1=69/10;
  hMinv_ep_em_clean_scale_extrapolated->Scale(scale_1);
  hMinv_ep_em_clean_scale_extrapolated->Draw();
  hMinv_ep_em_clean_scale->Draw("same");
  
  cCS->Write();
  cCounts->Write();
  cBeta->Write();
  cDeDx->Write();
  cBetaP->Write();
  cEpEm->Write();
  cPPimEpEm->Write();
  cPPimEpEmClean->Write();
  cEpEmClean->Write();
  cExtrapolated->Write();
  
  out2->Write();
  //out2->Close();
  //delete out2;
}
