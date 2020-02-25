
#define createHistos_cxx
#include "createHistos.h"
#include <TCutG.h>
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TMath.h>
#include <TGraphErrors.h>
#include <TH1F.h>
#include <TFile.h>
#include <iostream>
#include <TF1.h>
#include <TLine.h>

using namespace std;

void createHistos::Loop(char* output)
{
  //   In a ROOT session, you can do:
  //      root> .L createHistos.C
  //      root> createHistos t
  //      root> t.GetEntry(12); // Fill t data members with entry number 12
  //      root> t.Show();       // Show values of entry 12
  //      root> t.Show(16);     // Read and show values of entry 16
  //      root> t.Loop();       // Loop on all entries
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
  gStyle->SetOptStat(0);

  if (fChain == 0) return;
  const int bin=200;
  const int xmin=1000;
  const int xmax=2000;
  const int nsignal=20;
  double sidebandmin=10;
  double sidebandmax=27;
  TLine* line1=new TLine(1116-sidebandmax,0,1116-sidebandmax,120);
  TLine* line2=new TLine(1116-sidebandmin,0,1116-sidebandmin,120);
  TLine* line3=new TLine(1116+sidebandmin,0,1116+sidebandmin,120);
  TLine* line4=new TLine(1116+sidebandmax,0,1116+sidebandmax,120);
  
  int step;
  TH1F* signal=new TH1F("signal","signal simulated from gaus",bin,xmin,xmax);
  TH1F* background=new TH1F("background","background from side-band;M^{inv}_{p #pi- e^{+} e^{-}}[MeV]",bin,xmin,xmax);
  TH1F* data=new TH1F("data","data from experiment;M^{inv}_{p #pi- e^{+} e^{-}}[MeV]",bin,xmin,xmax);
  TH1F* orginal_spectrum=new TH1F("orginal_spectrum","orginal spectrum for side-band;M^{inv}_{p #pi-}[MeV]",bin*2,xmin,xmax);
  TH1F* hMEpEm_signal=new TH1F("hMEpEm_signal","M^{inv}_{#pi^{+} #pi^{-}} from #Lambda(1520)",120,0,600);
  TH1F* hMEpEm_background=new TH1F("hMEpEm_background","M^{inv}_{#pi^{+} #pi^{-}} from #Lambda(1520)",120,0,600);
  
  TGraphErrors* resi=new TGraphErrors(bin);
  TF1* background_fit=new TF1("background_fit","pol2(0)",1000,1200);
  TF1* K0_fit=new TF1("K0_fit","[0]*TMath::Voigt(x-[1],[2],[3])+pol2(4)",483,506);
  TF1* K0_signal=new TF1("K0_signal","[0]*TMath::Voigt(x-[1],[2],[3])",400,600);
  TF1* L1116_fit=new TF1("L1116_fit","[0]*TMath::Voigt(x-[1],[2],[3])+pol2(4)",1100,1125);
  TF1* L1116_signal=new TF1("L1116_signal","[0]*TMath::Voigt(x-[1],[2],[3])",1000,1200);
  
  TH1F* missing_mass_K0_L=new TH1F("missing_mass_K0_L","missing mass for #Lambda K^{0} candidates",1000,600,1600);
  TH2F* dedx_lambda=new TH2F("dedx_lambda","de/dx for #Lambda events",250,0,2000,250,0,18);
  TH2F* miss_m_vs_ep_em=new TH2F("miss_m_vs_ep_em","M^{miss} vs. M_{e^{+} e^{-}}",50,1340,1650,50,200,450);
  background->Sumw2();
  data->Sumw2();
  orginal_spectrum->Sumw2();

  int dM=1;
  int Lmin=1000;
  int Lmax=1500;
  int LdM=(Lmax-Lmin)/dM;
  int Kmin=250;
  int Kmax=950;
  int KdM=(Kmax-Kmin)/dM;
  //Histograms for all stages of analysis
  TH1F* hMPPim_start=new TH1F("hMPPim_start","M^{inv}_{p #pi^{-}} after identification cuts; M^{inv}_{p #pi^{-}} [MeV];N",LdM,Lmin,Lmax);
  TH1F* hMEpEm_start=new TH1F("hMEpEm_start","M^{inv}_{e^{+} e^{-}} after identification cuts; M^{inv}_{e^{+} e^{-}} [MeV];N",KdM,Kmin,Kmax);
  TH2F* miss_m_vs_ep_p_start=new TH2F("miss_m_vs_ep_p_start","M^{miss}_{p #pi^{-} e^{+} e^{-}} vs. M_{p #pi^{+}};M^{miss}_{p #pi^{-} e^{+} e^{-}}[MeV];M^{inv}_{#pi+ p}[MeV];N",90,500,1400,50,1100,1600);
  TH2F* p_pim_vs_ep_em_start=new TH2F("p_pim_vs_ep_em_start","M_{p #pi-} vs. M_{e^{+} e^{-}};M^{inv}_{p #pi^{-}}[MeV];M^{inv}_{e^{+} e^{-}}[MeV];N",200,1050,1450,200,250,700);
  TH1F* hMPPim_TMVA=new TH1F("hMPPim_TMVA","M^{inv}_{p #pi^{-}} after MLP; M^{inv}_{p #pi^{-}} [MeV];N",LdM,Lmin,Lmax);
  TH1F* hMEpEm_TMVA=new TH1F("hMEpEm_TMVA","M^{inv}_{e^{+} e^{-}} after MLP; M^{inv}_{e^{+} e^{-}} [MeV];N",KdM,Kmin,Kmax);
  TH2F* miss_m_vs_ep_p_TMVA=new TH2F("miss_m_vs_ep_p_TMVA","M^{miss}_{p #pi^{-} e^{+} e^{-}} vs. M_{#pi^{+} p};M^{miss}_{p #pi^{-} e^{+} e^{-}}[MeV];M^{inv}_{#pi^{+} p}[MeV];N",90,500,1400,50,1100,1600);
  TH2F* p_pim_vs_ep_em_TMVA=new TH2F("p_pim_vs_ep_em_TMVA","M_{p #pi-} vs. M_{e^{+} e^{-}};M^{inv}_{p #pi^{-}}[MeV];M^{inv}_{e^{+} e^{-}}[MeV];N",200,1050,1400,200,250,700);
  TH1F* hMPPim_TMVA_K0mass=new TH1F("hMPPim_TMVA_K0mass","M^{inv}_{p #pi^{-}} after MLP and a gate for K^{0}; M^{inv}_{p #pi^{-}} [MeV];N",LdM,Lmin,Lmax);
  TH1F* hMEpEm_TMVA_Lmass=new TH1F("hMEpEm_TMVA_Lmass","M^{inv}_{e^{+} e^{-}} after MLP and a gate for #Lambda; M^{inv}_{e^{+} e^{-}} [MeV];N",KdM,Kmin,Kmax);
  TH1F* hMPPim_TMVAMass=new TH1F("hMPPim_TMVAMass","M^{inv}_{p #pi^{-}} after MLP and a #Delta^{++} mass cut; M^{inv}_{p #pi^{-}} [MeV];N",LdM,Lmin,Lmax);
  TH1F* hMEpEm_TMVAMass=new TH1F("hMEpEm_TMVAMass","M^{inv}_{e^{+} e^{-}} after MLP and a #Delta^{++} mass cut; M^{inv}_{e^{+} e^{-}} [MeV];N",KdM,Kmin,Kmax);
  
  hMPPim_TMVA_K0mass->Sumw2();
  hMEpEm_TMVA_Lmass->Sumw2();
  
  
  TFile *cutFile=new TFile("/lustre/hades/user/knowakow/PP/FAT/PPIMPIPPIM_sim/TMVAeval_DD/cut_miss_mass_vs_pip_pim.root","READ");
  //TFile *cutFile=new TFile("/lustre/hades/user/knowakow/PP/FAT/PPIMPIPPIM_sim/TMVAeval_DD/cut_miss_pip_pim_tight.root","READ");
  TCutG *graph_cut=0;
  cutFile->GetObject("CUTG",graph_cut);
  cutFile->Close();

  double mlp_cut=0.57;
  TFile *MyFile = new TFile(output,"recreate");
 
  Long64_t nentries = fChain->GetEntries();
  Long64_t nbytes = 0, nb = 0;
  step =(int)nentries/15;

  for (Long64_t jentry=0; jentry<nentries;jentry++)
    {
      Long64_t ientry = LoadTree(jentry);
      if(jentry%step==0)
	cout<<"Progresss: "<<(double)jentry/nentries<<endl;
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;

      if(isBest==1)
	{
	  hMPPim_start->Fill(m_inv_p_pim);
	  hMEpEm_start->Fill(m_inv_ep_em);
	  //miss_m_vs_ep_p_start->Fill(miss_mass_kp,m_inv_p_ep);
	  p_pim_vs_ep_em_start->Fill(m_inv_p_pim,m_inv_ep_em);

	  if(mlp_output>mlp_cut)
	    {
	      hMPPim_TMVA->Fill(m_inv_p_pim);
	      hMEpEm_TMVA->Fill(m_inv_ep_em);
	      //miss_m_vs_ep_p_TMVA->Fill(miss_mass_kp,m_inv_p_ep);
	      p_pim_vs_ep_em_TMVA->Fill(m_inv_p_pim,m_inv_ep_em);
	      if(m_inv_p_pim<1120 && m_inv_p_pim>1110 && miss_mass_kp>1077)
		hMEpEm_TMVA_Lmass->Fill(m_inv_ep_em);
	      if(m_inv_ep_em<500 && m_inv_ep_em>480 && miss_mass_kp>1077)
		hMPPim_TMVA_K0mass->Fill(m_inv_p_pim);
	      if(graph_cut->IsInside(miss_mass_kp,m_inv_ep_em))
		{
		  hMPPim_TMVAMass->Fill(m_inv_p_pim);
		  hMEpEm_TMVAMass->Fill(m_inv_ep_em);
		}
	    }
	}
      
      
      if(
	 m_inv_ep_em<500
	 && m_inv_ep_em>480
	 && m_inv_p_pim<1126
	 && m_inv_p_pim>1106
	 && mlp_output<mlp_cut
	 && miss_mass_kp>1077
	 )//K0 and L(1116)
	{
	  dedx_lambda->Fill(ep_p,ep_dedx);
	  missing_mass_K0_L->Fill(miss_mass_kp);
	}


      //all events for final pictures
      if(isBest!=1
	 ||mlp_output<mlp_cut
	 //||miss_mass_kp<1432 //replaced by graphical cut
	 //||m_inv_ep_em>410 //replaced by graphical cut
	 ||dist_ver_to_ver<5
	 ||(oa_lambda>20)
	 ||!(graph_cut->IsInside(miss_mass_kp,m_inv_ep_em))
	 //||p_theta>20 //to clean up proton sample
	 //||dist_ep_em>5
	 //||dist_ep_em>150
	 //||ver_p_pim_z<-5
	 //||dist_ver_to_ver<14
	 )
	continue;
      orginal_spectrum->Fill(m_inv_p_pim);
      
      if(m_inv_p_pim<1116+sidebandmin && m_inv_p_pim>1116-sidebandmin)
	{
	  data->Fill(m_inv_p_pim_ep_em);
	  miss_m_vs_ep_em->Fill(miss_mass_kp,m_inv_ep_em);	  
	  if(m_inv_p_pim_ep_em>1440 && m_inv_p_pim_ep_em<1600)
	    hMEpEm_signal->Fill(m_inv_ep_em);
	}

      if(m_inv_p_pim<1116.-sidebandmin && m_inv_p_pim>1116.-sidebandmax)
	{
	  background->Fill(m_inv_p_pim_ep_em);
	  if(m_inv_p_pim_ep_em>1440 && m_inv_p_pim_ep_em<1600)
	    hMEpEm_background->Fill(m_inv_ep_em);
	}
      if(m_inv_p_pim>1116.+sidebandmin && m_inv_p_pim<1116.+sidebandmax)
	{
	  background->Fill(m_inv_p_pim_ep_em);
	  if(m_inv_p_pim_ep_em>1440 && m_inv_p_pim_ep_em<1600)
	    hMEpEm_background->Fill(m_inv_ep_em);
	
	}
    }

  //normalize background to signal
  TCanvas* cFit1116=new TCanvas("cFit1116");
  cFit1116->cd();
     
  TF1* fVoigt_bg= new TF1("fVoigt_bg","[0]*TMath::Voigt(x-[1],[2],[3])+pol5(4)",1090.00,1156.67);
  TF1* fVoigt= new TF1("fVoigt","[0]*TMath::Voigt(x-[1],[2],[3])",1090.00,1156.67);
  TF1* fbg= new TF1("fbg","pol5(0)",1090.00,1156.67);

  fVoigt_bg->SetParameters(527.3,1114,2.73,1,-9266,1.7,0.0061,5.51379e-6,1.23803e-9,-5.64175e-12);
  fVoigt_bg->SetParLimits(3,0,2);
  fVoigt_bg->SetParLimits(1,1112,1117);
  fVoigt_bg->SetRange(1106,1121);
  orginal_spectrum->Fit(fVoigt_bg,"R");
  fVoigt_bg->SetRange(1100,1126);
  orginal_spectrum->Fit(fVoigt_bg,"R");
  fVoigt_bg->SetRange(1092,1137);
  orginal_spectrum->Fit(fVoigt_bg,"R");
  fVoigt_bg->SetRange(1088,1141); 
  orginal_spectrum->Fit(fVoigt_bg,"R");
  fVoigt_bg->SetRange(1082,1141); 
  orginal_spectrum->Fit(fVoigt_bg,"R");


  orginal_spectrum->Draw();
  fbg->SetParameters(fVoigt_bg->GetParameter(4),fVoigt_bg->GetParameter(5),fVoigt_bg->GetParameter(6),fVoigt_bg->GetParameter(7),fVoigt_bg->GetParameter(8),fVoigt_bg->GetParameter(9));
  fVoigt->SetParameters(fVoigt_bg->GetParameter(0),fVoigt_bg->GetParameter(1),fVoigt_bg->GetParameter(2),fVoigt_bg->GetParameter(3));
  fVoigt->Draw("same");
  fVoigt->SetLineColor(kGreen);
  fbg->Draw("same");
  fbg->SetLineColor(kBlue);
   
  double intS=fVoigt->Integral(1116-sidebandmin,1116+sidebandmin)/orginal_spectrum->GetBinWidth(1);
  double intB=fbg->Integral(1116-sidebandmin,1116+sidebandmin)/orginal_spectrum->GetBinWidth(1);
  double intsideband=(fbg->Integral(1116-sidebandmax,1116-sidebandmin)+fbg->Integral(1116+sidebandmin,1116+sidebandmax))/orginal_spectrum->GetBinWidth(1);
  double intAll=fVoigt_bg->Integral(1116-sidebandmin,1116+sidebandmin)/orginal_spectrum->GetBinWidth(1);
  
  cout<<"signal integral: "<<intS<<endl<<"beckground integral: "<<intB<<endl<<"sideband integral: "<<intsideband<<endl;
  cout<<"all in signal range: "<<intAll<<endl;

  background->Scale(intB/intsideband);
  hMEpEm_background->Scale(intB/intsideband);
  //Fill random signal
  //TF1* L1520Spectral=new TF1("L1520Spectral","100*exp(-0.5*((x-1520)/16)**2)",xmin,xmax);
  TF1* L1520Spectral=new TF1("L1520Spectral","TMath::BreitWigner(x,1519.5,15.6)",xmin,xmax);
  signal->FillRandom("L1520Spectral",10000);

  signal->Scale((double)nsignal/10000);
  //Calculate residuals
  for(int i=0;i<bin;i++)
    {
      resi->SetPoint(i,background->GetBinCenter(i),data->GetBinContent(i)-background->GetBinContent(i));
      resi->SetPointError(i,background->GetBinWidth(i),TMath::Sqrt(TMath::Power(data->GetBinError(i),2)+TMath::Power(background->GetBinError(i),2)));
    }

  //Fit K0 i L1116
  double r_f=2; //r_f rebin fit histograms
  L1116_fit->SetParameters(1845*r_f*r_f*r_f,1115,0.059,7.5,-24741*r_f,44*r_f,-0.019*r_f);
  hMPPim_TMVA_K0mass->Rebin(r_f);
  hMPPim_TMVA_K0mass->Fit(L1116_fit,"R0");
  L1116_fit->SetRange(1088,1147);
  hMPPim_TMVA_K0mass->Fit(L1116_fit,"R0");
  L1116_signal->SetParameters(L1116_fit->GetParameter(0),L1116_fit->GetParameter(1),L1116_fit->GetParameter(2),L1116_fit->GetParameter(3));
  L1116_signal->SetLineColor(kGreen-2);
  //L1116_signal->Draw("same");

  K0_fit->SetParameters(1528*r_f*r_f,492,0.03,11,1207*r_f,-3.3*r_f,0.002*r_f);
  hMEpEm_TMVA_Lmass->Rebin(r_f);
  hMEpEm_TMVA_Lmass->Fit(K0_fit,"R0");
  K0_fit->SetRange(476,510);
  hMEpEm_TMVA_Lmass->Fit(K0_fit,"R0");
  K0_fit->SetRange(466,523);
  hMEpEm_TMVA_Lmass->Fit(K0_fit,"R0");
  K0_fit->SetRange(450,536);
  hMEpEm_TMVA_Lmass->Fit(K0_fit,"R0");
  K0_signal->SetParameters(K0_fit->GetParameter(0),K0_fit->GetParameter(1),K0_fit->GetParameter(2),K0_fit->GetParameter(3));
  K0_signal->SetLineColor(kGreen-2);
  //K0_signal->Draw("same");

  
  dedx_lambda->Write();
  cFit1116->Write();
  fVoigt->Write();
  fbg->Write();
  fVoigt_bg->Write();
  resi->Write();
  signal->Write();
  background->Write();
  data->Write();
  hMEpEm_signal->Write();
  hMEpEm_background->Write();
  orginal_spectrum->Write();
  missing_mass_K0_L->Write();
  miss_m_vs_ep_em->Write();
  graph_cut->Write();

  styleTH1(hMPPim_start);
  styleTH1(hMEpEm_start);
  styleTH1(hMPPim_TMVA);
  styleTH1(hMEpEm_TMVA);
  styleTH1(hMPPim_TMVA_K0mass);
  styleTH1(hMEpEm_TMVA_Lmass);
  styleTH1(hMPPim_TMVAMass);
  styleTH1(hMEpEm_TMVAMass);
 
  
  hMPPim_start->Write();
  hMEpEm_start->Write();
  miss_m_vs_ep_p_start->Write();
  p_pim_vs_ep_em_start->Write();
  hMPPim_TMVA->Write();
  hMEpEm_TMVA->Write();
  miss_m_vs_ep_p_TMVA->Write();
  p_pim_vs_ep_em_TMVA->Write();
  hMPPim_TMVA_K0mass->Write();
  hMEpEm_TMVA_Lmass->Write();
  hMPPim_TMVAMass->Write();
  hMEpEm_TMVAMass->Write(); 

  K0_fit->Write();
  K0_signal->Write();
  L1116_fit->Write();
  L1116_signal->Write();
   
  line1->Write("line1");
  line2->Write("line2");
  line3->Write("line3");
  line4->Write("line4");

  line1->Delete();
  line2->Delete();
  line3->Delete();
  line4->Delete();
    
  hMPPim_start->Delete();
  hMEpEm_start->Delete();
  miss_m_vs_ep_p_start->Delete();
  p_pim_vs_ep_em_start->Delete();
  hMPPim_TMVA->Delete();
  hMEpEm_TMVA->Delete();
  miss_m_vs_ep_p_TMVA->Delete();
  p_pim_vs_ep_em_TMVA->Delete();
  hMPPim_TMVA_K0mass->Delete();
  hMEpEm_TMVA_Lmass->Delete();
  hMPPim_TMVAMass->Delete();
  hMEpEm_TMVAMass->Delete();

  K0_fit->Delete();
  K0_signal->Delete();
  L1116_fit->Delete();
  L1116_signal->Delete();
 
  dedx_lambda->Delete();
  //cFit1116->Delete();
  fVoigt->Delete();
  fbg->Delete();
  fVoigt_bg->Delete();
  resi->Delete();
  signal->Delete();
  background->Delete();
  data->Delete();
  orginal_spectrum->Delete();
  missing_mass_K0_L->Delete();
  miss_m_vs_ep_em->Delete();
  graph_cut->Delete();


  MyFile->Close();
}
