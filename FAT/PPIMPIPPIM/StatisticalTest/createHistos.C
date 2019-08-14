#define createHistos_cxx
#include "createHistos.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TMath.h>
#include <TGraphErrors.h>
#include <TH1F.h>
#include <TFile.h>
#include <iostream>
#include <TF1.h>

using namespace std;

void createHistos::Loop()
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
   if (fChain == 0) return;
   const int bin=200;
   const int xmin=1000;
   const int xmax=2000;
   const int nsignal=20;
   double sidebandmin=16;
   double sidebandmax=26;
   int step;
   TH1F* signal=new TH1F("signal","signal simulated from gaus",bin,xmin,xmax);
   TH1F* background=new TH1F("background","background from side-band",bin,xmin,xmax);
   TH1F* data=new TH1F("data","data from experiment",bin,xmin,xmax);
   TH1F* oryginal_spectrum=new TH1F("oryginal_spectrum","oryginal spectrum for side-band",bin*6,xmin,xmax);
   TGraphErrors* resi=new TGraphErrors(bin);
   TF1* background_fit=new TF1("background_fit","pol2(0)",1000,1200);   
   TFile *MyFile = new TFile("Event.root","recreate");
 
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
      if(isBest_new!=1
	 ||mlp_output<0.58
	 ||miss_mass_kp<1350
	 //||m_inv_pip_pim>410
	 //||dist_pip_pim>150
	 //||ver_p_pim_z<-5
	 //||dist_ver_to_ver<14
	 )
	continue;

      oryginal_spectrum->Fill(m_inv_p_pim);
      
      if(m_inv_p_pim<1120 && m_inv_p_pim>1110)
	data->Fill(m_inv_p_pim_pip_pim);
      if((m_inv_p_pim<1116-sidebandmin && m_inv_p_pim>1110-sidebandmax)
	 ||(m_inv_p_pim>1116+sidebandmin && m_inv_p_pim<1116+sidebandmax))
	background->Fill(m_inv_p_pim_pip_pim);
   }

   //normalize background to signal
   TCanvas* cFit1116=new TCanvas("cFit1116");
   cFit1116->cd();

   oryginal_spectrum->Sumw2();
   
   TF1* fVoigt_bg= new TF1("fVoigt_bg","[0]*TMath::Voigt(x-[1],[2],[3])+pol5(4)",1090.00,1156.67);
   TF1* fVoigt= new TF1("fVoigt","[0]*TMath::Voigt(x-[1],[2],[3])",1090.00,1156.67);
   TF1* fbg= new TF1("fbg","pol5(4)",1090.00,1156.67);

   fVoigt_bg->SetParameters(13406.8,1114.78,0.1819,15.63,-53756.9,9.35043,0.0342083,3.07071e-5,7.00426e-9,-3.02549e-11);
   oryginal_spectrum->Fit(fVoigt_bg,"R");
   fVoigt_bg->SetRange(1080,1165);
   //oryginal_spectrum->Fit(fVoigt_bg,"R");
   oryginal_spectrum->Draw();
   fbg->SetParameters(fVoigt_bg->GetParameter(4),fVoigt_bg->GetParameter(5),fVoigt_bg->GetParameter(6),fVoigt_bg->GetParameter(7),fVoigt_bg->GetParameter(8),fVoigt_bg->GetParameter(9));
   fVoigt->SetParameters(fVoigt_bg->GetParameter(0),fVoigt_bg->GetParameter(1),fVoigt_bg->GetParameter(2),fVoigt_bg->GetParameter(3));
   fVoigt->Draw("same");
   fVoigt->SetLineColor(kGreen);
   fbg->Draw("same");
   fbg->SetLineColor(kBlue);
   
   double intS=fVoigt->Integral(1116-sidebandmin,1116+sidebandmin);
   double intB=fbg->Integral(1116-sidebandmin,1116+sidebandmin);
   double intsideband=fbg->Integral(1116-sidebandmax,1116-sidebandmin)+fbg->Integral(1116+sidebandmin,1116+sidebandmax);
   background->Scale(intsideband/intB);

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

   cFit1116->Write();
   fVoigt->Write();
   fbg->Write();
   fVoigt_bg->Write();
   resi->Write();
   signal->Write();
   background->Write();
   data->Write();
   oryginal_spectrum->Write();
   MyFile->Close();
}