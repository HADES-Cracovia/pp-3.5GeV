#define createHistos_cxx
#include "createHistos.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

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
   const int bin=250;
   const int xmin=1000;
   const int xmax=2000;
   const int nsignal=20;
   double sideband=10;
   int step;
   TH1F* signal=new TH1F("signal","signal simulated from gaus",bin,xmin,xmax);
   TH1F* background=new TH1F("background","background from side-band",bin,xmin,xmax);
   TH1F* data=new TH1F("data","data from experiment",bin,xmin,xmax);
   TGraphErrors* resi=new TGraphErrors(bin);
   
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
	 ||m_inv_pip_pim>410
	 ||dist_pip_pim>15
	 )
	continue;

      if(m_inv_p_pim<1120 && m_inv_p_pim>1110)
	data->Fill(m_inv_p_pim_pip_pim);
      if((m_inv_p_pim<1110 && m_inv_p_pim>1110-sideband)
	 ||(m_inv_p_pim>1120 && m_inv_p_pim<1120+sideband))
	background->Fill(m_inv_p_pim_pip_pim);
   }

   //normalize background to signal
   double intS=data->Integral();
   double intB=background->Integral();
   background->Scale(intS/intB);

   //Fill random signal
   TF1* L1520Spectral=new TF1("L1520Spectral","100*exp(-0.5*((x-1520)/16)**2)",xmin,xmax);
   signal->FillRandom("L1520Spectral",10000);

   signal->Scale((double)nsignal/10000);
   //Calculate residuals
   for(int i=0;i<bin;i++)
     {
       resi->SetPoint(i,background->GetBinCenter(i),data->GetBinContent(i)-background->GetBinContent(i));
       resi->SetPointError(i,background->GetBinWidth(i),TMath::Sqrt(TMath::Power(data->GetBinError(i),2)+TMath::Power(background->GetBinError(i),2)));
     }

   resi->Write();
   signal->Write();
   background->Write();
   data->Write();
   MyFile->Close();
}
