#define ppimpippim_cxx

#include "ppimpippim.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <iostream>
#include <TKey.h>
#include <TMath.h>
#include <TObject.h>
#include <TF1.h>
#include <TGraph2D.h>

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
  TFile *MyFile = new TFile("result.root","RECREATE");
  
  
  if ( MyFile->IsOpen() )
    printf("File opened successfully\n");
  else
    printf("Cann't open a file\n");
  //MyFile->cd();
  
  const int ndist=50;
  const int nangle=30;

  int nmass=500;
  int minmass=1000;
  int maxmass=2000;
  
  double maxdist=100;
  double mindist=-50;
  double maxangle=90;

  char hist_name[100];
  char sig_name[100];
  char bg_name[100];
  char sig_bg_name[100];
  char hist_title[1000];

  TH1F* m_L1520[ndist][nangle];
  TH1F* signal_sim[ndist][nangle];
  TH1F* background_sim[ndist][nangle];

  TF1* sig[ndist][nangle];
  TF1* bg[ndist][nangle];
  TF1* sig_bg[ndist][nangle];

  double dist;
  double int_min=1504;
  double int_max=1535;
  double int_sig[ndist][nangle];
  double int_sig_sim[ndist][nangle];
  double int_bg[ndist][nangle];
  double int_bg_sim[ndist][nangle];
  double sig_to_bg[ndist][nangle];
  double sig_to_bg_sim[ndist][nangle];
  double signif[ndist][nangle];
  double signif_sim[ndist][nangle];

  TGraph2D* gsignificance=new TGraph2D(ndist*nangle);
  TGraph2D* gsig_to_bg=new TGraph2D(ndist*nangle);
  gsignificance->SetName("gSignificance");
  gsignificance->SetTitle("significance");
  gsig_to_bg->SetName("gSig_To_Bg");
  gsig_to_bg->SetTitle("S/B");

  TGraph2D* gsignificance_sim=new TGraph2D(ndist*nangle);
  TGraph2D* gsig_to_bg_sim=new TGraph2D(ndist*nangle);
  gsignificance_sim->SetName("gsignificance_sim");
  gsignificance_sim->SetTitle("significance");
  gsig_to_bg_sim->SetName("gsig_to_bg_sim");
  gsig_to_bg_sim->SetTitle("S/B");

  //create histograms
  for(int i=0; i<ndist;i++)
    for(int j=0; j<nangle;j++)
      {
	double dist_cut=mindist+(double)i/ndist*(maxdist-mindist);
	double angle_cut=(double)j/nangle*maxangle;
	
	sprintf(hist_name,"L_1520_dist_%.2f_angle_%.2f",dist_cut,angle_cut);
	sprintf(hist_title,"#Lambda (1520) dist_{min}=%.2f angle_{max}=%.2f",dist_cut,angle_cut);
	sprintf(sig_name,"sig_dist_%.2f_angle_%.2f",dist_cut,angle_cut);
	sprintf(bg_name,"bg_dist_%.2f_angle_%.2f",dist_cut,angle_cut);
	sprintf(sig_bg_name,"sig_bg_dist_%.2f_angle_%.2f",dist_cut,angle_cut);
	m_L1520[i][j]=new TH1F(hist_name,hist_title,nmass,minmass,maxmass);
	sig[i][j]=new TF1(sig_name,"[0]*TMath::Voigt(x-[1],[2],[3])",1500,1540);
	bg[i][j]=new TF1(bg_name,"pol4(0)",1400,1640);
	sig_bg[i][j]=new TF1(sig_bg_name,"[0]*TMath::Voigt(x-[1],[2],[3])+pol4(4)",1480,1560);

	sprintf(hist_name,"L_1520_dist_%.2f_angle_%.2f_sim_signal",dist_cut,angle_cut);
	sprintf(hist_title,"#Lambda (1520) dist_{min}=%.2f angle_{max}=%.2f signal-id from sim",dist_cut,angle_cut);
	signal_sim[i][j]=new TH1F(hist_name,hist_title,nmass,minmass,maxmass);

	sprintf(hist_name,"L_1520_dist_%.2f_angle_%.2f_sim_background",dist_cut,angle_cut);
	sprintf(hist_title,"#Lambda (1520) dist_{min}=%.2f angle_{max}=%.2f background-id from sim",dist_cut,angle_cut);
	background_sim[i][j]=new TH1F(hist_name,hist_title,nmass,minmass,maxmass);
      }


  
  if (fChain == 0) return;

  Long64_t nentries = fChain->GetEntries();
  Long64_t nbytes = 0, nb = 0;

  //Fill histograms
  for (Long64_t jentry=0; jentry<nentries;jentry++)
    {
      Long64_t ientry = LoadTree(jentry);

      //dist=dist_ver_to_ver;
      dist=ver_p_pim_z-ver_pip_pim_z;
      
      if (ientry < 0)
	break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      if(jentry%5000==0)
	std::cout<<"progress "<<(double)jentry/nentries*100<<" %"<<endl;
      //simulation id part
      if(
	 mlp_output<0.58
	 ||miss_mass_kp<1430
	 ||isBest_new<1
	 )
	continue;

      

      if(p_sim_parentid==18
	 &&pim_sim_parentid==18
	 &&pip_sim_parentid==0
	 &&p_sim_id==14
	 )
	for(int i=0; i<ndist;i++)
	  for(int j=0; j<nangle;j++)
	    {
	      double dist_cut=mindist+(double)i/ndist*(maxdist-mindist);
	      double angle_cut=(double)j/nangle*maxangle;
	      if(dist>dist_cut && (oa_lambda<angle_cut || oa_lambda>180-angle_cut))
		signal_sim[i][j]->Fill(m_inv_p_pim_pip_pim);
	    }
      else
	for(int i=0; i<ndist;i++)
	  for(int j=0; j<nangle;j++)
	    {
	      double dist_cut=mindist+(double)i/ndist*(maxdist-mindist);
	      double angle_cut=(double)j/nangle*maxangle;
	      if(dist>dist_cut && (oa_lambda<angle_cut || oa_lambda>180-angle_cut))
		background_sim[i][j]->Fill(m_inv_p_pim_pip_pim);
	    }

      // analysis-like part
      // if (Cut(ientry) < 0) continue;
      for(int i=0; i<ndist;i++)
	for(int j=0; j<nangle;j++)
	  {
	    double dist_cut=mindist+(double)i/ndist*(maxdist-mindist);
	    double angle_cut=(double)j/nangle*maxangle;
	    if(dist>dist_cut && (oa_lambda<angle_cut || oa_lambda>180-angle_cut))
	      m_L1520[i][j]->Fill(m_inv_p_pim_pip_pim);
	  }
    }

  //fit all histograms
  int counter=0;
  for(int i=0; i<ndist;i++)
    for(int j=0; j<nangle;j++)
      {
	if(i==0 && j==0)
	  counter=0;
	else
	  counter++;

	double dist_cut=mindist+(double)i/ndist*(maxdist-mindist);
	double angle_cut=(double)j/nangle*maxangle;
	//calculate integrals from sim-id part
	int_bg_sim[i][j]=background_sim[i][j]->Integral(background_sim[i][j]->FindBin(int_min),background_sim[i][j]->FindBin(int_max));
	int_sig_sim[i][j]=signal_sim[i][j]->Integral(signal_sim[i][j]->FindBin(int_min),signal_sim[i][j]->FindBin(int_max));
	if(int_bg_sim[i][j]!=0)
	  {
	    sig_to_bg_sim[i][j]=int_sig_sim[i][j]/int_bg_sim[i][j];
	    signif_sim[i][j]=int_sig_sim[i][j]/TMath::Sqrt(int_bg_sim[i][j]+int_sig_sim[i][j]);
	  }
	else
	  {
	    sig_to_bg_sim[i][j]=-1;
	    signif_sim[i][j]=-1;
	  }
	//fitting part
	m_L1520[i][j]->Sumw2();
	
	sig[i][j]->SetParameters(m_L1520[i][j]->GetMaximum()*10,1520,5,16);
	m_L1520[i][j]->Fit(sig[i][j],"R");

	for(int k=0;k<=3;k++)
	  sig_bg[i][j]->SetParameter(k,sig[i][j]->GetParameter(k));
	sig_bg[i][j]->SetLineColor(kGreen);
	m_L1520[i][j]->Fit(sig_bg[i][j],"R");
	sig_bg[i][j]->SetRange(1470,1580);
	m_L1520[i][j]->Fit(sig_bg[i][j],"R");
	sig_bg[i][j]->SetRange(1450,1600);
	m_L1520[i][j]->Fit(sig_bg[i][j],"R");
	sig_bg[i][j]->SetRange(1450,1640);
	m_L1520[i][j]->Fit(sig_bg[i][j],"R");

	for(int k=4;k<sig_bg[i][j]->GetNpar();k++)
	  bg[i][j]->SetParameter(k-4,sig_bg[i][j]->GetParameter(k));
	bg[i][j]->SetLineColor(kRed);
	bg[i][j]->Draw("same");

	for(int k=0;k<4;k++)
	  sig[i][j]->SetParameter(k,sig_bg[i][j]->GetParameter(k));
	sig[i][j]->SetRange(1460,1600);
	sig[i][j]->SetLineColor(kBlue);

	//calculate integrals
	double ch2=sig_bg[i][j]->GetChisquare()/sig_bg[i][j]->GetNDF();
	  if(ch2<2 && ch2>0.5)//check fit quality
	  {
	    int_bg[i][j]=bg[i][j]->Integral(int_min,int_max)/m_L1520[i][j]->GetBinWidth(1);
	    int_sig[i][j]=sig[i][j]->Integral(int_min,int_max)/m_L1520[i][j]->GetBinWidth(1);
	    sig_to_bg[i][j]=int_sig[i][j]/int_bg[i][j];
	    signif[i][j]=int_sig[i][j]/TMath::Sqrt(int_bg[i][j]+int_sig[i][j]);
	  }
	else
	  {
	    int_bg[i][j]=0;
	    int_sig[i][j]=0;
	    sig_to_bg[i][j]=0;
	    signif[i][j]=0;
	  }
	gsig_to_bg->SetPoint(counter,dist_cut,angle_cut,sig_to_bg[i][j]);
	gsignificance->SetPoint(counter,dist_cut,angle_cut,signif[i][j]);
	gsig_to_bg_sim->SetPoint(counter,dist_cut,angle_cut,sig_to_bg_sim[i][j]);
	gsignificance_sim->SetPoint(counter,dist_cut,angle_cut,signif_sim[i][j]);
      }
  //write all histograms

  gStyle->SetPalette(1);

  TCanvas* cSignif=new TCanvas("cSignif");
  gsignificance->Draw("surf1");
  gsignificance->Write();
  cSignif->Write();

  TCanvas* cSig_To_Bg=new TCanvas("cSig_To_Bg");
  gsig_to_bg->Draw("surf1");
  gsig_to_bg->Write();
  cSig_To_Bg->Write();

  TCanvas* cSignif_sim=new TCanvas("cSignif_sim");
  gsignificance_sim->Draw("surf1");
  gsignificance_sim->Write();
  cSignif_sim->Write();

  TCanvas* cSig_To_Bg_sim=new TCanvas("cSig_To_Bg_sim");
  gsig_to_bg_sim->Draw("surf1");
  gsig_to_bg_sim->Write();
  cSig_To_Bg_sim->Write();

 
  //TCanvas* cHistos=new TCanvas("cHistos");
  for(int i=0; i<ndist;i++)
    for(int j=0; j<nangle;j++)
      {
      m_L1520[i][j]->Write();
      sig[i][j]->Write();
      bg[i][j]->Write();
      sig_bg[i][j]->Write();
      signal_sim[i][j]->Write();
      background_sim[i][j]->Write();
      }
  
  MyFile->Close();
}
