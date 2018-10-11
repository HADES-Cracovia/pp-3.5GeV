
int fit_DZ()
{
  TFile *f = new TFile("pp_ppim_all.root");
  TCanvas* c1=new TCanvas("c1","c1");
  TCanvas* c2=new TCanvas("c2","c2");

  c1->Divide(5,5);
  c2->Divide(2);
  
  TH1F *hists[25][25];
  f->ls();
  char hname[40];
  char new_name[40];
  char f_bg[40];
  char f_sg[40];
  char f_al[40];
  char f_sc[40];
  char f_bc[40];

  TF1 *bg[25][25];
  TF1 *signal[25][25];
  TF1 *all[25][25];
  TF1 *signal_cp[25][25];
  TF1 *bg_cp[25][25];
  //Load histograms
  
  double d_cut[25];
  double z_cut[25];
  double d_signal[25][25];
  double d_bcg[25][25];
  double d_sig_to_bg[25][25];
  double d_signif[25][25];

  
  for(int n=1;n<=25;n++)
    {
      d_cut[n-1]=2*n;
      z_cut[n-1]=(n*2)-10;
      for(int z=1;z<=25;z++)
	{
	  sprintf(hname,"DZ_p_pim_mass_%d_%d",n*2,(z*2)-10);
	  sprintf(new_name,"mass_%d_%d",n*2,(z*2)-10);
	  hists[n-1][z-1]=new TH1F("new_name","new_name",2000,500,2500);
	  hists[n-1][z-1]=(TH1F*)f->Get(hname);
	}
    }
  //fit histograms
  for(int j=0;j<25;j++)
    for(int k=0; k<25; k++)
    {
      //c1->cd(j+1);
      hists[j][k]->GetXaxis()->SetRange(hists[j][k]->FindFixBin(1070),hists[j][k]->FindFixBin(1200));
      //hists[j][k]->Draw();
      sprintf(f_bg,"background_dist_%d_z_%d",(j+1)*2,(k+1)*2-10);
      sprintf(f_sg,"signal_dist_%d_z_%d",(j+1)*2,(k+1)*2-10);
      sprintf(f_al,"all_dist_%d_z_%d",(j+1)*2,(k+1)*2-10);
      sprintf(f_sc,"signal_copy_dist_%d_z_%d",(j+1)*2,(k+1)*2-10);
      sprintf(f_bc,"background_copy_dist_%d_z_%d",(j+1)*2,(k+1)*2-10);
      
      bg[j][k]=new TF1(f_bg,"pol4",1080,1180);
      bg[j][k]->SetLineColor(kBlue);

      signal[j][k]=new TF1(f_sg,"gaus(0)+pol1(4)",1110,1120);
      signal[j][k]->SetLineColor(kMagenta);

      all[j][k]=new TF1(f_al,"pol4(0)+gaus(5)",1080,1150);
      all[j][k]->SetLineColor(kGreen);

      signal_cp[j][k]=new TF1(f_sc,"gaus(0)",1110,1120);
      signal_cp[j][k]->SetLineColor(kRed);

      bg_cp[j][k]=new TF1(f_bc,"pol4",1080,1180);
      bg_cp[j][k]->SetLineColor(kBlue);
      
      hists[j][k]->Fit(bg[j][k],"R");
      double start_base=0.5*(hists[j][k]->GetBinContent(hists[j][k]->FindBin(1110))+hists[j][k]->GetBinContent(hists[j][k]->FindBin(1120)));
      double start_high=hists[j][k]->GetBinContent(hists[j][k]->FindBin(1115))-start_base;
      cout<<"starting parameters "<<start_base<<" "<<start_high<<endl;
      signal[j][k]->SetParameters(start_high,1115,5,start_base);
      hists[j][k]->Fit(signal[j],"R");
      
      all[j][k]->SetParameter(0,bg[j][k]->GetParameter(0));
      all[j][k]->SetParameter(1,bg[j][k]->GetParameter(1));
      all[j][k]->SetParameter(2,bg[j][k]->GetParameter(2));
      all[j][k]->SetParameter(3,bg[j][k]->GetParameter(3));
      all[j][k]->SetParameter(4,bg[j][k]->GetParameter(4));
      all[j][k]->SetParameter(5,signal[j][k]->GetParameter(0));
      all[j][k]->SetParameter(6,signal[j][k]->GetParameter(1));
      all[j][k]->SetParameter(7,signal[j][k]->GetParameter(2));

      hists[j][k]->Fit(all[j][k],"R");
      
      signal_cp[j][k]->SetParameter(0,all[j][k]->GetParameter(5));
      signal_cp[j][k]->SetParameter(1,all[j][k]->GetParameter(6));
      signal_cp[j][k]->SetParameter(2,all[j][k]->GetParameter(7));

      bg_cp[j][k]->SetParameter(0,all[j][k]->GetParameter(0));
      bg_cp[j][k]->SetParameter(1,all[j][k]->GetParameter(1));
      bg_cp[j][k]->SetParameter(2,all[j][k]->GetParameter(2));
      bg_cp[j][k]->SetParameter(3,all[j][k]->GetParameter(3));
      bg_cp[j][k]->SetParameter(4,all[j][k]->GetParameter(4));
      
      //signal_cp[j][k]->Draw("same");
      //signal[j][k]->Draw("same");
      //bg_cp[j][k]->Draw("same");

      d_signal[j][k]=signal_cp[j][k]->Integral(1110,1120);
      d_bcg[j][k]=bg_cp[j][k]->Integral(1110,1120);
      d_sig_to_bg[j][k]=d_signal[j][k]/d_bcg[j][k];
      d_signif[j][k]=d_signal[j][k]/TMath::Sqrt(d_signal[j][k]+d_bcg[j][k]);
    }

  
  TGraph2D *sig_to_bacg=new TGraph2D();
  TGraph2D *signifi=new TGraph2D();
  int npoint=0;

  for(int l=0;l<25;l++)
    for(int j=0;j<25;j++)
      {
	sig_to_bacg->SetPoint(npoint,d_cut[l],z_cut[j],d_sig_to_bg[l][j]);
	signifi->SetPoint(npoint,d_cut[l],z_cut[j],d_signif[l][j]);
      }

  
  c2->cd(1);
  sig_to_bacg->Draw();
  c2->cd(2);
  signifi->Draw();
  
  return 0;
}
