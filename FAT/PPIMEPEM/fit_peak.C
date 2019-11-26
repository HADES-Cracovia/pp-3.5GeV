int fit_peak(void)
{
  TFile *file = new TFile("output_full.root");
  if(file->IsOpen())
    cout<<"File open"<<endl;
  else
    cout<<"file opening problem"<<endl;
  //TFile* file = TFile::Open("data.root", "RECREATE");  // my output file
  TH1F* h1 = (TH1F*)file->Get("miss_mass_RF_OK");
  int sig_min=950;
  int sig_max=1050;
  int bcg_min=1000;
  int bcg_max=1200;
  TF1* signal=new TF1("signal","gaus(0)",sig_min,sig_max);
  signal->SetLineColor(kGreen);
  signal->SetRange(sig_min,sig_max);
  signal->SetParameter(0,h1->GetBinContent(h1->GetBin((sig_max+sig_min)/2)));
  signal->SetParameter(1,(sig_max+sig_min)/2);
  signal->SetParLimits(1,sig_min,sig_max);
  signal->SetParameter(2,(sig_max-sig_min)/10);
  signal->SetParLimits(2,0,(sig_max-sig_min));
  //signal->SetParameter(3,(h1->GetBinContent(h1->GetBin(sig_min))+h1->GetBinContent(h1->GetBin(sig_max)))/2);
  //signal->SetParLimits(2,0,(sig_max-sig_min)/2);
  
  
  TF1* bcg=new TF1("bcg","pol1(0)",bcg_min,bcg_max);
  bcg->SetLineColor(kViolet);
  bcg->SetRange(bcg_min,bcg_max);
  //bcg->SetParameter(0,h1->GetMaximumStored());
  //bcg->SetParLimits(0,h1->GetMaximumStored()/3,h1->GetMaximumStored()*2);
  //bcg->SetParameter(1,h1->GetBinCenter(h1->GetMaximumBin()));
  //bcg->SetParLimits(1,bcg_min,bcg_max);
  
  h1->Fit(signal);
  h1->Fit(bcg);

  TF1* sig_bc=new TF1("sig_bc","gaus(0)+pol1(3)",0.9*sig_min,1.1*sig_max);
 
  sig_bc->SetRange(0.9*sig_min,1.1*sig_max);
  sig_bc->SetParameter(0,signal->GetParameter(0));
  sig_bc->SetParameter(1,signal->GetParameter(1));
  sig_bc->SetParameter(2,signal->GetParameter(2));
  sig_bc->SetParameter(3,bcg->GetParameter(0));
  sig_bc->SetParameter(4,bcg->GetParameter(1));

  //h1->Fit(sig_bc);


  signal->SetParameter(0,sig_bc->GetParameter(0));
  signal->SetParameter(1,sig_bc->GetParameter(1));
  signal->SetParameter(2,sig_bc->GetParameter(2));
  bcg->SetParameter(0,sig_bc->GetParameter(3));
  bcg->SetParameter(1,sig_bc->GetParameter(4));

  
 
  
  h1->Draw();
  signal->Draw("same");
  bcg->Draw("same");
  sig_bc->Draw("same");
  
}
