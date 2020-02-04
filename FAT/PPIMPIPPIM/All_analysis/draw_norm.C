void setHistogramStyleData(TH1* hist)
{
  hist->SetLineWidth(2);
  
  hist->SetMarkerColor(hist->GetLineColor());
  hist->SetMarkerSize(2);
  hist->SetMarkerStyle(8);
}

void setHistogramStyleSimul(TH1* hist)
{
  hist->SetLineWidth(2);
}

double hist_error(TH1* hist, double x1=2, double x2=1)
{
  int nbin_min;
  int nbin_max;
  double err_sum=0;

  cout<<endl<<"***calculating sum of errors in range***"<<endl;
  cout<<"min value= "<<x1<<endl;
  cout<<"max value= "<<x2<<endl;
  
  if(x1>x2)
    {
      nbin_min=1;
      nbin_max=hist->GetNbinsX();
    }
  else
    {
      nbin_min=hist->FindBin(x1);
      nbin_max=hist->FindBin(x2);
    }
  cout<<"min bin= "<<nbin_min<<endl;
  cout<<"max bin= "<<nbin_max<<endl;
  
  for(int i=nbin_min;i<=nbin_max;i++)
    {
      cout<<"bin number: "<<i<<" bin error: "<<hist->GetBinError(i)<<endl;
      err_sum=err_sum+hist->GetBinError(i)*hist->GetBinError(i); 
    }

  cout<<"***end of hist_error function***"<<endl<<endl;
  return TMath::Sqrt(err_sum);
}

void scale_error(TH1* hist, double err)
{
  cout<<endl<<"***scaling the histogram errors according to one relative error***"<<endl;
  cout<<"hist name"<<hist->GetName()<<endl;
  cout<<"scaling error"<<err<<endl;
  
  int nbin_min=1;
  int nbin_max=hist->GetNbinsX();
  for(int i=nbin_min;i<=nbin_max;i++)
    {
      double cont=hist->GetBinContent(i);
      double error=hist->GetBinContent(i)*err;
      cout<<"bin number: "<<i<<" bin contetnt: "<<cont<<endl;
      cout<<"                       bin error:  "<<error<<endl;
      hist->SetBinError(i,error); 
    }

  cout<<"***end of scale_error function***"<<endl<<endl;
  
}


void normalize(TH1* hist)
{
  for (Int_t j=1; j<hist->GetNbinsX()+1; ++j)
    {
      double scale=1.0/(3.13 * TMath::Power(10,8)) *1000; /*to get micro barn*/
      double binErr=hist->GetBinError(j);
      hist->SetBinContent(j, hist->GetBinContent(j) / hist->GetBinWidth(j) *scale);
      //hist->SetBinError( j, TMath::Sqrt( hist->GetBinContent(j) ) );
      hist->SetBinError( j,  binErr/ hist->GetBinWidth(j) * scale );
      hist->GetYaxis()->SetTitle("#frac{dN}{dE} [#frac{#mu b}{MeV}]");
    }

}  

int draw_norm(void)
{

  TF1 *voigt=new TF1("signal_voit","[0]*TMath::Voigt(x-[1],[2],[3],4)",1380,1750);
  TFile *fileS1385 = new TFile("SB_sim_S1385pK0.root","READ");
  TFile *fileSDpp = new TFile("SB_sim_SDppK0.root","READ");
  TFile *fileLDpp = new TFile("SB_sim_LDppK0.root","READ");
  TFile *fileL1520= new TFile("SB_sim_L1520pippim.root","READ");
  TFile *fileExp= new TFile("SB_experiment.root","READ");
  
  //TFile *output= new TFile("pictures.root","RECREATE");
  TH1F *hS1385_data = (TH1F*)fileS1385->Get("data");
  hS1385_data->SetName("hS1385_data");
  hS1385_data->Sumw2(kFALSE);
  TH1F *hSDpp_data = (TH1F*)fileSDpp->Get("data");
  hSDpp_data->SetName("hSDpp_data");
  hSDpp_data->Sumw2(kFALSE);;
  TH1F *hLDpp_data = (TH1F*)fileLDpp->Get("data");
  hLDpp_data->SetName("hLDpp_data");
  hLDpp_data->Sumw2(kFALSE);
  TH1F *hexperiment_data=(TH1F*)fileExp->Get("data");
  hexperiment_data->SetName("hexperiment_data");
  //hexperiment_data->Sumw2();
  TH1F *hL1520_data=(TH1F*)fileL1520->Get("data");
  hL1520_data->SetName("hL1520_data");
  hL1520_data->Sumw2(kFALSE);
  TH1F *hsum_data=(TH1F*)hS1385_data->Clone("hsum_data");
  hsum_data->Reset();
  
  TH1F *hS1385_background = (TH1F*)fileS1385->Get("background");
  hS1385_background->SetName("hS1385_background");
  hS1385_background->Sumw2(kFALSE);
  TH1F *hSDpp_background = (TH1F*)fileSDpp->Get("background");
  hSDpp_background->SetName("hSDpp_background");
  hSDpp_background->Sumw2(kFALSE);
  TH1F *hLDpp_background = (TH1F*)fileLDpp->Get("background");
  hLDpp_background->SetName("hLDpp_background");
  hLDpp_background->Sumw2(kFALSE);
  TH1F *hexperiment_background=(TH1F*)fileExp->Get("background");
  hexperiment_background->SetName("hexperiment_background");
  //hexperiment_background->Sumw2();
  TH1F *hL1520_background=(TH1F*)fileL1520->Get("background");
  hL1520_background->SetName("hL1520_background");
  hL1520_background->Sumw2(kFALSE);

  TH1F *hexperiment_SB_spectrum=(TH1F*)fileExp->Get("oryginal_spectrum");
  
  
  TH1F *hsum_background=(TH1F*)hS1385_background->Clone("hsum_background");
  TH1F *hclean_background=(TH1F*)hS1385_background->Clone("hclean_background");
  TH1F *hclean_experiment=(TH1F*)hexperiment_background->Clone("hclean_experiment");
  TH1F *hclean_L1520=(TH1F*)hL1520_background->Clone("hclean_L1520");
  TH1F *hclean_sum=(TH1F*)hL1520_background->Clone("hclean_sum");
  TH1F *hclean_L1520_ren=(TH1F*)hL1520_background->Clone("hclean_L1520_ren");
  TH1F *hclean_sum_ren=(TH1F*)hL1520_background->Clone("hclean_sum_ren");
  TH1F *hpure_signal=(TH1F*)hL1520_background->Clone("hpure_signal");
  hsum_background->Reset();
  hclean_background->Reset();
  hclean_experiment->Reset();
  hclean_L1520->Reset();
  hclean_sum->Reset();
  hclean_L1520_ren->Reset();
  hclean_sum_ren->Reset();
  hpure_signal->Reset();

  
  //scale according to CS
  //double nsim=40*TMath::Power(10,6);//number of simulated events
  double nsim=120*TMath::Power(10,6);
  double scale=3.13*TMath::Power(10,8);
  double downscale=3;//trigger downscale for simulated events
  double cs[4]=
    {14.05/1000*scale/(nsim*downscale),//S1385
     9.26/1000*scale/(nsim*downscale),//SDpp
     29.45/1000*scale/(nsim*downscale),//LDpp
     5.6/1000*scale/(100*100000*downscale)//L(1520)pK+->Lpi+pi-pK+
    };
  double err[4]=
    {2.25/14.05,//S1385
     1.47/9.26,//SDpp
     2.55/29.45,//LDpp
     0/5.6//L(1520)pK+->Lpi+pi-pK+
    };
  double cs_sig;
  // cs in \mu barns, have to me re-calculated to mb!!

  hS1385_background->Scale(cs[0]);
  hSDpp_background->Scale(cs[1]);
  hLDpp_background->Scale(cs[2]);
  hL1520_background->Scale(cs[3]);

  /*  hS1385_background->Sumw2();
  hSDpp_background->Sumw2();
  hLDpp_background->Sumw2();
  hL1520_background->Sumw2();
  */
  scale_error(hS1385_background,err[0]);
  scale_error(hSDpp_background,err[1]);
  scale_error(hLDpp_background,err[2]);
  scale_error(hL1520_background,err[3]);
  
  hS1385_data->Scale(cs[0]);
  hSDpp_data->Scale(cs[1]);
  hLDpp_data->Scale(cs[2]);
  hL1520_data->Scale(cs[3]);
  /*
  hS1385_data->Sumw2();
  hSDpp_data->Sumw2();
  hLDpp_data->Sumw2();
  hL1520_data->Sumw2();
  */
  scale_error(hS1385_data,err[0]);
  scale_error(hSDpp_data,err[1]);
  scale_error(hLDpp_data,err[2]);
  scale_error(hL1520_data,err[3]);
    
  hsum_background->Add(hS1385_background);
  hsum_background->Add(hSDpp_background);
  hsum_background->Add(hLDpp_background);

  hsum_data->Add(hS1385_data);
  hsum_data->Add(hSDpp_data);
  hsum_data->Add(hLDpp_data);

  hclean_background->Add(hsum_data,hsum_background,1,-1);
  hclean_experiment->Add(hexperiment_data,hexperiment_background,1,-1);
  hclean_L1520->Add(hL1520_data,hL1520_background,1,-1);
  hclean_sum->Add(hclean_L1520,hclean_background,1,1);
  //cs_sig=1/(hclean_L1520->Integral())*20;
  //hclean_L1520->Scale(cs_sig);    

  //scale signal to difference between signal and background
  double int_min=1410;
  double int_max=1600;
  double err_sum;
  
  double sig_int=hclean_L1520->Integral(hclean_L1520->FindBin(int_min),hclean_L1520->FindBin(int_max));
  double backgroud_int=hclean_background->Integral(hclean_background->FindBin(int_min),hclean_background->FindBin(int_max));
  double experiment_int=hclean_experiment->Integral(hclean_experiment->FindBin(int_min),hclean_experiment->FindBin(int_max));

  hclean_L1520_ren->Add(hclean_L1520,1);
  hclean_L1520_ren->Scale((experiment_int-backgroud_int)/sig_int);
  hclean_sum_ren->Add(hclean_L1520_ren,1);
  hclean_sum_ren->Add(hclean_background,1);

  
  TCanvas *cRes=new TCanvas("cRes","cRes");
  cRes->Divide(2,2);
  cRes->cd(1);
  hS1385_data->Draw();
  hS1385_background->SetLineColor(kRed);
  hS1385_background->Draw("same");
  cRes->cd(2);
  hSDpp_data->Draw();
  hSDpp_background->SetLineColor(kRed);
  hSDpp_background->Draw("same");
  cRes->cd(3);
  hLDpp_data->Draw();
  hLDpp_background->SetLineColor(kRed);
  hLDpp_background->Draw("same");
  cRes->cd(4);
  hL1520_data->Draw();
  hL1520_background->SetLineColor(kRed);
  hL1520_background->Draw("same");
  
  
  TCanvas *cSum=new TCanvas("cSum","cSum");
  int rebin2=2; //wrong error propagation for simul events
  hexperiment_data->Rebin(rebin2);
  hexperiment_data->Draw();
  setHistogramStyleData(hexperiment_data);

  hexperiment_background->Rebin(rebin2);
  hexperiment_background->SetLineColor(kRed);
  hexperiment_background->Draw("same");
  setHistogramStyleData(hexperiment_background);

  hsum_data->Rebin(rebin2);
  hsum_data->Draw("same");
  setHistogramStyleSimul(hsum_data);
  
  hsum_background->Rebin(rebin2);
  hsum_background->SetLineColor(kRed);
  setHistogramStyleSimul(hsum_background);
  hsum_background->Draw("Same");

  int rebin=4;
  TCanvas *cClean=new TCanvas("cClean","cClean");
  hclean_experiment->Draw();
  hclean_experiment->Rebin(rebin);
  
  hclean_background->SetLineColor(kRed);
  hclean_background->SetFillColor(kRed);
  hclean_background->Rebin(rebin);
  hclean_background->Draw("samee2B");
  setHistogramStyleSimul(hclean_background);
  
  hclean_L1520->SetLineColor(kGreen);
  hclean_L1520->Rebin(rebin);
  hclean_L1520->Draw("same");
  setHistogramStyleSimul(hclean_L1520);
  
  hclean_sum->Rebin(rebin);
  hclean_sum->SetLineColor(kMagenta);
  hclean_sum->SetFillColor(kMagenta);
  hclean_sum->Draw("samee2B");
  setHistogramStyleSimul(hclean_sum);

  TCanvas *cClean_ren=new TCanvas("cClean_ren","cClean_ren");
  cClean_ren->Divide(2);
  cClean_ren->cd(1);
  hclean_experiment->Draw();
  hclean_experiment->GetXaxis()->SetRangeUser(1360,1780);
  setHistogramStyleData(hclean_experiment);
  
  //hclean_experiment->Rebin(rebin);
  //hclean_background->SetLineColor(kRed);
  //hclean_background->Rebin(rebin);
  hclean_background->Draw("samee2B");
  setHistogramStyleSimul(hclean_background);
  hclean_L1520_ren->SetLineColor(kGreen);
  hclean_L1520_ren->Rebin(rebin);
  hclean_L1520_ren->Draw("same");
  setHistogramStyleSimul(hclean_L1520_ren);

  hclean_sum_ren->Rebin(rebin);
  hclean_sum_ren->SetLineColor(kMagenta);
  hclean_sum_ren->SetFillColor(kMagenta);
  hclean_sum_ren->Draw("samee2B");
  setHistogramStyleSimul(hclean_sum_ren);

  cClean_ren->cd(2);
  hpure_signal->Rebin(rebin);
  hpure_signal->GetXaxis()->SetRangeUser(1360,1780);
  hpure_signal->Draw();
    
  //fit Voigt to data
  hpure_signal->Add(hclean_experiment,hclean_background,1,-1);
  setHistogramStyleData(hpure_signal);
  voigt->SetParameter(0,2412);
  voigt->SetParameter(1,1500);
  voigt->SetParameter(2,5);
  voigt->SetParameter(3,50);
  hpure_signal->Fit(voigt,"RL");
  hpure_signal->Fit(voigt,"RL");

  
  TCanvas *cSB=new TCanvas("cSB","Spectrum for side-band");
  hexperiment_SB_spectrum->Draw();

  err_sum=hist_error(hpure_signal,int_min,int_max);
  
  cout<<"Integral for pK0L(1520) (CS from Laura paper):"<<endl;
  cout<<hclean_L1520->Integral()<<endl;
  cout<<"Integral for inclusive L(1520) production:"<<endl;
  cout<<hclean_L1520_ren->Integral()<<endl;
  cout<<"C-S for pp->pK0L(1520):"<<endl;
  cout<<"5.6 \mu b:"<<endl;
  cout<<"inclusive L(1520) production C-S:"<<endl;
  cout<<5.6*(experiment_int-backgroud_int)/sig_int<<endl;
  cout<<"a scaling factor"<<endl;
  cout<<(experiment_int-backgroud_int)/sig_int<<endl;
  cout<<"****************error estimation****************"<<endl;
  cout<<"error sum= "<<err_sum<<endl;
  cout<<endl<<endl;
  
  //save all
  TFile* output=new TFile("final_output.root","recreate");

  hS1385_data->Write();
  hSDpp_data->Write();
  hLDpp_data->Write();
  hexperiment_data->Write();
  hL1520_data->Write();
  hsum_data->Write();
  
  hS1385_background->Write();
  hSDpp_background->Write();
  hLDpp_background->Write();
  hexperiment_background->Write();
  hL1520_background->Write();
  
  hsum_background->Write();
  hclean_background->Write();
  hclean_experiment->Write();
  hclean_L1520->Write();
  hclean_sum->Write();
  hexperiment_SB_spectrum->Write();

  hclean_L1520_ren->Write();
  hclean_sum_ren->Write();
  hpure_signal->Write();
  
  cClean_ren->Write();
  cRes->Write();
  cClean->Write();
  cSum->Write();
  cSB->Write();
}

