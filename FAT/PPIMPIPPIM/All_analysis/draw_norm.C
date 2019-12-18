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
  double nsim=40*TMath::Power(10,6);//number of simulated events
  double scale=3.13*TMath::Power(10,8);
  double downscale=3;//trigger downscale for simulated events
  double cs[4]={14.05/1000*scale/(nsim*downscale),//S1385
		9.26/1000*scale/(nsim*downscale),//SDpp
		29.45/1000*scale/(nsim*downscale),//LDpp
		35.26*0.06/1000*scale/(100*100000*downscale)//L(1520)pK+->Lpi+pi-pK+
  };
  double cs_sig;
  // cs in \mu barns, have to me re-calculated to mb!!

  hS1385_background->Scale(cs[0]);
  hSDpp_background->Scale(cs[1]);
  hLDpp_background->Scale(cs[2]);
  hL1520_background->Scale(cs[3]);
  
  hS1385_data->Scale(cs[0]);
  hSDpp_data->Scale(cs[1]);
  hLDpp_data->Scale(cs[2]);
  hL1520_data->Scale(cs[3]);
  
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
  double int_max=1590;
  
  double sig_int=hclean_L1520->Integral(hclean_L1520->FindBin(int_min),hclean_L1520->FindBin(int_max));
  double backgroud_int=hclean_background->Integral(hclean_background->FindBin(int_min),hclean_background->FindBin(int_max));
  double experiment_int=hclean_experiment->Integral(hclean_experiment->FindBin(int_min),hclean_experiment->FindBin(int_max));

  hclean_L1520_ren->Add(hclean_L1520,1);
  hclean_L1520_ren->Scale((experiment_int-backgroud_int)/sig_int);
  hclean_sum_ren->Add(hclean_L1520_ren,1);
  hclean_sum_ren->Add(hclean_background,1);
    
  cout<<"Integral for pK0L(1520):"<<endl;
  cout<<hclean_L1520->Integral()<<endl;
  cout<<"Integral for inclusive L(1520) production:"<<endl;
  cout<<hclean_L1520_ren->Integral()<<endl;
  cout<<"C-S for pp->pK0L(1520):"<<endl;
  cout<<"35.26 \mu b:"<<endl;
  cout<<"inclusive L(1520) production C-S:"<<endl;
  cout<<35.26*(experiment_int-backgroud_int)/sig_int<<endl;
  cout<<"a scaling factor"<<endl;
  cout<<(experiment_int-backgroud_int)/sig_int<<endl;

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
  hexperiment_data->Rebin(2);
  hexperiment_data->Draw();
  hexperiment_background->Rebin(2);
  hexperiment_background->SetLineColor(kRed);
  hexperiment_background->Draw("same");  
  hsum_data->Rebin(2);
  hsum_data->Draw("same");
  hsum_background->Rebin(2);
  hsum_background->SetLineColor(kRed);
  hsum_background->Draw("Same");

  int rebin=4;
  TCanvas *cClean=new TCanvas("cClean","cClean");
  hclean_experiment->Draw();
  hclean_experiment->Rebin(rebin);
  hclean_background->SetLineColor(kRed);
  hclean_background->Rebin(rebin);
  hclean_background->Draw("same");
  hclean_L1520->SetLineColor(kGreen);
  hclean_L1520->Rebin(rebin);
  hclean_L1520->Draw("same");
  hclean_sum->Rebin(rebin);
  hclean_sum->SetLineColor(kMagenta);
  hclean_sum->Draw("same");

  TCanvas *cClean_ren=new TCanvas("cClean_ren","cClean_ren");
  cClean_ren->Divide(2);
  cClean_ren->cd(1);
  hclean_experiment->Draw();
  //hclean_experiment->Rebin(rebin);
  //hclean_background->SetLineColor(kRed);
  //hclean_background->Rebin(rebin);
  hclean_background->Draw("same");
  hclean_L1520_ren->SetLineColor(kGreen);
  hclean_L1520_ren->Rebin(rebin);
  hclean_L1520_ren->Draw("same");
  hclean_L1520_ren->GetXaxis()->SetRange(hclean_L1520_ren->FindBin(1350),hclean_L1520_ren->FindBin(1850));
  hclean_sum_ren->Rebin(rebin);
  hclean_sum_ren->SetLineColor(kMagenta);
  hclean_sum_ren->Draw("same");
  cClean_ren->cd(2);
  hpure_signal->Rebin(rebin);
  hpure_signal->GetXaxis()->SetRange(hpure_signal->FindBin(1350),hpure_signal->FindBin(1850));
  hpure_signal->Draw();
  
  TCanvas *cSB=new TCanvas("cSB","Spectrum for side-band");
  hexperiment_SB_spectrum->Draw();


  //fit Voigt to data
  hpure_signal->Add(hclean_experiment,hclean_background,1,-1);
  voigt->SetParameter(0,2412);
  voigt->SetParameter(1,1500);
  voigt->SetParameter(2,5);
  voigt->SetParameter(3,50);
  hpure_signal->Fit(voigt,"RL");
  hpure_signal->Fit(voigt,"RL");

  
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
  
  cRes->Write();
  cClean->Write();
  cSum->Write();
  cSB->Write();
}

