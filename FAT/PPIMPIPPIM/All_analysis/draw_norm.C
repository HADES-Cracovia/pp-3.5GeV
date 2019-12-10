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
  TFile *fileS1385 = new TFile("SB_sim_S1385pK0.root","READ");
  TFile *fileSDpp = new TFile("SB_sim_SDppK0.root","READ");
  TFile *fileLDpp = new TFile("SB_sim_LDppK0.root","READ");
  TFile *fileL1520= new TFile("SB_sim_L1520pippim.root","READ");
  TFile *fileExp= new TFile("SB_experiment.root","READ");
  
  TFile *output= new TFile("pictures.root","RECREATE");

  TH1F *hS1385_data = (TH1F*)fileS1385->Get("data");
  TH1F *hSDpp_data = (TH1F*)fileSDpp->Get("data");
  TH1F *hLDpp_data = (TH1F*)fileLDpp->Get("data");
  TH1F *hexperiment_data=(TH1F*)fileExp->Get("data");
  TH1F *hL1520_data=(TH1F*)fileL1520->Get("data");
  TH1F *hsum_data=(TH1F*)hS1385_data->Clone("hsum_data");
  hsum_data->Reset();
  
  TH1F *hS1385_background = (TH1F*)fileS1385->Get("background");
  TH1F *hSDpp_background = (TH1F*)fileSDpp->Get("background");
  TH1F *hLDpp_background = (TH1F*)fileLDpp->Get("background");
  TH1F *hexperiment_background=(TH1F*)fileExp->Get("background");
  TH1F *hL1520_background=(TH1F*)fileL1520->Get("background");
  
  TH1F *hsum_background=(TH1F*)hS1385_background->Clone("hsum_background");
  TH1F *hclean_background=(TH1F*)hS1385_background->Clone("hclean_background");
  TH1F *hclean_experiment=(TH1F*)hexperiment_background->Clone("hclean_experiment");
  TH1F *hclean_L1520=(TH1F*)hL1520_background->Clone("hclean_L1520");
  hsum_background->Reset();
  hclean_background->Reset();
  hclean_experiment->Reset();
  hclean_L1520->Reset();
  
  //scale according to CS
  double nsim=40*TMath::Power(10,6);//number of simulated events
  double scale=3.13*TMath::Power(10,8);
  double cs[3]={14.05/1000*scale/nsim,//S1385
		9.26/1000*scale/nsim,//SDpp
		29.45/1000*scale/nsim};//LDpp

  double cs_sig;
  // cs in \mu barns, have to me re-calculated to mb!!

  hS1385_background->Scale(cs[0]);
  hSDpp_background->Scale(cs[1]);
  hLDpp_background->Scale(cs[2]);

  hS1385_data->Scale(cs[0]);
  hSDpp_data->Scale(cs[1]);
  hLDpp_data->Scale(cs[2]);

  hsum_background->Add(hS1385_background);
  hsum_background->Add(hSDpp_background);
  hsum_background->Add(hLDpp_background);

  hsum_data->Add(hS1385_data);
  hsum_data->Add(hSDpp_data);
  hsum_data->Add(hLDpp_data);

  hclean_background->Add(hsum_data,hsum_background,1,-1);
  hclean_experiment->Add(hexperiment_data,hexperiment_background,1,-1);
  hclean_L1520->Add(hL1520_data,hL1520_background,1,-1);

  cs_sig=1/(hclean_L1520->Integral())*20;
  hclean_L1520->Scale(cs_sig);    
  
  TCanvas *res=new TCanvas("res","res");
  res->Divide(2,2);
  res->cd(1);
  hS1385_data->Draw();
  hS1385_background->SetLineColor(kRed);
  hS1385_background->Draw("same");
  res->cd(2);
  hSDpp_data->Draw();
  hSDpp_background->SetLineColor(kRed);
  hSDpp_background->Draw("same");
  res->cd(3);
  hLDpp_data->Draw();
  hLDpp_background->SetLineColor(kRed);
  hLDpp_background->Draw("same");
  res->cd(4);
  hL1520_data->Draw();
  hL1520_background->SetLineColor(kRed);
  hL1520_background->Draw("same");
  
  
  TCanvas *cSum=new TCanvas("cSum","cSum");
  hexperiment_data->Draw();
  hexperiment_background->SetLineColor(kRed);
  hexperiment_background->Draw("same");
  hsum_data->Draw("same");
  hsum_background->SetLineColor(kRed);
  hsum_background->Draw("Same");

  int rebin=2;
  TCanvas *cClean=new TCanvas("cClean","cClean");
  hclean_experiment->Draw();
  hclean_experiment->Rebin(rebin);
  hclean_background->SetLineColor(kRed);
  hclean_background->Rebin(rebin);
  hclean_background->Draw("same");
  hclean_L1520->SetLineColor(kGreen);
  hclean_L1520->Rebin(rebin);
  hclean_L1520->Draw("same");
}
