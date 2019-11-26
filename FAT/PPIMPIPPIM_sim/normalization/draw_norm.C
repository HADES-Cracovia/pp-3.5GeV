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
  TFile *fileS1385 = new TFile("./../StatisticalTest/Rafal_sim_S1385pK0.root","READ");
  TFile *fileSDpp = new TFile("./../StatisticalTest/Rafal_sim_SDppK0.root","READ");
  TFile *fileLDpp = new TFile("./../StatisticalTest/Rafal_sim_LDppK0.root","READ");

  TFile *output= new TFile("output.root","RECREATE");

  TH1F *hS1385_data = (TH1F*)fileS1385->Get("data");
  TH1F *hSDpp_data = (TH1F*)fileSDpp->Get("data");
  TH1F *hLDpp_data = (TH1F*)fileLDpp->Get("data");
  TH1F *sum_data=(TH1F*)hS1385_data->Clone("sum_data");
  sum_data->Reset();
  
  TH1F *hS1385_background = (TH1F*)fileS1385->Get("background");
  TH1F *hSDpp_background = (TH1F*)fileSDpp->Get("background");
  TH1F *hLDpp_background = (TH1F*)fileLDpp->Get("background");
  TH1F *sum_background=(TH1F*)hS1385_background->Clone("sum_background");
  sum_background->Reset();
  //scale according to CS
  double nsim=40*TMath::Power(10,6);//number of simulated events
  double scale=3.13*TMath::Power(10,8)/(13.6*1000)/nsim;
  double cs[3]={14.05*scale,//S1385
		9.26*scale,//SDpp
		29.45*scale};//LDpp
  // cs in \mu barns!!

  hS1385_background->Scale(cs[0]);
  hSDpp_background->Scale(cs[1]);
  hLDpp_background->Scale(cs[2]);

  hS1385_data->Scale(cs[0]);
  hSDpp_data->Scale(cs[1]);
  hLDpp_data->Scale(cs[2]);

  sum_background->Add(hS1385_background);
  sum_background->Add(hSDpp_background);
  sum_background->Add(hLDpp_background);

  sum_data->Add(hS1385_data);
  sum_data->Add(hSDpp_data);
  sum_data->Add(hLDpp_data);
  
  TCanvas *res=new TCanvas("res","res");
  res->Divide(3);
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
  
  TCanvas *cSum=new TCanvas("cSum","cSum");
  sum_data->Draw();
  sum_background->SetLineColor(kRed);
  sum_background->Draw("Same");
}
