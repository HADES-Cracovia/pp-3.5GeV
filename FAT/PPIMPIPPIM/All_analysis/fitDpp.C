int fitDpp(void)
{
  TFile* f_dpp=new TFile("h_dpp2.root","READ");
  TFile* f_exp=new TFile("h_exp.root","READ");

  TH1F* h_dpp=(TH1F*)f_dpp->Get("f1");
  TH1F* h_exp=(TH1F*)f_exp->Get("f2");
  TH1F* h_dpp_scaled=(TH1F*)h_dpp->Clone("h_dpp_scaled");
  TH1F* h_bg=(TH1F*)h_exp->Clone("h_bg");
  
  //scaling a simulated signal to estimate background
  double scale=2;
  h_dpp_scaled->Scale(scale);
  h_bg->Add(h_dpp_scaled,-1);
  
  TCanvas* cData=new TCanvas("cData");
  h_dpp->Draw();
  h_exp->Draw("same");
  h_dpp->SetLineColor(kGreen);

  TCanvas* cScale=new TCanvas("cScale");
  h_exp->Draw();
  h_dpp_scaled->Draw("same");
  h_dpp_scaled->SetLineColor(kGreen);
  h_bg->Draw("same");
  h_bg->SetLineColor(kRed);

  return 0;
}
