int fitDpp(void)
{
  TFile* f_dpp=new TFile("SB_sim_LDppK0.root","READ");
  TFile* f_exp=new TFile("SB_experiment.root","READ");
  TFile* f_ideal=new TFile("/lustre/hades/user/knowakow/PP/FAT/PPIMPIPPIM_sim/Dpp_ideal.root","read");
  TFile* out=new TFile("Dpp_out.root","recreate");
  
  TH1F* h_dpp=(TH1F*)f_dpp->Get("miss_m_start");
  h_dpp->SetName("h_dpp");
  TH1F* h_exp=(TH1F*)f_exp->Get("miss_m_start");
  h_exp->SetName("h_exp");
  TH1F* h_dpp_ideal=(TH1F*)f_ideal->Get("h2");

  int rebin =2;
  h_dpp->Rebin(rebin);
  h_exp->Rebin(rebin);
  h_dpp_ideal->Rebin(rebin);
  
  //scaling dpp and ideal dpp to the same area as data
  h_dpp->Scale(h_exp->GetMaximum()/h_dpp->GetMaximum());
  h_dpp_ideal->Scale(h_exp->GetMaximum()/h_dpp_ideal->GetMaximum());

  TH1F* h_dpp_scaled=(TH1F*)h_dpp->Clone("h_dpp_scaled");
  TH1F* h_bg=(TH1F*)h_exp->Clone("h_bg");

  
  //scaling a simulated signal to estimate background
  double scale=0.83;
  h_dpp_scaled->Scale(scale);
  h_bg->Add(h_dpp_scaled,-1);
 
  TCanvas* cData=new TCanvas("cData");
  h_dpp->Draw();
  h_exp->Draw("same");
  h_dpp->SetLineColor(kGreen);
  h_dpp_ideal->Draw("same");
  h_dpp_ideal->SetLineColor(kMagenta);
  
  TCanvas* cScale=new TCanvas("cScale");
  h_exp->Draw();
  h_dpp_scaled->Draw("same");
  h_dpp_scaled->SetLineColor(kGreen);
  h_bg->Draw("same");
  h_bg->SetLineColor(kRed);

  cScale->Write();
  cData->Write();
  out->Write();

  return 0;
}
