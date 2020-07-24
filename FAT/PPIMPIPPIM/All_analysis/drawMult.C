void norm(TH1* h)
{
  h->Scale(1/h->Integral());
  h->SetLineWidth(3);
}

int drawMult(void)
{
  TFile* f_poz=new TFile("/lustre/hades/user/knowakow/PP/PAT_1/FILES/ppimpippim/h_pos_mult.root","READ");
  TFile* f_neg=new TFile("/lustre/hades/user/knowakow/PP/PAT_1/FILES/ppimpippim/h_neg_mult.root","READ");

  TH1F* h_mult_poz=(TH1F*)f_poz->Get("hposmult");
  TH1F* h_mult_neg=(TH1F*)f_neg->Get("hnegmult");
  TH1F* h_mult_sum=(TH1F*)h_mult_neg->Clone("h_mult_sum");
  h_mult_sum->Add(h_mult_poz);

  norm(h_mult_poz);
  norm(h_mult_neg);
  norm(h_mult_sum);
  
  //draw results
  TCanvas* cRes=new TCanvas("cRes","particles multiplicity");
  h_mult_poz->Draw("same");
  h_mult_poz->SetLineColor(kRed);
  h_mult_neg->Draw("same");
  h_mult_neg->SetLineColor(kBlue);
  h_mult_sum->Draw("same");
  h_mult_sum->SetLineColor(kGreen-1);
  //save results
  TFile* out =new TFile("mult_out.root","recreate");
  h_mult_poz->Write();
  h_mult_neg->Write();
  h_mult_sum->Write();

  
  
  return 0;
}
