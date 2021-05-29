void norm(TH1* h)
{
  h->Scale(1/h->Integral());
  h->SetLineWidth(3);
  h->SetLineStyle(2);
  h->SetFillColor(h->GetLineColor());

}

int drawMult(void)
{
  TFile* f_poz=new TFile("/lustre/hades/user/knowakow/PP/PAT_1/FILES/ppimpippim/h_pos_mult.root","READ");
  TFile* f_neg=new TFile("/lustre/hades/user/knowakow/PP/PAT_1/FILES/ppimpippim/h_neg_mult.root","READ");
  TFile* f_sum=new TFile("/lustre/hades/user/knowakow/PP/PAT_1/FILES/ppimpippim/h_sum_mult.root","READ");
  
  TH1F* h_mult_poz=(TH1F*)f_poz->Get("hposmult");
  TH1F* h_mult_neg=(TH1F*)f_neg->Get("hnegmult");
  TH1F* h_mult_sum=(TH1F*)f_sum->Get("hmultsum");

  
  //draw results
  TCanvas* cRes=new TCanvas("cRes","particles multiplicity");
  h_mult_neg->Draw("same");
  h_mult_neg->SetLineColor(kBlue);
  h_mult_neg->SetFillStyle(3354);
  h_mult_poz->Draw("same");
  h_mult_poz->SetLineColor(kRed);
  h_mult_poz->SetFillStyle(3345);
  h_mult_sum->Draw("same");
  h_mult_sum->SetLineColor(kGreen-1);
  h_mult_sum->SetFillStyle(3015);
  
  TLegend* legend1=new TLegend(0.6,0.6,0.9,0.9);
  legend1->AddEntry(h_mult_sum,"all particles");
  legend1->AddEntry(h_mult_poz,"+1 charge");
  legend1->AddEntry(h_mult_neg,"-1 charge");
  legend1->Draw();

  norm(h_mult_poz);
  norm(h_mult_neg);
  norm(h_mult_sum);
  

  //save results
  TFile* out =new TFile("mult_out.root","recreate");
  cRes->Write();
  h_mult_poz->Write();
  h_mult_neg->Write();
  h_mult_sum->Write();

  out->Write();
  
  return 0;
}
