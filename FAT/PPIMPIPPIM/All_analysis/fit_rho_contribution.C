int fit_rho_contribution(void)
{
  TFile* fin=new TFile("final_output.root","READ");
  TH1F* rho_component=(TH1F*)fin->Get("hL1520_hMPipPim_rho");
  TH1F* bez_rho_component=(TH1F*)fin->Get("hclean_L1520_ren_PipPim");
  TH1F* bg_component=(TH1F*)fin->Get("hclean_background_PipPim");
  TH1F* data_component=(TH1F*)fin->Get("hclean_experiment_PipPim");
  TH1F* sum_component=(TH1F*)fin->Get("hclean_sum_ren_PipPim");
  
  const int steps=10000;
  double x[steps];
  double y[steps];
  double resi[steps];

  for(int i=0;i<steps;i++)
    {
      sum_component->Reset();
      sum_component->Add(rho_component,bez_rho_component,(double)i/(double)steps,(double)(steps-i)/(double)steps);
      sum_component->Add(bg_component);

      double chi=sum_component->Chi2Test(data_component,"CHI2");
      x[i]=(double)i/(double)steps;
      y[i]=chi;

      cout<<"n: "<<i<<" x: "<<x[i]<<" y: "<<y[i]<<endl;
    }

  TGraph* chi_test=new TGraph(steps,x,y);

  TCanvas *cChi2=new TCanvas("cChi2");
  chi_test->Draw("AC*");

  TCanvas *cOptimal=new TCanvas("cOptimal");
  double ww=0;
  sum_component->Add(rho_component,bez_rho_component,ww,1.0-ww);
  sum_component->Add(bg_component);
  data_component->Draw();
  sum_component->Draw("same");
  bg_component->Draw("same");
  rho_component->Scale(ww);
  rho_component->Draw("same");
  bez_rho_component->Scale(1.0-ww);
  bez_rho_component->Draw("same");
}
