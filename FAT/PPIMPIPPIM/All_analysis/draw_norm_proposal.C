void mult_error(TH1* hist, double err,bool verbose=0)
{
  cout<<endl<<"***scaling the histogram by const value***"<<endl;
  cout<<"hist name"<<hist->GetName()<<endl;
  cout<<"scaling error "<<err<<endl;
  
  int nbin_min=1;
  int nbin_max=hist->GetNbinsX();
  for(int i=nbin_min;i<=nbin_max;i++)
    {
      //double cont=hist->GetBinContent(i);
      double error=hist->GetBinError(i)*err;
      if(verbose)
	{
	  cout<<"bin number: "<<i<<" bin contetnt: "<<cont<<endl;
	  cout<<"                       bin error:  "<<error<<endl;
	}
      hist->SetBinError(i,error); 
    }

  cout<<"***end of scale_error function***"<<endl<<endl;
  
}


void set_Y_name(TH1* hist)
{
  char name[10000]; // enough to hold all numbers up to 64-bits
  sprintf(name, "#frac{counts}{%.1f MeV}", hist->GetBinWidth(2));
  cout<<"Y axis name: "<<name<<endl;
  hist->GetYaxis()->SetTitle(name);
}

void setHistogramStyleData(TH1* hist)
{
  hist->SetLineWidth(2);
  
  hist->SetMarkerColor(hist->GetLineColor());
  hist->SetMarkerSize(2);
  hist->SetMarkerStyle(8);
  set_Y_name(hist);

  //hist->GetXaxis()->SetLabelFont(42);
  hist->GetXaxis()->SetNdivisions(508);
  //hist->GetXaxis()->SetLabelSize(0.05);
  // hist->GetXaxis()->SetTitleSize(0.05);
  //hist->GetXaxis()->SetTitleOffset(1.1);
  //hist->GetXaxis()->SetTitleFont(42);

  hist->GetYaxis()->SetNdivisions(508);
  //hist->GetYaxis()->SetLabelFont(42);
  //hist->GetYaxis()->SetLabelSize(0.05);
  //hist->GetYaxis()->SetTitleOffset(0.8);
  //hist->GetYaxis()->SetTitleSize(0.05);
  //hist->GetYaxis()->SetTitleFont(42);  
  
}

void setHistogramStyleSimul(TH1* hist)
{
  hist->SetLineWidth(2);

  hist->SetMarkerColor(hist->GetLineColor());
  hist->SetMarkerSize(1);
  hist->SetMarkerStyle(8);
  hist->SetFillColor(hist->GetLineColor());
  set_Y_name(hist);

  //hist->GetXaxis()->SetLabelFont(42);
  hist->GetXaxis()->SetNdivisions(508);
  //hist->GetXaxis()->SetLabelSize(0.05);
  //hist->GetXaxis()->SetTitleSize(0.05);
  //hist->GetXaxis()->SetTitleOffset(1.1);
  //hist->GetXaxis()->SetTitleFont(42);

  hist->GetYaxis()->SetNdivisions(508);
  //hist->GetYaxis()->SetLabelFont(42);
  //hist->GetYaxis()->SetLabelSize(0.05);
  //hist->GetYaxis()->SetTitleOffset(0.8);
  //hist->GetYaxis()->SetTitleSize(0.05);
  //hist->GetYaxis()->SetTitleFont(42);  
  
}

void setLineStyle(TLine* line)
{
  line->SetLineWidth(4);
  line->SetLineStyle(9);
  line->SetLineColor(kRed-3);
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

void scale_error(TH1* hist, double err,bool verbose=0)
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
    if(verbose)
    {
    cout<<"bin number: "<<i<<" bin contetnt: "<<cont<<endl;
    cout<<"                       bin error:  "<<error<<endl;
    }
    hist->SetBinError(i,error); 
    }

    cout<<"***end of scale_error function***"<<endl<<endl;
    /*
  cout<<endl<<"***scaling the histogram errors according expected statistical errors***"<<endl;
  cout<<"hist name"<<hist->GetName()<<endl;
  cout<<"scaling error"<<err<<endl;
  hist->Sumw2();
  cout<<"***end of scale_error function***"<<endl<<endl;
    */
}

void scale(TH1* hist,double sc)
{
  for (Int_t j=1; j<hist->GetNbinsX()+1; ++j)
    {
      //double scale=1.0/(3.13 * TMath::Power(10,8)) *1000; /*to get micro barn*/
      double binErr=hist->GetBinError(j);
      hist->SetBinContent(j, hist->GetBinContent(j)*sc);
      //hist->SetBinError( j, TMath::Sqrt( hist->GetBinContent(j) ) );
      hist->SetBinError(j, binErr*sc);
      //hist->GetYaxis()->SetTitle("#frac{dN}{dE} [#frac{#mu b}{MeV}]");
    }

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


int draw_norm_proposal(void)
{

  TFile *fileS1385 = new TFile("SB_sim_S1385pK0.root","READ");
  TFile *fileSDpp = new TFile("SB_sim_SDppK0.root","READ");
  TFile *fileLDpp = new TFile("SB_sim_LDppK0.root","READ");
  TFile *fileL1520= new TFile("SB_sim_L1520pippim.root","READ");
  TFile *fileLK0=new TFile("SB_sim_LK0ppip.root","READ");
  TFile *fileExp= new TFile("SB_experiment.root","READ");
  TFile *fileRho=new TFile("z_RHO.root","READ");
  TFile *fileBezRho=new TFile("bez_RHO.root","READ");

  
  TH1F *hexperiment_L=(TH1F*)fileExp->Get("hMPPim_TMVA_K0mass");
  TH1F *hexperiemnt_K0=(TH1F*)fileExp->Get("hMPipPim_TMVA_Lmass");
  TH1F *hsim_L=(TH1F*)fileLK0->Get("hMPPim_TMVA_K0mass");
  TH1F *hsim_K0=(TH1F*)fileLK0->Get("hMPipPim_TMVA_Lmass");
  TF1 *fK0_experiment_fit=(TF1 *)fileExp->Get("K0_fit");
  TF1 *fK0_experiment_sig=(TF1 *)fileExp->Get("K0_signal");
  TF1 *fL1116_experiment_fit=(TF1 *)fileExp->Get("L1116_fit");
  TF1 *fL1116_experiment_sig=(TF1 *)fileExp->Get("L1116_signal");
  TH1F *hexperiment_SB_spectrum=(TH1F*)fileExp->Get("orginal_spectrum");

  
  TF1* fVoigt_bg=(TF1*)fileExp->Get("fVoigt_bg");
  TF1* fVoigt=(TF1*)fileExp->Get("fVoigt");
  TF1* fbg=(TF1*)fileExp->Get("fbg");
  TF1 *voigt=new TF1("signal_voit","[0]*TMath::Voigt(x-[1],[2],[3],4)",1380,1750); //only for pure signal

  TLine* line1=(TLine*)fileExp->Get("line1");
  TLine* line2=(TLine*)fileExp->Get("line2");
  TLine* line3=(TLine*)fileExp->Get("line3");
  TLine* line4=(TLine*)fileExp->Get("line4");

  
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

  TH1F *hS1385_hMPipPim_signal = (TH1F*)fileS1385->Get("hMPipPim_signal");
  hS1385_hMPipPim_signal->SetName("hS1385_hMPipPim_signal");
  hS1385_hMPipPim_signal->Sumw2(kFALSE);
  TH1F *hSDpp_hMPipPim_signal = (TH1F*)fileSDpp->Get("hMPipPim_signal");
  hSDpp_hMPipPim_signal->SetName("hSDpp_hMPipPim_signal");
  hSDpp_hMPipPim_signal->Sumw2(kFALSE);;
  TH1F *hLDpp_hMPipPim_signal = (TH1F*)fileLDpp->Get("hMPipPim_signal");
  hLDpp_hMPipPim_signal->SetName("hLDpp_hMPipPim_signal");
  hLDpp_hMPipPim_signal->Sumw2(kFALSE);
  TH1F *hexperiment_hMPipPim_signal=(TH1F*)fileExp->Get("hMPipPim_signal");
  hexperiment_hMPipPim_signal->SetName("hexperiment_hMPipPim_signal");
  //hexperiment_hMPipPim_signal->Sumw2();
  TH1F *hL1520_hMPipPim_signal=(TH1F*)fileL1520->Get("hMPipPim_signal");
  hL1520_hMPipPim_signal->SetName("hL1520_hMPipPim_signal");
  hL1520_hMPipPim_signal->Sumw2(kFALSE);
  TH1F *hsum_hMPipPim_signal=(TH1F*)hS1385_hMPipPim_signal->Clone("hsum_hMPipPim_signal");
  hsum_hMPipPim_signal->Reset();

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

  TH1F *hS1385_hMPipPim_background = (TH1F*)fileS1385->Get("hMPipPim_background");
  hS1385_hMPipPim_background->SetName("hS1385_hMPipPim_background");
  hS1385_hMPipPim_background->Sumw2(kFALSE);
  TH1F *hSDpp_hMPipPim_background = (TH1F*)fileSDpp->Get("hMPipPim_background");
  hSDpp_hMPipPim_background->SetName("hSDpp_hMPipPim_background");
  hSDpp_hMPipPim_background->Sumw2(kFALSE);
  TH1F *hLDpp_hMPipPim_background = (TH1F*)fileLDpp->Get("hMPipPim_background");
  hLDpp_hMPipPim_background->SetName("hLDpp_hMPipPim_background");
  hLDpp_hMPipPim_background->Sumw2(kFALSE);
  TH1F *hexperiment_hMPipPim_background=(TH1F*)fileExp->Get("hMPipPim_background");
  hexperiment_hMPipPim_background->SetName("hexperiment_hMPipPim_background");
  //hexperiment_hMPipPim_background->Sumw2();
  TH1F *hL1520_hMPipPim_background=(TH1F*)fileL1520->Get("hMPipPim_background");
  hL1520_hMPipPim_background->SetName("hL1520_hMPipPim_background");
  hL1520_hMPipPim_background->Sumw2(kFALSE);
  //RHO file
  TH1F *hL1520_hMPipPim_rho_distribution=(TH1F*)fileRho->Get("histo1");
  hL1520_hMPipPim_rho_distribution->SetName("hL1520_hMPipPim_rho_distribution");
  TH1F *hL1520_hMPipPim_bez_rho_distribution=(TH1F*)fileBezRho->Get("histo1");
  hL1520_hMPipPim_bez_rho_distribution->SetName("hL1520_hMPipPim_bez_rho_distribution");

  TH1F *hL1520_hMPipPim_rho=(TH1F*)hL1520_hMPipPim_signal->Clone("hL1520_hMPipPim_rho");
  hL1520_hMPipPim_rho->Reset();
  hL1520_hMPipPim_rho->SetNameTitle("hL1520_hMPipPim_rho","M^{inv}_{#pi^{+} #pi^{-}} for #Lambda(1520);M^{inv}_{#pi^{+} #pi^{-}}[MeV];counts");
  TH1F *hL1520_hMPipPim_bez_rho=(TH1F*)hL1520_hMPipPim_signal->Clone("hL1520_hMPipPim_bez_rho");
  hL1520_hMPipPim_bez_rho->Reset();
  hL1520_hMPipPim_bez_rho->SetNameTitle("hL1520_hMPipPim_bez_rho","M^{inv}_{#pi^{+} #pi^{-}} for #Lambda(1520);M^{inv}_{#pi^{+} #pi^{-}}[MeV];counts");

  TH1F *hL1520_hMPipPim_efficiency=(TH1F*)hL1520_hMPipPim_background->Clone("hL1520_hMPipPim_efficiency");
  hL1520_hMPipPim_efficiency->Reset();
  hL1520_hMPipPim_efficiency->SetNameTitle("hL1520_hMPipPim_efficiency","Shape distortion by the limited acceptance;M^{inv}_{#pi^{+} #pi^{-}}[MeV];#frac{reconstructed}{simulated in 4 #pi}");
  //end of HRO files 
  
  TH1F *hsum_background=(TH1F*)hS1385_background->Clone("hsum_background");
  TH1F *hclean_background=(TH1F*)hS1385_background->Clone("hclean_background");
  TH1F *hclean_experiment=(TH1F*)hexperiment_background->Clone("hclean_experiment");
  TH1F *hclean_L1520=(TH1F*)hL1520_background->Clone("hclean_L1520");
  TH1F *hclean_sum=(TH1F*)hL1520_background->Clone("hclean_sum");
  TH1F *hclean_L1520_ren=(TH1F*)hL1520_background->Clone("hclean_L1520_ren");
  TH1F *hclean_sum_ren=(TH1F*)hL1520_background->Clone("hclean_sum_ren");
  TH1F *hpure_signal=(TH1F*)hL1520_background->Clone("hpure_signal");

  TH1F *hsum_background_PipPim=(TH1F*)hS1385_hMPipPim_background->Clone("hsum_background_PipPim");
  TH1F *hclean_background_PipPim=(TH1F*)hS1385_hMPipPim_background->Clone("hclean_background_PipPim");
  TH1F *hclean_experiment_PipPim=(TH1F*)hexperiment_hMPipPim_background->Clone("hclean_experiment_PipPim");
  TH1F *hclean_L1520_PipPim=(TH1F*)hL1520_hMPipPim_background->Clone("hclean_L1520_PipPim");
  TH1F *hclean_sum_PipPim=(TH1F*)hL1520_hMPipPim_background->Clone("hclean_sum_PipPim");
  TH1F *hclean_L1520_ren_PipPim=(TH1F*)hL1520_hMPipPim_background->Clone("hclean_L1520_ren_PipPim");
  TH1F *hclean_sum_ren_PipPim=(TH1F*)hL1520_hMPipPim_background->Clone("hclean_sum_ren_PipPim");
  TH1F *hpure_signal_PipPim=(TH1F*)hL1520_hMPipPim_background->Clone("hpure_signal_PipPim");

  hsum_background->Reset();
  hclean_background->Reset();
  hclean_experiment->Reset();
  hclean_L1520->Reset();
  hclean_sum->Reset();
  hclean_L1520_ren->Reset();
  hclean_sum_ren->Reset();
  hpure_signal->Reset();

  
  hsum_background_PipPim->Reset();
  hclean_background_PipPim->Reset();
  hclean_experiment_PipPim->Reset();
  hclean_L1520_PipPim->Reset();
  hclean_sum_PipPim->Reset();
  hclean_L1520_ren_PipPim->Reset();
  hclean_sum_ren_PipPim->Reset();
  hpure_signal_PipPim->Reset();

  //scale according to CS
  //double nsim=40*TMath::Power(10,6);//number of simulated events
  double nsim=120*TMath::Power(10,6);
  double scale=3.13*TMath::Power(10,8);
  double downscale=3;//trigger downscale for simulated events comensate by sim_factor
  double sim_factor=2.64*10e4*28/1400*370/430;//fit to prdicted signal area
  double cs[5]=
    {14.05/1000*scale/(nsim*downscale)*sim_factor,//S1385
     9.26/1000*scale/(nsim*downscale)*sim_factor,//SDpp
     29.45/1000*scale/(nsim*downscale)*sim_factor,//LDpp
     5.6/1000*scale/(100*100000*downscale)*sim_factor,//L(1520)pK+->Lpi+pi-pK+
     (2.57+14.05+9.26+29.45+5.0+3.5+2.3+14)/1000*scale/(100*100000*downscale)*0.5*sim_factor//L K0 p pi+ (0.5 because of Ks i Kl)
    };
  double err[4]={0,0,0,0};
    /*{2.25/14.05,//S1385
     1.47/9.26,//SDpp
     2.55/29.45,//LDpp
     2/5.6//L(1520)pK+->Lpi+pi-pK+
     };*/
  double cs_sig;
  // cs in \mu barns, have to me re-calculated to mb!!

  
  
  scale(hS1385_background,cs[0]);
  scale(hSDpp_background,cs[1]);
  scale(hLDpp_background,cs[2]);
  scale(hL1520_background,cs[3]);
  scale(hS1385_hMPipPim_background,cs[0]);
  scale(hSDpp_hMPipPim_background,cs[1]);
  scale(hLDpp_hMPipPim_background,cs[2]);
  scale(hL1520_hMPipPim_background,cs[3]);
  scale(hsim_L,cs[4]);
  scale(hsim_K0,cs[4]);
  scale(hexperiment_background,sim_factor);
  scale(hexperiment_hMPipPim_background,sim_factor);
  /*
  scale_error(hS1385_background,err[0]);
  scale_error(hSDpp_background,err[1]);
  scale_error(hLDpp_background,err[2]);
  scale_error(hL1520_background,err[3]);
  scale_error(hS1385_hMPipPim_background,err[0]);
  scale_error(hSDpp_hMPipPim_background,err[1]);
  scale_error(hLDpp_hMPipPim_background,err[2]);
  scale_error(hL1520_hMPipPim_background,err[3]);
  */
  
  
  hS1385_background->Sumw2();
  hSDpp_background->Sumw2();
  hLDpp_background->Sumw2();
  hL1520_background->Sumw2();
  hS1385_hMPipPim_background->Sumw2();
  hSDpp_hMPipPim_background->Sumw2();
  hLDpp_hMPipPim_background->Sumw2();
  hL1520_hMPipPim_background->Sumw2();
  hexperiment_background->Sumw2();
  hexperiment_hMPipPim_background->Sumw2();
  //hsim_L->Sumw2();
  //hsim_K0->Sumw2();

  /*  hS1385_background->Sumw2();
      hSDpp_background->Sumw2();
      hLDpp_background->Sumw2();
      hL1520_background->Sumw2();
  */
  
  scale(hS1385_data,cs[0]);
  scale(hSDpp_data,cs[1]);
  scale(hLDpp_data,cs[2]);
  scale(hL1520_data,cs[3]);
  scale(hS1385_hMPipPim_signal,cs[0]);
  scale(hSDpp_hMPipPim_signal,cs[1]);
  scale(hLDpp_hMPipPim_signal,cs[2]);
  scale(hL1520_hMPipPim_signal,cs[3]);
  scale(hexperiment_data,sim_factor);
  scale(hexperiment_hMPipPim_signal,sim_factor);
  /*
  scale_error(hS1385_data,err[0]);
  scale_error(hSDpp_data,err[1]);
  scale_error(hLDpp_data,err[2]);
  scale_error(hL1520_data,err[3]);
  scale_error(hS1385_hMPipPim_signal,err[0]);
  scale_error(hSDpp_hMPipPim_signal,err[1]);
  scale_error(hLDpp_hMPipPim_signal,err[2]);
  scale_error(hL1520_hMPipPim_signal,err[3]);
  */
  
  hS1385_data->Sumw2();
  hSDpp_data->Sumw2();
  hLDpp_data->Sumw2();
  hL1520_data->Sumw2();
  hS1385_hMPipPim_signal->Sumw2();
  hSDpp_hMPipPim_signal->Sumw2();
  hLDpp_hMPipPim_signal->Sumw2();
  hL1520_hMPipPim_signal->Sumw2();
  hexperiment_data->Sumw2();
  hexperiment_hMPipPim_signal->Sumw2();

  /*
    hS1385_data->Sumw2();
    hSDpp_data->Sumw2();
    hLDpp_data->Sumw2();
    hL1520_data->Sumw2();
  */
  hsum_background->Add(hS1385_background);
  hsum_background->Add(hSDpp_background);
  hsum_background->Add(hLDpp_background);
  hsum_background_PipPim->Add(hS1385_hMPipPim_background);
  hsum_background_PipPim->Add(hSDpp_hMPipPim_background);
  hsum_background_PipPim->Add(hLDpp_hMPipPim_background);
 
  hsum_data->Add(hS1385_data);
  hsum_data->Add(hSDpp_data);
  hsum_data->Add(hLDpp_data);
  hsum_hMPipPim_signal->Add(hS1385_hMPipPim_signal);
  hsum_hMPipPim_signal->Add(hSDpp_hMPipPim_signal);
  hsum_hMPipPim_signal->Add(hLDpp_hMPipPim_signal);

  hclean_background->Add(hsum_data,hsum_background,1,-1);
  hclean_experiment->Add(hexperiment_data,hexperiment_background,1,-1);
  hclean_L1520->Add(hL1520_data,hL1520_background,1,-1);
  hclean_sum->Add(hclean_L1520,hclean_background,1,1);
  //hclean_sum->Add(hclean_experiment,hclean_background,1,-1);
  
  hclean_background_PipPim->Add(hsum_hMPipPim_signal,hsum_background_PipPim,1,-1);
  hclean_experiment_PipPim->Add(hexperiment_hMPipPim_signal,hexperiment_hMPipPim_background,1,-1);
  hclean_L1520_PipPim->Add(hL1520_hMPipPim_signal,hL1520_hMPipPim_background,1,-1);
  hclean_sum_PipPim->Add(hclean_L1520_PipPim,hclean_background_PipPim,1,1);
  //hclean_sum_PipPim->Add(hclean_experiment_PipPim,hclean_background_PipPim,1,-1);
  //cs_sig=1/(hclean_L1520->Integral())*20;
  //hclean_L1520->Scale(cs_sig);


  //scale signal to difference between signal and background
  double int_min=1410;
  double int_max=1600;
  double err_sum;
  
  double sig_int=hclean_L1520->Integral(hclean_L1520->FindBin(int_min),hclean_L1520->FindBin(int_max));
  double backgroud_int=hclean_background->Integral(hclean_background->FindBin(int_min),hclean_background->FindBin(int_max));
  double experiment_int=hclean_experiment->Integral(hclean_experiment->FindBin(int_min),hclean_experiment->FindBin(int_max));
  double scale_sig=5.6/7.2;//hard code from experiment part

  
  hclean_L1520_ren->Add(hclean_L1520,1);
  //hclean_L1520_ren->Scale(scale_sig);
  hclean_L1520_ren->Scale((experiment_int-backgroud_int)/sig_int);
  hclean_sum_ren->Add(hclean_L1520_ren,1);
  hclean_sum_ren->Add(hclean_background,1);
  //hclean_sum_ren->Add(hclean_experiment,hclean_background,1,-1);
  
  hclean_L1520_ren_PipPim->Add(hclean_L1520_PipPim,1);
  //hclean_L1520_ren_PipPim->Scale(scale_sig);
  hclean_L1520_ren_PipPim->Scale((experiment_int-backgroud_int)/sig_int);
  hclean_sum_ren_PipPim->Add(hclean_L1520_ren_PipPim,1);
  hclean_sum_ren_PipPim->Add(hclean_background_PipPim,1);
  //hclean_sum_ren_PipPim->Add(hclean_experiment_PipPim,hclean_background_PipPim,1,-1);
  
  int rebin_res=1;  
  TCanvas *cRes=new TCanvas("cRes","cRes");
  cRes->Divide(2,2);
  cRes->cd(1);
  hS1385_data->Rebin(rebin_res);
  set_Y_name(hS1385_data);
  hS1385_data->SetAxisRange(1350,1800);
  hS1385_data->Draw("e1");
  hS1385_background->Rebin(rebin_res);
  hS1385_background->SetLineColor(kRed);
  hS1385_background->Draw("samee1");

  cRes->cd(2);
  hSDpp_data->Rebin(rebin_res);
  set_Y_name(hSDpp_data);
  hSDpp_data->SetAxisRange(1350,1800);
  hSDpp_data->Draw("e1");
  hSDpp_background->SetLineColor(kRed);
  hSDpp_background->Rebin(rebin_res);
  hSDpp_background->SetAxisRange(1350,1800); 
  hSDpp_background->Draw("samee1");

  cRes->cd(3);
  hLDpp_data->Rebin(rebin_res);
  set_Y_name(hLDpp_data);
  hLDpp_data->SetAxisRange(1350,1800);
  hLDpp_data->Draw("e1");
  hLDpp_background->SetLineColor(kRed);
  hLDpp_background->Rebin(rebin_res);
  hLDpp_background->SetAxisRange(1350,1800);
  hLDpp_background->Draw("samee1");

  cRes->cd(4);
  hL1520_data->Rebin(rebin_res);
  set_Y_name(hL1520_data);
  hL1520_data->SetAxisRange(1350,1800);
  hL1520_data->Draw("e1");  
  hL1520_background->SetLineColor(kRed);
  hL1520_background->Rebin(rebin_res);
  hL1520_background->SetAxisRange(1350,1800);
  hL1520_background->Draw("samee1");

  TCanvas *cRes_PipPim=new TCanvas("cRes_PipPim","cRes_PipPim");
  cRes_PipPim->Divide(2,2);
  cRes_PipPim->cd(1);
  hS1385_hMPipPim_signal->Rebin(rebin_res);
  set_Y_name(hS1385_hMPipPim_signal);
  hS1385_hMPipPim_signal->SetAxisRange(250,450);
  hS1385_hMPipPim_signal->Draw("e1");
  hS1385_hMPipPim_background->Rebin(rebin_res);
  hS1385_hMPipPim_background->SetLineColor(kRed);
  hS1385_hMPipPim_background->Draw("samee1");

  cRes_PipPim->cd(2);
  hSDpp_hMPipPim_signal->Rebin(rebin_res);
  set_Y_name(hSDpp_hMPipPim_signal);
  hSDpp_hMPipPim_signal->SetAxisRange(250,450);
  hSDpp_hMPipPim_signal->Draw("e1");
  hSDpp_hMPipPim_background->SetLineColor(kRed);
  hSDpp_hMPipPim_background->Rebin(rebin_res);
  //hSDpp_hMPipPim_background->SetAxisRange(1350,1800); 
  hSDpp_hMPipPim_background->Draw("samee1");

  cRes_PipPim->cd(3);
  hLDpp_hMPipPim_signal->Rebin(rebin_res);
  set_Y_name(hLDpp_hMPipPim_signal);
  hLDpp_hMPipPim_signal->SetAxisRange(250,450);
  hLDpp_hMPipPim_signal->Draw("e1");
  hLDpp_hMPipPim_background->SetLineColor(kRed);
  hLDpp_hMPipPim_background->Rebin(rebin_res);
  //hLDpp_hMPipPim_background->SetAxisRange(1350,1800);
  hLDpp_hMPipPim_background->Draw("samee1");

  cRes_PipPim->cd(4);
  hL1520_hMPipPim_signal->Rebin(rebin_res);
  set_Y_name(hL1520_hMPipPim_signal);
  hL1520_hMPipPim_signal->SetAxisRange(250,450);
  hL1520_hMPipPim_signal->Draw("e1");  
  hL1520_hMPipPim_background->SetLineColor(kRed);
  hL1520_hMPipPim_background->Rebin(rebin_res);
  //hL1520_hMPipPim_background->SetAxisRange(1350,1800);
  hL1520_hMPipPim_background->Draw("samee1");  
  
  int rebin2=2; //wrong error propagation for simul events
  TCanvas *cSum=new TCanvas("cSum","cSum");  
  
  hexperiment_data->Rebin(rebin2);
  hexperiment_data->SetAxisRange(1300,1800);
  hexperiment_data->Draw("e1");
  setHistogramStyleData(hexperiment_data);
  
  hsum_data->Rebin(rebin2);
  hsum_data->SetAxisRange(1300,1800);
  hsum_data->Draw("samee2");
  setHistogramStyleSimul(hsum_data);
 
  
  hsum_background->Rebin(rebin2);
  hsum_background->SetLineColor(kRed);
  setHistogramStyleSimul(hsum_background);
  hsum_background->Draw("samee2");

  
  hexperiment_background->SetAxisRange(1300,1800);
  hexperiment_background->Rebin(rebin2);
  hexperiment_background->SetLineColor(kRed);
  hexperiment_background->Draw("samee2");
  setHistogramStyleData(hexperiment_background);
  
   
  
  TCanvas *cSum_PipPim=new TCanvas("cSum_PipPim","cSumPipPim");  
  
  hexperiment_hMPipPim_signal->Rebin(rebin2);
  hexperiment_hMPipPim_signal->SetAxisRange(250,450);
  hexperiment_hMPipPim_signal->Draw("e1");
  setHistogramStyleData(hexperiment_hMPipPim_signal);
    
  hexperiment_hMPipPim_background->Rebin(rebin2);
  hexperiment_hMPipPim_background->SetLineColor(kRed);
  hexperiment_hMPipPim_background->Draw("samee2");
  setHistogramStyleData(hexperiment_hMPipPim_background);
  
  hsum_hMPipPim_signal->Rebin(rebin2);
  hsum_hMPipPim_signal->SetAxisRange(250,450);
  hsum_hMPipPim_signal->Draw("samee2");
  setHistogramStyleSimul(hsum_hMPipPim_signal);
  
  hsum_background_PipPim->Rebin(rebin2);
  hsum_background_PipPim->SetLineColor(kRed);
  setHistogramStyleSimul(hsum_background_PipPim);
  hsum_background_PipPim->Draw("samee2");

  
  int rebin=4;
  TCanvas *cClean=new TCanvas("cClean","cClean");
  
  hclean_experiment->Draw("e1");
  hclean_experiment->Rebin(rebin);
  hclean_experiment->GetXaxis()->SetRangeUser(1360,1780);
  setHistogramStyleData(hclean_experiment);
  
  hclean_sum->Rebin(rebin);
  hclean_sum->GetXaxis()->SetRangeUser(1360,1780);
  hclean_sum->SetLineColor(kMagenta);
  hclean_sum->SetFillStyle(3145);
  hclean_sum->Draw("samee2");
  setHistogramStyleSimul(hclean_sum);


  hclean_background->SetLineColor(kRed);
  hclean_background->SetFillColor(kRed);
  hclean_background->Rebin(rebin);
  hclean_background->Draw("samee2");
  setHistogramStyleSimul(hclean_background);
  
  hclean_L1520->SetLineColor(kGreen+3);
  hclean_L1520->Rebin(rebin);
  hclean_L1520->Draw("samee2");
  hclean_L1520->SetFillStyle(3154);
  setHistogramStyleSimul(hclean_L1520);
  
  int rebin_pippim=4;
  TCanvas *cClean_PipPim=new TCanvas("cClean_PipPim","cClean_PipPim");
  
  hclean_experiment_PipPim->Draw("e1");
  hclean_experiment_PipPim->Rebin(rebin_pippim);
  hclean_experiment_PipPim->GetXaxis()->SetRangeUser(250,450);
  setHistogramStyleData(hclean_experiment_PipPim);
  
  hclean_sum_PipPim->Rebin(rebin_pippim);
  hclean_sum_PipPim->GetXaxis()->SetRangeUser(250,450);
  hclean_sum_PipPim->SetLineColor(kMagenta);
  hclean_sum_PipPim->SetFillStyle(3145);
  hclean_sum_PipPim->Draw("samee2");
  setHistogramStyleSimul(hclean_sum_PipPim);

  hclean_background_PipPim->SetLineColor(kRed);
  hclean_background_PipPim->SetFillColor(kRed);
  hclean_background_PipPim->Rebin(rebin_pippim);
  hclean_background_PipPim->Draw("samee2");
  setHistogramStyleSimul(hclean_background_PipPim);
  
  hclean_L1520_PipPim->SetLineColor(kGreen+3);
  hclean_L1520_PipPim->Rebin(rebin_pippim);
  hclean_L1520_PipPim->Draw("samee2");
  hclean_L1520_PipPim->SetFillStyle(3154);
  setHistogramStyleSimul(hclean_L1520_PipPim);

   
  TCanvas *cClean_ren=new TCanvas("cClean_ren","cClean_ren");
  cClean_ren->Divide(1);
  cClean_ren->cd(1);

  hclean_experiment->Draw("e1");

  hclean_sum_ren->Rebin(rebin);
  hclean_sum_ren->SetLineColor(kMagenta);
  hclean_sum_ren->SetFillColor(kMagenta);
  setHistogramStyleSimul(hclean_sum_ren);
  hclean_sum_ren->SetFillStyle(3145);
  hclean_sum_ren->Draw("samee2");
     
  //hclean_experiment->Rebin(rebin);
  //hclean_background->SetLineColor(kRed);
  //hclean_background->Rebin(rebin);
  hclean_background->Draw("samee2");
  setHistogramStyleSimul(hclean_background);
  hclean_background->SetFillStyle(3125);
  
  hclean_L1520_ren->SetLineColor(kGreen+3);
  hclean_L1520_ren->Rebin(rebin);
  setHistogramStyleSimul(hclean_L1520_ren);
  hclean_L1520_ren->SetFillStyle(3154);
  hclean_L1520_ren->Draw("samee2");
  

  
  /*
  cClean_ren->cd(2);
  hpure_signal->Rebin(rebin);
  hpure_signal->GetXaxis()->SetRangeUser(1360,1780);
  hpure_signal->Draw("e1");
    
  //fit Voigt to data
  hpure_signal->Add(hclean_experiment,hclean_background,1,-1);
  setHistogramStyleData(hpure_signal);
  voigt->SetParameter(0,2412);
  voigt->SetParameter(1,1500);
  voigt->SetParameter(2,5);
  voigt->SetParameter(3,50);
  hpure_signal->Fit(voigt,"RL");
  hpure_signal->Fit(voigt,"RL");

  TLatex *printFormula1 = new TLatex();
  double high=0.85;
  char text4[10000];
  char text5[10000];
  char text6[10000];
  char text7[10000];
  char text8[10000];
  sprintf(text4, "#sigma = %.2f MeV",voigt->GetParameter(2));
  sprintf(text5, "#Gamma = %.2f MeV",voigt->GetParameter(3));
  sprintf(text6, "#bar{M_{p #pi^{-} #pi^{+} #pi^{-}}} = %.1f MeV",voigt->GetParameter(1));
  sprintf(text7, "#int_{1400 MeV}^{1620 MeV} = %.2f ",voigt->Integral(1400,1620)/hpure_signal->GetBinWidth(2));
  sprintf(text8, "#sum_{1400 MeV}^{1620 MeV} = %.2f ",hpure_signal->Integral(hpure_signal->FindBin(1400),hpure_signal->FindBin(1620)));
  printFormula1->SetNDC();
  printFormula1->SetTextFont(32);
  printFormula1->SetTextColor(1);
  printFormula1->SetTextSize(0.04);
  printFormula1->SetTextAlign(13);
  printFormula1->DrawLatex(0.5,high,text4);
  printFormula1->DrawLatex(0.5,high-printFormula1->GetTextSize(),text5);
  printFormula1->DrawLatex(0.5,high-printFormula1->GetTextSize()*2,text6);
  printFormula1->DrawLatex(0.5,high-printFormula1->GetTextSize()*4,text7);
  printFormula1->DrawLatex(0.5,high-printFormula1->GetTextSize()*7,text8);
  */
  
  TCanvas *cClean_ren_PipPim=new TCanvas("cClean_ren_PipPim","cClean_ren_PipPim");
  hclean_experiment_PipPim->Draw("e1");
  
  hclean_sum_ren_PipPim->Rebin(rebin_pippim);
  hclean_sum_ren_PipPim->SetLineColor(kMagenta);
  hclean_sum_ren_PipPim->SetFillColor(kMagenta);
  setHistogramStyleSimul(hclean_sum_ren_PipPim);
  hclean_sum_ren_PipPim->SetFillStyle(3145);
  hclean_sum_ren_PipPim->Draw("samee2");
      
  //hclean_experiment->Rebin(rebin_pippim);
  //hclean_background->SetLineColor(kRed);
  //hclean_background->Rebin(rebin_pippim);
  hclean_background_PipPim->Draw("samee2");
  setHistogramStyleSimul(hclean_background_PipPim);
  hclean_background_PipPim->SetFillStyle(3125);
  
  hclean_L1520_ren_PipPim->SetLineColor(kGreen+3);
  hclean_L1520_ren_PipPim->Rebin(rebin_pippim);
  setHistogramStyleSimul(hclean_L1520_ren_PipPim);
  hclean_L1520_ren_PipPim->SetFillStyle(3154);
  hclean_L1520_ren_PipPim->Draw("samee2");
  
  TCanvas *cClean_ren_PipPim_rho=new TCanvas("cClean_ren_PipPim_rho","cClean_ren_PipPim_rho");
  hclean_experiment_PipPim->Draw("e1");
    
  //hclean_experiment->Rebin(rebin_pippim);
  //hclean_background->SetLineColor(kRed);
  //hclean_background->Rebin(rebin_pippim);
  
  hclean_sum_ren_PipPim->Draw("samee2");
  hclean_background_PipPim->Draw("samee2");
  hclean_L1520_ren_PipPim->Draw("samee2");
  
  int rr=100000;
  for(int r=0;r<rr;r++)
    {
      hL1520_hMPipPim_rho->Fill(hL1520_hMPipPim_rho_distribution->GetRandom()*1000);
      hL1520_hMPipPim_bez_rho->Fill(hL1520_hMPipPim_bez_rho_distribution->GetRandom()*1000);
    }

  //rebin the same as data
  hL1520_hMPipPim_rho->Rebin(rebin_pippim);
  hL1520_hMPipPim_bez_rho->Rebin(rebin_pippim);
  hL1520_hMPipPim_efficiency->Rebin(rebin_pippim);

  hL1520_hMPipPim_efficiency->Divide(hclean_L1520_ren_PipPim, hL1520_hMPipPim_bez_rho);
  hL1520_hMPipPim_rho->Multiply(hL1520_hMPipPim_efficiency);
  hL1520_hMPipPim_rho->Scale(hclean_L1520_ren_PipPim->Integral()/hL1520_hMPipPim_rho->Integral());
    

  hL1520_hMPipPim_rho->Sumw2();
  mult_error(hL1520_hMPipPim_rho,0.9);//scale up according to assumption that the same values has the same eror 
  
  hL1520_hMPipPim_rho->Draw("samee2");
  hL1520_hMPipPim_rho->SetFillStyle(3195);
  hL1520_hMPipPim_rho->SetFillColor(kGreen-4);
  hL1520_hMPipPim_rho->SetLineColor(kGreen-4);
  setHistogramStyleSimul(hL1520_hMPipPim_rho);
  
  
  TCanvas *cSB=new TCanvas("cSB","Spectrum for side-band");
  hexperiment_SB_spectrum->SetAxisRange(1050,1250);
  setHistogramStyleData(hexperiment_SB_spectrum);
  hexperiment_SB_spectrum->Draw();
  fVoigt->Draw("same");
  fVoigt->SetLineColor(kGreen);
  fbg->Draw("same");
  fbg->SetLineColor(kBlue);
  setLineStyle(line1);
  setLineStyle(line2);
  setLineStyle(line3);
  setLineStyle(line4);
  line1->Draw("same");
  line2->Draw("same");
  line3->Draw("same");
  line4->Draw("same");

  TLatex *printFormula = new TLatex();
  double high=0.8;
  char text1[10000];
  char text2[10000];
  char text3[10000];
  sprintf(text1, "#sigma = %.1f MeV",fVoigt->GetParameter(2));
  sprintf(text2, "#Gamma = %.1f MeV",fVoigt->GetParameter(3));
  sprintf(text3, "#bar{M_{p#pi^{-}}} = %.1f MeV",fVoigt->GetParameter(1));
  printFormula->SetNDC();
  printFormula->SetTextFont(32);
  printFormula->SetTextColor(1);
  printFormula->SetTextSize(0.05);
  printFormula->SetTextAlign(13);
  printFormula->DrawLatex(0.6,high, text1);
  printFormula->DrawLatex(0.6,high-printFormula->GetTextSize(), text2);
  printFormula->DrawLatex(0.6,high-printFormula->GetTextSize()*2 , text3);

  
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

  TCanvas* cLK0=new TCanvas("cLK0", "Signal for final state p #pi^{+} #L^{0} K^{0}");
  int npx=300;
  cLK0->Divide(2);

  cLK0->cd(1);
  hexperiment_L->SetAxisRange(1080,1200);
  hexperiment_L->SetMinimum(0);
  hexperiment_L->Draw("e1");
  hsim_L->Draw("samee1");
  hsim_L->SetLineColor(kMagenta+1);
  fL1116_experiment_fit->SetNpx(npx);
  fL1116_experiment_sig->SetNpx(npx);
  fL1116_experiment_fit->Draw("same");
  fL1116_experiment_sig->Draw("same");
  fL1116_experiment_sig->SetLineWidth(3);
  
  cLK0->cd(2);
  hexperiemnt_K0->SetAxisRange(300,650);
  hexperiemnt_K0->SetMinimum(0);
  hexperiemnt_K0->Draw("e1");  
  hsim_K0->Draw("samee1");
  hsim_K0->SetLineColor(kMagenta+1);
  fK0_experiment_sig->SetNpx(npx);
  fK0_experiment_fit->SetNpx(npx);
  fK0_experiment_fit->Draw("same");
  fK0_experiment_sig->Draw("same");
  fK0_experiment_sig->SetLineWidth(3);

  TCanvas* cL=new TCanvas("cL", "Signal for final state p #pi^{+} #L^{0} K^{0}");
  TF1* L_sim_bg=new TF1("L_sim_bg","pol1(0)",1095,1135);
  TF1* L_sim_sig_bg=new TF1("L_sim_sig_bg","[0]*TMath::Voigt(x-[1],[2],[3])+pol1(4)",1095,1135);
  TF1* fL1116_experiment_fit_bg=new TF1("fL1116_experiment_fit_bg","pol2(0)",1095,1135);
  fL1116_experiment_fit_bg->SetParameters(fL1116_experiment_fit->GetParameter(4),fL1116_experiment_fit->GetParameter(5),fL1116_experiment_fit->GetParameter(6));
  L_sim_sig_bg->SetParameters(1494,1115,0.0004,5.4,47.3,-0.033);
  hexperiment_L->SetAxisRange(1080,1200);
  hexperiment_L->SetMinimum(0);
  hexperiment_L->Draw("e1");
  hsim_L->Draw("samee1");
  hsim_L->SetLineColor(kMagenta+1);
  L_sim_bg->SetLineColor(kMagenta);
  L_sim_bg->SetLineWidth(3);
  L_sim_bg->SetLineStyle(2);
  L_sim_sig_bg->SetLineColor(kMagenta);
  L_sim_sig_bg->SetLineWidth(3);
  L_sim_sig_bg->SetLineStyle(1);
  L_sim_sig_bg->SetNpx(npx);

  hsim_L->Fit(L_sim_sig_bg,"R");
  L_sim_bg->SetParameters(L_sim_sig_bg->GetParameter(4),L_sim_sig_bg->GetParameter(5));
  L_sim_bg->Draw("same");

  //fL1116_experiment_sig->SetNpx(npx);
  fL1116_experiment_fit->SetLineWidth(3);
  fL1116_experiment_fit->SetLineColor(kBlue);
  fL1116_experiment_fit->Draw("same");
  //fL1116_experiment_sig->Draw("same");
  //fL1116_experiment_sig->SetLineWidth(3);

  fL1116_experiment_fit_bg->SetLineWidth(3);
  fL1116_experiment_fit_bg->SetLineColor(kBlue);
  fL1116_experiment_fit_bg->SetLineStyle(2);
  fL1116_experiment_fit_bg->Draw("same");

  TLatex *printFormula2 = new TLatex();
  double high2=0.90;
  char text9[10000];
  char text10[10000];
  char text11[10000];
  sprintf(text9, "I_{exp}=#int_{1095}^{1135} S_{exp} = %.1f",fL1116_experiment_sig->Integral(1095,1135)/hsim_L->GetBinWidth(3));
  sprintf(text10, "I_{simul}=#int_{1095}^{1135} S_{simul} = %.1f",(L_sim_sig_bg->Integral(1095,1135)-L_sim_bg->Integral(1095,1135))/hsim_L->GetBinWidth(3));
  sprintf(text11, "I_{exp}/I_{simul} = %.2f",fL1116_experiment_sig->Integral(1095,1135)/(L_sim_sig_bg->Integral(1095,1135)-L_sim_bg->Integral(1095,1135)));
  printFormula2->SetNDC();
  printFormula2->SetTextFont(32);
  printFormula2->SetTextColor(1);
  printFormula2->SetTextSize(0.05);
  printFormula2->SetTextAlign(13);
  printFormula2->DrawLatex(0.6,high2, text9);
  printFormula2->DrawLatex(0.6,high2-printFormula->GetTextSize()*4, text10);
  printFormula2->DrawLatex(0.6,high2-printFormula->GetTextSize()*8 , text11);

  
  TCanvas* cK0=new TCanvas("cK0", "Signal for final state p #pi^{+} #L^{0} K^{0}");
  TF1* K0_sim_bg=new TF1("K0_sim_bg","pol1(0)",450,550);
  TF1* K0_sim_sig_bg=new TF1("K0_sim_sig_bg","[0]*TMath::Voigt(x-[1],[2],[3])+pol1(4)",450,550);
  TF1* fK0_experiment_fit_bg=new TF1("fK0_experiment_fit_bg","pol2(0)",450,550);
  fK0_experiment_fit_bg->SetParameters(fK0_experiment_fit->GetParameter(4),fK0_experiment_fit->GetParameter(5),fK0_experiment_fit->GetParameter(6));
  K0_sim_sig_bg->SetParameters(1490,495,3.21,9.54,60.71,-0.1);
  hexperiemnt_K0->SetAxisRange(300,650);
  hexperiemnt_K0->SetMinimum(0);
  hexperiemnt_K0->Draw("e1");  
  hsim_K0->Draw("samee1");
  hsim_K0->SetLineColor(kMagenta+1);

  K0_sim_bg->SetLineColor(kMagenta);
  K0_sim_bg->SetLineWidth(4);
  K0_sim_bg->SetLineStyle(2);
  K0_sim_sig_bg->SetLineColor(kMagenta);
  K0_sim_sig_bg->SetLineWidth(4);
  K0_sim_sig_bg->SetLineStyle(1);
  K0_sim_sig_bg->SetNpx(npx);

  hsim_K0->Fit(K0_sim_sig_bg,"R");
  K0_sim_bg->SetParameters(K0_sim_sig_bg->GetParameter(4),K0_sim_sig_bg->GetParameter(5));
  K0_sim_bg->Draw("same");

  fK0_experiment_fit->SetLineWidth(4);
  fK0_experiment_fit->SetLineColor(kBlue);
  fK0_experiment_fit->Draw("same");
  //fK0_experiment_sig->Draw("same");
  //fK0_experiment_sig->SetLineWidth(3);

  fK0_experiment_fit_bg->SetLineWidth(4);
  fK0_experiment_fit_bg->SetLineColor(kBlue);
  fK0_experiment_fit_bg->SetLineStyle(2);
  fK0_experiment_fit_bg->Draw("same");

  TLatex *printFormula3 = new TLatex();
  double high3=0.90;
  char text12[10000];
  char text13[10000];
  char text14[10000];
  sprintf(text12, "I_{exp}=#int_{450}^{550} S_{exp} = %.1f",fK0_experiment_sig->Integral(450,550)/hsim_K0->GetBinWidth(3));
  sprintf(text13, "I_{simul}=#int_{450}^{550} S_{simul} = %.1f",(K0_sim_sig_bg->Integral(450,550)-K0_sim_bg->Integral(450,550))/hsim_L->GetBinWidth(3));
  sprintf(text14, "I_{exp}/I_{simul} = %.2f",fK0_experiment_sig->Integral(450,550)/(K0_sim_sig_bg->Integral(450,550)-K0_sim_bg->Integral(450,550)));
  printFormula3->SetNDC();
  printFormula3->SetTextFont(32);
  printFormula3->SetTextColor(1);
  printFormula3->SetTextSize(0.05);
  printFormula3->SetTextAlign(13);
  printFormula3->DrawLatex(0.6,high2, text12);
  printFormula3->DrawLatex(0.6,high2-printFormula->GetTextSize()*4, text13);
  printFormula3->DrawLatex(0.6,high2-printFormula->GetTextSize()*8 , text14);

  
  //save all
  TFile* output=new TFile("final_output_proposal.root","recreate");

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

  hexperiemnt_K0->Write();
  hexperiment_L->Write();
  hsim_K0->Write();
  hsim_L->Write();

  fK0_experiment_fit->Write();
  fK0_experiment_sig->Write();
  fL1116_experiment_fit->Write();
  fL1116_experiment_sig->Write();

  
  line1->Write();
  line2->Write();
  line3->Write();
  line4->Write();

  fbg->Write();
  fVoigt_bg->Write();
  fVoigt->Write();

  cClean_ren->Write();
  cClean_ren_PipPim->Write();
  cRes->Write();
  cRes_PipPim->Write();
  cClean->Write();
  cClean_PipPim->Write();
  cClean_ren_PipPim_rho->Write();
  cSum->Write();
  cSB->Write();
  cLK0->Write();
  cK0->Write();
  cL->Write();
}

