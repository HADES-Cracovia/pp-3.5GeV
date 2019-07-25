double calcchi2( TF1* fuc, TH1F* hist, double xmin=0, double xmax=-1)
{
  return calcchi2(hist, fuc, xmin, xmax);
}

double calcchi2(TH1F* hist, TF1* fuc, double xmin=0, double xmax=-1)
{
  int nbin=0;
  int binmin=0;
  int binmax=0;
  double chi2=0;
  if(xmax>xmin)
    {
      binmax=hist->GetBin(xmax);
      binmin=hist->GetBin(xmin);
      nbin=binmax-binmin+1;
    }
  else
    {
      minmin=0;
      binmax=hist->GetNBinsX()-1;
      nbin=hist->GetNBinsX()
    }

  for(int i=binmin;i<binmax;i++)
    {
      chi2=chi2+TMath::Pow(hist->GetBinContent(i)-fuc(fuc->GetBinCenter(i)),2);
    }
  return chi2/nbin;
}



void TLimit_test()
{
  TFile* infile=new TFile("Event.root","READ");
  infile->cd();
  double chi2;
  TH1* sh=(TH1*)infile->Get("signal");
  TH1* bh=(TH1*)infile->Get("background");
  TH1* dh=(TH1*)infile->Get("data");
  TLimitDataSource* mydatasource = new TLimitDataSource(sh,bh,dh);
  TConfidenceLevel *myconfidence = TLimit::ComputeLimit(mydatasource,500000);

  std::cout << "  CLs    : " << myconfidence->CLs()  << std::endl;
  std::cout << "  CLsb   : " << myconfidence->CLsb() << std::endl;
  std::cout << "  CLb    : " << myconfidence->CLb()  << std::endl;
  std::cout << "< CLs >  : " << myconfidence->GetExpectedCLs_b()  << std::endl;
  std::cout << "< CLsb > : " << myconfidence->GetExpectedCLsb_b() << std::endl;
  std::cout << "< CLb >  : " << myconfidence->GetExpectedCLb_b()  << std::endl;

  myconfidence->Draw();
  delete myconfidence;
  delete mydatasource;
  infile->Close();

  //chi2 test
  double fit_min=1450;
  double fit_max=1520;
  
  TF1* bg=new TF1("bg","pol3(0)",fit_min,fit_max);
  dh->Fit(bg,"r");

  chi2=calcchi2(dh,bg,fit_min,fit_max);
  cout<<"chi2 by my program= "<<chi2<<end;
}
