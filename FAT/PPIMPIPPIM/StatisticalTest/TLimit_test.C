
void TLimit_test()
{
TFile* infile=new TFile("Event.root","READ");
infile->cd();
TH1* sh=(TH1*)infile->Get("signal");
TH1* bh=(TH1*)infile->Get("background");
TH1* dh=(TH1*)infile->Get("data");
TLimitDataSource* mydatasource = new TLimitDataSource(sh,bh,dh);
TConfidenceLevel *myconfidence = TLimit::ComputeLimit(mydatasource,50000);
std::cout << "  CLs    : " << myconfidence->CLs()  << std::endl;
std::cout << "  CLsb   : " << myconfidence->CLsb() << std::endl;
std::cout << "  CLb    : " << myconfidence->CLb()  << std::endl;
std::cout << "< CLs >  : " << myconfidence->GetExpectedCLs_b()  << std::endl;
std::cout << "< CLsb > : " << myconfidence->GetExpectedCLsb_b() << std::endl;
std::cout << "< CLb >  : " << myconfidence->GetExpectedCLb_b()  << std::endl;
delete myconfidence;
delete mydatasource;
infile->Close();
}
