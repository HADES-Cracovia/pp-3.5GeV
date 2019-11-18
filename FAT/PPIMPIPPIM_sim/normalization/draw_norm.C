int draw_norm(void)
{
  TFile *fileS1385 = new TFile("./../TMVAeval_DD/pp_after_TMVA_DD_6(n+4)_pNb_newVertex_S1385pK0_Rafal.root","READ");
  TFile *fileSDpp = new TFile("./../TMVAeval_DD/pp_after_TMVA_DD_6(n+4)_pNb_newVertex_SDppK0_Rafal.root","READ");
  TFile *fileLDpp = new TFile("./../TMVAeval_DD/pp_after_TMVA_DD_6(n+4)_pNb_newVertex_LDppK0_Rafal.root","READ");

  TFile *output= new TFile("output.root","RECREATE");

  TH1F * hS1385_data = (TH1F*)f.Get("");
  h1->Draw();

  
}
