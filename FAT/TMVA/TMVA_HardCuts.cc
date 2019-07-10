
//TMVAClassification.cc

#include <cstdlib>
#include <iostream>

#include "TChain.h"
#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TObjString.h"
#include "TSystem.h"
#include "TROOT.h"

#include "TMVA/Factory.h"
//#include "TMVA/DataLoader.h"
#include "TMVA/Tools.h"
//#include "TMVA/TMVAGui.h"

//#include "ttextfile.h"

using namespace std;

//______________________________________________________________________________
void TMVA_HardCuts(TString extraSuffix = "_RecCuts", Long64_t DesEntries = -1) {
  TMVA::Tools::Instance();

  TString treeFile="tree";
  
  cout << "==> Start TMVAClassification" << endl;

  TString NameSuffix = (treeFile.Contains("_") ? treeFile(treeFile.First('_'), treeFile.Last('.') - treeFile.First('_')) : TString("New")) + extraSuffix;
  TFile* outputFile = TFile::Open("TMVATraining" + NameSuffix + ".root", "RECREATE");

  TMVA::Factory* factory = new TMVA::Factory("TMVAClassification_MLP_BDT_HC" + NameSuffix, outputFile, "!V:!Silent:Color:DrawProgressBar:AnalysisType=Classification");

  //factory->AddVariable("p_p",&p_beta);
  //factory->AddVariable("pip_p",&pip_p);
  //factory->AddVariable("dist_p_pim", "dist_p_pim","mm",'F',0,250);
  factory->AddVariable("dist_pip_pim", "dist_pip_pim","mm",'F',0,250);
  //factory->AddVariable("eVert_x", &eVert_x);
  //factory->AddVariable("eVert_y", &eVert_y);
  //factory->AddVariable("eVert_z", &eVert_z);
  //factory->AddVariable("ver_pip_pim_x","ver_pip_pim_x","mm",'F',-500,500);
  //factory->AddVariable("ver_pip_pim_y","ver_pip_pim_y","mm",'F',-500,500);
  factory->AddVariable("ver_pip_pim_z","ver_pip_pim_z","mm",'F',-500,500);
  //factory->AddVariable("ver_p_pim_x","ver_p_pim_x","mm",'F',-500,500);
  //factory->AddVariable("ver_p_pim_y","ver_p_pim_y","mm",'F',-500,500);
  factory->AddVariable("ver_p_pim_z","ver_p_pim_z","mm",'F',-500,500);
  //factory->AddVariable("oa_lambda", &oa_lambda);
  //factory->AddVariable("oa_pip_p",&oa_pip_p);
  //factory->AddVariable("lambda_mom_z",&lambda_mom_z);
  //factory->AddVariable("dist_p_eVert",&dist_p_eVert);
  //factory->AddVariable("dist_pim_eVert",&dist_pim_eVert);
  //factory->AddVariable("dist_lambda_eVert",&dist_lambda_eVert);
  factory->AddVariable("dist_lambda_ver_pip_pim","dist_lambda_ver_pip_pim","mm",'F',0,800);
  factory->AddVariable("dist_ver_to_ver","dist_ver_to_ver","mm",'F',0,1200);
  //factory->AddVariable("m_inv_pip_pim","m_inv_pip_pim","MeV",'F',0,1400);
  //factory->AddVariable("miss_mass_kp","miss_mass_kp","MeV",'F',0,1800);

  
  //#warning Momentum diabled!

  TFile* input    = TFile::Open("hardCuts_signal.root", "UPDATE");
  TFile* input2   = TFile::Open("hardCuts_background.root", "UPDATE");
  TTree* tSigAll  = (TTree*) input->Get("signal");
  TTree* tBackAll = (TTree*) input->Get("background");
  TTree* tBackData= (TTree*) input2->Get("background_data");
  if(tBackData==0)
    cout<<"uninicialized pointer!"<<endl;
  
  Long64_t MaxEntries = TMath::Min(TMath::Min(tSigAll->GetEntries(), tBackAll->GetEntries()),tBackData->GetEntries());
//if (DesEntries > 0 && DesEntries < MaxEntries)
//  MaxEntries = DesEntries;
  
  //TTree* tSig  = tSigAll ->CloneTree(MaxEntries);
  TTree* tSig  = tSigAll ->CloneTree();
  //TTree* tBack = tBackAll->CloneTree(MaxEntries);
  TTree* tBack = tBackAll->CloneTree();
  
  //TTree* tBackDataSel =tBackData-> CloneTree(MaxEntries);
  TTree* tBackDataSel =tBackData-> CloneTree();

  factory->AddSignalTree    (tSig,  1.);
  //factory->AddBackgroundTree(tBack, 1);
  factory->AddBackgroundTree(tBackDataSel, 1);
  //factory->SetSignalWeightExpression("Weight");
  //factory->SetBackgroundWeightExpression("Weight");

  factory->PrepareTrainingAndTestTree("", "", "!V:SplitMode=Random:SplitSeed=0:NormMode=EqualNumEvents");

  // Default
  //factory->BookMethod(factory, TMVA::Types::kMLP, "MLP", "H:!V:NCycles=500:HiddenLayers=N,N-1:NeuronType=sigmoid:EstimatorType=MSE:VarTransform=N:TestRate=10:!UseRegulator");
  // Shzymon
  //factory->BookMethod(factory, TMVA::Types::kMLP, "MLP", "H:!V:NCycles=500:HiddenLayers=N+1,N:NeuronType=sigmoid:EstimatorType=CE:VarTransform=N:TestRate=5:!UseRegulator");

  // MLP
  //factory->BookMethod(TMVA::Types::kMLP, "kMLP_pca_mse_500", "!H:!V:NCycles=500:HiddenLayers=N+1,N:NeuronType=sigmoid:NeuronInputType=sum:EstimatorType=MSE:TrainingMethod=BP:VarTransform=N,P:BPMode=sequential");

  //factory->BookMethod(TMVA::Types::kMLP, "kMLP_pca_mse_1000", "!H:!V:NCycles=1000:HiddenLayers=N+1,N:NeuronType=sigmoid:NeuronInputType=sum:EstimatorType=MSE:TrainingMethod=BP:VarTransform=N,P:BPMode=sequential");

  //factory->BookMethod(TMVA::Types::kMLP, "kMLP_pca_ce_1000_n4_no_ev", "!H:!V:NCycles=1000:HiddenLayers=N,N+5,N+5,N:NeuronType=sigmoid:NeuronInputType=sum:EstimatorType=CE:TrainingMethod=BP:VarTransform=N,P:BPMode=sequential");

  //factory->BookMethod(TMVA::Types::kMLP, "kMLP_ce_600_n2_no_ev", "!H:!V:NCycles=600:HiddenLayers=N,N:NeuronType=sigmoid:NeuronInputType=sum:EstimatorType=CE:TrainingMethod=BP:VarTransform=N:BPMode=sequential:CalculateErrors=True");

  //factory->BookMethod(TMVA::Types::kMLP, "kMLP_pca_ce_600_n_no_ev", "!H:!V:NCycles=600:HiddenLayers=N:NeuronType=sigmoid:NeuronInputType=sum:EstimatorType=CE:TrainingMethod=BP:VarTransform=N,P:BPMode=sequential:CalculateErrors=True");

  //factory->BookMethod( TMVA::Types::kLikelihood, "Likelihood", "VarTransform=N,P" );

  factory->BookMethod( TMVA::Types::kBDT, "BDT", "UseYesNoLeaf=False:nCuts=40:VarTransform=N" );
  
  factory->BookMethod( TMVA::Types::kCuts,"RecCuts_GA","!V:FitMethod=GA");
  factory->BookMethod( TMVA::Types::kCuts,"RecCuts_SA","!V:FitMethod=SA");
  factory->BookMethod( TMVA::Types::kCuts,"RecCuts_MC","!V:FitMethod=MC");

  factory->TrainAllMethods();
  factory->TestAllMethods();
  factory->EvaluateAllMethods();

  outputFile->cd();
  //(new TTextFile(Form("%s/TMVAClassification.cc", gSystem->pwd()), "TMVAClassification_cc"))->Write("TMVAClassification_cc", TObject::kOverwrite);

  outputFile->Close();

  cout << "==> Wrote root file: " << outputFile->GetName() << endl;
  cout << "==> TMVAClassification is done!" << endl;

  delete factory;
  //delete factory;
}
