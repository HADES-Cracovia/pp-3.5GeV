
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
#include "TMVA/DataLoader.h"
#include "TMVA/Tools.h"
#include "TMVA/TMVAGui.h"


using namespace std;

//______________________________________________________________________________
void TMVAClassification_data_driven_multiple(TString extraSuffix = "new_Vertex", Long64_t DesEntries = -1) {

  TMVA::Tools::Instance();
  cout << "==> Start TMVAClassification" << endl;

  TFile* outputFile = TFile::Open("TMVATraining_data_driven_miss_mass" + extraSuffix + ".root", "RECREATE");

  if(outputFile==0)
    cout<<"uninicialized output file"<<endl;
  else
    cout<<"output file set"<<endl;
  
  TMVA::Factory* factory = new TMVA::Factory("TMVAClassification_data_driven" + extraSuffix, outputFile, "!V:!Silent:Color:DrawProgressBar:Transformations=N:AnalysisType=Classification");
  cout<<"Factory set"<<endl;

  TMVA::DataLoader* dataloader = new TMVA::DataLoader("dataset");
  dataloader->AddVariable("dist_p_pim", "dist_p_pim","mm",'F',0,150);
  dataloader->AddVariable("dist_pip_pim", "dist_pip_pim","mm",'F',0,150);
  dataloader->AddVariable("ver_pip_pim_x","ver_pip_pim_x","mm",'F',-80,80);
  dataloader->AddVariable("ver_pip_pim_y","ver_pip_pim_y","mm",'F',-80,80);
  dataloader->AddVariable("ver_pip_pim_z","ver_pip_pim_z","mm",'F',-150,150);
  dataloader->AddVariable("ver_p_pim_x","ver_p_pim_x","mm",'F',-100,100);
  dataloader->AddVariable("ver_p_pim_y","ver_p_pim_y","mm",'F',-100,100);
  dataloader->AddVariable("ver_p_pim_z","ver_p_pim_z","mm",'F',-150,150);
  dataloader->AddVariable("oa_lambda", "oa_lambda","deg",'F',0,180);
  dataloader->AddVariable("dist_lambda_ver_pip_pim","dist_lambda_ver_pip_pim","mm",'F',0,140);
  dataloader->AddVariable("dist_ver_to_ver","dist_ver_to_ver","mm",'F',0,150);
  
  

  TFile* input2   = TFile::Open("input_from_data_miss_mass_4_new_vertex.root","UPDATE");
  cout<<"load input file"<<endl;

  TTree* tBackData= (TTree*) input2->Get("background_data");
  TTree* tSignalData= (TTree*) input2->Get("signal_data");

  if(tBackData==0)
    cout<<"uninicialized pointer!"<<endl;
  else
    cout<<"background data load"<<endl;

  Long64_t MaxEntries = TMath::Min(tSignalData->GetEntries(),tBackData->GetEntries());
  
  TTree* tSig  = tSignalData -> CloneTree();
  TTree* tBack = tBackData -> CloneTree();
  

  dataloader->AddSignalTree    (tSig,  1.);
  dataloader->AddBackgroundTree(tBack, 1);

  dataloader->PrepareTrainingAndTestTree("", "", "V:SplitMode=Random:SplitSeed=0:NormMode=EqualNumEvents");
  
  factory->BookMethod(dataloader,TMVA::Types::kMLP, "kMLP_pca_ce_600_n2_no_ev", "!H:!V:NCycles=600:HiddenLayers=N,N:NeuronType=sigmoid:NeuronInputType=sum:EstimatorType=CE:TrainingMethod=BP:VarTransform=N,P:BPMode=sequential:CalculateErrors=True");
  factory->BookMethod(dataloader,TMVA::Types::kLikelihood, "Likelihood", "VarTransform=N,P" );
  factory->BookMethod(dataloader,TMVA::Types::kBDT, "BDT", "UseYesNoLeaf=False:nCuts=40" );
  factory->BookMethod(dataloader,TMVA::Types::kCuts,"RecCuts","!V:FitMethod=GA");
  
  factory->TrainAllMethods();
  factory->TestAllMethods();
  factory->EvaluateAllMethods();

  outputFile->cd();
  outputFile->Close();

  cout << "==> Wrote root file: " << outputFile->GetName() << endl;
  cout << "==> TMVAClassification is done!" << endl;

  delete factory;
}
