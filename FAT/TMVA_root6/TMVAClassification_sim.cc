
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
void TMVAClassification_sim(TString extraSuffix = "new_Vertex", Long64_t DesEntries = -1) {
  TMVA::Tools::Instance();

  cout << "==> Start TMVAClassification" << endl;

  //TString NameSuffix = (treeFile.Contains("_") ? treeFile(treeFile.First('_'), treeFile.Last('.') - treeFile.First('_')) : TString("New")) + extraSuffix;
  TFile* outputFile = TFile::Open("Output_sim" + extraSuffix + ".root", "RECREATE");

  if(outputFile==0)
    cout<<"uninicialized output file"<<endl;
  else
    cout<<"output file set"<<endl;
  
  TMVA::Factory* factory = new TMVA::Factory("TMVAClassification_data_driven" + extraSuffix, outputFile, "!V:!Silent:Color:DrawProgressBar:Transformations=N:AnalysisType=Classification");
  cout<<"Factory set"<<endl;
  TMVA::DataLoader* dataloader = new TMVA::DataLoader("dataset");

  //dataloader->AddVariable("oa_lambda","oa_lambda","deg",'F',0,180);
  //dataloader->AddVariable("p_p","p_beta","mm",'F');
  //dataloader->AddVariable("pip_p","pip_p","mm",'F');
  dataloader->AddVariable("dist_p_pim", "dist_p_pim","mm",'F',0,150);
  dataloader->AddVariable("dist_pip_pim", "dist_pip_pim","mm",'F',0,150);
  //dataloader->AddVariable("eVert_x", "eVert_x","mm",'F',-60,60);
  //dataloader->AddVariable("eVert_y", "eVert_y","mm",'F',-60,60);
  //dataloader->AddVariable("eVert_z", "eVert_z","mm",'F',-100,40);
  dataloader->AddVariable("ver_pip_pim_x","ver_pip_pim_x","mm",'F',-80,80);
  dataloader->AddVariable("ver_pip_pim_y","ver_pip_pim_y","mm",'F',-80,80);
  dataloader->AddVariable("ver_pip_pim_z","ver_pip_pim_z","mm",'F',-150,150);
  dataloader->AddVariable("ver_p_pim_x","ver_p_pim_x","mm",'F',-100,100);
  dataloader->AddVariable("ver_p_pim_y","ver_p_pim_y","mm",'F',-100,100);
  dataloader->AddVariable("ver_p_pim_z","ver_p_pim_z","mm",'F',-150,150);
  dataloader->AddVariable("oa_lambda", "oa_lambda","deg",'F',0,180);
  //dataloader->AddVariable("oa_pip_p","oa_pip_p","mm",'F');
  //dataloader->AddVariable("lambda_mom_z","lambda_mom_z","mm",'F');
  //dataloader->AddVariable("dist_p_eVert","dist_p_eVert","mm",'F',0,140);
  //dataloader->AddVariable("dist_pim_eVert","dist_pim_eVert","mm",'F',0,140);
  //dataloader->AddVariable("dist_lambda_eVert","dist_lambda_eVert","mm",'F',0,140);
  dataloader->AddVariable("dist_lambda_ver_pip_pim","dist_lambda_ver_pip_pim","mm",'F',0,140);
  dataloader->AddVariable("dist_ver_to_ver","dist_ver_to_ver","mm",'F',0,150);
  
  
  //#warning Momentum diabled!

  TFile* input1   = TFile::Open("signal_from_sim.root","UPDATE");
  TFile* input2   = TFile::Open("background_from_sim.root","UPDATE");
  cout<<"load input file"<<endl;
  //TTree* tSigAll  = (TTree*) input->Get("signal");
  //TTree* tBackAll = (TTree*) input->Get("background");
  TTree* tBackData= (TTree*) input2->Get("background");
  TTree* tSignalData= (TTree*) input1->Get("signal");
  if(tBackData==0)
    cout<<"uninicialized background pointer!"<<endl;
  else
    cout<<"background data load"<<endl;

  if(tSignalData==0)
    cout<<"uninicialized signal pointer!"<<endl;
  else
    cout<<"signal data load"<<endl;
  
  
  Long64_t MaxEntries = TMath::Min(tSignalData->GetEntries(),tBackData->GetEntries());
//if (DesEntries > 0 && DesEntries < MaxEntries)
//  MaxEntries = DesEntries;
  
  TTree* tSig  = tSignalData -> CloneTree();
  TTree* tBack = tBackData -> CloneTree();
  
  
  dataloader->AddSignalTree    (tSig,  1.);
  dataloader->AddBackgroundTree(tBack, 1);
  
  dataloader->PrepareTrainingAndTestTree("", "", "!V:SplitMode=Random:SplitSeed=0:NormMode=EqualNumEvents");


  factory->BookMethod(dataloader,TMVA::Types::kCuts,"RecCuts","!V:VarTransform=N:FitMethod=GA");
  factory->BookMethod(dataloader,TMVA::Types::kKNN, "kNN_30", "nkNN=30:VarTransform=N" );
  factory->BookMethod(dataloader,TMVA::Types::kKNN, "kNN_20", "nkNN=20:VarTransform=N" );
  factory->BookMethod(dataloader,TMVA::Types::kKNN, "kNN_10", "nkNN=10:VarTransform=N" );
  //factory->BookMethod(dataloader,TMVA::Types::kPDERS, "PDERS", "VarTransform=N,P" );
  //factory->BookMethod(dataloader,TMVA::Types::kSVM, "SVM", "VarTransform=N" );

  //factory->BookMethod( TMVA::Types::kBDT, "BDT", "UseYesNoLeaf=False:nCuts=40" );
   
  
  //factory->BookMethod( TMVA::Types::kDNN, "DNN_d_0.1", "Layout:RELU|N,RELU|N,TANH:TrainingStrategy = LearningRate=1e-1, BatchSize=256| LearningRate=1e-2, BatchSize=256| LearningRate=1e-3, BatchSize=256:DropConfig=0.1:VarTransform=N,P");
  //factory->BookMethod( TMVA::Types::kDNN, "DNN_d_0.2", "Layout:RELU|N,RELU|N,TANH:TrainingStrategy = LearningRate=1e-1, BatchSize=256| LearningRate=1e-2, BatchSize=256| LearningRate=1e-3, BatchSize=256:DropConfig=0.2:VarTransform=N,P");
  //factory->BookMethod( TMVA::Types::kDNN, "DNN_d_0.0", "Layout:RELU|N,RELU|N,TANH:TrainingStrategy = LearningRate=1e-1, BatchSize=256| LearningRate=1e-2, BatchSize=256| LearningRate=1e-3, BatchSize=256:DropConfig=0.0:VarTransform=N,P");
  // Keras - TensorFlow
  //factory->BookMethod(factory, TMVA::Types::kPyKeras, "PyKeras", "H:!V:VarTransform=N:FilenameModel=model.h5:NumEpochs=20:BatchSize=32");

  // Comparison
  /*factory->BookMethod( factory, TMVA::Types::kPDEFoam, "PDEFoam",
    "!H:!V:SigBgSeparate=F:TailCut=0.001:VolFrac=0.0666:nActiveCells=500:nSampl=2000:nBin=5:Nmin=100:Kernel=None:Compress=T" );
    factory->BookMethod( factory, TMVA::Types::kPDEFoam, "PDEFoamBoost",
    "!H:!V:Boost_Num=30:Boost_Transform=linear:SigBgSeparate=F:MaxDepth=4:UseYesNoCell=T:DTLogic=MisClassificationError:FillFoamWithOrigWeights=F:TailCut=0:nActiveCells=500:nBin=20:Nmin=400:Kernel=None:Compress=T" );

    factory->BookMethod( factory, TMVA::Types::kMLP, "MLP", "H:!V:NeuronType=tanh:VarTransform=N:NCycles=600:HiddenLayers=N+5:TestRate=5:!UseRegulator" );
    factory->BookMethod( factory, TMVA::Types::kMLP, "MLPBFGS", "H:!V:NeuronType=tanh:VarTransform=N:NCycles=600:HiddenLayers=N+5:TestRate=5:TrainingMethod=BFGS:!UseRegulator" );
    factory->BookMethod( factory, TMVA::Types::kMLP, "MLPBNN", "H:!V:NeuronType=tanh:VarTransform=N:NCycles=60:HiddenLayers=N+5:TestRate=5:TrainingMethod=BFGS:UseRegulator" ); // BFGS training with bayesian regulators

    
    // Multi-core CPU implementation.
    TString cpuOptions = dnnOptions + ":Architecture=CPU";
    factory->BookMethod(factory, TMVA::Types::kDNN, "DNN_CPU", cpuOptions);

    factory->BookMethod( factory, TMVA::Types::kCFMlpANN, "CFMlpANN", "!H:!V:NCycles=200:HiddenLayers=N+1,N"  ); // n_cycles:#nodes:#nodes:...
    factory->BookMethod( factory, TMVA::Types::kTMlpANN, "TMlpANN", "!H:!V:NCycles=200:HiddenLayers=N+1,N:LearningMethod=BFGS:ValidationFraction=0.3"  ); // n_cycles:#nodes:#nodes:...

    factory->BookMethod( factory, TMVA::Types::kSVM, "SVM", "Gamma=0.25:Tol=0.001:VarTransform=Norm" );

    factory->BookMethod( factory, TMVA::Types::kBDT, "BDTG",
    "!H:!V:NTrees=1000:MinNodeSize=2.5%:BoostType=Grad:Shrinkage=0.10:UseBaggedBoost:BaggedSampleFraction=0.5:nCuts=20:MaxDepth=2" );
    factory->BookMethod( factory, TMVA::Types::kBDT, "BDT",
    "!H:!V:NTrees=850:MinNodeSize=2.5%:MaxDepth=3:BoostType=AdaBoost:AdaBoostBeta=0.5:UseBaggedBoost:BaggedSampleFraction=0.5:SeparationType=GiniIndex:nCuts=20" );
    factory->BookMethod( factory, TMVA::Types::kBDT, "BDTB",
    "!H:!V:NTrees=400:BoostType=Bagging:SeparationType=GiniIndex:nCuts=20" );
    factory->BookMethod( factory, TMVA::Types::kBDT, "BDTD",
    "!H:!V:NTrees=400:MinNodeSize=5%:MaxDepth=3:BoostType=AdaBoost:SeparationType=GiniIndex:nCuts=20:VarTransform=Decorrelate" );
    factory->BookMethod( factory, TMVA::Types::kBDT, "BDTF",
    "!H:!V:NTrees=50:MinNodeSize=2.5%:UseFisherCuts:MaxDepth=3:BoostType=AdaBoost:AdaBoostBeta=0.5:SeparationType=GiniIndex:nCuts=20" );

    factory->BookMethod( factory, TMVA::Types::kRuleFit, "RuleFit",
    "H:!V:RuleFitModule=RFTMVA:Model=ModRuleLinear:MinImp=0.001:RuleMinDist=0.001:NTrees=20:fEventsMin=0.01:fEventsMax=0.5:GDTau=-1.0:GDTauPrec=0.01:GDStep=0.01:GDNSteps=10000:GDErrScale=1.02" );
    */
  factory->TrainAllMethods();
  factory->TestAllMethods();
  factory->EvaluateAllMethods();

  outputFile->cd();
  //(new TTextFile(Form("%s/TMVAClassification.cc", gSystem->pwd()), "TMVAClassification_cc"))->Write("TMVAClassification_cc", TObject::kOverwrite);

  outputFile->Close();

  cout << "==> Wrote root file: " << outputFile->GetName() << endl;
  cout << "==> TMVAClassification is done!" << endl;

  delete factory;
  delete dataloader;
}
