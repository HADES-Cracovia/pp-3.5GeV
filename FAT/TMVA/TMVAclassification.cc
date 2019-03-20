TMVAClassification.cc

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

//#include "ttextfile.h"

using namespace std;

//______________________________________________________________________________
void TMVAClassification(TString treeFile, TString extraSuffix = "", Long64_t DesEntries = -1) {
  TMVA::Tools::Instance();

  cout << "==> Start TMVAClassification" << endl;

  TString NameSuffix = (treeFile.Contains("_") ? treeFile(treeFile.First('_'), treeFile.Last('.') - treeFile.First('_')) : TString("New")) + extraSuffix;
  TFile* outputFile = TFile::Open("TMVATraining" + NameSuffix + ".root", "RECREATE");

  TMVA::Factory* factory = new TMVA::Factory("TMVAClassification" + NameSuffix, outputFile, "!V:!Silent:Color:DrawProgressBar:Transformations=N:AnalysisType=Classification");
  TMVA::DataLoader* dataloader = new TMVA::DataLoader("dataset");

  dataloader->AddVariable("dist_p_pim", "dist_p_pim", "mm",    'F');
  dataloader->AddVariable("dist_pip_pim", "dist_pip_pim", "mm",    'F');
  dataloader->AddVariable("eVert_x",  "eVert_x",  "mm",    'F');
  dataloader->AddVariable("eVert_y",  "eVert_y",  "mm",    'F');
  dataloader->AddVariable("eVert_z",  "eVert_z",  "mm",    'F');
  dataloader->AddVariable("dist_ver_to_ver",     "dist_ver_to_ver",     "mm",    'F');
  dataloader->AddVariable("oa_lambda",     "oa_lambda",     "mm",    'F');
  //#warning Momentum diabled!

  TFile* input    = TFile::Open(train, "UPDATE");
  TTree* tSigAll  = (TTree*) input->Get("signal");
  TTree* tBackAll = (TTree*) input->Get("background");

  Long64_t MaxEntries = TMath::Min(tSigAll->GetEntries(), tBackAll->GetEntries());
  if (DesEntries > 0 && DesEntries < MaxEntries)
    MaxEntries = DesEntries;

  TTree* tSig  = tSigAll ->CloneTree(MaxEntries);
  TTree* tBack = tBackAll->CloneTree(MaxEntries);

  dataloader->AddSignalTree    (tSig,  1.);
  dataloader->AddBackgroundTree(tBack, 1.);

  dataloader->SetSignalWeightExpression("Weight");
  dataloader->SetBackgroundWeightExpression("Weight");

  dataloader->PrepareTrainingAndTestTree("", "", "!V:SplitMode=Alternate:NormMode=EqualNumEvents");

  // Default
  //factory->BookMethod(dataloader, TMVA::Types::kMLP, "MLP", "H:!V:NCycles=500:HiddenLayers=N,N-1:NeuronType=sigmoid:EstimatorType=MSE:VarTransform=N:TestRate=10:!UseRegulator");
  // Shzymon
  //factory->BookMethod(dataloader, TMVA::Types::kMLP, "MLP", "H:!V:NCycles=500:HiddenLayers=N+1,N:NeuronType=sigmoid:EstimatorType=CE:VarTransform=N:TestRate=5:!UseRegulator");

  // MLP
  factory->BookMethod(dataloader, TMVA::Types::kMLP, "MLP", "!H:!V:NCycles=1500:HiddenLayers=N+2,N:NeuronType=sigmoid:NeuronInputType=sum:EstimatorType=MSE:TrainingMethod=BP:VarTransform=N:BPMode=sequential");

  // Keras - TensorFlow
  //factory->BookMethod(dataloader, TMVA::Types::kPyKeras, "PyKeras", "H:!V:VarTransform=N:FilenameModel=model.h5:NumEpochs=20:BatchSize=32");

  // Comparison
  /*factory->BookMethod( dataloader, TMVA::Types::kPDEFoam, "PDEFoam",
    "!H:!V:SigBgSeparate=F:TailCut=0.001:VolFrac=0.0666:nActiveCells=500:nSampl=2000:nBin=5:Nmin=100:Kernel=None:Compress=T" );
    factory->BookMethod( dataloader, TMVA::Types::kPDEFoam, "PDEFoamBoost",
    "!H:!V:Boost_Num=30:Boost_Transform=linear:SigBgSeparate=F:MaxDepth=4:UseYesNoCell=T:DTLogic=MisClassificationError:FillFoamWithOrigWeights=F:TailCut=0:nActiveCells=500:nBin=20:Nmin=400:Kernel=None:Compress=T" );

    factory->BookMethod( dataloader, TMVA::Types::kMLP, "MLP", "H:!V:NeuronType=tanh:VarTransform=N:NCycles=600:HiddenLayers=N+5:TestRate=5:!UseRegulator" );
    factory->BookMethod( dataloader, TMVA::Types::kMLP, "MLPBFGS", "H:!V:NeuronType=tanh:VarTransform=N:NCycles=600:HiddenLayers=N+5:TestRate=5:TrainingMethod=BFGS:!UseRegulator" );
    factory->BookMethod( dataloader, TMVA::Types::kMLP, "MLPBNN", "H:!V:NeuronType=tanh:VarTransform=N:NCycles=60:HiddenLayers=N+5:TestRate=5:TrainingMethod=BFGS:UseRegulator" ); // BFGS training with bayesian regulators

    TString layoutString ("Layout=TANH|128,TANH|128,TANH|128,LINEAR");

    // Training strategies.
    TString training0("LearningRate=1e-1,Momentum=0.9,Repetitions=1,"
          "ConvergenceSteps=20,BatchSize=256,TestRepetitions=10,"
	        "WeightDecay=1e-4,Regularization=L2,"
		      "DropConfig=0.0+0.5+0.5+0.5, Multithreading=True");
    TString training1("LearningRate=1e-2,Momentum=0.9,Repetitions=1,"
          "ConvergenceSteps=20,BatchSize=256,TestRepetitions=10,"
	        "WeightDecay=1e-4,Regularization=L2,"
		      "DropConfig=0.0+0.0+0.0+0.0, Multithreading=True");
    TString training2("LearningRate=1e-3,Momentum=0.0,Repetitions=1,"
          "ConvergenceSteps=20,BatchSize=256,TestRepetitions=10,"
	        "WeightDecay=1e-4,Regularization=L2,"
		      "DropConfig=0.0+0.0+0.0+0.0, Multithreading=True");
    TString trainingStrategyString ("TrainingStrategy=");
    trainingStrategyString += training0 + "|" + training1 + "|" + training2;

    // General Options.
    TString dnnOptions ("!H:V:ErrorStrategy=CROSSENTROPY:VarTransform=N:"
    "WeightInitialization=XAVIERUNIFORM");
    dnnOptions.Append (":"); dnnOptions.Append (layoutString);
    dnnOptions.Append (":"); dnnOptions.Append (trainingStrategyString);

    // Multi-core CPU implementation.
    TString cpuOptions = dnnOptions + ":Architecture=CPU";
    factory->BookMethod(dataloader, TMVA::Types::kDNN, "DNN_CPU", cpuOptions);

    factory->BookMethod( dataloader, TMVA::Types::kCFMlpANN, "CFMlpANN", "!H:!V:NCycles=200:HiddenLayers=N+1,N"  ); // n_cycles:#nodes:#nodes:...
    factory->BookMethod( dataloader, TMVA::Types::kTMlpANN, "TMlpANN", "!H:!V:NCycles=200:HiddenLayers=N+1,N:LearningMethod=BFGS:ValidationFraction=0.3"  ); // n_cycles:#nodes:#nodes:...

    factory->BookMethod( dataloader, TMVA::Types::kSVM, "SVM", "Gamma=0.25:Tol=0.001:VarTransform=Norm" );

    factory->BookMethod( dataloader, TMVA::Types::kBDT, "BDTG",
    "!H:!V:NTrees=1000:MinNodeSize=2.5%:BoostType=Grad:Shrinkage=0.10:UseBaggedBoost:BaggedSampleFraction=0.5:nCuts=20:MaxDepth=2" );
    factory->BookMethod( dataloader, TMVA::Types::kBDT, "BDT",
    "!H:!V:NTrees=850:MinNodeSize=2.5%:MaxDepth=3:BoostType=AdaBoost:AdaBoostBeta=0.5:UseBaggedBoost:BaggedSampleFraction=0.5:SeparationType=GiniIndex:nCuts=20" );
    factory->BookMethod( dataloader, TMVA::Types::kBDT, "BDTB",
    "!H:!V:NTrees=400:BoostType=Bagging:SeparationType=GiniIndex:nCuts=20" );
    factory->BookMethod( dataloader, TMVA::Types::kBDT, "BDTD",
    "!H:!V:NTrees=400:MinNodeSize=5%:MaxDepth=3:BoostType=AdaBoost:SeparationType=GiniIndex:nCuts=20:VarTransform=Decorrelate" );
    factory->BookMethod( dataloader, TMVA::Types::kBDT, "BDTF",
    "!H:!V:NTrees=50:MinNodeSize=2.5%:UseFisherCuts:MaxDepth=3:BoostType=AdaBoost:AdaBoostBeta=0.5:SeparationType=GiniIndex:nCuts=20" );

    factory->BookMethod( dataloader, TMVA::Types::kRuleFit, "RuleFit",
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
