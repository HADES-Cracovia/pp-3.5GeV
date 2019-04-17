
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
void TMVAClassification_data(TString treeFile, TString extraSuffix = "", Long64_t DesEntries = -1) {
  TMVA::Tools::Instance();

  cout << "==> Start TMVAClassification" << endl;

  TString NameSuffix = (treeFile.Contains("_") ? treeFile(treeFile.First('_'), treeFile.Last('.') - treeFile.First('_')) : TString("New")) + extraSuffix;
  TFile* outputFile = TFile::Open("TMVATraining" + NameSuffix + ".root", "RECREATE");

  TMVA::Factory* factory = new TMVA::Factory("TMVAClassification" + NameSuffix, outputFile, "!V:!Silent:Color:DrawProgressBar:Transformations=N:AnalysisType=Classification");
  //TMVA::DataLoader* dataloader = new TMVA::DataLoader("dataset");


  factory->AddVariable("oa_lambda","oa_lambda","deg",'F',0,180);
  factory->AddVariable("oa_pip_p","oa_pip_p","deg",'F',0,180);
  factory->AddVariable("p_p","p_p","MeV",'F');
  factory->AddVariable("pip_p","pip_p","MeV",'F');
  //factory->AddVariable("lambda_mom_z","lambda_mom_z","MeV",'F');
  factory->AddVariable("eVert_x",  "eVert_x",  "mm",    'F',-200,200);
  factory->AddVariable("eVert_y",  "eVert_y",  "mm",    'F',-200,200);
  factory->AddVariable("eVert_z",  "eVert_z",  "mm",    'F',-200,200);
  factory->AddVariable("ver_pip_pim_x","ver_pip_pim_x","mm",'F',-200,200);
  factory->AddVariable("ver_pip_pim_y","ver_pip_pim_y","mm",'F',-200,200);
  factory->AddVariable("ver_pip_pim_z","ver_pip_pim_z","mm",'F',-200,200);
  factory->AddVariable("ver_p_pim_x","ver_p_pim_x","mm",'F',-200,200);
  factory->AddVariable("ver_p_pim_y","ver_p_pim_y","mm",'F',-200,200);
  factory->AddVariable("ver_p_pim_z","ver_p_pim_z","mm",'F',-200,200);
  factory->AddVariable("dist_p_pim", "dist_p_pim", "mm",'F',0,200);
  factory->AddVariable("dist_pip_pim", "dist_pip_pim", "mm",'F',0,200);
  factory->AddVariable("dist_p_eVert", "dist_p_eVert","mm",'F',0,200);
  factory->AddVariable("dist_pim_eVert", "dist_pim_eVert",     "mm",'F',0,200);
  factory->AddVariable("dist_lambda_eVert","dist_lambda_eVert",     "mm",    'F',0,200);
  factory->AddVariable("dist_lambda_ver_pip_pim","dist_lambda_ver_pip_pim","mm",'F',0,200);
  //factory->AddVariable("pip_sim_vertex_x","pip_sim_vertex_x","mm",'F');
  //factory->AddVariable("pip_sim_vertex_y","pip_sim_vertex_y","mm",'F');
  //factory->AddVariable("pip_sim_vertex_z","pip_sim_vertex_z","mm",'F');
  //factory->AddVariable("p_sim_vertex_x","p_sim_vertex_x","mm",'F');
  //factory->AddVariable("p_sim_vertex_y","p_sim_vertex_y","mm",'F');
  //factory->AddVariable("p_sim_vertex_z","p_sim_vertex_z","mm",'F');

  
  //#warning Momentum diabled!

  TFile* input    = TFile::Open(treeFile, "UPDATE");
  TTree* tSigAll  = (TTree*) input->Get("signal");
  TTree* tBackAll = (TTree*) input->Get("background");

  /*
Long64_t MaxEntries = TMath::Min(tSigAll->GetEntries(), tBackAll->GetEntries());
  if (DesEntries > 0 && DesEntries < MaxEntries)
    MaxEntries = DesEntries;
  
  TTree* tSig  = tSigAll ->CloneTree(MaxEntries);
  TTree* tBack = tBackAll->CloneTree(MaxEntries);
  */
  TTree* tSig  = tSigAll ->CloneTree();
  TTree* tBack = tBackAll->CloneTree();

  factory->AddSignalTree    (tSig,  1.);
  factory->AddBackgroundTree(tBack, 1);

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

  factory->BookMethod(TMVA::Types::kMLP, "kMLP_pca_ce_600_n2_no_ev", "!H:!V:NCycles=600:HiddenLayers=N+4,N+4:NeuronType=sigmoid:NeuronInputType=sum:EstimatorType=CE:TrainingMethod=BP:VarTransform=N,P:BPMode=sequential:CalculateErrors=True");

  //factory->PrepareTrainingAndTestTree("", "", "!V:SplitMode=Random:SplitSeed=0:NormMode=EqualNumEvents");

  factory->BookMethod(TMVA::Types::kMLP, "kMLP_pca_ce_600_n4_no_ev", "!H:!V:NCycles=600:HiddenLayers=N+4,N+4,N+4,N+4:NeuronType=sigmoid:NeuronInputType=sum:EstimatorType=CE:TrainingMethod=BP:VarTransform=N,P:BPMode=sequential:UseRegulator=True");

  //factory->BookMethod(TMVA::Types::kMLP, "kMLP_pca_ce_1000_n5_eq_ev", "!H:!V:NCycles=1000:HiddenLayers=N,N,N,N,N:NeuronType=sigmoid:NeuronInputType=sum:EstimatorType=CE:TrainingMethod=BP:VarTransform=N,P:BPMode=sequential:UseRegulator=True");

  //factory->BookMethod(TMVA::Types::kMLP, "kMLP_pca_ce_600_n2_eq_ev", "!H:!V:NCycles=300:HiddenLayers=N,N+1:NeuronType=sigmoid:NeuronInputType=sum:EstimatorType=CE:TrainingMethod=BP:VarTransform=N,P:BPMode=sequential:CalculateErrors=True");


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
  //delete factory;
}
