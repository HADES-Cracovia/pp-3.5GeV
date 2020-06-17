#define TMVAeval_cxx
#include "TMVAeval.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TF1.h>
#include <TH1F.h>
#include "hntuple.h"
#include "TMVA/Reader.h"
#include <TGraph.h>
#include <TMath.h>
#include <iostream>

void TMVAeval::Loop()
{
  //   In a ROOT session, you can do:
  //      Root > .L TMVAeval.C
  //      Root > TMVAeval t
  //      Root > t.GetEntry(12); // Fill t data members with entry number 12
  //      Root > t.Show();       // Show values of entry 12
  //      Root > t.Show(16);     // Read and show values of entry 16
  //      Root > t.Loop();       // Loop on all entries
  //

  //     This is the loop skeleton where:
  //    jentry is the global entry number in the chain
  //    ientry is the entry number in the current Tree
  //  Note that the argument to GetEntry must be:
  //    jentry for TChain::GetEntry
  //    ientry for TTree::GetEntry and TBranch::GetEntry
  //
  //       To read only selected branches, Insert statements like:
  // METHOD1:
  //    fChain->SetBranchStatus("*",0);  // disable all branches
  //    fChain->SetBranchStatus("branchname",1);  // activate branchname
  // METHOD2: replace line
  //    fChain->GetEntry(jentry);       //read all branches
  //by  b_branchname->GetEntry(ientry); //read only this branch
  if (fChain == 0) return;

  TFile* outFileData = new TFile("TMVAoutput.root","recreate");
  if(outFileData!=0)
    std::cout<<"Output file created: "<<"TMVAoutput.root"<<endl;
  
  HNtuple *tlo = new HNtuple("TMVAeval","TMVAeval_after TMVA");
  tlo->setFile( outFileData );

   
  Long64_t nentries = fChain->GetEntries();
  Long64_t nbytes = 0, nb = 0;

  TMVA::Reader *reader = new TMVA::Reader();
  reader->AddVariable("dist_p_pim", &dist_p_pim);
  reader->AddVariable("dist_pip_pim", &dist_ep_em);
  reader->AddVariable("eVert_x",  &eVert_x);
  reader->AddVariable("eVert_y",  &eVert_y);
  reader->AddVariable("eVert_z",  &eVert_z);
  reader->AddVariable("ver_pip_pim_x", &ver_ep_em_x);
  reader->AddVariable("ver_pip_pim_y", &ver_ep_em_y);
  reader->AddVariable("ver_pip_pim_z", &ver_ep_em_z);
  reader->AddVariable("ver_p_pim_x",&ver_p_pim_x);
  reader->AddVariable("ver_p_pim_y",&ver_p_pim_y);
  reader->AddVariable("ver_p_pim_z",&ver_p_pim_z);
  reader->AddVariable("oa_lambda",&oa_lambda);
  reader->AddVariable("dist_p_eVert", &dist_p_eVert);
  reader->AddVariable("dist_pim_eVert", &dist_pim_eVert);
  reader->AddVariable("dist_lambda_eVert",&dist_lambda_eVert);
  reader->AddVariable("dist_lambda_ver_pip_pim",&dist_lambda_ver_ep_em);
  reader->AddVariable("dist_ver_to_ver",&dist_ver_to_ver);

  reader->BookMVA("kMLP","/lustre/hades/user/knowakow/PNB/FAT/TMVA/weights/TMVAClassification_data_driven_kMLP_pca_ce_600_(n6+4)_no_ev.weights.xml");
  

  for (Long64_t jentry=0; jentry<nentries;jentry++)
    {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;
      Double_t mlp_output=reader->EvaluateMVA("kMLP");


      //save all important variables
      (*tlo)["isBest"]=isBest;
      //(*tlo)["isBest_new"]=isBest_new;
      (*tlo)["mlp_output"]=mlp_output;
      (*tlo)["event"]=event;
      (*tlo)["hneg_mult"]=hneg_mult;
      (*tlo)["hpos_mult"]=hpos_mult;
      (*tlo)["eVert_x"]=eVert_x;
      (*tlo)["eVert_y"]=eVert_y;
      (*tlo)["eVert_z"]=eVert_z;
      (*tlo)["totalmult"]=totalmult;
      (*tlo)["trigdownscaleflag"]=trigdownscaleflag;
      (*tlo)["trigdownscale"]=trigdownscale;
      (*tlo)["event_mult"]=event_mult;
      //(*tlo)["hypothesis"]=pim_no;
      //(*tlo)["hypothesis_quality"]=quality;
  
      (*tlo)["p_p"]=p_p;
      (*tlo)["p_theta"] = p_theta;
      (*tlo)["p_phi"] = p_phi;
      (*tlo)["p_beta"] = p_beta;
      (*tlo)["p_m"] = p_m;
      (*tlo)["p_dedx"]=p_dedx;
      (*tlo)["p_q"]=p_q;
	  
      //(*tlo)["p_sim_p"]=p_sim_p;
      //(*tlo)["p_sim_id"]=p_sim_id;
      //(*tlo)["p_sim_parentid"]=p_sim_parentid;
      //(*tlo)["p_sim_vertex_x"]=p_sim_vertexx;
      //(*tlo)["p_sim_vertex_y"]=p_sim_vertexy;
      //(*tlo)["p_sim_vertex_z"]=p_sim_vertexz;
	  
      (*tlo)["ep_p"]=ep_p;
      (*tlo)["ep_theta"] = ep_theta;
      (*tlo)["ep_phi"] = ep_phi;
      (*tlo)["ep_beta"] = ep_beta;
      (*tlo)["ep_m"] = ep_m;
      (*tlo)["ep_dedx"]=ep_dedx;
      (*tlo)["ep_q"]=ep_q;
	  
      //(*tlo)["ep_sim_p"]=ep_sim_p;
      //(*tlo)["ep_sim_id"]=ep_sim_id;
      //(*tlo)["ep_sim_parentid"]=ep_sim_parentid;
      //(*tlo)["ep_sim_vertex_x"]=ep_sim_vertexx;
      //(*tlo)["ep_sim_vertex_y"]=ep_sim_vertexy;
      //(*tlo)["ep_sim_vertex_z"]=ep_sim_vertexz;
	  
	  
      (*tlo)["pim_p"]=pim_p;
      (*tlo)["pim_theta"] = pim_theta;
      (*tlo)["pim_phi"] = pim_phi;
      (*tlo)["pim_beta"] = pim_beta;
      (*tlo)["pim_m"] = pim_m;
      (*tlo)["pim_dedx"]=pim_dedx;
      (*tlo)["pim_q"]=pim_q;
	  
      //(*tlo)["pim_sim_p"]=pim_sim_p;
      //(*tlo)["pim_sim_id"]=pim_sim_id;
      //(*tlo)["pim_sim_parentid"]=pim_sim_parentid;
      //(*tlo)["pim_sim_vertex_x"]=pim_sim_vertexx;
      //(*tlo)["pim_sim_vertex_y"]=pim_sim_vertexy;
      //(*tlo)["pim_sim_vertex_z"]=pim_sim_vertexz;
	  
	  
      (*tlo)["em_p"]=em_p;
      (*tlo)["em_theta"] = em_theta;
      (*tlo)["em_phi"] = em_phi;
      (*tlo)["em_beta"] = em_beta;
      (*tlo)["em_m"] = em_m;
      (*tlo)["em_dedx"]=em_dedx;
      (*tlo)["em_q"]=em_q;
	  
      //(*tlo)["em_sim_p"]=em_sim_p;
      //(*tlo)["em_sim_id"]=em_sim_id;
      //(*tlo)["em_sim_parentid"]=em_sim_parentid;
      //(*tlo)["em_sim_vertex_x"]=em_sim_vertexx;
      //(*tlo)["em_sim_vertex_y"]=em_sim_vertexy;
      //(*tlo)["em_sim_vertex_z"]=em_sim_vertexz;
	  	  
      //(*tlo)["pim_sim_id"]=pim_sim_id;
      //(*tlo)["pim_sim_parentid"]=pim_sim_parentid;

      (*tlo)["dist_ep_pim"]=dist_ep_pim;
      (*tlo)["dist_ep_em"] = dist_ep_em;
      (*tlo)["dist_ep_pim"] = dist_ep_pim;
      (*tlo)["dist_p_pim"] = dist_p_pim;
      (*tlo)["dist_p_em"] = dist_p_em;
      (*tlo)["dist_p_pim"] = dist_p_pim;
      (*tlo)["dist_lambda_em"] = dist_lambda_em;
      //(*tlo)["dist_lambda1_ep"] = dist_lambda1_ep;
      //(*tlo)["dist_lambda_pim"] = dist_lambda_pim;
      (*tlo)["dist_lambda_ep"] = dist_lambda_ep;
      (*tlo)["dist_ver_to_ver"]=dist_ver_to_ver;
      //(*tlo)["dist_ver_to_ver_2"]=dist_ver_to_ver_2;
      //(*tlo)["dist_ver_to_ver"]=dist_ver_to_ver;
      (*tlo)["dist_lambda_eVert"]=dist_lambda_eVert;
      (*tlo)["dist_lambda_ver_ep_em"]=dist_lambda_ver_ep_em;
      //(*tlo)["dist_lambda2_eVert"]=dist_lambda2_eVert;
      //(*tlo)["dist_lambda2_ver_ep_pim"]=dist_lambda2_ver_ep_pim;
      //(*tlo)["dist_lambda_eVert"]=dist_lambda_eVert;
      //(*tlo)["dist_lambda_ver_ep_pim"]=dist_lambda_ver_ep_pim;
      (*tlo)["dist_p_eVert"]=dist_p_eVert;
      (*tlo)["dist_pim_eVert"]=dist_pim_eVert;
      //(*tlo)["dist_p2_eVert"]=dist_p2_eVert;
      //(*tlo)["dist_em_eVert"]=dist_em_eVert;
      //(*tlo)["dist_p_eVert"]=dist_p_eVert;
      //(*tlo)["dist_pim_eVert"]=dist_pim_eVert;
  

  
      (*tlo)["m_inv_p_pim"] = m_inv_p_pim;
      (*tlo)["m_inv_p_em"] = m_inv_p_em;
	  
      //(*tlo)["m_inv_p_pim_em"]=m_inv_ppimem;
      //(*tlo)["m_inv_p_pim_ep"]=m_inv_ppimep;

      (*tlo)["m_inv_ep_pim"] = m_inv_ep_pim;
      (*tlo)["m_inv_ep_em"] = m_inv_ep_em;
      //(*tlo)["m_inv_ep_pim"] = m_inv_eppim;
      (*tlo)["m_inv_p_pim_ep_em"] = m_inv_p_pim_ep_em;
      //(*tlo)["m_inv_p_ep"] = m_inv_pep;
  
      (*tlo)["ver_p_pim_x"]=ver_p_pim_x;
      (*tlo)["ver_p_pim_y"]=ver_p_pim_y;
      (*tlo)["ver_p_pim_z"]=ver_p_pim_z;

      (*tlo)["ver_p_em_x"]=ver_p_em_x;
      (*tlo)["ver_p_em_y"]=ver_p_em_y;
      (*tlo)["ver_p_em_z"]=ver_p_em_z;

      (*tlo)["ver_ep_pim_x"]=ver_ep_pim_x;
      (*tlo)["ver_ep_pim_y"]=ver_ep_pim_y;
      (*tlo)["ver_ep_pim_z"]=ver_ep_pim_z;

      (*tlo)["ver_ep_em_x"]=ver_ep_em_x;
      (*tlo)["ver_ep_em_y"]=ver_ep_em_y;
      (*tlo)["ver_ep_em_z"]=ver_ep_em_z;

      (*tlo)["oa_lambda"]=oa_lambda;
      //(*tlo)["oa_lambda_2"]=oa_lambda_2;
      //(*tlo)["oa_lambda"]=oa_lambda;
      (*tlo)["oa_pim_p"]=oa_pim_p;
      (*tlo)["oa_em_p"]=oa_em_p;
      //(*tlo)["oa_pim_p"]=oa_p_pim;
      (*tlo)["oa_ep_p"]=oa_ep_p;
      (*tlo)["oa_pim_em"]=oa_pim_em;
      (*tlo)["oa_pim_ep"]=oa_pim_ep;
      (*tlo)["oa_em_ep"]=oa_em_ep;
	 
      (*tlo)["lambda_mom_z"]=lambda_mom_z;
      (*tlo)["simon_cuts"]=simon_cuts;
      (*tlo)["miss_mass_kp"]=miss_mass_kp;

      (*tlo)["lambda_pt"]=lambda_pt;
      (*tlo)["lambda_w"]=lambda_w;
      (*tlo)["lambda_p"]=lambda_p;
      (*tlo)["lambda_theta"]=lambda_theta;
  
      (*tlo)["k0_pt"]=k0_pt;
      (*tlo)["k0_w"]=k0_w;
      (*tlo)["k0_p"]=k0_p;
      (*tlo)["k0_theta"]=k0_theta;
    
      tlo->fill();
    }
  cout<<"Writing the files"<<endl;
  
  tlo->Write();
  outFileData->Close();
}
