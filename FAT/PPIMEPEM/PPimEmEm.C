#include "PPimEmEm.h"
//#include "PPimEmEm_buffer.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include "data.h"
#include <iostream>
#include <TLorentzVector.h>
#include <TVector3.h>
#include "hntuple.h"


using namespace std;
using namespace PATData;

void PPimEmEm::Loop()
{
  int licznik = 0;
  if (fChain == 0) return;

  Long64_t nentries = fChain->GetEntries();

  double nbytes = 0, nb = 0;
  int event_number=-1;
  int event_mult=-1;

  int dif_events=1;
  int real_lambdas=0;
    
  for(Long64_t jentry=0; jentry<nentries;jentry++)
    {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0)
	{
	  cout<<"jentry out of range!!!"<<endl;
	  break;
	}
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      cout.precision(10);
      if(jentry%1000==0 && isBest!=-1)
	{
	  cout << "netry no. "<< jentry<<" from "<<nentries;
	  //cout<<" isBest: "<< isBest<<" event: "<<event;
	  cout<<endl;
	}

      if(isBest!=-1 /*&& trigdownscaleflag==1*/)
	{
	  TVector3 eVert (eVert_x,eVert_y,eVert_z);
	  
	  dif_events++;
	  double min_value=10000000;
	  
	  double F = 1.006;
	  //double F=1;
	  TVector3 v1, v2, v3, v4, v5;
	  v2.SetXYZ(F*p_p*sin(D2R*p_theta)*cos(D2R*p_phi),F*p_p*sin(D2R*p_theta)*sin(D2R*p_phi),F*p_p*cos(D2R*p_theta));
	  v3.SetXYZ(F*pim_p*sin(D2R*pim_theta)*cos(D2R*pim_phi),F*pim_p*sin(D2R*pim_theta)*sin(D2R*pim_phi),F*pim_p*cos(D2R*pim_theta));
	  v4.SetXYZ(F*em1_p*sin(D2R*em1_theta)*cos(D2R*em1_phi),F*em1_p*sin(D2R*em1_theta)*sin(D2R*em1_phi),F*em1_p*cos(D2R*em1_theta));
	  v5.SetXYZ(F*em2_p*sin(D2R*em2_theta)*cos(D2R*em2_phi),F*em2_p*sin(D2R*em2_theta)*sin(D2R*em2_phi),F*em2_p*cos(D2R*em2_theta));

	  p->SetVectM( v2, 938.272013 );
	  pim->SetVectM( v3, 139.57018 );
	  em1->SetVectM( v4, 0.51099894 );
	  em2->SetVectM( v5, 0.51099894 );

	  *gammapem1 = *p + *em1;
	  *gammappim = *p + *pim;
	  *gammapem2 = *p + *em2;
	  *gammapimem1= *pim + *em1;
	  *gammaem2em1= *em2 + *em1;
	  *gammappimem1em2=*pim +*em2 + *em1 + *p;
	  *miss=*beam-*gammappimem1em2;

	  
	  Float_t m_inv_ppim = gammappim->M();
	  Float_t m_inv_pem2 = gammapem2->M();
	  Float_t m_inv_em1pim = gammapimem1->M();
	  Float_t m_inv_em1em2 = gammaem2em1->M();
	  Float_t m_inv_ppimem1em2 = gammappimem1em2->M();
	  Float_t oa = R2D * openingangle(*p, *pim);
	  //Float_t oa_rich = R2D * openingangle(r1, r2);

	  Float_t p_mass = p_p*p_p * (  1. / (p_beta*p_beta)  - 1. ) ;
	  Float_t pi_mass = pim_p*pim_p * (  1. / (pim_beta*pim_beta_new)  - 1. ) ;
	  Float_t em1_mass = em1_p*em1_p * (  1. / (em1_beta*em1_beta_new)  - 1. ) ;
	  Float_t pim_mass = pim_p*pim_p * (  1. / (pim_beta*pim_beta_new)  - 1. ) ;
	  Float_t em2_mass = em2_p*em2_p * (  1. / (em2_beta*em2_beta_new)  - 1. ) ;

	  TVector3 ver_p_pim=vertex(p_r,p_z,*p,pim_r,pim_z,*pim);
	  TVector3 ver_p_em2=vertex(p_r,p_z,*p,em2_r,em2_z,*em2);
	  TVector3 ver_em1_pim=vertex(em1_r,em1_z,*em1,pim_r,pim_z,*pim);
	  TVector3 ver_em1_em2=vertex(em1_r,em1_z,*em1,em2_r,em2_z,*em2);

	  TVector3 ver_to_ver=ver_p_pim-ver_em1_em2;
	  //TVector3 ver_to_ver_2=ver_p_em2-ver_em1_pim;

	  Float_t oa_pim_p=R2D*openingangle(pim->Vect(),p->Vect());
	  Float_t oa_em2_p=R2D*openingangle(em2->Vect(),p->Vect());
	  Float_t oa_em1_p=R2D*openingangle(em1->Vect(),p->Vect());
	  Float_t oa_pim_em2=R2D*openingangle(pim->Vect(),em2->Vect());
	  Float_t oa_pim_em1=R2D*openingangle(pim->Vect(),em1->Vect());
	  Float_t oa_em2_em1=R2D*openingangle(em2->Vect(),em1->Vect());
                  
	  Float_t oa_lambda=R2D*openingangle(ver_to_ver,gammappim->Vect());
	  //Float_t oa_lambda_2=R2D*openingangle(ver_to_ver_2,gammapem2->Vect());
      
	  Float_t dist_p_pim=trackDistance(p_r,p_z,*p,pim_r,pim_z,*pim);
	  Float_t dist_p_em2=trackDistance(p_r,p_z,*p,em2_r,em2_z,*em2);
	  Float_t dist_em1_pim=trackDistance(em1_r,em1_z,*em1,pim_r,pim_z,*pim);
	  Float_t dist_em1_em2=trackDistance(em1_r,em1_z,*em1,em2_r,em2_z,*em2);
	  Float_t dist_lambda_em1=trackDistance(em1_r,em1_z,*em1,getR(ver_p_pim),ver_p_pim.Z(),*gammappim);
	  //Float_t dist_lambda2_em1=trackDistance(em1_r,em1_z,*em1,ver_p_em2.Z(),getR(ver_p_em2),*gammapem2);
	  Float_t dist_lambda_em2=trackDistance(em2_r,em2_z,*em2,getR(ver_p_pim),ver_p_pim.Z(),*gammappim);
	  //Float_t dist_lambda2_pim=trackDistance(pim_r,pim_z,*pim,ver_p_em2.Z(),getR(ver_p_em2),*gammapem2);
	  Float_t dist_ver_to_ver=ver_to_ver.Mag();
	  //Float_t dist_ver_to_ver_2=ver_to_ver_2.Mag();
      

	  //cout<<"p pim dist from main part:"<<dist_p_pim<<endl;
	  Float_t dist_lambda_eVert=trackToPoint(ver_p_pim,gammappim->Vect(),eVert);
	  //Float_t dist_lambda2_eVert=trackToPoint(ver_p_em2,gammapem2->Vect(),eVert);;
	  Float_t dist_lambda_ver_em1_em2=trackToPoint(ver_p_pim,gammappim->Vect(),ver_em1_em2);;
	  //Float_t dist_lambda2_ver_em1_pim=trackToPoint(ver_p_em2,gammapem2->Vect(),ver_em1_pim);;

	  Float_t dist_p_eVert=trackToPoint(ver_p_pim,p->Vect(),eVert);
	  //Float_t dist_p2_eVert=trackToPoint(ver_p_em2,p->Vect(),eVert);
	  Float_t dist_pim_eVert=trackToPoint(ver_p_pim,pim->Vect(),eVert);
	  //Float_t dist_em2_eVert=trackToPoint(ver_p_em2,em2->Vect(),eVert);
  	    
	  Float_t lambda_mom_z;
	  TLorentzVector lorentz_lambda1115;
	  TLorentzVector lorentz_k0;

	  lorentz_lambda1115=*gammappim;
	  lorentz_k0=*gammaem2em1;
	  dist_lambda_eVert=dist_lambda_eVert;
	  lambda_mom_z=gammappim->Z();
	  
	  bool simon_cut=(oa_pim_p > 15
			  && dist_p_pim < 10
			  && dist_lambda_eVert > 50
			  && dist_p_eVert > 5
			  && dist_pim_eVert > 15
			  );

	  Float_t lambda_pt=lorentz_lambda1115.Pt();
	  Float_t lambda_w=lorentz_lambda1115.Rapidity();
	  Float_t k0_pt=lorentz_k0.Pt();
	  Float_t k0_w=lorentz_k0.Rapidity();
  
	  //save all important variables
	  (*tlo_emem)["isBest"]=isBest;
	  //(*tlo_emem)["isBest_new"]=isBest_new;
	  (*tlo_emem)["event"]=event;
	  (*tlo_emem)["hneg_mult"]=hneg_mult;
	  (*tlo_emem)["hpos_mult"]=hpos_mult;
	  (*tlo_emem)["eVert_x"]=eVert_x;
	  (*tlo_emem)["eVert_y"]=eVert_y;
	  (*tlo_emem)["eVert_z"]=eVert_z;
	  (*tlo_emem)["totalmult"]=totalmult;
	  (*tlo_emem)["trigdownscaleflag"]=trigdownscaleflag;
	  (*tlo_emem)["trigdownscale"]=trigdownscale;
	  (*tlo_emem)["trigbit"]=trigbit;
	  (*tlo_emem)["trigdec"]=trigdec;
	  (*tlo_emem)["event_mult"]=event_mult;
	  //(*tlo_emem)["hypothesis"]=pim_no;
	  //(*tlo_emem)["hypothesis_quality"]=quality;
  
	  (*tlo_emem)["p_p"]=p_p;
	  (*tlo_emem)["p_theta"] = p_theta;
	  (*tlo_emem)["p_phi"] = p_phi;
	  (*tlo_emem)["p_beta"] = p_beta_new;
	  (*tlo_emem)["p_m"] = p_mass;
	  (*tlo_emem)["p_dedx"]=p_dedx_mdc;
	  (*tlo_emem)["p_q"]=p_q;
	  
	  //(*tlo_emem)["p_sim_p"]=p_sim_p;
	  //(*tlo_emem)["p_sim_id"]=p_sim_id;
	  //(*tlo_emem)["p_sim_parentid"]=p_sim_parentid;
	  //(*tlo_emem)["p_sim_vertex_x"]=p_sim_vertexx;
	  //(*tlo_emem)["p_sim_vertex_y"]=p_sim_vertexy;
	  //(*tlo_emem)["p_sim_vertex_z"]=p_sim_vertexz;
	  
	  (*tlo_emem)["em1_p"]=em1_p;
	  (*tlo_emem)["em1_theta"] = em1_theta;
	  (*tlo_emem)["em1_phi"] = em1_phi;
	  (*tlo_emem)["em1_beta"] = em1_beta_new;
	  (*tlo_emem)["em1_m"] = em1_mass;
	  (*tlo_emem)["em1_dedx"]=em1_dedx_mdc;
	  (*tlo_emem)["em1_q"]=em1_q;
	  
	  //(*tlo_emem)["em1_sim_p"]=em1_sim_p;
	  //(*tlo_emem)["em1_sim_id"]=em1_sim_id;
	  //(*tlo_emem)["em1_sim_parentid"]=em1_sim_parentid;
	  //(*tlo_emem)["em1_sim_vertex_x"]=em1_sim_vertexx;
	  //(*tlo_emem)["em1_sim_vertex_y"]=em1_sim_vertexy;
	  //(*tlo_emem)["em1_sim_vertex_z"]=em1_sim_vertexz;
	  
	  
	  (*tlo_emem)["pim_p"]=pim_p;
	  (*tlo_emem)["pim_theta"] = pim_theta;
	  (*tlo_emem)["pim_phi"] = pim_phi;
	  (*tlo_emem)["pim_beta"] = pim_beta_new;
	  (*tlo_emem)["pim_m"] = pim_mass;
	  (*tlo_emem)["pim_dedx"]=pim_dedx_mdc;
	  (*tlo_emem)["pim_q"]=pim_q;
	  
	  //(*tlo_emem)["pim_sim_p"]=pim_sim_p;
	  //(*tlo_emem)["pim_sim_id"]=pim_sim_id;
	  //(*tlo_emem)["pim_sim_parentid"]=pim_sim_parentid;
	  //(*tlo_emem)["pim_sim_vertex_x"]=pim_sim_vertexx;
	  //(*tlo_emem)["pim_sim_vertex_y"]=pim_sim_vertexy;
	  //(*tlo_emem)["pim_sim_vertex_z"]=pim_sim_vertexz;
	  
	  
	  (*tlo_emem)["em2_p"]=em2_p;
	  (*tlo_emem)["em2_theta"] = em2_theta;
	  (*tlo_emem)["em2_phi"] = em2_phi;
	  (*tlo_emem)["em2_beta"] = em2_beta_new;
	  (*tlo_emem)["em2_m"] = em2_mass;
	  (*tlo_emem)["em2_dedx"]=em2_dedx_mdc;
	  (*tlo_emem)["em2_q"]=em2_q;
	  
	  //(*tlo_emem)["em2_sim_p"]=em2_sim_p;
	  //(*tlo_emem)["em2_sim_id"]=em2_sim_id;
	  //(*tlo_emem)["em2_sim_parentid"]=em2_sim_parentid;
	  //(*tlo_emem)["em2_sim_vertex_x"]=em2_sim_vertexx;
	  //(*tlo_emem)["em2_sim_vertex_y"]=em2_sim_vertexy;
	  //(*tlo_emem)["em2_sim_vertex_z"]=em2_sim_vertexz;
	  	  
	  //(*tlo_emem)["pim_sim_id"]=pim_sim_id;
	  //(*tlo_emem)["pim_sim_parentid"]=pim_sim_parentid;

	  (*tlo_emem)["dist_em1_pim"]=dist_em1_pim;
	  (*tlo_emem)["dist_em1_em2"] = dist_em1_em2;
	  (*tlo_emem)["dist_em1_pim"] = dist_em1_pim;
	  (*tlo_emem)["dist_p_pim"] = dist_p_pim;
	  (*tlo_emem)["dist_p_em2"] = dist_p_em2;
	  (*tlo_emem)["dist_p_pim"] = dist_p_pim;
	  (*tlo_emem)["dist_lambda_em2"] = dist_lambda_em2;
	  //(*tlo_emem)["dist_lambda1_em1"] = dist_lambda1_em1;
	  //(*tlo_emem)["dist_lambda_pim"] = dist_lambda_pim;
	  (*tlo_emem)["dist_lambda_em1"] = dist_lambda_em1;
	  (*tlo_emem)["dist_ver_to_ver"]=dist_ver_to_ver;
	  //(*tlo_emem)["dist_ver_to_ver_2"]=dist_ver_to_ver_2;
	  //(*tlo_emem)["dist_ver_to_ver"]=dist_ver_to_ver;
	  (*tlo_emem)["dist_lambda_eVert"]=dist_lambda_eVert;
	  (*tlo_emem)["dist_lambda_ver_em1_em2"]=dist_lambda_ver_em1_em2;
	  //(*tlo_emem)["dist_lambda_ver_em1_pim"]=dist_lambda_ver_em1_pim;
	  //(*tlo_emem)["dist_lambda2_eVert"]=dist_lambda2_eVert;
	  //(*tlo_emem)["dist_lambda2_ver_em1_pim"]=dist_lambda2_ver_em1_pim;
	  //(*tlo_emem)["dist_lambda_eVert"]=dist_lambda_eVert;
	  //(*tlo_emem)["dist_lambda_ver_em1_pim"]=dist_lambda_ver_em1_pim;
	  (*tlo_emem)["dist_p_eVert"]=dist_p_eVert;
	  (*tlo_emem)["dist_pim_eVert"]=dist_pim_eVert;
	  //(*tlo_emem)["dist_p2_eVert"]=dist_p2_eVert;
	  //(*tlo_emem)["dist_em2_eVert"]=dist_em2_eVert;
	  //(*tlo_emem)["dist_p_eVert"]=dist_p_eVert;
	  //(*tlo_emem)["dist_pim_eVert"]=dist_pim_eVert;
  

  
	  (*tlo_emem)["m_inv_p_pim"] = m_inv_ppim;
	  (*tlo_emem)["m_inv_p_em2"] = m_inv_pem2;
	  
	  //(*tlo_emem)["m_inv_p_pim_em2"]=m_inv_ppimem2;
	  //(*tlo_emem)["m_inv_p_pim_em1"]=m_inv_ppimem1;

	  (*tlo_emem)["m_inv_em1_pim"] = m_inv_em1pim;
	  (*tlo_emem)["m_inv_em1_em2"] = m_inv_em1em2;
	  //(*tlo_emem)["m_inv_em1_pim"] = m_inv_em1pim;
	  (*tlo_emem)["m_inv_p_pim_em1_em2"] = m_inv_ppimem1em2;
	  //(*tlo_emem)["m_inv_p_em1"] = m_inv_pem1;
  
	  (*tlo_emem)["ver_p_pim_x"]=ver_p_pim.X();
	  (*tlo_emem)["ver_p_pim_y"]=ver_p_pim.Y();
	  (*tlo_emem)["ver_p_pim_z"]=ver_p_pim.Z();

	  (*tlo_emem)["ver_p_em2_x"]=ver_p_em2.X();
	  (*tlo_emem)["ver_p_em2_y"]=ver_p_em2.Y();
	  (*tlo_emem)["ver_p_em2_z"]=ver_p_em2.Z();

	  (*tlo_emem)["ver_em1_pim_x"]=ver_em1_pim.X();
	  (*tlo_emem)["ver_em1_pim_y"]=ver_em1_pim.Y();
	  (*tlo_emem)["ver_em1_pim_z"]=ver_em1_pim.Z();

	  (*tlo_emem)["ver_em1_em2_x"]=ver_em1_em2.X();
	  (*tlo_emem)["ver_em1_em2_y"]=ver_em1_em2.Y();
	  (*tlo_emem)["ver_em1_em2_z"]=ver_em1_em2.Z();

	  (*tlo_emem)["oa_lambda"]=oa_lambda;
	  //(*tlo_emem)["oa_lambda_2"]=oa_lambda_2;
	  //(*tlo_emem)["oa_lambda"]=oa_lambda;
	  (*tlo_emem)["oa_pim_p"]=oa_pim_p;
	  (*tlo_emem)["oa_em2_p"]=oa_em2_p;
	  //(*tlo_emem)["oa_pim_p"]=oa_p_pim;
	  (*tlo_emem)["oa_em1_p"]=oa_em1_p;
	  (*tlo_emem)["oa_pim_em2"]=oa_pim_em2;
	  (*tlo_emem)["oa_pim_em1"]=oa_pim_em1;
	  (*tlo_emem)["oa_em2_em1"]=oa_em2_em1;
	 
	  (*tlo_emem)["lambda_mom_z"]=lambda_mom_z;
	  (*tlo_emem)["simon_cuts"]=simon_cut;
	  (*tlo_emem)["miss_mass_kp"]=miss->M();

	  (*tlo_emem)["lambda_pt"]=lambda_pt;
	  (*tlo_emem)["lambda_w"]=lambda_w;
	  (*tlo_emem)["lambda_p"]=lorentz_lambda1115.P();
	  lorentz_lambda1115.Boost(-1*(beam->Vect()));
	  (*tlo_emem)["lambda_theta"]=lorentz_lambda1115.Theta();
  
	  (*tlo_emem)["k0_pt"]=k0_pt;
	  (*tlo_emem)["k0_w"]=k0_w;
	  (*tlo_emem)["k0_p"]=lorentz_k0.P();
	  lorentz_k0.Boost(-1*(beam->Vect()));
	  (*tlo_emem)["k0_theta"]=lorentz_k0.Theta();
    
	  tlo_emem->fill();	
	  
	}
    }
}
PPimEmEm::PPimEmEm(TTree *tree)
{
  // if parameter tree is not specified (or zero), connect the file
  // used to generate this class and read the Tree.
  /*if (tree == 0) {
    TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("be08280235056_dst_gen1_sep08_hadron_out.root");
    if (!f || !f->IsOpen()) {
    f = new TFile("be08280235056_dst_gen1_sep08_hadron_out.root");
    }
    f->GetObject("PPimEpPim",tree);
    }
    Init(tree);*/
  if (tree == 0)
    {
      TChain * chain = new TChain("PPimEmEm_ID","");
      
      chain->Add("/lustre/nyx/hades/user/knowakow/PP/PAT_1/FILES/ppimepem/lepton00.root/PPimEmEm_ID");
      chain->Add("/lustre/nyx/hades/user/knowakow/PP/PAT_1/FILES/ppimepem/lepton02.root/PPimEmEm_ID");
      chain->Add("/lustre/nyx/hades/user/knowakow/PP/PAT_1/FILES/ppimepem/lepton03.root/PPimEmEm_ID");
      chain->Add("/lustre/nyx/hades/user/knowakow/PP/PAT_1/FILES/ppimepem/lepton04.root/PPimEmEm_ID");
      chain->Add("/lustre/nyx/hades/user/knowakow/PP/PAT_1/FILES/ppimepem/lepton05.root/PPimEmEm_ID");
      chain->Add("/lustre/nyx/hades/user/knowakow/PP/PAT_1/FILES/ppimepem/lepton06.root/PPimEmEm_ID");
      chain->Add("/lustre/nyx/hades/user/knowakow/PP/PAT_1/FILES/ppimepem/lepton07.root/PPimEmEm_ID");
      chain->Add("/lustre/nyx/hades/user/knowakow/PP/PAT_1/FILES/ppimepem/lepton08.root/PPimEmEm_ID");
      chain->Add("/lustre/nyx/hades/user/knowakow/PP/PAT_1/FILES/ppimepem/lepton09.root/PPimEmEm_ID");
      chain->Add("/lustre/nyx/hades/user/knowakow/PP/PAT_1/FILES/ppimepem/lepton10.root/PPimEmEm_ID");
      chain->Add("/lustre/nyx/hades/user/knowakow/PP/PAT_1/FILES/ppimepem/lepton11.root/PPimEmEm_ID");
      //chain->Add("/lustre/nyx/hades/user/knowakow/PP/PAT_1/FILES/ppimepem/lepton12.root/PPimEmEm_ID");
      chain->Add("/lustre/nyx/hades/user/knowakow/PP/PAT_1/FILES/ppimepem/lepton01.root/PPimEmEm_ID");
      
      
      tree = chain;
    }

  Init(tree);
}

PPimEmEm::~PPimEmEm()
{
  if (!fChain) return;
  delete fChain->GetCurrentFile();
}

Int_t PPimEmEm::GetEntry(Long64_t entry)
{
  // Read contents of entry.
  if (!fChain) return 0;
  return fChain->GetEntry(entry);
}
Long64_t PPimEmEm::LoadTree(Long64_t entry)
{
  // Set the environment to read one entry
  if (!fChain) return -5;
  Long64_t centry = fChain->LoadTree(entry);
  if (centry < 0) return centry;
  if (fChain->GetTreeNumber() != fCurrent) {
    fCurrent = fChain->GetTreeNumber();
    Notify();
  }
  return centry;
}

void PPimEmEm::Init(TTree *tree)
{
  // The Init() function is called when the selector needs to initialize
  // a new tree or chain. Typically here the branch addresses and branch
  // pointers of the tree will be set.
  // It is normally not necessary to make changes to the generated
  // code, but the routine can be extended by the user if needed.
  // Init() will be called many times when running on PROOF
  // (once per file to be processed).

  // Set branch addresses and branch pointers
  if (!tree) return;
  fChain = tree;
  fCurrent = -1;
  fChain->SetMakeClass(1);
  
  fChain->SetBranchAddress("eVert_chi2", &eVert_chi2, &b_eVert_chi2);
  fChain->SetBranchAddress("eVert_x", &eVert_x, &b_eVert_x);
  fChain->SetBranchAddress("eVert_y", &eVert_y, &b_eVert_y);
  fChain->SetBranchAddress("eVert_z", &eVert_z, &b_eVert_z);
  fChain->SetBranchAddress("event", &event, &b_event);
  fChain->SetBranchAddress("hneg_mult", &hneg_mult, &b_hneg_mult);
  fChain->SetBranchAddress("hpos_mult", &hpos_mult, &b_hpos_mult);
  fChain->SetBranchAddress("isBest", &isBest, &b_isBest);
  fChain->SetBranchAddress("lneg_mult", &lneg_mult, &b_lneg_mult);
  fChain->SetBranchAddress("lpos_mult", &lpos_mult, &b_lpos_mult);
  fChain->SetBranchAddress("p_beta", &p_beta, &b_p_beta);
  fChain->SetBranchAddress("p_beta_new", &p_beta_new, &b_p_beta_new);
  fChain->SetBranchAddress("p_dedx_in", &p_dedx_in, &b_p_dedx_in);
  fChain->SetBranchAddress("p_dedx_in_sigma", &p_dedx_in_sigma, &b_p_dedx_in_sigma);
  fChain->SetBranchAddress("p_dedx_mdc", &p_dedx_mdc, &b_p_dedx_mdc);
  fChain->SetBranchAddress("p_dedx_mdc_sigma", &p_dedx_mdc_sigma, &b_p_dedx_mdc_sigma);
  fChain->SetBranchAddress("p_dedx_out", &p_dedx_out, &b_p_dedx_out);
  fChain->SetBranchAddress("p_dedx_out_sigma", &p_dedx_out_sigma, &b_p_dedx_out_sigma);
  fChain->SetBranchAddress("p_dedx_tof", &p_dedx_tof, &b_p_dedx_tof);
  fChain->SetBranchAddress("p_id", &p_id, &b_p_id);
  fChain->SetBranchAddress("p_isring", &p_isring, &b_p_isring);
  fChain->SetBranchAddress("p_kIsLepton", &p_kIsLepton, &b_p_kIsLepton);
  fChain->SetBranchAddress("p_kIsUsed", &p_kIsUsed, &b_p_kIsUsed);
  fChain->SetBranchAddress("p_mdcchi2", &p_mdcchi2, &b_p_mdcchi2);
  fChain->SetBranchAddress("p_oa_hadr", &p_oa_hadr, &b_p_oa_hadr);
  fChain->SetBranchAddress("p_oa_lept", &p_oa_lept, &b_p_oa_lept);
  fChain->SetBranchAddress("p_p", &p_p, &b_p_p);
  fChain->SetBranchAddress("p_phi", &p_phi, &b_p_phi);
  fChain->SetBranchAddress("p_q", &p_q, &b_p_q);
  fChain->SetBranchAddress("p_r", &p_r, &b_p_r);
  fChain->SetBranchAddress("p_resolution", &p_resolution, &b_p_resolution);
  fChain->SetBranchAddress("p_rkchi2", &p_rkchi2, &b_p_rkchi2);
  fChain->SetBranchAddress("p_sector", &p_sector, &b_p_sector);
  fChain->SetBranchAddress("p_shw_sum0", &p_shw_sum0, &b_p_shw_sum0);
  fChain->SetBranchAddress("p_shw_sum1", &p_shw_sum1, &b_p_shw_sum1);
  fChain->SetBranchAddress("p_shw_sum2", &p_shw_sum2, &b_p_shw_sum2);
  /*
    fChain->SetBranchAddress("p_sim_corrflag", &p_sim_corrflag, &b_p_sim_corrflag);
    fChain->SetBranchAddress("p_sim_geninfo", &p_sim_geninfo, &b_p_sim_geninfo);
    fChain->SetBranchAddress("p_sim_geninfo1", &p_sim_geninfo1, &b_p_sim_geninfo1);
    fChain->SetBranchAddress("p_sim_geninfo2", &p_sim_geninfo2, &b_p_sim_geninfo2);
    fChain->SetBranchAddress("p_sim_genweight", &p_sim_genweight, &b_p_sim_genweight);
    fChain->SetBranchAddress("p_sim_id", &p_sim_id, &b_p_sim_id);
    fChain->SetBranchAddress("p_sim_iscommon", &p_sim_iscommon, &b_p_sim_iscommon);
    fChain->SetBranchAddress("p_sim_mediumid", &p_sim_mediumid, &b_p_sim_mediumid);
    fChain->SetBranchAddress("p_sim_p", &p_sim_p, &b_p_sim_p);
    fChain->SetBranchAddress("p_sim_parentid", &p_sim_parentid, &b_p_sim_parentid);
    fChain->SetBranchAddress("p_sim_primaryflag", &p_sim_primaryflag, &b_p_sim_primaryflag);
    fChain->SetBranchAddress("p_sim_processid", &p_sim_processid, &b_p_sim_processid);
    fChain->SetBranchAddress("p_sim_px", &p_sim_px, &b_p_sim_px);
    fChain->SetBranchAddress("p_sim_py", &p_sim_py, &b_p_sim_py);
    fChain->SetBranchAddress("p_sim_pz", &p_sim_pz, &b_p_sim_pz);
    fChain->SetBranchAddress("p_sim_vertexx", &p_sim_vertexx, &b_p_sim_vertexx);
    fChain->SetBranchAddress("p_sim_vertexy", &p_sim_vertexy, &b_p_sim_vertexy);
    fChain->SetBranchAddress("p_sim_vertexz", &p_sim_vertexz, &b_p_sim_vertexz);
  */  
  fChain->SetBranchAddress("p_system", &p_system, &b_p_system);
  fChain->SetBranchAddress("p_theta", &p_theta, &b_p_theta);
  fChain->SetBranchAddress("p_tof_exp", &p_tof_exp, &b_p_tof_exp);
  fChain->SetBranchAddress("p_tof_mom", &p_tof_mom, &b_p_tof_mom);
  fChain->SetBranchAddress("p_tof_new", &p_tof_new, &b_p_tof_new);
  fChain->SetBranchAddress("p_tofino_mult", &p_tofino_mult, &b_p_tofino_mult);
  fChain->SetBranchAddress("p_track_length", &p_track_length, &b_p_track_length);
  fChain->SetBranchAddress("p_z", &p_z, &b_p_z);
  fChain->SetBranchAddress("pim_beta", &pim_beta, &b_pim_beta);
  fChain->SetBranchAddress("pim_beta_new", &pim_beta_new, &b_pim_beta_new);
  fChain->SetBranchAddress("pim_dedx_in", &pim_dedx_in, &b_pim_dedx_in);
  fChain->SetBranchAddress("pim_dedx_in_sigma", &pim_dedx_in_sigma, &b_pim_dedx_in_sigma);
  fChain->SetBranchAddress("pim_dedx_mdc", &pim_dedx_mdc, &b_pim_dedx_mdc);
  fChain->SetBranchAddress("pim_dedx_mdc_sigma", &pim_dedx_mdc_sigma, &b_pim_dedx_mdc_sigma);
  fChain->SetBranchAddress("pim_dedx_out", &pim_dedx_out, &b_pim_dedx_out);
  fChain->SetBranchAddress("pim_dedx_out_sigma", &pim_dedx_out_sigma, &b_pim_dedx_out_sigma);
  fChain->SetBranchAddress("pim_dedx_tof", &pim_dedx_tof, &b_pim_dedx_tof);
  fChain->SetBranchAddress("pim_id", &pim_id, &b_pim_id);
  fChain->SetBranchAddress("pim_isring", &pim_isring, &b_pim_isring);
  fChain->SetBranchAddress("pim_kIsLepton", &pim_kIsLepton, &b_pim_kIsLepton);
  fChain->SetBranchAddress("pim_kIsUsed", &pim_kIsUsed, &b_pim_kIsUsed);
  fChain->SetBranchAddress("pim_mdcchi2", &pim_mdcchi2, &b_pim_mdcchi2);
  fChain->SetBranchAddress("pim_oa_hadr", &pim_oa_hadr, &b_pim_oa_hadr);
  fChain->SetBranchAddress("pim_oa_lept", &pim_oa_lept, &b_pim_oa_lept);
  fChain->SetBranchAddress("pim_p", &pim_p, &b_pim_p);
  fChain->SetBranchAddress("pim_phi", &pim_phi, &b_pim_phi);
  fChain->SetBranchAddress("pim_q", &pim_q, &b_pim_q);
  fChain->SetBranchAddress("pim_r", &pim_r, &b_pim_r);
  fChain->SetBranchAddress("pim_resolution", &pim_resolution, &b_pim_resolution);
  fChain->SetBranchAddress("pim_rkchi2", &pim_rkchi2, &b_pim_rkchi2);
  fChain->SetBranchAddress("pim_sector", &pim_sector, &b_pim_sector);
  fChain->SetBranchAddress("pim_shw_sum0", &pim_shw_sum0, &b_pim_shw_sum0);
  fChain->SetBranchAddress("pim_shw_sum1", &pim_shw_sum1, &b_pim_shw_sum1);
  fChain->SetBranchAddress("pim_shw_sum2", &pim_shw_sum2, &b_pim_shw_sum2);
  /*
    fChain->SetBranchAddress("pim_sim_corrflag", &pim_sim_corrflag, &b_pim_sim_corrflag);
    fChain->SetBranchAddress("pim_sim_geninfo", &pim_sim_geninfo, &b_pim_sim_geninfo);
    fChain->SetBranchAddress("pim_sim_geninfo1", &pim_sim_geninfo1, &b_pim_sim_geninfo1);
    fChain->SetBranchAddress("pim_sim_geninfo2", &pim_sim_geninfo2, &b_pim_sim_geninfo2);
    fChain->SetBranchAddress("pim_sim_genweight", &pim_sim_genweight, &b_pim_sim_genweight);
    fChain->SetBranchAddress("pim_sim_id", &pim_sim_id, &b_pim_sim_id);
    fChain->SetBranchAddress("pim_sim_iscommon", &pim_sim_iscommon, &b_pim_sim_iscommon);
    fChain->SetBranchAddress("pim_sim_mediumid", &pim_sim_mediumid, &b_pim_sim_mediumid);
    fChain->SetBranchAddress("pim_sim_p", &pim_sim_p, &b_pim_sim_p);
    fChain->SetBranchAddress("pim_sim_parentid", &pim_sim_parentid, &b_pim_sim_parentid);
    fChain->SetBranchAddress("pim_sim_primaryflag", &pim_sim_primaryflag, &b_pim_sim_primaryflag);
    fChain->SetBranchAddress("pim_sim_processid", &pim_sim_processid, &b_pim_sim_processid);
    fChain->SetBranchAddress("pim_sim_px", &pim_sim_px, &b_pim_sim_px);
    fChain->SetBranchAddress("pim_sim_py", &pim_sim_py, &b_pim_sim_py);
    fChain->SetBranchAddress("pim_sim_pz", &pim_sim_pz, &b_pim_sim_pz);
    fChain->SetBranchAddress("pim_sim_vertexx", &pim_sim_vertexx, &b_pim_sim_vertexx);
    fChain->SetBranchAddress("pim_sim_vertexy", &pim_sim_vertexy, &b_pim_sim_vertexy);
    fChain->SetBranchAddress("pim_sim_vertexz", &pim_sim_vertexz, &b_pim_sim_vertexz);
  */  
  fChain->SetBranchAddress("pim_system", &pim_system, &b_pim_system);
  fChain->SetBranchAddress("pim_theta", &pim_theta, &b_pim_theta);
  fChain->SetBranchAddress("pim_tof_exp", &pim_tof_exp, &b_pim_tof_exp);
  fChain->SetBranchAddress("pim_tof_mom", &pim_tof_mom, &b_pim_tof_mom);
  fChain->SetBranchAddress("pim_tof_new", &pim_tof_new, &b_pim_tof_new);
  fChain->SetBranchAddress("pim_tofino_mult", &pim_tofino_mult, &b_pim_tofino_mult);
  fChain->SetBranchAddress("pim_track_length", &pim_track_length, &b_pim_track_length);
  fChain->SetBranchAddress("pim_z", &pim_z, &b_pim_z);
  fChain->SetBranchAddress("em2_beta", &em2_beta, &b_em2_beta);
  fChain->SetBranchAddress("em2_beta_new", &em2_beta_new, &b_em2_beta_new);
  fChain->SetBranchAddress("em2_dedx_in", &em2_dedx_in, &b_em2_dedx_in);
  fChain->SetBranchAddress("em2_dedx_in_sigma", &em2_dedx_in_sigma, &b_em2_dedx_in_sigma);
  fChain->SetBranchAddress("em2_dedx_mdc", &em2_dedx_mdc, &b_em2_dedx_mdc);
  fChain->SetBranchAddress("em2_dedx_mdc_sigma", &em2_dedx_mdc_sigma, &b_em2_dedx_mdc_sigma);
  fChain->SetBranchAddress("em2_dedx_out", &em2_dedx_out, &b_em2_dedx_out);
  fChain->SetBranchAddress("em2_dedx_out_sigma", &em2_dedx_out_sigma, &b_em2_dedx_out_sigma);
  fChain->SetBranchAddress("em2_dedx_tof", &em2_dedx_tof, &b_em2_dedx_tof);
  fChain->SetBranchAddress("em2_id", &em2_id, &b_em2_id);
  fChain->SetBranchAddress("em2_isring", &em2_isring, &b_em2_isring);
  fChain->SetBranchAddress("em2_kIsLepton", &em2_kIsLepton, &b_em2_kIsLepton);
  fChain->SetBranchAddress("em2_kIsUsed", &em2_kIsUsed, &b_em2_kIsUsed);
  fChain->SetBranchAddress("em2_mdcchi2", &em2_mdcchi2, &b_em2_mdcchi2);
  fChain->SetBranchAddress("em2_oa_hadr", &em2_oa_hadr, &b_em2_oa_hadr);
  fChain->SetBranchAddress("em2_oa_lept", &em2_oa_lept, &b_em2_oa_lept);
  fChain->SetBranchAddress("em2_p", &em2_p, &b_em2_p);
  fChain->SetBranchAddress("em2_phi", &em2_phi, &b_em2_phi);
  fChain->SetBranchAddress("em2_q", &em2_q, &b_em2_q);
  fChain->SetBranchAddress("em2_r", &em2_r, &b_em2_r);
  fChain->SetBranchAddress("em2_resolution", &em2_resolution, &b_em2_resolution);
  fChain->SetBranchAddress("em2_rkchi2", &em2_rkchi2, &b_em2_rkchi2);
  fChain->SetBranchAddress("em2_sector", &em2_sector, &b_em2_sector);
  fChain->SetBranchAddress("em2_shw_sum0", &em2_shw_sum0, &b_em2_shw_sum0);
  fChain->SetBranchAddress("em2_shw_sum1", &em2_shw_sum1, &b_em2_shw_sum1);
  fChain->SetBranchAddress("em2_shw_sum2", &em2_shw_sum2, &b_em2_shw_sum2);
  /*
    fChain->SetBranchAddress("em2_sim_corrflag", &em2_sim_corrflag, &b_em2_sim_corrflag);
    fChain->SetBranchAddress("em2_sim_geninfo", &em2_sim_geninfo, &b_em2_sim_geninfo);
    fChain->SetBranchAddress("em2_sim_geninfo1", &em2_sim_geninfo1, &b_em2_sim_geninfo1);
    fChain->SetBranchAddress("em2_sim_geninfo2", &em2_sim_geninfo2, &b_em2_sim_geninfo2);
    fChain->SetBranchAddress("em2_sim_genweight", &em2_sim_genweight, &b_em2_sim_genweight);
    fChain->SetBranchAddress("em2_sim_id", &em2_sim_id, &b_em2_sim_id);
    fChain->SetBranchAddress("em2_sim_iscommon", &em2_sim_iscommon, &b_em2_sim_iscommon);
    fChain->SetBranchAddress("em2_sim_mediumid", &em2_sim_mediumid, &b_em2_sim_mediumid);
    fChain->SetBranchAddress("em2_sim_p", &em2_sim_p, &b_em2_sim_p);
    fChain->SetBranchAddress("em2_sim_parentid", &em2_sim_parentid, &b_em2_sim_parentid);
    fChain->SetBranchAddress("em2_sim_primaryflag", &em2_sim_primaryflag, &b_em2_sim_primaryflag);
    fChain->SetBranchAddress("em2_sim_processid", &em2_sim_processid, &b_em2_sim_processid);
    fChain->SetBranchAddress("em2_sim_px", &em2_sim_px, &b_em2_sim_px);
    fChain->SetBranchAddress("em2_sim_py", &em2_sim_py, &b_em2_sim_py);
    fChain->SetBranchAddress("em2_sim_pz", &em2_sim_pz, &b_em2_sim_pz);
    fChain->SetBranchAddress("em2_sim_vertexx", &em2_sim_vertexx, &b_em2_sim_vertexx);
    fChain->SetBranchAddress("em2_sim_vertexy", &em2_sim_vertexy, &b_em2_sim_vertexy);
    fChain->SetBranchAddress("em2_sim_vertexz", &em2_sim_vertexz, &b_em2_sim_vertexz);
  */  
  fChain->SetBranchAddress("em2_system", &em2_system, &b_em2_system);
  fChain->SetBranchAddress("em2_theta", &em2_theta, &b_em2_theta);
  fChain->SetBranchAddress("em2_tof_exp", &em2_tof_exp, &b_em2_tof_exp);
  fChain->SetBranchAddress("em2_tof_mom", &em2_tof_mom, &b_em2_tof_mom);
  fChain->SetBranchAddress("em2_tof_new", &em2_tof_new, &b_em2_tof_new);
  fChain->SetBranchAddress("em2_tofino_mult", &em2_tofino_mult, &b_em2_tofino_mult);
  fChain->SetBranchAddress("em2_track_length", &em2_track_length, &b_em2_track_length);
  fChain->SetBranchAddress("em2_z", &em2_z, &b_em2_z);
  fChain->SetBranchAddress("em1_beta", &em1_beta, &b_em1_beta);
  fChain->SetBranchAddress("em1_beta_new", &em1_beta_new, &b_em1_beta_new);
  fChain->SetBranchAddress("em1_dedx_in", &em1_dedx_in, &b_em1_dedx_in);
  fChain->SetBranchAddress("em1_dedx_in_sigma", &em1_dedx_in_sigma, &b_em1_dedx_in_sigma);
  fChain->SetBranchAddress("em1_dedx_mdc", &em1_dedx_mdc, &b_em1_dedx_mdc);
  fChain->SetBranchAddress("em1_dedx_mdc_sigma", &em1_dedx_mdc_sigma, &b_em1_dedx_mdc_sigma);
  fChain->SetBranchAddress("em1_dedx_out", &em1_dedx_out, &b_em1_dedx_out);
  fChain->SetBranchAddress("em1_dedx_out_sigma", &em1_dedx_out_sigma, &b_em1_dedx_out_sigma);
  fChain->SetBranchAddress("em1_dedx_tof", &em1_dedx_tof, &b_em1_dedx_tof);
  fChain->SetBranchAddress("em1_id", &em1_id, &b_em1_id);
  fChain->SetBranchAddress("em1_isring", &em1_isring, &b_em1_isring);
  fChain->SetBranchAddress("em1_kIsLepton", &em1_kIsLepton, &b_em1_kIsLepton);
  fChain->SetBranchAddress("em1_kIsUsed", &em1_kIsUsed, &b_em1_kIsUsed);
  fChain->SetBranchAddress("em1_mdcchi2", &em1_mdcchi2, &b_em1_mdcchi2);
  fChain->SetBranchAddress("em1_oa_hadr", &em1_oa_hadr, &b_em1_oa_hadr);
  fChain->SetBranchAddress("em1_oa_lept", &em1_oa_lept, &b_em1_oa_lept);
  fChain->SetBranchAddress("em1_p", &em1_p, &b_em1_p);
  fChain->SetBranchAddress("em1_phi", &em1_phi, &b_em1_phi);
  fChain->SetBranchAddress("em1_q", &em1_q, &b_em1_q);
  fChain->SetBranchAddress("em1_r", &em1_r, &b_em1_r);
  fChain->SetBranchAddress("em1_resolution", &em1_resolution, &b_em1_resolution);
  fChain->SetBranchAddress("em1_rkchi2", &em1_rkchi2, &b_em1_rkchi2);
  fChain->SetBranchAddress("em1_sector", &em1_sector, &b_em1_sector);
  fChain->SetBranchAddress("em1_shw_sum0", &em1_shw_sum0, &b_em1_shw_sum0);
  fChain->SetBranchAddress("em1_shw_sum1", &em1_shw_sum1, &b_em1_shw_sum1);
  fChain->SetBranchAddress("em1_shw_sum2", &em1_shw_sum2, &b_em1_shw_sum2);
  /*
    fChain->SetBranchAddress("em1_sim_corrflag", &em1_sim_corrflag, &b_em1_sim_corrflag);
    fChain->SetBranchAddress("em1_sim_geninfo", &em1_sim_geninfo, &b_em1_sim_geninfo);
    fChain->SetBranchAddress("em1_sim_geninfo1", &em1_sim_geninfo1, &b_em1_sim_geninfo1);
    fChain->SetBranchAddress("em1_sim_geninfo2", &em1_sim_geninfo2, &b_em1_sim_geninfo2);
    fChain->SetBranchAddress("em1_sim_genweight", &em1_sim_genweight, &b_em1_sim_genweight);
    fChain->SetBranchAddress("em1_sim_id", &em1_sim_id, &b_em1_sim_id);
    fChain->SetBranchAddress("em1_sim_iscommon", &em1_sim_iscommon, &b_em1_sim_iscommon);
    fChain->SetBranchAddress("em1_sim_mediumid", &em1_sim_mediumid, &b_em1_sim_mediumid);
    fChain->SetBranchAddress("em1_sim_p", &em1_sim_p, &b_em1_sim_p);
    fChain->SetBranchAddress("em1_sim_parentid", &em1_sim_parentid, &b_em1_sim_parentid);
    fChain->SetBranchAddress("em1_sim_primaryflag", &em1_sim_primaryflag, &b_em1_sim_primaryflag);
    fChain->SetBranchAddress("em1_sim_processid", &em1_sim_processid, &b_em1_sim_processid);
    fChain->SetBranchAddress("em1_sim_px", &em1_sim_px, &b_em1_sim_px);
    fChain->SetBranchAddress("em1_sim_py", &em1_sim_py, &b_em1_sim_py);
    fChain->SetBranchAddress("em1_sim_pz", &em1_sim_pz, &b_em1_sim_pz);
    fChain->SetBranchAddress("em1_sim_vertexx", &em1_sim_vertexx, &b_em1_sim_vertexx);
    fChain->SetBranchAddress("em1_sim_vertexy", &em1_sim_vertexy, &b_em1_sim_vertexy);
    fChain->SetBranchAddress("em1_sim_vertexz", &em1_sim_vertexz, &b_em1_sim_vertexz);
  */  
  fChain->SetBranchAddress("em1_system", &em1_system, &b_em1_system);
  fChain->SetBranchAddress("em1_theta", &em1_theta, &b_em1_theta);
  fChain->SetBranchAddress("em1_tof_exp", &em1_tof_exp, &b_em1_tof_exp);
  fChain->SetBranchAddress("em1_tof_mom", &em1_tof_mom, &b_em1_tof_mom);
  fChain->SetBranchAddress("em1_tof_new", &em1_tof_new, &b_em1_tof_new);
  fChain->SetBranchAddress("em1_tofino_mult", &em1_tofino_mult, &b_em1_tofino_mult);
  fChain->SetBranchAddress("em1_track_length", &em1_track_length, &b_em1_track_length);
  fChain->SetBranchAddress("em1_z", &em1_z, &b_em1_z);
  fChain->SetBranchAddress("runnumber", &runnumber, &b_runnumber);
  fChain->SetBranchAddress("totalmult", &totalmult, &b_totalmult);
  fChain->SetBranchAddress("trigbit", &trigbit, &b_trigbit);
  fChain->SetBranchAddress("trigdec", &trigdec, &b_trigdec);
  fChain->SetBranchAddress("trigdownscale", &trigdownscale, &b_trigdownscale);
  fChain->SetBranchAddress("trigdownscaleflag", &trigdownscaleflag, &b_trigdownscaleflag);
  Notify();
}

Bool_t PPimEmEm::Notify()
{
  // The Notify() function is called when a new file is opened. This
  // can be either for a new TTree in a TChain or when when a new TTree
  // is started when using PROOF. It is normally not necessary to make changes
  // to the generated code, but the routine can be extended by the
  // user if needed. The return value is currently not used.

  return kTRUE;
}

void PPimEmEm::Show(Long64_t entry)
{
  // Print contents of entry.
  // If entry is not specified, print current entry
  if (!fChain) return;
  fChain->Show(entry);
}
Int_t PPimEmEm::Cut(Long64_t entry)
{
  // This function may be called from Loop.
  // returns  1 if entry is accepted.
  // returns -1 otherwise.
  return 1;
}
