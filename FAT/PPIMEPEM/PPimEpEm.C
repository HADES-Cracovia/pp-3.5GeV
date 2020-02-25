#include "PPimEpEm.h"
//#include "PPimEpEm_buffer.h"
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

void PPimEpEm::Loop()
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
	  v4.SetXYZ(F*ep_p*sin(D2R*ep_theta)*cos(D2R*ep_phi),F*ep_p*sin(D2R*ep_theta)*sin(D2R*ep_phi),F*ep_p*cos(D2R*ep_theta));
	  v5.SetXYZ(F*em_p*sin(D2R*em_theta)*cos(D2R*em_phi),F*em_p*sin(D2R*em_theta)*sin(D2R*em_phi),F*em_p*cos(D2R*em_theta));

	  p->SetVectM( v2, 938.272013 );
	  pim->SetVectM( v3, 139.57018 );
	  ep->SetVectM( v4, 0.51099894 );
	  em->SetVectM( v5, 0.51099894 );

	  *gammapep = *p + *ep;
	  *gammappim = *p + *pim;
	  *gammapem = *p + *em;
	  *gammapimep= *pim + *ep;
	  *gammaemep= *em + *ep;
	  *gammappimepem=*pim +*em + *ep + *p;
	  *miss=*beam-*gammappimepem;

	  
	  Float_t m_inv_ppim = gammappim->M();
	  Float_t m_inv_pem = gammapem->M();
	  Float_t m_inv_eppim = gammapimep->M();
	  Float_t m_inv_epem = gammaemep->M();
	  Float_t m_inv_ppimepem = gammappimepem->M();
	  Float_t oa = R2D * openingangle(*p, *pim);
	  //Float_t oa_rich = R2D * openingangle(r1, r2);

	  Float_t p_mass = p_p*p_p * (  1. / (p_beta*p_beta)  - 1. ) ;
	  Float_t pi_mass = pim_p*pim_p * (  1. / (pim_beta*pim_beta_new)  - 1. ) ;
	  Float_t ep_mass = ep_p*ep_p * (  1. / (ep_beta*ep_beta_new)  - 1. ) ;
	  Float_t pim_mass = pim_p*pim_p * (  1. / (pim_beta*pim_beta_new)  - 1. ) ;
	  Float_t em_mass = em_p*em_p * (  1. / (em_beta*em_beta_new)  - 1. ) ;

	  TVector3 ver_p_pim=vertex(p_r,p_z,*p,pim_r,pim_z,*pim);
	  TVector3 ver_p_em=vertex(p_r,p_z,*p,em_r,em_z,*em);
	  TVector3 ver_ep_pim=vertex(ep_r,ep_z,*ep,pim_r,pim_z,*pim);
	  TVector3 ver_ep_em=vertex(ep_r,ep_z,*ep,em_r,em_z,*em);

	  TVector3 ver_to_ver=ver_p_pim-ver_ep_em;
	  //TVector3 ver_to_ver_2=ver_p_em-ver_ep_pim;

	  Float_t oa_pim_p=R2D*openingangle(pim->Vect(),p->Vect());
	  Float_t oa_em_p=R2D*openingangle(em->Vect(),p->Vect());
	  Float_t oa_ep_p=R2D*openingangle(ep->Vect(),p->Vect());
	  Float_t oa_pim_em=R2D*openingangle(pim->Vect(),em->Vect());
	  Float_t oa_pim_ep=R2D*openingangle(pim->Vect(),ep->Vect());
	  Float_t oa_em_ep=R2D*openingangle(em->Vect(),ep->Vect());
                  
	  Float_t oa_lambda=R2D*openingangle(ver_to_ver,gammappim->Vect());
	  //Float_t oa_lambda_2=R2D*openingangle(ver_to_ver_2,gammapem->Vect());
      
	  Float_t dist_p_pim=trackDistance(p_r,p_z,*p,pim_r,pim_z,*pim);
	  Float_t dist_p_em=trackDistance(p_r,p_z,*p,em_r,em_z,*em);
	  Float_t dist_ep_pim=trackDistance(ep_r,ep_z,*ep,pim_r,pim_z,*pim);
	  Float_t dist_ep_em=trackDistance(ep_r,ep_z,*ep,em_r,em_z,*em);
	  Float_t dist_lambda_ep=trackDistance(ep_r,ep_z,*ep,getR(ver_p_pim),ver_p_pim.Z(),*gammappim);
	  //Float_t dist_lambda2_ep=trackDistance(ep_r,ep_z,*ep,ver_p_em.Z(),getR(ver_p_em),*gammapem);
	  Float_t dist_lambda_em=trackDistance(em_r,em_z,*em,getR(ver_p_pim),ver_p_pim.Z(),*gammappim);
	  //Float_t dist_lambda2_pim=trackDistance(pim_r,pim_z,*pim,ver_p_em.Z(),getR(ver_p_em),*gammapem);
	  Float_t dist_ver_to_ver=ver_to_ver.Mag();
	  //Float_t dist_ver_to_ver_2=ver_to_ver_2.Mag();
      

	  //cout<<"p pim dist from main part:"<<dist_p_pim<<endl;
	  Float_t dist_lambda_eVert=trackToPoint(ver_p_pim,gammappim->Vect(),eVert);
	  //Float_t dist_lambda2_eVert=trackToPoint(ver_p_em,gammapem->Vect(),eVert);;
	  Float_t dist_lambda_ver_ep_em=trackToPoint(ver_p_pim,gammappim->Vect(),ver_ep_em);;
	  //Float_t dist_lambda2_ver_ep_pim=trackToPoint(ver_p_em,gammapem->Vect(),ver_ep_pim);;

	  Float_t dist_p_eVert=trackToPoint(ver_p_pim,p->Vect(),eVert);
	  //Float_t dist_p2_eVert=trackToPoint(ver_p_em,p->Vect(),eVert);
	  Float_t dist_pim_eVert=trackToPoint(ver_p_pim,pim->Vect(),eVert);
	  //Float_t dist_em_eVert=trackToPoint(ver_p_em,em->Vect(),eVert);
  	    
	  Float_t lambda_mom_z;
	  TLorentzVector lorentz_lambda1115;
	  TLorentzVector lorentz_k0;

	  lorentz_lambda1115=*gammappim;
	  lorentz_k0=*gammaemep;
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
	  (*tlo)["isBest"]=isBest;
	  //(*tlo)["isBest_new"]=isBest_new;
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
	  (*tlo)["p_beta"] = p_beta_new;
	  (*tlo)["p_m"] = p_mass;
	  (*tlo)["p_dedx"]=p_dedx_mdc;
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
	  (*tlo)["ep_beta"] = ep_beta_new;
	  (*tlo)["ep_m"] = ep_mass;
	  (*tlo)["ep_dedx"]=ep_dedx_mdc;
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
	  (*tlo)["pim_beta"] = pim_beta_new;
	  (*tlo)["pim_m"] = pim_mass;
	  (*tlo)["pim_dedx"]=pim_dedx_mdc;
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
	  (*tlo)["em_beta"] = em_beta_new;
	  (*tlo)["em_m"] = em_mass;
	  (*tlo)["em_dedx"]=em_dedx_mdc;
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
	  //(*tlo)["dist_lambda_ver_ep_pim"]=dist_lambda_ver_ep_pim;
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
  

  
	  (*tlo)["m_inv_p_pim"] = m_inv_ppim;
	  (*tlo)["m_inv_p_em"] = m_inv_pem;
	  
	  //(*tlo)["m_inv_p_pim_em"]=m_inv_ppimem;
	  //(*tlo)["m_inv_p_pim_ep"]=m_inv_ppimep;

	  (*tlo)["m_inv_ep_pim"] = m_inv_eppim;
	  (*tlo)["m_inv_ep_em"] = m_inv_epem;
	  //(*tlo)["m_inv_ep_pim"] = m_inv_eppim;
	  (*tlo)["m_inv_p_pim_ep_em"] = m_inv_ppimepem;
	  //(*tlo)["m_inv_p_ep"] = m_inv_pep;
  
	  (*tlo)["ver_p_pim_x"]=ver_p_pim.X();
	  (*tlo)["ver_p_pim_y"]=ver_p_pim.Y();
	  (*tlo)["ver_p_pim_z"]=ver_p_pim.Z();

	  (*tlo)["ver_p_em_x"]=ver_p_em.X();
	  (*tlo)["ver_p_em_y"]=ver_p_em.Y();
	  (*tlo)["ver_p_em_z"]=ver_p_em.Z();

	  (*tlo)["ver_ep_pim_x"]=ver_ep_pim.X();
	  (*tlo)["ver_ep_pim_y"]=ver_ep_pim.Y();
	  (*tlo)["ver_ep_pim_z"]=ver_ep_pim.Z();

	  (*tlo)["ver_ep_em_x"]=ver_ep_em.X();
	  (*tlo)["ver_ep_em_y"]=ver_ep_em.Y();
	  (*tlo)["ver_ep_em_z"]=ver_ep_em.Z();

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
	  (*tlo)["simon_cuts"]=simon_cut;
	  (*tlo)["miss_mass_kp"]=miss->M();

	  (*tlo)["lambda_pt"]=lambda_pt;
	  (*tlo)["lambda_w"]=lambda_w;
	  (*tlo)["lambda_p"]=lorentz_lambda1115.P();
	  lorentz_lambda1115.Boost(-1*(beam->Vect()));
	  (*tlo)["lambda_theta"]=lorentz_lambda1115.Theta();
  
	  (*tlo)["k0_pt"]=k0_pt;
	  (*tlo)["k0_w"]=k0_w;
	  (*tlo)["k0_p"]=lorentz_k0.P();
	  lorentz_k0.Boost(-1*(beam->Vect()));
	  (*tlo)["k0_theta"]=lorentz_k0.Theta();
    
	  tlo->fill();	
	  
	}
    }
}
PPimEpEm::PPimEpEm(TTree *tree)
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
      TChain * chain = new TChain("PPimEpEm_ID","");
      
      chain->Add("/lustre/nyx/hades/user/knowakow/PP/PAT_1/FILES/ppimepem/lepton00.root/PPimEpEm_ID");
      chain->Add("/lustre/nyx/hades/user/knowakow/PP/PAT_1/FILES/ppimepem/lepton02.root/PPimEpEm_ID");
      chain->Add("/lustre/nyx/hades/user/knowakow/PP/PAT_1/FILES/ppimepem/lepton03.root/PPimEpEm_ID");
      chain->Add("/lustre/nyx/hades/user/knowakow/PP/PAT_1/FILES/ppimepem/lepton04.root/PPimEpEm_ID");
      chain->Add("/lustre/nyx/hades/user/knowakow/PP/PAT_1/FILES/ppimepem/lepton05.root/PPimEpEm_ID");
      chain->Add("/lustre/nyx/hades/user/knowakow/PP/PAT_1/FILES/ppimepem/lepton06.root/PPimEpEm_ID");
      chain->Add("/lustre/nyx/hades/user/knowakow/PP/PAT_1/FILES/ppimepem/lepton07.root/PPimEpEm_ID");
      chain->Add("/lustre/nyx/hades/user/knowakow/PP/PAT_1/FILES/ppimepem/lepton08.root/PPimEpEm_ID");
      chain->Add("/lustre/nyx/hades/user/knowakow/PP/PAT_1/FILES/ppimepem/lepton09.root/PPimEpEm_ID");
      chain->Add("/lustre/nyx/hades/user/knowakow/PP/PAT_1/FILES/ppimepem/lepton10.root/PPimEpEm_ID");
      chain->Add("/lustre/nyx/hades/user/knowakow/PP/PAT_1/FILES/ppimepem/lepton11.root/PPimEpEm_ID");
      chain->Add("/lustre/nyx/hades/user/knowakow/PP/PAT_1/FILES/ppimepem/lepton12.root/PPimEpEm_ID");
      chain->Add("/lustre/nyx/hades/user/knowakow/PP/PAT_1/FILES/ppimepem/lepton01.root/PPimEpEm_ID");
      
      
      tree = chain;
    }

  Init(tree);
}

PPimEpEm::~PPimEpEm()
{
  if (!fChain) return;
  delete fChain->GetCurrentFile();
}

Int_t PPimEpEm::GetEntry(Long64_t entry)
{
  // Read contents of entry.
  if (!fChain) return 0;
  return fChain->GetEntry(entry);
}
Long64_t PPimEpEm::LoadTree(Long64_t entry)
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

void PPimEpEm::Init(TTree *tree)
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
  fChain->SetBranchAddress("em_beta", &em_beta, &b_em_beta);
  fChain->SetBranchAddress("em_beta_new", &em_beta_new, &b_em_beta_new);
  fChain->SetBranchAddress("em_dedx_in", &em_dedx_in, &b_em_dedx_in);
  fChain->SetBranchAddress("em_dedx_in_sigma", &em_dedx_in_sigma, &b_em_dedx_in_sigma);
  fChain->SetBranchAddress("em_dedx_mdc", &em_dedx_mdc, &b_em_dedx_mdc);
  fChain->SetBranchAddress("em_dedx_mdc_sigma", &em_dedx_mdc_sigma, &b_em_dedx_mdc_sigma);
  fChain->SetBranchAddress("em_dedx_out", &em_dedx_out, &b_em_dedx_out);
  fChain->SetBranchAddress("em_dedx_out_sigma", &em_dedx_out_sigma, &b_em_dedx_out_sigma);
  fChain->SetBranchAddress("em_dedx_tof", &em_dedx_tof, &b_em_dedx_tof);
  fChain->SetBranchAddress("em_id", &em_id, &b_em_id);
  fChain->SetBranchAddress("em_isring", &em_isring, &b_em_isring);
  fChain->SetBranchAddress("em_kIsLepton", &em_kIsLepton, &b_em_kIsLepton);
  fChain->SetBranchAddress("em_kIsUsed", &em_kIsUsed, &b_em_kIsUsed);
  fChain->SetBranchAddress("em_mdcchi2", &em_mdcchi2, &b_em_mdcchi2);
  fChain->SetBranchAddress("em_oa_hadr", &em_oa_hadr, &b_em_oa_hadr);
  fChain->SetBranchAddress("em_oa_lept", &em_oa_lept, &b_em_oa_lept);
  fChain->SetBranchAddress("em_p", &em_p, &b_em_p);
  fChain->SetBranchAddress("em_phi", &em_phi, &b_em_phi);
  fChain->SetBranchAddress("em_q", &em_q, &b_em_q);
  fChain->SetBranchAddress("em_r", &em_r, &b_em_r);
  fChain->SetBranchAddress("em_resolution", &em_resolution, &b_em_resolution);
  fChain->SetBranchAddress("em_rkchi2", &em_rkchi2, &b_em_rkchi2);
  fChain->SetBranchAddress("em_sector", &em_sector, &b_em_sector);
  fChain->SetBranchAddress("em_shw_sum0", &em_shw_sum0, &b_em_shw_sum0);
  fChain->SetBranchAddress("em_shw_sum1", &em_shw_sum1, &b_em_shw_sum1);
  fChain->SetBranchAddress("em_shw_sum2", &em_shw_sum2, &b_em_shw_sum2);
  /*
    fChain->SetBranchAddress("em_sim_corrflag", &em_sim_corrflag, &b_em_sim_corrflag);
    fChain->SetBranchAddress("em_sim_geninfo", &em_sim_geninfo, &b_em_sim_geninfo);
    fChain->SetBranchAddress("em_sim_geninfo1", &em_sim_geninfo1, &b_em_sim_geninfo1);
    fChain->SetBranchAddress("em_sim_geninfo2", &em_sim_geninfo2, &b_em_sim_geninfo2);
    fChain->SetBranchAddress("em_sim_genweight", &em_sim_genweight, &b_em_sim_genweight);
    fChain->SetBranchAddress("em_sim_id", &em_sim_id, &b_em_sim_id);
    fChain->SetBranchAddress("em_sim_iscommon", &em_sim_iscommon, &b_em_sim_iscommon);
    fChain->SetBranchAddress("em_sim_mediumid", &em_sim_mediumid, &b_em_sim_mediumid);
    fChain->SetBranchAddress("em_sim_p", &em_sim_p, &b_em_sim_p);
    fChain->SetBranchAddress("em_sim_parentid", &em_sim_parentid, &b_em_sim_parentid);
    fChain->SetBranchAddress("em_sim_primaryflag", &em_sim_primaryflag, &b_em_sim_primaryflag);
    fChain->SetBranchAddress("em_sim_processid", &em_sim_processid, &b_em_sim_processid);
    fChain->SetBranchAddress("em_sim_px", &em_sim_px, &b_em_sim_px);
    fChain->SetBranchAddress("em_sim_py", &em_sim_py, &b_em_sim_py);
    fChain->SetBranchAddress("em_sim_pz", &em_sim_pz, &b_em_sim_pz);
    fChain->SetBranchAddress("em_sim_vertexx", &em_sim_vertexx, &b_em_sim_vertexx);
    fChain->SetBranchAddress("em_sim_vertexy", &em_sim_vertexy, &b_em_sim_vertexy);
    fChain->SetBranchAddress("em_sim_vertexz", &em_sim_vertexz, &b_em_sim_vertexz);
  */  
  fChain->SetBranchAddress("em_system", &em_system, &b_em_system);
  fChain->SetBranchAddress("em_theta", &em_theta, &b_em_theta);
  fChain->SetBranchAddress("em_tof_exp", &em_tof_exp, &b_em_tof_exp);
  fChain->SetBranchAddress("em_tof_mom", &em_tof_mom, &b_em_tof_mom);
  fChain->SetBranchAddress("em_tof_new", &em_tof_new, &b_em_tof_new);
  fChain->SetBranchAddress("em_tofino_mult", &em_tofino_mult, &b_em_tofino_mult);
  fChain->SetBranchAddress("em_track_length", &em_track_length, &b_em_track_length);
  fChain->SetBranchAddress("em_z", &em_z, &b_em_z);
  fChain->SetBranchAddress("ep_beta", &ep_beta, &b_ep_beta);
  fChain->SetBranchAddress("ep_beta_new", &ep_beta_new, &b_ep_beta_new);
  fChain->SetBranchAddress("ep_dedx_in", &ep_dedx_in, &b_ep_dedx_in);
  fChain->SetBranchAddress("ep_dedx_in_sigma", &ep_dedx_in_sigma, &b_ep_dedx_in_sigma);
  fChain->SetBranchAddress("ep_dedx_mdc", &ep_dedx_mdc, &b_ep_dedx_mdc);
  fChain->SetBranchAddress("ep_dedx_mdc_sigma", &ep_dedx_mdc_sigma, &b_ep_dedx_mdc_sigma);
  fChain->SetBranchAddress("ep_dedx_out", &ep_dedx_out, &b_ep_dedx_out);
  fChain->SetBranchAddress("ep_dedx_out_sigma", &ep_dedx_out_sigma, &b_ep_dedx_out_sigma);
  fChain->SetBranchAddress("ep_dedx_tof", &ep_dedx_tof, &b_ep_dedx_tof);
  fChain->SetBranchAddress("ep_id", &ep_id, &b_ep_id);
  fChain->SetBranchAddress("ep_isring", &ep_isring, &b_ep_isring);
  fChain->SetBranchAddress("ep_kIsLepton", &ep_kIsLepton, &b_ep_kIsLepton);
  fChain->SetBranchAddress("ep_kIsUsed", &ep_kIsUsed, &b_ep_kIsUsed);
  fChain->SetBranchAddress("ep_mdcchi2", &ep_mdcchi2, &b_ep_mdcchi2);
  fChain->SetBranchAddress("ep_oa_hadr", &ep_oa_hadr, &b_ep_oa_hadr);
  fChain->SetBranchAddress("ep_oa_lept", &ep_oa_lept, &b_ep_oa_lept);
  fChain->SetBranchAddress("ep_p", &ep_p, &b_ep_p);
  fChain->SetBranchAddress("ep_phi", &ep_phi, &b_ep_phi);
  fChain->SetBranchAddress("ep_q", &ep_q, &b_ep_q);
  fChain->SetBranchAddress("ep_r", &ep_r, &b_ep_r);
  fChain->SetBranchAddress("ep_resolution", &ep_resolution, &b_ep_resolution);
  fChain->SetBranchAddress("ep_rkchi2", &ep_rkchi2, &b_ep_rkchi2);
  fChain->SetBranchAddress("ep_sector", &ep_sector, &b_ep_sector);
  fChain->SetBranchAddress("ep_shw_sum0", &ep_shw_sum0, &b_ep_shw_sum0);
  fChain->SetBranchAddress("ep_shw_sum1", &ep_shw_sum1, &b_ep_shw_sum1);
  fChain->SetBranchAddress("ep_shw_sum2", &ep_shw_sum2, &b_ep_shw_sum2);
  /*
    fChain->SetBranchAddress("ep_sim_corrflag", &ep_sim_corrflag, &b_ep_sim_corrflag);
    fChain->SetBranchAddress("ep_sim_geninfo", &ep_sim_geninfo, &b_ep_sim_geninfo);
    fChain->SetBranchAddress("ep_sim_geninfo1", &ep_sim_geninfo1, &b_ep_sim_geninfo1);
    fChain->SetBranchAddress("ep_sim_geninfo2", &ep_sim_geninfo2, &b_ep_sim_geninfo2);
    fChain->SetBranchAddress("ep_sim_genweight", &ep_sim_genweight, &b_ep_sim_genweight);
    fChain->SetBranchAddress("ep_sim_id", &ep_sim_id, &b_ep_sim_id);
    fChain->SetBranchAddress("ep_sim_iscommon", &ep_sim_iscommon, &b_ep_sim_iscommon);
    fChain->SetBranchAddress("ep_sim_mediumid", &ep_sim_mediumid, &b_ep_sim_mediumid);
    fChain->SetBranchAddress("ep_sim_p", &ep_sim_p, &b_ep_sim_p);
    fChain->SetBranchAddress("ep_sim_parentid", &ep_sim_parentid, &b_ep_sim_parentid);
    fChain->SetBranchAddress("ep_sim_primaryflag", &ep_sim_primaryflag, &b_ep_sim_primaryflag);
    fChain->SetBranchAddress("ep_sim_processid", &ep_sim_processid, &b_ep_sim_processid);
    fChain->SetBranchAddress("ep_sim_px", &ep_sim_px, &b_ep_sim_px);
    fChain->SetBranchAddress("ep_sim_py", &ep_sim_py, &b_ep_sim_py);
    fChain->SetBranchAddress("ep_sim_pz", &ep_sim_pz, &b_ep_sim_pz);
    fChain->SetBranchAddress("ep_sim_vertexx", &ep_sim_vertexx, &b_ep_sim_vertexx);
    fChain->SetBranchAddress("ep_sim_vertexy", &ep_sim_vertexy, &b_ep_sim_vertexy);
    fChain->SetBranchAddress("ep_sim_vertexz", &ep_sim_vertexz, &b_ep_sim_vertexz);
  */  
  fChain->SetBranchAddress("ep_system", &ep_system, &b_ep_system);
  fChain->SetBranchAddress("ep_theta", &ep_theta, &b_ep_theta);
  fChain->SetBranchAddress("ep_tof_exp", &ep_tof_exp, &b_ep_tof_exp);
  fChain->SetBranchAddress("ep_tof_mom", &ep_tof_mom, &b_ep_tof_mom);
  fChain->SetBranchAddress("ep_tof_new", &ep_tof_new, &b_ep_tof_new);
  fChain->SetBranchAddress("ep_tofino_mult", &ep_tofino_mult, &b_ep_tofino_mult);
  fChain->SetBranchAddress("ep_track_length", &ep_track_length, &b_ep_track_length);
  fChain->SetBranchAddress("ep_z", &ep_z, &b_ep_z);
  fChain->SetBranchAddress("runnumber", &runnumber, &b_runnumber);
  fChain->SetBranchAddress("totalmult", &totalmult, &b_totalmult);
  fChain->SetBranchAddress("trigbit", &trigbit, &b_trigbit);
  fChain->SetBranchAddress("trigdec", &trigdec, &b_trigdec);
  fChain->SetBranchAddress("trigdownscale", &trigdownscale, &b_trigdownscale);
  fChain->SetBranchAddress("trigdownscaleflag", &trigdownscaleflag, &b_trigdownscaleflag);
  Notify();
}

Bool_t PPimEpEm::Notify()
{
  // The Notify() function is called when a new file is opened. This
  // can be either for a new TTree in a TChain or when when a new TTree
  // is started when using PROOF. It is normally not necessary to make changes
  // to the generated code, but the routine can be extended by the
  // user if needed. The return value is currently not used.

  return kTRUE;
}

void PPimEpEm::Show(Long64_t entry)
{
  // Print contents of entry.
  // If entry is not specified, print current entry
  if (!fChain) return;
  fChain->Show(entry);
}
Int_t PPimEpEm::Cut(Long64_t entry)
{
  // This function may be called from Loop.
  // returns  1 if entry is accepted.
  // returns -1 otherwise.
  return 1;
}
