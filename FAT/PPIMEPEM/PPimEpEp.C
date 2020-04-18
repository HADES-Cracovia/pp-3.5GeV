#include "PPimEpEp.h"
//#include "PPimEpEp_buffer.h"
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

void PPimEpEp::Loop()
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
	  v4.SetXYZ(F*ep1_p*sin(D2R*ep1_theta)*cos(D2R*ep1_phi),F*ep1_p*sin(D2R*ep1_theta)*sin(D2R*ep1_phi),F*ep1_p*cos(D2R*ep1_theta));
	  v5.SetXYZ(F*ep2_p*sin(D2R*ep2_theta)*cos(D2R*ep2_phi),F*ep2_p*sin(D2R*ep2_theta)*sin(D2R*ep2_phi),F*ep2_p*cos(D2R*ep2_theta));

	  p->SetVectM( v2, 938.272013 );
	  pim->SetVectM( v3, 139.57018 );
	  ep1->SetVectM( v4, 0.51099894 );
	  ep2->SetVectM( v5, 0.51099894 );

	  *gammapep1 = *p + *ep1;
	  *gammappim = *p + *pim;
	  *gammapep2 = *p + *ep2;
	  *gammapimep1= *pim + *ep1;
	  *gammaep2ep1= *ep2 + *ep1;
	  *gammappimep1ep2=*pim +*ep2 + *ep1 + *p;
	  *miss=*beam-*gammappimep1ep2;

	  
	  Float_t m_inv_ppim = gammappim->M();
	  Float_t m_inv_pep2 = gammapep2->M();
	  Float_t m_inv_ep1pim = gammapimep1->M();
	  Float_t m_inv_ep1ep2 = gammaep2ep1->M();
	  Float_t m_inv_ppimep1ep2 = gammappimep1ep2->M();
	  Float_t oa = R2D * openingangle(*p, *pim);
	  //Float_t oa_rich = R2D * openingangle(r1, r2);

	  Float_t p_mass = p_p*p_p * (  1. / (p_beta*p_beta)  - 1. ) ;
	  Float_t pi_mass = pim_p*pim_p * (  1. / (pim_beta*pim_beta_new)  - 1. ) ;
	  Float_t ep1_mass = ep1_p*ep1_p * (  1. / (ep1_beta*ep1_beta_new)  - 1. ) ;
	  Float_t pim_mass = pim_p*pim_p * (  1. / (pim_beta*pim_beta_new)  - 1. ) ;
	  Float_t ep2_mass = ep2_p*ep2_p * (  1. / (ep2_beta*ep2_beta_new)  - 1. ) ;

	  TVector3 ver_p_pim=vertex(p_r,p_z,*p,pim_r,pim_z,*pim);
	  TVector3 ver_p_ep2=vertex(p_r,p_z,*p,ep2_r,ep2_z,*ep2);
	  TVector3 ver_ep1_pim=vertex(ep1_r,ep1_z,*ep1,pim_r,pim_z,*pim);
	  TVector3 ver_ep1_ep2=vertex(ep1_r,ep1_z,*ep1,ep2_r,ep2_z,*ep2);

	  TVector3 ver_to_ver=ver_p_pim-ver_ep1_ep2;
	  //TVector3 ver_to_ver_2=ver_p_ep2-ver_ep1_pim;

	  Float_t oa_pim_p=R2D*openingangle(pim->Vect(),p->Vect());
	  Float_t oa_ep2_p=R2D*openingangle(ep2->Vect(),p->Vect());
	  Float_t oa_ep1_p=R2D*openingangle(ep1->Vect(),p->Vect());
	  Float_t oa_pim_ep2=R2D*openingangle(pim->Vect(),ep2->Vect());
	  Float_t oa_pim_ep1=R2D*openingangle(pim->Vect(),ep1->Vect());
	  Float_t oa_ep2_ep1=R2D*openingangle(ep2->Vect(),ep1->Vect());
                  
	  Float_t oa_lambda=R2D*openingangle(ver_to_ver,gammappim->Vect());
	  //Float_t oa_lambda_2=R2D*openingangle(ver_to_ver_2,gammapep2->Vect());
      
	  Float_t dist_p_pim=trackDistance(p_r,p_z,*p,pim_r,pim_z,*pim);
	  Float_t dist_p_ep2=trackDistance(p_r,p_z,*p,ep2_r,ep2_z,*ep2);
	  Float_t dist_ep1_pim=trackDistance(ep1_r,ep1_z,*ep1,pim_r,pim_z,*pim);
	  Float_t dist_ep1_ep2=trackDistance(ep1_r,ep1_z,*ep1,ep2_r,ep2_z,*ep2);
	  Float_t dist_lambda_ep1=trackDistance(ep1_r,ep1_z,*ep1,getR(ver_p_pim),ver_p_pim.Z(),*gammappim);
	  //Float_t dist_lambda2_ep1=trackDistance(ep1_r,ep1_z,*ep1,ver_p_ep2.Z(),getR(ver_p_ep2),*gammapep2);
	  Float_t dist_lambda_ep2=trackDistance(ep2_r,ep2_z,*ep2,getR(ver_p_pim),ver_p_pim.Z(),*gammappim);
	  //Float_t dist_lambda2_pim=trackDistance(pim_r,pim_z,*pim,ver_p_ep2.Z(),getR(ver_p_ep2),*gammapep2);
	  Float_t dist_ver_to_ver=ver_to_ver.Mag();
	  //Float_t dist_ver_to_ver_2=ver_to_ver_2.Mag();
      

	  //cout<<"p pim dist from main part:"<<dist_p_pim<<endl;
	  Float_t dist_lambda_eVert=trackToPoint(ver_p_pim,gammappim->Vect(),eVert);
	  //Float_t dist_lambda2_eVert=trackToPoint(ver_p_ep2,gammapep2->Vect(),eVert);;
	  Float_t dist_lambda_ver_ep1_ep2=trackToPoint(ver_p_pim,gammappim->Vect(),ver_ep1_ep2);;
	  //Float_t dist_lambda2_ver_ep1_pim=trackToPoint(ver_p_ep2,gammapep2->Vect(),ver_ep1_pim);;

	  Float_t dist_p_eVert=trackToPoint(ver_p_pim,p->Vect(),eVert);
	  //Float_t dist_p2_eVert=trackToPoint(ver_p_ep2,p->Vect(),eVert);
	  Float_t dist_pim_eVert=trackToPoint(ver_p_pim,pim->Vect(),eVert);
	  //Float_t dist_ep2_eVert=trackToPoint(ver_p_ep2,ep2->Vect(),eVert);
  	    
	  Float_t lambda_mom_z;
	  TLorentzVector lorentz_lambda1115;
	  TLorentzVector lorentz_k0;

	  lorentz_lambda1115=*gammappim;
	  lorentz_k0=*gammaep2ep1;
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
	  (*tlo)["trigbit"]=trigbit;
	  (*tlo)["trigdec"]=trigdec;
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
	  
	  (*tlo)["ep1_p"]=ep1_p;
	  (*tlo)["ep1_theta"] = ep1_theta;
	  (*tlo)["ep1_phi"] = ep1_phi;
	  (*tlo)["ep1_beta"] = ep1_beta_new;
	  (*tlo)["ep1_m"] = ep1_mass;
	  (*tlo)["ep1_dedx"]=ep1_dedx_mdc;
	  (*tlo)["ep1_q"]=ep1_q;
	  
	  //(*tlo)["ep1_sim_p"]=ep1_sim_p;
	  //(*tlo)["ep1_sim_id"]=ep1_sim_id;
	  //(*tlo)["ep1_sim_parentid"]=ep1_sim_parentid;
	  //(*tlo)["ep1_sim_vertex_x"]=ep1_sim_vertexx;
	  //(*tlo)["ep1_sim_vertex_y"]=ep1_sim_vertexy;
	  //(*tlo)["ep1_sim_vertex_z"]=ep1_sim_vertexz;
	  
	  
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
	  
	  
	  (*tlo)["ep2_p"]=ep2_p;
	  (*tlo)["ep2_theta"] = ep2_theta;
	  (*tlo)["ep2_phi"] = ep2_phi;
	  (*tlo)["ep2_beta"] = ep2_beta_new;
	  (*tlo)["ep2_m"] = ep2_mass;
	  (*tlo)["ep2_dedx"]=ep2_dedx_mdc;
	  (*tlo)["ep2_q"]=ep2_q;
	  
	  //(*tlo)["ep2_sim_p"]=ep2_sim_p;
	  //(*tlo)["ep2_sim_id"]=ep2_sim_id;
	  //(*tlo)["ep2_sim_parentid"]=ep2_sim_parentid;
	  //(*tlo)["ep2_sim_vertex_x"]=ep2_sim_vertexx;
	  //(*tlo)["ep2_sim_vertex_y"]=ep2_sim_vertexy;
	  //(*tlo)["ep2_sim_vertex_z"]=ep2_sim_vertexz;
	  	  
	  //(*tlo)["pim_sim_id"]=pim_sim_id;
	  //(*tlo)["pim_sim_parentid"]=pim_sim_parentid;

	  (*tlo)["dist_ep1_pim"]=dist_ep1_pim;
	  (*tlo)["dist_ep1_ep2"] = dist_ep1_ep2;
	  (*tlo)["dist_ep1_pim"] = dist_ep1_pim;
	  (*tlo)["dist_p_pim"] = dist_p_pim;
	  (*tlo)["dist_p_ep2"] = dist_p_ep2;
	  (*tlo)["dist_p_pim"] = dist_p_pim;
	  (*tlo)["dist_lambda_ep2"] = dist_lambda_ep2;
	  //(*tlo)["dist_lambda1_ep1"] = dist_lambda1_ep1;
	  //(*tlo)["dist_lambda_pim"] = dist_lambda_pim;
	  (*tlo)["dist_lambda_ep1"] = dist_lambda_ep1;
	  (*tlo)["dist_ver_to_ver"]=dist_ver_to_ver;
	  //(*tlo)["dist_ver_to_ver_2"]=dist_ver_to_ver_2;
	  //(*tlo)["dist_ver_to_ver"]=dist_ver_to_ver;
	  (*tlo)["dist_lambda_eVert"]=dist_lambda_eVert;
	  (*tlo)["dist_lambda_ver_ep1_ep2"]=dist_lambda_ver_ep1_ep2;
	  //(*tlo)["dist_lambda_ver_ep1_pim"]=dist_lambda_ver_ep1_pim;
	  //(*tlo)["dist_lambda2_eVert"]=dist_lambda2_eVert;
	  //(*tlo)["dist_lambda2_ver_ep1_pim"]=dist_lambda2_ver_ep1_pim;
	  //(*tlo)["dist_lambda_eVert"]=dist_lambda_eVert;
	  //(*tlo)["dist_lambda_ver_ep1_pim"]=dist_lambda_ver_ep1_pim;
	  (*tlo)["dist_p_eVert"]=dist_p_eVert;
	  (*tlo)["dist_pim_eVert"]=dist_pim_eVert;
	  //(*tlo)["dist_p2_eVert"]=dist_p2_eVert;
	  //(*tlo)["dist_ep2_eVert"]=dist_ep2_eVert;
	  //(*tlo)["dist_p_eVert"]=dist_p_eVert;
	  //(*tlo)["dist_pim_eVert"]=dist_pim_eVert;
  

  
	  (*tlo)["m_inv_p_pim"] = m_inv_ppim;
	  (*tlo)["m_inv_p_ep2"] = m_inv_pep2;
	  
	  //(*tlo)["m_inv_p_pim_ep2"]=m_inv_ppimep2;
	  //(*tlo)["m_inv_p_pim_ep1"]=m_inv_ppimep1;

	  (*tlo)["m_inv_ep1_pim"] = m_inv_ep1pim;
	  (*tlo)["m_inv_ep1_ep2"] = m_inv_ep1ep2;
	  //(*tlo)["m_inv_ep1_pim"] = m_inv_ep1pim;
	  (*tlo)["m_inv_p_pim_ep1_ep2"] = m_inv_ppimep1ep2;
	  //(*tlo)["m_inv_p_ep1"] = m_inv_pep1;
  
	  (*tlo)["ver_p_pim_x"]=ver_p_pim.X();
	  (*tlo)["ver_p_pim_y"]=ver_p_pim.Y();
	  (*tlo)["ver_p_pim_z"]=ver_p_pim.Z();

	  (*tlo)["ver_p_ep2_x"]=ver_p_ep2.X();
	  (*tlo)["ver_p_ep2_y"]=ver_p_ep2.Y();
	  (*tlo)["ver_p_ep2_z"]=ver_p_ep2.Z();

	  (*tlo)["ver_ep1_pim_x"]=ver_ep1_pim.X();
	  (*tlo)["ver_ep1_pim_y"]=ver_ep1_pim.Y();
	  (*tlo)["ver_ep1_pim_z"]=ver_ep1_pim.Z();

	  (*tlo)["ver_ep1_ep2_x"]=ver_ep1_ep2.X();
	  (*tlo)["ver_ep1_ep2_y"]=ver_ep1_ep2.Y();
	  (*tlo)["ver_ep1_ep2_z"]=ver_ep1_ep2.Z();

	  (*tlo)["oa_lambda"]=oa_lambda;
	  //(*tlo)["oa_lambda_2"]=oa_lambda_2;
	  //(*tlo)["oa_lambda"]=oa_lambda;
	  (*tlo)["oa_pim_p"]=oa_pim_p;
	  (*tlo)["oa_ep2_p"]=oa_ep2_p;
	  //(*tlo)["oa_pim_p"]=oa_p_pim;
	  (*tlo)["oa_ep1_p"]=oa_ep1_p;
	  (*tlo)["oa_pim_ep2"]=oa_pim_ep2;
	  (*tlo)["oa_pim_ep1"]=oa_pim_ep1;
	  (*tlo)["oa_ep2_ep1"]=oa_ep2_ep1;
	 
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
PPimEpEp::PPimEpEp(TTree *tree)
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
      TChain * chain = new TChain("PPimEpEp_ID","");
      
      chain->Add("/lustre/nyx/hades/user/knowakow/PP/PAT_1/FILES/ppimepem/lepton00.root/PPimEpEp_ID");
      chain->Add("/lustre/nyx/hades/user/knowakow/PP/PAT_1/FILES/ppimepem/lepton02.root/PPimEpEp_ID");
      chain->Add("/lustre/nyx/hades/user/knowakow/PP/PAT_1/FILES/ppimepem/lepton03.root/PPimEpEp_ID");
      chain->Add("/lustre/nyx/hades/user/knowakow/PP/PAT_1/FILES/ppimepem/lepton04.root/PPimEpEp_ID");
      chain->Add("/lustre/nyx/hades/user/knowakow/PP/PAT_1/FILES/ppimepem/lepton05.root/PPimEpEp_ID");
      chain->Add("/lustre/nyx/hades/user/knowakow/PP/PAT_1/FILES/ppimepem/lepton06.root/PPimEpEp_ID");
      chain->Add("/lustre/nyx/hades/user/knowakow/PP/PAT_1/FILES/ppimepem/lepton07.root/PPimEpEp_ID");
      chain->Add("/lustre/nyx/hades/user/knowakow/PP/PAT_1/FILES/ppimepem/lepton08.root/PPimEpEp_ID");
      chain->Add("/lustre/nyx/hades/user/knowakow/PP/PAT_1/FILES/ppimepem/lepton09.root/PPimEpEp_ID");
      chain->Add("/lustre/nyx/hades/user/knowakow/PP/PAT_1/FILES/ppimepem/lepton10.root/PPimEpEp_ID");
      chain->Add("/lustre/nyx/hades/user/knowakow/PP/PAT_1/FILES/ppimepem/lepton11.root/PPimEpEp_ID");
      chain->Add("/lustre/nyx/hades/user/knowakow/PP/PAT_1/FILES/ppimepem/lepton12.root/PPimEpEp_ID");
      chain->Add("/lustre/nyx/hades/user/knowakow/PP/PAT_1/FILES/ppimepem/lepton01.root/PPimEpEp_ID");
      
      
      tree = chain;
    }

  Init(tree);
}

PPimEpEp::~PPimEpEp()
{
  if (!fChain) return;
  delete fChain->GetCurrentFile();
}

Int_t PPimEpEp::GetEntry(Long64_t entry)
{
  // Read contents of entry.
  if (!fChain) return 0;
  return fChain->GetEntry(entry);
}
Long64_t PPimEpEp::LoadTree(Long64_t entry)
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

void PPimEpEp::Init(TTree *tree)
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
  fChain->SetBranchAddress("ep2_beta", &ep2_beta, &b_ep2_beta);
  fChain->SetBranchAddress("ep2_beta_new", &ep2_beta_new, &b_ep2_beta_new);
  fChain->SetBranchAddress("ep2_dedx_in", &ep2_dedx_in, &b_ep2_dedx_in);
  fChain->SetBranchAddress("ep2_dedx_in_sigma", &ep2_dedx_in_sigma, &b_ep2_dedx_in_sigma);
  fChain->SetBranchAddress("ep2_dedx_mdc", &ep2_dedx_mdc, &b_ep2_dedx_mdc);
  fChain->SetBranchAddress("ep2_dedx_mdc_sigma", &ep2_dedx_mdc_sigma, &b_ep2_dedx_mdc_sigma);
  fChain->SetBranchAddress("ep2_dedx_out", &ep2_dedx_out, &b_ep2_dedx_out);
  fChain->SetBranchAddress("ep2_dedx_out_sigma", &ep2_dedx_out_sigma, &b_ep2_dedx_out_sigma);
  fChain->SetBranchAddress("ep2_dedx_tof", &ep2_dedx_tof, &b_ep2_dedx_tof);
  fChain->SetBranchAddress("ep2_id", &ep2_id, &b_ep2_id);
  fChain->SetBranchAddress("ep2_isring", &ep2_isring, &b_ep2_isring);
  fChain->SetBranchAddress("ep2_kIsLepton", &ep2_kIsLepton, &b_ep2_kIsLepton);
  fChain->SetBranchAddress("ep2_kIsUsed", &ep2_kIsUsed, &b_ep2_kIsUsed);
  fChain->SetBranchAddress("ep2_mdcchi2", &ep2_mdcchi2, &b_ep2_mdcchi2);
  fChain->SetBranchAddress("ep2_oa_hadr", &ep2_oa_hadr, &b_ep2_oa_hadr);
  fChain->SetBranchAddress("ep2_oa_lept", &ep2_oa_lept, &b_ep2_oa_lept);
  fChain->SetBranchAddress("ep2_p", &ep2_p, &b_ep2_p);
  fChain->SetBranchAddress("ep2_phi", &ep2_phi, &b_ep2_phi);
  fChain->SetBranchAddress("ep2_q", &ep2_q, &b_ep2_q);
  fChain->SetBranchAddress("ep2_r", &ep2_r, &b_ep2_r);
  fChain->SetBranchAddress("ep2_resolution", &ep2_resolution, &b_ep2_resolution);
  fChain->SetBranchAddress("ep2_rkchi2", &ep2_rkchi2, &b_ep2_rkchi2);
  fChain->SetBranchAddress("ep2_sector", &ep2_sector, &b_ep2_sector);
  fChain->SetBranchAddress("ep2_shw_sum0", &ep2_shw_sum0, &b_ep2_shw_sum0);
  fChain->SetBranchAddress("ep2_shw_sum1", &ep2_shw_sum1, &b_ep2_shw_sum1);
  fChain->SetBranchAddress("ep2_shw_sum2", &ep2_shw_sum2, &b_ep2_shw_sum2);
  /*
    fChain->SetBranchAddress("ep2_sim_corrflag", &ep2_sim_corrflag, &b_ep2_sim_corrflag);
    fChain->SetBranchAddress("ep2_sim_geninfo", &ep2_sim_geninfo, &b_ep2_sim_geninfo);
    fChain->SetBranchAddress("ep2_sim_geninfo1", &ep2_sim_geninfo1, &b_ep2_sim_geninfo1);
    fChain->SetBranchAddress("ep2_sim_geninfo2", &ep2_sim_geninfo2, &b_ep2_sim_geninfo2);
    fChain->SetBranchAddress("ep2_sim_genweight", &ep2_sim_genweight, &b_ep2_sim_genweight);
    fChain->SetBranchAddress("ep2_sim_id", &ep2_sim_id, &b_ep2_sim_id);
    fChain->SetBranchAddress("ep2_sim_iscommon", &ep2_sim_iscommon, &b_ep2_sim_iscommon);
    fChain->SetBranchAddress("ep2_sim_mediumid", &ep2_sim_mediumid, &b_ep2_sim_mediumid);
    fChain->SetBranchAddress("ep2_sim_p", &ep2_sim_p, &b_ep2_sim_p);
    fChain->SetBranchAddress("ep2_sim_parentid", &ep2_sim_parentid, &b_ep2_sim_parentid);
    fChain->SetBranchAddress("ep2_sim_primaryflag", &ep2_sim_primaryflag, &b_ep2_sim_primaryflag);
    fChain->SetBranchAddress("ep2_sim_processid", &ep2_sim_processid, &b_ep2_sim_processid);
    fChain->SetBranchAddress("ep2_sim_px", &ep2_sim_px, &b_ep2_sim_px);
    fChain->SetBranchAddress("ep2_sim_py", &ep2_sim_py, &b_ep2_sim_py);
    fChain->SetBranchAddress("ep2_sim_pz", &ep2_sim_pz, &b_ep2_sim_pz);
    fChain->SetBranchAddress("ep2_sim_vertexx", &ep2_sim_vertexx, &b_ep2_sim_vertexx);
    fChain->SetBranchAddress("ep2_sim_vertexy", &ep2_sim_vertexy, &b_ep2_sim_vertexy);
    fChain->SetBranchAddress("ep2_sim_vertexz", &ep2_sim_vertexz, &b_ep2_sim_vertexz);
  */  
  fChain->SetBranchAddress("ep2_system", &ep2_system, &b_ep2_system);
  fChain->SetBranchAddress("ep2_theta", &ep2_theta, &b_ep2_theta);
  fChain->SetBranchAddress("ep2_tof_exp", &ep2_tof_exp, &b_ep2_tof_exp);
  fChain->SetBranchAddress("ep2_tof_mom", &ep2_tof_mom, &b_ep2_tof_mom);
  fChain->SetBranchAddress("ep2_tof_new", &ep2_tof_new, &b_ep2_tof_new);
  fChain->SetBranchAddress("ep2_tofino_mult", &ep2_tofino_mult, &b_ep2_tofino_mult);
  fChain->SetBranchAddress("ep2_track_length", &ep2_track_length, &b_ep2_track_length);
  fChain->SetBranchAddress("ep2_z", &ep2_z, &b_ep2_z);
  fChain->SetBranchAddress("ep1_beta", &ep1_beta, &b_ep1_beta);
  fChain->SetBranchAddress("ep1_beta_new", &ep1_beta_new, &b_ep1_beta_new);
  fChain->SetBranchAddress("ep1_dedx_in", &ep1_dedx_in, &b_ep1_dedx_in);
  fChain->SetBranchAddress("ep1_dedx_in_sigma", &ep1_dedx_in_sigma, &b_ep1_dedx_in_sigma);
  fChain->SetBranchAddress("ep1_dedx_mdc", &ep1_dedx_mdc, &b_ep1_dedx_mdc);
  fChain->SetBranchAddress("ep1_dedx_mdc_sigma", &ep1_dedx_mdc_sigma, &b_ep1_dedx_mdc_sigma);
  fChain->SetBranchAddress("ep1_dedx_out", &ep1_dedx_out, &b_ep1_dedx_out);
  fChain->SetBranchAddress("ep1_dedx_out_sigma", &ep1_dedx_out_sigma, &b_ep1_dedx_out_sigma);
  fChain->SetBranchAddress("ep1_dedx_tof", &ep1_dedx_tof, &b_ep1_dedx_tof);
  fChain->SetBranchAddress("ep1_id", &ep1_id, &b_ep1_id);
  fChain->SetBranchAddress("ep1_isring", &ep1_isring, &b_ep1_isring);
  fChain->SetBranchAddress("ep1_kIsLepton", &ep1_kIsLepton, &b_ep1_kIsLepton);
  fChain->SetBranchAddress("ep1_kIsUsed", &ep1_kIsUsed, &b_ep1_kIsUsed);
  fChain->SetBranchAddress("ep1_mdcchi2", &ep1_mdcchi2, &b_ep1_mdcchi2);
  fChain->SetBranchAddress("ep1_oa_hadr", &ep1_oa_hadr, &b_ep1_oa_hadr);
  fChain->SetBranchAddress("ep1_oa_lept", &ep1_oa_lept, &b_ep1_oa_lept);
  fChain->SetBranchAddress("ep1_p", &ep1_p, &b_ep1_p);
  fChain->SetBranchAddress("ep1_phi", &ep1_phi, &b_ep1_phi);
  fChain->SetBranchAddress("ep1_q", &ep1_q, &b_ep1_q);
  fChain->SetBranchAddress("ep1_r", &ep1_r, &b_ep1_r);
  fChain->SetBranchAddress("ep1_resolution", &ep1_resolution, &b_ep1_resolution);
  fChain->SetBranchAddress("ep1_rkchi2", &ep1_rkchi2, &b_ep1_rkchi2);
  fChain->SetBranchAddress("ep1_sector", &ep1_sector, &b_ep1_sector);
  fChain->SetBranchAddress("ep1_shw_sum0", &ep1_shw_sum0, &b_ep1_shw_sum0);
  fChain->SetBranchAddress("ep1_shw_sum1", &ep1_shw_sum1, &b_ep1_shw_sum1);
  fChain->SetBranchAddress("ep1_shw_sum2", &ep1_shw_sum2, &b_ep1_shw_sum2);
  /*
    fChain->SetBranchAddress("ep1_sim_corrflag", &ep1_sim_corrflag, &b_ep1_sim_corrflag);
    fChain->SetBranchAddress("ep1_sim_geninfo", &ep1_sim_geninfo, &b_ep1_sim_geninfo);
    fChain->SetBranchAddress("ep1_sim_geninfo1", &ep1_sim_geninfo1, &b_ep1_sim_geninfo1);
    fChain->SetBranchAddress("ep1_sim_geninfo2", &ep1_sim_geninfo2, &b_ep1_sim_geninfo2);
    fChain->SetBranchAddress("ep1_sim_genweight", &ep1_sim_genweight, &b_ep1_sim_genweight);
    fChain->SetBranchAddress("ep1_sim_id", &ep1_sim_id, &b_ep1_sim_id);
    fChain->SetBranchAddress("ep1_sim_iscommon", &ep1_sim_iscommon, &b_ep1_sim_iscommon);
    fChain->SetBranchAddress("ep1_sim_mediumid", &ep1_sim_mediumid, &b_ep1_sim_mediumid);
    fChain->SetBranchAddress("ep1_sim_p", &ep1_sim_p, &b_ep1_sim_p);
    fChain->SetBranchAddress("ep1_sim_parentid", &ep1_sim_parentid, &b_ep1_sim_parentid);
    fChain->SetBranchAddress("ep1_sim_primaryflag", &ep1_sim_primaryflag, &b_ep1_sim_primaryflag);
    fChain->SetBranchAddress("ep1_sim_processid", &ep1_sim_processid, &b_ep1_sim_processid);
    fChain->SetBranchAddress("ep1_sim_px", &ep1_sim_px, &b_ep1_sim_px);
    fChain->SetBranchAddress("ep1_sim_py", &ep1_sim_py, &b_ep1_sim_py);
    fChain->SetBranchAddress("ep1_sim_pz", &ep1_sim_pz, &b_ep1_sim_pz);
    fChain->SetBranchAddress("ep1_sim_vertexx", &ep1_sim_vertexx, &b_ep1_sim_vertexx);
    fChain->SetBranchAddress("ep1_sim_vertexy", &ep1_sim_vertexy, &b_ep1_sim_vertexy);
    fChain->SetBranchAddress("ep1_sim_vertexz", &ep1_sim_vertexz, &b_ep1_sim_vertexz);
  */  
  fChain->SetBranchAddress("ep1_system", &ep1_system, &b_ep1_system);
  fChain->SetBranchAddress("ep1_theta", &ep1_theta, &b_ep1_theta);
  fChain->SetBranchAddress("ep1_tof_exp", &ep1_tof_exp, &b_ep1_tof_exp);
  fChain->SetBranchAddress("ep1_tof_mom", &ep1_tof_mom, &b_ep1_tof_mom);
  fChain->SetBranchAddress("ep1_tof_new", &ep1_tof_new, &b_ep1_tof_new);
  fChain->SetBranchAddress("ep1_tofino_mult", &ep1_tofino_mult, &b_ep1_tofino_mult);
  fChain->SetBranchAddress("ep1_track_length", &ep1_track_length, &b_ep1_track_length);
  fChain->SetBranchAddress("ep1_z", &ep1_z, &b_ep1_z);
  fChain->SetBranchAddress("runnumber", &runnumber, &b_runnumber);
  fChain->SetBranchAddress("totalmult", &totalmult, &b_totalmult);
  fChain->SetBranchAddress("trigbit", &trigbit, &b_trigbit);
  fChain->SetBranchAddress("trigdec", &trigdec, &b_trigdec);
  fChain->SetBranchAddress("trigdownscale", &trigdownscale, &b_trigdownscale);
  fChain->SetBranchAddress("trigdownscaleflag", &trigdownscaleflag, &b_trigdownscaleflag);
  Notify();
}

Bool_t PPimEpEp::Notify()
{
  // The Notify() function is called when a new file is opened. This
  // can be either for a new TTree in a TChain or when when a new TTree
  // is started when using PROOF. It is normally not necessary to make changes
  // to the generated code, but the routine can be extended by the
  // user if needed. The return value is currently not used.

  return kTRUE;
}

void PPimEpEp::Show(Long64_t entry)
{
  // Print contents of entry.
  // If entry is not specified, print current entry
  if (!fChain) return;
  fChain->Show(entry);
}
Int_t PPimEpEp::Cut(Long64_t entry)
{
  // This function may be called from Loop.
  // returns  1 if entry is accepted.
  // returns -1 otherwise.
  return 1;
}
