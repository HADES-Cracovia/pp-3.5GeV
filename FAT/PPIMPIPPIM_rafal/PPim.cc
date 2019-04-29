#include "PPim.h"
#include "data.h"
#include <iostream>
#include <TLorentzVector.h>
#include <TVector3.h>
#include "hntuple.h"
#include "hparticletool.h"
#include <hgeomvector.h>

using namespace std;
using namespace PATData;

TVector3 vertex_new(Double_t z1,Double_t r1,TLorentzVector v1, Double_t z2,Double_t r2,TLorentzVector v2)
{
  TVector3 out;
  HGeomVector ver;
  HGeomVector base_1, base_2, dir_1, dir_2;
  HParticleTool p_tool_1;
  double phi1, phi2;
  //r1=TMath::Abs(r1);
  //r2=TMath::Abs(r2);

  phi1=(p_tool_1.getLabPhiDeg(v1))*D2R;
  phi2=(p_tool_1.getLabPhiDeg(v2))*D2R;

  cout<<"przy liczeniu odleglosci :"<<endl;
  cout<<"phi1 :"<<phi1*R2D<<" theta1 "<<v1.Theta()*R2D<<endl;
  cout<<"phi2 :"<<phi2*R2D<<" theta2 "<<v2.Theta()*R2D<<endl;
  cout<<"r1: "<<r1<<" z1: "<<z1<<endl;
  cout<<"r2: "<<r2<<" z2: "<<z2<<endl;
  
  p_tool_1.calcSegVector(z1,r1,phi1,v1.Theta(),base_1,dir_1);
  p_tool_1.calcSegVector(z2,r2,phi2,v2.Theta(),base_2,dir_2);

  cout<<"base_1 "<<endl;
  cout<<base_1.getX()<<" "<<base_1.getY()<<" "<<base_1.getZ()<<endl;
  cout<<"dir_1 "<<endl;
  cout<<dir_1.getX()<<" "<<dir_1.getY()<<" "<<dir_1.getZ()<<endl;
  cout<<"base_2 "<<endl;
  cout<<base_2.getX()<<" "<<base_2.getY()<<" "<<base_2.getZ()<<endl;
  cout<<"dir_2 "<<endl;
  cout<<dir_2.getX()<<" "<<dir_2.getY()<<" "<<dir_2.getZ()<<endl;
 
  ver=p_tool_1.calcVertexAnalytical(base_1,dir_1,base_2,dir_2);
  cout<<"vertex in function: "<<ver.getX()<<" " <<ver.getY()<<" "<<ver.getZ()<<endl;
  
  out.SetXYZ(ver.getX(),ver.getY(),ver.getZ());
  //cout<<"out vector"; out.Print(); cout<<endl;
  return out; 
}


void PPim::Loop()
{
  static long licznik = 0;

  if (fChain == 0) return;

  Long64_t nentries = fChain->GetEntries();

  Long64_t nbytes = 0, nb = 0;
  for (Long64_t jentry=0; jentry<nentries;jentry++)
    {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;

      ++licznik;
      //      if ((licznik % 100000)==0) cout << "Events: " << licznik << endl;

      if(jentry%10000==0 && isBest!=-1)
	{
	  cout << "netry no. "<< jentry<<" from "<<nentries;
	  //cout<<" isBest: "<< isBest<<" event: "<<event;
	  cout<<endl;
	}
      if(isBest!=1)
	continue;

      //double F = 1.006;
      double F=1.0;
      TVector3 v1, v2, v3;
      v2.SetXYZ(F*p_p*sin(D2R*p_theta)*cos(D2R*p_phi),F*p_p*sin(D2R*p_theta)*sin(D2R*p_phi),F*p_p*cos(D2R*p_theta));
      v3.SetXYZ(F*pim_p*sin(D2R*pim_theta)*cos(D2R*pim_phi),F*pim_p*sin(D2R*pim_theta)*sin(D2R*pim_phi),F*pim_p*cos(D2R*pim_theta));


      //cout<<"v2.X(),Y(),Z() :"<<v2.X()<<" "<<v2.Y()<<" "<<v2.Z()<<endl;
      //cout<<"px,py,pz       :"<<p_sim_px<<" "<<p_sim_py<<" "<<p_sim_pz<<endl<<endl;
      //TVector3 r1, r2;
      //      r1.SetXYZ(sin(D2R*p_theta_rich)*cos(D2R*p_phi_rich),sin(D2R*p_theta_rich)*sin(D2R*p_phi_rich),cos(D2R*p_theta_rich));
      //r2.SetXYZ(sin(D2R*pim_theta_rich)*cos(D2R*pim_phi_rich),sin(D2R*pim_theta_rich)*sin(D2R*pim_phi_rich),cos(D2R*pim_theta_rich));

      p->SetVectM( v2, 938.272013 );
      pi->SetVectM( v3, 139.57018 );

      *gammappi = *p + *pi;
      *ppi = *p + *pi;
      *p_delta = *p;
      *pi_delta = *pi;
      *ppi_miss = *beam - *p - *pi;

      double m2_inv_ppi = gammappi->M2();
      double m_inv_ppim = gammappi->M();
      double oa = R2D * openingangle(*p, *pi);
      //double oa_rich = R2D * openingangle(r1, r2);

      double p_mass = p_p*p_p * (  1. / (p_beta*p_beta)  - 1. ) ;
      double pi_mass = pim_p*pim_p * (  1. / (pim_beta*pim_beta)  - 1. ) ;

      TVector3 ver_p_pim=vertex(p_r,p_z,*p,pim_r,pim_z,*pi);
      double dist_p_pim=trackDistance(p_r,p_z,*p,pim_r,pim_z,*pi);
      //	  cout << "opening angle = " << oa << endl;
      /*
	cout<<"main part"<<endl;
	HParticleTool p_tool;  
	cout<<"p_phi: "<<p_phi<<endl;
	cout<<"p_phi by 3 vector: "<<v2.Phi()*R2D<<endl;
	cout<<"p_phi by 4 vector: "<<p->Phi()*R2D<<endl;
	cout<<"p_phi by getLabPhiDeg: "<<p_tool.getLabPhiDeg(*p)<<endl;
	cout<<"p_phi by phiLabToPhiSecDeg: "<<p_tool.phiLabToPhiSecDeg(v2.Phi()*R2D)<<endl;
	cout<<"p_phi by Atan(x/y): "<<TMath::ATan(v2.Y()/v2.X())*R2D<<endl;
	cout<<"p_phi by Atan(px_sim/py_sim): "<<TMath::ATan(p_sim_py/p_sim_px)*R2D<<endl;

	cout<<"p_theta: "<<p_theta<<endl;
	cout<<"p_theta by 3 vector: "<<v2.Theta()*R2D<<endl;
	cout<<"p_theta by 4 vector: "<<p->Theta()*R2D<<endl;
	//cout<<"p_theta by getLabThetaDeg: "<<p_tool.getLabThetaDeg(*p)<<endl;
	//cout<<"p_theta by thetaLabToThetaSecDeg: "<<p_tool.thetaLabToThetaSecDeg(v2.Theta()*R2D)<<endl;
	//cout<<"p_theta by Atan(x/y): "<<TMath::ATan(v2.Y()/v2.X())*R2D<<endl;
	v2.Print();
	cout<<endl;
      */
      ACC = 1.;
      EFF = 1.;


      gammappi->Boost(0., 0., -(beam->Beta()));
      p_delta->Boost(0., 0., -(beam->Beta()));
      pi_delta->Boost(0., 0., -(beam->Beta()));

      p_delta->Boost( -gammappi->Px()/gammappi->E(), -gammappi->Py()/gammappi->E(), -gammappi->Pz()/gammappi->E());
      pi_delta->Boost( -gammappi->Px()/gammappi->E(), -gammappi->Py()/gammappi->E(), -gammappi->Pz()/gammappi->E());

      //cout << "Poczatek obliczen..." << endl;

      //double ang_cut = 0.;
      double ang_cut = 9.;

      double close_cut = 9.;
      double nonfit_close_cut = -4.;
      //double close_cut = 0.;
      //double nonfit_close_cut = 0.;
      //double close_cut = 4.;
      /*if(event==1018 ||
	 event==1048 ||
	 event==1111 ||
	 event==1943 ||
	 event==1988 ||
	 event==4464 ||
	 event==4471 ||
	 event==4487 ||
	 event==4498 */
	 /*
	 p_sim_parentid==18 && pim_sim_parentid==18 && p_sim_id==14 && pim_sim_id==9)
	{
	  cout<<"!!!!!!!!!!!!!!!!!!!!!"<<endl;
	  cout<<"event number: "<<event<<endl;
	  cout<<"p phi: "<<p_phi<<" p theta: "<<p_theta<<endl;
	  cout<<"p z: "<<p_z<<" p r: "<<p_r<<endl;
	  cout<<"pim phi: "<<pim_phi<<" pim theta: "<<pim_theta<<endl;
	  cout<<"pim z: "<<pim_z<<" pim r: "<<pim_r<<endl;
	  cout<<"vertex rec"; ver_p_pim.Print();
	  cout<<"distance: "<<dist_p_pim<<endl;
	  cout<<"ideal vertex p: "<<p_sim_vertexx<<" "<<p_sim_vertexy<<" "<<p_sim_vertexz<<endl;
	  cout<<"ideal vertex pim: "<<pim_sim_vertexx<<" "<<pim_sim_vertexy<<" "<<pim_sim_vertexz<<endl<<endl;
	  }*/
#ifdef FLANCH
      //insidePimS0 = (pPimS0 == 0) ? 0 : pPimS0->IsInside(pim_z,pim_theta);
      //insidePimS1 = (pPimS1 == 0) ? 0 : pPimS1->IsInside(pim_z,pim_theta);
      //insideEpS0 = (pPS0 == 0) ? 0 : pPS0->IsInside(p_z,p_theta);
      //insidePS1 = (pPS1 == 0) ? 0 : pPS1->IsInside(p_z,p_theta);
      //insidePimS0 = (pPimS0 == 0) ? 0 : pPimS0->IsInside(eVert_z,pim_theta);
      //insidePimS1 = (pPimS1 == 0) ? 0 : pPimS1->IsInside(eVert_z,pim_theta);
      //insidePS0 = (pPS0 == 0) ? 0 : pPS0->IsInside(eVert_z,p_theta);
      //insideEpS1 = (pPS1 == 0) ? 0 : pPS1->IsInside(eVert_z,p_theta);
#endif

      insideTarget = 1;

#ifdef RECTANG
      //insidePimS0 = (pim_theta > 50 && pim_z < -50 /* && pim_p<200.*/) ? 1 : 0;
      //insidePimS1 = (pim_theta > 50 && pim_z < -50 /* && pim_p<200.*/) ? 1 : 0;
      //insidePS0 = (p_theta > 50 && p_z < -50 /* && p_p<200.*/) ? 1 : 0;
      //insidePS1 = (p_theta > 50 && p_z < -50 /* && p_p<200.*/) ? 1 : 0;
#endif

      //#ifdef NOCUT
      //insidePimS0 = 0;
      //insidePimS1 = 0;
      //insidePS0 = 0;
      //insidePS1 = 0;
      //#endif


      NoLeptonP = !((p_oa_lept< close_cut&&p_oa_lept>0.0) &&p_oa_lept>nonfit_close_cut );
      NoHadronP = !(p_oa_hadr< close_cut &&p_oa_hadr>nonfit_close_cut );
      NoLeptonPI = !((pim_oa_lept< close_cut&&pim_oa_lept>0.0) &&pim_oa_lept>nonfit_close_cut );
      NoHadronPI = !(pim_oa_hadr< close_cut &&pim_oa_hadr>nonfit_close_cut );
      NoHadronP = 1;
      NoHadronPI = 1;

      /*
	NoLeptonP = 1;
	NoHadronP = 1;
	NoLeptonPI = 1;
	NoHadronPI = 1;
      */


      (*n_ppim)["isBest"]=isBest;
      //(*n_ppim)["isBest_new"]=isBest_new;
      (*n_ppim)["event"]=event;
      (*n_ppim)["hneg_mult"]=hneg_mult;
      (*n_ppim)["hpos_mult"]=hpos_mult;
      (*n_ppim)["eVert_x"]=eVert_x;
      (*n_ppim)["eVert_y"]=eVert_y;
      (*n_ppim)["eVert_z"]=eVert_z;
      (*n_ppim)["totalmult"]=totalmult;
      (*n_ppim)["trigdownscaleflag"]=trigdownscaleflag;
      (*n_ppim)["trigdownscale"]=trigdownscale;
      //(*n_ppim)["event_mult"]=event_mult;
      //(*n_ppim)["hypothesis"]=pim_no;
      //(*n_ppim)["hypothesis_quality"]=quality;
  
      (*n_ppim)["p_p"]=p_p;
      (*n_ppim)["p_theta"] = p_theta;
      (*n_ppim)["p_phi"] = p_phi;
      (*n_ppim)["p_beta"] = p_beta_new;
      (*n_ppim)["p_m"] = p_mass;
      (*n_ppim)["p_dedx"]=p_dedx_mdc;
      (*n_ppim)["p_q"]=p_q;
	  
      (*n_ppim)["p_sim_p"]=p_sim_p;
      (*n_ppim)["p_sim_iscommon"]=p_sim_iscommon;
      (*n_ppim)["p_sim_id"]=p_sim_id;
      (*n_ppim)["p_sim_parentid"]=p_sim_parentid;
      (*n_ppim)["p_sim_vertex_x"]=p_sim_vertexx;
      (*n_ppim)["p_sim_vertex_y"]=p_sim_vertexy;
      (*n_ppim)["p_sim_vertex_z"]=p_sim_vertexz;
	  
      (*n_ppim)["pim_p"]=pim_p;
      (*n_ppim)["pim_theta"] = pim_theta;
      (*n_ppim)["pim_phi"] = pim_phi;
      (*n_ppim)["pim_beta"] = pim_beta_new;
      (*n_ppim)["pim_m"] = pi_mass;
      (*n_ppim)["pim_dedx"]=pim_dedx_mdc;
      (*n_ppim)["pim_q"]=pim_q;
	  
      (*n_ppim)["pim_sim_p"]=pim_sim_p;
      (*n_ppim)["pim_sim_iscommon"]=pim_sim_iscommon;
      (*n_ppim)["pim_sim_id"]=pim_sim_id;
      (*n_ppim)["pim_sim_parentid"]=pim_sim_parentid;
      (*n_ppim)["pim_sim_vertex_x"]=pim_sim_vertexx;
      (*n_ppim)["pim_sim_vertex_y"]=pim_sim_vertexy;
      (*n_ppim)["pim_sim_vertex_z"]=pim_sim_vertexz;
	  	  
  
      (*n_ppim)["m_inv_p_pim"] = m_inv_ppim;
  
      (*n_ppim)["ver_p_pim_x"]=ver_p_pim.X();
      (*n_ppim)["ver_p_pim_y"]=ver_p_pim.Y();
      (*n_ppim)["ver_p_pim_z"]=ver_p_pim.Z();
  
      //  (*n_ppim)["mlp_output"]=mlp_output;
      //(*n_ppim)["mlp_response"]=mlp_response;  
      n_ppim->fill();
    } // end of main loop
} // eof Loop 



// -------------------------------------------------------------------------------------------------
PPim::PPim(TTree *tree) 
{
  // if parameter tree is not specified (or zero), connect the file
  // used to generate this class and read the Tree.

  //pim_acc = pim_acc_err = p_acc = p_acc_err = 0.;
  //pim_eff = pim_eff_err = p_eff = p_eff_err = 0.;

  if (tree == 0) {
	  
    TChain * chain = new TChain("PPim_ID","");
    //chain->Add("/lustre/nyx/hades/user/knowakow/PNB/PAT/FILES/sep08_all/list5/sum5.root/PPim_ID");
    //chain->Add("/lustre/nyx/hades/user/knowakow/PNB/PAT/FILES/sep08_all/list4/sum4.root/PPim_ID");
    //chain->Add("/lustre/nyx/hades/user/knowakow/PNB/PAT/FILES/sep08_all/list3/sum3.root/PPim_ID");
    //chain->Add("/lustre/nyx/hades/user/knowakow/PNB/PAT/FILES/sep08_all/list2/sum2.root/PPim_ID");
    //chain->Add("/lustre/nyx/hades/user/knowakow/PNB/PAT_ppim/FILES/day280/hadron.root/PPim_ID");
    /*
      chain->Add("/lustre/nyx/hades/user/knowakow/PP/PAT_1/FILES/ppimpippim_dedx/hadron00.root/PPim_ID");
      chain->Add("/lustre/nyx/hades/user/knowakow/PP/PAT_1/FILES/ppimpippim_dedx/hadron01.root/PPim_ID");
      chain->Add("/lustre/nyx/hades/user/knowakow/PP/PAT_1/FILES/ppimpippim_dedx/hadron02.root/PPim_ID");
      chain->Add("/lustre/nyx/hades/user/knowakow/PP/PAT_1/FILES/ppimpippim_dedx/hadron03.root/PPim_ID");
      chain->Add("/lustre/nyx/hades/user/knowakow/PP/PAT_1/FILES/ppimpippim_dedx/hadron04.root/PPim_ID");
      chain->Add("/lustre/nyx/hades/user/knowakow/PP/PAT_1/FILES/ppimpippim_dedx/hadron05.root/PPim_ID");
      chain->Add("/lustre/nyx/hades/user/knowakow/PP/PAT_1/FILES/ppimpippim_dedx/hadron06.root/PPim_ID");
      chain->Add("/lustre/nyx/hades/user/knowakow/PP/PAT_1/FILES/ppimpippim_dedx/hadron07.root/PPim_ID");
      chain->Add("/lustre/nyx/hades/user/knowakow/PP/PAT_1/FILES/ppimpippim_dedx/hadron08.root/PPim_ID");
      chain->Add("/lustre/nyx/hades/user/knowakow/PP/PAT_1/FILES/ppimpippim_dedx/hadron09.root/PPim_ID");
      chain->Add("/lustre/nyx/hades/user/knowakow/PP/PAT_1/FILES/ppimpippim_dedx/hadron10.root/PPim_ID");
      chain->Add("/lustre/nyx/hades/user/knowakow/PP/PAT_1/FILES/ppimpippim_dedx/hadron11.root/PPim_ID");
    */
    chain->Add("/lustre/nyx/hades/user/knowakow/PP/PAT_35_Rafal/FILES/all.root/PPim_ID");
    //chain->Add("/lustre/nyx/hades/user/knowakow/PP/PAT_35_Rafal/FILES/pp35Sigmap_chan_001_evt_50000_nfile_007_hgeant1_dst_hadron_out.root/PPim_ID");
    //chain->Add("/lustre/nyx/hades/user/knowakow/PP/PAT_1/FILES/ppimpippim/hadron12.root/PPim_ID");
    //chain->Add("/lustre/nyx/hades/user/knowakow/PP/PAT_1/FILES/ppimpippim/hadron13.root/PPim_ID");
    
    tree = chain; 
  }

  Init(tree);
}


// -------------------------------------------------------------------------------------------------
PPim::~PPim()
{
  if (!fChain) return;
  delete fChain->GetCurrentFile();
}

// -------------------------------------------------------------------------------------------------
Int_t PPim::GetEntry(Long64_t entry)
{
  // Read contents of entry.
  if (!fChain) return 0;
  return fChain->GetEntry(entry);
}

// -------------------------------------------------------------------------------------------------
Long64_t PPim::LoadTree(Long64_t entry)
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

// -------------------------------------------------------------------------------------------------
void PPim::Init(TTree *tree)
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
  fChain->SetBranchAddress("pim_system", &pim_system, &b_pim_system);
  fChain->SetBranchAddress("pim_theta", &pim_theta, &b_pim_theta);
  fChain->SetBranchAddress("pim_tof_exp", &pim_tof_exp, &b_pim_tof_exp);
  fChain->SetBranchAddress("pim_tof_mom", &pim_tof_mom, &b_pim_tof_mom);
  fChain->SetBranchAddress("pim_tof_new", &pim_tof_new, &b_pim_tof_new);
  fChain->SetBranchAddress("pim_tofino_mult", &pim_tofino_mult, &b_pim_tofino_mult);
  fChain->SetBranchAddress("pim_track_length", &pim_track_length, &b_pim_track_length);
  fChain->SetBranchAddress("pim_z", &pim_z, &b_pim_z);
  fChain->SetBranchAddress("runnumber", &runnumber, &b_runnumber);
  fChain->SetBranchAddress("totalmult", &totalmult, &b_totalmult);
  fChain->SetBranchAddress("trigbit", &trigbit, &b_trigbit);
  fChain->SetBranchAddress("trigdec", &trigdec, &b_trigdec);
  fChain->SetBranchAddress("trigdownscale", &trigdownscale, &b_trigdownscale);
  fChain->SetBranchAddress("trigdownscaleflag", &trigdownscaleflag, &b_trigdownscaleflag);

  Notify();
}

// -------------------------------------------------------------------------------------------------
Bool_t PPim::Notify()
{
  // The Notify() function is called when a new file is opened. This
  // can be either for a new TTree in a TChain or when when a new TTree
  // is started when using PROOF. It is normally not necessary to make changes
  // to the generated code, but the routine can be extended by the
  // user if needed. The return value is currently not used.

  return kTRUE;
}

// -------------------------------------------------------------------------------------------------
void PPim::Show(Long64_t entry)
{
  // Print contents of entry.
  // If entry is not specified, print current entry
  if (!fChain) return;
  fChain->Show(entry);
}

// -------------------------------------------------------------------------------------------------
Int_t PPim::Cut(Long64_t entry)
{
  // This function may be called from Loop.
  // returns  1 if entry is accepted.
  // returns -1 otherwise.
  return 1;
}
