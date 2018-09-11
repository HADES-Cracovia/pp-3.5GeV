#include "PPipPim.h"
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

void PPipPim::Loop()
{
  static long licznik = 0;
  if (fChain == 0) return;

  Long64_t nentries = fChain->GetEntriesFast();

  Long64_t nbytes = 0, nb = 0;
  for (Long64_t jentry=0; jentry<nentries;jentry++)
    {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;
      
      ++licznik;
      if ((licznik % 100000)==0) cout << "Events: " << licznik << endl;


      double F = 1.006;
      TVector3 v1, v2, v3, v4;// v5;
      v2.SetXYZ(F*p_p*sin(D2R*p_theta)*cos(D2R*p_phi),F*p_p*sin(D2R*p_theta)*sin(D2R*p_phi),F*p_p*cos(D2R*p_theta));
      v3.SetXYZ(F*pim_p*sin(D2R*pim_theta)*cos(D2R*pim_phi),F*pim_p*sin(D2R*pim_theta)*sin(D2R*pim_phi),F*pim_p*cos(D2R*pim_theta));
      v4.SetXYZ(F*pip_p*sin(D2R*pip_theta)*cos(D2R*pip_phi),F*pip_p*sin(D2R*pip_theta)*sin(D2R*pip_phi),F*pip_p*cos(D2R*pip_theta));
      //v5.SetXYZ(F*pim2_p*sin(D2R*pim2_theta)*cos(D2R*pim2_phi),F*pim2_p*sin(D2R*pim2_theta)*sin(D2R*pim2_phi),F*pim2_p*cos(D2R*pim2_theta));
      
      TVector3 r1, r2, r3,r4;
      r1.SetXYZ(sin(D2R*p_theta_rich)*cos(D2R*p_phi_rich),sin(D2R*p_theta_rich)*sin(D2R*p_phi_rich),cos(D2R*p_theta_rich));
      r2.SetXYZ(sin(D2R*pim_theta_rich)*cos(D2R*pim_phi_rich),sin(D2R*pim_theta_rich)*sin(D2R*pim_phi_rich),cos(D2R*pim_theta_rich));
      r3.SetXYZ(sin(D2R*pip_theta_rich)*cos(D2R*pip_phi_rich),sin(D2R*pip_theta_rich)*sin(D2R*pip_phi_rich),cos(D2R*pip_theta_rich));
      //r4.SetXYZ(sin(D2R*pim2_theta_rich)*cos(D2R*pim2_phi_rich),sin(D2R*pim2_theta_rich)*sin(D2R*pim2_phi_rich),cos(D2R*pim2_theta_rich));
      p->SetVectM( v2, 938.272013 );
      pim->SetVectM( v3, 139.57018 );
      pip->SetVectM( v4, 139.57018 );
      //pim2->SetVectM( v5, 139.57018 );
     
      *gammappim = *p + *pim;
      //*gammappim2 = *p + *pim2;
      *gammapimpip= *pim + *pip;
      //*gammapim2pip= *pim2 + *pip;
      *gammappimpip=*pim + *pip + *p;
      
      //*ppim = *p + *pim;
      //*p_delta = *p;
      //*pim_delta = *pim;
      //*ppim_miss = *beam - *p - *pim;

      //double m2_inv_ppim = gammappim->M2();
      double m_inv_ppim = gammappim->M();
      //double m_inv_ppim2 = gammappim2->M();
      double m_inv_pippim = gammapimpip->M();
      //double m_inv_pippim2 = gammapim2pip->M();
      double m_inv_ppimpip = gammappimpip->M();
      double oa = R2D * openingangle(*p, *pim);
      double oa_rich = R2D * openingangle(r1, r2);

      double p_mass = p_p*p_p * (  1. / (p_beta*p_beta)  - 1. ) ;
      double pi_mass = pim_p*pim_p * (  1. / (pim_beta*pim_beta)  - 1. ) ;

      double missing_energy=beam->E()-gammappimpip->E();
      
      double d_dist_p_pim=trackDistance(p_r,p_z,v2,pim_r,pim_z,v3);
      //double d_dist_p_pim2=trackDistance(p_r,p_z,v2,pim2_r,pim2_z,v5);
      double d_dist_pip_pim=trackDistance(pip_r,pip_z,v4,pim_r,pim_z,v3);
      //double d_dist_pip_pim2=trackDistance(pip_r,pip_z,v4,pim2_r,pim2_z,v5);

      TVector3 ver_p_pim=vertex(p_r,p_z,v2,pim_r,pim_z,v3);
      //TVector3 ver_p_pim2=vertex(p_r,p_z,v2,pim2_r,pim2_z,v5);
      TVector3 ver_pip_pim=vertex(pip_r,pip_z,v4,pim_r,pim_z,v3);
      //TVector3 ver_pip_pim2=vertex(pip_r,pip_z,v4,pim2_r,pim2_z,v5);
      

      //	  cout << "opening angle = " << oa << endl;

      ACC = 1.;
      EFF = 1.;

      /*
	gammappi->Boost(0., 0., -(beam->Beta()));
	p_delta->Boost(0., 0., -(beam->Beta()));
	pi_delta->Boost(0., 0., -(beam->Beta()));

	p_delta->Boost( -gammappi->Px()/gammappi->E(), -gammappi->Py()/gammappi->E(), -gammappi->Pz()/gammappi->E());
	pi_delta->Boost( -gammappi->Px()/gammappi->E(), -gammappi->Py()/gammappi->E(), -gammappi->Pz()/gammappi->E());
      */
      //cout << "Poczatek obliczen..." << endl;

      //double ang_cut = 0.;
      //double ang_cut = 9.;

      //double close_cut = 9.;
      //double nonfit_close_cut = -4.;
      //double close_cut = 0.;
      //double nonfit_close_cut = 0.;
      //double close_cut = 4.;


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


      //NoLeptonP = !((p_oa_lept< close_cut&&p_oa_lept>0.0) &&p_oa_lept>nonfit_close_cut );
      //NoHadronP = !(p_oa_hadr< close_cut &&p_oa_hadr>nonfit_close_cut );
      //NoLeptonPI = !((pim_oa_lept< close_cut&&pim_oa_lept>0.0) &&pim_oa_lept>nonfit_close_cut );
      //NoHadronPI = !(pim_oa_hadr< close_cut &&pim_oa_hadr>nonfit_close_cut );
      //NoHadronP = 1;
      //NoHadronPI = 1;

      /*
	NoLeptonP = 1;
	NoHadronP = 1;
	NoLeptonPI = 1;
	NoHadronPI = 1;
      */
      bool c_mass=(m_inv_ppim<1120 &&m_inv_ppim>1110);
      //bool c_mass2=(m_inv_ppim2<1120 &&m_inv_ppim2>1110);

      bool dist=(d_dist_p_pim<15 && ver_p_pim.Z()>20);
      //bool dist2=(d_dist_p_pim2<15 && ver_p_pim2.Z()>20);
      
      if(isBest==1)
	{
	  p_p_beta->Fill(p_p,p_beta_new);
	  pim_p_beta->Fill(pim_p,pim_beta_new);
	  //pim_p_beta->Fill(pim_p,pim_beta_new);
	  pip_p_beta->Fill(pip_p,pip_beta_new);
	  
	  p_pim_mass->Fill(m_inv_ppim);
	  //p_pim_mass->Fill(m_inv_ppim2);
	  p_pim_mass->Fill(m_inv_ppim);
	  //p_pim2_mass->Fill(m_inv_ppim2);

	  //pim_pip_mass->Fill(m_inv_pippim2);
	  pim_pip_mass->Fill(m_inv_pippim);
	  //pim2_pip_mass->Fill(m_inv_pippim2);
	  //pim_pip_mass->Fill(m_inv_pippim);

	  p_pim_pip_mass->Fill(m_inv_ppimpip);

	  //dist_p_pim_pim_pip->Fill(d_dist_p_pim,d_dist_pip_pim);
	  //dist_p_pim_pim_pip->Fill(d_dist_p_pim2,d_dist_pip_pim2);
	  //dist_p_pim->Fill(d_dist_p_pim2);
	  dist_p_pim->Fill(d_dist_p_pim);
	  //dist_pim_pip->Fill(d_dist_pip_pim2);
	  dist_pim_pip->Fill(d_dist_pip_pim);

	  miss_energy->Fill(missing_energy);
	}
      if(isBest==1 && dist)
	{
	  DL_p_pim_mass->Fill(m_inv_ppim);
	  //DL_p_pim_mass->Fill(m_inv_ppim);
	  DL_pim_mass->Fill(m_inv_pippim);
	  //DL_pim2_pip_mass->Fill(m_inv_pippim2);
	  DL_p_pim_pip_mass->Fill(m_inv_ppimpip);
	  //DL_dist_p_pim_pim_pip->Fill(d_dist_p_pim,d_dist_pip_pim2);
	  DL_dist_p_pim->Fill(d_dist_p_pim);
	  DL_dist_pim_pip->Fill(d_dist_pip_pim);
	  DL_miss_energy->Fill(missing_energy);
	}
      if(isBest==1 && dist && c_mass)
	{
	  DML_p_pim_mass->Fill(m_inv_ppim);
	  //DML_p_pim_mass->Fill(m_inv_ppim);
	  DML_pim_pip_mass->Fill(m_inv_pippim);
	  //L_pim2_pip_mass->Fill(m_inv_pippim2);
	  DML_p_pim_pip_mass->Fill(m_inv_ppimpip);
	  //DML_dist_p_pim_pim_pip->Fill(d_dist_p_pim,d_dist_pip_pim2);
	  DML_dist_p_pim->Fill(d_dist_p_pim);
	  DML_dist_pim_pip->Fill(d_dist_pip_pim);
	  DML_miss_energy->Fill(missing_energy);
	}
    }
}


PPipPim::PPipPim(TTree *tree)  
{
  // if parameter tree is not specified (or zero), connect the file
  // used to generate this class and read the Tree.
  /*if (tree == 0) {
    TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("be08280235056_dst_gen1_sep08_hadron_out.root");
    if (!f || !f->IsOpen()) {
    f = new TFile("be08280235056_dst_gen1_sep08_hadron_out.root");
    }
    f->GetObject("PPipPim",tree);

    }
    Init(tree);*/
  if (tree == 0) {
	  
    TChain * chain = new TChain("PPipPim_ID","");
    //chain->Add("/lustre/nyx/hades/user/knowakow/PNB/PAT_ppim/FILES/day280/hadron.root/PPipPim_ID");
    chain->Add("/lustre/nyx/hades/user/knowakow/PP/PAT_1/FILES/ppippim/hadron00.root/PPimPip_ID");
    ///*
    chain->Add("/lustre/nyx/hades/user/knowakow/PP/PAT_1/FILES/ppippim/hadron01.root/PPimPip_ID");
    chain->Add("/lustre/nyx/hades/user/knowakow/PP/PAT_1/FILES/ppippim/hadron02.root/PPimPip_ID");
    chain->Add("/lustre/nyx/hades/user/knowakow/PP/PAT_1/FILES/ppippim/hadron03.root/PPimPip_ID");
    chain->Add("/lustre/nyx/hades/user/knowakow/PP/PAT_1/FILES/ppippim/hadron04.root/PPimPip_ID");
    chain->Add("/lustre/nyx/hades/user/knowakow/PP/PAT_1/FILES/ppippim/hadron05.root/PPimPip_ID");
    //chain->Add("/lustre/nyx/hades/user/knowakow/PP/PAT_1/FILES/ppippim/hadron06.root/PPimPip_ID");
    chain->Add("/lustre/nyx/hades/user/knowakow/PP/PAT_1/FILES/ppippim/hadron07.root/PPimPip_ID");
    chain->Add("/lustre/nyx/hades/user/knowakow/PP/PAT_1/FILES/ppippim/hadron08.root/PPimPip_ID");
    chain->Add("/lustre/nyx/hades/user/knowakow/PP/PAT_1/FILES/ppippim/hadron09.root/PPimPip_ID");
    chain->Add("/lustre/nyx/hades/user/knowakow/PP/PAT_1/FILES/ppippim/hadron10.root/PPimPip_ID");
    //chain->Add("/lustre/nyx/hades/user/knowakow/PP/PAT_1/FILES/ppippim/hadron11.root/PPimPip_ID");
    //*/  
    
    tree = chain; 
  }

  Init(tree);
}

PPipPim::~PPipPim()
{
  if (!fChain) return;
  delete fChain->GetCurrentFile();
}

Int_t PPipPim::GetEntry(Long64_t entry)
{
  // Read contents of entry.
  if (!fChain) return 0;
  return fChain->GetEntry(entry);
}
Long64_t PPipPim::LoadTree(Long64_t entry)
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

void PPipPim::Init(TTree *tree)
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

  fChain->SetBranchAddress("isBest", &isBest, &b_isBest);
  fChain->SetBranchAddress("p_beta", &p_beta, &b_p_beta);
  fChain->SetBranchAddress("p_beta_new", &p_beta_new, &b_p_beta_new);
  fChain->SetBranchAddress("p_btChargeRing", &p_btChargeRing, &b_p_btChargeRing);
  fChain->SetBranchAddress("p_btChargeSum", &p_btChargeSum, &b_p_btChargeSum);
  fChain->SetBranchAddress("p_btChi2", &p_btChi2, &b_p_btChi2);
  fChain->SetBranchAddress("p_btClusters", &p_btClusters, &b_p_btClusters);
  fChain->SetBranchAddress("p_btMaxima", &p_btMaxima, &b_p_btMaxima);
  fChain->SetBranchAddress("p_btMaximaCharge", &p_btMaximaCharge, &b_p_btMaximaCharge);
  fChain->SetBranchAddress("p_btMaximaChargeShared", &p_btMaximaChargeShared, &b_p_btMaximaChargeShared);
  fChain->SetBranchAddress("p_btMaximaChargeSharedFragment", &p_btMaximaChargeSharedFragment, &b_p_btMaximaChargeSharedFragment);
  fChain->SetBranchAddress("p_btMaximaShared", &p_btMaximaShared, &b_p_btMaximaShared);
  fChain->SetBranchAddress("p_btMaximaSharedFragment", &p_btMaximaSharedFragment, &b_p_btMaximaSharedFragment);
  fChain->SetBranchAddress("p_btMeanDist", &p_btMeanDist, &b_p_btMeanDist);
  fChain->SetBranchAddress("p_btNearbyMaxima", &p_btNearbyMaxima, &b_p_btNearbyMaxima);
  fChain->SetBranchAddress("p_btNearbyMaximaShared", &p_btNearbyMaximaShared, &b_p_btNearbyMaximaShared);
  fChain->SetBranchAddress("p_btPadsClus", &p_btPadsClus, &b_p_btPadsClus);
  fChain->SetBranchAddress("p_btPadsRing", &p_btPadsRing, &b_p_btPadsRing);
  fChain->SetBranchAddress("p_btRingMatrix", &p_btRingMatrix, &b_p_btRingMatrix);
  fChain->SetBranchAddress("p_dedx_mdc", &p_dedx_mdc, &b_p_dedx_mdc);
  fChain->SetBranchAddress("p_dedx_tof", &p_dedx_tof, &b_p_dedx_tof);
  fChain->SetBranchAddress("p_id", &p_id, &b_p_id);
  fChain->SetBranchAddress("p_isBT", &p_isBT, &b_p_isBT);
  fChain->SetBranchAddress("p_isOffVertexClust", &p_isOffVertexClust, &b_p_isOffVertexClust);
  fChain->SetBranchAddress("p_isPrimaryVertex", &p_isPrimaryVertex, &b_p_isPrimaryVertex);
  fChain->SetBranchAddress("p_isUsedVertex", &p_isUsedVertex, &b_p_isUsedVertex);
  fChain->SetBranchAddress("p_isring", &p_isring, &b_p_isring);
  fChain->SetBranchAddress("p_isringmdc", &p_isringmdc, &b_p_isringmdc);
  fChain->SetBranchAddress("p_isringnomatch", &p_isringnomatch, &b_p_isringnomatch);
  fChain->SetBranchAddress("p_isringtrack", &p_isringtrack, &b_p_isringtrack);
  fChain->SetBranchAddress("p_kIsLepton", &p_kIsLepton, &b_p_kIsLepton);
  fChain->SetBranchAddress("p_kIsUsed", &p_kIsUsed, &b_p_kIsUsed);
  fChain->SetBranchAddress("p_mdcinnerchi2", &p_mdcinnerchi2, &b_p_mdcinnerchi2);
  fChain->SetBranchAddress("p_mdcouterchi2", &p_mdcouterchi2, &b_p_mdcouterchi2);
  fChain->SetBranchAddress("p_oa_hadr", &p_oa_hadr, &b_p_oa_hadr);
  fChain->SetBranchAddress("p_oa_lept", &p_oa_lept, &b_p_oa_lept);
  fChain->SetBranchAddress("p_p", &p_p, &b_p_p);
  fChain->SetBranchAddress("p_p_corr_em", &p_p_corr_em, &b_p_p_corr_em);
  fChain->SetBranchAddress("p_p_corr_ep", &p_p_corr_ep, &b_p_p_corr_ep);
  fChain->SetBranchAddress("p_p_corr_p", &p_p_corr_p, &b_p_p_corr_p);
  fChain->SetBranchAddress("p_p_corr_pim", &p_p_corr_pim, &b_p_p_corr_pim);
  fChain->SetBranchAddress("p_p_corr_pip", &p_p_corr_pip, &b_p_p_corr_pip);
  fChain->SetBranchAddress("p_phi", &p_phi, &b_p_phi);
  fChain->SetBranchAddress("p_phi_rich", &p_phi_rich, &b_p_phi_rich);
  fChain->SetBranchAddress("p_pid", &p_pid, &b_p_pid);
  fChain->SetBranchAddress("p_q", &p_q, &b_p_q);
  fChain->SetBranchAddress("p_r", &p_r, &b_p_r);
  fChain->SetBranchAddress("p_resolution", &p_resolution, &b_p_resolution);
  fChain->SetBranchAddress("p_resoultion", &p_resoultion, &b_p_resoultion);
  fChain->SetBranchAddress("p_rich_amp", &p_rich_amp, &b_p_rich_amp);
  fChain->SetBranchAddress("p_rich_centr", &p_rich_centr, &b_p_rich_centr);
  fChain->SetBranchAddress("p_rich_houtra", &p_rich_houtra, &b_p_rich_houtra);
  fChain->SetBranchAddress("p_rich_padnum", &p_rich_padnum, &b_p_rich_padnum);
  fChain->SetBranchAddress("p_rich_patmat", &p_rich_patmat, &b_p_rich_patmat);
  fChain->SetBranchAddress("p_rkchi2", &p_rkchi2, &b_p_rkchi2);
  fChain->SetBranchAddress("p_sector", &p_sector, &b_p_sector);
  fChain->SetBranchAddress("p_shw_sum0", &p_shw_sum0, &b_p_shw_sum0);
  fChain->SetBranchAddress("p_shw_sum1", &p_shw_sum1, &b_p_shw_sum1);
  fChain->SetBranchAddress("p_shw_sum2", &p_shw_sum2, &b_p_shw_sum2);
  fChain->SetBranchAddress("p_system", &p_system, &b_p_system);
  fChain->SetBranchAddress("p_theta", &p_theta, &b_p_theta);
  fChain->SetBranchAddress("p_theta_rich", &p_theta_rich, &b_p_theta_rich);
  fChain->SetBranchAddress("p_tof_mom", &p_tof_mom, &b_p_tof_mom);
  fChain->SetBranchAddress("p_tof_new", &p_tof_new, &b_p_tof_new);
  fChain->SetBranchAddress("p_tof_rec", &p_tof_rec, &b_p_tof_rec);
  fChain->SetBranchAddress("p_track_length", &p_track_length, &b_p_track_length);
  fChain->SetBranchAddress("p_tracklength", &p_tracklength, &b_p_tracklength);
  fChain->SetBranchAddress("p_z", &p_z, &b_p_z);
  fChain->SetBranchAddress("pim_beta", &pim_beta, &b_pim_beta);
  fChain->SetBranchAddress("pim_beta_new", &pim_beta_new, &b_pim_beta_new);
  fChain->SetBranchAddress("pim_btChargeRing", &pim_btChargeRing, &b_pim_btChargeRing);
  fChain->SetBranchAddress("pim_btChargeSum", &pim_btChargeSum, &b_pim_btChargeSum);
  fChain->SetBranchAddress("pim_btChi2", &pim_btChi2, &b_pim_btChi2);
  fChain->SetBranchAddress("pim_btClusters", &pim_btClusters, &b_pim_btClusters);
  fChain->SetBranchAddress("pim_btMaxima", &pim_btMaxima, &b_pim_btMaxima);
  fChain->SetBranchAddress("pim_btMaximaCharge", &pim_btMaximaCharge, &b_pim_btMaximaCharge);
  fChain->SetBranchAddress("pim_btMaximaChargeShared", &pim_btMaximaChargeShared, &b_pim_btMaximaChargeShared);
  fChain->SetBranchAddress("pim_btMaximaChargeSharedFragment", &pim_btMaximaChargeSharedFragment, &b_pim_btMaximaChargeSharedFragment);
  fChain->SetBranchAddress("pim_btMaximaShared", &pim_btMaximaShared, &b_pim_btMaximaShared);
  fChain->SetBranchAddress("pim_btMaximaSharedFragment", &pim_btMaximaSharedFragment, &b_pim_btMaximaSharedFragment);
  fChain->SetBranchAddress("pim_btMeanDist", &pim_btMeanDist, &b_pim_btMeanDist);
  fChain->SetBranchAddress("pim_btNearbyMaxima", &pim_btNearbyMaxima, &b_pim_btNearbyMaxima);
  fChain->SetBranchAddress("pim_btNearbyMaximaShared", &pim_btNearbyMaximaShared, &b_pim_btNearbyMaximaShared);
  fChain->SetBranchAddress("pim_btPadsClus", &pim_btPadsClus, &b_pim_btPadsClus);
  fChain->SetBranchAddress("pim_btPadsRing", &pim_btPadsRing, &b_pim_btPadsRing);
  fChain->SetBranchAddress("pim_btRingMatrix", &pim_btRingMatrix, &b_pim_btRingMatrix);
  fChain->SetBranchAddress("pim_dedx_mdc", &pim_dedx_mdc, &b_pim_dedx_mdc);
  fChain->SetBranchAddress("pim_dedx_tof", &pim_dedx_tof, &b_pim_dedx_tof);
  fChain->SetBranchAddress("pim_id", &pim_id, &b_pim_id);
  fChain->SetBranchAddress("pim_isBT", &pim_isBT, &b_pim_isBT);
  fChain->SetBranchAddress("pim_isOffVertexClust", &pim_isOffVertexClust, &b_pim_isOffVertexClust);
  fChain->SetBranchAddress("pim_isPrimaryVertex", &pim_isPrimaryVertex, &b_pim_isPrimaryVertex);
  fChain->SetBranchAddress("pim_isUsedVertex", &pim_isUsedVertex, &b_pim_isUsedVertex);
  fChain->SetBranchAddress("pim_isring", &pim_isring, &b_pim_isring);
  fChain->SetBranchAddress("pim_isringmdc", &pim_isringmdc, &b_pim_isringmdc);
  fChain->SetBranchAddress("pim_isringnomatch", &pim_isringnomatch, &b_pim_isringnomatch);
  fChain->SetBranchAddress("pim_isringtrack", &pim_isringtrack, &b_pim_isringtrack);
  fChain->SetBranchAddress("pim_kIsLepton", &pim_kIsLepton, &b_pim_kIsLepton);
  fChain->SetBranchAddress("pim_kIsUsed", &pim_kIsUsed, &b_pim_kIsUsed);
  fChain->SetBranchAddress("pim_mdcinnerchi2", &pim_mdcinnerchi2, &b_pim_mdcinnerchi2);
  fChain->SetBranchAddress("pim_mdcouterchi2", &pim_mdcouterchi2, &b_pim_mdcouterchi2);
  fChain->SetBranchAddress("pim_oa_hadr", &pim_oa_hadr, &b_pim_oa_hadr);
  fChain->SetBranchAddress("pim_oa_lept", &pim_oa_lept, &b_pim_oa_lept);
  fChain->SetBranchAddress("pim_p", &pim_p, &b_pim_p);
  fChain->SetBranchAddress("pim_p_corr_em", &pim_p_corr_em, &b_pim_p_corr_em);
  fChain->SetBranchAddress("pim_p_corr_ep", &pim_p_corr_ep, &b_pim_p_corr_ep);
  fChain->SetBranchAddress("pim_p_corr_p", &pim_p_corr_p, &b_pim_p_corr_p);
  fChain->SetBranchAddress("pim_p_corr_pim", &pim_p_corr_pim, &b_pim_p_corr_pim);
  fChain->SetBranchAddress("pim_p_corr_pip", &pim_p_corr_pip, &b_pim_p_corr_pip);
  fChain->SetBranchAddress("pim_phi", &pim_phi, &b_pim_phi);
  fChain->SetBranchAddress("pim_phi_rich", &pim_phi_rich, &b_pim_phi_rich);
  fChain->SetBranchAddress("pim_pid", &pim_pid, &b_pim_pid);
  fChain->SetBranchAddress("pim_q", &pim_q, &b_pim_q);
  fChain->SetBranchAddress("pim_r", &pim_r, &b_pim_r);
  fChain->SetBranchAddress("pim_resolution", &pim_resolution, &b_pim_resolution);
  fChain->SetBranchAddress("pim_resoultion", &pim_resoultion, &b_pim_resoultion);
  fChain->SetBranchAddress("pim_rich_amp", &pim_rich_amp, &b_pim_rich_amp);
  fChain->SetBranchAddress("pim_rich_centr", &pim_rich_centr, &b_pim_rich_centr);
  fChain->SetBranchAddress("pim_rich_houtra", &pim_rich_houtra, &b_pim_rich_houtra);
  fChain->SetBranchAddress("pim_rich_padnum", &pim_rich_padnum, &b_pim_rich_padnum);
  fChain->SetBranchAddress("pim_rich_patmat", &pim_rich_patmat, &b_pim_rich_patmat);
  fChain->SetBranchAddress("pim_rkchi2", &pim_rkchi2, &b_pim_rkchi2);
  fChain->SetBranchAddress("pim_sector", &pim_sector, &b_pim_sector);
  fChain->SetBranchAddress("pim_shw_sum0", &pim_shw_sum0, &b_pim_shw_sum0);
  fChain->SetBranchAddress("pim_shw_sum1", &pim_shw_sum1, &b_pim_shw_sum1);
  fChain->SetBranchAddress("pim_shw_sum2", &pim_shw_sum2, &b_pim_shw_sum2);
  fChain->SetBranchAddress("pim_system", &pim_system, &b_pim_system);
  fChain->SetBranchAddress("pim_theta", &pim_theta, &b_pim_theta);
  fChain->SetBranchAddress("pim_theta_rich", &pim_theta_rich, &b_pim_theta_rich);
  fChain->SetBranchAddress("pim_tof_mom", &pim_tof_mom, &b_pim_tof_mom);
  fChain->SetBranchAddress("pim_tof_new", &pim_tof_new, &b_pim_tof_new);
  fChain->SetBranchAddress("pim_tof_rec", &pim_tof_rec, &b_pim_tof_rec);
  fChain->SetBranchAddress("pim_track_length", &pim_track_length, &b_pim_track_length);
  fChain->SetBranchAddress("pim_tracklength", &pim_tracklength, &b_pim_tracklength);
  fChain->SetBranchAddress("pim_z", &pim_z, &b_pim_z);
  fChain->SetBranchAddress("pip_beta", &pip_beta, &b_pip_beta);
  fChain->SetBranchAddress("pip_beta_new", &pip_beta_new, &b_pip_beta_new);
  fChain->SetBranchAddress("pip_btChargeRing", &pip_btChargeRing, &b_pip_btChargeRing);
  fChain->SetBranchAddress("pip_btChargeSum", &pip_btChargeSum, &b_pip_btChargeSum);
  fChain->SetBranchAddress("pip_btChi2", &pip_btChi2, &b_pip_btChi2);
  fChain->SetBranchAddress("pip_btClusters", &pip_btClusters, &b_pip_btClusters);
  fChain->SetBranchAddress("pip_btMaxima", &pip_btMaxima, &b_pip_btMaxima);
  fChain->SetBranchAddress("pip_btMaximaCharge", &pip_btMaximaCharge, &b_pip_btMaximaCharge);
  fChain->SetBranchAddress("pip_btMaximaChargeShared", &pip_btMaximaChargeShared, &b_pip_btMaximaChargeShared);
  fChain->SetBranchAddress("pip_btMaximaChargeSharedFragment", &pip_btMaximaChargeSharedFragment, &b_pip_btMaximaChargeSharedFragment);
  fChain->SetBranchAddress("pip_btMaximaShared", &pip_btMaximaShared, &b_pip_btMaximaShared);
  fChain->SetBranchAddress("pip_btMaximaSharedFragment", &pip_btMaximaSharedFragment, &b_pip_btMaximaSharedFragment);
  fChain->SetBranchAddress("pip_btMeanDist", &pip_btMeanDist, &b_pip_btMeanDist);
  fChain->SetBranchAddress("pip_btNearbyMaxima", &pip_btNearbyMaxima, &b_pip_btNearbyMaxima);
  fChain->SetBranchAddress("pip_btNearbyMaximaShared", &pip_btNearbyMaximaShared, &b_pip_btNearbyMaximaShared);
  fChain->SetBranchAddress("pip_btPadsClus", &pip_btPadsClus, &b_pip_btPadsClus);
  fChain->SetBranchAddress("pip_btPadsRing", &pip_btPadsRing, &b_pip_btPadsRing);
  fChain->SetBranchAddress("pip_btRingMatrix", &pip_btRingMatrix, &b_pip_btRingMatrix);
  fChain->SetBranchAddress("pip_dedx_mdc", &pip_dedx_mdc, &b_pip_dedx_mdc);
  fChain->SetBranchAddress("pip_dedx_tof", &pip_dedx_tof, &b_pip_dedx_tof);
  fChain->SetBranchAddress("pip_id", &pip_id, &b_pip_id);
  fChain->SetBranchAddress("pip_isBT", &pip_isBT, &b_pip_isBT);
  fChain->SetBranchAddress("pip_isOffVertexClust", &pip_isOffVertexClust, &b_pip_isOffVertexClust);
  fChain->SetBranchAddress("pip_isPrimaryVertex", &pip_isPrimaryVertex, &b_pip_isPrimaryVertex);
  fChain->SetBranchAddress("pip_isUsedVertex", &pip_isUsedVertex, &b_pip_isUsedVertex);
  fChain->SetBranchAddress("pip_isring", &pip_isring, &b_pip_isring);
  fChain->SetBranchAddress("pip_isringmdc", &pip_isringmdc, &b_pip_isringmdc);
  fChain->SetBranchAddress("pip_isringnomatch", &pip_isringnomatch, &b_pip_isringnomatch);
  fChain->SetBranchAddress("pip_isringtrack", &pip_isringtrack, &b_pip_isringtrack);
  fChain->SetBranchAddress("pip_kIsLepton", &pip_kIsLepton, &b_pip_kIsLepton);
  fChain->SetBranchAddress("pip_kIsUsed", &pip_kIsUsed, &b_pip_kIsUsed);
  fChain->SetBranchAddress("pip_mdcinnerchi2", &pip_mdcinnerchi2, &b_pip_mdcinnerchi2);
  fChain->SetBranchAddress("pip_mdcouterchi2", &pip_mdcouterchi2, &b_pip_mdcouterchi2);
  fChain->SetBranchAddress("pip_oa_hadr", &pip_oa_hadr, &b_pip_oa_hadr);
  fChain->SetBranchAddress("pip_oa_lept", &pip_oa_lept, &b_pip_oa_lept);
  fChain->SetBranchAddress("pip_p", &pip_p, &b_pip_p);
  fChain->SetBranchAddress("pip_p_corr_em", &pip_p_corr_em, &b_pip_p_corr_em);
  fChain->SetBranchAddress("pip_p_corr_ep", &pip_p_corr_ep, &b_pip_p_corr_ep);
  fChain->SetBranchAddress("pip_p_corr_p", &pip_p_corr_p, &b_pip_p_corr_p);
  fChain->SetBranchAddress("pip_p_corr_pim", &pip_p_corr_pim, &b_pip_p_corr_pim);
  fChain->SetBranchAddress("pip_p_corr_pip", &pip_p_corr_pip, &b_pip_p_corr_pip);
  fChain->SetBranchAddress("pip_phi", &pip_phi, &b_pip_phi);
  fChain->SetBranchAddress("pip_phi_rich", &pip_phi_rich, &b_pip_phi_rich);
  fChain->SetBranchAddress("pip_pid", &pip_pid, &b_pip_pid);
  fChain->SetBranchAddress("pip_q", &pip_q, &b_pip_q);
  fChain->SetBranchAddress("pip_r", &pip_r, &b_pip_r);
  fChain->SetBranchAddress("pip_resolution", &pip_resolution, &b_pip_resolution);
  fChain->SetBranchAddress("pip_resoultion", &pip_resoultion, &b_pip_resoultion);
  fChain->SetBranchAddress("pip_rich_amp", &pip_rich_amp, &b_pip_rich_amp);
  fChain->SetBranchAddress("pip_rich_centr", &pip_rich_centr, &b_pip_rich_centr);
  fChain->SetBranchAddress("pip_rich_houtra", &pip_rich_houtra, &b_pip_rich_houtra);
  fChain->SetBranchAddress("pip_rich_padnum", &pip_rich_padnum, &b_pip_rich_padnum);
  fChain->SetBranchAddress("pip_rich_patmat", &pip_rich_patmat, &b_pip_rich_patmat);
  fChain->SetBranchAddress("pip_rkchi2", &pip_rkchi2, &b_pip_rkchi2);
  fChain->SetBranchAddress("pip_sector", &pip_sector, &b_pip_sector);
  fChain->SetBranchAddress("pip_shw_sum0", &pip_shw_sum0, &b_pip_shw_sum0);
  fChain->SetBranchAddress("pip_shw_sum1", &pip_shw_sum1, &b_pip_shw_sum1);
  fChain->SetBranchAddress("pip_shw_sum2", &pip_shw_sum2, &b_pip_shw_sum2);
  fChain->SetBranchAddress("pip_system", &pip_system, &b_pip_system);
  fChain->SetBranchAddress("pip_theta", &pip_theta, &b_pip_theta);
  fChain->SetBranchAddress("pip_theta_rich", &pip_theta_rich, &b_pip_theta_rich);
  fChain->SetBranchAddress("pip_tof_mom", &pip_tof_mom, &b_pip_tof_mom);
  fChain->SetBranchAddress("pip_tof_new", &pip_tof_new, &b_pip_tof_new);
  fChain->SetBranchAddress("pip_tof_rec", &pip_tof_rec, &b_pip_tof_rec);
  fChain->SetBranchAddress("pip_track_length", &pip_track_length, &b_pip_track_length);
  fChain->SetBranchAddress("pip_tracklength", &pip_tracklength, &b_pip_tracklength);
  fChain->SetBranchAddress("pip_z", &pip_z, &b_pip_z);
  Notify();
}

Bool_t PPipPim::Notify()
{
  // The Notify() function is called when a new file is opened. This
  // can be either for a new TTree in a TChain or when when a new TTree
  // is started when using PROOF. It is normally not necessary to make changes
  // to the generated code, but the routine can be extended by the
  // user if needed. The return value is currently not used.

  return kTRUE;
}

void PPipPim::Show(Long64_t entry)
{
  // Print contents of entry.
  // If entry is not specified, print current entry
  if (!fChain) return;
  fChain->Show(entry);
}
Int_t PPipPim::Cut(Long64_t entry)
{
  // This function may be called from Loop.
  // returns  1 if entry is accepted.
  // returns -1 otherwise.
  return 1;
}

