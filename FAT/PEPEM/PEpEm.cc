#include "PEpEm.h"
#include "data.h"
#include <iostream>
#include <TLorentzVector.h>
#include <TVector3.h>
#include "hntuple.h"


using namespace std;
using namespace PATData;

void PEpEm::Loop()
{
  static long licznik = 0;

  if (fChain == 0) return;
  /*
    (*tlo)["em_btChargeRing"] = 0;
    (*tlo)["em_btChargeSum"] = 0;
    (*tlo)["em_btChi2"] = 0;
    (*tlo)["em_btClusters"] = 0;
    (*tlo)["em_btMaxima"] = 0;
    (*tlo)["em_btMaximaCharge"] = 0;
    (*tlo)["em_btMaximaChargeShared"] = 0;
    (*tlo)["em_btMaximaChargeSharedFragment"] = 0;
    (*tlo)["em_btMaximaShared"] = 0;
    (*tlo)["em_btMaximaSharedFragment"] = 0;
    (*tlo)["em_btMeanDist"] = 0;
    (*tlo)["em_btNearbyMaxima"] = 0;
    (*tlo)["em_btNearbyMaximaShared"] = 0;
    (*tlo)["em_btPadsClus"] = 0;
    (*tlo)["em_btPadsRing"] = 0;
    (*tlo)["em_btRingMatrix"] = 0;
    (*tlo)["ep_btChargeRing"] = 0;
    (*tlo)["ep_btChargeSum"] = 0;
    (*tlo)["ep_btChi2"] = 0;
    (*tlo)["ep_btClusters"] = 0;
    (*tlo)["ep_btMaxima"] = 0;
    (*tlo)["ep_btMaximaCharge"] = 0;
    (*tlo)["ep_btMaximaChargeShared"] = 0;
    (*tlo)["ep_btMaximaChargeSharedFragment"] = 0;
    (*tlo)["ep_btMaximaShared"] = 0;
    (*tlo)["ep_btMaximaSharedFragment"] = 0;
    (*tlo)["ep_btMeanDist"] = 0;
    (*tlo)["ep_btNearbyMaxima"] = 0;
    (*tlo)["ep_btNearbyMaximaShared"] = 0;
    (*tlo)["ep_btPadsClus"] = 0;
    (*tlo)["ep_btPadsRing"] = 0;
    (*tlo)["ep_btRingMatrix"] = 0;


    (*tlo)["ep_mom"] = 0;
    (*tlo)["ep_theta"] = 0;
    (*tlo)["ep_theta_rich"] = 0;
    (*tlo)["ep_phi"] = 0;
    (*tlo)["ep_phi_rich"] = 0;
    (*tlo)["ep_beta"] = 0;
    (*tlo)["em_mom"] = 0;
    (*tlo)["em_theta"] = 0;
    (*tlo)["em_theta_rich"] = 0;
    (*tlo)["em_phi"] = 0;
    (*tlo)["em_phi_rich"] = 0;
    (*tlo)["em_beta"] = 0;
    (*tlo)["oa"] = 0;
    (*tlo)["oa_rich"] = 0;
    (*tlo)["sig"] = 0;
    (*tlo)["ep_m"] = 0;
    (*tlo)["em_m"] = 0;
    (*tlo)["epem_inv_mass"] = 0;
    (*tlo)["epem_inv_mass2"] = 0;
    (*tlo)["epem_miss_mass"] = 0;
    (*tlo)["epem_miss_mass2"] = 0;
    (*tlo)["epem_y"] = 0;
    (*tlo)["epem_pt"] = 0;
    (*tlo)["ep_rich_amp"] = 0;
    (*tlo)["ep_rich_centr"] = 0;
    (*tlo)["ep_rich_padnum"] = 0;
    (*tlo)["ep_rich_patmat"] = 0;
    (*tlo)["ep_rich_houtra"] = 0;
    (*tlo)["em_rich_amp"] = 0;
    (*tlo)["em_rich_centr"] = 0;
    (*tlo)["em_rich_padnum"] = 0;
    (*tlo)["em_rich_patmat"] = 0;
    (*tlo)["em_rich_houtra"] = 0;

    (*tlo)["eVert_x"] = 0;
    (*tlo)["eVert_y"] = 0;
    (*tlo)["eVert_z"] = 0;

    (*tlo)["eVertReco_z"] = -1000.;
    (*tlo)["eVertReco_x"] = -1000.;
    (*tlo)["eVertReco_y"] = -1000.;
    (*tlo)["evtPileupMeta"] = 0.;
    (*tlo)["evtPileupStart"] = 0.;
    (*tlo)["ep_isOffVertexClust"] = 0.;
    (*tlo)["ep_p_corr_ep"] = 0.;
    (*tlo)["em_isOffVertexClust"] = 0.;
    (*tlo)["em_p_corr_em"] = 0.;

    tlo->fill();
  */

  Long64_t nentries = fChain->GetEntries();

  Long64_t nbytes = 0, nb = 0;
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;

    ++licznik;
    if ((licznik % 100000)==0) cout << "Events: " << licznik << endl;


    double F = 1.006;
    TVector3 v1, v2, v3, v4;
    v2.SetXYZ(F*ep_p*sin(D2R*ep_theta)*cos(D2R*ep_phi),F*ep_p*sin(D2R*ep_theta)*sin(D2R*ep_phi),F*ep_p*cos(D2R*ep_theta));
    v3.SetXYZ(F*em_p*sin(D2R*em_theta)*cos(D2R*em_phi),F*em_p*sin(D2R*em_theta)*sin(D2R*em_phi),F*em_p*cos(D2R*em_theta));
    v4.SetXYZ(F*p_p*sin(D2R*p_theta)*cos(D2R*p_phi),F*p_p*sin(D2R*p_theta)*sin(D2R*p_phi),F*p_p*cos(D2R*p_theta));

      
    TVector3 r1, r2;
    r1.SetXYZ(sin(D2R*ep_theta_rich)*cos(D2R*ep_phi_rich),sin(D2R*ep_theta_rich)*sin(D2R*ep_phi_rich),cos(D2R*ep_theta_rich));
    r2.SetXYZ(sin(D2R*em_theta_rich)*cos(D2R*em_phi_rich),sin(D2R*em_theta_rich)*sin(D2R*em_phi_rich),cos(D2R*em_theta_rich));

    e1->SetVectM( v2, 0.51099906 );
    e2->SetVectM( v3, 0.51099906 );
    p->SetVectM( v4, 938.27231);

    *gammae1e2 = *e1 + *e2;
    *e1e2 = *e1 + *e2;
    *e1_delta = *e1;
    *e2_delta = *e2;
    *e1e2_miss = *beam - *e1 - *e2;
    *pe1e2_miss= *beam - (*e1 + *e2 + *p);

    double m2_inv_e1e2 = gammae1e2->M2();
    double m_inv_e1e2 = gammae1e2->M();
    double oa = R2D * openingangle(*e1, *e2);
    double oa_rich = R2D * openingangle(r1, r2);

    double e1_mass = ep_p*ep_p * (  1. / (ep_beta_new*ep_beta_new)  - 1. ) ;
    double e2_mass = em_p*em_p * (  1. / (em_beta_new*em_beta_new)  - 1. ) ;

    //	  cout << "opening angle = " << oa << endl;

    ACC = 1.;
    EFF = 1.;


    gammae1e2->Boost(0., 0., -(beam->Beta()));
    e1_delta->Boost(0., 0., -(beam->Beta()));
    e2_delta->Boost(0., 0., -(beam->Beta()));

    e1_delta->Boost( -gammae1e2->Px()/gammae1e2->E(), -gammae1e2->Py()/gammae1e2->E(), -gammae1e2->Pz()/gammae1e2->E());
    e2_delta->Boost( -gammae1e2->Px()/gammae1e2->E(), -gammae1e2->Py()/gammae1e2->E(), -gammae1e2->Pz()/gammae1e2->E());

    //cout << "Poczatek obliczen..." << endl;

    //double ang_cut = 0.;
    double ang_cut = 9.;

    double close_cut = 4.;
    double nonfit_close_cut = -4.;
    //double close_cut = 0.;
    //double nonfit_close_cut = 0.;
    //double close_cut = 4.;


#ifdef FLANCH
    insideEmS0 = (pEmS0 == 0) ? 0 : pEmS0->IsInside(em_z,em_theta);
    insideEmS1 = (pEmS1 == 0) ? 0 : pEmS1->IsInside(em_z,em_theta);
    insideEpS0 = (pEpS0 == 0) ? 0 : pEpS0->IsInside(ep_z,ep_theta);
    insideEpS1 = (pEpS1 == 0) ? 0 : pEpS1->IsInside(ep_z,ep_theta);
    //insideEmS0 = (pEmS0 == 0) ? 0 : pEmS0->IsInside(eVert_z,em_theta);
    //insideEmS1 = (pEmS1 == 0) ? 0 : pEmS1->IsInside(eVert_z,em_theta);
    //insideEpS0 = (pEpS0 == 0) ? 0 : pEpS0->IsInside(eVert_z,ep_theta);
    //insideEpS1 = (pEpS1 == 0) ? 0 : pEpS1->IsInside(eVert_z,ep_theta);
#endif

    insideTarget = 1;

#ifdef RECTANG
    insideEmS0 = (em_theta > 50 && em_z < -50 /* && em_p<200.*/) ? 1 : 0;
    insideEmS1 = (em_theta > 50 && em_z < -50 /* && em_p<200.*/) ? 1 : 0;
    insideEpS0 = (ep_theta > 50 && ep_z < -50 /* && ep_p<200.*/) ? 1 : 0;
    insideEpS1 = (ep_theta > 50 && ep_z < -50 /* && ep_p<200.*/) ? 1 : 0;
#endif

    //#ifdef NOCUT
    insideEmS0 = 0;
    insideEmS1 = 0;
    insideEpS0 = 0;
    insideEpS1 = 0;
    //#endif


    NoLeptonE1 = !((ep_oa_lept< close_cut&&ep_oa_lept>0.0) &&ep_oa_lept>nonfit_close_cut );
    NoHadronE1 = 1;//(ep_oa_hadr< close_cut &&ep_oa_hadr>nonfit_close_cut );
    NoLeptonE2 = !((em_oa_lept< close_cut&&em_oa_lept>0.0) &&em_oa_lept>nonfit_close_cut );
    NoHadronE2 = 1;//(em_oa_hadr< close_cut &&em_oa_hadr>nonfit_close_cut );
    NoHadronE1 = 1;
    NoHadronE2 = 1;

    /*
      NoLeptonE1 = 1;
      NoHadronE1 = 1;
      NoLeptonE2 = 1;
      NoHadronE2 = 1;
    */

    Positron = (((ep_system==0&&insideEpS0==0)||(ep_system==1&&insideEpS1==0)));
    Electron = (((em_system==0&&insideEmS0==0)||(em_system==1&&insideEmS1==0)));

    ElectronPositron = Positron && NoLeptonE1 && NoHadronE1  &&  Electron && NoLeptonE2 && NoHadronE2  &&  insideTarget;


      
           
    bool bt_em_condition=(em_isBT!=-1
			  //&& em_btMaxima>=2
			  && em_btPadsRing>=2
			  );
    bool bt_ep_condition=(ep_isBT!=-1
			  //&& ep_btMaxima>=2
			  && ep_btPadsRing>=2
			  );
    bool bt_condition=(bt_em_condition && bt_ep_condition);
    bool pre_shower= (ep_system==0?(ep_shw_sum1+ep_shw_sum2-ep_shw_sum0) > (parametrization(ep_p)):true)
      &&(em_system==0?(em_shw_sum1+em_shw_sum2-em_shw_sum0) > (parametrization(em_p)):true);
    //bool flanch=!(ep_theta>65 && eVert_z<-55) && !(em_theta>65 && eVert_z<-55); 
    bool flanch = !((ep_theta>65 || em_theta >65) && p_z<-55);
    bool mass_condition=(ep_p>50 && em_p>50 && ep_p<2000. && em_p<2000.
			 //&& ep_p>200 && em_p>200
			 //&&(ep_system==0?ep_beta_new>0.95:ep_beta_new>0.92)&&(em_system==0?em_beta_new>0.95:em_beta_new>0.92)
			 //&& em_beta_new<1.1 && ep_beta_new<1.1
			 && em_beta_new>0.8 && ep_beta_new>0.8
			 && pre_shower
			 && flanch
			 && m_inv_e1e2>140
			 //&& TMath::Sqrt(p_p*p_p*(1/p_beta_new*1/p_beta_new-1))>750 //mass cut for protons
			 );
    bool proton_condition=(isBest==1)
      //&& TMath::Sqrt(p_p*p_p*(1/p_beta_new*1/p_beta_new-1))<1100
      ;
    int i_array=0;

    //proton part****************************
    if(proton_condition)
      {
	proton_p_beta->Fill(p_p*p_q,p_beta_new);
	proton_p->Fill(p_p);
	proton_E->Fill(TMath::Sqrt(p_p*p_p+938.27*938.27));
	proton_m->Fill(TMath::Sqrt(p_p*p_p*(1/p_beta_new*1/p_beta_new-1)));
      }
    //end of proton part********************* 
      
    if (ElectronPositron && /*(((int)trigbit)&16) && trigdec>0 &&*/  isBest>=0 && oa > ang_cut /*&& eVertReco_z>-500 */ && mass_condition && proton_condition) 
      {
	if((ep_isring>0 || bt_ep_condition) && (em_isring>0 || bt_em_condition))// aLL DI-LEPTONS
	  {
	    miss_mass_all->Fill(pe1e2_miss->M());
	      
	  }
	    
	if(ep_isring>0 && em_isring>0)//RF signal
	  { 
	    miss_mass_RF->Fill(pe1e2_miss->M());
	    sig_all->Fill(m_inv_e1e2);
	    sig_all_var->Fill(m_inv_e1e2);
	    momentum_spectrum->Fill(ep_p);
	    momentum_spectrum->Fill(-1*em_p);
	    q_vs_p_leptons_RF->Fill(ep_p,ep_shw_sum1+ep_shw_sum2-ep_shw_sum0);
	    q_vs_p_leptons_RF->Fill(em_p,em_shw_sum1+em_shw_sum2-em_shw_sum0);
	    ep_beta_mom->Fill(ep_beta_new,ep_p);
	    em_beta_mom->Fill(em_beta_new,em_p);
	  }
	if( (bt_ep_condition||ep_isring) && (bt_em_condition||em_isring) && !(em_isring && ep_isring))//BT profit
	  {
	    miss_mass_profit->Fill(pe1e2_miss->M());
	  }
	if(bt_condition)//bt signal  
	  {
	    miss_mass_BT->Fill(pe1e2_miss->M());
	    sig_all_bt->Fill(m_inv_e1e2);
	    sig_all_var_bt->Fill(m_inv_e1e2);
	    momentum_spectrum_bt->Fill(ep_p);
	    momentum_spectrum_bt->Fill(-1*em_p);
	    q_vs_p_leptons_BT->Fill(ep_p,ep_shw_sum1+ep_shw_sum2-ep_shw_sum0);
	    q_vs_p_leptons_BT->Fill(em_p,em_shw_sum1+em_shw_sum2-em_shw_sum0);
	    ep_beta_mom_bt->Fill(ep_beta_new,ep_p);
	    em_beta_mom_bt->Fill(em_beta_new,em_p);
	  }
	      
      }

    //tlo->fill();

    // if (Cut(ientry) < 0) continue;
  } // end of main loop
} // eof Loop 



// -------------------------------------------------------------------------------------------------
PEpEm::PEpEm(TTree *tree) 
{
  // if parameter tree is not specified (or zero), connect the file
  // used to generate this class and read the Tree.

  em_acc = em_acc_err = ep_acc = ep_acc_err = 0.;
  em_eff = em_eff_err = ep_eff = ep_eff_err = 0.;

  if (tree == 0) {
	  
    TChain * chain = new TChain("EpEm_ID","");

    // -- NOWA EPOKA ANALIZY ----------------------------------------
    //chain->Add("/lustre/nyx/hades/user/knowakow/PNB/PAT/FILES/288/day288.root/EpEm_ID");
    //chain->Add("/lustre/nyx/hades/user/knowakow/PNB/PAT/FILES/288_new/lepton288new.root/EpEm_ID");
    //chain->Add("/lustre/nyx/hades/user/knowakow/PNB/PAT/FILES/288_new_pcut/lepton288new_p.root/EpEm_ID");
    //chain->Add("/lustre/nyx/hades/user/knowakow/PNB/PAT/FILES/sep08_all/list5/sum5.root/EpEm_ID");
    //chain->Add("/lustre/nyx/hades/user/knowakow/PNB/PAT/FILES/sep08_all/list4/sum4.root/EpEm_ID");
    //chain->Add("/lustre/nyx/hades/user/knowakow/PNB/PAT/FILES/sep08_all/list3/sum3.root/EpEm_ID");
    //chain->Add("/lustre/nyx/hades/user/knowakow/PNB/PAT/FILES/sep08_all/list2/sum2.root/EpEm_ID");
    //chain->Add("/lustre/nyx/hades/user/knowakow/PNB/PAT/FILES/sep08_all/list1/sum1.root/EpEm_ID");
    //chain->Add("/lustre/nyx/hades/user/knowakow/PP/PAT/FILES/115_combinatorics/lepton.root/EpEm_ID");
    //chain->Add("/lustre/nyx/hades/user/knowakow/PP/PAT/FILES/full_stat/leptons.root/EpEm_ID");
      
    chain->Add("/lustre/nyx/hades/user/knowakow/PP/PAT_1/FILES/bt/lepton00.root/PEpEm_ID");
    chain->Add("/lustre/nyx/hades/user/knowakow/PP/PAT_1/FILES/bt/lepton01.root/PEpEm_ID");
    chain->Add("/lustre/nyx/hades/user/knowakow/PP/PAT_1/FILES/bt/lepton02.root/PEpEm_ID");
    chain->Add("/lustre/nyx/hades/user/knowakow/PP/PAT_1/FILES/bt/lepton03.root/PEpEm_ID");
    chain->Add("/lustre/nyx/hades/user/knowakow/PP/PAT_1/FILES/bt/lepton04.root/PEpEm_ID");
    chain->Add("/lustre/nyx/hades/user/knowakow/PP/PAT_1/FILES/bt/lepton05.root/PEpEm_ID");
    chain->Add("/lustre/nyx/hades/user/knowakow/PP/PAT_1/FILES/bt/lepton06.root/PEpEm_ID");
    chain->Add("/lustre/nyx/hades/user/knowakow/PP/PAT_1/FILES/bt/lepton07.root/PEpEm_ID");
    chain->Add("/lustre/nyx/hades/user/knowakow/PP/PAT_1/FILES/bt/lepton08.root/PEpEm_ID");
    chain->Add("/lustre/nyx/hades/user/knowakow/PP/PAT_1/FILES/bt/lepton09.root/PEpEm_ID");
    chain->Add("/lustre/nyx/hades/user/knowakow/PP/PAT_1/FILES/bt/lepton10.root/PEpEm_ID");
    chain->Add("/lustre/nyx/hades/user/knowakow/PP/PAT_1/FILES/bt/lepton11.root/PEpEm_ID");
    tree = chain; 
  }

  Init(tree);
}


// -------------------------------------------------------------------------------------------------
PEpEm::~PEpEm()
{
  if (!fChain) return;
  delete fChain->GetCurrentFile();
}

// -------------------------------------------------------------------------------------------------
Int_t PEpEm::GetEntry(Long64_t entry)
{
  // Read contents of entry.
  if (!fChain) return 0;
  return fChain->GetEntry(entry);
}

// -------------------------------------------------------------------------------------------------
Long64_t PEpEm::LoadTree(Long64_t entry)
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
void PEpEm::Init(TTree *tree)
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

  fChain->SetBranchAddress("eVertClust_chi2", &eVertClust_chi2, &b_eVertClust_chi2);
  fChain->SetBranchAddress("eVertClust_x", &eVertClust_x, &b_eVertClust_x);
  fChain->SetBranchAddress("eVertClust_y", &eVertClust_y, &b_eVertClust_y);
  fChain->SetBranchAddress("eVertClust_z", &eVertClust_z, &b_eVertClust_z);
  fChain->SetBranchAddress("eVertReco_chi2", &eVertReco_chi2, &b_eVertReco_chi2);
  fChain->SetBranchAddress("eVertReco_x", &eVertReco_x, &b_eVertReco_x);
  fChain->SetBranchAddress("eVertReco_y", &eVertReco_y, &b_eVertReco_y);
  fChain->SetBranchAddress("eVertReco_z", &eVertReco_z, &b_eVertReco_z);
  fChain->SetBranchAddress("eVert_chi2", &eVert_chi2, &b_eVert_chi2);
  fChain->SetBranchAddress("eVert_x", &eVert_x, &b_eVert_x);
  fChain->SetBranchAddress("eVert_y", &eVert_y, &b_eVert_y);
  fChain->SetBranchAddress("eVert_z", &eVert_z, &b_eVert_z);
  fChain->SetBranchAddress("em_beta", &em_beta, &b_em_beta);
  fChain->SetBranchAddress("em_beta_new", &em_beta_new, &b_em_beta_new);
  fChain->SetBranchAddress("em_btChargeRing", &em_btChargeRing, &b_em_btChargeRing);
  fChain->SetBranchAddress("em_btChargeSum", &em_btChargeSum, &b_em_btChargeSum);
  fChain->SetBranchAddress("em_btChi2", &em_btChi2, &b_em_btChi2);
  fChain->SetBranchAddress("em_btClusters", &em_btClusters, &b_em_btClusters);
  fChain->SetBranchAddress("em_btMaxima", &em_btMaxima, &b_em_btMaxima);
  fChain->SetBranchAddress("em_btMaximaCharge", &em_btMaximaCharge, &b_em_btMaximaCharge);
  fChain->SetBranchAddress("em_btMaximaChargeShared", &em_btMaximaChargeShared, &b_em_btMaximaChargeShared);
  fChain->SetBranchAddress("em_btMaximaChargeSharedFragment", &em_btMaximaChargeSharedFragment, &b_em_btMaximaChargeSharedFragment);
  fChain->SetBranchAddress("em_btMaximaShared", &em_btMaximaShared, &b_em_btMaximaShared);
  fChain->SetBranchAddress("em_btMaximaSharedFragment", &em_btMaximaSharedFragment, &b_em_btMaximaSharedFragment);
  fChain->SetBranchAddress("em_btMeanDist", &em_btMeanDist, &b_em_btMeanDist);
  fChain->SetBranchAddress("em_btNearbyMaxima", &em_btNearbyMaxima, &b_em_btNearbyMaxima);
  fChain->SetBranchAddress("em_btNearbyMaximaShared", &em_btNearbyMaximaShared, &b_em_btNearbyMaximaShared);
  fChain->SetBranchAddress("em_btPadsClus", &em_btPadsClus, &b_em_btPadsClus);
  fChain->SetBranchAddress("em_btPadsRing", &em_btPadsRing, &b_em_btPadsRing);
  fChain->SetBranchAddress("em_btRingMatrix", &em_btRingMatrix, &b_em_btRingMatrix);
  fChain->SetBranchAddress("em_dedx_mdc", &em_dedx_mdc, &b_em_dedx_mdc);
  fChain->SetBranchAddress("em_dedx_tof", &em_dedx_tof, &b_em_dedx_tof);
  fChain->SetBranchAddress("em_id", &em_id, &b_em_id);
  fChain->SetBranchAddress("em_isBT", &em_isBT, &b_em_isBT);
  fChain->SetBranchAddress("em_isOffVertexClust", &em_isOffVertexClust, &b_em_isOffVertexClust);
  fChain->SetBranchAddress("em_isPrimaryVertex", &em_isPrimaryVertex, &b_em_isPrimaryVertex);
  fChain->SetBranchAddress("em_isUsedVertex", &em_isUsedVertex, &b_em_isUsedVertex);
  fChain->SetBranchAddress("em_isring", &em_isring, &b_em_isring);
  fChain->SetBranchAddress("em_isringmdc", &em_isringmdc, &b_em_isringmdc);
  fChain->SetBranchAddress("em_isringnomatch", &em_isringnomatch, &b_em_isringnomatch);
  fChain->SetBranchAddress("em_isringtrack", &em_isringtrack, &b_em_isringtrack);
  fChain->SetBranchAddress("em_kIsLepton", &em_kIsLepton, &b_em_kIsLepton);
  fChain->SetBranchAddress("em_kIsUsed", &em_kIsUsed, &b_em_kIsUsed);
  fChain->SetBranchAddress("em_mdcinnerchi2", &em_mdcinnerchi2, &b_em_mdcinnerchi2);
  fChain->SetBranchAddress("em_mdcouterchi2", &em_mdcouterchi2, &b_em_mdcouterchi2);
  fChain->SetBranchAddress("em_oa_hadr", &em_oa_hadr, &b_em_oa_hadr);
  fChain->SetBranchAddress("em_oa_lept", &em_oa_lept, &b_em_oa_lept);
  fChain->SetBranchAddress("em_p", &em_p, &b_em_p);
  fChain->SetBranchAddress("em_p_corr_em", &em_p_corr_em, &b_em_p_corr_em);
  fChain->SetBranchAddress("em_p_corr_ep", &em_p_corr_ep, &b_em_p_corr_ep);
  fChain->SetBranchAddress("em_p_corr_p", &em_p_corr_p, &b_em_p_corr_p);
  fChain->SetBranchAddress("em_p_corr_pim", &em_p_corr_pim, &b_em_p_corr_pim);
  fChain->SetBranchAddress("em_p_corr_pip", &em_p_corr_pip, &b_em_p_corr_pip);
  fChain->SetBranchAddress("em_phi", &em_phi, &b_em_phi);
  fChain->SetBranchAddress("em_phi_rich", &em_phi_rich, &b_em_phi_rich);
  fChain->SetBranchAddress("em_pid", &em_pid, &b_em_pid);
  fChain->SetBranchAddress("em_q", &em_q, &b_em_q);
  fChain->SetBranchAddress("em_r", &em_r, &b_em_r);
  fChain->SetBranchAddress("em_resolution", &em_resolution, &b_em_resolution);
  fChain->SetBranchAddress("em_resoultion", &em_resoultion, &b_em_resoultion);
  fChain->SetBranchAddress("em_rich_amp", &em_rich_amp, &b_em_rich_amp);
  fChain->SetBranchAddress("em_rich_centr", &em_rich_centr, &b_em_rich_centr);
  fChain->SetBranchAddress("em_rich_houtra", &em_rich_houtra, &b_em_rich_houtra);
  fChain->SetBranchAddress("em_rich_padnum", &em_rich_padnum, &b_em_rich_padnum);
  fChain->SetBranchAddress("em_rich_patmat", &em_rich_patmat, &b_em_rich_patmat);
  fChain->SetBranchAddress("em_rkchi2", &em_rkchi2, &b_em_rkchi2);
  fChain->SetBranchAddress("em_sector", &em_sector, &b_em_sector);
  fChain->SetBranchAddress("em_shw_sum0", &em_shw_sum0, &b_em_shw_sum0);
  fChain->SetBranchAddress("em_shw_sum1", &em_shw_sum1, &b_em_shw_sum1);
  fChain->SetBranchAddress("em_shw_sum2", &em_shw_sum2, &b_em_shw_sum2);
  fChain->SetBranchAddress("em_system", &em_system, &b_em_system);
  fChain->SetBranchAddress("em_theta", &em_theta, &b_em_theta);
  fChain->SetBranchAddress("em_theta_rich", &em_theta_rich, &b_em_theta_rich);
  fChain->SetBranchAddress("em_tof_mom", &em_tof_mom, &b_em_tof_mom);
  fChain->SetBranchAddress("em_tof_new", &em_tof_new, &b_em_tof_new);
  fChain->SetBranchAddress("em_tof_rec", &em_tof_rec, &b_em_tof_rec);
  fChain->SetBranchAddress("em_track_length", &em_track_length, &b_em_track_length);
  fChain->SetBranchAddress("em_tracklength", &em_tracklength, &b_em_tracklength);
  fChain->SetBranchAddress("em_z", &em_z, &b_em_z);
  fChain->SetBranchAddress("ep_beta", &ep_beta, &b_ep_beta);
  fChain->SetBranchAddress("ep_beta_new", &ep_beta_new, &b_ep_beta_new);
  fChain->SetBranchAddress("ep_btChargeRing", &ep_btChargeRing, &b_ep_btChargeRing);
  fChain->SetBranchAddress("ep_btChargeSum", &ep_btChargeSum, &b_ep_btChargeSum);
  fChain->SetBranchAddress("ep_btChi2", &ep_btChi2, &b_ep_btChi2);
  fChain->SetBranchAddress("ep_btClusters", &ep_btClusters, &b_ep_btClusters);
  fChain->SetBranchAddress("ep_btMaxima", &ep_btMaxima, &b_ep_btMaxima);
  fChain->SetBranchAddress("ep_btMaximaCharge", &ep_btMaximaCharge, &b_ep_btMaximaCharge);
  fChain->SetBranchAddress("ep_btMaximaChargeShared", &ep_btMaximaChargeShared, &b_ep_btMaximaChargeShared);
  fChain->SetBranchAddress("ep_btMaximaChargeSharedFragment", &ep_btMaximaChargeSharedFragment, &b_ep_btMaximaChargeSharedFragment);
  fChain->SetBranchAddress("ep_btMaximaShared", &ep_btMaximaShared, &b_ep_btMaximaShared);
  fChain->SetBranchAddress("ep_btMaximaSharedFragment", &ep_btMaximaSharedFragment, &b_ep_btMaximaSharedFragment);
  fChain->SetBranchAddress("ep_btMeanDist", &ep_btMeanDist, &b_ep_btMeanDist);
  fChain->SetBranchAddress("ep_btNearbyMaxima", &ep_btNearbyMaxima, &b_ep_btNearbyMaxima);
  fChain->SetBranchAddress("ep_btNearbyMaximaShared", &ep_btNearbyMaximaShared, &b_ep_btNearbyMaximaShared);
  fChain->SetBranchAddress("ep_btPadsClus", &ep_btPadsClus, &b_ep_btPadsClus);
  fChain->SetBranchAddress("ep_btPadsRing", &ep_btPadsRing, &b_ep_btPadsRing);
  fChain->SetBranchAddress("ep_btRingMatrix", &ep_btRingMatrix, &b_ep_btRingMatrix);
  fChain->SetBranchAddress("ep_dedx_mdc", &ep_dedx_mdc, &b_ep_dedx_mdc);
  fChain->SetBranchAddress("ep_dedx_tof", &ep_dedx_tof, &b_ep_dedx_tof);
  fChain->SetBranchAddress("ep_id", &ep_id, &b_ep_id);
  fChain->SetBranchAddress("ep_isBT", &ep_isBT, &b_ep_isBT);
  fChain->SetBranchAddress("ep_isOffVertexClust", &ep_isOffVertexClust, &b_ep_isOffVertexClust);
  fChain->SetBranchAddress("ep_isPrimaryVertex", &ep_isPrimaryVertex, &b_ep_isPrimaryVertex);
  fChain->SetBranchAddress("ep_isUsedVertex", &ep_isUsedVertex, &b_ep_isUsedVertex);
  fChain->SetBranchAddress("ep_isring", &ep_isring, &b_ep_isring);
  fChain->SetBranchAddress("ep_isringmdc", &ep_isringmdc, &b_ep_isringmdc);
  fChain->SetBranchAddress("ep_isringnomatch", &ep_isringnomatch, &b_ep_isringnomatch);
  fChain->SetBranchAddress("ep_isringtrack", &ep_isringtrack, &b_ep_isringtrack);
  fChain->SetBranchAddress("ep_kIsLepton", &ep_kIsLepton, &b_ep_kIsLepton);
  fChain->SetBranchAddress("ep_kIsUsed", &ep_kIsUsed, &b_ep_kIsUsed);
  fChain->SetBranchAddress("ep_mdcinnerchi2", &ep_mdcinnerchi2, &b_ep_mdcinnerchi2);
  fChain->SetBranchAddress("ep_mdcouterchi2", &ep_mdcouterchi2, &b_ep_mdcouterchi2);
  fChain->SetBranchAddress("ep_oa_hadr", &ep_oa_hadr, &b_ep_oa_hadr);
  fChain->SetBranchAddress("ep_oa_lept", &ep_oa_lept, &b_ep_oa_lept);
  fChain->SetBranchAddress("ep_p", &ep_p, &b_ep_p);
  fChain->SetBranchAddress("ep_p_corr_em", &ep_p_corr_em, &b_ep_p_corr_em);
  fChain->SetBranchAddress("ep_p_corr_ep", &ep_p_corr_ep, &b_ep_p_corr_ep);
  fChain->SetBranchAddress("ep_p_corr_p", &ep_p_corr_p, &b_ep_p_corr_p);
  fChain->SetBranchAddress("ep_p_corr_pim", &ep_p_corr_pim, &b_ep_p_corr_pim);
  fChain->SetBranchAddress("ep_p_corr_pip", &ep_p_corr_pip, &b_ep_p_corr_pip);
  fChain->SetBranchAddress("ep_phi", &ep_phi, &b_ep_phi);
  fChain->SetBranchAddress("ep_phi_rich", &ep_phi_rich, &b_ep_phi_rich);
  fChain->SetBranchAddress("ep_pid", &ep_pid, &b_ep_pid);
  fChain->SetBranchAddress("ep_q", &ep_q, &b_ep_q);
  fChain->SetBranchAddress("ep_r", &ep_r, &b_ep_r);
  fChain->SetBranchAddress("ep_resolution", &ep_resolution, &b_ep_resolution);
  fChain->SetBranchAddress("ep_resoultion", &ep_resoultion, &b_ep_resoultion);
  fChain->SetBranchAddress("ep_rich_amp", &ep_rich_amp, &b_ep_rich_amp);
  fChain->SetBranchAddress("ep_rich_centr", &ep_rich_centr, &b_ep_rich_centr);
  fChain->SetBranchAddress("ep_rich_houtra", &ep_rich_houtra, &b_ep_rich_houtra);
  fChain->SetBranchAddress("ep_rich_padnum", &ep_rich_padnum, &b_ep_rich_padnum);
  fChain->SetBranchAddress("ep_rich_patmat", &ep_rich_patmat, &b_ep_rich_patmat);
  fChain->SetBranchAddress("ep_rkchi2", &ep_rkchi2, &b_ep_rkchi2);
  fChain->SetBranchAddress("ep_sector", &ep_sector, &b_ep_sector);
  fChain->SetBranchAddress("ep_shw_sum0", &ep_shw_sum0, &b_ep_shw_sum0);
  fChain->SetBranchAddress("ep_shw_sum1", &ep_shw_sum1, &b_ep_shw_sum1);
  fChain->SetBranchAddress("ep_shw_sum2", &ep_shw_sum2, &b_ep_shw_sum2);
  fChain->SetBranchAddress("ep_system", &ep_system, &b_ep_system);
  fChain->SetBranchAddress("ep_theta", &ep_theta, &b_ep_theta);
  fChain->SetBranchAddress("ep_theta_rich", &ep_theta_rich, &b_ep_theta_rich);
  fChain->SetBranchAddress("ep_tof_mom", &ep_tof_mom, &b_ep_tof_mom);
  fChain->SetBranchAddress("ep_tof_new", &ep_tof_new, &b_ep_tof_new);
  fChain->SetBranchAddress("ep_tof_rec", &ep_tof_rec, &b_ep_tof_rec);
  fChain->SetBranchAddress("ep_track_length", &ep_track_length, &b_ep_track_length);
  fChain->SetBranchAddress("ep_tracklength", &ep_tracklength, &b_ep_tracklength);
  fChain->SetBranchAddress("ep_z", &ep_z, &b_ep_z);
  fChain->SetBranchAddress("event", &event, &b_event);
  fChain->SetBranchAddress("fw_beta_1", &fw_beta_1, &b_fw_beta_1);
  fChain->SetBranchAddress("fw_beta_2", &fw_beta_2, &b_fw_beta_2);
  fChain->SetBranchAddress("fw_beta_3", &fw_beta_3, &b_fw_beta_3);
  fChain->SetBranchAddress("fw_charge_1", &fw_charge_1, &b_fw_charge_1);
  fChain->SetBranchAddress("fw_charge_2", &fw_charge_2, &b_fw_charge_2);
  fChain->SetBranchAddress("fw_charge_3", &fw_charge_3, &b_fw_charge_3);
  fChain->SetBranchAddress("fw_cluster_mult", &fw_cluster_mult, &b_fw_cluster_mult);
  fChain->SetBranchAddress("fw_distance_1", &fw_distance_1, &b_fw_distance_1);
  fChain->SetBranchAddress("fw_distance_2", &fw_distance_2, &b_fw_distance_2);
  fChain->SetBranchAddress("fw_distance_3", &fw_distance_3, &b_fw_distance_3);
  fChain->SetBranchAddress("fw_mult", &fw_mult, &b_fw_mult);
  fChain->SetBranchAddress("fw_p_1", &fw_p_1, &b_fw_p_1);
  fChain->SetBranchAddress("fw_p_2", &fw_p_2, &b_fw_p_2);
  fChain->SetBranchAddress("fw_p_3", &fw_p_3, &b_fw_p_3);
  fChain->SetBranchAddress("fw_phi_1", &fw_phi_1, &b_fw_phi_1);
  fChain->SetBranchAddress("fw_phi_2", &fw_phi_2, &b_fw_phi_2);
  fChain->SetBranchAddress("fw_phi_3", &fw_phi_3, &b_fw_phi_3);
  fChain->SetBranchAddress("fw_size_1", &fw_size_1, &b_fw_size_1);
  fChain->SetBranchAddress("fw_size_2", &fw_size_2, &b_fw_size_2);
  fChain->SetBranchAddress("fw_size_3", &fw_size_3, &b_fw_size_3);
  fChain->SetBranchAddress("fw_spectator_1", &fw_spectator_1, &b_fw_spectator_1);
  fChain->SetBranchAddress("fw_spectator_2", &fw_spectator_2, &b_fw_spectator_2);
  fChain->SetBranchAddress("fw_spectator_3", &fw_spectator_3, &b_fw_spectator_3);
  fChain->SetBranchAddress("fw_theta_1", &fw_theta_1, &b_fw_theta_1);
  fChain->SetBranchAddress("fw_theta_2", &fw_theta_2, &b_fw_theta_2);
  fChain->SetBranchAddress("fw_theta_3", &fw_theta_3, &b_fw_theta_3);
  fChain->SetBranchAddress("fw_time_1", &fw_time_1, &b_fw_time_1);
  fChain->SetBranchAddress("fw_time_2", &fw_time_2, &b_fw_time_2);
  fChain->SetBranchAddress("fw_time_3", &fw_time_3, &b_fw_time_3);
  fChain->SetBranchAddress("fw_time_min_1", &fw_time_min_1, &b_fw_time_min_1);
  fChain->SetBranchAddress("fw_time_min_2", &fw_time_min_2, &b_fw_time_min_2);
  fChain->SetBranchAddress("fw_time_min_3", &fw_time_min_3, &b_fw_time_min_3);
  fChain->SetBranchAddress("fw_x_lab_1", &fw_x_lab_1, &b_fw_x_lab_1);
  fChain->SetBranchAddress("fw_x_lab_2", &fw_x_lab_2, &b_fw_x_lab_2);
  fChain->SetBranchAddress("fw_x_lab_3", &fw_x_lab_3, &b_fw_x_lab_3);
  fChain->SetBranchAddress("fw_y_lab_1", &fw_y_lab_1, &b_fw_y_lab_1);
  fChain->SetBranchAddress("fw_y_lab_2", &fw_y_lab_2, &b_fw_y_lab_2);
  fChain->SetBranchAddress("fw_y_lab_3", &fw_y_lab_3, &b_fw_y_lab_3);
  fChain->SetBranchAddress("fw_z_lab_1", &fw_z_lab_1, &b_fw_z_lab_1);
  fChain->SetBranchAddress("fw_z_lab_2", &fw_z_lab_2, &b_fw_z_lab_2);
  fChain->SetBranchAddress("fw_z_lab_3", &fw_z_lab_3, &b_fw_z_lab_3);
  fChain->SetBranchAddress("hneg_mult", &hneg_mult, &b_hneg_mult);
  fChain->SetBranchAddress("hpos_mult", &hpos_mult, &b_hpos_mult);
  fChain->SetBranchAddress("isBest", &isBest, &b_isBest);
  fChain->SetBranchAddress("lneg_mult", &lneg_mult, &b_lneg_mult);
  fChain->SetBranchAddress("lpos_mult", &lpos_mult, &b_lpos_mult);
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
  fChain->SetBranchAddress("runnumber", &runnumber, &b_runnumber);
  fChain->SetBranchAddress("totalmult", &totalmult, &b_totalmult);
  fChain->SetBranchAddress("trigbit", &trigbit, &b_trigbit);
  fChain->SetBranchAddress("trigdec", &trigdec, &b_trigdec);
  fChain->SetBranchAddress("trigdownscale", &trigdownscale, &b_trigdownscale);
  fChain->SetBranchAddress("trigdownscaleflag", &trigdownscaleflag, &b_trigdownscaleflag);
  Notify();
}

// -------------------------------------------------------------------------------------------------
Bool_t PEpEm::Notify()
{
  // The Notify() function is called when a new file is opened. This
  // can be either for a new TTree in a TChain or when when a new TTree
  // is started when using PROOF. It is normally not necessary to make changes
  // to the generated code, but the routine can be extended by the
  // user if needed. The return value is currently not used.

  return kTRUE;
}

// -------------------------------------------------------------------------------------------------
void PEpEm::Show(Long64_t entry)
{
  // Print contents of entry.
  // If entry is not specified, print current entry
  if (!fChain) return;
  fChain->Show(entry);
}

// -------------------------------------------------------------------------------------------------
Int_t PEpEm::Cut(Long64_t entry)
{
  // This function may be called from Loop.
  // returns  1 if entry is accepted.
  // returns -1 otherwise.
  return 1;
}
