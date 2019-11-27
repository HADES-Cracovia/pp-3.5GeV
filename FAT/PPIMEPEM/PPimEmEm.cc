#include "PPimEmEm.h"
#include "data.h"
#include <iostream>
#include <TLorentzVector.h>
#include <TVector3.h>
#include "hntuple.h"


using namespace std;
using namespace PATData;

void PPimEmEm::Loop()
{
    static long licznik = 0;

   if (fChain == 0) return;
/*
        (*tlo)["em2_btChargeRing"] = 0;
        (*tlo)["em2_btChargeSum"] = 0;
        (*tlo)["em2_btChi2"] = 0;
        (*tlo)["em2_btClusters"] = 0;
        (*tlo)["em2_btMaxima"] = 0;
        (*tlo)["em2_btMaximaCharge"] = 0;
        (*tlo)["em2_btMaximaChargeShared"] = 0;
        (*tlo)["em2_btMaximaChargeSharedFragment"] = 0;
        (*tlo)["em2_btMaximaShared"] = 0;
        (*tlo)["em2_btMaximaSharedFragment"] = 0;
        (*tlo)["em2_btMeanDist"] = 0;
        (*tlo)["em2_btNearbyMaxima"] = 0;
        (*tlo)["em2_btNearbyMaximaShared"] = 0;
        (*tlo)["em2_btPadsClus"] = 0;
        (*tlo)["em2_btPadsRing"] = 0;
        (*tlo)["em2_btRingMatrix"] = 0;
        (*tlo)["em1_btChargeRing"] = 0;
        (*tlo)["em1_btChargeSum"] = 0;
        (*tlo)["em1_btChi2"] = 0;
        (*tlo)["em1_btClusters"] = 0;
        (*tlo)["em1_btMaxima"] = 0;
        (*tlo)["em1_btMaximaCharge"] = 0;
        (*tlo)["em1_btMaximaChargeShared"] = 0;
        (*tlo)["em1_btMaximaChargeSharedFragment"] = 0;
        (*tlo)["em1_btMaximaShared"] = 0;
        (*tlo)["em1_btMaximaSharedFragment"] = 0;
        (*tlo)["em1_btMeanDist"] = 0;
        (*tlo)["em1_btNearbyMaxima"] = 0;
        (*tlo)["em1_btNearbyMaximaShared"] = 0;
        (*tlo)["em1_btPadsClus"] = 0;
        (*tlo)["em1_btPadsRing"] = 0;
        (*tlo)["em1_btRingMatrix"] = 0;


        (*tlo)["em1_mom"] = 0;
        (*tlo)["em1_theta"] = 0;
        (*tlo)["em1_theta_rich"] = 0;
        (*tlo)["em1_phi"] = 0;
        (*tlo)["em1_phi_rich"] = 0;
        (*tlo)["em1_beta"] = 0;
        (*tlo)["em2_mom"] = 0;
        (*tlo)["em2_theta"] = 0;
        (*tlo)["em2_theta_rich"] = 0;
        (*tlo)["em2_phi"] = 0;
        (*tlo)["em2_phi_rich"] = 0;
        (*tlo)["em2_beta"] = 0;
        (*tlo)["oa"] = 0;
        (*tlo)["oa_rich"] = 0;
        (*tlo)["sig"] = 0;
        (*tlo)["em1_m"] = 0;
        (*tlo)["em2_m"] = 0;
        (*tlo)["em1em2_inv_mass"] = 0;
        (*tlo)["em1em2_inv_mass2"] = 0;
        (*tlo)["em1em2_miss_mass"] = 0;
        (*tlo)["em1em2_miss_mass2"] = 0;
        (*tlo)["em1em2_y"] = 0;
        (*tlo)["em1em2_pt"] = 0;
        (*tlo)["em1_rich_amp"] = 0;
        (*tlo)["em1_rich_centr"] = 0;
        (*tlo)["em1_rich_padnum"] = 0;
        (*tlo)["em1_rich_patmat"] = 0;
        (*tlo)["em1_rich_houtra"] = 0;
        (*tlo)["em2_rich_amp"] = 0;
        (*tlo)["em2_rich_centr"] = 0;
        (*tlo)["em2_rich_padnum"] = 0;
        (*tlo)["em2_rich_patmat"] = 0;
        (*tlo)["em2_rich_houtra"] = 0;

	    (*tlo)["eVert_x"] = 0;
	    (*tlo)["eVert_y"] = 0;
	    (*tlo)["eVert_z"] = 0;

        (*tlo)["eVertReco_z"] = -1000.;
        (*tlo)["eVertReco_x"] = -1000.;
        (*tlo)["eVertReco_y"] = -1000.;
        (*tlo)["evtPileupMeta"] = 0.;
        (*tlo)["evtPileupStart"] = 0.;
        (*tlo)["em1_isOffVertexClust"] = 0.;
        (*tlo)["em1_p_corr_em1"] = 0.;
        (*tlo)["em2_isOffVertexClust"] = 0.;
        (*tlo)["em2_p_corr_em2"] = 0.;

        tlo->fill();
*/

   Long64_t nentries = fChain->GetEntries();

   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;

      ++licznik;
	  if ((licznik % 10000)==0) cout << "Events: " << licznik << endl;


      double F = 1.006;
      TVector3 v1, v2, v3, v4;
      v2.SetXYZ(F*em1_p*sin(D2R*em1_theta)*cos(D2R*em1_phi),F*em1_p*sin(D2R*em1_theta)*sin(D2R*em1_phi),F*em1_p*cos(D2R*em1_theta));
      v3.SetXYZ(F*em2_p*sin(D2R*em2_theta)*cos(D2R*em2_phi),F*em2_p*sin(D2R*em2_theta)*sin(D2R*em2_phi),F*em2_p*cos(D2R*em2_theta));
      v4.SetXYZ(F*p_p*sin(D2R*p_theta)*cos(D2R*p_phi),F*p_p*sin(D2R*p_theta)*sin(D2R*p_phi),F*p_p*cos(D2R*p_theta));

      
      TVector3 r1, r2;
      r1.SetXYZ(sin(D2R*em1_theta_rich)*cos(D2R*em1_phi_rich),sin(D2R*em1_theta_rich)*sin(D2R*em1_phi_rich),cos(D2R*em1_theta_rich));
      r2.SetXYZ(sin(D2R*em2_theta_rich)*cos(D2R*em2_phi_rich),sin(D2R*em2_theta_rich)*sin(D2R*em2_phi_rich),cos(D2R*em2_theta_rich));

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

      double e1_mass = em1_p*em1_p * (  1. / (em1_beta_new*em1_beta_new)  - 1. ) ;
      double e2_mass = em2_p*em2_p * (  1. / (em2_beta_new*em2_beta_new)  - 1. ) ;

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
      insideEm2S0 = (pEm2S0 == 0) ? 0 : pEm2S0->IsInside(em2_z,em2_theta);
      insideEm2S1 = (pEm2S1 == 0) ? 0 : pEm2S1->IsInside(em2_z,em2_theta);
      insideEm1S0 = (pEm1S0 == 0) ? 0 : pEm1S0->IsInside(em1_z,em1_theta);
      insideEm1S1 = (pEm1S1 == 0) ? 0 : pEm1S1->IsInside(em1_z,em1_theta);
      //insideEm2S0 = (pEm2S0 == 0) ? 0 : pEm2S0->IsInside(eVert_z,em2_theta);
      //insideEm2S1 = (pEm2S1 == 0) ? 0 : pEm2S1->IsInside(eVert_z,em2_theta);
      //insideEm1S0 = (pEm1S0 == 0) ? 0 : pEm1S0->IsInside(eVert_z,em1_theta);
      //insideEm1S1 = (pEm1S1 == 0) ? 0 : pEm1S1->IsInside(eVert_z,em1_theta);
#endif

      insideTarget = 1;

#ifdef RECTANG
      insideEm2S0 = (em2_theta > 50 && em2_z < -50 /* && em2_p<200.*/) ? 1 : 0;
      insideEm2S1 = (em2_theta > 50 && em2_z < -50 /* && em2_p<200.*/) ? 1 : 0;
      insideEm1S0 = (em1_theta > 50 && em1_z < -50 /* && em1_p<200.*/) ? 1 : 0;
      insideEm1S1 = (em1_theta > 50 && em1_z < -50 /* && em1_p<200.*/) ? 1 : 0;
#endif

//#ifdef NOCUT
      insideEm2S0 = 0;
      insideEm2S1 = 0;
      insideEm1S0 = 0;
      insideEm1S1 = 0;
//#endif


      NoLeptonE1 = !((em1_oa_lept< close_cut&&em1_oa_lept>0.0) &&em1_oa_lept>nonfit_close_cut );
      NoHadronE1 = 1;//(em1_oa_hadr< close_cut &&em1_oa_hadr>nonfit_close_cut );
      NoLeptonE2 = !((em2_oa_lept< close_cut&&em2_oa_lept>0.0) &&em2_oa_lept>nonfit_close_cut );
      NoHadronE2 = 1;//(em2_oa_hadr< close_cut &&em2_oa_hadr>nonfit_close_cut );
      NoHadronE1 = 1;
      NoHadronE2 = 1;

/*
      NoLem1tonE1 = 1;
      NoHadronE1 = 1;
      NoLem1tonE2 = 1;
      NoHadronE2 = 1;
*/

      Electron1 = (((em1_system==0&&insideEm1S0==0)||(em1_system==1&&insideEm1S1==0)));
      Electron2 = (((em2_system==0&&insideEm2S0==0)||(em2_system==1&&insideEm2S1==0)));

      ElectronElectron = Electron1 && NoLeptonE1 && NoHadronE1  &&  Electron2 && NoLeptonE2 && NoHadronE2 && insideTarget &&
                         ( e1_mass < 5000. && e2_mass < 5000. );

           
      bool bt_em2_condition=(em2_isBT!=-1
			    //&& em2_btMaxima>=2
			    && em2_btPadsRing>=2
			    );
      bool bt_em1_condition=(em1_isBT!=-1
			    //&& em1_btMaxima>=2
			    && em1_btPadsRing>=2
			    );
      bool bt_condition=(bt_em2_condition && bt_em1_condition);
      bool pre_shower= (em1_system==0?(em1_shw_sum1+em1_shw_sum2-em1_shw_sum0) > (parametrization(em1_p)):true)
	             &&(em2_system==0?(em2_shw_sum1+em2_shw_sum2-em2_shw_sum0) > (parametrization(em2_p)):true);
      //bool flanch=!(em1_theta>65 && eVert_z<-55) && !(em2_theta>65 && eVert_z<-55); 
      bool flanch = !((em1_theta>65 || em2_theta >65)&& p_z<-55);
      bool mass_condition=(em1_p>50 && em2_p>50 && em1_p<2000. && em2_p<2000.
			   //&& em1_p>200 && em2_p>200
			   //&&(em1_system2==0?em1_beta_new>0.95:em1_beta_new>0.92)&&(em2_system2==0?em2_beta_new>0.95:em2_beta_new>0.92)
			   //&& em2_beta_new<1.1 && em1_beta_new<1.1
			   && em2_beta_new>0.8 && em1_beta_new>0.8
			   && pre_shower
			   && flanch
			   && m_inv_e1e2>140
			   );
      bool proton_condition=(isBest==1)
	//&& TMath::Sqrt(p_p*p_p*(1/p_beta_new*1/p_beta_new-1))<1100
			  ;
      int i_array=0;

      
      //proton part****************************
      /*
        if(proton_condition)
	{
	  proton_p_beta->Fill(p_p*p_q,p_beta_new);
	  proton_p->Fill(p_p);
	  proton_E->Fill(TMath::Sqrt(p_p*p_p+938.27*938.27));
	  proton_m->Fill(TMath::Sqrt(p_p*p_p*(1/p_beta_new*1/p_beta_new-1)));
	}
      */
      //end of proton part********************* 
      
      if (ElectronElectron && /*(((int)trigbit)&16) && trigdec>0 &&*/  isBest>=0 && oa > ang_cut /*&& eVertReco_z>-500 */ && mass_condition && proton_condition) 
	{
	  //cout<<"Event info: \n"<<"bt_condition:"<<bt_condition<<"\nPreShower: "<<pre_shower<<"\nFlanche: "<<flanch<<"\n Mass condition:" <<mass_condition<<"\nproton condition: "<<proton_condition<<endl;
	  if((em1_isring>0 || bt_em1_condition) && (em2_isring>0 || bt_em2_condition))// aLL DI-LEM1TONS
	    {
	      miss_mass_all_bcg2->Fill(pe1e2_miss->M());
	    }
	    
	  if(em1_isring>0 && em2_isring>0)//RF signal
	    { 
	      miss_mass_RF_bcg2->Fill(pe1e2_miss->M());
	      sig_all_back2->Fill(m_inv_e1e2);
	      sig_all_var_back2->Fill(m_inv_e1e2);
	    }
	  if( (bt_em1_condition||em1_isring) && (bt_em2_condition||em2_isring) && !(em2_isring && em1_isring))//BT profit
	    {
	      miss_mass_profit_bcg2->Fill(pe1e2_miss->M());
	    }
	  if(bt_condition)//bt signal  
	    {
	      miss_mass_BT_bcg2->Fill(pe1e2_miss->M());
	      sig_all_bt_back2->Fill(m_inv_e1e2);
	      sig_all_var_bt_back2->Fill(m_inv_e1e2);
	    }
	      
	}

      //tlo->fill();

      // if (Cut(ientry) < 0) continue;
   } // end of main loop
} // eof Loop 



// -------------------------------------------------------------------------------------------------
PPimEmEm::PPimEmEm(TTree *tree) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.

  em2_acc = em2_acc_err = em1_acc = em1_acc_err = 0.;
  em2_eff = em2_eff_err = em1_eff = em1_eff_err = 0.;

  if (tree == 0) {
	  
    TChain * chain = new TChain("PPimEmEm_ID","");

    //chain->Add("/lustre/nyx/hades/user/knowakow/PP/PAT/FILES/full_stat/hadron.root/PPimEmEm_ID");
    chain->Add("/lustre/nyx/hades/user/knowakow/PP/PAT_1/FILES/bt/lepton00.root/PPimEmEm_ID");
    chain->Add("/lustre/nyx/hades/user/knowakow/PP/PAT_1/FILES/bt/lepton01.root/PPimEmEm_ID");
    chain->Add("/lustre/nyx/hades/user/knowakow/PP/PAT_1/FILES/bt/lepton02.root/PPimEmEm_ID");
    chain->Add("/lustre/nyx/hades/user/knowakow/PP/PAT_1/FILES/bt/lepton03.root/PPimEmEm_ID");
    chain->Add("/lustre/nyx/hades/user/knowakow/PP/PAT_1/FILES/bt/lepton04.root/PPimEmEm_ID");
    chain->Add("/lustre/nyx/hades/user/knowakow/PP/PAT_1/FILES/bt/lepton05.root/PPimEmEm_ID");
    chain->Add("/lustre/nyx/hades/user/knowakow/PP/PAT_1/FILES/bt/lepton06.root/PPimEmEm_ID");
    chain->Add("/lustre/nyx/hades/user/knowakow/PP/PAT_1/FILES/bt/lepton07.root/PPimEmEm_ID");
    chain->Add("/lustre/nyx/hades/user/knowakow/PP/PAT_1/FILES/bt/lepton08.root/PPimEmEm_ID");
    chain->Add("/lustre/nyx/hades/user/knowakow/PP/PAT_1/FILES/bt/lepton09.root/PPimEmEm_ID");
    chain->Add("/lustre/nyx/hades/user/knowakow/PP/PAT_1/FILES/bt/lepton10.root/PPimEmEm_ID");
    chain->Add("/lustre/nyx/hades/user/knowakow/PP/PAT_1/FILES/bt/lepton11.root/PPimEmEm_ID");
   
    tree = chain; 
  }

  Init(tree);
}


// -------------------------------------------------------------------------------------------------
PPimEmEm::~PPimEmEm()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

// -------------------------------------------------------------------------------------------------
Int_t PPimEmEm::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}

// -------------------------------------------------------------------------------------------------
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

// -------------------------------------------------------------------------------------------------
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
   fChain->SetBranchAddress("em2_beta", &em2_beta, &b_em2_beta);
   fChain->SetBranchAddress("em2_beta_new", &em2_beta_new, &b_em2_beta_new);
   fChain->SetBranchAddress("em2_btChargeRing", &em2_btChargeRing, &b_em2_btChargeRing);
   fChain->SetBranchAddress("em2_btChargeSum", &em2_btChargeSum, &b_em2_btChargeSum);
   fChain->SetBranchAddress("em2_btChi2", &em2_btChi2, &b_em2_btChi2);
   fChain->SetBranchAddress("em2_btClusters", &em2_btClusters, &b_em2_btClusters);
   fChain->SetBranchAddress("em2_btMaxima", &em2_btMaxima, &b_em2_btMaxima);
   fChain->SetBranchAddress("em2_btMaximaCharge", &em2_btMaximaCharge, &b_em2_btMaximaCharge);
   fChain->SetBranchAddress("em2_btMaximaChargeShared", &em2_btMaximaChargeShared, &b_em2_btMaximaChargeShared);
   fChain->SetBranchAddress("em2_btMaximaChargeSharedFragment", &em2_btMaximaChargeSharedFragment, &b_em2_btMaximaChargeSharedFragment);
   fChain->SetBranchAddress("em2_btMaximaShared", &em2_btMaximaShared, &b_em2_btMaximaShared);
   fChain->SetBranchAddress("em2_btMaximaSharedFragment", &em2_btMaximaSharedFragment, &b_em2_btMaximaSharedFragment);
   fChain->SetBranchAddress("em2_btMeanDist", &em2_btMeanDist, &b_em2_btMeanDist);
   fChain->SetBranchAddress("em2_btNearbyMaxima", &em2_btNearbyMaxima, &b_em2_btNearbyMaxima);
   fChain->SetBranchAddress("em2_btNearbyMaximaShared", &em2_btNearbyMaximaShared, &b_em2_btNearbyMaximaShared);
   fChain->SetBranchAddress("em2_btPadsClus", &em2_btPadsClus, &b_em2_btPadsClus);
   fChain->SetBranchAddress("em2_btPadsRing", &em2_btPadsRing, &b_em2_btPadsRing);
   fChain->SetBranchAddress("em2_btRingMatrix", &em2_btRingMatrix, &b_em2_btRingMatrix);
   fChain->SetBranchAddress("em2_dedx_mdc", &em2_dedx_mdc, &b_em2_dedx_mdc);
   fChain->SetBranchAddress("em2_dedx_tof", &em2_dedx_tof, &b_em2_dedx_tof);
   fChain->SetBranchAddress("em2_id", &em2_id, &b_em2_id);
   fChain->SetBranchAddress("em2_isBT", &em2_isBT, &b_em2_isBT);
   fChain->SetBranchAddress("em2_isOffVertexClust", &em2_isOffVertexClust, &b_em2_isOffVertexClust);
   fChain->SetBranchAddress("em2_isPrimaryVertex", &em2_isPrimaryVertex, &b_em2_isPrimaryVertex);
   fChain->SetBranchAddress("em2_isUsedVertex", &em2_isUsedVertex, &b_em2_isUsedVertex);
   fChain->SetBranchAddress("em2_isring", &em2_isring, &b_em2_isring);
   fChain->SetBranchAddress("em2_isringmdc", &em2_isringmdc, &b_em2_isringmdc);
   fChain->SetBranchAddress("em2_isringnomatch", &em2_isringnomatch, &b_em2_isringnomatch);
   fChain->SetBranchAddress("em2_isringtrack", &em2_isringtrack, &b_em2_isringtrack);
   fChain->SetBranchAddress("em2_kIsLepton", &em2_kIsLepton, &b_em2_kIsLepton);
   fChain->SetBranchAddress("em2_kIsUsed", &em2_kIsUsed, &b_em2_kIsUsed);
   fChain->SetBranchAddress("em2_mdcinnerchi2", &em2_mdcinnerchi2, &b_em2_mdcinnerchi2);
   fChain->SetBranchAddress("em2_mdcouterchi2", &em2_mdcouterchi2, &b_em2_mdcouterchi2);
   fChain->SetBranchAddress("em2_oa_hadr", &em2_oa_hadr, &b_em2_oa_hadr);
   fChain->SetBranchAddress("em2_oa_lept", &em2_oa_lept, &b_em2_oa_lept);
   fChain->SetBranchAddress("em2_p", &em2_p, &b_em2_p);
   fChain->SetBranchAddress("em2_p_corr_em2", &em2_p_corr_em2, &b_em2_p_corr_em2);
   fChain->SetBranchAddress("em2_p_corr_em1", &em2_p_corr_em1, &b_em2_p_corr_em1);
   fChain->SetBranchAddress("em2_p_corr_p", &em2_p_corr_p, &b_em2_p_corr_p);
   fChain->SetBranchAddress("em2_p_corr_pim", &em2_p_corr_pim, &b_em2_p_corr_pim);
   fChain->SetBranchAddress("em2_p_corr_pip", &em2_p_corr_pip, &b_em2_p_corr_pip);
   fChain->SetBranchAddress("em2_phi", &em2_phi, &b_em2_phi);
   fChain->SetBranchAddress("em2_phi_rich", &em2_phi_rich, &b_em2_phi_rich);
   fChain->SetBranchAddress("em2_pid", &em2_pid, &b_em2_pid);
   fChain->SetBranchAddress("em2_q", &em2_q, &b_em2_q);
   fChain->SetBranchAddress("em2_r", &em2_r, &b_em2_r);
   fChain->SetBranchAddress("em2_resolution", &em2_resolution, &b_em2_resolution);
   fChain->SetBranchAddress("em2_resoultion", &em2_resoultion, &b_em2_resoultion);
   fChain->SetBranchAddress("em2_rich_amp", &em2_rich_amp, &b_em2_rich_amp);
   fChain->SetBranchAddress("em2_rich_centr", &em2_rich_centr, &b_em2_rich_centr);
   fChain->SetBranchAddress("em2_rich_houtra", &em2_rich_houtra, &b_em2_rich_houtra);
   fChain->SetBranchAddress("em2_rich_padnum", &em2_rich_padnum, &b_em2_rich_padnum);
   fChain->SetBranchAddress("em2_rich_patmat", &em2_rich_patmat, &b_em2_rich_patmat);
   fChain->SetBranchAddress("em2_rkchi2", &em2_rkchi2, &b_em2_rkchi2);
   fChain->SetBranchAddress("em2_sector", &em2_sector, &b_em2_sector);
   fChain->SetBranchAddress("em2_shw_sum0", &em2_shw_sum0, &b_em2_shw_sum0);
   fChain->SetBranchAddress("em2_shw_sum1", &em2_shw_sum1, &b_em2_shw_sum1);
   fChain->SetBranchAddress("em2_shw_sum2", &em2_shw_sum2, &b_em2_shw_sum2);
   fChain->SetBranchAddress("em2_system", &em2_system, &b_em2_system);
   fChain->SetBranchAddress("em2_theta", &em2_theta, &b_em2_theta);
   fChain->SetBranchAddress("em2_theta_rich", &em2_theta_rich, &b_em2_theta_rich);
   fChain->SetBranchAddress("em2_tof_mom", &em2_tof_mom, &b_em2_tof_mom);
   fChain->SetBranchAddress("em2_tof_new", &em2_tof_new, &b_em2_tof_new);
   fChain->SetBranchAddress("em2_tof_rec", &em2_tof_rec, &b_em2_tof_rec);
   fChain->SetBranchAddress("em2_track_length", &em2_track_length, &b_em2_track_length);
   fChain->SetBranchAddress("em2_tracklength", &em2_tracklength, &b_em2_tracklength);
   fChain->SetBranchAddress("em2_z", &em2_z, &b_em2_z);
   fChain->SetBranchAddress("em1_beta", &em1_beta, &b_em1_beta);
   fChain->SetBranchAddress("em1_beta_new", &em1_beta_new, &b_em1_beta_new);
   fChain->SetBranchAddress("em1_btChargeRing", &em1_btChargeRing, &b_em1_btChargeRing);
   fChain->SetBranchAddress("em1_btChargeSum", &em1_btChargeSum, &b_em1_btChargeSum);
   fChain->SetBranchAddress("em1_btChi2", &em1_btChi2, &b_em1_btChi2);
   fChain->SetBranchAddress("em1_btClusters", &em1_btClusters, &b_em1_btClusters);
   fChain->SetBranchAddress("em1_btMaxima", &em1_btMaxima, &b_em1_btMaxima);
   fChain->SetBranchAddress("em1_btMaximaCharge", &em1_btMaximaCharge, &b_em1_btMaximaCharge);
   fChain->SetBranchAddress("em1_btMaximaChargeShared", &em1_btMaximaChargeShared, &b_em1_btMaximaChargeShared);
   fChain->SetBranchAddress("em1_btMaximaChargeSharedFragment", &em1_btMaximaChargeSharedFragment, &b_em1_btMaximaChargeSharedFragment);
   fChain->SetBranchAddress("em1_btMaximaShared", &em1_btMaximaShared, &b_em1_btMaximaShared);
   fChain->SetBranchAddress("em1_btMaximaSharedFragment", &em1_btMaximaSharedFragment, &b_em1_btMaximaSharedFragment);
   fChain->SetBranchAddress("em1_btMeanDist", &em1_btMeanDist, &b_em1_btMeanDist);
   fChain->SetBranchAddress("em1_btNearbyMaxima", &em1_btNearbyMaxima, &b_em1_btNearbyMaxima);
   fChain->SetBranchAddress("em1_btNearbyMaximaShared", &em1_btNearbyMaximaShared, &b_em1_btNearbyMaximaShared);
   fChain->SetBranchAddress("em1_btPadsClus", &em1_btPadsClus, &b_em1_btPadsClus);
   fChain->SetBranchAddress("em1_btPadsRing", &em1_btPadsRing, &b_em1_btPadsRing);
   fChain->SetBranchAddress("em1_btRingMatrix", &em1_btRingMatrix, &b_em1_btRingMatrix);
   fChain->SetBranchAddress("em1_dedx_mdc", &em1_dedx_mdc, &b_em1_dedx_mdc);
   fChain->SetBranchAddress("em1_dedx_tof", &em1_dedx_tof, &b_em1_dedx_tof);
   fChain->SetBranchAddress("em1_id", &em1_id, &b_em1_id);
   fChain->SetBranchAddress("em1_isBT", &em1_isBT, &b_em1_isBT);
   fChain->SetBranchAddress("em1_isOffVertexClust", &em1_isOffVertexClust, &b_em1_isOffVertexClust);
   fChain->SetBranchAddress("em1_isPrimaryVertex", &em1_isPrimaryVertex, &b_em1_isPrimaryVertex);
   fChain->SetBranchAddress("em1_isUsedVertex", &em1_isUsedVertex, &b_em1_isUsedVertex);
   fChain->SetBranchAddress("em1_isring", &em1_isring, &b_em1_isring);
   fChain->SetBranchAddress("em1_isringmdc", &em1_isringmdc, &b_em1_isringmdc);
   fChain->SetBranchAddress("em1_isringnomatch", &em1_isringnomatch, &b_em1_isringnomatch);
   fChain->SetBranchAddress("em1_isringtrack", &em1_isringtrack, &b_em1_isringtrack);
   fChain->SetBranchAddress("em1_kIsLepton", &em1_kIsLepton, &b_em1_kIsLepton);
   fChain->SetBranchAddress("em1_kIsUsed", &em1_kIsUsed, &b_em1_kIsUsed);
   fChain->SetBranchAddress("em1_mdcinnerchi2", &em1_mdcinnerchi2, &b_em1_mdcinnerchi2);
   fChain->SetBranchAddress("em1_mdcouterchi2", &em1_mdcouterchi2, &b_em1_mdcouterchi2);
   fChain->SetBranchAddress("em1_oa_hadr", &em1_oa_hadr, &b_em1_oa_hadr);
   fChain->SetBranchAddress("em1_oa_lem1t", &em1_oa_lept, &b_em1_oa_lept);
   fChain->SetBranchAddress("em1_p", &em1_p, &b_em1_p);
   fChain->SetBranchAddress("em1_p_corr_em2", &em1_p_corr_em2, &b_em1_p_corr_em2);
   fChain->SetBranchAddress("em1_p_corr_em1", &em1_p_corr_em1, &b_em1_p_corr_em1);
   fChain->SetBranchAddress("em1_p_corr_p", &em1_p_corr_p, &b_em1_p_corr_p);
   fChain->SetBranchAddress("em1_p_corr_pim", &em1_p_corr_pim, &b_em1_p_corr_pim);
   fChain->SetBranchAddress("em1_p_corr_pip", &em1_p_corr_pip, &b_em1_p_corr_pip);
   fChain->SetBranchAddress("em1_phi", &em1_phi, &b_em1_phi);
   fChain->SetBranchAddress("em1_phi_rich", &em1_phi_rich, &b_em1_phi_rich);
   fChain->SetBranchAddress("em1_pid", &em1_pid, &b_em1_pid);
   fChain->SetBranchAddress("em1_q", &em1_q, &b_em1_q);
   fChain->SetBranchAddress("em1_r", &em1_r, &b_em1_r);
   fChain->SetBranchAddress("em1_resolution", &em1_resolution, &b_em1_resolution);
   fChain->SetBranchAddress("em1_resoultion", &em1_resoultion, &b_em1_resoultion);
   fChain->SetBranchAddress("em1_rich_amp", &em1_rich_amp, &b_em1_rich_amp);
   fChain->SetBranchAddress("em1_rich_centr", &em1_rich_centr, &b_em1_rich_centr);
   fChain->SetBranchAddress("em1_rich_houtra", &em1_rich_houtra, &b_em1_rich_houtra);
   fChain->SetBranchAddress("em1_rich_padnum", &em1_rich_padnum, &b_em1_rich_padnum);
   fChain->SetBranchAddress("em1_rich_patmat", &em1_rich_patmat, &b_em1_rich_patmat);
   fChain->SetBranchAddress("em1_rkchi2", &em1_rkchi2, &b_em1_rkchi2);
   fChain->SetBranchAddress("em1_sector", &em1_sector, &b_em1_sector);
   fChain->SetBranchAddress("em1_shw_sum0", &em1_shw_sum0, &b_em1_shw_sum0);
   fChain->SetBranchAddress("em1_shw_sum1", &em1_shw_sum1, &b_em1_shw_sum1);
   fChain->SetBranchAddress("em1_shw_sum2", &em1_shw_sum2, &b_em1_shw_sum2);
   fChain->SetBranchAddress("em1_system", &em1_system, &b_em1_system);
   fChain->SetBranchAddress("em1_theta", &em1_theta, &b_em1_theta);
   fChain->SetBranchAddress("em1_theta_rich", &em1_theta_rich, &b_em1_theta_rich);
   fChain->SetBranchAddress("em1_tof_mom", &em1_tof_mom, &b_em1_tof_mom);
   fChain->SetBranchAddress("em1_tof_new", &em1_tof_new, &b_em1_tof_new);
   fChain->SetBranchAddress("em1_tof_rec", &em1_tof_rec, &b_em1_tof_rec);
   fChain->SetBranchAddress("em1_track_length", &em1_track_length, &b_em1_track_length);
   fChain->SetBranchAddress("em1_tracklength", &em1_tracklength, &b_em1_tracklength);
   fChain->SetBranchAddress("em1_z", &em1_z, &b_em1_z);
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
   fChain->SetBranchAddress("p_p_corr_em2", &p_p_corr_em2, &b_p_p_corr_em2);
   fChain->SetBranchAddress("p_p_corr_em1", &p_p_corr_em1, &b_p_p_corr_em1);
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
   fChain->SetBranchAddress("pim_p_corr_em2", &pim_p_corr_em2, &b_pim_p_corr_em2);
   fChain->SetBranchAddress("pim_p_corr_em1", &pim_p_corr_em1, &b_pim_p_corr_em1);
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
   fChain->SetBranchAddress("runnumber", &runnumber, &b_runnumber);
   fChain->SetBranchAddress("totalmult", &totalmult, &b_totalmult);
   fChain->SetBranchAddress("trigbit", &trigbit, &b_trigbit);
   fChain->SetBranchAddress("trigdec", &trigdec, &b_trigdec);
   fChain->SetBranchAddress("trigdownscale", &trigdownscale, &b_trigdownscale);
   fChain->SetBranchAddress("trigdownscaleflag", &trigdownscaleflag, &b_trigdownscaleflag);
   Notify();
}

// -------------------------------------------------------------------------------------------------
Bool_t PPimEmEm::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

// -------------------------------------------------------------------------------------------------
void PPimEmEm::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}

// -------------------------------------------------------------------------------------------------
Int_t PPimEmEm::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is acceted.
// returns -1 otherwise.
   return 1;
}
