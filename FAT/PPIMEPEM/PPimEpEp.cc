#include "PPimEpEp.h"
#include "data.h"
#include <iostream>
#include <TLorentzVector.h>
#include <TVector3.h>
#include "hntuple.h"


using namespace std;
using namespace PATData;

void PPimEpEp::Loop()
{
    static long licznik = 0;

   if (fChain == 0) return;
/*
        (*tlo)["ep2_btChargeRing"] = 0;
        (*tlo)["ep2_btChargeSum"] = 0;
        (*tlo)["ep2_btChi2"] = 0;
        (*tlo)["ep2_btClusters"] = 0;
        (*tlo)["ep2_btMaxima"] = 0;
        (*tlo)["ep2_btMaximaCharge"] = 0;
        (*tlo)["ep2_btMaximaChargeShared"] = 0;
        (*tlo)["ep2_btMaximaChargeSharedFragment"] = 0;
        (*tlo)["ep2_btMaximaShared"] = 0;
        (*tlo)["ep2_btMaximaSharedFragment"] = 0;
        (*tlo)["ep2_btMeanDist"] = 0;
        (*tlo)["ep2_btNearbyMaxima"] = 0;
        (*tlo)["ep2_btNearbyMaximaShared"] = 0;
        (*tlo)["ep2_btPadsClus"] = 0;
        (*tlo)["ep2_btPadsRing"] = 0;
        (*tlo)["ep2_btRingMatrix"] = 0;
        (*tlo)["ep1_btChargeRing"] = 0;
        (*tlo)["ep1_btChargeSum"] = 0;
        (*tlo)["ep1_btChi2"] = 0;
        (*tlo)["ep1_btClusters"] = 0;
        (*tlo)["ep1_btMaxima"] = 0;
        (*tlo)["ep1_btMaximaCharge"] = 0;
        (*tlo)["ep1_btMaximaChargeShared"] = 0;
        (*tlo)["ep1_btMaximaChargeSharedFragment"] = 0;
        (*tlo)["ep1_btMaximaShared"] = 0;
        (*tlo)["ep1_btMaximaSharedFragment"] = 0;
        (*tlo)["ep1_btMeanDist"] = 0;
        (*tlo)["ep1_btNearbyMaxima"] = 0;
        (*tlo)["ep1_btNearbyMaximaShared"] = 0;
        (*tlo)["ep1_btPadsClus"] = 0;
        (*tlo)["ep1_btPadsRing"] = 0;
        (*tlo)["ep1_btRingMatrix"] = 0;


        (*tlo)["ep1_mom"] = 0;
        (*tlo)["ep1_theta"] = 0;
        (*tlo)["ep1_theta_rich"] = 0;
        (*tlo)["ep1_phi"] = 0;
        (*tlo)["ep1_phi_rich"] = 0;
        (*tlo)["ep1_beta"] = 0;
        (*tlo)["ep2_mom"] = 0;
        (*tlo)["ep2_theta"] = 0;
        (*tlo)["ep2_theta_rich"] = 0;
        (*tlo)["ep2_phi"] = 0;
        (*tlo)["ep2_phi_rich"] = 0;
        (*tlo)["ep2_beta"] = 0;
        (*tlo)["oa"] = 0;
        (*tlo)["oa_rich"] = 0;
        (*tlo)["sig"] = 0;
        (*tlo)["ep1_m"] = 0;
        (*tlo)["ep2_m"] = 0;
        (*tlo)["ep1ep2_inv_mass"] = 0;
        (*tlo)["ep1ep2_inv_mass2"] = 0;
        (*tlo)["ep1ep2_miss_mass"] = 0;
        (*tlo)["ep1ep2_miss_mass2"] = 0;
        (*tlo)["ep1ep2_y"] = 0;
        (*tlo)["ep1ep2_pt"] = 0;
        (*tlo)["ep1_rich_amp"] = 0;
        (*tlo)["ep1_rich_centr"] = 0;
        (*tlo)["ep1_rich_padnum"] = 0;
        (*tlo)["ep1_rich_patmat"] = 0;
        (*tlo)["ep1_rich_houtra"] = 0;
        (*tlo)["ep2_rich_amp"] = 0;
        (*tlo)["ep2_rich_centr"] = 0;
        (*tlo)["ep2_rich_padnum"] = 0;
        (*tlo)["ep2_rich_patmat"] = 0;
        (*tlo)["ep2_rich_houtra"] = 0;

	    (*tlo)["eVert_x"] = 0;
	    (*tlo)["eVert_y"] = 0;
	    (*tlo)["eVert_z"] = 0;

        (*tlo)["eVertReco_z"] = -1000.;
        (*tlo)["eVertReco_x"] = -1000.;
        (*tlo)["eVertReco_y"] = -1000.;
        (*tlo)["evtPileupMeta"] = 0.;
        (*tlo)["evtPileupStart"] = 0.;
        (*tlo)["ep1_isOffVertexClust"] = 0.;
        (*tlo)["ep1_p_corr_ep1"] = 0.;
        (*tlo)["ep2_isOffVertexClust"] = 0.;
        (*tlo)["ep2_p_corr_ep2"] = 0.;

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
      v2.SetXYZ(F*ep1_p*sin(D2R*ep1_theta)*cos(D2R*ep1_phi),F*ep1_p*sin(D2R*ep1_theta)*sin(D2R*ep1_phi),F*ep1_p*cos(D2R*ep1_theta));
      v3.SetXYZ(F*ep2_p*sin(D2R*ep2_theta)*cos(D2R*ep2_phi),F*ep2_p*sin(D2R*ep2_theta)*sin(D2R*ep2_phi),F*ep2_p*cos(D2R*ep2_theta));
      v4.SetXYZ(F*p_p*sin(D2R*p_theta)*cos(D2R*p_phi),F*p_p*sin(D2R*p_theta)*sin(D2R*p_phi),F*p_p*cos(D2R*p_theta));

      
      TVector3 r1, r2;
      r1.SetXYZ(sin(D2R*ep1_theta_rich)*cos(D2R*ep1_phi_rich),sin(D2R*ep1_theta_rich)*sin(D2R*ep1_phi_rich),cos(D2R*ep1_theta_rich));
      r2.SetXYZ(sin(D2R*ep2_theta_rich)*cos(D2R*ep2_phi_rich),sin(D2R*ep2_theta_rich)*sin(D2R*ep2_phi_rich),cos(D2R*ep2_theta_rich));

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

      double e1_mass = ep1_p*ep1_p * (  1. / (ep1_beta_new*ep1_beta_new)  - 1. ) ;
      double e2_mass = ep2_p*ep2_p * (  1. / (ep2_beta_new*ep2_beta_new)  - 1. ) ;

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
      insideEp2S0 = (pEp2S0 == 0) ? 0 : pEp2S0->IsInside(ep2_z,ep2_theta);
      insideEp2S1 = (pEp2S1 == 0) ? 0 : pEp2S1->IsInside(ep2_z,ep2_theta);
      insideEp1S0 = (pEp1S0 == 0) ? 0 : pEp1S0->IsInside(ep1_z,ep1_theta);
      insideEp1S1 = (pEp1S1 == 0) ? 0 : pEp1S1->IsInside(ep1_z,ep1_theta);
      //insideEp2S0 = (pEp2S0 == 0) ? 0 : pEp2S0->IsInside(eVert_z,ep2_theta);
      //insideEp2S1 = (pEp2S1 == 0) ? 0 : pEp2S1->IsInside(eVert_z,ep2_theta);
      //insideEp1S0 = (pEp1S0 == 0) ? 0 : pEp1S0->IsInside(eVert_z,ep1_theta);
      //insideEp1S1 = (pEp1S1 == 0) ? 0 : pEp1S1->IsInside(eVert_z,ep1_theta);
#endif

      insideTarget = 1;

#ifdef RECTANG
      insideEp2S0 = (ep2_theta > 50 && ep2_z < -50 /* && ep2_p<200.*/) ? 1 : 0;
      insideEp2S1 = (ep2_theta > 50 && ep2_z < -50 /* && ep2_p<200.*/) ? 1 : 0;
      insideEp1S0 = (ep1_theta > 50 && ep1_z < -50 /* && ep1_p<200.*/) ? 1 : 0;
       insideEp1S1 = (ep1_theta > 50 && ep1_z < -50 /* && ep1_p<200.*/) ? 1 : 0;
#endif

//#ifdef NOCUT
      insideEp2S0 = 0;
      insideEp2S1 = 0;
      insideEp1S0 = 0;
      insideEp1S1 = 0;
//#endif


      NoLeptonE1 = !((ep1_oa_lept< close_cut&&ep1_oa_lept>0.0) &&ep1_oa_lept>nonfit_close_cut );
      NoHadronE1 = 1;//(ep1_oa_hadr< close_cut &&ep1_oa_hadr>nonfit_close_cut );
      NoLeptonE2 = !((ep2_oa_lept< close_cut&&ep2_oa_lept>0.0) &&ep2_oa_lept>nonfit_close_cut );
      NoHadronE2 = 1;//(ep2_oa_hadr< close_cut &&ep2_oa_hadr>nonfit_close_cut );
      NoHadronE1 = 1;
      NoHadronE2 = 1;

/*
      NoLep1tonE1 = 1;
      NoHadronE1 = 1;
      NoLep1tonE2 = 1;
      NoHadronE2 = 1;
*/

      Positron1 = (((ep1_system==0&&insideEp1S0==0)||(ep1_system==1&&insideEp1S1==0)));
      Positron2 = (((ep2_system==0&&insideEp2S0==0)||(ep2_system==1&&insideEp2S1==0)));

      PositronPositron = Positron1 && NoLeptonE1 && NoHadronE1  &&  Positron2 && NoLeptonE2 && NoHadronE2  && insideTarget &&
                         ( e1_mass < 5000. && e2_mass < 5000. );

           
      bool bt_ep2_condition=(ep2_isBT!=-1
			    //&& ep2_btMaxima>=2
			    && ep2_btPadsRing>=2
			    );
      bool bt_ep1_condition=(ep1_isBT!=-1
			    //&& ep1_btMaxima>=2
			    && ep1_btPadsRing>=2
			    );
      bool bt_condition=(bt_ep2_condition && bt_ep1_condition);
      bool pre_shower= (ep1_system==0?(ep1_shw_sum1+ep1_shw_sum2-ep1_shw_sum0) > (parametrization(ep1_p)):true)
	             &&(ep2_system==0?(ep2_shw_sum1+ep2_shw_sum2-ep2_shw_sum0) > (parametrization(ep2_p)):true);
      //bool flanch=!(ep1_theta>65 && eVert_z<-55) && !(ep2_theta>65 && eVert_z<-55); 
      bool flanch = !((ep1_theta>65 || ep2_theta >65)&& p_z<-55);
      bool mass_condition=(ep1_p>50 && ep2_p>50 && ep1_p<2000. && ep2_p<2000.
			   //&& ep1_p>200 && ep2_p>200
			   //&&(ep1_systep2==0?ep1_beta_new>0.95:ep1_beta_new>0.92)&&(ep2_systep2==0?ep2_beta_new>0.95:ep2_beta_new>0.92)
			   //&& ep2_beta_new<1.1 && ep1_beta_new<1.1
			   && ep2_beta_new>0.8 && ep1_beta_new>0.8
			   && pre_shower
			   && flanch
			   && m_inv_e1e2>140
			   //&& TMath::Sqrt(p_p*p_p*(1/p_beta_new*1/p_beta_new-1))>750 //mass cut for protons
			   );
      bool proton_condition=(isBest==1)
	//&& TMath::Sqrt(p_p*p_p*(1/p_beta_new*1/p_beta_new-1))<1100
			     ;
      int i_array=0;
      //cout<<"PepEp inside loop no. "<<licznik<<endl;
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
      
      if (PositronPositron && /*(((int)trigbit)&16) && trigdec>0 &&*/  isBest>=0 && oa > ang_cut /*&& eVertReco_z>-500 */ && mass_condition && proton_condition) 
	{
	  if((ep1_isring>0 || bt_ep1_condition) && (ep2_isring>0 || bt_ep2_condition))// aLL DI-LEP1TONS
	    {
	      miss_mass_all_bcg1->Fill(pe1e2_miss->M());
	    }
	    
	  if(ep1_isring>0 && ep2_isring>0)//RF signal
	    { 
	      miss_mass_RF_bcg1->Fill(pe1e2_miss->M());
	      sig_all_back1->Fill(m_inv_e1e2);
	      sig_all_var_back1->Fill(m_inv_e1e2);
	    }
	  if( (bt_ep1_condition||ep1_isring) && (bt_ep2_condition||ep2_isring) && !(ep2_isring && ep1_isring))//BT profit
	    {
	      miss_mass_profit_bcg1->Fill(pe1e2_miss->M());
	    }
	  if(bt_condition)//bt signal  
	    {
	      miss_mass_BT_bcg1->Fill(pe1e2_miss->M());
	      sig_all_bt_back1->Fill(m_inv_e1e2);
	      sig_all_var_bt_back1->Fill(m_inv_e1e2);
	    }
	      
	}

      //tlo->fill();

      // if (Cut(ientry) < 0) continue;
   } // end of main loop
} // eof Loop 



// -------------------------------------------------------------------------------------------------
PPimEpEp::PPimEpEp(TTree *tree) 
{
  // if parameter tree is not specified (or zero), connect the file
  // used to generate this class and read the Tree.

  ep2_acc = ep2_acc_err = ep1_acc = ep1_acc_err = 0.;
  ep2_eff = ep2_eff_err = ep1_eff = ep1_eff_err = 0.;

  if (tree == 0) {
	  
    TChain * chain = new TChain("PPimEpEp_ID","");

    //chain->Add("/lustre/nyx/hades/user/knowakow/PP/PAT/FILES/full_stat/hadron.root/PPimEpEp_ID");
    chain->Add("/lustre/nyx/hades/user/knowakow/PP/PAT_1/FILES/ppimepem/lepton00.root/PPimEpEm_ID");
    chain->Add("/lustre/nyx/hades/user/knowakow/PP/PAT_1/FILES/ppimepem/lepton01.root/PPimEpEm_ID");
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

    tree = chain; 
  }

  Init(tree);
}


// -------------------------------------------------------------------------------------------------
PPimEpEp::~PPimEpEp()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

// -------------------------------------------------------------------------------------------------
Int_t PPimEpEp::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}

// -------------------------------------------------------------------------------------------------
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

// -------------------------------------------------------------------------------------------------
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
   fChain->SetBranchAddress("ep2_beta", &ep2_beta, &b_ep2_beta);
   fChain->SetBranchAddress("ep2_beta_new", &ep2_beta_new, &b_ep2_beta_new);
   fChain->SetBranchAddress("ep2_btChargeRing", &ep2_btChargeRing, &b_ep2_btChargeRing);
   fChain->SetBranchAddress("ep2_btChargeSum", &ep2_btChargeSum, &b_ep2_btChargeSum);
   fChain->SetBranchAddress("ep2_btChi2", &ep2_btChi2, &b_ep2_btChi2);
   fChain->SetBranchAddress("ep2_btClusters", &ep2_btClusters, &b_ep2_btClusters);
   fChain->SetBranchAddress("ep2_btMaxima", &ep2_btMaxima, &b_ep2_btMaxima);
   fChain->SetBranchAddress("ep2_btMaximaCharge", &ep2_btMaximaCharge, &b_ep2_btMaximaCharge);
   fChain->SetBranchAddress("ep2_btMaximaChargeShared", &ep2_btMaximaChargeShared, &b_ep2_btMaximaChargeShared);
   fChain->SetBranchAddress("ep2_btMaximaChargeSharedFragment", &ep2_btMaximaChargeSharedFragment, &b_ep2_btMaximaChargeSharedFragment);
   fChain->SetBranchAddress("ep2_btMaximaShared", &ep2_btMaximaShared, &b_ep2_btMaximaShared);
   fChain->SetBranchAddress("ep2_btMaximaSharedFragment", &ep2_btMaximaSharedFragment, &b_ep2_btMaximaSharedFragment);
   fChain->SetBranchAddress("ep2_btMeanDist", &ep2_btMeanDist, &b_ep2_btMeanDist);
   fChain->SetBranchAddress("ep2_btNearbyMaxima", &ep2_btNearbyMaxima, &b_ep2_btNearbyMaxima);
   fChain->SetBranchAddress("ep2_btNearbyMaximaShared", &ep2_btNearbyMaximaShared, &b_ep2_btNearbyMaximaShared);
   fChain->SetBranchAddress("ep2_btPadsClus", &ep2_btPadsClus, &b_ep2_btPadsClus);
   fChain->SetBranchAddress("ep2_btPadsRing", &ep2_btPadsRing, &b_ep2_btPadsRing);
   fChain->SetBranchAddress("ep2_btRingMatrix", &ep2_btRingMatrix, &b_ep2_btRingMatrix);
   fChain->SetBranchAddress("ep2_dedx_mdc", &ep2_dedx_mdc, &b_ep2_dedx_mdc);
   fChain->SetBranchAddress("ep2_dedx_tof", &ep2_dedx_tof, &b_ep2_dedx_tof);
   fChain->SetBranchAddress("ep2_id", &ep2_id, &b_ep2_id);
   fChain->SetBranchAddress("ep2_isBT", &ep2_isBT, &b_ep2_isBT);
   fChain->SetBranchAddress("ep2_isOffVertexClust", &ep2_isOffVertexClust, &b_ep2_isOffVertexClust);
   fChain->SetBranchAddress("ep2_isPrimaryVertex", &ep2_isPrimaryVertex, &b_ep2_isPrimaryVertex);
   fChain->SetBranchAddress("ep2_isUsedVertex", &ep2_isUsedVertex, &b_ep2_isUsedVertex);
   fChain->SetBranchAddress("ep2_isring", &ep2_isring, &b_ep2_isring);
   fChain->SetBranchAddress("ep2_isringmdc", &ep2_isringmdc, &b_ep2_isringmdc);
   fChain->SetBranchAddress("ep2_isringnomatch", &ep2_isringnomatch, &b_ep2_isringnomatch);
   fChain->SetBranchAddress("ep2_isringtrack", &ep2_isringtrack, &b_ep2_isringtrack);
   fChain->SetBranchAddress("ep2_kIsLepton", &ep2_kIsLepton, &b_ep2_kIsLepton);
   fChain->SetBranchAddress("ep2_kIsUsed", &ep2_kIsUsed, &b_ep2_kIsUsed);
   fChain->SetBranchAddress("ep2_mdcinnerchi2", &ep2_mdcinnerchi2, &b_ep2_mdcinnerchi2);
   fChain->SetBranchAddress("ep2_mdcouterchi2", &ep2_mdcouterchi2, &b_ep2_mdcouterchi2);
   fChain->SetBranchAddress("ep2_oa_hadr", &ep2_oa_hadr, &b_ep2_oa_hadr);
   fChain->SetBranchAddress("ep2_oa_lept", &ep2_oa_lept, &b_ep2_oa_lept);
   fChain->SetBranchAddress("ep2_p", &ep2_p, &b_ep2_p);
   fChain->SetBranchAddress("ep2_p_corr_ep2", &ep2_p_corr_ep2, &b_ep2_p_corr_ep2);
   fChain->SetBranchAddress("ep2_p_corr_ep1", &ep2_p_corr_ep1, &b_ep2_p_corr_ep1);
   fChain->SetBranchAddress("ep2_p_corr_p", &ep2_p_corr_p, &b_ep2_p_corr_p);
   fChain->SetBranchAddress("ep2_p_corr_pim", &ep2_p_corr_pim, &b_ep2_p_corr_pim);
   fChain->SetBranchAddress("ep2_p_corr_pip", &ep2_p_corr_pip, &b_ep2_p_corr_pip);
   fChain->SetBranchAddress("ep2_phi", &ep2_phi, &b_ep2_phi);
   fChain->SetBranchAddress("ep2_phi_rich", &ep2_phi_rich, &b_ep2_phi_rich);
   fChain->SetBranchAddress("ep2_pid", &ep2_pid, &b_ep2_pid);
   fChain->SetBranchAddress("ep2_q", &ep2_q, &b_ep2_q);
   fChain->SetBranchAddress("ep2_r", &ep2_r, &b_ep2_r);
   fChain->SetBranchAddress("ep2_resolution", &ep2_resolution, &b_ep2_resolution);
   fChain->SetBranchAddress("ep2_resoultion", &ep2_resoultion, &b_ep2_resoultion);
   fChain->SetBranchAddress("ep2_rich_amp", &ep2_rich_amp, &b_ep2_rich_amp);
   fChain->SetBranchAddress("ep2_rich_centr", &ep2_rich_centr, &b_ep2_rich_centr);
   fChain->SetBranchAddress("ep2_rich_houtra", &ep2_rich_houtra, &b_ep2_rich_houtra);
   fChain->SetBranchAddress("ep2_rich_padnum", &ep2_rich_padnum, &b_ep2_rich_padnum);
   fChain->SetBranchAddress("ep2_rich_patmat", &ep2_rich_patmat, &b_ep2_rich_patmat);
   fChain->SetBranchAddress("ep2_rkchi2", &ep2_rkchi2, &b_ep2_rkchi2);
   fChain->SetBranchAddress("ep2_sector", &ep2_sector, &b_ep2_sector);
   fChain->SetBranchAddress("ep2_shw_sum0", &ep2_shw_sum0, &b_ep2_shw_sum0);
   fChain->SetBranchAddress("ep2_shw_sum1", &ep2_shw_sum1, &b_ep2_shw_sum1);
   fChain->SetBranchAddress("ep2_shw_sum2", &ep2_shw_sum2, &b_ep2_shw_sum2);
   fChain->SetBranchAddress("ep2_system", &ep2_system, &b_ep2_system);
   fChain->SetBranchAddress("ep2_theta", &ep2_theta, &b_ep2_theta);
   fChain->SetBranchAddress("ep2_theta_rich", &ep2_theta_rich, &b_ep2_theta_rich);
   fChain->SetBranchAddress("ep2_tof_mom", &ep2_tof_mom, &b_ep2_tof_mom);
   fChain->SetBranchAddress("ep2_tof_new", &ep2_tof_new, &b_ep2_tof_new);
   fChain->SetBranchAddress("ep2_tof_rec", &ep2_tof_rec, &b_ep2_tof_rec);
   fChain->SetBranchAddress("ep2_track_length", &ep2_track_length, &b_ep2_track_length);
   fChain->SetBranchAddress("ep2_tracklength", &ep2_tracklength, &b_ep2_tracklength);
   fChain->SetBranchAddress("ep2_z", &ep2_z, &b_ep2_z);
   fChain->SetBranchAddress("ep1_beta", &ep1_beta, &b_ep1_beta);
   fChain->SetBranchAddress("ep1_beta_new", &ep1_beta_new, &b_ep1_beta_new);
   fChain->SetBranchAddress("ep1_btChargeRing", &ep1_btChargeRing, &b_ep1_btChargeRing);
   fChain->SetBranchAddress("ep1_btChargeSum", &ep1_btChargeSum, &b_ep1_btChargeSum);
   fChain->SetBranchAddress("ep1_btChi2", &ep1_btChi2, &b_ep1_btChi2);
   fChain->SetBranchAddress("ep1_btClusters", &ep1_btClusters, &b_ep1_btClusters);
   fChain->SetBranchAddress("ep1_btMaxima", &ep1_btMaxima, &b_ep1_btMaxima);
   fChain->SetBranchAddress("ep1_btMaximaCharge", &ep1_btMaximaCharge, &b_ep1_btMaximaCharge);
   fChain->SetBranchAddress("ep1_btMaximaChargeShared", &ep1_btMaximaChargeShared, &b_ep1_btMaximaChargeShared);
   fChain->SetBranchAddress("ep1_btMaximaChargeSharedFragment", &ep1_btMaximaChargeSharedFragment, &b_ep1_btMaximaChargeSharedFragment);
   fChain->SetBranchAddress("ep1_btMaximaShared", &ep1_btMaximaShared, &b_ep1_btMaximaShared);
   fChain->SetBranchAddress("ep1_btMaximaSharedFragment", &ep1_btMaximaSharedFragment, &b_ep1_btMaximaSharedFragment);
   fChain->SetBranchAddress("ep1_btMeanDist", &ep1_btMeanDist, &b_ep1_btMeanDist);
   fChain->SetBranchAddress("ep1_btNearbyMaxima", &ep1_btNearbyMaxima, &b_ep1_btNearbyMaxima);
   fChain->SetBranchAddress("ep1_btNearbyMaximaShared", &ep1_btNearbyMaximaShared, &b_ep1_btNearbyMaximaShared);
   fChain->SetBranchAddress("ep1_btPadsClus", &ep1_btPadsClus, &b_ep1_btPadsClus);
   fChain->SetBranchAddress("ep1_btPadsRing", &ep1_btPadsRing, &b_ep1_btPadsRing);
   fChain->SetBranchAddress("ep1_btRingMatrix", &ep1_btRingMatrix, &b_ep1_btRingMatrix);
   fChain->SetBranchAddress("ep1_dedx_mdc", &ep1_dedx_mdc, &b_ep1_dedx_mdc);
   fChain->SetBranchAddress("ep1_dedx_tof", &ep1_dedx_tof, &b_ep1_dedx_tof);
   fChain->SetBranchAddress("ep1_id", &ep1_id, &b_ep1_id);
   fChain->SetBranchAddress("ep1_isBT", &ep1_isBT, &b_ep1_isBT);
   fChain->SetBranchAddress("ep1_isOffVertexClust", &ep1_isOffVertexClust, &b_ep1_isOffVertexClust);
   fChain->SetBranchAddress("ep1_isPrimaryVertex", &ep1_isPrimaryVertex, &b_ep1_isPrimaryVertex);
   fChain->SetBranchAddress("ep1_isUsedVertex", &ep1_isUsedVertex, &b_ep1_isUsedVertex);
   fChain->SetBranchAddress("ep1_isring", &ep1_isring, &b_ep1_isring);
   fChain->SetBranchAddress("ep1_isringmdc", &ep1_isringmdc, &b_ep1_isringmdc);
   fChain->SetBranchAddress("ep1_isringnomatch", &ep1_isringnomatch, &b_ep1_isringnomatch);
   fChain->SetBranchAddress("ep1_isringtrack", &ep1_isringtrack, &b_ep1_isringtrack);
   fChain->SetBranchAddress("ep1_kIsLepton", &ep1_kIsLepton, &b_ep1_kIsLepton);
   fChain->SetBranchAddress("ep1_kIsUsed", &ep1_kIsUsed, &b_ep1_kIsUsed);
   fChain->SetBranchAddress("ep1_mdcinnerchi2", &ep1_mdcinnerchi2, &b_ep1_mdcinnerchi2);
   fChain->SetBranchAddress("ep1_mdcouterchi2", &ep1_mdcouterchi2, &b_ep1_mdcouterchi2);
   fChain->SetBranchAddress("ep1_oa_hadr", &ep1_oa_hadr, &b_ep1_oa_hadr);
   fChain->SetBranchAddress("ep1_oa_lep1t", &ep1_oa_lept, &b_ep1_oa_lept);
   fChain->SetBranchAddress("ep1_p", &ep1_p, &b_ep1_p);
   fChain->SetBranchAddress("ep1_p_corr_ep2", &ep1_p_corr_ep2, &b_ep1_p_corr_ep2);
   fChain->SetBranchAddress("ep1_p_corr_ep1", &ep1_p_corr_ep1, &b_ep1_p_corr_ep1);
   fChain->SetBranchAddress("ep1_p_corr_p", &ep1_p_corr_p, &b_ep1_p_corr_p);
   fChain->SetBranchAddress("ep1_p_corr_pim", &ep1_p_corr_pim, &b_ep1_p_corr_pim);
   fChain->SetBranchAddress("ep1_p_corr_pip", &ep1_p_corr_pip, &b_ep1_p_corr_pip);
   fChain->SetBranchAddress("ep1_phi", &ep1_phi, &b_ep1_phi);
   fChain->SetBranchAddress("ep1_phi_rich", &ep1_phi_rich, &b_ep1_phi_rich);
   fChain->SetBranchAddress("ep1_pid", &ep1_pid, &b_ep1_pid);
   fChain->SetBranchAddress("ep1_q", &ep1_q, &b_ep1_q);
   fChain->SetBranchAddress("ep1_r", &ep1_r, &b_ep1_r);
   fChain->SetBranchAddress("ep1_resolution", &ep1_resolution, &b_ep1_resolution);
   fChain->SetBranchAddress("ep1_resoultion", &ep1_resoultion, &b_ep1_resoultion);
   fChain->SetBranchAddress("ep1_rich_amp", &ep1_rich_amp, &b_ep1_rich_amp);
   fChain->SetBranchAddress("ep1_rich_centr", &ep1_rich_centr, &b_ep1_rich_centr);
   fChain->SetBranchAddress("ep1_rich_houtra", &ep1_rich_houtra, &b_ep1_rich_houtra);
   fChain->SetBranchAddress("ep1_rich_padnum", &ep1_rich_padnum, &b_ep1_rich_padnum);
   fChain->SetBranchAddress("ep1_rich_patmat", &ep1_rich_patmat, &b_ep1_rich_patmat);
   fChain->SetBranchAddress("ep1_rkchi2", &ep1_rkchi2, &b_ep1_rkchi2);
   fChain->SetBranchAddress("ep1_sector", &ep1_sector, &b_ep1_sector);
   fChain->SetBranchAddress("ep1_shw_sum0", &ep1_shw_sum0, &b_ep1_shw_sum0);
   fChain->SetBranchAddress("ep1_shw_sum1", &ep1_shw_sum1, &b_ep1_shw_sum1);
   fChain->SetBranchAddress("ep1_shw_sum2", &ep1_shw_sum2, &b_ep1_shw_sum2);
   fChain->SetBranchAddress("ep1_system", &ep1_system, &b_ep1_system);
   fChain->SetBranchAddress("ep1_theta", &ep1_theta, &b_ep1_theta);
   fChain->SetBranchAddress("ep1_theta_rich", &ep1_theta_rich, &b_ep1_theta_rich);
   fChain->SetBranchAddress("ep1_tof_mom", &ep1_tof_mom, &b_ep1_tof_mom);
   fChain->SetBranchAddress("ep1_tof_new", &ep1_tof_new, &b_ep1_tof_new);
   fChain->SetBranchAddress("ep1_tof_rec", &ep1_tof_rec, &b_ep1_tof_rec);
   fChain->SetBranchAddress("ep1_track_length", &ep1_track_length, &b_ep1_track_length);
   fChain->SetBranchAddress("ep1_tracklength", &ep1_tracklength, &b_ep1_tracklength);
   fChain->SetBranchAddress("ep1_z", &ep1_z, &b_ep1_z);
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
   fChain->SetBranchAddress("p_p_corr_ep2", &p_p_corr_ep2, &b_p_p_corr_ep2);
   fChain->SetBranchAddress("p_p_corr_ep1", &p_p_corr_ep1, &b_p_p_corr_ep1);
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
   fChain->SetBranchAddress("pim_p_corr_ep2", &pim_p_corr_ep2, &b_pim_p_corr_ep2);
   fChain->SetBranchAddress("pim_p_corr_ep1", &pim_p_corr_ep1, &b_pim_p_corr_ep1);
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
Bool_t PPimEpEp::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

// -------------------------------------------------------------------------------------------------
void PPimEpEp::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}

// -------------------------------------------------------------------------------------------------
Int_t PPimEpEp::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is acceted.
// returns -1 otherwise.
   return 1;
}
