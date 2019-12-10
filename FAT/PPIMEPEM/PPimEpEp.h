#ifndef PPimEpEp_h
#define PPimEpEp_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

using namespace std;

/*********************************************************************************************/
class PPimEpEp {

public :

   float ep2_acc, ep2_acc_err, ep1_acc, ep1_acc_err;
   float ep2_eff, ep2_eff_err, ep1_eff, ep1_eff_err;

   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain
   // Declaration of leaf types
   // Declaration of leaf types
   Float_t         eVertClust_chi2;
   Float_t         eVertClust_x;
   Float_t         eVertClust_y;
   Float_t         eVertClust_z;
   Float_t         eVertReco_chi2;
   Float_t         eVertReco_x;
   Float_t         eVertReco_y;
   Float_t         eVertReco_z;
   Float_t         eVert_chi2;
   Float_t         eVert_x;
   Float_t         eVert_y;
   Float_t         eVert_z;
   Float_t         ep2_beta;
   Float_t         ep2_beta_new;
   Float_t         ep2_btChargeRing;
   Float_t         ep2_btChargeSum;
   Float_t         ep2_btChi2;
   Float_t         ep2_btClusters;
   Float_t         ep2_btMaxima;
   Float_t         ep2_btMaximaCharge;
   Float_t         ep2_btMaximaChargeShared;
   Float_t         ep2_btMaximaChargeSharedFragment;
   Float_t         ep2_btMaximaShared;
   Float_t         ep2_btMaximaSharedFragment;
   Float_t         ep2_btMeanDist;
   Float_t         ep2_btNearbyMaxima;
   Float_t         ep2_btNearbyMaximaShared;
   Float_t         ep2_btPadsClus;
   Float_t         ep2_btPadsRing;
   Float_t         ep2_btRingMatrix;
   Float_t         ep2_dedx_mdc;
   Float_t         ep2_dedx_tof;
   Float_t         ep2_id;
   Float_t         ep2_isBT;
   Float_t         ep2_isOffVertexClust;
   Float_t         ep2_isPrimaryVertex;
   Float_t         ep2_isUsedVertex;
   Float_t         ep2_isring;
   Float_t         ep2_isringmdc;
   Float_t         ep2_isringnomatch;
   Float_t         ep2_isringtrack;
   Float_t         ep2_kIsLepton;
   Float_t         ep2_kIsUsed;
   Float_t         ep2_mdcinnerchi2;
   Float_t         ep2_mdcouterchi2;
   Float_t         ep2_oa_hadr;
   Float_t         ep2_oa_lept;
   Float_t         ep2_p;
   Float_t         ep2_p_corr_ep2;
   Float_t         ep2_p_corr_ep1;
   Float_t         ep2_p_corr_p;
   Float_t         ep2_p_corr_pim;
   Float_t         ep2_p_corr_pip;
   Float_t         ep2_phi;
   Float_t         ep2_phi_rich;
   Float_t         ep2_pid;
   Float_t         ep2_q;
   Float_t         ep2_r;
   Float_t         ep2_resolution;
   Float_t         ep2_resoultion;
   Float_t         ep2_rich_amp;
   Float_t         ep2_rich_centr;
   Float_t         ep2_rich_houtra;
   Float_t         ep2_rich_padnum;
   Float_t         ep2_rich_patmat;
   Float_t         ep2_rkchi2;
   Float_t         ep2_sector;
   Float_t         ep2_shw_sum0;
   Float_t         ep2_shw_sum1;
   Float_t         ep2_shw_sum2;
   Float_t         ep2_system;
   Float_t         ep2_theta;
   Float_t         ep2_theta_rich;
   Float_t         ep2_tof_mom;
   Float_t         ep2_tof_new;
   Float_t         ep2_tof_rec;
   Float_t         ep2_track_length;
   Float_t         ep2_tracklength;
   Float_t         ep2_z;
   Float_t         ep1_beta;
   Float_t         ep1_beta_new;
   Float_t         ep1_btChargeRing;
   Float_t         ep1_btChargeSum;
   Float_t         ep1_btChi2;
   Float_t         ep1_btClusters;
   Float_t         ep1_btMaxima;
   Float_t         ep1_btMaximaCharge;
   Float_t         ep1_btMaximaChargeShared;
   Float_t         ep1_btMaximaChargeSharedFragment;
   Float_t         ep1_btMaximaShared;
   Float_t         ep1_btMaximaSharedFragment;
   Float_t         ep1_btMeanDist;
   Float_t         ep1_btNearbyMaxima;
   Float_t         ep1_btNearbyMaximaShared;
   Float_t         ep1_btPadsClus;
   Float_t         ep1_btPadsRing;
   Float_t         ep1_btRingMatrix;
   Float_t         ep1_dedx_mdc;
   Float_t         ep1_dedx_tof;
   Float_t         ep1_id;
   Float_t         ep1_isBT;
   Float_t         ep1_isOffVertexClust;
   Float_t         ep1_isPrimaryVertex;
   Float_t         ep1_isUsedVertex;
   Float_t         ep1_isring;
   Float_t         ep1_isringmdc;
   Float_t         ep1_isringnomatch;
   Float_t         ep1_isringtrack;
   Float_t         ep1_kIsLepton;
   Float_t         ep1_kIsUsed;
   Float_t         ep1_mdcinnerchi2;
   Float_t         ep1_mdcouterchi2;
   Float_t         ep1_oa_hadr;
   Float_t         ep1_oa_lept;
   Float_t         ep1_p;
   Float_t         ep1_p_corr_ep2;
   Float_t         ep1_p_corr_ep1;
   Float_t         ep1_p_corr_p;
   Float_t         ep1_p_corr_pim;
   Float_t         ep1_p_corr_pip;
   Float_t         ep1_phi;
   Float_t         ep1_phi_rich;
   Float_t         ep1_pid;
   Float_t         ep1_q;
   Float_t         ep1_r;
   Float_t         ep1_resolution;
   Float_t         ep1_resoultion;
   Float_t         ep1_rich_amp;
   Float_t         ep1_rich_centr;
   Float_t         ep1_rich_houtra;
   Float_t         ep1_rich_padnum;
   Float_t         ep1_rich_patmat;
   Float_t         ep1_rkchi2;
   Float_t         ep1_sector;
   Float_t         ep1_shw_sum0;
   Float_t         ep1_shw_sum1;
   Float_t         ep1_shw_sum2;
   Float_t         ep1_system;
   Float_t         ep1_theta;
   Float_t         ep1_theta_rich;
   Float_t         ep1_tof_mom;
   Float_t         ep1_tof_new;
   Float_t         ep1_tof_rec;
   Float_t         ep1_track_length;
   Float_t         ep1_tracklength;
   Float_t         ep1_z;
   Float_t         event;
   Float_t         fw_beta_1;
   Float_t         fw_beta_2;
   Float_t         fw_beta_3;
   Float_t         fw_charge_1;
   Float_t         fw_charge_2;
   Float_t         fw_charge_3;
   Float_t         fw_cluster_mult;
   Float_t         fw_distance_1;
   Float_t         fw_distance_2;
   Float_t         fw_distance_3;
   Float_t         fw_mult;
   Float_t         fw_p_1;
   Float_t         fw_p_2;
   Float_t         fw_p_3;
   Float_t         fw_phi_1;
   Float_t         fw_phi_2;
   Float_t         fw_phi_3;
   Float_t         fw_size_1;
   Float_t         fw_size_2;
   Float_t         fw_size_3;
   Float_t         fw_spectator_1;
   Float_t         fw_spectator_2;
   Float_t         fw_spectator_3;
   Float_t         fw_theta_1;
   Float_t         fw_theta_2;
   Float_t         fw_theta_3;
   Float_t         fw_time_1;
   Float_t         fw_time_2;
   Float_t         fw_time_3;
   Float_t         fw_time_min_1;
   Float_t         fw_time_min_2;
   Float_t         fw_time_min_3;
   Float_t         fw_x_lab_1;
   Float_t         fw_x_lab_2;
   Float_t         fw_x_lab_3;
   Float_t         fw_y_lab_1;
   Float_t         fw_y_lab_2;
   Float_t         fw_y_lab_3;
   Float_t         fw_z_lab_1;
   Float_t         fw_z_lab_2;
   Float_t         fw_z_lab_3;
   Float_t         hneg_mult;
   Float_t         hpos_mult;
   Float_t         isBest;
   Float_t         lneg_mult;
   Float_t         lpos_mult;
   Float_t         p_beta;
   Float_t         p_beta_new;
   Float_t         p_btChargeRing;
   Float_t         p_btChargeSum;
   Float_t         p_btChi2;
   Float_t         p_btClusters;
   Float_t         p_btMaxima;
   Float_t         p_btMaximaCharge;
   Float_t         p_btMaximaChargeShared;
   Float_t         p_btMaximaChargeSharedFragment;
   Float_t         p_btMaximaShared;
   Float_t         p_btMaximaSharedFragment;
   Float_t         p_btMeanDist;
   Float_t         p_btNearbyMaxima;
   Float_t         p_btNearbyMaximaShared;
   Float_t         p_btPadsClus;
   Float_t         p_btPadsRing;
   Float_t         p_btRingMatrix;
   Float_t         p_dedx_mdc;
   Float_t         p_dedx_tof;
   Float_t         p_id;
   Float_t         p_isBT;
   Float_t         p_isOffVertexClust;
   Float_t         p_isPrimaryVertex;
   Float_t         p_isUsedVertex;
   Float_t         p_isring;
   Float_t         p_isringmdc;
   Float_t         p_isringnomatch;
   Float_t         p_isringtrack;
   Float_t         p_kIsLepton;
   Float_t         p_kIsUsed;
   Float_t         p_mdcinnerchi2;
   Float_t         p_mdcouterchi2;
   Float_t         p_oa_hadr;
   Float_t         p_oa_lept;
   Float_t         p_p;
   Float_t         p_p_corr_ep2;
   Float_t         p_p_corr_ep1;
   Float_t         p_p_corr_p;
   Float_t         p_p_corr_pim;
   Float_t         p_p_corr_pip;
   Float_t         p_phi;
   Float_t         p_phi_rich;
   Float_t         p_pid;
   Float_t         p_q;
   Float_t         p_r;
   Float_t         p_resolution;
   Float_t         p_resoultion;
   Float_t         p_rich_amp;
   Float_t         p_rich_centr;
   Float_t         p_rich_houtra;
   Float_t         p_rich_padnum;
   Float_t         p_rich_patmat;
   Float_t         p_rkchi2;
   Float_t         p_sector;
   Float_t         p_shw_sum0;
   Float_t         p_shw_sum1;
   Float_t         p_shw_sum2;
   Float_t         p_system;
   Float_t         p_theta;
   Float_t         p_theta_rich;
   Float_t         p_tof_mom;
   Float_t         p_tof_new;
   Float_t         p_tof_rec;
   Float_t         p_track_length;
   Float_t         p_tracklength;
   Float_t         p_z;

  Float_t         pim_beta;
   Float_t         pim_beta_new;
   Float_t         pim_btChargeRing;
   Float_t         pim_btChargeSum;
   Float_t         pim_btChi2;
   Float_t         pim_btClusters;
   Float_t         pim_btMaxima;
   Float_t         pim_btMaximaCharge;
   Float_t         pim_btMaximaChargeShared;
   Float_t         pim_btMaximaChargeSharedFragment;
   Float_t         pim_btMaximaShared;
   Float_t         pim_btMaximaSharedFragment;
   Float_t         pim_btMeanDist;
   Float_t         pim_btNearbyMaxima;
   Float_t         pim_btNearbyMaximaShared;
   Float_t         pim_btPadsClus;
   Float_t         pim_btPadsRing;
   Float_t         pim_btRingMatrix;
   Float_t         pim_dedx_mdc;
   Float_t         pim_dedx_tof;
   Float_t         pim_id;
   Float_t         pim_isBT;
   Float_t         pim_isOffVertexClust;
   Float_t         pim_isPrimaryVertex;
   Float_t         pim_isUsedVertex;
   Float_t         pim_isring;
   Float_t         pim_isringmdc;
   Float_t         pim_isringnomatch;
   Float_t         pim_isringtrack;
   Float_t         pim_kIsLepton;
   Float_t         pim_kIsUsed;
   Float_t         pim_mdcinnerchi2;
   Float_t         pim_mdcouterchi2;
   Float_t         pim_oa_hadr;
   Float_t         pim_oa_lept;
   Float_t         pim_p;
   Float_t         pim_p_corr_ep2;
   Float_t         pim_p_corr_ep1;
   Float_t         pim_p_corr_p;
   Float_t         pim_p_corr_pim;
   Float_t         pim_p_corr_pip;
   Float_t         pim_phi;
   Float_t         pim_phi_rich;
   Float_t         pim_pid;
   Float_t         pim_q;
   Float_t         pim_r;
   Float_t         pim_resolution;
   Float_t         pim_resoultion;
   Float_t         pim_rich_amp;
   Float_t         pim_rich_centr;
   Float_t         pim_rich_houtra;
   Float_t         pim_rich_padnum;
   Float_t         pim_rich_patmat;
   Float_t         pim_rkchi2;
   Float_t         pim_sector;
   Float_t         pim_shw_sum0;
   Float_t         pim_shw_sum1;
   Float_t         pim_shw_sum2;
   Float_t         pim_system;
   Float_t         pim_theta;
   Float_t         pim_theta_rich;
   Float_t         pim_tof_mom;
   Float_t         pim_tof_new;
   Float_t         pim_tof_rec;
   Float_t         pim_track_length;
   Float_t         pim_tracklength;
   Float_t         pim_z;
  
  Float_t         runnumber;
   Float_t         totalmult;
   Float_t         trigbit;
   Float_t         trigdec;
   Float_t         trigdownscale;
   Float_t         trigdownscaleflag;

   // List of branches
   TBranch        *b_eVertClust_chi2;   //!
   TBranch        *b_eVertClust_x;   //!
   TBranch        *b_eVertClust_y;   //!
   TBranch        *b_eVertClust_z;   //!
   TBranch        *b_eVertReco_chi2;   //!
   TBranch        *b_eVertReco_x;   //!
   TBranch        *b_eVertReco_y;   //!
   TBranch        *b_eVertReco_z;   //!
   TBranch        *b_eVert_chi2;   //!
   TBranch        *b_eVert_x;   //!
   TBranch        *b_eVert_y;   //!
   TBranch        *b_eVert_z;   //!
   TBranch        *b_ep2_beta;   //!
   TBranch        *b_ep2_beta_new;   //!
   TBranch        *b_ep2_btChargeRing;   //!
   TBranch        *b_ep2_btChargeSum;   //!
   TBranch        *b_ep2_btChi2;   //!
   TBranch        *b_ep2_btClusters;   //!
   TBranch        *b_ep2_btMaxima;   //!
   TBranch        *b_ep2_btMaximaCharge;   //!
   TBranch        *b_ep2_btMaximaChargeShared;   //!
   TBranch        *b_ep2_btMaximaChargeSharedFragment;   //!
   TBranch        *b_ep2_btMaximaShared;   //!
   TBranch        *b_ep2_btMaximaSharedFragment;   //!
   TBranch        *b_ep2_btMeanDist;   //!
   TBranch        *b_ep2_btNearbyMaxima;   //!
   TBranch        *b_ep2_btNearbyMaximaShared;   //!
   TBranch        *b_ep2_btPadsClus;   //!
   TBranch        *b_ep2_btPadsRing;   //!
   TBranch        *b_ep2_btRingMatrix;   //!
   TBranch        *b_ep2_dedx_mdc;   //!
   TBranch        *b_ep2_dedx_tof;   //!
   TBranch        *b_ep2_id;   //!
   TBranch        *b_ep2_isBT;   //!
   TBranch        *b_ep2_isOffVertexClust;   //!
   TBranch        *b_ep2_isPrimaryVertex;   //!
   TBranch        *b_ep2_isUsedVertex;   //!
   TBranch        *b_ep2_isring;   //!
   TBranch        *b_ep2_isringmdc;   //!
   TBranch        *b_ep2_isringnomatch;   //!
   TBranch        *b_ep2_isringtrack;   //!
   TBranch        *b_ep2_kIsLepton;   //!
   TBranch        *b_ep2_kIsUsed;   //!
   TBranch        *b_ep2_mdcinnerchi2;   //!
   TBranch        *b_ep2_mdcouterchi2;   //!
   TBranch        *b_ep2_oa_hadr;   //!
   TBranch        *b_ep2_oa_lept;   //!
   TBranch        *b_ep2_p;   //!
   TBranch        *b_ep2_p_corr_ep2;   //!
   TBranch        *b_ep2_p_corr_ep1;   //!
   TBranch        *b_ep2_p_corr_p;   //!
   TBranch        *b_ep2_p_corr_pim;   //!
   TBranch        *b_ep2_p_corr_pip;   //!
   TBranch        *b_ep2_phi;   //!
   TBranch        *b_ep2_phi_rich;   //!
   TBranch        *b_ep2_pid;   //!
   TBranch        *b_ep2_q;   //!
   TBranch        *b_ep2_r;   //!
   TBranch        *b_ep2_resolution;   //!
   TBranch        *b_ep2_resoultion;   //!
   TBranch        *b_ep2_rich_amp;   //!
   TBranch        *b_ep2_rich_centr;   //!
   TBranch        *b_ep2_rich_houtra;   //!
   TBranch        *b_ep2_rich_padnum;   //!
   TBranch        *b_ep2_rich_patmat;   //!
   TBranch        *b_ep2_rkchi2;   //!
   TBranch        *b_ep2_sector;   //!
   TBranch        *b_ep2_shw_sum0;   //!
   TBranch        *b_ep2_shw_sum1;   //!
   TBranch        *b_ep2_shw_sum2;   //!
   TBranch        *b_ep2_system;   //!
   TBranch        *b_ep2_theta;   //!
   TBranch        *b_ep2_theta_rich;   //!
   TBranch        *b_ep2_tof_mom;   //!
   TBranch        *b_ep2_tof_new;   //!
   TBranch        *b_ep2_tof_rec;   //!
   TBranch        *b_ep2_track_length;   //!
   TBranch        *b_ep2_tracklength;   //!
   TBranch        *b_ep2_z;   //!
   TBranch        *b_ep1_beta;   //!
   TBranch        *b_ep1_beta_new;   //!
   TBranch        *b_ep1_btChargeRing;   //!
   TBranch        *b_ep1_btChargeSum;   //!
   TBranch        *b_ep1_btChi2;   //!
   TBranch        *b_ep1_btClusters;   //!
   TBranch        *b_ep1_btMaxima;   //!
   TBranch        *b_ep1_btMaximaCharge;   //!
   TBranch        *b_ep1_btMaximaChargeShared;   //!
   TBranch        *b_ep1_btMaximaChargeSharedFragment;   //!
   TBranch        *b_ep1_btMaximaShared;   //!
   TBranch        *b_ep1_btMaximaSharedFragment;   //!
   TBranch        *b_ep1_btMeanDist;   //!
   TBranch        *b_ep1_btNearbyMaxima;   //!
   TBranch        *b_ep1_btNearbyMaximaShared;   //!
   TBranch        *b_ep1_btPadsClus;   //!
   TBranch        *b_ep1_btPadsRing;   //!
   TBranch        *b_ep1_btRingMatrix;   //!
   TBranch        *b_ep1_dedx_mdc;   //!
   TBranch        *b_ep1_dedx_tof;   //!
   TBranch        *b_ep1_id;   //!
   TBranch        *b_ep1_isBT;   //!
   TBranch        *b_ep1_isOffVertexClust;   //!
   TBranch        *b_ep1_isPrimaryVertex;   //!
   TBranch        *b_ep1_isUsedVertex;   //!
   TBranch        *b_ep1_isring;   //!
   TBranch        *b_ep1_isringmdc;   //!
   TBranch        *b_ep1_isringnomatch;   //!
   TBranch        *b_ep1_isringtrack;   //!
   TBranch        *b_ep1_kIsLepton;   //!
   TBranch        *b_ep1_kIsUsed;   //!
   TBranch        *b_ep1_mdcinnerchi2;   //!
   TBranch        *b_ep1_mdcouterchi2;   //!
   TBranch        *b_ep1_oa_hadr;   //!
   TBranch        *b_ep1_oa_lept;   //!
   TBranch        *b_ep1_p;   //!
   TBranch        *b_ep1_p_corr_ep2;   //!
   TBranch        *b_ep1_p_corr_ep1;   //!
   TBranch        *b_ep1_p_corr_p;   //!
   TBranch        *b_ep1_p_corr_pim;   //!
   TBranch        *b_ep1_p_corr_pip;   //!
   TBranch        *b_ep1_phi;   //!
   TBranch        *b_ep1_phi_rich;   //!
   TBranch        *b_ep1_pid;   //!
   TBranch        *b_ep1_q;   //!
   TBranch        *b_ep1_r;   //!
   TBranch        *b_ep1_resolution;   //!
   TBranch        *b_ep1_resoultion;   //!
   TBranch        *b_ep1_rich_amp;   //!
   TBranch        *b_ep1_rich_centr;   //!
   TBranch        *b_ep1_rich_houtra;   //!
   TBranch        *b_ep1_rich_padnum;   //!
   TBranch        *b_ep1_rich_patmat;   //!
   TBranch        *b_ep1_rkchi2;   //!
   TBranch        *b_ep1_sector;   //!
   TBranch        *b_ep1_shw_sum0;   //!
   TBranch        *b_ep1_shw_sum1;   //!
   TBranch        *b_ep1_shw_sum2;   //!
   TBranch        *b_ep1_system;   //!
   TBranch        *b_ep1_theta;   //!
   TBranch        *b_ep1_theta_rich;   //!
   TBranch        *b_ep1_tof_mom;   //!
   TBranch        *b_ep1_tof_new;   //!
   TBranch        *b_ep1_tof_rec;   //!
   TBranch        *b_ep1_track_length;   //!
   TBranch        *b_ep1_tracklength;   //!
   TBranch        *b_ep1_z;   //!
   TBranch        *b_event;   //!
   TBranch        *b_fw_beta_1;   //!
   TBranch        *b_fw_beta_2;   //!
   TBranch        *b_fw_beta_3;   //!
   TBranch        *b_fw_charge_1;   //!
   TBranch        *b_fw_charge_2;   //!
   TBranch        *b_fw_charge_3;   //!
   TBranch        *b_fw_cluster_mult;   //!
   TBranch        *b_fw_distance_1;   //!
   TBranch        *b_fw_distance_2;   //!
   TBranch        *b_fw_distance_3;   //!
   TBranch        *b_fw_mult;   //!
   TBranch        *b_fw_p_1;   //!
   TBranch        *b_fw_p_2;   //!
   TBranch        *b_fw_p_3;   //!
   TBranch        *b_fw_phi_1;   //!
   TBranch        *b_fw_phi_2;   //!
   TBranch        *b_fw_phi_3;   //!
   TBranch        *b_fw_size_1;   //!
   TBranch        *b_fw_size_2;   //!
   TBranch        *b_fw_size_3;   //!
   TBranch        *b_fw_spectator_1;   //!
   TBranch        *b_fw_spectator_2;   //!
   TBranch        *b_fw_spectator_3;   //!
   TBranch        *b_fw_theta_1;   //!
   TBranch        *b_fw_theta_2;   //!
   TBranch        *b_fw_theta_3;   //!
   TBranch        *b_fw_time_1;   //!
   TBranch        *b_fw_time_2;   //!
   TBranch        *b_fw_time_3;   //!
   TBranch        *b_fw_time_min_1;   //!
   TBranch        *b_fw_time_min_2;   //!
   TBranch        *b_fw_time_min_3;   //!
   TBranch        *b_fw_x_lab_1;   //!
   TBranch        *b_fw_x_lab_2;   //!
   TBranch        *b_fw_x_lab_3;   //!
   TBranch        *b_fw_y_lab_1;   //!
   TBranch        *b_fw_y_lab_2;   //!
   TBranch        *b_fw_y_lab_3;   //!
   TBranch        *b_fw_z_lab_1;   //!
   TBranch        *b_fw_z_lab_2;   //!
   TBranch        *b_fw_z_lab_3;   //!
   TBranch        *b_hneg_mult;   //!
   TBranch        *b_hpos_mult;   //!
   TBranch        *b_isBest;   //!
   TBranch        *b_lneg_mult;   //!
   TBranch        *b_lpos_mult;   //!
   TBranch        *b_p_beta;   //!
   TBranch        *b_p_beta_new;   //!
   TBranch        *b_p_btChargeRing;   //!
   TBranch        *b_p_btChargeSum;   //!
   TBranch        *b_p_btChi2;   //!
   TBranch        *b_p_btClusters;   //!
   TBranch        *b_p_btMaxima;   //!
   TBranch        *b_p_btMaximaCharge;   //!
   TBranch        *b_p_btMaximaChargeShared;   //!
   TBranch        *b_p_btMaximaChargeSharedFragment;   //!
   TBranch        *b_p_btMaximaShared;   //!
   TBranch        *b_p_btMaximaSharedFragment;   //!
   TBranch        *b_p_btMeanDist;   //!
   TBranch        *b_p_btNearbyMaxima;   //!
   TBranch        *b_p_btNearbyMaximaShared;   //!
   TBranch        *b_p_btPadsClus;   //!
   TBranch        *b_p_btPadsRing;   //!
   TBranch        *b_p_btRingMatrix;   //!
   TBranch        *b_p_dedx_mdc;   //!
   TBranch        *b_p_dedx_tof;   //!
   TBranch        *b_p_id;   //!
   TBranch        *b_p_isBT;   //!
   TBranch        *b_p_isOffVertexClust;   //!
   TBranch        *b_p_isPrimaryVertex;   //!
   TBranch        *b_p_isUsedVertex;   //!
   TBranch        *b_p_isring;   //!
   TBranch        *b_p_isringmdc;   //!
   TBranch        *b_p_isringnomatch;   //!
   TBranch        *b_p_isringtrack;   //!
   TBranch        *b_p_kIsLepton;   //!
   TBranch        *b_p_kIsUsed;   //!
   TBranch        *b_p_mdcinnerchi2;   //!
   TBranch        *b_p_mdcouterchi2;   //!
   TBranch        *b_p_oa_hadr;   //!
   TBranch        *b_p_oa_lept;   //!
   TBranch        *b_p_p;   //!
   TBranch        *b_p_p_corr_ep2;   //!
   TBranch        *b_p_p_corr_ep1;   //!
   TBranch        *b_p_p_corr_p;   //!
   TBranch        *b_p_p_corr_pim;   //!
   TBranch        *b_p_p_corr_pip;   //!
   TBranch        *b_p_phi;   //!
   TBranch        *b_p_phi_rich;   //!
   TBranch        *b_p_pid;   //!
   TBranch        *b_p_q;   //!
   TBranch        *b_p_r;   //!
   TBranch        *b_p_resolution;   //!
   TBranch        *b_p_resoultion;   //!
   TBranch        *b_p_rich_amp;   //!
   TBranch        *b_p_rich_centr;   //!
   TBranch        *b_p_rich_houtra;   //!
   TBranch        *b_p_rich_padnum;   //!
   TBranch        *b_p_rich_patmat;   //!
   TBranch        *b_p_rkchi2;   //!
   TBranch        *b_p_sector;   //!
   TBranch        *b_p_shw_sum0;   //!
   TBranch        *b_p_shw_sum1;   //!
   TBranch        *b_p_shw_sum2;   //!
   TBranch        *b_p_system;   //!
   TBranch        *b_p_theta;   //!
   TBranch        *b_p_theta_rich;   //!
   TBranch        *b_p_tof_mom;   //!
   TBranch        *b_p_tof_new;   //!
   TBranch        *b_p_tof_rec;   //!
   TBranch        *b_p_track_length;   //!
   TBranch        *b_p_tracklength;   //!
   TBranch        *b_p_z;   //!

TBranch        *b_pim_beta;   //!
   TBranch        *b_pim_beta_new;   //!
   TBranch        *b_pim_btChargeRing;   //!
   TBranch        *b_pim_btChargeSum;   //!
   TBranch        *b_pim_btChi2;   //!
   TBranch        *b_pim_btClusters;   //!
   TBranch        *b_pim_btMaxima;   //!
   TBranch        *b_pim_btMaximaCharge;   //!
   TBranch        *b_pim_btMaximaChargeShared;   //!
   TBranch        *b_pim_btMaximaChargeSharedFragment;   //!
   TBranch        *b_pim_btMaximaShared;   //!
   TBranch        *b_pim_btMaximaSharedFragment;   //!
   TBranch        *b_pim_btMeanDist;   //!
   TBranch        *b_pim_btNearbyMaxima;   //!
   TBranch        *b_pim_btNearbyMaximaShared;   //!
   TBranch        *b_pim_btPadsClus;   //!
   TBranch        *b_pim_btPadsRing;   //!
   TBranch        *b_pim_btRingMatrix;   //!
   TBranch        *b_pim_dedx_mdc;   //!
   TBranch        *b_pim_dedx_tof;   //!
   TBranch        *b_pim_id;   //!
   TBranch        *b_pim_isBT;   //!
   TBranch        *b_pim_isOffVertexClust;   //!
   TBranch        *b_pim_isPrimaryVertex;   //!
   TBranch        *b_pim_isUsedVertex;   //!
   TBranch        *b_pim_isring;   //!
   TBranch        *b_pim_isringmdc;   //!
   TBranch        *b_pim_isringnomatch;   //!
   TBranch        *b_pim_isringtrack;   //!
   TBranch        *b_pim_kIsLepton;   //!
   TBranch        *b_pim_kIsUsed;   //!
   TBranch        *b_pim_mdcinnerchi2;   //!
   TBranch        *b_pim_mdcouterchi2;   //!
   TBranch        *b_pim_oa_hadr;   //!
   TBranch        *b_pim_oa_lept;   //!
   TBranch        *b_pim_p;   //!
   TBranch        *b_pim_p_corr_ep2;   //!
   TBranch        *b_pim_p_corr_ep1;   //!
   TBranch        *b_pim_p_corr_p;   //!
   TBranch        *b_pim_p_corr_pim;   //!
   TBranch        *b_pim_p_corr_pip;   //!
   TBranch        *b_pim_phi;   //!
   TBranch        *b_pim_phi_rich;   //!
   TBranch        *b_pim_pid;   //!
   TBranch        *b_pim_q;   //!
   TBranch        *b_pim_r;   //!
   TBranch        *b_pim_resolution;   //!
   TBranch        *b_pim_resoultion;   //!
   TBranch        *b_pim_rich_amp;   //!
   TBranch        *b_pim_rich_centr;   //!
   TBranch        *b_pim_rich_houtra;   //!
   TBranch        *b_pim_rich_padnum;   //!
   TBranch        *b_pim_rich_patmat;   //!
   TBranch        *b_pim_rkchi2;   //!
   TBranch        *b_pim_sector;   //!
   TBranch        *b_pim_shw_sum0;   //!
   TBranch        *b_pim_shw_sum1;   //!
   TBranch        *b_pim_shw_sum2;   //!
   TBranch        *b_pim_system;   //!
   TBranch        *b_pim_theta;   //!
   TBranch        *b_pim_theta_rich;   //!
   TBranch        *b_pim_tof_mom;   //!
   TBranch        *b_pim_tof_new;   //!
   TBranch        *b_pim_tof_rec;   //!
   TBranch        *b_pim_track_length;   //!
   TBranch        *b_pim_tracklength;   //!
   TBranch        *b_pim_z;   //!
  
  TBranch        *b_runnumber;   //!
   TBranch        *b_totalmult;   //!
   TBranch        *b_trigbit;   //!
   TBranch        *b_trigdec;   //!
   TBranch        *b_trigdownscale;   //!
   TBranch        *b_trigdownscaleflag;   //!

   PPimEpEp(TTree *tree=0);
   virtual ~PPimEpEp();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif
