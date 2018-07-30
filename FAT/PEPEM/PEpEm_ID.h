//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Thu Apr  5 12:15:37 2018 by ROOT version 5.34/34
// from TTree PEpEm_ID/PEpEm_ID
// found on file: /lustre/nyx/hades/user/knowakow/PP/PAT/FILES/full_stat/hadrons.root
//////////////////////////////////////////////////////////

#ifndef PEpEm_ID_h
#define PEpEm_ID_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.

// Fixed size dimensions of array or collections stored in the TTree if any.

class PEpEm_ID {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

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
   Float_t         em_beta;
   Float_t         em_beta_new;
   Float_t         em_btChargeRing;
   Float_t         em_btChargeSum;
   Float_t         em_btChi2;
   Float_t         em_btClusters;
   Float_t         em_btMaxima;
   Float_t         em_btMaximaCharge;
   Float_t         em_btMaximaChargeShared;
   Float_t         em_btMaximaChargeSharedFragment;
   Float_t         em_btMaximaShared;
   Float_t         em_btMaximaSharedFragment;
   Float_t         em_btMeanDist;
   Float_t         em_btNearbyMaxima;
   Float_t         em_btNearbyMaximaShared;
   Float_t         em_btPadsClus;
   Float_t         em_btPadsRing;
   Float_t         em_btRingMatrix;
   Float_t         em_dedx_mdc;
   Float_t         em_dedx_tof;
   Float_t         em_id;
   Float_t         em_isBT;
   Float_t         em_isOffVertexClust;
   Float_t         em_isPrimaryVertex;
   Float_t         em_isUsedVertex;
   Float_t         em_isring;
   Float_t         em_isringmdc;
   Float_t         em_isringnomatch;
   Float_t         em_isringtrack;
   Float_t         em_kIsLepton;
   Float_t         em_kIsUsed;
   Float_t         em_mdcinnerchi2;
   Float_t         em_mdcouterchi2;
   Float_t         em_oa_hadr;
   Float_t         em_oa_lept;
   Float_t         em_p;
   Float_t         em_p_corr_em;
   Float_t         em_p_corr_ep;
   Float_t         em_p_corr_p;
   Float_t         em_p_corr_pim;
   Float_t         em_p_corr_pip;
   Float_t         em_phi;
   Float_t         em_phi_rich;
   Float_t         em_pid;
   Float_t         em_q;
   Float_t         em_r;
   Float_t         em_resolution;
   Float_t         em_resoultion;
   Float_t         em_rich_amp;
   Float_t         em_rich_centr;
   Float_t         em_rich_houtra;
   Float_t         em_rich_padnum;
   Float_t         em_rich_patmat;
   Float_t         em_rkchi2;
   Float_t         em_sector;
   Float_t         em_shw_sum0;
   Float_t         em_shw_sum1;
   Float_t         em_shw_sum2;
   Float_t         em_system;
   Float_t         em_theta;
   Float_t         em_theta_rich;
   Float_t         em_tof_mom;
   Float_t         em_tof_new;
   Float_t         em_tof_rec;
   Float_t         em_track_length;
   Float_t         em_tracklength;
   Float_t         em_z;
   Float_t         ep_beta;
   Float_t         ep_beta_new;
   Float_t         ep_btChargeRing;
   Float_t         ep_btChargeSum;
   Float_t         ep_btChi2;
   Float_t         ep_btClusters;
   Float_t         ep_btMaxima;
   Float_t         ep_btMaximaCharge;
   Float_t         ep_btMaximaChargeShared;
   Float_t         ep_btMaximaChargeSharedFragment;
   Float_t         ep_btMaximaShared;
   Float_t         ep_btMaximaSharedFragment;
   Float_t         ep_btMeanDist;
   Float_t         ep_btNearbyMaxima;
   Float_t         ep_btNearbyMaximaShared;
   Float_t         ep_btPadsClus;
   Float_t         ep_btPadsRing;
   Float_t         ep_btRingMatrix;
   Float_t         ep_dedx_mdc;
   Float_t         ep_dedx_tof;
   Float_t         ep_id;
   Float_t         ep_isBT;
   Float_t         ep_isOffVertexClust;
   Float_t         ep_isPrimaryVertex;
   Float_t         ep_isUsedVertex;
   Float_t         ep_isring;
   Float_t         ep_isringmdc;
   Float_t         ep_isringnomatch;
   Float_t         ep_isringtrack;
   Float_t         ep_kIsLepton;
   Float_t         ep_kIsUsed;
   Float_t         ep_mdcinnerchi2;
   Float_t         ep_mdcouterchi2;
   Float_t         ep_oa_hadr;
   Float_t         ep_oa_lept;
   Float_t         ep_p;
   Float_t         ep_p_corr_em;
   Float_t         ep_p_corr_ep;
   Float_t         ep_p_corr_p;
   Float_t         ep_p_corr_pim;
   Float_t         ep_p_corr_pip;
   Float_t         ep_phi;
   Float_t         ep_phi_rich;
   Float_t         ep_pid;
   Float_t         ep_q;
   Float_t         ep_r;
   Float_t         ep_resolution;
   Float_t         ep_resoultion;
   Float_t         ep_rich_amp;
   Float_t         ep_rich_centr;
   Float_t         ep_rich_houtra;
   Float_t         ep_rich_padnum;
   Float_t         ep_rich_patmat;
   Float_t         ep_rkchi2;
   Float_t         ep_sector;
   Float_t         ep_shw_sum0;
   Float_t         ep_shw_sum1;
   Float_t         ep_shw_sum2;
   Float_t         ep_system;
   Float_t         ep_theta;
   Float_t         ep_theta_rich;
   Float_t         ep_tof_mom;
   Float_t         ep_tof_new;
   Float_t         ep_tof_rec;
   Float_t         ep_track_length;
   Float_t         ep_tracklength;
   Float_t         ep_z;
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
   Float_t         p_p_corr_em;
   Float_t         p_p_corr_ep;
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
   TBranch        *b_em_beta;   //!
   TBranch        *b_em_beta_new;   //!
   TBranch        *b_em_btChargeRing;   //!
   TBranch        *b_em_btChargeSum;   //!
   TBranch        *b_em_btChi2;   //!
   TBranch        *b_em_btClusters;   //!
   TBranch        *b_em_btMaxima;   //!
   TBranch        *b_em_btMaximaCharge;   //!
   TBranch        *b_em_btMaximaChargeShared;   //!
   TBranch        *b_em_btMaximaChargeSharedFragment;   //!
   TBranch        *b_em_btMaximaShared;   //!
   TBranch        *b_em_btMaximaSharedFragment;   //!
   TBranch        *b_em_btMeanDist;   //!
   TBranch        *b_em_btNearbyMaxima;   //!
   TBranch        *b_em_btNearbyMaximaShared;   //!
   TBranch        *b_em_btPadsClus;   //!
   TBranch        *b_em_btPadsRing;   //!
   TBranch        *b_em_btRingMatrix;   //!
   TBranch        *b_em_dedx_mdc;   //!
   TBranch        *b_em_dedx_tof;   //!
   TBranch        *b_em_id;   //!
   TBranch        *b_em_isBT;   //!
   TBranch        *b_em_isOffVertexClust;   //!
   TBranch        *b_em_isPrimaryVertex;   //!
   TBranch        *b_em_isUsedVertex;   //!
   TBranch        *b_em_isring;   //!
   TBranch        *b_em_isringmdc;   //!
   TBranch        *b_em_isringnomatch;   //!
   TBranch        *b_em_isringtrack;   //!
   TBranch        *b_em_kIsLepton;   //!
   TBranch        *b_em_kIsUsed;   //!
   TBranch        *b_em_mdcinnerchi2;   //!
   TBranch        *b_em_mdcouterchi2;   //!
   TBranch        *b_em_oa_hadr;   //!
   TBranch        *b_em_oa_lept;   //!
   TBranch        *b_em_p;   //!
   TBranch        *b_em_p_corr_em;   //!
   TBranch        *b_em_p_corr_ep;   //!
   TBranch        *b_em_p_corr_p;   //!
   TBranch        *b_em_p_corr_pim;   //!
   TBranch        *b_em_p_corr_pip;   //!
   TBranch        *b_em_phi;   //!
   TBranch        *b_em_phi_rich;   //!
   TBranch        *b_em_pid;   //!
   TBranch        *b_em_q;   //!
   TBranch        *b_em_r;   //!
   TBranch        *b_em_resolution;   //!
   TBranch        *b_em_resoultion;   //!
   TBranch        *b_em_rich_amp;   //!
   TBranch        *b_em_rich_centr;   //!
   TBranch        *b_em_rich_houtra;   //!
   TBranch        *b_em_rich_padnum;   //!
   TBranch        *b_em_rich_patmat;   //!
   TBranch        *b_em_rkchi2;   //!
   TBranch        *b_em_sector;   //!
   TBranch        *b_em_shw_sum0;   //!
   TBranch        *b_em_shw_sum1;   //!
   TBranch        *b_em_shw_sum2;   //!
   TBranch        *b_em_system;   //!
   TBranch        *b_em_theta;   //!
   TBranch        *b_em_theta_rich;   //!
   TBranch        *b_em_tof_mom;   //!
   TBranch        *b_em_tof_new;   //!
   TBranch        *b_em_tof_rec;   //!
   TBranch        *b_em_track_length;   //!
   TBranch        *b_em_tracklength;   //!
   TBranch        *b_em_z;   //!
   TBranch        *b_ep_beta;   //!
   TBranch        *b_ep_beta_new;   //!
   TBranch        *b_ep_btChargeRing;   //!
   TBranch        *b_ep_btChargeSum;   //!
   TBranch        *b_ep_btChi2;   //!
   TBranch        *b_ep_btClusters;   //!
   TBranch        *b_ep_btMaxima;   //!
   TBranch        *b_ep_btMaximaCharge;   //!
   TBranch        *b_ep_btMaximaChargeShared;   //!
   TBranch        *b_ep_btMaximaChargeSharedFragment;   //!
   TBranch        *b_ep_btMaximaShared;   //!
   TBranch        *b_ep_btMaximaSharedFragment;   //!
   TBranch        *b_ep_btMeanDist;   //!
   TBranch        *b_ep_btNearbyMaxima;   //!
   TBranch        *b_ep_btNearbyMaximaShared;   //!
   TBranch        *b_ep_btPadsClus;   //!
   TBranch        *b_ep_btPadsRing;   //!
   TBranch        *b_ep_btRingMatrix;   //!
   TBranch        *b_ep_dedx_mdc;   //!
   TBranch        *b_ep_dedx_tof;   //!
   TBranch        *b_ep_id;   //!
   TBranch        *b_ep_isBT;   //!
   TBranch        *b_ep_isOffVertexClust;   //!
   TBranch        *b_ep_isPrimaryVertex;   //!
   TBranch        *b_ep_isUsedVertex;   //!
   TBranch        *b_ep_isring;   //!
   TBranch        *b_ep_isringmdc;   //!
   TBranch        *b_ep_isringnomatch;   //!
   TBranch        *b_ep_isringtrack;   //!
   TBranch        *b_ep_kIsLepton;   //!
   TBranch        *b_ep_kIsUsed;   //!
   TBranch        *b_ep_mdcinnerchi2;   //!
   TBranch        *b_ep_mdcouterchi2;   //!
   TBranch        *b_ep_oa_hadr;   //!
   TBranch        *b_ep_oa_lept;   //!
   TBranch        *b_ep_p;   //!
   TBranch        *b_ep_p_corr_em;   //!
   TBranch        *b_ep_p_corr_ep;   //!
   TBranch        *b_ep_p_corr_p;   //!
   TBranch        *b_ep_p_corr_pim;   //!
   TBranch        *b_ep_p_corr_pip;   //!
   TBranch        *b_ep_phi;   //!
   TBranch        *b_ep_phi_rich;   //!
   TBranch        *b_ep_pid;   //!
   TBranch        *b_ep_q;   //!
   TBranch        *b_ep_r;   //!
   TBranch        *b_ep_resolution;   //!
   TBranch        *b_ep_resoultion;   //!
   TBranch        *b_ep_rich_amp;   //!
   TBranch        *b_ep_rich_centr;   //!
   TBranch        *b_ep_rich_houtra;   //!
   TBranch        *b_ep_rich_padnum;   //!
   TBranch        *b_ep_rich_patmat;   //!
   TBranch        *b_ep_rkchi2;   //!
   TBranch        *b_ep_sector;   //!
   TBranch        *b_ep_shw_sum0;   //!
   TBranch        *b_ep_shw_sum1;   //!
   TBranch        *b_ep_shw_sum2;   //!
   TBranch        *b_ep_system;   //!
   TBranch        *b_ep_theta;   //!
   TBranch        *b_ep_theta_rich;   //!
   TBranch        *b_ep_tof_mom;   //!
   TBranch        *b_ep_tof_new;   //!
   TBranch        *b_ep_tof_rec;   //!
   TBranch        *b_ep_track_length;   //!
   TBranch        *b_ep_tracklength;   //!
   TBranch        *b_ep_z;   //!
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
   TBranch        *b_p_p_corr_em;   //!
   TBranch        *b_p_p_corr_ep;   //!
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
   TBranch        *b_runnumber;   //!
   TBranch        *b_totalmult;   //!
   TBranch        *b_trigbit;   //!
   TBranch        *b_trigdec;   //!
   TBranch        *b_trigdownscale;   //!
   TBranch        *b_trigdownscaleflag;   //!

   PEpEm_ID(TTree *tree=0);
   virtual ~PEpEm_ID();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef PEpEm_ID_cxx
PEpEm_ID::PEpEm_ID(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("/lustre/nyx/hades/user/knowakow/PP/PAT/FILES/full_stat/hadrons.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("/lustre/nyx/hades/user/knowakow/PP/PAT/FILES/full_stat/hadrons.root");
      }
      f->GetObject("PEpEm_ID",tree);

   }
   Init(tree);
}

PEpEm_ID::~PEpEm_ID()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t PEpEm_ID::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t PEpEm_ID::LoadTree(Long64_t entry)
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

void PEpEm_ID::Init(TTree *tree)
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

Bool_t PEpEm_ID::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void PEpEm_ID::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t PEpEm_ID::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef PEpEm_ID_cxx
