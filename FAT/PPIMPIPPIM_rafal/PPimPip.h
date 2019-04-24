//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Wed Apr 24 14:57:29 2019 by ROOT version 5.34/34
// from TTree PPimPip/PPimPip
// found on file: all.root
//////////////////////////////////////////////////////////

#ifndef PPimPip_h
#define PPimPip_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.

// Fixed size dimensions of array or collections stored in the TTree if any.

class PPimPip {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types
   Float_t         eVert_chi2;
   Float_t         eVert_x;
   Float_t         eVert_y;
   Float_t         eVert_z;
   Float_t         event;
   Float_t         hneg_mult;
   Float_t         hpos_mult;
   Float_t         isBest;
   Float_t         lneg_mult;
   Float_t         lpos_mult;
   Float_t         p_beta;
   Float_t         p_beta_new;
   Float_t         p_dedx_in;
   Float_t         p_dedx_in_sigma;
   Float_t         p_dedx_mdc;
   Float_t         p_dedx_mdc_sigma;
   Float_t         p_dedx_out;
   Float_t         p_dedx_out_sigma;
   Float_t         p_dedx_tof;
   Float_t         p_id;
   Float_t         p_isring;
   Float_t         p_kIsLepton;
   Float_t         p_kIsUsed;
   Float_t         p_mdcchi2;
   Float_t         p_oa_hadr;
   Float_t         p_oa_lept;
   Float_t         p_p;
   Float_t         p_phi;
   Float_t         p_q;
   Float_t         p_r;
   Float_t         p_resolution;
   Float_t         p_rkchi2;
   Float_t         p_sector;
   Float_t         p_shw_sum0;
   Float_t         p_shw_sum1;
   Float_t         p_shw_sum2;
   Float_t         p_sim_corrflag;
   Float_t         p_sim_geninfo;
   Float_t         p_sim_geninfo1;
   Float_t         p_sim_geninfo2;
   Float_t         p_sim_genweight;
   Float_t         p_sim_id;
   Float_t         p_sim_iscommon;
   Float_t         p_sim_mediumid;
   Float_t         p_sim_p;
   Float_t         p_sim_parentid;
   Float_t         p_sim_primaryflag;
   Float_t         p_sim_processid;
   Float_t         p_sim_px;
   Float_t         p_sim_py;
   Float_t         p_sim_pz;
   Float_t         p_sim_vertexx;
   Float_t         p_sim_vertexy;
   Float_t         p_sim_vertexz;
   Float_t         p_system;
   Float_t         p_theta;
   Float_t         p_tof_exp;
   Float_t         p_tof_mom;
   Float_t         p_tof_new;
   Float_t         p_tofino_mult;
   Float_t         p_track_length;
   Float_t         p_z;
   Float_t         pim_beta;
   Float_t         pim_beta_new;
   Float_t         pim_dedx_in;
   Float_t         pim_dedx_in_sigma;
   Float_t         pim_dedx_mdc;
   Float_t         pim_dedx_mdc_sigma;
   Float_t         pim_dedx_out;
   Float_t         pim_dedx_out_sigma;
   Float_t         pim_dedx_tof;
   Float_t         pim_id;
   Float_t         pim_isring;
   Float_t         pim_kIsLepton;
   Float_t         pim_kIsUsed;
   Float_t         pim_mdcchi2;
   Float_t         pim_oa_hadr;
   Float_t         pim_oa_lept;
   Float_t         pim_p;
   Float_t         pim_phi;
   Float_t         pim_q;
   Float_t         pim_r;
   Float_t         pim_resolution;
   Float_t         pim_rkchi2;
   Float_t         pim_sector;
   Float_t         pim_shw_sum0;
   Float_t         pim_shw_sum1;
   Float_t         pim_shw_sum2;
   Float_t         pim_sim_corrflag;
   Float_t         pim_sim_geninfo;
   Float_t         pim_sim_geninfo1;
   Float_t         pim_sim_geninfo2;
   Float_t         pim_sim_genweight;
   Float_t         pim_sim_id;
   Float_t         pim_sim_iscommon;
   Float_t         pim_sim_mediumid;
   Float_t         pim_sim_p;
   Float_t         pim_sim_parentid;
   Float_t         pim_sim_primaryflag;
   Float_t         pim_sim_processid;
   Float_t         pim_sim_px;
   Float_t         pim_sim_py;
   Float_t         pim_sim_pz;
   Float_t         pim_sim_vertexx;
   Float_t         pim_sim_vertexy;
   Float_t         pim_sim_vertexz;
   Float_t         pim_system;
   Float_t         pim_theta;
   Float_t         pim_tof_exp;
   Float_t         pim_tof_mom;
   Float_t         pim_tof_new;
   Float_t         pim_tofino_mult;
   Float_t         pim_track_length;
   Float_t         pim_z;
   Float_t         pip_beta;
   Float_t         pip_beta_new;
   Float_t         pip_dedx_in;
   Float_t         pip_dedx_in_sigma;
   Float_t         pip_dedx_mdc;
   Float_t         pip_dedx_mdc_sigma;
   Float_t         pip_dedx_out;
   Float_t         pip_dedx_out_sigma;
   Float_t         pip_dedx_tof;
   Float_t         pip_id;
   Float_t         pip_isring;
   Float_t         pip_kIsLepton;
   Float_t         pip_kIsUsed;
   Float_t         pip_mdcchi2;
   Float_t         pip_oa_hadr;
   Float_t         pip_oa_lept;
   Float_t         pip_p;
   Float_t         pip_phi;
   Float_t         pip_q;
   Float_t         pip_r;
   Float_t         pip_resolution;
   Float_t         pip_rkchi2;
   Float_t         pip_sector;
   Float_t         pip_shw_sum0;
   Float_t         pip_shw_sum1;
   Float_t         pip_shw_sum2;
   Float_t         pip_sim_corrflag;
   Float_t         pip_sim_geninfo;
   Float_t         pip_sim_geninfo1;
   Float_t         pip_sim_geninfo2;
   Float_t         pip_sim_genweight;
   Float_t         pip_sim_id;
   Float_t         pip_sim_iscommon;
   Float_t         pip_sim_mediumid;
   Float_t         pip_sim_p;
   Float_t         pip_sim_parentid;
   Float_t         pip_sim_primaryflag;
   Float_t         pip_sim_processid;
   Float_t         pip_sim_px;
   Float_t         pip_sim_py;
   Float_t         pip_sim_pz;
   Float_t         pip_sim_vertexx;
   Float_t         pip_sim_vertexy;
   Float_t         pip_sim_vertexz;
   Float_t         pip_system;
   Float_t         pip_theta;
   Float_t         pip_tof_exp;
   Float_t         pip_tof_mom;
   Float_t         pip_tof_new;
   Float_t         pip_tofino_mult;
   Float_t         pip_track_length;
   Float_t         pip_z;
   Float_t         runnumber;
   Float_t         totalmult;
   Float_t         trigbit;
   Float_t         trigdec;
   Float_t         trigdownscale;
   Float_t         trigdownscaleflag;

   // List of branches
   TBranch        *b_eVert_chi2;   //!
   TBranch        *b_eVert_x;   //!
   TBranch        *b_eVert_y;   //!
   TBranch        *b_eVert_z;   //!
   TBranch        *b_event;   //!
   TBranch        *b_hneg_mult;   //!
   TBranch        *b_hpos_mult;   //!
   TBranch        *b_isBest;   //!
   TBranch        *b_lneg_mult;   //!
   TBranch        *b_lpos_mult;   //!
   TBranch        *b_p_beta;   //!
   TBranch        *b_p_beta_new;   //!
   TBranch        *b_p_dedx_in;   //!
   TBranch        *b_p_dedx_in_sigma;   //!
   TBranch        *b_p_dedx_mdc;   //!
   TBranch        *b_p_dedx_mdc_sigma;   //!
   TBranch        *b_p_dedx_out;   //!
   TBranch        *b_p_dedx_out_sigma;   //!
   TBranch        *b_p_dedx_tof;   //!
   TBranch        *b_p_id;   //!
   TBranch        *b_p_isring;   //!
   TBranch        *b_p_kIsLepton;   //!
   TBranch        *b_p_kIsUsed;   //!
   TBranch        *b_p_mdcchi2;   //!
   TBranch        *b_p_oa_hadr;   //!
   TBranch        *b_p_oa_lept;   //!
   TBranch        *b_p_p;   //!
   TBranch        *b_p_phi;   //!
   TBranch        *b_p_q;   //!
   TBranch        *b_p_r;   //!
   TBranch        *b_p_resolution;   //!
   TBranch        *b_p_rkchi2;   //!
   TBranch        *b_p_sector;   //!
   TBranch        *b_p_shw_sum0;   //!
   TBranch        *b_p_shw_sum1;   //!
   TBranch        *b_p_shw_sum2;   //!
   TBranch        *b_p_sim_corrflag;   //!
   TBranch        *b_p_sim_geninfo;   //!
   TBranch        *b_p_sim_geninfo1;   //!
   TBranch        *b_p_sim_geninfo2;   //!
   TBranch        *b_p_sim_genweight;   //!
   TBranch        *b_p_sim_id;   //!
   TBranch        *b_p_sim_iscommon;   //!
   TBranch        *b_p_sim_mediumid;   //!
   TBranch        *b_p_sim_p;   //!
   TBranch        *b_p_sim_parentid;   //!
   TBranch        *b_p_sim_primaryflag;   //!
   TBranch        *b_p_sim_processid;   //!
   TBranch        *b_p_sim_px;   //!
   TBranch        *b_p_sim_py;   //!
   TBranch        *b_p_sim_pz;   //!
   TBranch        *b_p_sim_vertexx;   //!
   TBranch        *b_p_sim_vertexy;   //!
   TBranch        *b_p_sim_vertexz;   //!
   TBranch        *b_p_system;   //!
   TBranch        *b_p_theta;   //!
   TBranch        *b_p_tof_exp;   //!
   TBranch        *b_p_tof_mom;   //!
   TBranch        *b_p_tof_new;   //!
   TBranch        *b_p_tofino_mult;   //!
   TBranch        *b_p_track_length;   //!
   TBranch        *b_p_z;   //!
   TBranch        *b_pim_beta;   //!
   TBranch        *b_pim_beta_new;   //!
   TBranch        *b_pim_dedx_in;   //!
   TBranch        *b_pim_dedx_in_sigma;   //!
   TBranch        *b_pim_dedx_mdc;   //!
   TBranch        *b_pim_dedx_mdc_sigma;   //!
   TBranch        *b_pim_dedx_out;   //!
   TBranch        *b_pim_dedx_out_sigma;   //!
   TBranch        *b_pim_dedx_tof;   //!
   TBranch        *b_pim_id;   //!
   TBranch        *b_pim_isring;   //!
   TBranch        *b_pim_kIsLepton;   //!
   TBranch        *b_pim_kIsUsed;   //!
   TBranch        *b_pim_mdcchi2;   //!
   TBranch        *b_pim_oa_hadr;   //!
   TBranch        *b_pim_oa_lept;   //!
   TBranch        *b_pim_p;   //!
   TBranch        *b_pim_phi;   //!
   TBranch        *b_pim_q;   //!
   TBranch        *b_pim_r;   //!
   TBranch        *b_pim_resolution;   //!
   TBranch        *b_pim_rkchi2;   //!
   TBranch        *b_pim_sector;   //!
   TBranch        *b_pim_shw_sum0;   //!
   TBranch        *b_pim_shw_sum1;   //!
   TBranch        *b_pim_shw_sum2;   //!
   TBranch        *b_pim_sim_corrflag;   //!
   TBranch        *b_pim_sim_geninfo;   //!
   TBranch        *b_pim_sim_geninfo1;   //!
   TBranch        *b_pim_sim_geninfo2;   //!
   TBranch        *b_pim_sim_genweight;   //!
   TBranch        *b_pim_sim_id;   //!
   TBranch        *b_pim_sim_iscommon;   //!
   TBranch        *b_pim_sim_mediumid;   //!
   TBranch        *b_pim_sim_p;   //!
   TBranch        *b_pim_sim_parentid;   //!
   TBranch        *b_pim_sim_primaryflag;   //!
   TBranch        *b_pim_sim_processid;   //!
   TBranch        *b_pim_sim_px;   //!
   TBranch        *b_pim_sim_py;   //!
   TBranch        *b_pim_sim_pz;   //!
   TBranch        *b_pim_sim_vertexx;   //!
   TBranch        *b_pim_sim_vertexy;   //!
   TBranch        *b_pim_sim_vertexz;   //!
   TBranch        *b_pim_system;   //!
   TBranch        *b_pim_theta;   //!
   TBranch        *b_pim_tof_exp;   //!
   TBranch        *b_pim_tof_mom;   //!
   TBranch        *b_pim_tof_new;   //!
   TBranch        *b_pim_tofino_mult;   //!
   TBranch        *b_pim_track_length;   //!
   TBranch        *b_pim_z;   //!
   TBranch        *b_pip_beta;   //!
   TBranch        *b_pip_beta_new;   //!
   TBranch        *b_pip_dedx_in;   //!
   TBranch        *b_pip_dedx_in_sigma;   //!
   TBranch        *b_pip_dedx_mdc;   //!
   TBranch        *b_pip_dedx_mdc_sigma;   //!
   TBranch        *b_pip_dedx_out;   //!
   TBranch        *b_pip_dedx_out_sigma;   //!
   TBranch        *b_pip_dedx_tof;   //!
   TBranch        *b_pip_id;   //!
   TBranch        *b_pip_isring;   //!
   TBranch        *b_pip_kIsLepton;   //!
   TBranch        *b_pip_kIsUsed;   //!
   TBranch        *b_pip_mdcchi2;   //!
   TBranch        *b_pip_oa_hadr;   //!
   TBranch        *b_pip_oa_lept;   //!
   TBranch        *b_pip_p;   //!
   TBranch        *b_pip_phi;   //!
   TBranch        *b_pip_q;   //!
   TBranch        *b_pip_r;   //!
   TBranch        *b_pip_resolution;   //!
   TBranch        *b_pip_rkchi2;   //!
   TBranch        *b_pip_sector;   //!
   TBranch        *b_pip_shw_sum0;   //!
   TBranch        *b_pip_shw_sum1;   //!
   TBranch        *b_pip_shw_sum2;   //!
   TBranch        *b_pip_sim_corrflag;   //!
   TBranch        *b_pip_sim_geninfo;   //!
   TBranch        *b_pip_sim_geninfo1;   //!
   TBranch        *b_pip_sim_geninfo2;   //!
   TBranch        *b_pip_sim_genweight;   //!
   TBranch        *b_pip_sim_id;   //!
   TBranch        *b_pip_sim_iscommon;   //!
   TBranch        *b_pip_sim_mediumid;   //!
   TBranch        *b_pip_sim_p;   //!
   TBranch        *b_pip_sim_parentid;   //!
   TBranch        *b_pip_sim_primaryflag;   //!
   TBranch        *b_pip_sim_processid;   //!
   TBranch        *b_pip_sim_px;   //!
   TBranch        *b_pip_sim_py;   //!
   TBranch        *b_pip_sim_pz;   //!
   TBranch        *b_pip_sim_vertexx;   //!
   TBranch        *b_pip_sim_vertexy;   //!
   TBranch        *b_pip_sim_vertexz;   //!
   TBranch        *b_pip_system;   //!
   TBranch        *b_pip_theta;   //!
   TBranch        *b_pip_tof_exp;   //!
   TBranch        *b_pip_tof_mom;   //!
   TBranch        *b_pip_tof_new;   //!
   TBranch        *b_pip_tofino_mult;   //!
   TBranch        *b_pip_track_length;   //!
   TBranch        *b_pip_z;   //!
   TBranch        *b_runnumber;   //!
   TBranch        *b_totalmult;   //!
   TBranch        *b_trigbit;   //!
   TBranch        *b_trigdec;   //!
   TBranch        *b_trigdownscale;   //!
   TBranch        *b_trigdownscaleflag;   //!

   PPimPip(TTree *tree=0);
   virtual ~PPimPip();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

