#ifndef PPimPipPim_ID_buffer_h
#define PPipPimPim_ID_buffer_h

#include "PPimPipPim.h"

class PPimPipPim_ID_buffer
{
public :
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
  /*
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
  */  
  Float_t         p_system;
  Float_t         p_theta;
  Float_t         p_tof_exp;
  Float_t         p_tof_mom;
  Float_t         p_tof_new;
  Float_t         p_tofino_mult;
  Float_t         p_track_length;
  Float_t         p_z;
  Float_t         pim1_beta;
  Float_t         pim1_beta_new;
  Float_t         pim1_dedx_in;
  Float_t         pim1_dedx_in_sigma;
  Float_t         pim1_dedx_mdc;
  Float_t         pim1_dedx_mdc_sigma;
  Float_t         pim1_dedx_out;
  Float_t         pim1_dedx_out_sigma;
  Float_t         pim1_dedx_tof;
  Float_t         pim1_id;
  Float_t         pim1_isring;
  Float_t         pim1_kIsLepton;
  Float_t         pim1_kIsUsed;
  Float_t         pim1_mdcchi2;
  Float_t         pim1_oa_hadr;
  Float_t         pim1_oa_lept;
  Float_t         pim1_p;
  Float_t         pim1_phi;
  Float_t         pim1_q;
  Float_t         pim1_r;
  Float_t         pim1_resolution;
  Float_t         pim1_rkchi2;
  Float_t         pim1_sector;
  Float_t         pim1_shw_sum0;
  Float_t         pim1_shw_sum1;
  Float_t         pim1_shw_sum2;
  /*
    Float_t         pim1_sim_corrflag;
    Float_t         pim1_sim_geninfo;
    Float_t         pim1_sim_geninfo1;
    Float_t         pim1_sim_geninfo2;
    Float_t         pim1_sim_genweight;
    Float_t         pim1_sim_id;
    Float_t         pim1_sim_iscommon;
    Float_t         pim1_sim_mediumid;
    Float_t         pim1_sim_p;
    Float_t         pim1_sim_parentid;
    Float_t         pim1_sim_primaryflag;
    Float_t         pim1_sim_processid;
    Float_t         pim1_sim_px;
    Float_t         pim1_sim_py;
    Float_t         pim1_sim_pz;
    Float_t         pim1_sim_vertexx;
    Float_t         pim1_sim_vertexy;
    Float_t         pim1_sim_vertexz;
  */  
  Float_t         pim1_system;
  Float_t         pim1_theta;
  Float_t         pim1_tof_exp;
  Float_t         pim1_tof_mom;
  Float_t         pim1_tof_new;
  Float_t         pim1_tofino_mult;
  Float_t         pim1_track_length;
  Float_t         pim1_z;
  Float_t         pim2_beta;
  Float_t         pim2_beta_new;
  Float_t         pim2_dedx_in;
  Float_t         pim2_dedx_in_sigma;
  Float_t         pim2_dedx_mdc;
  Float_t         pim2_dedx_mdc_sigma;
  Float_t         pim2_dedx_out;
  Float_t         pim2_dedx_out_sigma;
  Float_t         pim2_dedx_tof;
  Float_t         pim2_id;
  Float_t         pim2_isring;
  Float_t         pim2_kIsLepton;
  Float_t         pim2_kIsUsed;
  Float_t         pim2_mdcchi2;
  Float_t         pim2_oa_hadr;
  Float_t         pim2_oa_lept;
  Float_t         pim2_p;
  Float_t         pim2_phi;
  Float_t         pim2_q;
  Float_t         pim2_r;
  Float_t         pim2_resolution;
  Float_t         pim2_rkchi2;
  Float_t         pim2_sector;
  Float_t         pim2_shw_sum0;
  Float_t         pim2_shw_sum1;
  Float_t         pim2_shw_sum2;
  /*
    Float_t         pim2_sim_corrflag;
    Float_t         pim2_sim_geninfo;
    Float_t         pim2_sim_geninfo1;
    Float_t         pim2_sim_geninfo2;
    Float_t         pim2_sim_genweight;
    Float_t         pim2_sim_id;
    Float_t         pim2_sim_iscommon;
    Float_t         pim2_sim_mediumid;
    Float_t         pim2_sim_p;
    Float_t         pim2_sim_parentid;
    Float_t         pim2_sim_primaryflag;
    Float_t         pim2_sim_processid;
    Float_t         pim2_sim_px;
    Float_t         pim2_sim_py;
    Float_t         pim2_sim_pz;
    Float_t         pim2_sim_vertexx;
    Float_t         pim2_sim_vertexy;
    Float_t         pim2_sim_vertexz;
  */  
  Float_t         pim2_system;
  Float_t         pim2_theta;
  Float_t         pim2_tof_exp;
  Float_t         pim2_tof_mom;
  Float_t         pim2_tof_new;
  Float_t         pim2_tofino_mult;
  Float_t         pim2_track_length;
  Float_t         pim2_z;
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
  /*
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
  */  
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

  PPimPipPim_ID_buffer(const PPimPipPim_ID_buffer& src)
  {
    eVert_chi2=src.eVert_chi2;
    eVert_x=src.eVert_x;
    eVert_y=src.eVert_y;
    eVert_z=src.eVert_z;
    event=src.event;
    hneg_mult=src.hneg_mult;
    hpos_mult=src.hpos_mult;
    isBest=src.isBest;
    lneg_mult=src.lneg_mult;
    lpos_mult=src.lpos_mult;
    p_beta=src.p_beta;
    p_beta_new=src.p_beta_new;
    p_dedx_in=src.p_dedx_in;
    p_dedx_in_sigma=src.p_dedx_in_sigma;
    p_dedx_mdc=src.p_dedx_mdc;
    p_dedx_mdc_sigma=src.p_dedx_mdc_sigma;
    p_dedx_out=src.p_dedx_out;
    p_dedx_out_sigma=src.p_dedx_out_sigma;
    p_dedx_tof=src.p_dedx_tof;
    p_id=src.p_id;
    p_isring=src.p_isring;
    p_kIsLepton=src.p_kIsLepton;
    p_kIsUsed=src.p_kIsUsed;
    p_mdcchi2=src.p_mdcchi2;
    p_oa_hadr=src.p_oa_hadr;
    p_oa_lept=src.p_oa_lept;
    p_p=src.p_p;
    p_phi=src.p_phi;
    p_q=src.p_q;
    p_r=src.p_r;
    p_resolution=src.p_resolution;
    p_rkchi2=src.p_rkchi2;
    p_sector=src.p_sector;
    p_shw_sum0=src.p_shw_sum0;
    p_shw_sum1=src.p_shw_sum1;
    p_shw_sum2=src.p_shw_sum2;
    /*
      p_sim_corrflag=src.p_sim_corrflag;
      p_sim_geninfo=src.p_sim_geninfo;
      p_sim_geninfo1=src.p_sim_geninfo1;
      p_sim_geninfo2=src.p_sim_geninfo2;
      p_sim_genweight=src.p_sim_genweight;
      p_sim_id=src.p_sim_id;
      p_sim_iscommon=src.p_sim_iscommon;
      p_sim_mediumid=src.p_sim_mediumid;
      p_sim_p=src.p_sim_p;
      p_sim_parentid=src.p_sim_parentid;
      p_sim_primaryflag=src.p_sim_primaryflag;
      p_sim_processid=src.p_sim_processid;
      p_sim_px=src.p_sim_px;
      p_sim_py=src.p_sim_py;
      p_sim_pz=src.p_sim_pz;
      p_sim_vertexx=src.p_sim_vertexx;
      p_sim_vertexy=src.p_sim_vertexy;
      p_sim_vertexz=src.p_sim_vertexz;
    */    
    p_system=src.p_system;
    p_theta=src.p_theta;
    p_tof_exp=src.p_tof_exp;
    p_tof_mom=src.p_tof_mom;
    p_tof_new=src.p_tof_new;
    p_tofino_mult=src.p_tofino_mult;
    p_track_length=src.p_track_length;
    p_z=src.p_z;
    pim1_beta=src.pim1_beta;
    pim1_beta_new=src.pim1_beta_new;
    pim1_dedx_in=src.pim1_dedx_in;
    pim1_dedx_in_sigma=src.pim1_dedx_in_sigma;
    pim1_dedx_mdc=src.pim1_dedx_mdc;
    pim1_dedx_mdc_sigma=src.pim1_dedx_mdc_sigma;
    pim1_dedx_out=src.pim1_dedx_out;
    pim1_dedx_out_sigma=src.pim1_dedx_out_sigma;
    pim1_dedx_tof=src.pim1_dedx_tof;
    pim1_id=src.pim1_id;
    pim1_isring=src.pim1_isring;
    pim1_kIsLepton=src.pim1_kIsLepton;
    pim1_kIsUsed=src.pim1_kIsUsed;
    pim1_mdcchi2=src.pim1_mdcchi2;
    pim1_oa_hadr=src.pim1_oa_hadr;
    pim1_oa_lept=src.pim1_oa_lept;
    pim1_p=src.pim1_p;
    pim1_phi=src.pim1_phi;
    pim1_q=src.pim1_q;
    pim1_r=src.pim1_r;
    pim1_resolution=src.pim1_resolution;
    pim1_rkchi2=src.pim1_rkchi2;
    pim1_sector=src.pim1_sector;
    pim1_shw_sum0=src.pim1_shw_sum0;
    pim1_shw_sum1=src.pim1_shw_sum1;
    pim1_shw_sum2=src.pim1_shw_sum2;
    /*
      pim1_sim_corrflag=src.pim1_sim_corrflag;
      pim1_sim_geninfo=src.pim1_sim_geninfo;
      pim1_sim_geninfo1=src.pim1_sim_geninfo1;
      pim1_sim_geninfo2=src.pim1_sim_geninfo2;
      pim1_sim_genweight=src.pim1_sim_genweight;
      pim1_sim_id=src.pim1_sim_id;
      pim1_sim_iscommon=src.pim1_sim_iscommon;
      pim1_sim_mediumid=src.pim1_sim_mediumid;
      pim1_sim_p=src.pim1_sim_p;
      pim1_sim_parentid=src.pim1_sim_parentid;
      pim1_sim_primaryflag=src.pim1_sim_primaryflag;
      pim1_sim_processid=src.pim1_sim_processid;
      pim1_sim_px=src.pim1_sim_px;
      pim1_sim_py=src.pim1_sim_py;
      pim1_sim_pz=src.pim1_sim_pz;
      pim1_sim_vertexx=src.pim1_sim_vertexx;
      pim1_sim_vertexy=src.pim1_sim_vertexy;
      pim1_sim_vertexz=src.pim1_sim_vertexz;
    */    
    pim1_system=src.pim1_system;
    pim1_theta=src.pim1_theta;
    pim1_tof_exp=src.pim1_tof_exp;
    pim1_tof_mom=src.pim1_tof_mom;
    pim1_tof_new=src.pim1_tof_new;
    pim1_tofino_mult=src.pim1_tofino_mult;
    pim1_track_length=src.pim1_track_length;
    pim1_z=src.pim1_z;
    pim2_beta=src.pim2_beta;
    pim2_beta_new=src.pim2_beta_new;
    pim2_dedx_in=src.pim2_dedx_in;
    pim2_dedx_in_sigma=src.pim2_dedx_in_sigma;
    pim2_dedx_mdc=src.pim2_dedx_mdc;
    pim2_dedx_mdc_sigma=src.pim2_dedx_mdc_sigma;
    pim2_dedx_out=src.pim2_dedx_out;
    pim2_dedx_out_sigma=src.pim2_dedx_out_sigma;
    pim2_dedx_tof=src.pim2_dedx_tof;
    pim2_id=src.pim2_id;
    pim2_isring=src.pim2_isring;
    pim2_kIsLepton=src.pim2_kIsLepton;
    pim2_kIsUsed=src.pim2_kIsUsed;
    pim2_mdcchi2=src.pim2_mdcchi2;
    pim2_oa_hadr=src.pim2_oa_hadr;
    pim2_oa_lept=src.pim2_oa_lept;
    pim2_p=src.pim2_p;
    pim2_phi=src.pim2_phi;
    pim2_q=src.pim2_q;
    pim2_r=src.pim2_r;
    pim2_resolution=src.pim2_resolution;
    pim2_rkchi2=src.pim2_rkchi2;
    pim2_sector=src.pim2_sector;
    pim2_shw_sum0=src.pim2_shw_sum0;
    pim2_shw_sum1=src.pim2_shw_sum1;
    pim2_shw_sum2=src.pim2_shw_sum2;
    /*
      pim2_sim_corrflag=src.pim2_sim_corrflag;
      pim2_sim_geninfo=src.pim2_sim_geninfo;
      pim2_sim_geninfo1=src.pim2_sim_geninfo1;
      pim2_sim_geninfo2=src.pim2_sim_geninfo2;
      pim2_sim_genweight=src.pim2_sim_genweight;
      pim2_sim_id=src.pim2_sim_id;
      pim2_sim_iscommon=src.pim2_sim_iscommon;
      pim2_sim_mediumid=src.pim2_sim_mediumid;
      pim2_sim_p=src.pim2_sim_p;
      pim2_sim_parentid=src.pim2_sim_parentid;
      pim2_sim_primaryflag=src.pim2_sim_primaryflag;
      pim2_sim_processid=src.pim2_sim_processid;
      pim2_sim_px=src.pim2_sim_px;
      pim2_sim_py=src.pim2_sim_py;
      pim2_sim_pz=src.pim2_sim_pz;
      pim2_sim_vertexx=src.pim2_sim_vertexx;
      pim2_sim_vertexy=src.pim2_sim_vertexy;
      pim2_sim_vertexz=src.pim2_sim_vertexz;
    */    
    pim2_system=src.pim2_system;
    pim2_theta=src.pim2_theta;
    pim2_tof_exp=src.pim2_tof_exp;
    pim2_tof_mom=src.pim2_tof_mom;
    pim2_tof_new=src.pim2_tof_new;
    pim2_tofino_mult=src.pim2_tofino_mult;
    pim2_track_length=src.pim2_track_length;
    pim2_z=src.pim2_z;
    pip_beta=src.pip_beta;
    pip_beta_new=src.pip_beta_new;
    pip_dedx_in=src.pip_dedx_in;
    pip_dedx_in_sigma=src.pip_dedx_in_sigma;
    pip_dedx_mdc=src.pip_dedx_mdc;
    pip_dedx_mdc_sigma=src.pip_dedx_mdc_sigma;
    pip_dedx_out=src.pip_dedx_out;
    pip_dedx_out_sigma=src.pip_dedx_out_sigma;
    pip_dedx_tof=src.pip_dedx_tof;
    pip_id=src.pip_id;
    pip_isring=src.pip_isring;
    pip_kIsLepton=src.pip_kIsLepton;
    pip_kIsUsed=src.pip_kIsUsed;
    pip_mdcchi2=src.pip_mdcchi2;
    pip_oa_hadr=src.pip_oa_hadr;
    pip_oa_lept=src.pip_oa_lept;
    pip_p=src.pip_p;
    pip_phi=src.pip_phi;
    pip_q=src.pip_q;
    pip_r=src.pip_r;
    pip_resolution=src.pip_resolution;
    pip_rkchi2=src.pip_rkchi2;
    pip_sector=src.pip_sector;
    pip_shw_sum0=src.pip_shw_sum0;
    pip_shw_sum1=src.pip_shw_sum1;
    pip_shw_sum2=src.pip_shw_sum2;
    /*
      pip_sim_corrflag=src.pip_sim_corrflag;
      pip_sim_geninfo=src.pip_sim_geninfo;
      pip_sim_geninfo1=src.pip_sim_geninfo1;
      pip_sim_geninfo2=src.pip_sim_geninfo2;
      pip_sim_genweight=src.pip_sim_genweight;
      pip_sim_id=src.pip_sim_id;
      pip_sim_iscommon=src.pip_sim_iscommon;
      pip_sim_mediumid=src.pip_sim_mediumid;
      pip_sim_p=src.pip_sim_p;
      pip_sim_parentid=src.pip_sim_parentid;
      pip_sim_primaryflag=src.pip_sim_primaryflag;
      pip_sim_processid=src.pip_sim_processid;
      pip_sim_px=src.pip_sim_px;
      pip_sim_py=src.pip_sim_py;
      pip_sim_pz=src.pip_sim_pz;
      pip_sim_vertexx=src.pip_sim_vertexx;
      pip_sim_vertexy=src.pip_sim_vertexy;
      pip_sim_vertexz=src.pip_sim_vertexz;
    */    
    pip_system=src.pip_system;
    pip_theta=src.pip_theta;
    pip_tof_exp=src.pip_tof_exp;
    pip_tof_mom=src.pip_tof_mom;
    pip_tof_new=src.pip_tof_new;
    pip_tofino_mult=src.pip_tofino_mult;
    pip_track_length=src.pip_track_length;
    pip_z=src.pip_z;
    runnumber=src.runnumber;
    totalmult=src.totalmult;
    trigbit=src.trigbit;
    trigdec=src.trigdec;
    trigdownscale=src.trigdownscale;
    trigdownscaleflag=src.trigdownscaleflag;
  }

  PPimPipPim_ID_buffer(PPimPipPim *src)
  {
    eVert_chi2=src->eVert_chi2;
    eVert_x=src->eVert_x;
    eVert_y=src->eVert_y;
    eVert_z=src->eVert_z;
    event=src->event;
    hneg_mult=src->hneg_mult;
    hpos_mult=src->hpos_mult;
    isBest=src->isBest;
    lneg_mult=src->lneg_mult;
    lpos_mult=src->lpos_mult;
    p_beta=src->p_beta;
    p_beta_new=src->p_beta_new;
    p_dedx_in=src->p_dedx_in;
    p_dedx_in_sigma=src->p_dedx_in_sigma;
    p_dedx_mdc=src->p_dedx_mdc;
    p_dedx_mdc_sigma=src->p_dedx_mdc_sigma;
    p_dedx_out=src->p_dedx_out;
    p_dedx_out_sigma=src->p_dedx_out_sigma;
    p_dedx_tof=src->p_dedx_tof;
    p_id=src->p_id;
    p_isring=src->p_isring;
    p_kIsLepton=src->p_kIsLepton;
    p_kIsUsed=src->p_kIsUsed;
    p_mdcchi2=src->p_mdcchi2;
    p_oa_hadr=src->p_oa_hadr;
    p_oa_lept=src->p_oa_lept;
    p_p=src->p_p;
    p_phi=src->p_phi;
    p_q=src->p_q;
    p_r=src->p_r;
    p_resolution=src->p_resolution;
    p_rkchi2=src->p_rkchi2;
    p_sector=src->p_sector;
    p_shw_sum0=src->p_shw_sum0;
    p_shw_sum1=src->p_shw_sum1;
    p_shw_sum2=src->p_shw_sum2;
    /*
      p_sim_corrflag=src->p_sim_corrflag;
      p_sim_geninfo=src->p_sim_geninfo;
      p_sim_geninfo1=src->p_sim_geninfo1;
      p_sim_geninfo2=src->p_sim_geninfo2;
      p_sim_genweight=src->p_sim_genweight;
      p_sim_id=src->p_sim_id;
      p_sim_iscommon=src->p_sim_iscommon;
      p_sim_mediumid=src->p_sim_mediumid;
      p_sim_p=src->p_sim_p;
      p_sim_parentid=src->p_sim_parentid;
      p_sim_primaryflag=src->p_sim_primaryflag;
      p_sim_processid=src->p_sim_processid;
      p_sim_px=src->p_sim_px;
      p_sim_py=src->p_sim_py;
      p_sim_pz=src->p_sim_pz;
      p_sim_vertexx=src->p_sim_vertexx;
      p_sim_vertexy=src->p_sim_vertexy;
      p_sim_vertexz=src->p_sim_vertexz;
    */    
    p_system=src->p_system;
    p_theta=src->p_theta;
    p_tof_exp=src->p_tof_exp;
    p_tof_mom=src->p_tof_mom;
    p_tof_new=src->p_tof_new;
    p_tofino_mult=src->p_tofino_mult;
    p_track_length=src->p_track_length;
    p_z=src->p_z;
    pim1_beta=src->pim1_beta;
    pim1_beta_new=src->pim1_beta_new;
    pim1_dedx_in=src->pim1_dedx_in;
    pim1_dedx_in_sigma=src->pim1_dedx_in_sigma;
    pim1_dedx_mdc=src->pim1_dedx_mdc;
    pim1_dedx_mdc_sigma=src->pim1_dedx_mdc_sigma;
    pim1_dedx_out=src->pim1_dedx_out;
    pim1_dedx_out_sigma=src->pim1_dedx_out_sigma;
    pim1_dedx_tof=src->pim1_dedx_tof;
    pim1_id=src->pim1_id;
    pim1_isring=src->pim1_isring;
    pim1_kIsLepton=src->pim1_kIsLepton;
    pim1_kIsUsed=src->pim1_kIsUsed;
    pim1_mdcchi2=src->pim1_mdcchi2;
    pim1_oa_hadr=src->pim1_oa_hadr;
    pim1_oa_lept=src->pim1_oa_lept;
    pim1_p=src->pim1_p;
    pim1_phi=src->pim1_phi;
    pim1_q=src->pim1_q;
    pim1_r=src->pim1_r;
    pim1_resolution=src->pim1_resolution;
    pim1_rkchi2=src->pim1_rkchi2;
    pim1_sector=src->pim1_sector;
    pim1_shw_sum0=src->pim1_shw_sum0;
    pim1_shw_sum1=src->pim1_shw_sum1;
    pim1_shw_sum2=src->pim1_shw_sum2;
    /*
      pim1_sim_corrflag=src->pim1_sim_corrflag;
      pim1_sim_geninfo=src->pim1_sim_geninfo;
      pim1_sim_geninfo1=src->pim1_sim_geninfo1;
      pim1_sim_geninfo2=src->pim1_sim_geninfo2;
      pim1_sim_genweight=src->pim1_sim_genweight;
      pim1_sim_id=src->pim1_sim_id;
      pim1_sim_iscommon=src->pim1_sim_iscommon;
      pim1_sim_mediumid=src->pim1_sim_mediumid;
      pim1_sim_p=src->pim1_sim_p;
      pim1_sim_parentid=src->pim1_sim_parentid;
      pim1_sim_primaryflag=src->pim1_sim_primaryflag;
      pim1_sim_processid=src->pim1_sim_processid;
      pim1_sim_px=src->pim1_sim_px;
      pim1_sim_py=src->pim1_sim_py;
      pim1_sim_pz=src->pim1_sim_pz;
      pim1_sim_vertexx=src->pim1_sim_vertexx;
      pim1_sim_vertexy=src->pim1_sim_vertexy;
      pim1_sim_vertexz=src->pim1_sim_vertexz;
    */    
    pim1_system=src->pim1_system;
    pim1_theta=src->pim1_theta;
    pim1_tof_exp=src->pim1_tof_exp;
    pim1_tof_mom=src->pim1_tof_mom;
    pim1_tof_new=src->pim1_tof_new;
    pim1_tofino_mult=src->pim1_tofino_mult;
    pim1_track_length=src->pim1_track_length;
    pim1_z=src->pim1_z;
    pim2_beta=src->pim2_beta;
    pim2_beta_new=src->pim2_beta_new;
    pim2_dedx_in=src->pim2_dedx_in;
    pim2_dedx_in_sigma=src->pim2_dedx_in_sigma;
    pim2_dedx_mdc=src->pim2_dedx_mdc;
    pim2_dedx_mdc_sigma=src->pim2_dedx_mdc_sigma;
    pim2_dedx_out=src->pim2_dedx_out;
    pim2_dedx_out_sigma=src->pim2_dedx_out_sigma;
    pim2_dedx_tof=src->pim2_dedx_tof;
    pim2_id=src->pim2_id;
    pim2_isring=src->pim2_isring;
    pim2_kIsLepton=src->pim2_kIsLepton;
    pim2_kIsUsed=src->pim2_kIsUsed;
    pim2_mdcchi2=src->pim2_mdcchi2;
    pim2_oa_hadr=src->pim2_oa_hadr;
    pim2_oa_lept=src->pim2_oa_lept;
    pim2_p=src->pim2_p;
    pim2_phi=src->pim2_phi;
    pim2_q=src->pim2_q;
    pim2_r=src->pim2_r;
    pim2_resolution=src->pim2_resolution;
    pim2_rkchi2=src->pim2_rkchi2;
    pim2_sector=src->pim2_sector;
    pim2_shw_sum0=src->pim2_shw_sum0;
    pim2_shw_sum1=src->pim2_shw_sum1;
    pim2_shw_sum2=src->pim2_shw_sum2;
    /*
      pim2_sim_corrflag=src->pim2_sim_corrflag;
      pim2_sim_geninfo=src->pim2_sim_geninfo;
      pim2_sim_geninfo1=src->pim2_sim_geninfo1;
      pim2_sim_geninfo2=src->pim2_sim_geninfo2;
      pim2_sim_genweight=src->pim2_sim_genweight;
      pim2_sim_id=src->pim2_sim_id;
      pim2_sim_iscommon=src->pim2_sim_iscommon;
      pim2_sim_mediumid=src->pim2_sim_mediumid;
      pim2_sim_p=src->pim2_sim_p;
      pim2_sim_parentid=src->pim2_sim_parentid;
      pim2_sim_primaryflag=src->pim2_sim_primaryflag;
      pim2_sim_processid=src->pim2_sim_processid;
      pim2_sim_px=src->pim2_sim_px;
      pim2_sim_py=src->pim2_sim_py;
      pim2_sim_pz=src->pim2_sim_pz;
      pim2_sim_vertexx=src->pim2_sim_vertexx;
      pim2_sim_vertexy=src->pim2_sim_vertexy;
      pim2_sim_vertexz=src->pim2_sim_vertexz;
    */    
    pim2_system=src->pim2_system;
    pim2_theta=src->pim2_theta;
    pim2_tof_exp=src->pim2_tof_exp;
    pim2_tof_mom=src->pim2_tof_mom;
    pim2_tof_new=src->pim2_tof_new;
    pim2_tofino_mult=src->pim2_tofino_mult;
    pim2_track_length=src->pim2_track_length;
    pim2_z=src->pim2_z;
    pip_beta=src->pip_beta;
    pip_beta_new=src->pip_beta_new;
    pip_dedx_in=src->pip_dedx_in;
    pip_dedx_in_sigma=src->pip_dedx_in_sigma;
    pip_dedx_mdc=src->pip_dedx_mdc;
    pip_dedx_mdc_sigma=src->pip_dedx_mdc_sigma;
    pip_dedx_out=src->pip_dedx_out;
    pip_dedx_out_sigma=src->pip_dedx_out_sigma;
    pip_dedx_tof=src->pip_dedx_tof;
    pip_id=src->pip_id;
    pip_isring=src->pip_isring;
    pip_kIsLepton=src->pip_kIsLepton;
    pip_kIsUsed=src->pip_kIsUsed;
    pip_mdcchi2=src->pip_mdcchi2;
    pip_oa_hadr=src->pip_oa_hadr;
    pip_oa_lept=src->pip_oa_lept;
    pip_p=src->pip_p;
    pip_phi=src->pip_phi;
    pip_q=src->pip_q;
    pip_r=src->pip_r;
    pip_resolution=src->pip_resolution;
    pip_rkchi2=src->pip_rkchi2;
    pip_sector=src->pip_sector;
    pip_shw_sum0=src->pip_shw_sum0;
    pip_shw_sum1=src->pip_shw_sum1;
    pip_shw_sum2=src->pip_shw_sum2;
    /*
      pip_sim_corrflag=src->pip_sim_corrflag;
      pip_sim_geninfo=src->pip_sim_geninfo;
      pip_sim_geninfo1=src->pip_sim_geninfo1;
      pip_sim_geninfo2=src->pip_sim_geninfo2;
      pip_sim_genweight=src->pip_sim_genweight;
      pip_sim_id=src->pip_sim_id;
      pip_sim_iscommon=src->pip_sim_iscommon;
      pip_sim_mediumid=src->pip_sim_mediumid;
      pip_sim_p=src->pip_sim_p;
      pip_sim_parentid=src->pip_sim_parentid;
      pip_sim_primaryflag=src->pip_sim_primaryflag;
      pip_sim_processid=src->pip_sim_processid;
      pip_sim_px=src->pip_sim_px;
      pip_sim_py=src->pip_sim_py;
      pip_sim_pz=src->pip_sim_pz;
      pip_sim_vertexx=src->pip_sim_vertexx;
      pip_sim_vertexy=src->pip_sim_vertexy;
      pip_sim_vertexz=src->pip_sim_vertexz;
    */    
    pip_system=src->pip_system;
    pip_theta=src->pip_theta;
    pip_tof_exp=src->pip_tof_exp;
    pip_tof_mom=src->pip_tof_mom;
    pip_tof_new=src->pip_tof_new;
    pip_tofino_mult=src->pip_tofino_mult;
    pip_track_length=src->pip_track_length;
    pip_z=src->pip_z;
    runnumber=src->runnumber;
    totalmult=src->totalmult;
    trigbit=src->trigbit;
    trigdec=src->trigdec;
    trigdownscale=src->trigdownscale;
    trigdownscaleflag=src->trigdownscaleflag;

  }
};
#endif
