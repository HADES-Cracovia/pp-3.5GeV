#include "EpEm.h"
#include "EpEp.h"
#include "EmEm.h"
#include "PEpEm.h"
#include "PEpEp.h"
#include "PEmEm.h"
#include "data.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <iostream>
#include "HFilter.h"
#include <fstream>
#include <TLine.h>
using namespace std;
using namespace PATData;

int main()
{
  //******************************************* B E A M   P A R A M E T E R S *****************************/
  /*
  //proj = new TLorentzVector(0,0,1700,1705.71972227); // PION BEAM 1.7 GeV/c momentum
  //targ = new TLorentzVector(0,0,0,938.27231); // PROTON
  beam = new TLorentzVector(0,0,0,0);
  proj = new TLorentzVector(0,0,700,713.7785); // PION BEAM 0.7 GeV/c momentum 
  targ = new TLorentzVector(0,0,0,938.27231); // PROTON
  *beam = *proj + *targ;
  */

  //PION BEAM SCANNER
  //   double pion_momentum = 612.0;

  //   double pion_momentum = 656.0;
  //   double pion_momentum = 690.0; 
  //   double pion_momentum = 748.0;
  //   double pion_momentum = 800.0;

  // ORIGIN: https://hades-wiki.gsi.de/foswiki/bin/view/PionBeam/WebHome

  double proton_energy = 3500+938.27231;//3.5 GeV
  double proton_momentum=TMath::Sqrt(proton_energy*proton_energy-938.27231*938.27231);

  proj = new TLorentzVector(0,0, proton_momentum, proton_energy); // PROTON BEAM momentum as above

  targ = new TLorentzVector(0,0,0,938.27231); // PROTON
  /*******************************************************************************************************/
  beam = new TLorentzVector(0,0,0,0);
  *beam = *proj + *targ;



  //*******************************************************************************************************/


  //filter = new HFilter;

  // *********************** FILES WITH CUTS ****************************
  insideEmS0 = -1;
  insideEmS1 = -1;
  insideEpS0 = -1;
  insideEpS1 = -1;

  insideEm1S0 = -1;
  insideEm1S1 = -1;
  insideEm2S0 = -1;
  insideEm2S1 = -1;

  insideEp1S0 = -1;
  insideEp1S1 = -1;
  insideEp2S0 = -1;
  insideEp2S1 = -1;

  e1 = new TLorentzVector(0,0,0,0);
  e2 = new TLorentzVector(0,0,0,0);
  p = new TLorentzVector(0,0,0,0);
  gammae1e2 = new TLorentzVector(0,0,0,0);
  e1e2 = new TLorentzVector(0,0,0,0);
  e1_delta = new TLorentzVector(0,0,0,0);
  e2_delta = new TLorentzVector(0,0,0,0);
  e1e2_miss = new TLorentzVector(0,0,0,0);
  pe1e2_miss = new TLorentzVector(0,0,0,0);
  /************************************** O U T P U T   F I L E ******************************************/
  outFileData = new TFile("output_bt.root","recreate");
  ofstream myfile;
  myfile.open ("raport.txt",ios::trunc);
  //outFileData = new TFile("ntuple_epem_656_C_gen1.root","recreate");
  //outFileData = new TFile("ntuple_epem_690_C_gen1.root","recreate");
  //outFileData = new TFile("ntuple_epem_690_PE_gen1_helicity.root","recreate");
  //outFileData = new TFile("ntuple_epem_690_C_gen1_helicity.root","recreate");
  //outFileData = new TFile("ntuple_epem_748_C_gen1.root","recreate");
  //outFileData = new TFile("ntuple_epem_800_C_gen1.root","recreate");
  /*******************************************************************************************************/

  /************************** control ntuple ***************************/
  tlo = new HNtuple("epem","epem");
  tlo->setFile( outFileData );
  /*********************************************************************/

  int signal_bins=24;
  double signal_min=0;
  double signal_max=1200;
  sig_all = new TH1F("sig_all","sig_all",signal_bins,signal_min,signal_max);
  sig_all->Sumw2();
  sig_all_bt = new TH1F("sig_all_bt","sig_all_bt",signal_bins,signal_min,signal_max);
  sig_all_bt->Sumw2();
  sig_all_back1 = new TH1F("sig_all_back1","sig_all_back1",signal_bins,signal_min,signal_max);
  sig_all_back1->Sumw2();
  sig_all_back2 = new TH1F("sig_all_back2","sig_all_back2",signal_bins,signal_min,signal_max);
  sig_all_back2->Sumw2();
  sig_all_bt_back1 = new TH1F("sig_all_bt_back1","sig_all_bt_back1",signal_bins,signal_min,signal_max);
  sig_all_bt_back1->Sumw2();
  sig_all_bt_back2 = new TH1F("sig_all_bt_back2","sig_all_bt_back2",signal_bins,signal_min,signal_max);
  sig_all_bt_back2->Sumw2();
  sig_rf_and_bt=new TH1F("sig_rf_and_bt","sig_rf_and_bt",signal_bins,signal_min,signal_max);
  sig_rf_and_bt_OK=new TH1F("sig_rf_and_bt_OK","sig_rf_and_bt_OK",signal_bins,signal_min,signal_max);
  sig_rf_and_bt_OK->Sumw2();
  sig_rf_and_bt_back1=new TH1F("sig_rf_and_bt_back1","sig_rf_and_bt_back1",signal_bins,signal_min,signal_max);
  sig_rf_and_bt_back1->Sumw2();
  sig_rf_and_bt_back2=new TH1F("sig_rf_and_bt_back2","sig_rf_and_bt_back2",signal_bins,signal_min,signal_max);
  sig_rf_and_bt_back2->Sumw2();
   
  sig_OK = 0;
  sig_bt_OK = 0;
   
  miss_all = new TH1F("miss_all","miss_all",80,0.6,1.4);
  miss_all->Sumw2();
  miss_all_back1 = new TH1F("miss_all_back1","miss_all_back1",80,0.6,1.4);
  miss_all_back1->Sumw2();
  miss_all_back2 = new TH1F("miss_all_back2","miss_all_back2",80,0.6,1.4);
  miss_all_back2->Sumw2();
  miss_OK = 0;

  //Float_t xbins[] = {0.00,0.04, 0.08, 0.12, 0.2,0.28,0.36,0.46, 0.56, 0.7, 0.9, 1.2};
  //Float_t xbins[] = {0.00,0.02,0.04,0.06,0.08,0.10,0.12,0.16,0.2,0.24,0.28,0.32,0.36,0.39,0.46,0.51,0.56,0.63,0.7, 0.9, 1.2};
  Float_t xbins[] = {0.00,0.02,0.04,0.06,0.08,0.10,0.12,0.16,0.2,0.24,0.28,0.32,0.36,0.40,0.44,0.48,0.52,0.56,0.60,0.64,0.68,0.70,0.75,0.8,0.85,0.9,0.95,1.0,1.1, 1.2};
  Int_t nbins = sizeof(xbins)/sizeof(Float_t);
  for(int nn=0;nn<nbins;nn++)//resize to MeV
    xbins[nn]=xbins[nn]*1000;

  sig_all_var = new TH1F("sig_all_var","sig_all_var",nbins-1,xbins);
  sig_all_var->Sumw2();
  sig_all_var_bt = new TH1F("sig_all_var_bt","sig_all_var_bt",nbins-1,xbins);
  sig_all_var_bt->Sumw2();
  sig_all_var_back1 = new TH1F("sig_all_var_back1","sig_all_var_back1",nbins-1,xbins);
  sig_all_var_back1->Sumw2();
  sig_all_var_back2 = new TH1F("sig_all_var_back2","sig_all_var_back2",nbins-1,xbins);
  sig_all_var_back2->Sumw2();
  sig_all_var_back = new TH1F("sig_all_var_back","sig_all_var_back",nbins-1,xbins);
  sig_all_var_back->Sumw2();
  sig_all_var_bt_back1 = new TH1F("sig_all_var_bt_back1","sig_all_var_bt_back1",nbins-1,xbins);
  sig_all_var_bt_back1->Sumw2();
  sig_all_var_bt_back2 = new TH1F("sig_all_var_bt_back2","sig_all_var_bt_back2",nbins-1,xbins);
  sig_all_var_bt_back2->Sumw2();
  sig_all_var_bt_back = new TH1F("sig_all_var_bt_back","sig_all_var_bt_back",nbins-1,xbins);
  sig_all_var_bt_back->Sumw2();
  sig_var_bt_OK = 0;
  sig_var_OK = 0;

  sig_all_var2 = new TH1F("sig_all_var_no","sig_all_var_no_norm",nbins-1,xbins);
  sig_all_var2->Sumw2();
  sig_all_var2_back1 = new TH1F("sig_all_var_back1_no","sig_all_var_back1_no_norm",nbins-1,xbins);
  sig_all_var2_back1->Sumw2();
  sig_all_var2_back2 = new TH1F("sig_all_var_back2_no","sig_all_var_back2_no_norm",nbins-1,xbins);
  sig_all_var2_back2->Sumw2();
  sig_var2_OK = 0;

   
  Float_t xbinsang[] =  {-1.0, -0.9, -0.8, -0.7, -0.6, -0.4, -0.2, 0, 0.2, 0.4, 0.6, 0.7, 0.8, 0.9, 1.0};
  Int_t nbinsang = sizeof(xbinsang)/sizeof(Float_t);

  cos_ep = new TH1F("cos_ep","cos_ep",10,-1,1); cos_ep->Sumw2();
  cos_em = new TH1F("cos_em","cos_em",10,-1,1); cos_em->Sumw2();
  cos_sum = new TH1F("cos_sum","cos_sum",10,-1,1); cos_sum->Sumw2();
  cos_ep_back1 = new TH1F("cos_ep_back1","cos_ep_back1",10,-1,1); cos_ep_back1->Sumw2();
  cos_em_back1 = new TH1F("cos_em_back1","cos_em_back1",10,-1,1); cos_em_back1->Sumw2();
  cos_sum_back1 = new TH1F("cos_sum_back1","cos_sum_back1",10,-1,1); cos_sum_back1->Sumw2();
  cos_ep_back2 = new TH1F("cos_ep_back2","cos_ep_back2",10,-1,1); cos_ep_back2->Sumw2();
  cos_em_back2 = new TH1F("cos_em_back2","cos_em_back2",10,-1,1); cos_em_back2->Sumw2();
  cos_sum_back2 = new TH1F("cos_sum_back2","cos_sum_back2",10,-1,1); cos_sum_back2->Sumw2();

  cos_ep_cm = new TH1F("cos_ep_cm","cos_ep_cm",10,-1,1); cos_ep_cm->Sumw2();
  cos_back1_cm = new TH1F("cos_back1_cm","cos_back1_cm",10,-1,1); cos_back1_cm->Sumw2();
  cos_back2_cm = new TH1F("cos_back2_cm","cos_back2_cm",10,-1,1); cos_back2_cm->Sumw2();
  cos_ep_cm_OK = new TH1F("cos_ep_cm_OK","cos_ep_cm_OK",10,-1,1); cos_ep_cm_OK->Sumw2();


  rapidity_all = new TH1F("rapidity_all","rapidity_all",50,0,2); rapidity_all->Sumw2();
  rapidity_back1 = new TH1F("rapidity_back1","rapidity_back1",50,0,2); rapidity_back1->Sumw2();
  rapidity_back2 = new TH1F("rapidity_back2","rapidity_back2",50,0,2); rapidity_back2->Sumw2();
  rapidity_OK = 0;

  rapidity_140_all = new TH1F("rapidity_140_all","rapidity_140_all",50,0,2); rapidity_140_all->Sumw2();
  rapidity_140_back1 = new TH1F("rapidity_140_back1","rapidity_140_back1",50,0,2); rapidity_140_back1->Sumw2();
  rapidity_140_back2 = new TH1F("rapidity_140_back2","rapidity_140_back2",50,0,2); rapidity_140_back2->Sumw2();
  rapidity_140_OK = 0;

  pt_all = new TH1F("pt_all","pt_all",50,0,1); pt_all->Sumw2();
  pt_back1 = new TH1F("pt_back1","pt_back1",50,0,1); pt_back1->Sumw2();
  pt_back2 = new TH1F("pt_back2","pt_back2",50,0,1); pt_back2->Sumw2();
  pt_OK = 0;

  pt_140_all = new TH1F("pt_140_all","pt_140_all",50,0,1); pt_140_all->Sumw2();
  pt_140_back1 = new TH1F("pt_140_back1","pt_140_back1",50,0,1); pt_140_back1->Sumw2();
  pt_140_back2 = new TH1F("pt_140_back2","pt_140_back2",50,0,1); pt_140_back2->Sumw2();
  pt_140_OK = 0;

  ep_beta_mom = new TH2F("ep_beta_mom","ep_beta_mom",50,0.6,1.3,100,0,1000); ep_beta_mom->Sumw2();
  em_beta_mom = new TH2F("em_beta_mom","em_beta_mom",50,0.6,1.3,100,0,1000); em_beta_mom->Sumw2();
  ep_beta_mom_bt = new TH2F("ep_beta_mom_bt","ep_beta_mom_bt",50,0.6,1.3,100,0,1000); ep_beta_mom_bt->Sumw2();
  em_beta_mom_bt = new TH2F("em_beta_mom_bt","em_beta_mom_bt",50,0.6,1.3,100,0,1000); em_beta_mom_bt->Sumw2();

  em_mom= new TH1F("em_mom","em_mom",250,0,1000);
  ep_mom= new TH1F("ep_mom","ep_mom",250,0,1000);
  em_mom_bt= new TH1F("em_mom_bt","em_mom_bt",250,0,1000);
  ep_mom_bt= new TH1F("ep_mom_bt","ep_mom_bt",250,0,1000);

  pureBT_signal=new TH1F("pureBT_signal","pureBT_signal",signal_bins,signal_min,signal_max);
  pureBT_signal_OK=new TH1F("pureBT_signal_OK","pureBT_signal_OK",signal_bins,signal_min,signal_max);
  pureBT_signal_back1=new TH1F("pureBT_signal_back1","pureBT_signal_back1",signal_bins,signal_min,signal_max);
  pureBT_signal_back2=new TH1F("pureBT_signal_back2","pureBT_signal_back2",signal_bins,signal_min,signal_max);

  pureBT_signal_var=new TH1F("pureBT_signal_var","pureBT_signal_var",nbins-1,xbins);
  pureBT_signal_OK_var=new TH1F("pureBT_signal_OK_var","pureBT_signal_OK_var",nbins-1,xbins);
  pureBT_signal_back1_var=new TH1F("pureBT_signal_back1_var","pureBT_signal_back1_var",nbins-1,xbins);
  pureBT_signal_back2_var=new TH1F("pureBT_signal_back2_var","pureBT_signal_back2_var",nbins-1,xbins);
  pureBT_signal_back_var=new TH1F("pureBT_signal_back_var","pureBT_signal_back_var",nbins-1,xbins);
   
  sig_sum=new TH1F("sig_sum","sig_sum",signal_bins,signal_min,signal_max);
  sig_sum_var=new TH1F("sig_sum_var","sig_sum_var",nbins-1,xbins);
   
  pureBT_beta_mom=new TH2F("pureBT_beta_mom","pureBT_beta_mom",50,0.7,1.3,100,0,1000);

  bt_rf_stat=new TH1I("bt_rf_stat","Statistic how many events heppend",23,-0.5,22.5);
  bt_rf_stat_pi=new TH1I("bt_rf_stat_pi","Statistic how many events heppend for 140 MeV<mass<700 MeV ",14,-0.5,13.5);
  bt_rf_stat_back1=new TH1I("bt_rf_stat_back1","Statistic how many events heppend",23,-0.5,22.5);
  bt_rf_stat_pi_back1=new TH1I("bt_rf_stat_pi_back1","Statistic how many events heppend for 140 MeV<mass<700 MeV ",14,-0.5,13.5);
  bt_rf_stat_back2=new TH1I("bt_rf_stat_back2","Statistic how many events heppend",23,-0.5,22.5);
  bt_rf_stat_pi_back2=new TH1I("bt_rf_stat_pi_back2","Statistic how many events heppend for 140 MeV<mass<700 MeV ",14,-0.5,13.5);
  bt_rf_stat_OK=new TH1I("bt_rf_stat_OK","Statistic how many events heppend",23,-0.5,22.5);
  bt_rf_stat_pi_OK=new TH1I("bt_rf_stat_pi_OK","Statistic how many events heppend for 140 MeV<mass<700 MeV ",14,-0.5,13.5);

  rf_freedom=new TH3F("rf_freedom","rf_freedom:d_theta:d_phi*sin(theta):mom",50,-2,2,50,-2,2,20,0,800);
  rf_f_dtheta=new TH2F("rf_f_dtheta","fr_f_dtheta",50,-4,4,20,0,800);
  rf_f_dphi=new TH2F("rf_f_dphi","fr_f_dphi",50,-4,4,20,0,800);

  momentum_spectrum=new TH1F("momentum_spectrum","leptons from RF momentum",400,-2000,2000);
  momentum_spectrum_bt=new TH1F("momentum_spectrum_bt","leptons from BT momentum",400,-2000,2000);
  momentum_spectrum_pureBT=new TH1F("momentum_spectrum_pureBT","leptons from BT profit",400,-2000,2000);

  sig_to_bg_var=new TH1F("sig_to_bg_var","signal to background ratio",nbins-1,xbins);
  sig_to_bg_bt_var=new TH1F("sig_to_bg_bt_var","signal to background ratio",nbins-1,xbins);
  sig_to_bg_pureBT_var=new TH1F("sig_to_bg_pureBT_var","signal to background ratio",nbins-1,xbins);

  q_vs_p_leptons_RF=new TH2F("q_vs_p_leptons_RF","charge in pre-shower for leptons from RF",600,0,1500,180,-100,200);
  q_vs_p_leptons_BT=new TH2F("q_vs_p_leptons_BT","charge in pre-shower for leptons from BT",600,0,1500,180,-100,200);

  int z_bins=30;
  int theta_bins=30;
  int zmin=-80;
  int zmax=0;
  int thetamin=0;
  int thetamax=90;
  z_theta_epep=new TH2F("z_theta_epep","Z vs. #theta for leptons from epep pairs",z_bins,zmin,zmax,theta_bins,thetamin,thetamax);
  z_theta_emem=new TH2F("z_theta_emem","Z vs. #theta for leptons from emem pairs",z_bins,zmin,zmax,theta_bins,thetamin,thetamax);
  z_theta_epem=new TH2F("z_theta_epem","Z vs. #theta for leptons from epem pairs",z_bins,zmin,zmax,theta_bins,thetamin,thetamax);
  z_theta_all=new TH2F("z_theta_all","Z vs. #theta for leptons from all pairs",z_bins,zmin,zmax,theta_bins,thetamin,thetamax);

  char hname[20];
  for(int h=0;h<9;h++)
    {
      sprintf(hname,"phi_theta_rich_%d",h);
      phi_theta_rich[h]=new TH2F(hname,"cut for RICH-MDC match",50,-3,3,50,-3,3);
    }

  proton_p_beta=new TH2F("proton_p_beta","Momentum versus beta for protons",800,-4000,4000,100,0,1.4);
  proton_p=new TH1F("proton_p","proton momentum",2000,0,4000);
  proton_E=new TH1F("proton_E","proton energy",720,900,4500);
  proton_m=new TH1F("proton_m","reconstructed proton mass",1000,500,1500);

  int miss_min=0;
  int miss_max=2000;
  int miss_n=100;
  miss_mass_all=new TH1F("miss_mass_all","Missing mass peak for all epem",miss_n,miss_min,miss_max);
  miss_mass_BT=new TH1F("miss_mass_bt","Missing mass peak for epem detected by BT",miss_n,miss_min,miss_max);
  miss_mass_RF=new TH1F("miss_mass_RF","Missing mass peak for epem detected by RF",miss_n,miss_min,miss_max);
  miss_mass_profit=new TH1F("miss_mass_profit","Missing mass peak for profit epem",miss_n,miss_min,miss_max);

  miss_mass_all_bcg1=new TH1F("miss_mass_all_bcg1","Missing mass peak for all epep",miss_n,miss_min,miss_max);
  miss_mass_BT_bcg1=new TH1F("miss_mass_bt_bcg1","Missing mass peak for epep detected by BT",miss_n,miss_min,miss_max);
  miss_mass_RF_bcg1=new TH1F("miss_mass_RF_bcg1","Missing mass peak for epep detected by RF",miss_n,miss_min,miss_max);
  miss_mass_profit_bcg1=new TH1F("miss_mass_profit_bcg1","Missing mass peak for profit epep",miss_n,miss_min,miss_max);

  miss_mass_all_bcg2=new TH1F("miss_mass_all_bcg2","Missing mass peak for all emem",miss_n,miss_min,miss_max);
  miss_mass_BT_bcg2=new TH1F("miss_mass_bt_bcg2","Missing mass peak for emem detected by BT",miss_n,miss_min,miss_max);
  miss_mass_RF_bcg2=new TH1F("miss_mass_RF_bcg2","Missing mass peak for emem  detected by RF",miss_n,miss_min,miss_max);
  miss_mass_profit_bcg2=new TH1F("miss_mass_profit_bcg2","Missing mass peak for profit emem",miss_n,miss_min,miss_max);

  miss_mass_all_bcg=new TH1F("miss_mass_all_bcg","Missing mass peak for all background pairs",miss_n,miss_min,miss_max);
  miss_mass_BT_bcg=new TH1F("miss_mass_bt_bcg","Missing mass peak for background pairs detected by BT",miss_n,miss_min,miss_max);
  miss_mass_RF_bcg=new TH1F("miss_mass_RF_bcg","Missing mass peak for background pairs  detected by RF",miss_n,miss_min,miss_max);
  miss_mass_profit_bcg=new TH1F("miss_mass_profit_bcg","Missing mass peak for profit background pairs",miss_n,miss_min,miss_max);

  miss_mass_all_OK=new TH1F("miss_mass_all_OK","Missing mass peak for all epem",miss_n,miss_min,miss_max);
  miss_mass_BT_OK=new TH1F("miss_mass_bt_OK","Missing mass peak for epem detected by BT",miss_n,miss_min,miss_max);
  miss_mass_RF_OK=new TH1F("miss_mass_RF_OK","Missing mass peak for epem detected by RF",miss_n,miss_min,miss_max);
  miss_mass_profit_OK=new TH1F("miss_mass_profit_OK","Missing mass peak for profit epem",miss_n,miss_min,miss_max);

  /**************************** M A I N   P A R T ****************************************/
      
  PEpEm pepem;
  cout << "START PEPEM!"<< endl;
  pepem.Loop();
   
  PEpEp pepem_back1;
  cout << "START PEPEP!"<< endl;
  pepem_back1.Loop();
   
  PEmEm pepem_back2;
  cout << "START PEMEM!"<< endl;
  pepem_back2.Loop();
  /*
    EpEm epem;
    cout << "START EPEM!" << endl;
    epem.Loop();
   
    EpEp epem_back1;
    cout << "START EPEP!" << endl;
    epem_back1.Loop();

    EmEm epem_back2;
    cout << "START EMEM!" << endl;
    epem_back2.Loop();   
  */
   
  cout << "FINALIZING!" << endl;


  /***************************** F I N A L     C A L C U L A T I O N S ********************/

  sig_OK = (TH1F*)signal("sig_OK", sig_all, sig_all_back1, sig_all_back2);
  sig_bt_OK = (TH1F*)signal("sig_bt_OK", sig_all_bt, sig_all_bt_back1, sig_all_bt_back2);
  pureBT_signal_OK=(TH1F*)signal("pureBT_signal_OK",pureBT_signal,pureBT_signal_back1,pureBT_signal_back2);
  sig_rf_and_bt_OK=(TH1F*)signal("sig_rf_and_bt_OK",sig_rf_and_bt,sig_rf_and_bt_back1,sig_rf_and_bt_back2);
  pureBT_signal_OK_var=(TH1F*)signal("pureBT_signal_OK_var",pureBT_signal_var,pureBT_signal_back1_var,pureBT_signal_back2_var);
  bt_rf_stat_OK=(TH1I*)signal("bt_rf_stat_OK",bt_rf_stat,bt_rf_stat_back1,bt_rf_stat_back2);
  bt_rf_stat_pi_OK=(TH1I*)signal("bt_rf_stat_pi_OK",bt_rf_stat_pi,bt_rf_stat_pi_back1,bt_rf_stat_pi_back2);
  miss_mass_all_OK=(TH1F*)signal("miss_mass_all_OK",miss_mass_all,miss_mass_all_bcg1,miss_mass_all_bcg2);
  miss_mass_BT_OK=(TH1F*)signal("miss_mass_BT_OK",miss_mass_BT,miss_mass_BT_bcg1,miss_mass_BT_bcg2);
  miss_mass_RF_OK=(TH1F*)signal("miss_mass_RF_OK",miss_mass_RF,miss_mass_RF_bcg1,miss_mass_RF_bcg2);
  miss_mass_profit_OK=(TH1F*)signal("miss_mass_profit_OK",miss_mass_profit,miss_mass_profit_bcg1,miss_mass_profit_bcg2);

  //normalize( sig_all );
  //normalize( sig_all_back1 );
  //normalize( sig_all_back2 );
  //normalize( sig_OK );
  miss_OK = (TH1F*)signal("miss_OK", miss_all, miss_all_back1, miss_all_back2);
  //normalize( miss_all );
  //normalize( miss_all_back1 );
  //normalize( miss_all_back2 );
  //normalize( miss_OK );
         
  sig_var_OK = (TH1F*)signal("sig_var_OK", sig_all_var, sig_all_var_back1, sig_all_var_back2);
  sig_var_bt_OK = (TH1F*)signal("sig_var_bt_OK", sig_all_var_bt, sig_all_var_bt_back1, sig_all_var_bt_back2);
  sig_var2_OK = (TH1F*)signal("sig_var2_OK", sig_all_var2, sig_all_var2_back1, sig_all_var2_back2);

   
  normalize( sig_all_var );
  normalize( sig_all_var_bt );
   
  normalize( sig_all_var_back1 );
  normalize( sig_all_var_back2 );
  normalize( sig_var_OK );

  normalize( sig_all_var_bt_back1 );
  normalize( sig_all_var_bt_back2 );
  normalize( sig_var_bt_OK );

  normalize( pureBT_signal_var );
  normalize( pureBT_signal_back1_var );
  normalize( pureBT_signal_back2_var );
  normalize( pureBT_signal_OK_var );
   

  sig_all_var_bt_back->Add(sig_all_var_bt_back1, sig_all_var_bt_back2);
  sig_all_var_back->Add(sig_all_var_back1, sig_all_var_back2);
  pureBT_signal_back_var->Add(pureBT_signal_back1_var,pureBT_signal_back2_var);

  sig_to_bg_var->Divide(sig_var_OK,sig_all_var_back);
  sig_to_bg_bt_var->Divide(sig_var_bt_OK,sig_all_var_bt_back);
  sig_to_bg_pureBT_var->Divide(pureBT_signal_OK_var,pureBT_signal_back_var);
   
  sig_sum->Add(sig_OK, pureBT_signal_OK,1,1);
  sig_sum_var->Add(sig_var_OK, pureBT_signal_OK_var,1,1);


  miss_mass_all_bcg->Add(miss_mass_all_bcg1,miss_mass_all_bcg2);
  miss_mass_BT_bcg->Add(miss_mass_BT_bcg1,miss_mass_BT_bcg2);
  miss_mass_RF_bcg->Add(miss_mass_RF_bcg1, miss_mass_RF_bcg2);
  miss_mass_profit_bcg->Add(miss_mass_profit_bcg1,miss_mass_profit_bcg2);
   
  SIGNAL = (TH1F*)sig_var2_OK->Clone("SIGNAL");
  SIGNAL->SetNameTitle("SIGNAL","SIGNAL");
  CB = (TH1F*)sig_all_var2_back1->Clone("CB");
  CB->Add(sig_all_var2_back2);
  CB->SetNameTitle("CB","CB");

  SIGNIFICANCE = (TH1F*)CB->Clone("SIGNIFICANCE");
  SIGNIFICANCE->Divide( SIGNAL );
  SIGNIFICANCE->Scale(2.);
  for (Int_t j=1; j<SIGNIFICANCE->GetNbinsX()+1; ++j)
    {
      SIGNIFICANCE->SetBinContent( j, SIGNIFICANCE->GetBinContent(j) + 1. );
    }
  TH1F *temp = (TH1F*)SIGNAL->Clone("temp");
  temp->Divide( SIGNIFICANCE );
  SIGNIFICANCE = temp;
  SIGNIFICANCE->SetNameTitle("SIGNIFICANCE", "SIGNIFICANCE");

  normalize( SIGNIFICANCE );
  normalize( SIGNAL );
  normalize( CB );

  cos_ep_OK = (TH1F*)signal("cos_ep_OK", cos_ep, cos_ep_back1, cos_ep_back2);
  //   normalize( cos_ep );
  //   normalize( cos_ep_back1 );
  //   normalize( cos_ep_back2 );
  //   normalize( cos_ep_OK );

  cos_em_OK = (TH1F*)signal("cos_em_OK", cos_em, cos_em_back1, cos_em_back2);
  //   normalize( cos_em );
  //   normalize( cos_em_back1 );
  //   normalize( cos_em_back2 );
  //   normalize( cos_em_OK );

  cos_sum_OK = (TH1F*)signal("cos_sum_OK", cos_sum, cos_sum_back1, cos_sum_back2);
  //   normalize( cos_sum );
  //   normalize( cos_sum_back1 );
  //   normalize( cos_sum_back2 );
  //   normalize( cos_sum_OK );

  cos_ep_cm_OK = (TH1F*)signal("cos_ep_cm_OK", cos_ep_cm, cos_back1_cm, cos_back2_cm);
  //   normalize( cos_ep_cm_OK );
  //   normalize( cos_back1_cm );
  //   normalize( cos_back2_cm );
  //   normalize( cos_ep_cm );
   
  rapidity_OK = (TH1F*)signal("rapidity_OK", rapidity_all, rapidity_back1, rapidity_back2);
  //   normalize( rapidity_all );
  //   normalize( rapidity_back1 );
  //   normalize( rapidity_back2 );
  //   normalize( rapidity_OK );

  rapidity_140_OK = (TH1F*)signal("rapidity_140_OK", rapidity_140_all, rapidity_140_back1, rapidity_140_back2);
  //   normalize( rapidity_140_all );
  //  normalize( rapidity_140_back1 );
  //   normalize( rapidity_140_back2 );
  //   normalize( rapidity_140_OK );

  pt_OK = (TH1F*)signal("pt_OK", pt_all, pt_back1, pt_back2);
  //   normalize( pt_all );
  //   normalize( pt_back1 );
  //   normalize( pt_back2 );
  //   normalize( pt_OK );

  pt_140_OK = (TH1F*)signal("pt_140_OK", pt_140_all, pt_140_back1, pt_140_back2);
  //   normalize( pt_140_all );
  //   normalize( pt_140_back1 );
  //   normalize( pt_140_back2 );
  //   normalize( pt_140_OK );

  /****************************************************************************************/

  outFileData->cd();

  // tlo->Write();
  
  TCanvas* cSpectra=new TCanvas("cSpectra","cSpectra");
  cSpectra->Divide(2);
  cSpectra->cd(1);
  gPad->SetLogy();
  sig_all_var->Draw();
  sig_all_var_bt->SetLineColor(kRed);
  sig_all_var_bt->Draw("same");
  cSpectra->cd(2);
  gPad->SetLogy();
  sig_var_OK->SetLineColor(kBlue);
  sig_var_OK->Draw();
  sig_var_bt_OK->SetLineColor(kRed);
  sig_var_bt_OK->Draw("SAME");

  TCanvas* cMomenta=new TCanvas("cMomenta","cMomenta");
  cMomenta->Divide(2);
  cMomenta->cd(1);
  ep_mom_bt->SetLineColor(kRed);
  ep_mom_bt->Draw();
  ep_mom->Draw("same");
  cMomenta->cd(2);
  em_mom_bt->SetLineColor(kRed);
  em_mom_bt->Draw();
  em_mom->Draw("same");
   
  TCanvas* cBetaMom=new TCanvas("cBetaMom","cBetaMom");
  cBetaMom->Divide(2,2);
  cBetaMom->cd(1);
  ep_beta_mom->Draw("colz");
  cBetaMom->cd(2);
  em_beta_mom->Draw("colz");
  cBetaMom->cd(3);
  ep_beta_mom_bt->Draw("colz");
  cBetaMom->cd(4);
  em_beta_mom_bt->Draw("colz");

  TCanvas* cPureBT=new TCanvas("cPureBT","cPureBT");
  cPureBT->Divide(2);
  cPureBT->cd(1);
  gPad->SetLogy();
  sig_var_OK->Draw();
  pureBT_signal_OK_var->SetLineColor(kRed);
  pureBT_signal_OK_var->Draw("same");
  cPureBT->cd(2);
  pureBT_beta_mom->Draw("colz");

  TCanvas* cSiggAll=new TCanvas("cSiggAll","cSiggAll");
  cSiggAll->SetLogy();
  sig_var_OK->Draw();
  format(sig_var_OK);
  pureBT_signal_OK_var->Draw("same");
  format(pureBT_signal_OK_var);
  sig_sum_var->SetLineColor(kGreen);
  format(sig_sum_var);
  sig_sum_var->Draw("same");

  TCanvas* cBackground=new TCanvas("cBackground","cBackground");
  cBackground->Divide(3,2);
  cBackground->cd(1);
  gPad->SetLogy();
  sig_var_OK->SetLineColor(kBlue+1);
  sig_var_OK->Draw();
  sig_all_var_back->Draw("same");
  format(sig_var_OK);
  format(sig_all_var_back);
  cBackground->cd(2);
  sig_var_bt_OK->SetLineColor(kRed+1);
  sig_var_bt_OK->Draw();
  sig_all_var_bt_back->Draw("same");
  sig_all_var_bt_back->SetLineColor(kRed);
  gPad->SetLogy();
  format(sig_var_bt_OK);
  format(sig_all_var_bt_back);
  cBackground->cd(3);
  pureBT_signal_OK_var->Draw();
  pureBT_signal_back_var->Draw("SAME");
  gPad->SetLogy();
  format(pureBT_signal_OK_var);
  format(pureBT_signal_back_var);
  cBackground->cd(4);
  sig_to_bg_var->Draw();
  gPad->SetLogy();
  cBackground->cd(5);
  sig_to_bg_bt_var->Draw();
  gPad->SetLogy();
  cBackground->cd(6);
  sig_to_bg_pureBT_var->Draw();
  gPad->SetLogy();

  TCanvas* cBackground_normal= new TCanvas("cBackground_normal","cBackground_normal");
  cBackground_normal->Divide(2);
  cBackground_normal->cd(1);
  gPad->SetLogy();
  sig_all->Draw();
  sig_OK->SetLineColor(kBlue+2);
  sig_OK->Draw("same");
  cBackground_normal->cd(2);
  sig_all_bt->Draw();
  sig_all_bt->SetLineColor(kRed);
  sig_bt_OK->SetLineColor(kRed+2);
  sig_bt_OK->Draw("same");
  gPad->SetLogy();

  TCanvas* cLeptonMom=new TCanvas("cLeptonMom","cLeptonMom");
  momentum_spectrum->Draw();
  //momentum_spectrum->SetFillColor(kBlue);
  momentum_spectrum_bt->Draw("same");
  momentum_spectrum_bt->SetLineColor(kRed);
  //momentum_spectrum_bt->SetFillColor(kRed);
  momentum_spectrum_pureBT->Draw("same");
  momentum_spectrum_pureBT->SetLineColor(kGreen);

  TCanvas* cp_q=new TCanvas("cp_q","cp_q");
  double linex[1000];
  double liney[1000];
  for(int v=0; v<1000;v++)
    {
      linex[v]=v;
      liney[v]=parametrization(v);
    }
  TGraph *l2=new TGraph(1000,linex,liney);
  //TLine *l1=new TLine(0,parametrization(0),1400,parametrization(1400));
  //l1->SetLineColor(kRed);
  //l1->SetLineWidth(3);
  cp_q->Divide(2);
  cp_q->cd(1);
  q_vs_p_leptons_RF->Draw("colz");
  //l1->Draw("same");
  l2->Draw("Lsame");
  gPad->SetLogz();
  cp_q->cd(2);
  q_vs_p_leptons_BT->Draw("colz");
  //l1->Draw("same");
  l2->Draw("Lsame");
  gPad->SetLogz();

  TCanvas* cStatistics=new TCanvas("cStatistics","cStatistics");
  cStatistics->Divide(2);
  cStatistics->cd(1);
  bt_rf_stat_OK->GetXaxis()->SetBinLabel(11,"BT");
  bt_rf_stat_OK->GetXaxis()->SetBinLabel(12,"RF");
  bt_rf_stat_OK->GetXaxis()->SetBinLabel(13,"RF && RF");
  bt_rf_stat_OK->GetXaxis()->SetBinLabel(14,"Profit");
  bt_rf_stat_OK->GetXaxis()->SetBinLabel(15,"BT <140 MeV");
  bt_rf_stat_OK->GetXaxis()->SetBinLabel(16,"RF <140 MeV");
  bt_rf_stat_OK->GetXaxis()->SetBinLabel(17,"profit <140 MeV");
  bt_rf_stat_OK->GetXaxis()->SetBinLabel(18,"BT >700 MeV");
  bt_rf_stat_OK->GetXaxis()->SetBinLabel(19,"RF >700 MeV");
  bt_rf_stat_OK->GetXaxis()->SetBinLabel(20,"profit >700 MeV");
  bt_rf_stat_OK->GetXaxis()->SetBinLabel(21,"BT <700 >140 MeV");
  bt_rf_stat_OK->GetXaxis()->SetBinLabel(22,"RF <700 >140 MeV");
  bt_rf_stat_OK->GetXaxis()->SetBinLabel(23,"profit <700 MeV >140 MeV");
  bt_rf_stat_OK->Draw();
   
  cStatistics->cd(2);
  bt_rf_stat_pi_OK->GetXaxis()->SetBinLabel(11,"BT");
  bt_rf_stat_pi_OK->GetXaxis()->SetBinLabel(12,"RF");
  bt_rf_stat_pi_OK->GetXaxis()->SetBinLabel(13,"RF && RF");
  bt_rf_stat_pi_OK->GetXaxis()->SetBinLabel(14,"Profit");
  bt_rf_stat_pi_OK->Draw();

  TCanvas* cThetaPhi=new TCanvas("cThetaPhi","cThetaPhi");
  cThetaPhi->Divide(3,3);
   
  for(int hh=0;hh<9;hh++)
    {
      cThetaPhi->cd(hh+1);
      phi_theta_rich[hh]->Draw("colz");
    }

  TCanvas* cSigAndBg=new TCanvas("cSigAndBg","cSigAndBg");
  cSigAndBg->Divide(2);
  cSigAndBg->cd(1);
  sig_all_var->Draw();
  sig_all_var_back->Draw("same");
  gPad->SetLogy();
  cSigAndBg->cd(2);
  sig_all_var_bt->Draw();
  sig_all_var_bt_back->Draw("Same");
  gPad->SetLogy();

  TCanvas* cZvertexTheta= new TCanvas("cZvertexTheta","cZvertexTheta");
  cZvertexTheta->Divide(2,2);
  cZvertexTheta->cd(1);
  z_theta_epem->Draw("COLZ");
  cZvertexTheta->cd(2);
  z_theta_epep->Draw("COLZ");
  cZvertexTheta->cd(3);
  z_theta_emem->Draw("COLZ");
  cZvertexTheta->cd(4);
  z_theta_all->Draw("COLZ");

  TCanvas* cProton=new TCanvas("cProton","proton properties");
  cProton->Divide(2,2);
  cProton->cd(1);
  proton_E->Draw();
  cProton->cd(2);
  proton_p->Draw();
  cProton->cd(3);
  proton_p_beta->Draw("COLZ");
  cProton->cd(4);
  proton_m->Draw();

  TCanvas* cMissMass=new TCanvas("cMissMass","missing mass peaks");
  cMissMass->Divide(2,2);
  cMissMass->cd(1);
  miss_mass_all->Draw();
  cMissMass->cd(2);
  miss_mass_RF->Draw();
  cMissMass->cd(3);
  miss_mass_BT->Draw();
  cMissMass->cd(4);
  miss_mass_profit->Draw();

  TCanvas* cMissMassBg=new TCanvas("cMissMassBg","missing mass for background events");
  cMissMassBg->Divide(2,2);

  cMissMassBg->cd(1);
  miss_mass_all_bcg1->Draw();
  miss_mass_all_bcg2->SetLineColor(kRed);
  miss_mass_all_bcg2->Draw("same");

  cMissMassBg->cd(2);
  miss_mass_RF_bcg1->Draw();
  miss_mass_RF_bcg2->SetLineColor(kRed);
  miss_mass_RF_bcg2->Draw("same");

  cMissMassBg->cd(3);
  miss_mass_BT_bcg1->Draw();
  miss_mass_BT_bcg2->SetLineColor(kRed);
  miss_mass_BT_bcg2->Draw("same");

  cMissMassBg->cd(4);
  miss_mass_profit_bcg1->Draw();
  miss_mass_profit_bcg2->SetLineColor(kRed);
  miss_mass_profit_bcg2->Draw("same");

  TCanvas* cMissMassOK=new TCanvas("cMissMassOK","missing mass signal");
  cMissMassOK->Divide(2,2);
  cMissMassOK->cd(1);
  miss_mass_all_OK->Draw();
  cMissMassOK->cd(2);
  miss_mass_BT_OK->Draw();
  cMissMassOK->cd(3);
  miss_mass_RF_OK->Draw();
  cMissMassOK->cd(4);
  miss_mass_profit_OK->Draw();
   
  cBackground_normal->Write();
  cBetaMom->Write();
  cSpectra->Write();
  cMomenta->Write();
  cPureBT->Write();
  cSiggAll->Write();
  cBackground->Write();
  cLeptonMom->Write();
  cp_q->Write();
  cStatistics->Write();
  cThetaPhi->Write();
  cSigAndBg->Write();
  cZvertexTheta->Write();
  cProton->Write();
  cMissMass->Write();
  cMissMassBg->Write();
  cMissMassOK->Write();
   
  ep_beta_mom->Write();
  em_beta_mom->Write();
  ep_beta_mom_bt->Write();
  em_beta_mom_bt->Write();

  sig_all->Write();
  sig_all_bt->Write();
  sig_all_back1->Write();
  sig_all_back2->Write();
  sig_OK->Write();
  sig_all_bt_back1->Write();
  sig_all_bt_back2->Write();
  sig_bt_OK->Write();
  sig_rf_and_bt->Write();
  sig_rf_and_bt_OK->Write();

  miss_all->Write();
  miss_all_back1->Write();
  miss_all_back2->Write();
  miss_OK->Write();

  sig_all_var->Write();
  sig_all_var_bt->Write();
  sig_all_var_back1->Write();
  sig_all_var_back2->Write();
  sig_all_var_back->Write();
  sig_var_OK->Write();
  sig_all_var_bt_back1->Write();
  sig_all_var_bt_back2->Write();
  sig_all_var_bt_back->Write();
  sig_var_bt_OK->Write();
   
  SIGNAL->Write();
  CB->Write();
  SIGNIFICANCE->Write();

  sig_all_var2->Write();
  sig_all_var2_back1->Write();
  sig_all_var2_back2->Write();
  sig_var2_OK->Write();

  cos_ep->Write();
  cos_ep_back1->Write();
  cos_ep_back2->Write();
  cos_ep_OK->Write();
  cos_em->Write();
  cos_em_back1->Write();
  cos_em_back2->Write();
  cos_em_OK->Write();
  cos_sum->Write();
  cos_sum_back1->Write();
  cos_sum_back2->Write();
  cos_sum_OK->Write();
 
  cos_ep_cm_OK->Write();
  cos_ep_cm->Write();
  cos_back1_cm->Write();
  cos_back2_cm->Write();

  rapidity_all->Write();
  rapidity_back1->Write();
  rapidity_back2->Write();
  rapidity_OK->Write();
   
  rapidity_140_all->Write();
  rapidity_140_back1->Write();
  rapidity_140_back2->Write();
  rapidity_140_OK->Write();
   
  pt_all->Write();
  pt_back1->Write();
  pt_back2->Write();
  pt_OK->Write();
   
  pt_140_all->Write();
  pt_140_back1->Write();
  pt_140_back2->Write();
  pt_140_OK->Write();

  ep_mom->Write();
  em_mom->Write();
  ep_mom_bt->Write();
  em_mom_bt->Write();

  pureBT_beta_mom->Write();
  pureBT_signal->Write();
  pureBT_signal_OK->Write();
  pureBT_signal_back1->Write();
  pureBT_signal_back2->Write();

  pureBT_beta_mom->Write();
  pureBT_signal_var->Write();
  pureBT_signal_OK_var->Write();
  pureBT_signal_back1_var->Write();
  pureBT_signal_back2_var->Write();
  pureBT_signal_back_var->Write();
  sig_sum->Write();
  sig_sum_var->Write();

  bt_rf_stat->Write();
  bt_rf_stat_pi->Write();
  bt_rf_stat_back1->Write();
  bt_rf_stat_back2->Write();
  bt_rf_stat_pi_back2->Write();
  bt_rf_stat_pi_back1->Write();
  bt_rf_stat_pi_OK->Write();
  bt_rf_stat_OK->Write();
   
  rf_freedom->Write();
  rf_f_dphi->Write();
  rf_f_dtheta->Write();

  momentum_spectrum->Write();
  momentum_spectrum_bt->Write();
  momentum_spectrum_pureBT->Write();

  sig_to_bg_var->Write();
  sig_to_bg_bt_var->Write();
  sig_to_bg_pureBT_var->Write();

  q_vs_p_leptons_BT->Write();
  q_vs_p_leptons_RF->Write();

  z_theta_epem->Write();
  z_theta_emem->Write();
  z_theta_epem->Write();

  miss_mass_profit->Write();
  miss_mass_BT->Write();
  miss_mass_RF->Write();
  miss_mass_all->Write();

  miss_mass_profit_bcg2->Write();
  miss_mass_BT_bcg2->Write();
  miss_mass_RF_bcg2->Write();
  miss_mass_all_bcg2->Write();

  miss_mass_profit_bcg1->Write();
  miss_mass_BT_bcg1->Write();
  miss_mass_RF_bcg1->Write();
  miss_mass_all_bcg1->Write();

  miss_mass_profit_OK->Write();
  miss_mass_BT_OK->Write();
  miss_mass_RF_OK->Write();
  miss_mass_all_OK->Write();
   
  //Save histopgrams parameters into text file
   
  myfile <<"Di-lepton parameters\n";
  myfile<< "Number of all di-leptons reconstructed in Pi0 mass region\n";
  myfile<<"  RF: "<<sig_OK->Integral(0,43)<<" BT: "<<sig_bt_OK->Integral(0,43)<<"\n";
  myfile<< "Number of di-leptons reconstructed E>mPi0\n";
  myfile<<"  RF: "<<sig_OK->Integral(44,300)<<" BT: "<<sig_bt_OK->Integral(44,300)<<"\n";
  myfile<<" Number of di-lepton recnstructed by BT and RF together in Pi0 mass region: \n";
  myfile<<"  = "<<sig_rf_and_bt_OK->Integral(0,43)<<endl;
  myfile<<" Number of di-lepton recnstructed by BT and RF together E>mPi0: \n";
  myfile<<"  = "<<sig_rf_and_bt_OK->Integral(44,300)<<endl;
  myfile<< "Number of di-leptons reconstructed by BT and missed by RF in Pi0 mass region:\n";
  myfile<<"  = "<<pureBT_signal_OK->Integral(0,43)<<"\n";
  myfile<< "Number of di-leptons reconstructed by BT and missed by RF E>mPi0:\n";
  myfile<<"  = "<<pureBT_signal_OK->Integral(44,300)<<"\n";
  myfile.close();

  outFileData->Close();
}

