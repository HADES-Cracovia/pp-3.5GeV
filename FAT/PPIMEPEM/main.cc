//#include "PPim.h"
#include "PPimEpEm.h"
#include "PPimEpEp.h"
#include "PPimEmEm.h"
//#include "EpEp.h"
//#include "EmEm.h"
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
  //double pion_momentum = 690.0; 
  //   double pion_momentum = 748.0;
  //   double pion_momentum = 800.0;
  double proton_momentum=TMath::Sqrt((3.5+0.938)*(3.5+0.938)-0.938*0.938)*1000;
  // ORIGIN: https://hades-wiki.gsi.de/foswiki/bin/view/PionBeam/WebHome
  //double pion_energy = sqrt( pion_momentum*pion_momentum + 139.56995*139.56995 );
  double proton_energy=sqrt(proton_momentum*proton_momentum + 938.27231*938.27231);
  
  proj = new TLorentzVector(0,0, proton_momentum, proton_energy); // PION BEAM momentum as above
  targ = new TLorentzVector(0,0,0,938.27231); // PROTON

  /*******************************************************************************************************/
  beam = new TLorentzVector(0,0,0,0);
  miss = new TLorentzVector(0,0,0,0);
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

  p = new TLorentzVector(0,0,0,0);
  pim = new TLorentzVector(0,0,0,0);
  em = new TLorentzVector(0,0,0,0);
  ep = new TLorentzVector(0,0,0,0);
  em1 = new TLorentzVector(0,0,0,0);
  ep1 = new TLorentzVector(0,0,0,0);
  em2 = new TLorentzVector(0,0,0,0);
  ep2 = new TLorentzVector(0,0,0,0);
  pi = new TLorentzVector(0,0,0,0);
  gammapep = new TLorentzVector(0,0,0,0);
  gammapem = new TLorentzVector(0,0,0,0);
  gammapep1 = new TLorentzVector(0,0,0,0);
  gammapem1 = new TLorentzVector(0,0,0,0);
  gammapep2 = new TLorentzVector(0,0,0,0);
  gammapem2 = new TLorentzVector(0,0,0,0);
  gammappim = new TLorentzVector(0,0,0,0);
  gammapimep = new TLorentzVector(0,0,0,0);
  gammapimep1 = new TLorentzVector(0,0,0,0);
  gammapimem1 = new TLorentzVector(0,0,0,0);
  gammaemep = new TLorentzVector(0,0,0,0);
  gammaem2em1 = new TLorentzVector(0,0,0,0);
  gammaep2ep1 = new TLorentzVector(0,0,0,0);
  gammappimepem = new TLorentzVector(0,0,0,0);
  gammappimep1ep2 = new TLorentzVector(0,0,0,0);
  gammappimem1em2 = new TLorentzVector(0,0,0,0);
  gammappimem=new TLorentzVector(0,0,0,0);
  gammappimep=new TLorentzVector(0,0,0,0);
  gammapemep=new TLorentzVector(0,0,0,0);
  
  ppi = new TLorentzVector(0,0,0,0);
  p_delta = new TLorentzVector(0,0,0,0);
  pi_delta = new TLorentzVector(0,0,0,0);
  ppi_miss = new TLorentzVector(0,0,0,0);

  /************************************** O U T P U T   F I L E ******************************************/
  //outFileData = new TFile("pp_ppimeppim_full_stat_dedx_extended_2_new_miss_mass.root","recreate");
  outFileData = new TFile("pp_fullstat_ppimepep.root","recreate");
  //ofstream myfile;
  //myfile.open ("raport.txt",ios::trunc);
  //outFileData = new TFile("ntuple_epem_656_C_gen1.root","recreate");
  //outFileData = new TFile("ntuple_epem_690_C_gen1.root","recreate");
  //outFileData = new TFile("ntuple_epem_690_PE_gen1_helicity.root","recreate");
  //outFileData = new TFile("ntuple_epem_690_C_gen1_helicity.root","recreate");
  //outFileData = new TFile("ntuple_epem_748_C_gen1.root","recreate");
  //outFileData = new TFile("ntuple_epem_800_C_gen1.root","recreate");
  /*******************************************************************************************************/

  /************************** control ntuple ***************************/
  tlo = new HNtuple("ppimepem","ppimepem");
  tlo_epep = new HNtuple("ppimepep","ppimepep");
  tlo_emem = new HNtuple("ppimemem","ppimemem");
  tlo->setFile( outFileData );
  tlo_epep->setFile( outFileData );
  tlo_emem->setFile( outFileData );
  /*********************************************************************/

  int beta_min=0;
  int beta_max=1.2;
  int beta_n=240;
  int p_min=0;
  int p_max=3000;
  int p_n=3000;
  /*  
  p_p_beta=new TH2F("p_p_beta","Momentum vs. beta for protons",p_n,p_min,p_max,beta_n,beta_min,beta_max);
  pim_p_beta=new TH2F("pim_p_beta","Momentum vs. beta for #pi^{-}",p_n,p_min,p_max,beta_n,beta_min,beta_max);
  ep_p_beta=new TH2F("ep_p_beta","Momentum vs. beta for #pi^{+}",p_n,p_min,p_max,beta_n,beta_min,beta_max);
  p_pim_mass=new TH1F("p_pim_mass","Invariant mass #pi^{-} p",2000,500,2500);									       
  p_pim_mass=new TH1F("p_pim_mass","Invariant mass #pi_{1}^{-} p",2000,500,2500);
  p_em_mass=new TH1F("p_em_mass","Invariant mass #pi_{2}^{-} p",2000,500,2500);
  pim_ep_mass=new TH1F("pim_ep_mass","Invariant mass #pi^{-} #pi^{+} ",2000,200,1500);
  em_ep_mass=new TH1F("em_ep_mass","Invariant mass #pi_{2}^{-} #pi^{+} ",2000,200,1500);
  pim_ep_mass=new TH1F("pim_ep_mass","Invariant mass #pi_{1}^{-} #pi^{+} ",2000,200,1500);
  p_pim_ep_pim_mass=new TH1F("p_pim_ep_pim_mass","Invariant mass #pi^{-} #pi^{+} p #pi^{-}",2000,1000,2000);

  dist_p_pim_pim_ep=new TH2F("dist_p_pim_pim_ep","dist_p_pim_pim_ep",300,0,300,300,0,300);
  dist_p_pim=new TH1F("dist_p_pim","dist_p_pim",1000,0,300);
  dist_pim_ep=new TH1F("dist_ep_pim","dist_ep_pim",1000,0,300);

  
  DL_p_pim_mass=new TH1F("DL_p_pim_mass","Invariant mass #pi^{-} p",2000,500,2500);						       
  DL_p_pim_mass=new TH1F("DL_p_pim_mass","DL_Invariant mass #pi_{1}^{-} p",2000,500,2500);
  DL_p_em_mass=new TH1F("DL_p_em_mass","DL_Invariant mass #pi_{2}^{-} p",2000,500,2500);
  DL_pim_ep_mass=new TH1F("DL_pim_ep_mass","DL_Invariant mass #pi^{-} #pi^{+} ",2000,200,1500);
  DL_em_ep_mass=new TH1F("DL_em_ep_mass","DL_Invariant mass #pi_{2}^{-} #pi^{+} ",2000,200,1500);
  DL_pim_ep_mass=new TH1F("DL_pim_ep_mass","DL_Invariant mass #pi_{1}^{-} #pi^{+} ",2000,200,1500);
  DL_p_pim_ep_pim_mass=new TH1F("DL_p_pim_ep_pim_mass","DL_Invariant mass #pi^{-} #pi^{+} p #pi^{-}",2000,1000,2000);

  DL_dist_p_pim_pim_ep=new TH2F("DL_dist_p_pim_pim_ep","DL_dist_p_pim_pim_ep",300,0,300,300,0,300);
  DL_dist_p_pim=new TH1F("DL_dist_p_pim","DL_dist_p_pim",1000,0,300);
  DL_dist_pim_ep=new TH1F("DL_dist_ep_pim","DL_dist_ep_pim",1000,0,300);

  DML_p_pim_mass=new TH1F("DML_p_pim_mass","Invariant mass #pi^{-} p",2000,500,2500);						       
  DML_p_pim_mass=new TH1F("DML_p_pim_mass","DML_Invariant mass #pi_{1}^{-} p",2000,500,2500);
  DML_p_em_mass=new TH1F("DML_p_em_mass","DML_Invariant mass #pi_{2}^{-} p",2000,500,2500);
  DML_pim_ep_mass=new TH1F("DML_pim_ep_mass","DML_Invariant mass #pi^{-} #pi^{+} ",2000,200,1500);
  DML_em_ep_mass=new TH1F("DML_em_ep_mass","DML_Invariant mass #pi_{2}^{-} #pi^{+} ",2000,200,1500);
  DML_pim_ep_mass=new TH1F("DML_pim_ep_mass","DML_Invariant mass #pi_{1}^{-} #pi^{+} ",2000,200,1500);
  DML_p_pim_ep_pim_mass=new TH1F("DML_p_pim_ep_pim_mass","DML_Invariant mass #pi^{-} #pi^{+} p #pi^{-}",2000,1000,2000);

  DML_dist_p_pim_pim_ep=new TH2F("DML_dist_p_pim_pim_ep","DML_dist_p_pim_pim_ep",300,0,300,300,0,300);
  DML_dist_p_pim=new TH1F("DML_dist_p_pim","DML_dist_p_pim",1000,0,300);
  DML_dist_pim_ep=new TH1F("DML_dist_ep_pim","DML_dist_ep_pim",1000,0,300);

  //miss_energy=new TH1F("miss_energy","missing energy",2000,0,6000);
  //DL_miss_energy=new TH1F("DL_miss_energy","missing energy",2000,0,6000);
  //DML_miss_energy=new TH1F("DLM_miss_energy","missing energy",2000,0,6000);

  //ppim_eppim_mass=new TH2F("ppim_eppim_mass","ppim_eppim_mass;M_{p #pi^{-}}[MeV];M_{#pi^{+} #pi^{-}}[MeV]",600,1000,2200,400,200,1000);
  //dist_z_ppim_eppim_mass=new TH2F("dist_z_ppim_eppim_mass","dist_z_ppim_eppim_mass;M_{p #pi^{-}}[MeV];M_{#pi^{+} #pi^{-}}[MeV]",600,1000,2200,400,200,1000);

  sum_dist_1=new TH1F("sum_dist_1","Sum of all distances in hyp1",1000,0,1000);
  sum_dist_2=new TH1F("sum_dist_2","Sum of all distances in hyp2",1000,0,1000);
  sum_dist_diff=new TH1F("sum_dist_diff","Difference between hypothesis 1 and 2",600,0,600);

  chi_p_pim_mass=new TH1F("chi_p_pim_mass","p #pi^{-} mass after chi cut; M_{#pi^{-} p}[MeV]",600,1000,1600);
  chi_ep_pim_mass=new TH1F("chi_ep_pim_mass","#pi^{+} #pi^{-} mass after chi cut; M_{#pi^{-} #pi^{+}}[MeV]",900,200,1100);
  chi_lambda_vertex=new TH2F("chi_lambda_vertex","Reconstructed vertex of a #Lambda (1116);Z_{vertex}[mm];R_{vertex}[mm]",600,-100,200,400,0,200);
  chi_final_mass=new TH1F("chi_final_mass","#pi^{+} #pi^{-} P #Pi^{-} mass after chi cut; M_{#pi^{-} #pi^{+}  p #pi^{-}}[MeV]",1000,1300,2300);

  
  LM_chi_p_pim_mass=new TH1F("LM_chi_p_pim_mass","p #pi^{-} mass after chi cut and #Lambda mass; M_{#pi^{-} p}[MeV]",600,1000,1600);
  LM_chi_ep_pim_mass=new TH1F("LM_chi_ep_pim_mass","#pi^{+} #pi^{-} mass after chi cut and #Lambda mass; M_{#pi^{-} #pi^{+}}[MeV]",900,200,1100);
  LM_chi_lambda_vertex=new TH2F("LM_chi_lambda_vertex","Reconstructed vertex of a #Lambda (1116),after chi and #Lambda mass cut;Z_{vertex}[mm];R_{vertex}[mm]",600,-100,200,400,0,200);
  LM_chi_final_mass=new TH1F("LM_chi_final_mass","#pi^{+} #pi^{-} P #Pi^{-} mass after chi and mass cut; M_{#pi^{-} #pi^{+}  p #pi^{-}}[MeV]",1000,1300,2300);

  int chi_step=50;
  int dist_step=2;

  for(int i=0;i<10;i++)
    for(int j=0;j<10;j++)
      {
	char hname[20];
	char htitle[40];
	sprintf(htitle,"#pi^{+} #pi^{-} p #pi^{-} for chi < %d and vertex_distance > %d",250+chi_step*i,dist_step*j);
	sprintf(hname,"chi_%d_min_distance_%d",250+chi_step*i,dist_step*j);
	signal_fit[i][j]=new TH1F(hname,htitle,500,1300,2300);
      }

  */
  /**************************** M A I N   P A R T ****************************************/

  //PPim t;
  PPimEpEm t2;
  PPimEpEp t3;
  PPimEmEm t4;
  //cout << "START PPIM!" << endl;
  //t.Loop();
  //cout << "STOP PPIM!" << endl;

  cout << "START PPimEpEm!" << endl;
  t2.Loop();
  cout << "STOP PPimEpEm!!" << endl;
  cout << "START PPimEpEp!" << endl;
  t3.Loop();
  cout << "STOP PPimEpEp!!" << endl;
  cout << "START PPimEmEm!" << endl;
  t4.Loop();
  cout << "STOP PPimEmEm!!" << endl;
  
  /*EpEp t_back1;
    t_back1.Loop();
    cout << "START EPEP!" << endl;
    EmEm t_back2;
    t_back2.Loop();   
    cout << "START EMEM!" << endl;*/
  cout << "FINALIZING!" << endl;


  /***************************** F I N A L     C A L C U L A T I O N S ********************/


  /****************************************************************************************/

  outFileData->cd();

  tlo->Write();
  tlo_epep->Write();
  tlo_emem->Write();
  /*
  p_p_beta->Write();
  pim_p_beta->Write();
  ep_p_beta->Write();
  p_pim_mass->Write();
  p_pim_mass->Write();
  p_em_mass->Write();
  pim_ep_mass->Write();
  em_ep_mass->Write();
  pim_ep_mass->Write();
  p_pim_ep_pim_mass->Write();
  //miss_energy->Write();
  
  dist_p_pim_pim_ep->Write();
  dist_pim_ep->Write();
  dist_p_pim->Write();

  DL_p_pim_mass->Write();
  DL_p_pim_mass->Write();
  DL_p_em_mass->Write();
  DL_pim_ep_mass->Write();
  DL_em_ep_mass->Write();
  DL_pim_ep_mass->Write();
  DL_p_pim_ep_pim_mass->Write();

  DL_dist_p_pim_pim_ep->Write();
  DL_dist_pim_ep->Write();
  DL_dist_p_pim->Write();
  //DL_miss_energy->Write();

  DML_p_pim_mass->Write();
  DML_p_pim_mass->Write();
  DML_p_em_mass->Write();
  DML_pim_ep_mass->Write();
  DML_em_ep_mass->Write();
  DML_pim_ep_mass->Write();
  DML_p_pim_ep_pim_mass->Write();

  DML_dist_p_pim_pim_ep->Write();
  DML_dist_pim_ep->Write();
  DML_dist_p_pim->Write();
  //DML_miss_energy->Write();

  //dist_z_ppim_eppim_mass->Write();
  //ppim_eppim_mass->Write();

  //p_mass->Write();
  //pim_mass->Write();
 
  //myfile.close();

  
  sum_dist_diff->Write();
  sum_dist_2->Write();
  sum_dist_1->Write();

  chi_ep_pim_mass->Write();
  chi_p_pim_mass->Write();
  chi_lambda_vertex->Write();
  chi_final_mass->Write();
  
  LM_chi_ep_pim_mass->Write();
  LM_chi_p_pim_mass->Write();
  LM_chi_lambda_vertex->Write();
  LM_chi_final_mass->Write();
  //myfile.close();

  for(int i=0;i<10;i++)
    for(int j=0;j<10;j++)
      {
	signal_fit[i][j]->Write();
      }
  */
  outFileData->Close();
}

