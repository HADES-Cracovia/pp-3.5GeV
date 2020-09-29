#include "PPim.h"
#include "PPimPipPim.h"
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
  double proton_momentum=TMath::Sqrt((3.5+0.93827)*(3.5+0.93827)-0.93827*0.93827)*1000;
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
  pim1 = new TLorentzVector(0,0,0,0);
  pim2 = new TLorentzVector(0,0,0,0);
  pip = new TLorentzVector(0,0,0,0);
  pi = new TLorentzVector(0,0,0,0);
  gammappi = new TLorentzVector(0,0,0,0);
  gammappip = new TLorentzVector(0,0,0,0);
  gammappim1 = new TLorentzVector(0,0,0,0);
  gammappim2 = new TLorentzVector(0,0,0,0);
  gammapim1pip = new TLorentzVector(0,0,0,0);
  gammapim2pip = new TLorentzVector(0,0,0,0);
  gammappim1pip= new TLorentzVector(0,0,0,0);
  gammappim2pip= new TLorentzVector(0,0,0,0);
  //gammappimpip= new TLorentzVector(0,0,0,0);
  gammappim1pim2=new TLorentzVector(0,0,0,0);
  gammappim1pippim2 = new TLorentzVector(0,0,0,0);

  ppi = new TLorentzVector(0,0,0,0);
  p_delta = new TLorentzVector(0,0,0,0);
  pi_delta = new TLorentzVector(0,0,0,0);
  ppi_miss = new TLorentzVector(0,0,0,0);

  //reader->BookMVA("kMLP","/lustre/nyx/hades/user/knowakow/PP/FAT/TMVA/weights/TMVAClassification_from_simplus_rec_cuts_kMLP_ce_600_n4_no_ev.weights.xml" );
  /************************************** O U T P U T   F I L E ******************************************/
  //outFileData = new TFile("S1385pK0_Rafal_part2.root","recreate");
  //outFileData = new TFile("SDppK0_Rafal_part2.root","recreate");
  //outFileData = new TFile("LDppK0_Rafal_part2.root","recreate");
  outFileData=new TFile("pp_Lpippim_ver4_new_vertex.root","recreate");
  //outFileData=new TFile("pp_pK0Lpip_ver2.root","recreate");
  //outFileData=new TFile("DppPPimPipPim.root","recreate");
  //outFileData=new TFile("PPipPPimPipPim.root","recreate");
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
  n_out = new HNtuple("ppimpippim","ppimpippim");
  //n_ppim = new HNtuple("ppim","ppim");
  n_out->setFile( outFileData );
  //n_ppim->setFile( outFileData );
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
  pip_p_beta=new TH2F("pip_p_beta","Momentum vs. beta for #pi^{+}",p_n,p_min,p_max,beta_n,beta_min,beta_max);
  p_pim_mass=new TH1F("p_pim_mass","Invariant mass #pi^{-} p",2000,500,2500);									       
  p_pim1_mass=new TH1F("p_pim1_mass","Invariant mass #pi_{1}^{-} p",2000,500,2500);
  p_pim2_mass=new TH1F("p_pim2_mass","Invariant mass #pi_{2}^{-} p",2000,500,2500);
  pim_pip_mass=new TH1F("pim_pip_mass","Invariant mass #pi^{-} #pi^{+} ",2000,200,1500);
  pim2_pip_mass=new TH1F("pim2_pip_mass","Invariant mass #pi_{2}^{-} #pi^{+} ",2000,200,1500);
  pim1_pip_mass=new TH1F("pim1_pip_mass","Invariant mass #pi_{1}^{-} #pi^{+} ",2000,200,1500);
  p_pim_pip_pim_mass=new TH1F("p_pim_pip_pim_mass","Invariant mass #pi^{-} #pi^{+} p #pi^{-}",2000,1000,2000);

  dist_p_pim_pim_pip=new TH2F("dist_p_pim_pim_pip","dist_p_pim_pim_pip",300,0,300,300,0,300);
  dist_p_pim=new TH1F("dist_p_pim","dist_p_pim",1000,0,300);
  dist_pim_pip=new TH1F("dist_pip_pim","dist_pip_pim",1000,0,300);

  
  DL_p_pim_mass=new TH1F("DL_p_pim_mass","Invariant mass #pi^{-} p",2000,500,2500);						       
  DL_p_pim1_mass=new TH1F("DL_p_pim1_mass","DL_Invariant mass #pi_{1}^{-} p",2000,500,2500);
  DL_p_pim2_mass=new TH1F("DL_p_pim2_mass","DL_Invariant mass #pi_{2}^{-} p",2000,500,2500);
  DL_pim_pip_mass=new TH1F("DL_pim_pip_mass","DL_Invariant mass #pi^{-} #pi^{+} ",2000,200,1500);
  DL_pim2_pip_mass=new TH1F("DL_pim2_pip_mass","DL_Invariant mass #pi_{2}^{-} #pi^{+} ",2000,200,1500);
  DL_pim1_pip_mass=new TH1F("DL_pim1_pip_mass","DL_Invariant mass #pi_{1}^{-} #pi^{+} ",2000,200,1500);
  DL_p_pim_pip_pim_mass=new TH1F("DL_p_pim_pip_pim_mass","DL_Invariant mass #pi^{-} #pi^{+} p #pi^{-}",2000,1000,2000);

  DL_dist_p_pim_pim_pip=new TH2F("DL_dist_p_pim_pim_pip","DL_dist_p_pim_pim_pip",300,0,300,300,0,300);
  DL_dist_p_pim=new TH1F("DL_dist_p_pim","DL_dist_p_pim",1000,0,300);
  DL_dist_pim_pip=new TH1F("DL_dist_pip_pim","DL_dist_pip_pim",1000,0,300);

  DML_p_pim_mass=new TH1F("DML_p_pim_mass","Invariant mass #pi^{-} p",2000,500,2500);						       
  DML_p_pim1_mass=new TH1F("DML_p_pim1_mass","DML_Invariant mass #pi_{1}^{-} p",2000,500,2500);
  DML_p_pim2_mass=new TH1F("DML_p_pim2_mass","DML_Invariant mass #pi_{2}^{-} p",2000,500,2500);
  DML_pim_pip_mass=new TH1F("DML_pim_pip_mass","DML_Invariant mass #pi^{-} #pi^{+} ",2000,200,1500);
  DML_pim2_pip_mass=new TH1F("DML_pim2_pip_mass","DML_Invariant mass #pi_{2}^{-} #pi^{+} ",2000,200,1500);
  DML_pim1_pip_mass=new TH1F("DML_pim1_pip_mass","DML_Invariant mass #pi_{1}^{-} #pi^{+} ",2000,200,1500);
  DML_p_pim_pip_pim_mass=new TH1F("DML_p_pim_pip_pim_mass","DML_Invariant mass #pi^{-} #pi^{+} p #pi^{-}",2000,1000,2000);

  DML_dist_p_pim_pim_pip=new TH2F("DML_dist_p_pim_pim_pip","DML_dist_p_pim_pim_pip",300,0,300,300,0,300);
  DML_dist_p_pim=new TH1F("DML_dist_p_pim","DML_dist_p_pim",1000,0,300);
  DML_dist_pim_pip=new TH1F("DML_dist_pip_pim","DML_dist_pip_pim",1000,0,300);

  //miss_energy=new TH1F("miss_energy","missing energy",2000,0,6000);
  //DL_miss_energy=new TH1F("DL_miss_energy","missing energy",2000,0,6000);
  //DML_miss_energy=new TH1F("DLM_miss_energy","missing energy",2000,0,6000);

  //ppim_pippim_mass=new TH2F("ppim_pippim_mass","ppim_pippim_mass;M_{p #pi^{-}}[MeV];M_{#pi^{+} #pi^{-}}[MeV]",600,1000,2200,400,200,1000);
  //dist_z_ppim_pippim_mass=new TH2F("dist_z_ppim_pippim_mass","dist_z_ppim_pippim_mass;M_{p #pi^{-}}[MeV];M_{#pi^{+} #pi^{-}}[MeV]",600,1000,2200,400,200,1000);

  sum_dist_1=new TH1F("sum_dist_1","Sum of all distances in hyp1",1000,0,1000);
  sum_dist_2=new TH1F("sum_dist_2","Sum of all distances in hyp2",1000,0,1000);
  sum_dist_diff=new TH1F("sum_dist_diff","Difference between hypothesis 1 and 2",600,0,600);

  chi_p_pim_mass=new TH1F("chi_p_pim_mass","p #pi^{-} mass after chi cut; M_{#pi^{-} p}[MeV]",600,1000,1600);
  chi_pip_pim_mass=new TH1F("chi_pip_pim_mass","#pi^{+} #pi^{-} mass after chi cut; M_{#pi^{-} #pi^{+}}[MeV]",900,200,1100);
  chi_lambda_vertex=new TH2F("chi_lambda_vertex","Reconstructed vertex of a #Lambda (1116);Z_{vertex}[mm];R_{vertex}[mm]",600,-100,200,400,0,200);
  chi_final_mass=new TH1F("chi_final_mass","#pi^{+} #pi^{-} P #Pi^{-} mass after chi cut; M_{#pi^{-} #pi^{+}  p #pi^{-}}[MeV]",1000,1300,2300);

  
  LM_chi_p_pim_mass=new TH1F("LM_chi_p_pim_mass","p #pi^{-} mass after chi cut and #Lambda mass; M_{#pi^{-} p}[MeV]",600,1000,1600);
  LM_chi_pip_pim_mass=new TH1F("LM_chi_pip_pim_mass","#pi^{+} #pi^{-} mass after chi cut and #Lambda mass; M_{#pi^{-} #pi^{+}}[MeV]",900,200,1100);
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
  PPimPipPim t2;
  //cout << "START PPIM!" << endl;
  //t.Loop();
  //cout << "STOP PPIM!" << endl;

  cout << "START PPimPipPim!" << endl;
  t2.Loop();
  cout << "STOP PPimPipPim!!" << endl;
  
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

  if(n_out!=0)
    n_out->Write();
  else
    cout<<"n_out pointer empty"<<endl;

  if(n_ppim!=0)
    n_ppim->Write();
  else
    cout<<"n_out pointer empty"<<endl;
  /*
  p_p_beta->Write();
  pim_p_beta->Write();
  pip_p_beta->Write();
  p_pim_mass->Write();
  p_pim1_mass->Write();
  p_pim2_mass->Write();
  pim_pip_mass->Write();
  pim2_pip_mass->Write();
  pim1_pip_mass->Write();
  p_pim_pip_pim_mass->Write();
  //miss_energy->Write();
  
  dist_p_pim_pim_pip->Write();
  dist_pim_pip->Write();
  dist_p_pim->Write();

  DL_p_pim_mass->Write();
  DL_p_pim1_mass->Write();
  DL_p_pim2_mass->Write();
  DL_pim_pip_mass->Write();
  DL_pim2_pip_mass->Write();
  DL_pim1_pip_mass->Write();
  DL_p_pim_pip_pim_mass->Write();

  DL_dist_p_pim_pim_pip->Write();
  DL_dist_pim_pip->Write();
  DL_dist_p_pim->Write();
  //DL_miss_energy->Write();

  DML_p_pim_mass->Write();
  DML_p_pim1_mass->Write();
  DML_p_pim2_mass->Write();
  DML_pim_pip_mass->Write();
  DML_pim2_pip_mass->Write();
  DML_pim1_pip_mass->Write();
  DML_p_pim_pip_pim_mass->Write();

  DML_dist_p_pim_pim_pip->Write();
  DML_dist_pim_pip->Write();
  DML_dist_p_pim->Write();
  //DML_miss_energy->Write();

  //dist_z_ppim_pippim_mass->Write();
  //ppim_pippim_mass->Write();

  //p_mass->Write();
  //pim_mass->Write();
 
  //myfile.close();

  
  sum_dist_diff->Write();
  sum_dist_2->Write();
  sum_dist_1->Write();

  chi_pip_pim_mass->Write();
  chi_p_pim_mass->Write();
  chi_lambda_vertex->Write();
  chi_final_mass->Write();
  
  LM_chi_pip_pim_mass->Write();
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

