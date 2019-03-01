#include "PPimPipPim.h"
#include "PPimPipPim_buffer.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include "data.h"
#include <iostream>
#include <TLorentzVector.h>
#include <TVector3.h>
#include "hntuple.h"


using namespace std;
using namespace PATData;

void PPimPipPim::Loop()
{
  int licznik = 0;
  if (fChain == 0) return;

  int nentries = fChain->GetEntriesFast();

  double nbytes = 0, nb = 0;

  std::vector< PPimPipPim_ID_buffer >  buffer;
  event_number=-1;
  event_mult=1;

  std::vector<double> the_best;
    
  for (Long64_t jentry=0; jentry<nentries;jentry++)
    {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;

      ++licznik;
      if ((licznik % 100000)==0) cout << "Events: " << licznik << " "<<(1.0*licznik)/(1.0*nentries)*100<<" %"<< endl;


      if(isBest>=0 /*&& trigdownscaleflag==1*/)
	{
	  //reset event multiplisity in case of a new event 
	  if(event!=event_number)
	    {
	      if(event_number!=-1)
		{
		  for(int k=0;k<buffer.size();k++)
		    {
		      double min_value=1000;
		      int best_hipo;
		      for(int j=0;j<the_best.size();j++)
			{
			  cout<<the_best[j]<<" ";
			  if(the_best[j]<=min_value)
			    {
			      min_value=the_best[j];
			      best_hipo=j;
			    }
			}
		      cout<<endl;
		      cout<<"best hipo "<<best_hipo<<endl;
		      if(k==best_hipo)
			filler(buffer[k],event_mult,1,1);
		      else
			filler(buffer[k],event_mult,1,0);
			
		      //cout<<"end of ev. "<<event_number<<" mult "<<event_mult<<endl;
		    }
		  buffer.clear();
		  the_best.clear();
		}
	      event_number=event;
	      event_mult=1;
	    }
	  else
	    {
	      event_mult++;
	    }

	  double F = 1.006;
	  TVector3 v1, v2, v3, v4, v5;
	  v2.SetXYZ(F*p_p*sin(D2R*p_theta)*cos(D2R*p_phi),F*p_p*sin(D2R*p_theta)*sin(D2R*p_phi),F*p_p*cos(D2R*p_theta));
	  v3.SetXYZ(F*pim1_p*sin(D2R*pim1_theta)*cos(D2R*pim1_phi),F*pim1_p*sin(D2R*pim1_theta)*sin(D2R*pim1_phi),F*pim1_p*cos(D2R*pim1_theta));
	  v4.SetXYZ(F*pip_p*sin(D2R*pip_theta)*cos(D2R*pip_phi),F*pip_p*sin(D2R*pip_theta)*sin(D2R*pip_phi),F*pip_p*cos(D2R*pip_theta));
	  v5.SetXYZ(F*pim2_p*sin(D2R*pim2_theta)*cos(D2R*pim2_phi),F*pim2_p*sin(D2R*pim2_theta)*sin(D2R*pim2_phi),F*pim2_p*cos(D2R*pim2_theta));

	  p->SetVectM( v2, 938.272013 );
	  pim1->SetVectM( v3, 139.57018 );
	  pip->SetVectM( v4, 139.57018 );
	  pim2->SetVectM( v5, 139.57018 );

	  
	  double m_inv_ppim1 = gammappim1->M();
	  double m_inv_ppim2 = gammappim2->M();
	  double m_inv_pippim1 = gammapim1pip->M();
	  double m_inv_pippim2 = gammapim2pip->M();
	  double m_inv_ppimpippim = gammappim1pippim2->M();
	  double oa = R2D * openingangle(*p, *pim1);
	  //double oa_rich = R2D * openingangle(r1, r2);

	  double p_mass = p_p*p_p * (  1. / (p_beta*p_beta)  - 1. ) ;
	  double pi_mass = pim1_p*pim1_p * (  1. / (pim1_beta*pim1_beta_new)  - 1. ) ;
	  double pip_mass = pip_p*pip_p * (  1. / (pip_beta*pip_beta_new)  - 1. ) ;
	  double pim1_mass = pim1_p*pim1_p * (  1. / (pim1_beta*pim1_beta_new)  - 1. ) ;
	  double pim2_mass = pim2_p*pim2_p * (  1. / (pim2_beta*pim2_beta_new)  - 1. ) ;

	  TVector3 ver_p_pim1=vertex(p_r,p_z,v2,pim1_r,pim1_z,v3);
	  TVector3 ver_p_pim2=vertex(p_r,p_z,v2,pim2_r,pim2_z,v5);
	  TVector3 ver_pip_pim1=vertex(pip_r,pip_z,v4,pim1_r,pim1_z,v3);
	  TVector3 ver_pip_pim2=vertex(pip_r,pip_z,v4,pim2_r,pim2_z,v5);

	  TVector3 ver_to_ver_1=ver_p_pim1-ver_pip_pim2;
	  TVector3 ver_to_ver_2=ver_p_pim2-ver_pip_pim1;

	  double oa_pim1_p=R2D*openingangle(pim1->Vect(),p->Vect());
	  double oa_pim2_p=R2D*openingangle(pim2->Vect(),p->Vect());
	  double oa_pip_p=R2D*openingangle(pip->Vect(),p->Vect());
	  double oa_pim1_pim2=R2D*openingangle(pim1->Vect(),pim2->Vect());
	  double oa_pim1_pip=R2D*openingangle(pim1->Vect(),pip->Vect());
	  double oa_pim2_pip=R2D*openingangle(pim2->Vect(),pip->Vect());
                  
	  double oa_lambda_1=R2D*openingangle(ver_to_ver_1,gammappim1->Vect());
	  double oa_lambda_2=R2D*openingangle(ver_to_ver_2,gammappim2->Vect());
      
	  double dist_p_pim1=trackDistance(p_r,p_z,v2,pim1_r,pim1_z,v3);
	  double dist_p_pim2=trackDistance(p_r,p_z,v2,pim2_r,pim2_z,v5);
	  double dist_pip_pim1=trackDistance(pip_r,pip_z,v4,pim1_r,pim1_z,v3);
	  double dist_pip_pim2=trackDistance(pip_r,pip_z,v4,pim2_r,pim2_z,v5);
	  double dist_lambda1_pip=trackDistance(pip_r,pip_z,v4,ver_pip_pim1.Z(),getR(ver_pip_pim1),gammappim1->Vect());
	  double dist_lambda2_pip=trackDistance(pip_r,pip_z,v4,ver_pip_pim2.Z(),getR(ver_pip_pim2),gammappim2->Vect());
	  double dist_lambda1_pim2=trackDistance(pim2_r,pim2_z,v5,ver_pip_pim1.Z(),getR(ver_pip_pim1),gammappim1->Vect());
	  double dist_lambda2_pim1=trackDistance(pim1_r,pim1_z,v3,ver_pip_pim2.Z(),getR(ver_pip_pim2),gammappim2->Vect());
	  double dist_ver_to_ver_1=ver_to_ver_1.Mag();
	  double dist_ver_to_ver_2=ver_to_ver_2.Mag();
      
	  //double quality=trackDistance(p_r,p_z,v2,pim1_r,pim1_z,v3);
	  double quality1=dist_p_pim1*dist_p_pim1+dist_pip_pim2*dist_pip_pim2;
	  double quality2=dist_p_pim2*dist_p_pim2+dist_pip_pim1*dist_pip_pim1;

	  double quality=std::min(quality1,quality2);
	  
	  //add hypothesis to buffer and quality measure
	  buffer.push_back( PPimPipPim_ID_buffer( this ));
	  the_best.push_back(quality);
	  	
	}
	  
    }
}

PPimPipPim::PPimPipPim(TTree *tree)
{
  // if parameter tree is not specified (or zero), connect the file
  // used to generate this class and read the Tree.
  /*if (tree == 0) {
    TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("be08280235056_dst_gen1_sep08_hadron_out.root");
    if (!f || !f->IsOpen()) {
    f = new TFile("be08280235056_dst_gen1_sep08_hadron_out.root");
    }
    f->GetObject("PPimPipPim",tree);
    }
    Init(tree);*/
  if (tree == 0)
    {
      TChain * chain = new TChain("PPimPipPim_ID","");
      //chain->Add("/lustre/nyx/hades/user/knowakow/PP/PAT_sim/FILES/pip_pim/all.root");
      //chain->Add("/lustre/nyx/hades/user/knowakow/PP/PAT_sim/FILES/pip_pim_ver2/all.root/PPimPipPim_ID");
      chain->Add("/lustre/nyx/hades/user/knowakow/PP/PAT_sim/FILES/pippimL/all.root/PPimPipPim_ID");
      /*
	chain->Add("/lustre/nyx/hades/user/knowakow/PP/PAT_sim/FILES/lambda1520_100k1_dst_hadron_out.root/PPimPipPim_ID");
	chain->Add("/lustre/nyx/hades/user/knowakow/PP/PAT_sim/FILES/lambda1520_100k2_dst_hadron_out.root/PPimPipPim_ID");
	chain->Add("/lustre/nyx/hades/user/knowakow/PP/PAT_sim/FILES/lambda1520_100k3_dst_hadron_out.root/PPimPipPim_ID");
	chain->Add("/lustre/nyx/hades/user/knowakow/PP/PAT_sim/FILES/lambda1520_100k4_dst_hadron_out.root/PPimPipPim_ID");
	chain->Add("/lustre/nyx/hades/user/knowakow/PP/PAT_sim/FILES/lambda1520_100k5_dst_hadron_out.root/PPimPipPim_ID");
	chain->Add("/lustre/nyx/hades/user/knowakow/PP/PAT_sim/FILES/lambda1520_100k6_dst_hadron_out.root/PPimPipPim_ID");
	chain->Add("/lustre/nyx/hades/user/knowakow/PP/PAT_sim/FILES/lambda1520_100k7_dst_hadron_out.root/PPimPipPim_ID");
	chain->Add("/lustre/nyx/hades/user/knowakow/PP/PAT_sim/FILES/lambda1520_100k8_dst_hadron_out.root/PPimPipPim_ID");
	chain->Add("/lustre/nyx/hades/user/knowakow/PP/PAT_sim/FILES/lambda1520_100k9_dst_hadron_out.root/PPimPipPim_ID");
	/*
	chain->Add("/lustre/nyx/hades/user/knowakow/PP/PAT_1/FILES/ppimpippim_dedx_3/hadron01.root/PPimPipPim_ID");
	chain->Add("/lustre/nyx/hades/user/knowakow/PP/PAT_1/FILES/ppimpippim_dedx_3/hadron02.root/PPimPipPim_ID");
	chain->Add("/lustre/nyx/hades/user/knowakow/PP/PAT_1/FILES/ppimpippim_dedx_3/hadron03.root/PPimPipPim_ID");
	chain->Add("/lustre/nyx/hades/user/knowakow/PP/PAT_1/FILES/ppimpippim_dedx_3/hadron04.root/PPimPipPim_ID");
	chain->Add("/lustre/nyx/hades/user/knowakow/PP/PAT_1/FILES/ppimpippim_dedx_3/hadron05.root/PPimPipPim_ID");
	chain->Add("/lustre/nyx/hades/user/knowakow/PP/PAT_1/FILES/ppimpippim_dedx_3/hadron06.root/PPimPipPim_ID");
	chain->Add("/lustre/nyx/hades/user/knowakow/PP/PAT_1/FILES/ppimpippim_dedx_3/hadron07.root/PPimPipPim_ID");
	chain->Add("/lustre/nyx/hades/user/knowakow/PP/PAT_1/FILES/ppimpippim_dedx_3/hadron08.root/PPimPipPim_ID");
	chain->Add("/lustre/nyx/hades/user/knowakow/PP/PAT_1/FILES/ppimpippim_dedx_3/hadron09.root/PPimPipPim_ID");
	chain->Add("/lustre/nyx/hades/user/knowakow/PP/PAT_1/FILES/ppimpippim_dedx_3/hadron10.root/PPimPipPim_ID");
	chain->Add("/lustre/nyx/hades/user/knowakow/PP/PAT_1/FILES/ppimpippim_dedx_3/hadron11.root/PPimPipPim_ID");    //
      */
      tree = chain;
    }

  Init(tree);
}

PPimPipPim::~PPimPipPim()
{
  if (!fChain) return;
  delete fChain->GetCurrentFile();
}

void PPimPipPim::filler( const PPimPipPim_ID_buffer& s,int event_mult, double WEIGHT,int is_best)
{
  
  double F = 1.006;
  TVector3 v1, v2, v3, v4, v5;
  v2.SetXYZ(F*s.p_p*sin(D2R*s.p_theta)*cos(D2R*s.p_phi),F*s.p_p*sin(D2R*s.p_theta)*sin(D2R*s.p_phi),F*s.p_p*cos(D2R*s.p_theta));
  v3.SetXYZ(F*s.pim1_p*sin(D2R*s.pim1_theta)*cos(D2R*s.pim1_phi),F*s.pim1_p*sin(D2R*s.pim1_theta)*sin(D2R*s.pim1_phi),F*s.pim1_p*cos(D2R*s.pim1_theta));
  v4.SetXYZ(F*s.pip_p*sin(D2R*s.pip_theta)*cos(D2R*s.pip_phi),F*s.pip_p*sin(D2R*s.pip_theta)*sin(D2R*s.pip_phi),F*s.pip_p*cos(D2R*s.pip_theta));
  v5.SetXYZ(F*s.pim2_p*sin(D2R*s.pim2_theta)*cos(D2R*s.pim2_phi),F*s.pim2_p*sin(D2R*s.pim2_theta)*sin(D2R*s.pim2_phi),F*s.pim2_p*cos(D2R*s.pim2_theta));

  /*TVector3 r1, r2, r3,r4;
    r1.SetXYZ(sin(D2R*p_theta_rich)*cos(D2R*p_phi_rich),sin(D2R*p_theta_rich)*sin(D2R*p_phi_rich),cos(D2R*p_theta_rich));
    r2.SetXYZ(sin(D2R*pim1_theta_rich)*cos(D2R*pim1_phi_rich),sin(D2R*pim1_theta_rich)*sin(D2R*pim1_phi_rich),cos(D2R*pim1_theta_rich));
    r3.SetXYZ(sin(D2R*pip_theta_rich)*cos(D2R*pip_phi_rich),sin(D2R*pip_theta_rich)*sin(D2R*pip_phi_rich),cos(D2R*pip_theta_rich));
    r4.SetXYZ(sin(D2R*pim2_theta_rich)*cos(D2R*pim2_phi_rich),sin(D2R*pim2_theta_rich)*sin(D2R*pim2_phi_rich),cos(D2R*pim2_theta_rich));
  */      
  p->SetVectM( v2, 938.272013 );
  pim1->SetVectM( v3, 139.57018 );
  pip->SetVectM( v4, 139.57018 );
  pim2->SetVectM( v5, 139.57018 );

  *gammappim1 = *p + *pim1;
  *gammappim2 = *p + *pim2;
  *gammapim1pip= *pim1 + *pip;
  *gammapim2pip= *pim2 + *pip;
  *gammappim1pippim2=*pim1 +*pim2 + *pip + *p;
  *miss=*beam-*gammappim1pippim2;
      
  //*ppim1 = *p + *pim1;
  //*p_delta = *p;
  //*pim1_delta = *pim1;
  //*ppim1_miss = *beam - *p - *pim1;

  double m_inv_ppim1 = gammappim1->M();
  double m_inv_ppim2 = gammappim2->M();
  double m_inv_pippim1 = gammapim1pip->M();
  double m_inv_pippim2 = gammapim2pip->M();
  double m_inv_ppimpippim = gammappim1pippim2->M();
  double oa = R2D * openingangle(*p, *pim1);
  //double oa_rich = R2D * openingangle(r1, r2);

  double p_mass = s.p_p*s.p_p * (  1. / (s.p_beta*s.p_beta)  - 1. ) ;
  double pi_mass = s.pim1_p*s.pim1_p * (  1. / (s.pim1_beta*s.pim1_beta_new)  - 1. ) ;
  double pip_mass = s.pip_p*s.pip_p * (  1. / (s.pip_beta*s.pip_beta_new)  - 1. ) ;
  double pim1_mass = s.pim1_p*s.pim1_p * (  1. / (s.pim1_beta*s.pim1_beta_new)  - 1. ) ;
  double pim2_mass = s.pim2_p*s.pim2_p * (  1. / (s.pim2_beta*s.pim2_beta_new)  - 1. ) ;

  TVector3 ver_p_pim1=vertex(s.p_r,s.p_z,v2,s.pim1_r,s.pim1_z,v3);
  TVector3 ver_p_pim2=vertex(s.p_r,s.p_z,v2,s.pim2_r,s.pim2_z,v5);
  TVector3 ver_pip_pim1=vertex(s.pip_r,s.pip_z,v4,s.pim1_r,s.pim1_z,v3);
  TVector3 ver_pip_pim2=vertex(s.pip_r,s.pip_z,v4,s.pim2_r,s.pim2_z,v5);

  TVector3 ver_to_ver_1=ver_p_pim1-ver_pip_pim2;
  TVector3 ver_to_ver_2=ver_p_pim2-ver_pip_pim1;

  double oa_pim1_p=R2D*openingangle(pim1->Vect(),p->Vect());
  double oa_pim2_p=R2D*openingangle(pim2->Vect(),p->Vect());
  double oa_pip_p=R2D*openingangle(pip->Vect(),p->Vect());
  double oa_pim1_pim2=R2D*openingangle(pim1->Vect(),pim2->Vect());
  double oa_pim1_pip=R2D*openingangle(pim1->Vect(),pip->Vect());
  double oa_pim2_pip=R2D*openingangle(pim2->Vect(),pip->Vect());
                  
  double oa_lambda_1=R2D*openingangle(ver_to_ver_1,gammappim1->Vect());
  double oa_lambda_2=R2D*openingangle(ver_to_ver_2,gammappim2->Vect());
      
  double dist_p_pim1=trackDistance(s.p_r,s.p_z,v2,s.pim1_r,s.pim1_z,v3);
  double dist_p_pim2=trackDistance(s.p_r,s.p_z,v2,s.pim2_r,s.pim2_z,v5);
  double dist_pip_pim1=trackDistance(s.pip_r,s.pip_z,v4,s.pim1_r,s.pim1_z,v3);
  double dist_pip_pim2=trackDistance(s.pip_r,s.pip_z,v4,s.pim2_r,s.pim2_z,v5);
  double dist_lambda1_pip=trackDistance(s.pip_r,s.pip_z,v4,ver_pip_pim1.Z(),getR(ver_pip_pim1),gammappim1->Vect());
  double dist_lambda2_pip=trackDistance(s.pip_r,s.pip_z,v4,ver_pip_pim2.Z(),getR(ver_pip_pim2),gammappim2->Vect());
  double dist_lambda1_pim2=trackDistance(s.pim2_r,s.pim2_z,v5,ver_pip_pim1.Z(),getR(ver_pip_pim1),gammappim1->Vect());
  double dist_lambda2_pim1=trackDistance(s.pim1_r,s.pim1_z,v3,ver_pip_pim2.Z(),getR(ver_pip_pim2),gammappim2->Vect());
  double dist_ver_to_ver_1=ver_to_ver_1.Mag();
  double dist_ver_to_ver_2=ver_to_ver_2.Mag();

  //  cout << "opening angle = " << oa << endl;

  ACC = 1.;
  EFF = 1.;

  /*
    gammappi->Boost(0., 0., -(beam->Beta()));
    p_delta->Boost(0., 0., -(beam->Beta()));
    pi_delta->Boost(0., 0., -(beam->Beta()));
    p_delta->Boost( -gammappi->Px()/gammappi->E(), -gammappi->Py()/gammappi->E(), -gammappi->Pz()/gammappi->E());
    pi_delta->Boost( -gammappi->Px()/gammappi->E(), -gammappi->Py()/gammappi->E(), -gammappi->Pz()/gammappi->E());
  */
  //cout << "Poczatek obliczen..." << endl;

  //double ang_cut = 0.;
  //double ang_cut = 9.;

  //double close_cut = 9.;
  //double nonfit_close_cut = -4.;
  //double close_cut = 0.;
  //double nonfit_close_cut = 0.;
  //double close_cut = 4.;


#ifdef FLANCH
  //insidePim1S0 = (pPim1S0 == 0) ? 0 : pPim1S0->IsInside(pim1_z,pim1_theta);
  //insidePim1S1 = (pPim1S1 == 0) ? 0 : pPim1S1->IsInside(pim1_z,pim1_theta);
  //insideEpS0 = (pPS0 == 0) ? 0 : pPS0->IsInside(p_z,p_theta);
  //insidePS1 = (pPS1 == 0) ? 0 : pPS1->IsInside(p_z,p_theta);
  //insidePim1S0 = (pPim1S0 == 0) ? 0 : pPim1S0->IsInside(eVert_z,pim1_theta);
  //insidePim1S1 = (pPim1S1 == 0) ? 0 : pPim1S1->IsInside(eVert_z,pim1_theta);
  //insidePS0 = (pPS0 == 0) ? 0 : pPS0->IsInside(eVert_z,p_theta);
  //insideEpS1 = (pPS1 == 0) ? 0 : pPS1->IsInside(eVert_z,p_theta);
#endif

  insideTarget = 1;

#ifdef RECTANG
  //insidePim1S0 = (pim1_theta > 50 && pim1_z < -50 /* && pim1_p<200.*/) ? 1 : 0;
  //insidePim1S1 = (pim1_theta > 50 && pim1_z < -50 /* && pim1_p<200.*/) ? 1 : 0;
  //insidePS0 = (p_theta > 50 && p_z < -50 /* && p_p<200.*/) ? 1 : 0;
  //insidePS1 = (p_theta > 50 && p_z < -50 /* && p_p<200.*/) ? 1 : 0;
#endif

  //#ifdef NOCUT
  //insidePim1S0 = 0;
  //insidePim1S1 = 0;
  //insidePS0 = 0;
  //insidePS1 = 0;
  //#endif


  //NoLeptonP = !((p_oa_lept< close_cut&&p_oa_lept>0.0) &&p_oa_lept>nonfit_close_cut );
  //NoHadronP = !(p_oa_hadr< close_cut &&p_oa_hadr>nonfit_close_cut );
  //NoLeptonPI = !((pim1_oa_lept< close_cut&&pim1_oa_lept>0.0) &&pim1_oa_lept>nonfit_close_cut );
  //NoHadronPI = !(pim1_oa_hadr< close_cut &&pim1_oa_hadr>nonfit_close_cut );
  //NoHadronP = 1;
  //NoHadronPI = 1;

  /*
    NoLeptonP = 1;
    NoHadronP = 1;
    NoLeptonPI = 1;
    NoHadronPI = 1;
  */
  double chi_max=180;
  double sum1=dist_p_pim1*dist_p_pim1+dist_pip_pim2*dist_pip_pim2+dist_lambda1_pip*dist_lambda1_pip+dist_lambda1_pim2*dist_lambda1_pim2;
  double sum2=dist_p_pim2*dist_p_pim2+dist_pip_pim1*dist_pip_pim1+dist_lambda2_pip*dist_lambda2_pip+dist_lambda2_pim1*dist_lambda2_pim1;

  double sum1_1=dist_p_pim1+dist_pip_pim2+dist_lambda1_pip+dist_lambda1_pim2;
  double sum2_1=dist_p_pim2+dist_pip_pim1+dist_lambda2_pip+dist_lambda2_pim1;

      

  //save all important variables
  (*n_out)["isBest"]=s.isBest;
  (*n_out)["isBest_new"]=is_best;
  (*n_out)["event"]=event;
  (*n_out)["hneg_mult"]=hneg_mult;
  (*n_out)["hpos_mult"]=hpos_mult;
  (*n_out)["eVert_x"]=eVert_x;
  (*n_out)["eVert_y"]=eVert_y;
  (*n_out)["eVert_z"]=eVert_z;
  (*n_out)["totalmult"]=totalmult;
  (*n_out)["trigdownscaleflag"]=trigdownscaleflag;
  (*n_out)["trigdownscale"]=trigdownscale;
  (*n_out)["event_mult"]=event_mult;
	  
  (*n_out)["p_p"]=s.p_p;
  (*n_out)["p_theta"] = s.p_theta;
  (*n_out)["p_phi"] = s.p_phi;
  (*n_out)["p_beta"] = s.p_beta_new;
  (*n_out)["p_m"] = p_mass;
  (*n_out)["p_dedx"]=s.p_dedx_mdc;
  (*n_out)["p_q"]=s.p_q;
	  
  (*n_out)["p_sim_p"]=s.p_sim_p;
  (*n_out)["p_sim_id"]=s.p_sim_id;
  (*n_out)["p_sim_parentid"]=s.p_sim_parentid;
  (*n_out)["p_sim_vertex_x"]=s.p_sim_vertexx;
  (*n_out)["p_sim_vertex_y"]=s.p_sim_vertexy;
  (*n_out)["p_sim_vertex_z"]=s.p_sim_vertexz;
	  
  (*n_out)["pip_p"]=s.pip_p;
  (*n_out)["pip_theta"] = s.pip_theta;
  (*n_out)["pip_phi"] = s.pip_phi;
  (*n_out)["pip_beta"] = s.pip_beta_new;
  (*n_out)["pip_m"] = pip_mass;
  (*n_out)["pip_dedx"]=s.pip_dedx_mdc;
  (*n_out)["pip_q"]=s.pip_q;
	  
  (*n_out)["pip_sim_p"]=s.pip_sim_p;
  (*n_out)["pip_sim_id"]=s.pip_sim_id;
  (*n_out)["pip_sim_parentid"]=s.pip_sim_parentid;
  (*n_out)["pip_sim_vertex_x"]=s.pip_sim_vertexx;
  (*n_out)["pip_sim_vertex_y"]=s.pip_sim_vertexy;
  (*n_out)["pip_sim_vertex_z"]=s.pip_sim_vertexz;
	  
	  
  (*n_out)["pim1_p"]=s.pim1_p;
  (*n_out)["pim1_theta"] = s.pim1_theta;
  (*n_out)["pim1_phi"] = s.pim1_phi;
  (*n_out)["pim1_beta"] = s.pim1_beta_new;
  (*n_out)["pim1_m"] = pim1_mass;
  (*n_out)["pim1_dedx"]=s.pim1_dedx_mdc;
  (*n_out)["pim1_q"]=s.pim1_q;
	  
  (*n_out)["pim1_sim_p"]=s.pim1_sim_p;
  (*n_out)["pim1_sim_id"]=s.pim1_sim_id;
  (*n_out)["pim1_sim_parentid"]=s.pim1_sim_parentid;
  (*n_out)["pim1_sim_vertex_x"]=s.pim1_sim_vertexx;
  (*n_out)["pim1_sim_vertex_y"]=s.pim1_sim_vertexy;
  (*n_out)["pim1_sim_vertex_z"]=s.pim1_sim_vertexz;
	  
	  
  (*n_out)["pim2_p"]=s.pim2_p;
  (*n_out)["pim2_theta"] = s.pim2_theta;
  (*n_out)["pim2_phi"] = s.pim2_phi;
  (*n_out)["pim2_beta"] = s.pim2_beta_new;
  (*n_out)["pim2_m"] = pim2_mass;
  (*n_out)["pim2_dedx"]=s.pim2_dedx_mdc;
  (*n_out)["pim2_q"]=s.pim2_q;
	  
  (*n_out)["pim2_sim_p"]=s.pim2_sim_p;
  (*n_out)["pim2_sim_id"]=s.pim2_sim_id;
  (*n_out)["pim2_sim_parentid"]=s.pim2_sim_parentid;
  (*n_out)["pim2_sim_vertex_x"]=s.pim2_sim_vertexx;
  (*n_out)["pim2_sim_vertex_y"]=s.pim2_sim_vertexy;
  (*n_out)["pim2_sim_vertex_z"]=s.pim2_sim_vertexz;
	  	  
  (*n_out)["dist_pip_pim1"]=dist_pip_pim1;
  (*n_out)["dist_pip_pim2"] = dist_pip_pim2;
  (*n_out)["dist_p_pim1"] = dist_p_pim1;
  (*n_out)["dist_p_pim2"] = dist_p_pim2;
  (*n_out)["dist_lambda1_pim2"] = dist_lambda1_pim2;
  (*n_out)["dist_lambda1_pip"] = dist_lambda1_pip;
  (*n_out)["dist_lambda2_pim1"] = dist_lambda2_pim1;
  (*n_out)["dist_lambda2_pip"] = dist_lambda2_pip;
  (*n_out)["dist_ver_to_ver_1"]=dist_ver_to_ver_1;
  (*n_out)["dist_ver_to_ver_2"]=dist_ver_to_ver_2;
	  
  (*n_out)["m_inv_p_pim1"] = m_inv_ppim1;
  (*n_out)["m_inv_p_pim2"] = m_inv_ppim2;
  (*n_out)["m_inv_pip_pim1"] = m_inv_pippim1;
  (*n_out)["m_inv_pip_pim2"] = m_inv_pippim2;
  (*n_out)["m_inv_p_pim_pip_pim"] = m_inv_ppimpippim;

  (*n_out)["ver_p_pim1_x"]=ver_p_pim1.X();
  (*n_out)["ver_p_pim1_y"]=ver_p_pim1.Y();
  (*n_out)["ver_p_pim1_z"]=ver_p_pim1.Z();

  (*n_out)["ver_p_pim2_x"]=ver_p_pim2.X();
  (*n_out)["ver_p_pim2_y"]=ver_p_pim2.Y();
  (*n_out)["ver_p_pim2_z"]=ver_p_pim2.Z();

  (*n_out)["ver_pip_pim1_x"]=ver_pip_pim1.X();
  (*n_out)["ver_pip_pim1_y"]=ver_pip_pim1.Y();
  (*n_out)["ver_pip_pim1_z"]=ver_pip_pim1.Z();

  (*n_out)["ver_pip_pim2_x"]=ver_pip_pim2.X();
  (*n_out)["ver_pip_pim2_y"]=ver_pip_pim2.Y();
  (*n_out)["ver_pip_pim2_z"]=ver_pip_pim2.Z();

  (*n_out)["sum_dist_1"]=sum1_1;
  (*n_out)["sum_dist_2"]=sum2_1;
	  
  (*n_out)["sum_dist2_1"]=sum1;
  (*n_out)["sum_dist2_2"]=sum2;

  (*n_out)["oa_lambda_1"]=oa_lambda_1;
  (*n_out)["oa_lambda_2"]=oa_lambda_2;
  (*n_out)["oa_pim1_p"]=oa_pim1_p;
  (*n_out)["oa_pim2_p"]=oa_pim2_p;
  (*n_out)["oa_pip_p"]=oa_pip_p;
  (*n_out)["oa_pim1_pim2"]=oa_pim1_pim2;
  (*n_out)["oa_pim1_pip"]=oa_pim1_pip;
  (*n_out)["oa_pim2_pip"]=oa_pim2_pip;
	 
  (*n_out)["miss_mass_kp"]=miss->M();
	  	  
  n_out->fill();
}

Int_t PPimPipPim::GetEntry(Long64_t entry)
{
  // Read contents of entry.
  if (!fChain) return 0;
  return fChain->GetEntry(entry);
}
Long64_t PPimPipPim::LoadTree(Long64_t entry)
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

void PPimPipPim::Init(TTree *tree)
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
  
  fChain->SetBranchAddress("eVert_chi2", &eVert_chi2, &b_eVert_chi2);
  fChain->SetBranchAddress("eVert_x", &eVert_x, &b_eVert_x);
  fChain->SetBranchAddress("eVert_y", &eVert_y, &b_eVert_y);
  fChain->SetBranchAddress("eVert_z", &eVert_z, &b_eVert_z);
  fChain->SetBranchAddress("event", &event, &b_event);
  fChain->SetBranchAddress("hneg_mult", &hneg_mult, &b_hneg_mult);
  fChain->SetBranchAddress("hpos_mult", &hpos_mult, &b_hpos_mult);
  fChain->SetBranchAddress("isBest", &isBest, &b_isBest);
  fChain->SetBranchAddress("lneg_mult", &lneg_mult, &b_lneg_mult);
  fChain->SetBranchAddress("lpos_mult", &lpos_mult, &b_lpos_mult);
  fChain->SetBranchAddress("p_beta", &p_beta, &b_p_beta);
  fChain->SetBranchAddress("p_beta_new", &p_beta_new, &b_p_beta_new);
  fChain->SetBranchAddress("p_dedx_in", &p_dedx_in, &b_p_dedx_in);
  fChain->SetBranchAddress("p_dedx_in_sigma", &p_dedx_in_sigma, &b_p_dedx_in_sigma);
  fChain->SetBranchAddress("p_dedx_mdc", &p_dedx_mdc, &b_p_dedx_mdc);
  fChain->SetBranchAddress("p_dedx_mdc_sigma", &p_dedx_mdc_sigma, &b_p_dedx_mdc_sigma);
  fChain->SetBranchAddress("p_dedx_out", &p_dedx_out, &b_p_dedx_out);
  fChain->SetBranchAddress("p_dedx_out_sigma", &p_dedx_out_sigma, &b_p_dedx_out_sigma);
  fChain->SetBranchAddress("p_dedx_tof", &p_dedx_tof, &b_p_dedx_tof);
  fChain->SetBranchAddress("p_id", &p_id, &b_p_id);
  fChain->SetBranchAddress("p_isring", &p_isring, &b_p_isring);
  fChain->SetBranchAddress("p_kIsLepton", &p_kIsLepton, &b_p_kIsLepton);
  fChain->SetBranchAddress("p_kIsUsed", &p_kIsUsed, &b_p_kIsUsed);
  fChain->SetBranchAddress("p_mdcchi2", &p_mdcchi2, &b_p_mdcchi2);
  fChain->SetBranchAddress("p_oa_hadr", &p_oa_hadr, &b_p_oa_hadr);
  fChain->SetBranchAddress("p_oa_lept", &p_oa_lept, &b_p_oa_lept);
  fChain->SetBranchAddress("p_p", &p_p, &b_p_p);
  fChain->SetBranchAddress("p_phi", &p_phi, &b_p_phi);
  fChain->SetBranchAddress("p_q", &p_q, &b_p_q);
  fChain->SetBranchAddress("p_r", &p_r, &b_p_r);
  fChain->SetBranchAddress("p_resolution", &p_resolution, &b_p_resolution);
  fChain->SetBranchAddress("p_rkchi2", &p_rkchi2, &b_p_rkchi2);
  fChain->SetBranchAddress("p_sector", &p_sector, &b_p_sector);
  fChain->SetBranchAddress("p_shw_sum0", &p_shw_sum0, &b_p_shw_sum0);
  fChain->SetBranchAddress("p_shw_sum1", &p_shw_sum1, &b_p_shw_sum1);
  fChain->SetBranchAddress("p_shw_sum2", &p_shw_sum2, &b_p_shw_sum2);
  fChain->SetBranchAddress("p_sim_corrflag", &p_sim_corrflag, &b_p_sim_corrflag);
  fChain->SetBranchAddress("p_sim_geninfo", &p_sim_geninfo, &b_p_sim_geninfo);
  fChain->SetBranchAddress("p_sim_geninfo1", &p_sim_geninfo1, &b_p_sim_geninfo1);
  fChain->SetBranchAddress("p_sim_geninfo2", &p_sim_geninfo2, &b_p_sim_geninfo2);
  fChain->SetBranchAddress("p_sim_genweight", &p_sim_genweight, &b_p_sim_genweight);
  fChain->SetBranchAddress("p_sim_id", &p_sim_id, &b_p_sim_id);
  fChain->SetBranchAddress("p_sim_iscommon", &p_sim_iscommon, &b_p_sim_iscommon);
  fChain->SetBranchAddress("p_sim_mediumid", &p_sim_mediumid, &b_p_sim_mediumid);
  fChain->SetBranchAddress("p_sim_p", &p_sim_p, &b_p_sim_p);
  fChain->SetBranchAddress("p_sim_parentid", &p_sim_parentid, &b_p_sim_parentid);
  fChain->SetBranchAddress("p_sim_primaryflag", &p_sim_primaryflag, &b_p_sim_primaryflag);
  fChain->SetBranchAddress("p_sim_processid", &p_sim_processid, &b_p_sim_processid);
  fChain->SetBranchAddress("p_sim_px", &p_sim_px, &b_p_sim_px);
  fChain->SetBranchAddress("p_sim_py", &p_sim_py, &b_p_sim_py);
  fChain->SetBranchAddress("p_sim_pz", &p_sim_pz, &b_p_sim_pz);
  fChain->SetBranchAddress("p_sim_vertexx", &p_sim_vertexx, &b_p_sim_vertexx);
  fChain->SetBranchAddress("p_sim_vertexy", &p_sim_vertexy, &b_p_sim_vertexy);
  fChain->SetBranchAddress("p_sim_vertexz", &p_sim_vertexz, &b_p_sim_vertexz);
  fChain->SetBranchAddress("p_system", &p_system, &b_p_system);
  fChain->SetBranchAddress("p_theta", &p_theta, &b_p_theta);
  fChain->SetBranchAddress("p_tof_exp", &p_tof_exp, &b_p_tof_exp);
  fChain->SetBranchAddress("p_tof_mom", &p_tof_mom, &b_p_tof_mom);
  fChain->SetBranchAddress("p_tof_new", &p_tof_new, &b_p_tof_new);
  fChain->SetBranchAddress("p_tofino_mult", &p_tofino_mult, &b_p_tofino_mult);
  fChain->SetBranchAddress("p_track_length", &p_track_length, &b_p_track_length);
  fChain->SetBranchAddress("p_z", &p_z, &b_p_z);
  fChain->SetBranchAddress("pim1_beta", &pim1_beta, &b_pim1_beta);
  fChain->SetBranchAddress("pim1_beta_new", &pim1_beta_new, &b_pim1_beta_new);
  fChain->SetBranchAddress("pim1_dedx_in", &pim1_dedx_in, &b_pim1_dedx_in);
  fChain->SetBranchAddress("pim1_dedx_in_sigma", &pim1_dedx_in_sigma, &b_pim1_dedx_in_sigma);
  fChain->SetBranchAddress("pim1_dedx_mdc", &pim1_dedx_mdc, &b_pim1_dedx_mdc);
  fChain->SetBranchAddress("pim1_dedx_mdc_sigma", &pim1_dedx_mdc_sigma, &b_pim1_dedx_mdc_sigma);
  fChain->SetBranchAddress("pim1_dedx_out", &pim1_dedx_out, &b_pim1_dedx_out);
  fChain->SetBranchAddress("pim1_dedx_out_sigma", &pim1_dedx_out_sigma, &b_pim1_dedx_out_sigma);
  fChain->SetBranchAddress("pim1_dedx_tof", &pim1_dedx_tof, &b_pim1_dedx_tof);
  fChain->SetBranchAddress("pim1_id", &pim1_id, &b_pim1_id);
  fChain->SetBranchAddress("pim1_isring", &pim1_isring, &b_pim1_isring);
  fChain->SetBranchAddress("pim1_kIsLepton", &pim1_kIsLepton, &b_pim1_kIsLepton);
  fChain->SetBranchAddress("pim1_kIsUsed", &pim1_kIsUsed, &b_pim1_kIsUsed);
  fChain->SetBranchAddress("pim1_mdcchi2", &pim1_mdcchi2, &b_pim1_mdcchi2);
  fChain->SetBranchAddress("pim1_oa_hadr", &pim1_oa_hadr, &b_pim1_oa_hadr);
  fChain->SetBranchAddress("pim1_oa_lept", &pim1_oa_lept, &b_pim1_oa_lept);
  fChain->SetBranchAddress("pim1_p", &pim1_p, &b_pim1_p);
  fChain->SetBranchAddress("pim1_phi", &pim1_phi, &b_pim1_phi);
  fChain->SetBranchAddress("pim1_q", &pim1_q, &b_pim1_q);
  fChain->SetBranchAddress("pim1_r", &pim1_r, &b_pim1_r);
  fChain->SetBranchAddress("pim1_resolution", &pim1_resolution, &b_pim1_resolution);
  fChain->SetBranchAddress("pim1_rkchi2", &pim1_rkchi2, &b_pim1_rkchi2);
  fChain->SetBranchAddress("pim1_sector", &pim1_sector, &b_pim1_sector);
  fChain->SetBranchAddress("pim1_shw_sum0", &pim1_shw_sum0, &b_pim1_shw_sum0);
  fChain->SetBranchAddress("pim1_shw_sum1", &pim1_shw_sum1, &b_pim1_shw_sum1);
  fChain->SetBranchAddress("pim1_shw_sum2", &pim1_shw_sum2, &b_pim1_shw_sum2);
  fChain->SetBranchAddress("pim1_sim_corrflag", &pim1_sim_corrflag, &b_pim1_sim_corrflag);
  fChain->SetBranchAddress("pim1_sim_geninfo", &pim1_sim_geninfo, &b_pim1_sim_geninfo);
  fChain->SetBranchAddress("pim1_sim_geninfo1", &pim1_sim_geninfo1, &b_pim1_sim_geninfo1);
  fChain->SetBranchAddress("pim1_sim_geninfo2", &pim1_sim_geninfo2, &b_pim1_sim_geninfo2);
  fChain->SetBranchAddress("pim1_sim_genweight", &pim1_sim_genweight, &b_pim1_sim_genweight);
  fChain->SetBranchAddress("pim1_sim_id", &pim1_sim_id, &b_pim1_sim_id);
  fChain->SetBranchAddress("pim1_sim_iscommon", &pim1_sim_iscommon, &b_pim1_sim_iscommon);
  fChain->SetBranchAddress("pim1_sim_mediumid", &pim1_sim_mediumid, &b_pim1_sim_mediumid);
  fChain->SetBranchAddress("pim1_sim_p", &pim1_sim_p, &b_pim1_sim_p);
  fChain->SetBranchAddress("pim1_sim_parentid", &pim1_sim_parentid, &b_pim1_sim_parentid);
  fChain->SetBranchAddress("pim1_sim_primaryflag", &pim1_sim_primaryflag, &b_pim1_sim_primaryflag);
  fChain->SetBranchAddress("pim1_sim_processid", &pim1_sim_processid, &b_pim1_sim_processid);
  fChain->SetBranchAddress("pim1_sim_px", &pim1_sim_px, &b_pim1_sim_px);
  fChain->SetBranchAddress("pim1_sim_py", &pim1_sim_py, &b_pim1_sim_py);
  fChain->SetBranchAddress("pim1_sim_pz", &pim1_sim_pz, &b_pim1_sim_pz);
  fChain->SetBranchAddress("pim1_sim_vertexx", &pim1_sim_vertexx, &b_pim1_sim_vertexx);
  fChain->SetBranchAddress("pim1_sim_vertexy", &pim1_sim_vertexy, &b_pim1_sim_vertexy);
  fChain->SetBranchAddress("pim1_sim_vertexz", &pim1_sim_vertexz, &b_pim1_sim_vertexz);
  fChain->SetBranchAddress("pim1_system", &pim1_system, &b_pim1_system);
  fChain->SetBranchAddress("pim1_theta", &pim1_theta, &b_pim1_theta);
  fChain->SetBranchAddress("pim1_tof_exp", &pim1_tof_exp, &b_pim1_tof_exp);
  fChain->SetBranchAddress("pim1_tof_mom", &pim1_tof_mom, &b_pim1_tof_mom);
  fChain->SetBranchAddress("pim1_tof_new", &pim1_tof_new, &b_pim1_tof_new);
  fChain->SetBranchAddress("pim1_tofino_mult", &pim1_tofino_mult, &b_pim1_tofino_mult);
  fChain->SetBranchAddress("pim1_track_length", &pim1_track_length, &b_pim1_track_length);
  fChain->SetBranchAddress("pim1_z", &pim1_z, &b_pim1_z);
  fChain->SetBranchAddress("pim2_beta", &pim2_beta, &b_pim2_beta);
  fChain->SetBranchAddress("pim2_beta_new", &pim2_beta_new, &b_pim2_beta_new);
  fChain->SetBranchAddress("pim2_dedx_in", &pim2_dedx_in, &b_pim2_dedx_in);
  fChain->SetBranchAddress("pim2_dedx_in_sigma", &pim2_dedx_in_sigma, &b_pim2_dedx_in_sigma);
  fChain->SetBranchAddress("pim2_dedx_mdc", &pim2_dedx_mdc, &b_pim2_dedx_mdc);
  fChain->SetBranchAddress("pim2_dedx_mdc_sigma", &pim2_dedx_mdc_sigma, &b_pim2_dedx_mdc_sigma);
  fChain->SetBranchAddress("pim2_dedx_out", &pim2_dedx_out, &b_pim2_dedx_out);
  fChain->SetBranchAddress("pim2_dedx_out_sigma", &pim2_dedx_out_sigma, &b_pim2_dedx_out_sigma);
  fChain->SetBranchAddress("pim2_dedx_tof", &pim2_dedx_tof, &b_pim2_dedx_tof);
  fChain->SetBranchAddress("pim2_id", &pim2_id, &b_pim2_id);
  fChain->SetBranchAddress("pim2_isring", &pim2_isring, &b_pim2_isring);
  fChain->SetBranchAddress("pim2_kIsLepton", &pim2_kIsLepton, &b_pim2_kIsLepton);
  fChain->SetBranchAddress("pim2_kIsUsed", &pim2_kIsUsed, &b_pim2_kIsUsed);
  fChain->SetBranchAddress("pim2_mdcchi2", &pim2_mdcchi2, &b_pim2_mdcchi2);
  fChain->SetBranchAddress("pim2_oa_hadr", &pim2_oa_hadr, &b_pim2_oa_hadr);
  fChain->SetBranchAddress("pim2_oa_lept", &pim2_oa_lept, &b_pim2_oa_lept);
  fChain->SetBranchAddress("pim2_p", &pim2_p, &b_pim2_p);
  fChain->SetBranchAddress("pim2_phi", &pim2_phi, &b_pim2_phi);
  fChain->SetBranchAddress("pim2_q", &pim2_q, &b_pim2_q);
  fChain->SetBranchAddress("pim2_r", &pim2_r, &b_pim2_r);
  fChain->SetBranchAddress("pim2_resolution", &pim2_resolution, &b_pim2_resolution);
  fChain->SetBranchAddress("pim2_rkchi2", &pim2_rkchi2, &b_pim2_rkchi2);
  fChain->SetBranchAddress("pim2_sector", &pim2_sector, &b_pim2_sector);
  fChain->SetBranchAddress("pim2_shw_sum0", &pim2_shw_sum0, &b_pim2_shw_sum0);
  fChain->SetBranchAddress("pim2_shw_sum1", &pim2_shw_sum1, &b_pim2_shw_sum1);
  fChain->SetBranchAddress("pim2_shw_sum2", &pim2_shw_sum2, &b_pim2_shw_sum2);
  fChain->SetBranchAddress("pim2_sim_corrflag", &pim2_sim_corrflag, &b_pim2_sim_corrflag);
  fChain->SetBranchAddress("pim2_sim_geninfo", &pim2_sim_geninfo, &b_pim2_sim_geninfo);
  fChain->SetBranchAddress("pim2_sim_geninfo1", &pim2_sim_geninfo1, &b_pim2_sim_geninfo1);
  fChain->SetBranchAddress("pim2_sim_geninfo2", &pim2_sim_geninfo2, &b_pim2_sim_geninfo2);
  fChain->SetBranchAddress("pim2_sim_genweight", &pim2_sim_genweight, &b_pim2_sim_genweight);
  fChain->SetBranchAddress("pim2_sim_id", &pim2_sim_id, &b_pim2_sim_id);
  fChain->SetBranchAddress("pim2_sim_iscommon", &pim2_sim_iscommon, &b_pim2_sim_iscommon);
  fChain->SetBranchAddress("pim2_sim_mediumid", &pim2_sim_mediumid, &b_pim2_sim_mediumid);
  fChain->SetBranchAddress("pim2_sim_p", &pim2_sim_p, &b_pim2_sim_p);
  fChain->SetBranchAddress("pim2_sim_parentid", &pim2_sim_parentid, &b_pim2_sim_parentid);
  fChain->SetBranchAddress("pim2_sim_primaryflag", &pim2_sim_primaryflag, &b_pim2_sim_primaryflag);
  fChain->SetBranchAddress("pim2_sim_processid", &pim2_sim_processid, &b_pim2_sim_processid);
  fChain->SetBranchAddress("pim2_sim_px", &pim2_sim_px, &b_pim2_sim_px);
  fChain->SetBranchAddress("pim2_sim_py", &pim2_sim_py, &b_pim2_sim_py);
  fChain->SetBranchAddress("pim2_sim_pz", &pim2_sim_pz, &b_pim2_sim_pz);
  fChain->SetBranchAddress("pim2_sim_vertexx", &pim2_sim_vertexx, &b_pim2_sim_vertexx);
  fChain->SetBranchAddress("pim2_sim_vertexy", &pim2_sim_vertexy, &b_pim2_sim_vertexy);
  fChain->SetBranchAddress("pim2_sim_vertexz", &pim2_sim_vertexz, &b_pim2_sim_vertexz);
  fChain->SetBranchAddress("pim2_system", &pim2_system, &b_pim2_system);
  fChain->SetBranchAddress("pim2_theta", &pim2_theta, &b_pim2_theta);
  fChain->SetBranchAddress("pim2_tof_exp", &pim2_tof_exp, &b_pim2_tof_exp);
  fChain->SetBranchAddress("pim2_tof_mom", &pim2_tof_mom, &b_pim2_tof_mom);
  fChain->SetBranchAddress("pim2_tof_new", &pim2_tof_new, &b_pim2_tof_new);
  fChain->SetBranchAddress("pim2_tofino_mult", &pim2_tofino_mult, &b_pim2_tofino_mult);
  fChain->SetBranchAddress("pim2_track_length", &pim2_track_length, &b_pim2_track_length);
  fChain->SetBranchAddress("pim2_z", &pim2_z, &b_pim2_z);
  fChain->SetBranchAddress("pip_beta", &pip_beta, &b_pip_beta);
  fChain->SetBranchAddress("pip_beta_new", &pip_beta_new, &b_pip_beta_new);
  fChain->SetBranchAddress("pip_dedx_in", &pip_dedx_in, &b_pip_dedx_in);
  fChain->SetBranchAddress("pip_dedx_in_sigma", &pip_dedx_in_sigma, &b_pip_dedx_in_sigma);
  fChain->SetBranchAddress("pip_dedx_mdc", &pip_dedx_mdc, &b_pip_dedx_mdc);
  fChain->SetBranchAddress("pip_dedx_mdc_sigma", &pip_dedx_mdc_sigma, &b_pip_dedx_mdc_sigma);
  fChain->SetBranchAddress("pip_dedx_out", &pip_dedx_out, &b_pip_dedx_out);
  fChain->SetBranchAddress("pip_dedx_out_sigma", &pip_dedx_out_sigma, &b_pip_dedx_out_sigma);
  fChain->SetBranchAddress("pip_dedx_tof", &pip_dedx_tof, &b_pip_dedx_tof);
  fChain->SetBranchAddress("pip_id", &pip_id, &b_pip_id);
  fChain->SetBranchAddress("pip_isring", &pip_isring, &b_pip_isring);
  fChain->SetBranchAddress("pip_kIsLepton", &pip_kIsLepton, &b_pip_kIsLepton);
  fChain->SetBranchAddress("pip_kIsUsed", &pip_kIsUsed, &b_pip_kIsUsed);
  fChain->SetBranchAddress("pip_mdcchi2", &pip_mdcchi2, &b_pip_mdcchi2);
  fChain->SetBranchAddress("pip_oa_hadr", &pip_oa_hadr, &b_pip_oa_hadr);
  fChain->SetBranchAddress("pip_oa_lept", &pip_oa_lept, &b_pip_oa_lept);
  fChain->SetBranchAddress("pip_p", &pip_p, &b_pip_p);
  fChain->SetBranchAddress("pip_phi", &pip_phi, &b_pip_phi);
  fChain->SetBranchAddress("pip_q", &pip_q, &b_pip_q);
  fChain->SetBranchAddress("pip_r", &pip_r, &b_pip_r);
  fChain->SetBranchAddress("pip_resolution", &pip_resolution, &b_pip_resolution);
  fChain->SetBranchAddress("pip_rkchi2", &pip_rkchi2, &b_pip_rkchi2);
  fChain->SetBranchAddress("pip_sector", &pip_sector, &b_pip_sector);
  fChain->SetBranchAddress("pip_shw_sum0", &pip_shw_sum0, &b_pip_shw_sum0);
  fChain->SetBranchAddress("pip_shw_sum1", &pip_shw_sum1, &b_pip_shw_sum1);
  fChain->SetBranchAddress("pip_shw_sum2", &pip_shw_sum2, &b_pip_shw_sum2);
  fChain->SetBranchAddress("pip_sim_corrflag", &pip_sim_corrflag, &b_pip_sim_corrflag);
  fChain->SetBranchAddress("pip_sim_geninfo", &pip_sim_geninfo, &b_pip_sim_geninfo);
  fChain->SetBranchAddress("pip_sim_geninfo1", &pip_sim_geninfo1, &b_pip_sim_geninfo1);
  fChain->SetBranchAddress("pip_sim_geninfo2", &pip_sim_geninfo2, &b_pip_sim_geninfo2);
  fChain->SetBranchAddress("pip_sim_genweight", &pip_sim_genweight, &b_pip_sim_genweight);
  fChain->SetBranchAddress("pip_sim_id", &pip_sim_id, &b_pip_sim_id);
  fChain->SetBranchAddress("pip_sim_iscommon", &pip_sim_iscommon, &b_pip_sim_iscommon);
  fChain->SetBranchAddress("pip_sim_mediumid", &pip_sim_mediumid, &b_pip_sim_mediumid);
  fChain->SetBranchAddress("pip_sim_p", &pip_sim_p, &b_pip_sim_p);
  fChain->SetBranchAddress("pip_sim_parentid", &pip_sim_parentid, &b_pip_sim_parentid);
  fChain->SetBranchAddress("pip_sim_primaryflag", &pip_sim_primaryflag, &b_pip_sim_primaryflag);
  fChain->SetBranchAddress("pip_sim_processid", &pip_sim_processid, &b_pip_sim_processid);
  fChain->SetBranchAddress("pip_sim_px", &pip_sim_px, &b_pip_sim_px);
  fChain->SetBranchAddress("pip_sim_py", &pip_sim_py, &b_pip_sim_py);
  fChain->SetBranchAddress("pip_sim_pz", &pip_sim_pz, &b_pip_sim_pz);
  fChain->SetBranchAddress("pip_sim_vertexx", &pip_sim_vertexx, &b_pip_sim_vertexx);
  fChain->SetBranchAddress("pip_sim_vertexy", &pip_sim_vertexy, &b_pip_sim_vertexy);
  fChain->SetBranchAddress("pip_sim_vertexz", &pip_sim_vertexz, &b_pip_sim_vertexz);
  fChain->SetBranchAddress("pip_system", &pip_system, &b_pip_system);
  fChain->SetBranchAddress("pip_theta", &pip_theta, &b_pip_theta);
  fChain->SetBranchAddress("pip_tof_exp", &pip_tof_exp, &b_pip_tof_exp);
  fChain->SetBranchAddress("pip_tof_mom", &pip_tof_mom, &b_pip_tof_mom);
  fChain->SetBranchAddress("pip_tof_new", &pip_tof_new, &b_pip_tof_new);
  fChain->SetBranchAddress("pip_tofino_mult", &pip_tofino_mult, &b_pip_tofino_mult);
  fChain->SetBranchAddress("pip_track_length", &pip_track_length, &b_pip_track_length);
  fChain->SetBranchAddress("pip_z", &pip_z, &b_pip_z);
  fChain->SetBranchAddress("runnumber", &runnumber, &b_runnumber);
  fChain->SetBranchAddress("totalmult", &totalmult, &b_totalmult);
  fChain->SetBranchAddress("trigbit", &trigbit, &b_trigbit);
  fChain->SetBranchAddress("trigdec", &trigdec, &b_trigdec);
  fChain->SetBranchAddress("trigdownscale", &trigdownscale, &b_trigdownscale);
  fChain->SetBranchAddress("trigdownscaleflag", &trigdownscaleflag, &b_trigdownscaleflag);
  Notify();
}

Bool_t PPimPipPim::Notify()
{
  // The Notify() function is called when a new file is opened. This
  // can be either for a new TTree in a TChain or when when a new TTree
  // is started when using PROOF. It is normally not necessary to make changes
  // to the generated code, but the routine can be extended by the
  // user if needed. The return value is currently not used.

  return kTRUE;
}

void PPimPipPim::Show(Long64_t entry)
{
  // Print contents of entry.
  // If entry is not specified, print current entry
  if (!fChain) return;
  fChain->Show(entry);
}
Int_t PPimPipPim::Cut(Long64_t entry)
{
  // This function may be called from Loop.
  // returns  1 if entry is accepted.
  // returns -1 otherwise.
  return 1;
}
