#define T_cxx
#include "T.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

void T::Loop()
{
//   In a ROOT session, you can do:
//      Root > .L T.C
//      Root > T t
//      Root > t.GetEntry(12); // Fill t data members with entry number 12
//      Root > t.Show();       // Show values of entry 12
//      Root > t.Show(16);     // Read and show values of entry 16
//      Root > t.Loop();       // Loop on all entries
//

//     This is the loop skeleton where:
//    jentry is the global entry number in the chain
//    ientry is the entry number in the current Tree
//  Note that the argument to GetEntry must be:
//    jentry for TChain::GetEntry
//    ientry for TTree::GetEntry and TBranch::GetEntry
//
//       To read only selected branches, Insert statements like:
// METHOD1:
//    fChain->SetBranchStatus("*",0);  // disable all branches
//    fChain->SetBranchStatus("branchname",1);  // activate branchname
// METHOD2: replace line
//    fChain->GetEntry(jentry);       //read all branches
//by  b_branchname->GetEntry(ientry); //read only this branch
   if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntriesFast();
   int n_p=0;
   int n_pim=0;
   int n_p_pim=0;
   int n_p_pim_common=0;
   int n_p_pim_lambda=0;
   
   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++)
     {
      Long64_t ientry = LoadTree(jentry);
      if(ientry%100==0)
	cout<<"ientry = "<<ientry<<endl;

      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      
      //if (Cut(ientry) < 0) continue;
      
      //loop over all particles in kine
      for(int i=0; i<kMaxHGeantKine_fData;i++)
	{
	  if(HGeantKine_fData_particleID[i]==9)
	    n_pim++;
	  if(HGeantKine_fData_particleID[i]==14)
	    n_p++;

	for(int j=0; j<kMaxHGeantKine_fData;j++)
	  {
	    if(HGeantKine_fData_particleID[i]==14 && HGeantKine_fData_particleID[j] ==9)
	      {
		n_p_pim++;
		if(HGeantKine_fData_parentTrack[i]==HGeantKine_fData_parentTrack[j])
		  {
		    n_p_pim_common++;
		    //cout<<HGeantKine_fData_parentTrack[i]<<endl;
		    if( HGeantKine_fData_parentTrack[j] ==18)
		      n_p_pim_lambda++;
		  }
	      }
	  }
	}
   }
   cout<<"in file p    : "<<n_p<<"\n";
   cout<<"in file pim  : "<<n_pim<<"\n";
   cout<<"in file pim+p: "<<n_p_pim<<"\n";
   cout<<"pim+p from same source: "<<n_p_pim_common<<"\n";
   cout<<"pim+p from lambda: "<<n_p_pim_lambda<<"\n";
   cout<<"***"<<endl;
   cout<<"events       : "<<ientry<<"\n";
}
