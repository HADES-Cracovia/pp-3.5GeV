//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Wed Feb 19 17:04:27 2020 by ROOT version 5.34/34
// from TTree ppimepem/ppimepem
// found on file: pp_fullstat_ppimepep.root
//////////////////////////////////////////////////////////

#ifndef TMVAeval_h
#define TMVAeval_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.

// Fixed size dimensions of array or collections stored in the TTree if any.

class TMVAeval {
public :
  TTree          *fChain;   //!pointer to the analyzed TTree or TChain
  Int_t           fCurrent; //!current Tree number in a TChain

  // Declaration of leaf types
  Float_t         dist_ep_em;
  Float_t         dist_ep_pim;
  Float_t         dist_lambda_eVert;
  Float_t         dist_lambda_ver_ep_em;
  Float_t         dist_lambda_em;
  Float_t         dist_lambda_ep;
  Float_t         dist_p_eVert;
  Float_t         dist_p_em;
  Float_t         dist_p_pim;
  Float_t         dist_pim_eVert;
  Float_t         dist_ver_to_ver;
  Float_t         eVert_x;
  Float_t         eVert_y;
  Float_t         eVert_z;
  Float_t         em_beta;
  Float_t         em_dedx;
  Float_t         em_m;
  Float_t         em_p;
  Float_t         em_phi;
  Float_t         em_q;
  Float_t         em_theta;
  Float_t         ep_beta;
  Float_t         ep_dedx;
  Float_t         ep_m;
  Float_t         ep_p;
  Float_t         ep_phi;
  Float_t         ep_q;
  Float_t         ep_theta;
  Float_t         event;
  Float_t         event_mult;
  Float_t         hneg_mult;
  Float_t         hpos_mult;
  Float_t         isBest;
  Float_t         k0_p;
  Float_t         k0_pt;
  Float_t         k0_theta;
  Float_t         k0_w;
  Float_t         lambda_mom_z;
  Float_t         lambda_p;
  Float_t         lambda_pt;
  Float_t         lambda_theta;
  Float_t         lambda_w;
  Float_t         m_inv_ep_em;
  Float_t         m_inv_ep_pim;
  Float_t         m_inv_p_em;
  Float_t         m_inv_p_pim;
  Float_t         m_inv_p_pim_ep_em;
  Float_t         miss_mass_kp;
  Float_t         oa_em_ep;
  Float_t         oa_em_p;
  Float_t         oa_ep_p;
  Float_t         oa_lambda;
  Float_t         oa_pim_em;
  Float_t         oa_pim_ep;
  Float_t         oa_pim_p;
  Float_t         p_beta;
  Float_t         p_dedx;
  Float_t         p_m;
  Float_t         p_p;
  Float_t         p_phi;
  Float_t         p_q;
  Float_t         p_theta;
  Float_t         pim_beta;
  Float_t         pim_dedx;
  Float_t         pim_m;
  Float_t         pim_p;
  Float_t         pim_phi;
  Float_t         pim_q;
  Float_t         pim_theta;
  Float_t         simon_cuts;
  Float_t         totalmult;
  Float_t         trigdownscale;
  Float_t         trigdownscaleflag;
  Float_t         ver_ep_em_x;
  Float_t         ver_ep_em_y;
  Float_t         ver_ep_em_z;
  Float_t         ver_ep_pim_x;
  Float_t         ver_ep_pim_y;
  Float_t         ver_ep_pim_z;
  Float_t         ver_p_em_x;
  Float_t         ver_p_em_y;
  Float_t         ver_p_em_z;
  Float_t         ver_p_pim_x;
  Float_t         ver_p_pim_y;
  Float_t         ver_p_pim_z;

  // List of branches
  TBranch        *b_dist_ep_em;   //!
  TBranch        *b_dist_ep_pim;   //!
  TBranch        *b_dist_lambda_eVert;   //!
  TBranch        *b_dist_lambda_ver_ep_em;//!
  TBranch        *b_dist_lambda_em;   //!
  TBranch        *b_dist_lambda_ep;   //!
  TBranch        *b_dist_p_eVert;   //!
  TBranch        *b_dist_p_em;   //!
  TBranch        *b_dist_p_pim;   //!
  TBranch        *b_dist_pim_eVert;   //!
  TBranch        *b_dist_ver_to_ver;   //!
  TBranch        *b_eVert_x;   //!
  TBranch        *b_eVert_y;   //!
  TBranch        *b_eVert_z;   //!
  TBranch        *b_em_beta;   //!
  TBranch        *b_em_dedx;   //!
  TBranch        *b_em_m;   //!
  TBranch        *b_em_p;   //!
  TBranch        *b_em_phi;   //!
  TBranch        *b_em_q;   //!
  TBranch        *b_em_theta;   //!
  TBranch        *b_ep_beta;   //!
  TBranch        *b_ep_dedx;   //!
  TBranch        *b_ep_m;   //!
  TBranch        *b_ep_p;   //!
  TBranch        *b_ep_phi;   //!
  TBranch        *b_ep_q;   //!
  TBranch        *b_ep_theta;   //!
  TBranch        *b_event;   //!
  TBranch        *b_event_mult;   //!
  TBranch        *b_hneg_mult;   //!
  TBranch        *b_hpos_mult;   //!
  TBranch        *b_isBest;   //!
  TBranch        *b_k0_p;   //!
  TBranch        *b_k0_pt;   //!
  TBranch        *b_k0_theta;   //!
  TBranch        *b_k0_w;   //!
  TBranch        *b_lambda_mom_z;   //!
  TBranch        *b_lambda_p;   //!
  TBranch        *b_lambda_pt;   //!
  TBranch        *b_lambda_theta;   //!
  TBranch        *b_lambda_w;   //!
  TBranch        *b_m_inv_ep_em;   //!
  TBranch        *b_m_inv_ep_pim;   //!
  TBranch        *b_m_inv_p_em;   //!
  TBranch        *b_m_inv_p_pim;   //!
  TBranch        *b_m_inv_p_pim_ep_em;   //!
  TBranch        *b_miss_mass_kp;   //!
  TBranch        *b_oa_em_ep;   //!
  TBranch        *b_oa_em_p;   //!
  TBranch        *b_oa_ep_p;   //!
  TBranch        *b_oa_lambda;   //!
  TBranch        *b_oa_pim_em;   //!
  TBranch        *b_oa_pim_ep;   //!
  TBranch        *b_oa_pim_p;   //!
  TBranch        *b_p_beta;   //!
  TBranch        *b_p_dedx;   //!
  TBranch        *b_p_m;   //!
  TBranch        *b_p_p;   //!
  TBranch        *b_p_phi;   //!
  TBranch        *b_p_q;   //!
  TBranch        *b_p_theta;   //!
  TBranch        *b_pim_beta;   //!
  TBranch        *b_pim_dedx;   //!
  TBranch        *b_pim_m;   //!
  TBranch        *b_pim_p;   //!
  TBranch        *b_pim_phi;   //!
  TBranch        *b_pim_q;   //!
  TBranch        *b_pim_theta;   //!
  TBranch        *b_simon_cuts;   //!
  TBranch        *b_totalmult;   //!
  TBranch        *b_trigdownscale;   //!
  TBranch        *b_trigdownscaleflag;   //!
  TBranch        *b_ver_ep_em_x;   //!
  TBranch        *b_ver_ep_em_y;   //!
  TBranch        *b_ver_ep_em_z;   //!
  TBranch        *b_ver_ep_pim_x;   //!
  TBranch        *b_ver_ep_pim_y;   //!
  TBranch        *b_ver_ep_pim_z;   //!
  TBranch        *b_ver_p_em_x;   //!
  TBranch        *b_ver_p_em_y;   //!
  TBranch        *b_ver_p_em_z;   //!
  TBranch        *b_ver_p_pim_x;   //!
  TBranch        *b_ver_p_pim_y;   //!
  TBranch        *b_ver_p_pim_z;   //!

  TMVAeval(TTree *tree=0);
  virtual ~TMVAeval();
  virtual Int_t    Cut(Long64_t entry);
  virtual Int_t    GetEntry(Long64_t entry);
  virtual Long64_t LoadTree(Long64_t entry);
  virtual void     Init(TTree *tree);
  virtual void     Loop(char*  output);
  virtual Bool_t   Notify();
  virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef TMVAeval_cxx
TMVAeval::TMVAeval(TTree *tree) : fChain(0) 
{
  // if parameter tree is not specified (or zero), connect the file
  // used to generate this class and read the Tree.
  if (tree == 0) {
    TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("pp_fullstat_ppimepep.root");
    if (!f || !f->IsOpen()) {
      f = new TFile("pp_fullstat_ppimepep.root");
    }
    f->GetObject("ppimepem",tree);

  }
  Init(tree);
}

TMVAeval::~TMVAeval()
{
  if (!fChain) return;
  delete fChain->GetCurrentFile();
}

Int_t TMVAeval::GetEntry(Long64_t entry)
{
  // Read contents of entry.
  if (!fChain) return 0;
  return fChain->GetEntry(entry);
}
Long64_t TMVAeval::LoadTree(Long64_t entry)
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

void TMVAeval::Init(TTree *tree)
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

  fChain->SetBranchAddress("dist_ep_em", &dist_ep_em, &b_dist_ep_em);
  fChain->SetBranchAddress("dist_ep_pim", &dist_ep_pim, &b_dist_ep_pim);
  fChain->SetBranchAddress("dist_lambda_eVert", &dist_lambda_eVert, &b_dist_lambda_eVert);
  fChain->SetBranchAddress("dist_lambda_ver_ep_em",&dist_lambda_ver_ep_em, &b_dist_lambda_ver_ep_em);
  fChain->SetBranchAddress("dist_lambda_em", &dist_lambda_em, &b_dist_lambda_em);
  fChain->SetBranchAddress("dist_lambda_ep", &dist_lambda_ep, &b_dist_lambda_ep);
  fChain->SetBranchAddress("dist_p_eVert", &dist_p_eVert, &b_dist_p_eVert);
  fChain->SetBranchAddress("dist_p_em", &dist_p_em, &b_dist_p_em);
  fChain->SetBranchAddress("dist_p_pim", &dist_p_pim, &b_dist_p_pim);
  fChain->SetBranchAddress("dist_pim_eVert", &dist_pim_eVert, &b_dist_pim_eVert);
  fChain->SetBranchAddress("dist_ver_to_ver", &dist_ver_to_ver, &b_dist_ver_to_ver);
  fChain->SetBranchAddress("eVert_x", &eVert_x, &b_eVert_x);
  fChain->SetBranchAddress("eVert_y", &eVert_y, &b_eVert_y);
  fChain->SetBranchAddress("eVert_z", &eVert_z, &b_eVert_z);
  fChain->SetBranchAddress("em_beta", &em_beta, &b_em_beta);
  fChain->SetBranchAddress("em_dedx", &em_dedx, &b_em_dedx);
  fChain->SetBranchAddress("em_m", &em_m, &b_em_m);
  fChain->SetBranchAddress("em_p", &em_p, &b_em_p);
  fChain->SetBranchAddress("em_phi", &em_phi, &b_em_phi);
  fChain->SetBranchAddress("em_q", &em_q, &b_em_q);
  fChain->SetBranchAddress("em_theta", &em_theta, &b_em_theta);
  fChain->SetBranchAddress("ep_beta", &ep_beta, &b_ep_beta);
  fChain->SetBranchAddress("ep_dedx", &ep_dedx, &b_ep_dedx);
  fChain->SetBranchAddress("ep_m", &ep_m, &b_ep_m);
  fChain->SetBranchAddress("ep_p", &ep_p, &b_ep_p);
  fChain->SetBranchAddress("ep_phi", &ep_phi, &b_ep_phi);
  fChain->SetBranchAddress("ep_q", &ep_q, &b_ep_q);
  fChain->SetBranchAddress("ep_theta", &ep_theta, &b_ep_theta);
  fChain->SetBranchAddress("event", &event, &b_event);
  fChain->SetBranchAddress("event_mult", &event_mult, &b_event_mult);
  fChain->SetBranchAddress("hneg_mult", &hneg_mult, &b_hneg_mult);
  fChain->SetBranchAddress("hpos_mult", &hpos_mult, &b_hpos_mult);
  fChain->SetBranchAddress("isBest", &isBest, &b_isBest);
  fChain->SetBranchAddress("k0_p", &k0_p, &b_k0_p);
  fChain->SetBranchAddress("k0_pt", &k0_pt, &b_k0_pt);
  fChain->SetBranchAddress("k0_theta", &k0_theta, &b_k0_theta);
  fChain->SetBranchAddress("k0_w", &k0_w, &b_k0_w);
  fChain->SetBranchAddress("lambda_mom_z", &lambda_mom_z, &b_lambda_mom_z);
  fChain->SetBranchAddress("lambda_p", &lambda_p, &b_lambda_p);
  fChain->SetBranchAddress("lambda_pt", &lambda_pt, &b_lambda_pt);
  fChain->SetBranchAddress("lambda_theta", &lambda_theta, &b_lambda_theta);
  fChain->SetBranchAddress("lambda_w", &lambda_w, &b_lambda_w);
  fChain->SetBranchAddress("m_inv_ep_em", &m_inv_ep_em, &b_m_inv_ep_em);
  fChain->SetBranchAddress("m_inv_ep_pim", &m_inv_ep_pim, &b_m_inv_ep_pim);
  fChain->SetBranchAddress("m_inv_p_em", &m_inv_p_em, &b_m_inv_p_em);
  fChain->SetBranchAddress("m_inv_p_pim", &m_inv_p_pim, &b_m_inv_p_pim);
  fChain->SetBranchAddress("m_inv_p_pim_ep_em", &m_inv_p_pim_ep_em, &b_m_inv_p_pim_ep_em);
  fChain->SetBranchAddress("miss_mass_kp", &miss_mass_kp, &b_miss_mass_kp);
  fChain->SetBranchAddress("oa_em_ep", &oa_em_ep, &b_oa_em_ep);
  fChain->SetBranchAddress("oa_em_p", &oa_em_p, &b_oa_em_p);
  fChain->SetBranchAddress("oa_ep_p", &oa_ep_p, &b_oa_ep_p);
  fChain->SetBranchAddress("oa_lambda", &oa_lambda, &b_oa_lambda);
  fChain->SetBranchAddress("oa_pim_em", &oa_pim_em, &b_oa_pim_em);
  fChain->SetBranchAddress("oa_pim_ep", &oa_pim_ep, &b_oa_pim_ep);
  fChain->SetBranchAddress("oa_pim_p", &oa_pim_p, &b_oa_pim_p);
  fChain->SetBranchAddress("p_beta", &p_beta, &b_p_beta);
  fChain->SetBranchAddress("p_dedx", &p_dedx, &b_p_dedx);
  fChain->SetBranchAddress("p_m", &p_m, &b_p_m);
  fChain->SetBranchAddress("p_p", &p_p, &b_p_p);
  fChain->SetBranchAddress("p_phi", &p_phi, &b_p_phi);
  fChain->SetBranchAddress("p_q", &p_q, &b_p_q);
  fChain->SetBranchAddress("p_theta", &p_theta, &b_p_theta);
  fChain->SetBranchAddress("pim_beta", &pim_beta, &b_pim_beta);
  fChain->SetBranchAddress("pim_dedx", &pim_dedx, &b_pim_dedx);
  fChain->SetBranchAddress("pim_m", &pim_m, &b_pim_m);
  fChain->SetBranchAddress("pim_p", &pim_p, &b_pim_p);
  fChain->SetBranchAddress("pim_phi", &pim_phi, &b_pim_phi);
  fChain->SetBranchAddress("pim_q", &pim_q, &b_pim_q);
  fChain->SetBranchAddress("pim_theta", &pim_theta, &b_pim_theta);
  fChain->SetBranchAddress("simon_cuts", &simon_cuts, &b_simon_cuts);
  fChain->SetBranchAddress("totalmult", &totalmult, &b_totalmult);
  fChain->SetBranchAddress("trigdownscale", &trigdownscale, &b_trigdownscale);
  fChain->SetBranchAddress("trigdownscaleflag", &trigdownscaleflag, &b_trigdownscaleflag);
  fChain->SetBranchAddress("ver_ep_em_x", &ver_ep_em_x, &b_ver_ep_em_x);
  fChain->SetBranchAddress("ver_ep_em_y", &ver_ep_em_y, &b_ver_ep_em_y);
  fChain->SetBranchAddress("ver_ep_em_z", &ver_ep_em_z, &b_ver_ep_em_z);
  fChain->SetBranchAddress("ver_ep_pim_x", &ver_ep_pim_x, &b_ver_ep_pim_x);
  fChain->SetBranchAddress("ver_ep_pim_y", &ver_ep_pim_y, &b_ver_ep_pim_y);
  fChain->SetBranchAddress("ver_ep_pim_z", &ver_ep_pim_z, &b_ver_ep_pim_z);
  fChain->SetBranchAddress("ver_p_em_x", &ver_p_em_x, &b_ver_p_em_x);
  fChain->SetBranchAddress("ver_p_em_y", &ver_p_em_y, &b_ver_p_em_y);
  fChain->SetBranchAddress("ver_p_em_z", &ver_p_em_z, &b_ver_p_em_z);
  fChain->SetBranchAddress("ver_p_pim_x", &ver_p_pim_x, &b_ver_p_pim_x);
  fChain->SetBranchAddress("ver_p_pim_y", &ver_p_pim_y, &b_ver_p_pim_y);
  fChain->SetBranchAddress("ver_p_pim_z", &ver_p_pim_z, &b_ver_p_pim_z);
  Notify();
}

Bool_t TMVAeval::Notify()
{
  // The Notify() function is called when a new file is opened. This
  // can be either for a new TTree in a TChain or when when a new TTree
  // is started when using PROOF. It is normally not necessary to make changes
  // to the generated code, but the routine can be extended by the
  // user if needed. The return value is currently not used.

  return kTRUE;
}

void TMVAeval::Show(Long64_t entry)
{
  // Print contents of entry.
  // If entry is not specified, print current entry
  if (!fChain) return;
  fChain->Show(entry);
}
Int_t TMVAeval::Cut(Long64_t entry)
{
  // This function may be called from Loop.
  // returns  1 if entry is accepted.
  // returns -1 otherwise.
  return 1;
}
#endif // #ifdef TMVAeval_cxx
