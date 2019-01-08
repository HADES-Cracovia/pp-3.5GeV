#ifndef HGRAPHCUT_H
#define HGRAPHCUT_H

#include <TROOT.h>
#include <TFile.h>
#include <TH1.h>
#include <list>
#include <functional>
#include <algorithm>
#include <map>
#include <vector>
#include <TVector3.h>
#include <TLorentzVector.h>
#include "hhypcandidate.h"
#include "hcut.h"


// ****************************************************************************
class HGraphCut : public HCut
{
public:

  HGraphCut(const char* cutname, const char* filename=0, const char* opt="read") : HCut(cutname, filename, opt),
										   p_cut(0), d_cut(0), ep_cut(0), em_cut(0), pip_cut(0), pim_cut(0) 
  {
    //p_cut = getCut("pp_cut");
    p_cut = getCut("Mdc_dEdx_P_cut_mod_ChiiV1.root");
    if (p_cut) cout << "p_cut has been read" << endl;
    d_cut = getCut("d_cut");
    if (d_cut) cout << "d_cut has been read" << endl;
    ep_cut = getCut("ep_cut");
    if (ep_cut) cout << "ep_cut has been read" << endl;
    em_cut = getCut("em_cut");
    if (em_cut) cout << "em_cut has been read" << endl;
    pip_cut=getCut("Mdc_dEdx_PiP_cut_PID_mod_ChiiV2");
    //pip_cut = getCut("pip_cut");
    // **************** TEMPORARY 
    //pip_cut = getCut("p_cut"); // reverse cut
    // ****************************
    if (pip_cut) cout << "pip_cut has been read" << endl;
    pim_cut = getCut("Mdc_dEdx_PiP_cut_PID_mod_ChiiV2");
    if (pim_cut) cout << "pim_cut has been read" << endl;
  }
  virtual ~HGraphCut() {}

  Bool_t select(HReconstructor& rec);

private:

  void graphCut(HHypCandidate *pHyp);
  bool graphCut(HParticle *pPart);

  int id;
  double mom, beta, dedx_mdc;
  TCutG *p_cut;
  TCutG *d_cut;
  TCutG *ep_cut;
  TCutG *em_cut;
  TCutG *pip_cut;
  TCutG *pim_cut;


  ClassDef(HGraphCut, 0)
};
// ****************************************************************************

#endif // HGRAPHCUT_H

