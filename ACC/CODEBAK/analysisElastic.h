#ifndef HANALYSISDST_H
#define HANALYSISDST_H
      
#ifdef __CINT__
#define EXIT_SUCCESS  0
#define EXIT_FAILURE  1
#endif

#ifndef __CINT__
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <stdlib.h>
#include <TCanvas.h>
#include <TCutG.h>
#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <TLatex.h>
#include <TMath.h>
#include <TPad.h>
#include <TPaveText.h>
#include <TROOT.h>
//#include <TRFIOFile.h>
#include <TStopwatch.h>
#include <TString.h>
#include <TStyle.h>
#include <TNtuple.h>
#include <TTree.h>
#include <TUnixSystem.h>

#include "hades.h"

#include "haddef.h"
#include "hstartdef.h"
#include "richdef.h"
#include "hmdcdef.h"
#include "tofdef.h"
#include "rpcdef.h"
#include "showerdef.h"
#include "walldef.h"

#include "heventheader.h"
#include "hiterator.h"
#include "hmatrixcategory.h"
#include "hrecevent.h"
#include "hrootsource.h"
#include "htree.h"

//#include "hmdcraw.h"
//#include "htofraw.h"
//#include "hwallraw.h"



#include "hspectrometer.h"
#include "hgeantkine.h"
#include "hgeantmdc.h"
#include "hgeanttof.h"
#include "hgeantshower.h"
#include "hgeantrpc.h"
#include "hgeantwall.h"


#endif

#endif
