/********************************************************************
* ../build/pc/WallDict.h
* CAUTION: DON'T CHANGE THIS FILE. THIS FILE IS AUTOMATICALLY GENERATED
*          FROM HEADER FILES LISTED IN G__setup_cpp_environmentXXX().
*          CHANGE THOSE HEADER FILES AND REGENERATE THIS FILE.
********************************************************************/
#ifdef __CINT__
#error ../build/pc/WallDict.h/C is only for compilation. Abort cint.
#endif
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#define G__ANSIHEADER
#define G__DICTIONARY
#define G__PRIVATE_GVALUE
#include "cint/G__ci.h"
extern "C" {
extern void G__cpp_setup_tagtableWallDict();
extern void G__cpp_setup_inheritanceWallDict();
extern void G__cpp_setup_typetableWallDict();
extern void G__cpp_setup_memvarWallDict();
extern void G__cpp_setup_globalWallDict();
extern void G__cpp_setup_memfuncWallDict();
extern void G__cpp_setup_funcWallDict();
extern void G__set_cpp_environmentWallDict();
}


#include "TObject.h"
#include "TMemberInspector.h"
#include "./hwalldetector.h"
#include "./hwalltaskset.h"
#include "./hwallunpacker.h"
#include "./hwallraw.h"
#include "./hwallcontfact.h"
#include "./hwallparrootfileio.h"
#include "./hwallparasciifileio.h"
#include "./hwalllookup.h"
#include "./hwallcal.h"
#include "./hwallcalpar.h"
#include "./hwallcalibrater.h"
#include "./hwalldigitizer.h"
#include "./hwallrawsimfilter.h"
#include "./hwallrawsim.h"
#include "./hwallonehit.h"
#include "./hwallhit.h"
#include "./hwallhitsim.h"
#include "./hwallhitf.h"
#include "./hwallonehitf.h"
#include "./hwallrefwinpar.h"
#include "./hwalldigipar.h"
#include "./hwallhitfsim.h"
#include "./hwallgeompar.h"
#include "./hwalltrbunpacker.h"
#include "./walldef.h"
#include <algorithm>
namespace std { }
using namespace std;

#ifndef G__MEMFUNCBODY
#endif

extern G__linked_taginfo G__WallDictLN_TClass;
extern G__linked_taginfo G__WallDictLN_TBuffer;
extern G__linked_taginfo G__WallDictLN_TMemberInspector;
extern G__linked_taginfo G__WallDictLN_TObject;
extern G__linked_taginfo G__WallDictLN_TNamed;
extern G__linked_taginfo G__WallDictLN_basic_fstreamlEcharcOchar_traitslEchargRsPgR;
extern G__linked_taginfo G__WallDictLN_vectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgR;
extern G__linked_taginfo G__WallDictLN_reverse_iteratorlEvectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgRcLcLiteratorgR;
extern G__linked_taginfo G__WallDictLN_TObjArray;
extern G__linked_taginfo G__WallDictLN_TString;
extern G__linked_taginfo G__WallDictLN_HTask;
extern G__linked_taginfo G__WallDictLN_HCategory;
extern G__linked_taginfo G__WallDictLN_HParIo;
extern G__linked_taginfo G__WallDictLN_HDetector;
extern G__linked_taginfo G__WallDictLN_HWallDetector;
extern G__linked_taginfo G__WallDictLN_iteratorlEbidirectional_iterator_tagcOTObjectmUcOlongcOconstsPTObjectmUmUcOconstsPTObjectmUaNgR;
extern G__linked_taginfo G__WallDictLN_HTaskSet;
extern G__linked_taginfo G__WallDictLN_HWallTaskSet;
extern G__linked_taginfo G__WallDictLN_HldUnpack;
extern G__linked_taginfo G__WallDictLN_HLocation;
extern G__linked_taginfo G__WallDictLN_HWallLookup;
extern G__linked_taginfo G__WallDictLN_HWallUnpacker;
extern G__linked_taginfo G__WallDictLN_HLocatedDataObject;
extern G__linked_taginfo G__WallDictLN_HWallRaw;
extern G__linked_taginfo G__WallDictLN_HParSet;
extern G__linked_taginfo G__WallDictLN_HContainer;
extern G__linked_taginfo G__WallDictLN_HContFact;
extern G__linked_taginfo G__WallDictLN_HWallContFact;
extern G__linked_taginfo G__WallDictLN_HDetParIo;
extern G__linked_taginfo G__WallDictLN_HParRootFile;
extern G__linked_taginfo G__WallDictLN_HDetGeomPar;
extern G__linked_taginfo G__WallDictLN_HDetParRootFileIo;
extern G__linked_taginfo G__WallDictLN_HWallCalPar;
extern G__linked_taginfo G__WallDictLN_HWallDigiPar;
extern G__linked_taginfo G__WallDictLN_HWallRefWinPar;
extern G__linked_taginfo G__WallDictLN_HWallParRootFileIo;
extern G__linked_taginfo G__WallDictLN_HParCond;
extern G__linked_taginfo G__WallDictLN_HDetParAsciiFileIo;
extern G__linked_taginfo G__WallDictLN_HWallParAsciiFileIo;
extern G__linked_taginfo G__WallDictLN_HWallLookupChan;
extern G__linked_taginfo G__WallDictLN_HWallLookupSlot;
extern G__linked_taginfo G__WallDictLN_HWallLookupCrate;
extern G__linked_taginfo G__WallDictLN_HWallCal;
extern G__linked_taginfo G__WallDictLN_HWallCalParCell;
extern G__linked_taginfo G__WallDictLN_HReconstructor;
extern G__linked_taginfo G__WallDictLN_HIterator;
extern G__linked_taginfo G__WallDictLN_HWallCalibrater;
extern G__linked_taginfo G__WallDictLN_HWallGeomPar;
extern G__linked_taginfo G__WallDictLN_HWallDigitizer;
extern G__linked_taginfo G__WallDictLN_HFilter;
extern G__linked_taginfo G__WallDictLN_HWallRawSim;
extern G__linked_taginfo G__WallDictLN_HWallRawSimFilter;
extern G__linked_taginfo G__WallDictLN_HWallOneHit;
extern G__linked_taginfo G__WallDictLN_HWallHit;
extern G__linked_taginfo G__WallDictLN_HWallHitSim;
extern G__linked_taginfo G__WallDictLN_HParamList;
extern G__linked_taginfo G__WallDictLN_HSpecGeomPar;
extern G__linked_taginfo G__WallDictLN_HWallHitF;
extern G__linked_taginfo G__WallDictLN_HWallOneHitF;
extern G__linked_taginfo G__WallDictLN_HWallHitFSim;
extern G__linked_taginfo G__WallDictLN_HTrbBaseUnpacker;
extern G__linked_taginfo G__WallDictLN_HWallTrbUnpacker;

/* STUB derived class for protected member access */
