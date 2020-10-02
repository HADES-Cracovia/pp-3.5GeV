// Do NOT change. Changes will be lost next time file is generated

#define R__DICTIONARY_FILENAME TupleCint

/*******************************************************************/
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#define G__DICTIONARY
#include "RConfig.h"
#include "TClass.h"
#include "TDictAttributeMap.h"
#include "TInterpreter.h"
#include "TROOT.h"
#include "TBuffer.h"
#include "TMemberInspector.h"
#include "TInterpreter.h"
#include "TVirtualMutex.h"
#include "TError.h"

#ifndef G__ROOT
#define G__ROOT
#endif

#include "RtypesImp.h"
#include "TIsAProxy.h"
#include "TFileMergeInfo.h"
#include <algorithm>
#include "TCollectionProxyInfo.h"
/*******************************************************************/

#include "TDataMember.h"

// Since CINT ignores the std namespace, we need to do so in this file.
namespace std {} using namespace std;

// Header files passed as explicit arguments
#include "hntuple.h"

// Header files passed via #pragma extra_include

namespace ROOT {
   static void *new_HNtuple(void *p = 0);
   static void *newArray_HNtuple(Long_t size, void *p);
   static void delete_HNtuple(void *p);
   static void deleteArray_HNtuple(void *p);
   static void destruct_HNtuple(void *p);
   static void streamer_HNtuple(TBuffer &buf, void *obj);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::HNtuple*)
   {
      ::HNtuple *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::HNtuple >(0);
      static ::ROOT::TGenericClassInfo 
         instance("HNtuple", ::HNtuple::Class_Version(), "hntuple.h", 17,
                  typeid(::HNtuple), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::HNtuple::Dictionary, isa_proxy, 16,
                  sizeof(::HNtuple) );
      instance.SetNew(&new_HNtuple);
      instance.SetNewArray(&newArray_HNtuple);
      instance.SetDelete(&delete_HNtuple);
      instance.SetDeleteArray(&deleteArray_HNtuple);
      instance.SetDestructor(&destruct_HNtuple);
      instance.SetStreamerFunc(&streamer_HNtuple);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::HNtuple*)
   {
      return GenerateInitInstanceLocal((::HNtuple*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::HNtuple*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

//______________________________________________________________________________
atomic_TClass_ptr HNtuple::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *HNtuple::Class_Name()
{
   return "HNtuple";
}

//______________________________________________________________________________
const char *HNtuple::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::HNtuple*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int HNtuple::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::HNtuple*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *HNtuple::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::HNtuple*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *HNtuple::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::HNtuple*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
void HNtuple::Streamer(TBuffer &R__b)
{
   // Stream an object of class HNtuple.

   TObject::Streamer(R__b);
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_HNtuple(void *p) {
      return  p ? new(p) ::HNtuple : new ::HNtuple;
   }
   static void *newArray_HNtuple(Long_t nElements, void *p) {
      return p ? new(p) ::HNtuple[nElements] : new ::HNtuple[nElements];
   }
   // Wrapper around operator delete
   static void delete_HNtuple(void *p) {
      delete ((::HNtuple*)p);
   }
   static void deleteArray_HNtuple(void *p) {
      delete [] ((::HNtuple*)p);
   }
   static void destruct_HNtuple(void *p) {
      typedef ::HNtuple current_t;
      ((current_t*)p)->~current_t();
   }
   // Wrapper around a custom streamer member function.
   static void streamer_HNtuple(TBuffer &buf, void *obj) {
      ((::HNtuple*)obj)->::HNtuple::Streamer(buf);
   }
} // end of namespace ROOT for class ::HNtuple

namespace {
  void TriggerDictionaryInitialization_TupleCint_Impl() {
    static const char* headers[] = {
"hntuple.h",
0
    };
    static const char* includePaths[] = {
"/cvmfs/hades.gsi.de/install/root-6.12.06/include",
"/lustre/hades/user/knowakow/PP/FAT/PPIMPIPPIM/TMVAeval_DataDriven_root6/HNTUPLE/",
0
    };
    static const char* fwdDeclCode = R"DICTFWDDCLS(
#line 1 "TupleCint dictionary forward declarations' payload"
#pragma clang diagnostic ignored "-Wkeyword-compat"
#pragma clang diagnostic ignored "-Wignored-attributes"
#pragma clang diagnostic ignored "-Wreturn-type-c-linkage"
extern int __Cling_Autoloading_Map;
class __attribute__((annotate("$clingAutoload$hntuple.h")))  HNtuple;
)DICTFWDDCLS";
    static const char* payloadCode = R"DICTPAYLOAD(
#line 1 "TupleCint dictionary payload"

#ifndef G__VECTOR_HAS_CLASS_ITERATOR
  #define G__VECTOR_HAS_CLASS_ITERATOR 1
#endif

#define _BACKWARD_BACKWARD_WARNING_H
#include "hntuple.h"

#undef  _BACKWARD_BACKWARD_WARNING_H
)DICTPAYLOAD";
    static const char* classesHeaders[]={
"HNtuple", payloadCode, "@",
nullptr};

    static bool isInitialized = false;
    if (!isInitialized) {
      TROOT::RegisterModule("TupleCint",
        headers, includePaths, payloadCode, fwdDeclCode,
        TriggerDictionaryInitialization_TupleCint_Impl, {}, classesHeaders);
      isInitialized = true;
    }
  }
  static struct DictInit {
    DictInit() {
      TriggerDictionaryInitialization_TupleCint_Impl();
    }
  } __TheDictionaryInitializer;
}
void TriggerDictionaryInitialization_TupleCint() {
  TriggerDictionaryInitialization_TupleCint_Impl();
}
