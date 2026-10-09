// Do NOT change. Changes will be lost next time file is generated

#define R__DICTIONARY_FILENAME StRefMultCorr_Dict
#define R__NO_DEPRECATION

/*******************************************************************/
#include <cstddef>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cassert>
#define G__DICTIONARY
#include "ROOT/RConfig.hxx"
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

// Header files passed as explicit arguments
#include "StRefMultCorr.h"
#include "BadRun.h"
#include "Param.h"

// Header files passed via #pragma extra_include

// The generated code does not explicitly qualify STL entities
namespace std {} using namespace std;

namespace ROOT {
   static void *new_StRefMultCorr(void *p = nullptr);
   static void *newArray_StRefMultCorr(Long_t size, void *p);
   static void delete_StRefMultCorr(void *p);
   static void deleteArray_StRefMultCorr(void *p);
   static void destruct_StRefMultCorr(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::StRefMultCorr*)
   {
      ::StRefMultCorr *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::StRefMultCorr >(nullptr);
      static ::ROOT::TGenericClassInfo 
         instance("StRefMultCorr", ::StRefMultCorr::Class_Version(), "StRefMultCorr.h", 44,
                  typeid(::StRefMultCorr), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::StRefMultCorr::Dictionary, isa_proxy, 4,
                  sizeof(::StRefMultCorr) );
      instance.SetNew(&new_StRefMultCorr);
      instance.SetNewArray(&newArray_StRefMultCorr);
      instance.SetDelete(&delete_StRefMultCorr);
      instance.SetDeleteArray(&deleteArray_StRefMultCorr);
      instance.SetDestructor(&destruct_StRefMultCorr);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::StRefMultCorr*)
   {
      return GenerateInitInstanceLocal(static_cast<::StRefMultCorr*>(nullptr));
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal(static_cast<const ::StRefMultCorr*>(nullptr)); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

//______________________________________________________________________________
atomic_TClass_ptr StRefMultCorr::fgIsA(nullptr);  // static to hold class pointer

//______________________________________________________________________________
const char *StRefMultCorr::Class_Name()
{
   return "StRefMultCorr";
}

//______________________________________________________________________________
const char *StRefMultCorr::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::StRefMultCorr*)nullptr)->GetImplFileName();
}

//______________________________________________________________________________
int StRefMultCorr::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::StRefMultCorr*)nullptr)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *StRefMultCorr::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::StRefMultCorr*)nullptr)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *StRefMultCorr::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::StRefMultCorr*)nullptr)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
void StRefMultCorr::Streamer(TBuffer &R__b)
{
   // Stream an object of class StRefMultCorr.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(StRefMultCorr::Class(),this);
   } else {
      R__b.WriteClassBuffer(StRefMultCorr::Class(),this);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_StRefMultCorr(void *p) {
      return  p ? new(p) ::StRefMultCorr : new ::StRefMultCorr;
   }
   static void *newArray_StRefMultCorr(Long_t nElements, void *p) {
      return p ? new(p) ::StRefMultCorr[nElements] : new ::StRefMultCorr[nElements];
   }
   // Wrapper around operator delete
   static void delete_StRefMultCorr(void *p) {
      delete (static_cast<::StRefMultCorr*>(p));
   }
   static void deleteArray_StRefMultCorr(void *p) {
      delete [] (static_cast<::StRefMultCorr*>(p));
   }
   static void destruct_StRefMultCorr(void *p) {
      typedef ::StRefMultCorr current_t;
      (static_cast<current_t*>(p))->~current_t();
   }
} // end of namespace ROOT for class ::StRefMultCorr

namespace ROOT {
   // Registration Schema evolution read functions
   int RecordReadRules_StRefMultCorr_Dict() {
      return 0;
   }
   static int _R__UNIQUE_DICT_(ReadRules_StRefMultCorr_Dict) = RecordReadRules_StRefMultCorr_Dict();R__UseDummy(_R__UNIQUE_DICT_(ReadRules_StRefMultCorr_Dict));
} // namespace ROOT
namespace {
  void TriggerDictionaryInitialization_StRefMultCorr_Dict_Impl() {
    static const char* headers[] = {
"StRefMultCorr.h",
"BadRun.h",
"Param.h",
nullptr
    };
    static const char* includePaths[] = {
"/opt/homebrew/Cellar/root/6.38.04_1/include/root",
"/Users/kajihara/研究/ana_isobar/StRefMultCorr/",
nullptr
    };
    static const char* fwdDeclCode = R"DICTFWDDCLS(
#line 1 "StRefMultCorr_Dict dictionary forward declarations' payload"

#pragma diagnostic push
#pragma clang diagnostic ignored "-Wkeyword-compat"
#pragma clang diagnostic ignored "-Wignored-attributes"
#pragma clang diagnostic ignored "-Wreturn-type-c-linkage"
extern int __Cling_AutoLoading_Map;
class  __attribute__((annotate("$clingAutoload$StRefMultCorr.h")))  StRefMultCorr;
#pragma diagnostic pop
)DICTFWDDCLS";
    static const char* payloadCode = R"DICTPAYLOAD(
#line 1 "StRefMultCorr_Dict dictionary payload"

#ifndef ROOT_CINT
  #define ROOT_CINT 1
#endif
#ifndef __ROOT__
  #define __ROOT__ 1
#endif

#define _BACKWARD_BACKWARD_WARNING_H
// Inline headers
#include "StRefMultCorr.h"
#include "BadRun.h"
#include "Param.h"

#undef  _BACKWARD_BACKWARD_WARNING_H
)DICTPAYLOAD";
    static const char* classesHeaders[] = {
"StRefMultCorr", payloadCode, "@",
nullptr
};
    static bool isInitialized = false;
    if (!isInitialized) {
      TROOT::RegisterModule("StRefMultCorr_Dict",
        headers, includePaths, payloadCode, fwdDeclCode,
        TriggerDictionaryInitialization_StRefMultCorr_Dict_Impl, {}, classesHeaders, /*hasCxxModule*/false);
      isInitialized = true;
    }
  }
  static struct DictInit {
    DictInit() {
      TriggerDictionaryInitialization_StRefMultCorr_Dict_Impl();
    }
  } __TheDictionaryInitializer;
}
void TriggerDictionaryInitialization_StRefMultCorr_Dict() {
  TriggerDictionaryInitialization_StRefMultCorr_Dict_Impl();
}
