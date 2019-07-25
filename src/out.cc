// Do NOT change. Changes will be lost next time file is generated

#define R__DICTIONARY_FILENAME dOdIout

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
#include "../include/Hit_t.hh"
#include "../include/Track_t.hh"

// Header files passed via #pragma extra_include

namespace ROOT {
   static TClass *vectorlETrack_tgR_Dictionary();
   static void vectorlETrack_tgR_TClassManip(TClass*);
   static void *new_vectorlETrack_tgR(void *p = 0);
   static void *newArray_vectorlETrack_tgR(Long_t size, void *p);
   static void delete_vectorlETrack_tgR(void *p);
   static void deleteArray_vectorlETrack_tgR(void *p);
   static void destruct_vectorlETrack_tgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<Track_t>*)
   {
      vector<Track_t> *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<Track_t>));
      static ::ROOT::TGenericClassInfo 
         instance("vector<Track_t>", -2, "vector", 386,
                  typeid(vector<Track_t>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &vectorlETrack_tgR_Dictionary, isa_proxy, 4,
                  sizeof(vector<Track_t>) );
      instance.SetNew(&new_vectorlETrack_tgR);
      instance.SetNewArray(&newArray_vectorlETrack_tgR);
      instance.SetDelete(&delete_vectorlETrack_tgR);
      instance.SetDeleteArray(&deleteArray_vectorlETrack_tgR);
      instance.SetDestructor(&destruct_vectorlETrack_tgR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<Track_t> >()));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const vector<Track_t>*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlETrack_tgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const vector<Track_t>*)0x0)->GetClass();
      vectorlETrack_tgR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlETrack_tgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlETrack_tgR(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<Track_t> : new vector<Track_t>;
   }
   static void *newArray_vectorlETrack_tgR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<Track_t>[nElements] : new vector<Track_t>[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlETrack_tgR(void *p) {
      delete ((vector<Track_t>*)p);
   }
   static void deleteArray_vectorlETrack_tgR(void *p) {
      delete [] ((vector<Track_t>*)p);
   }
   static void destruct_vectorlETrack_tgR(void *p) {
      typedef vector<Track_t> current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class vector<Track_t>

namespace ROOT {
   static TClass *vectorlEHit_tgR_Dictionary();
   static void vectorlEHit_tgR_TClassManip(TClass*);
   static void *new_vectorlEHit_tgR(void *p = 0);
   static void *newArray_vectorlEHit_tgR(Long_t size, void *p);
   static void delete_vectorlEHit_tgR(void *p);
   static void deleteArray_vectorlEHit_tgR(void *p);
   static void destruct_vectorlEHit_tgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<Hit_t>*)
   {
      vector<Hit_t> *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<Hit_t>));
      static ::ROOT::TGenericClassInfo 
         instance("vector<Hit_t>", -2, "vector", 386,
                  typeid(vector<Hit_t>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &vectorlEHit_tgR_Dictionary, isa_proxy, 4,
                  sizeof(vector<Hit_t>) );
      instance.SetNew(&new_vectorlEHit_tgR);
      instance.SetNewArray(&newArray_vectorlEHit_tgR);
      instance.SetDelete(&delete_vectorlEHit_tgR);
      instance.SetDeleteArray(&deleteArray_vectorlEHit_tgR);
      instance.SetDestructor(&destruct_vectorlEHit_tgR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<Hit_t> >()));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const vector<Hit_t>*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlEHit_tgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const vector<Hit_t>*)0x0)->GetClass();
      vectorlEHit_tgR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlEHit_tgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlEHit_tgR(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<Hit_t> : new vector<Hit_t>;
   }
   static void *newArray_vectorlEHit_tgR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<Hit_t>[nElements] : new vector<Hit_t>[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlEHit_tgR(void *p) {
      delete ((vector<Hit_t>*)p);
   }
   static void deleteArray_vectorlEHit_tgR(void *p) {
      delete [] ((vector<Hit_t>*)p);
   }
   static void destruct_vectorlEHit_tgR(void *p) {
      typedef vector<Hit_t> current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class vector<Hit_t>

namespace {
  void TriggerDictionaryInitialization_out_Impl() {
    static const char* headers[] = {
"../include/Hit_t.hh",
"../include/Track_t.hh",
0
    };
    static const char* includePaths[] = {
"/data/sc_soft/root-6.14.06_install/include",
"/data/Projects_physics/LUXE/2019.07.23/Luxe_LuCaS/src/",
0
    };
    static const char* fwdDeclCode = R"DICTFWDDCLS(
#line 1 "out dictionary forward declarations' payload"
#pragma clang diagnostic ignored "-Wkeyword-compat"
#pragma clang diagnostic ignored "-Wignored-attributes"
#pragma clang diagnostic ignored "-Wreturn-type-c-linkage"
extern int __Cling_Autoloading_Map;
)DICTFWDDCLS";
    static const char* payloadCode = R"DICTPAYLOAD(
#line 1 "out dictionary payload"

#ifndef G__VECTOR_HAS_CLASS_ITERATOR
  #define G__VECTOR_HAS_CLASS_ITERATOR 1
#endif

#define _BACKWARD_BACKWARD_WARNING_H
#include "../include/Hit_t.hh"
#include "../include/Track_t.hh"

#undef  _BACKWARD_BACKWARD_WARNING_H
)DICTPAYLOAD";
    static const char* classesHeaders[]={
nullptr};

    static bool isInitialized = false;
    if (!isInitialized) {
      TROOT::RegisterModule("out",
        headers, includePaths, payloadCode, fwdDeclCode,
        TriggerDictionaryInitialization_out_Impl, {}, classesHeaders, /*has no C++ module*/false);
      isInitialized = true;
    }
  }
  static struct DictInit {
    DictInit() {
      TriggerDictionaryInitialization_out_Impl();
    }
  } __TheDictionaryInitializer;
}
void TriggerDictionaryInitialization_out() {
  TriggerDictionaryInitialization_out_Impl();
}
