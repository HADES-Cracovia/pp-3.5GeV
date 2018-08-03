#ifndef HPIDPOOL_H
#define HPIDPOOL_H

#include <list>
#include <string>
#include <memory>
#include <TVector3.h>
#include <TLorentzVector.h>
#include "hcommondef.h"
#include "hntuple.h"
#include "houtputfile.h"
#include "hparticlecandidate.h"
#include "hhypcandidate.h"
#include "hpool.h"
#include "hhypdatapool.h"
#include "hhyppool.h"
#include "heventpool.h"

using namespace CommonDefinitions;


struct DelPidCand {
  void operator()(MultiHypPair obj) {
     delete obj.second;
  }
};


class HPidPool : public HPool, public HHypDataPool {

public:

   HPidPool(HOutputFile *ptr = 0);
   virtual ~HPidPool();

   bool add(const char* oldname, const char* name, EParticle p1);
   bool add(const char* oldname, const char* name, EParticle p1, EParticle p2);
   bool add(const char* oldname, const char* name, EParticle p1, EParticle p2, EParticle p3);
   bool add(const char* oldname, const char* name, EParticle p1, EParticle p2, EParticle p3, EParticle p4);
   bool add(const char* oldname, const char* name, const char* p1);
   bool add(const char* oldname, const char* name, const char* p1, const char* p2);
   bool add(const char* oldname, const char* name, const char* p1, const char* p2, const char* p3);
   bool add(const char* oldname, const char* name, const char* p1, const char* p2, const char* p3, const char* p4);
   
   void reset();
   void loop(HHypPool& ref);
   void fill();

private:

   MultiHypPid hypPid; // name oldname link

ClassDef(HPidPool, 0) 
};

#endif // HPIDPOOL_H
