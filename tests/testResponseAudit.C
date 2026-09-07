#include "../ZJetResponseAudit.h"
#include <TMemFile.h>
#include <cassert>
#include <iostream>

void testResponseAudit() {
  using namespace ZJetResponseAudit;
  Scores s{.1,.1,.1,.9,.8,.1,.1,.1};
  auto c=categories(s);
  assert(c.at("deepjet_asstored")==1);
  assert(c.at("deepjet")==6);
  assert(c.at("hybrid_deep050")==6);
  assert(c.at("hybrid_pnet030")==6);
  s.qpnet=.3; assert(categories(s).at("hybrid_pnet030")==1);
  s.ucvb=-1.;assert(categories(s).at("hybrid_pnet030")==0);
  assert(block(123,9,3)==block(123,9,3));
  bool occupied[20]={false};for(int i=0;i<1000;++i)occupied[block(i,9,3)]=true;
  for(bool b:occupied) assert(b);
  TMemFile f("unit.root","RECREATE");double x[]={30.,100.};
  auto h=book(&f,1,x);
  fillRu(h,"native",50,6,6,.2,.1,.3,2.,3);
  fillRu(h,"native",50,6,6,.2,.1,.3,-.5,3);
  auto *m=h.moments.at("native_truth");
  assert(std::abs(m->GetBinContent(1,7,31)-1.5)<1.e-12);
  assert(std::abs(m->GetBinContent(1,7,35)-1.5*.02)<1.e-12);
  auto *copy=dynamic_cast<TH3D*>(m->Clone("copy")); copy->SetDirectory(nullptr);
  m->Add(copy); assert(std::abs(m->GetBinContent(1,7,31)-3.)<1.e-12); delete copy;
  std::cout<<"ResponseAudit unit tests passed\n";
}
