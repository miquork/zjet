#undef NDEBUG
#include "../writeCompactOutput.C"
#include <TProfile.h>
#include <TProfile2D.h>
#include <TProfile3D.h>
#include <TGraphErrors.h>
#include <cassert>

// Also callable on actual local smoke/merged outputs. Compare all histogram
// cells, variances, statistics and profile weight sums before and after export.
void checkCompactPair(const char *source,const char *compact) {
  TFile a(source),b(compact);assert(!a.IsZombie()&&!b.IsZombie());
  Long64_t checked=0;
  auto same=[](double x,double y){return x==y||(std::isnan(x)&&std::isnan(y));};
  std::function<void(TDirectory*,TDirectory*,std::string)> compare;
  compare=[&](TDirectory *x,TDirectory *y,std::string path) {
    TIter next(x->GetListOfKeys());
    while(auto *key=dynamic_cast<TKey*>(next())) {
      const std::string name=key->GetName(),full=path+name;
      if(full=="LegacyFlavor/events") {assert(!y->GetKey(name.c_str()));continue;}
      assert(y->GetKey(name.c_str()));
      if(TClass::GetClass(key->GetClassName())->InheritsFrom(TDirectory::Class())) {
        compare(x->GetDirectory(name.c_str()),y->GetDirectory(name.c_str()),full+"/");continue;
      }
      std::unique_ptr<TObject> aa(key->ReadObj()),bb(y->GetKey(name.c_str())->ReadObj());
      assert(aa&&bb&&aa->IsA()==bb->IsA());
      assert(std::string(aa->GetTitle())==bb->GetTitle());
      auto *h=dynamic_cast<TH1*>(aa.get()),*k=dynamic_cast<TH1*>(bb.get());
      if(h) {
        h->SetDirectory(nullptr);k->SetDirectory(nullptr);
        assert(h->GetNcells()==k->GetNcells()&&h->GetEntries()==k->GetEntries());
        assert(h->GetSumw2N()==k->GetSumw2N());
        for(int i=0;i<h->GetNcells();++i) {
          assert(same(h->GetBinContent(i),k->GetBinContent(i)));
          assert(same(h->GetBinError(i),k->GetBinError(i)));
          if(h->GetSumw2N())assert(same(h->GetSumw2()->At(i),k->GetSumw2()->At(i)));
          if(auto *p=dynamic_cast<TProfile*>(h))assert(p->GetBinEntries(i)==static_cast<TProfile*>(k)->GetBinEntries(i));
          if(auto *p=dynamic_cast<TProfile2D*>(h))assert(p->GetBinEntries(i)==static_cast<TProfile2D*>(k)->GetBinEntries(i));
          if(auto *p=dynamic_cast<TProfile3D*>(h))assert(p->GetBinEntries(i)==static_cast<TProfile3D*>(k)->GetBinEntries(i));
        }
        double sx[13]={},sy[13]={};h->GetStats(sx);k->GetStats(sy);
        for(int i=0;i<13;++i)assert(same(sx[i],sy[i]));
        ++checked;
      }
      if(auto *text=dynamic_cast<TObjString*>(aa.get()))assert(text->GetString()==static_cast<TObjString*>(bb.get())->GetString());
    }
  };
  compare(&a,&b,"");
  auto *tree=static_cast<TTree*>(a.Get("LegacyFlavor/events"));assert(tree);
  auto *n=dynamic_cast<TParameter<Long64_t>*>(b.Get("zjet_replay_entries"));assert(n&&n->GetVal()==tree->GetEntries());
  assert(b.Get("zjet_compact_definition")&&!b.Get("LegacyFlavor/events"));
  std::cout<<"Compact round-trip exact for "<<checked<<" histograms/profiles\n";
}

void testCompactOutput(const char *existing="") {
  const std::string dir="tmp/compact_test_"+std::to_string(gSystem->GetPid());
  assert(gSystem->mkdir(dir.c_str(),true)==0);
  const auto source=*existing?std::string(existing):dir+"/source.root";
  const auto compact=dir+"/compact.root";
  if(!*existing) {
    TFile f(source.c_str(),"RECREATE");TObjString("original provenance").Write("metadata");
    f.mkdir("controls")->cd();
    TH1D h("h","control",3,0,3);h.Sumw2();h.Fill(.5,2);h.Fill(.5,-.5);h.Write();
    TProfile3D p("p","weighted profile",2,0,2,2,0,2,2,0,2);
    p.Fill(.5,.5,.5,1.1,2);p.Fill(.5,.5,.5,.9,.5);p.Write();
    TGraphErrors g;g.SetName("graph");g.SetPoint(0,1,2);g.SetPointError(0,.1,.2);g.Write();
    f.mkdir("LegacyFlavor")->cd();TObjString("replay definition").Write("definition");
    TTree t("events","");float score=.4;t.Branch("score",&score);t.Fill();t.Write();
  }
  writeCompactOutput(source.c_str(),compact.c_str());checkCompactPair(source.c_str(),compact.c_str());
  bool caught=false;try{writeCompactOutput(source.c_str(),compact.c_str());}catch(const std::runtime_error&){caught=true;}
  assert(caught); // Existing output is never overwritten by the macro.
}
