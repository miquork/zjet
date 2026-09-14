#undef NDEBUG
#include "../inspectMergedFile.C"
#include <TSystem.h>
#include <cassert>
#include <vector>
#include <functional>
#include <unistd.h>

void testInspectMerged() {
  const std::string dir="tmp/inspect_merge_"+std::to_string(gSystem->GetPid());
  assert(gSystem->mkdir(dir.c_str(),true)==0);
  const std::string path=dir+"/valid.root", metadata=dir+"/metadata.json";
  {
    TFile f(path.c_str(),"RECREATE");
    TObjString("{\"campaign\":\"synthetic\"}").Write("zjet_campaign_metadata");
    TObjString("legacy").Write("zjet_analysis_mode");
    f.mkdir("configInfo")->cd();
    TH1D h("SkimCounter","",3,.5,3.5);h.SetBinContent(2,2);h.Write();
    f.mkdir("LegacyFlavor")->cd();
    TTree t("events","");Bool_t isMC=true;std::vector<float> otherB;
    t.Branch("isMC",&isMC);t.Branch("otherB",&otherB);
    for(int i=0;i<20;++i){otherB.assign(i,float(i)/20);t.Fill();}t.Write();
  }
  inspectMergedFile(path.c_str(),metadata.c_str(),true,2);
  std::ifstream in(metadata);std::string content;std::getline(in,content);
  assert(content=="{\"campaign\":\"synthetic\"}");
  auto fail=[](const std::function<void()> &call) {
    bool caught=false;try{call();}catch(const std::runtime_error &e){caught=true;std::cout<<"Expected: "<<e.what()<<"\n";}
    assert(caught);
  };
  fail([&]{inspectMergedFile(path.c_str(),metadata.c_str(),false,2);});
  fail([&]{inspectMergedFile(path.c_str(),metadata.c_str(),true,3);});
  const auto broken=dir+"/truncated.root";
  assert(gSystem->CopyFile(path.c_str(),broken.c_str())==0);
  assert(truncate(broken.c_str(),200)==0);
  fail([&]{inspectMergedFile(broken.c_str(),metadata.c_str(),true,2);});
  std::cout<<"Merged-file inspection tests passed\n";
}
