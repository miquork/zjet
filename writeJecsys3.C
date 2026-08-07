#include "TClass.h"
#include "TDirectory.h"
#include "TFile.h"
#include "TKey.h"
#include "TNamed.h"
#include "TObjString.h"
#include "TProfile2D.h"
#include "TROOT.h"

#include <cassert>

namespace {

void copyDirectory(TDirectory *source, TDirectory *target) {
  assert(source);
  assert(target);

  TIter next(source->GetListOfKeys());
  while (TKey *key = static_cast<TKey*>(next())) {
    TClass *objectClass = gROOT->GetClass(key->GetClassName());
    assert(objectClass);
    if (objectClass->InheritsFrom(TDirectory::Class())) {
      TDirectory *sourceChild = source->GetDirectory(key->GetName());
      TDirectory *targetChild = target->mkdir(key->GetName());
      copyDirectory(sourceChild,targetChild);
    }
    else {
      TObject *object = key->ReadObj();
      target->cd();
      object->Write(key->GetName());
      delete object;
    }
  }
}

void copySample(TFile *source, TFile *target, const char *sample) {
  TDirectory *sampleDirectory = target->mkdir(sample);
  assert(sampleDirectory);
  for (const char *name : {"l2res","l2res1"}) {
    TDirectory *sourceDirectory = source->GetDirectory(name);
    assert(sourceDirectory);
    TDirectory *targetDirectory = sampleDirectory->mkdir(name);
    copyDirectory(sourceDirectory,targetDirectory);
  }
}

} // namespace

// Create the raw-profile hierarchy expected by jecsys3/L2Res.C and L3Res.C.
// This is not the high-level eta-bin hierarchy consumed by reprocess.C.
void writeJecsys3(
  const char *dataFile="rootfiles/zjet_DATA.root",
  const char *mcFile="rootfiles/zjet_MC.root",
  const char *outputFile="rootfiles/zjet_JMENANO_compat.root") {
  TFile data(dataFile,"READ");
  TFile mc(mcFile,"READ");
  assert(!data.IsZombie());
  assert(!mc.IsZombie());

  for (TFile *source : {&data,&mc}) {
    TProfile2D *response = dynamic_cast<TProfile2D*>(source->Get("l2res/p2m0tc"));
    TProfile2D *residual = dynamic_cast<TProfile2D*>(source->Get("l2res/p2restc"));
    assert(response && response->GetEntries()>0);
    assert(residual && residual->GetEntries()>0);
    assert(source->Get("l2res/p2jes"));
    assert(source->Get("l2res/p2res"));
  }

  TFile output(outputFile,"RECREATE");
  assert(!output.IsZombie());
  copySample(&data,&output,"data");
  copySample(&mc,&output,"mc");
  output.cd();
  TNamed method("zjet_method",
                "all accepted Z-jet pairs; two transverse sidebands with half weight");
  method.Write();
  if (TObjString *metadata =
        dynamic_cast<TObjString*>(data.Get("zjet_campaign_metadata"))) {
    metadata->Write("zjet_campaign_metadata",TObject::kOverwrite);
  }
  output.Write();
  output.Close();
}
