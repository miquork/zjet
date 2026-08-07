#include "TFile.h"
#include "TObjString.h"

#include <cassert>
#include <fstream>
#include <iterator>
#include <string>

// Store the exact campaign provenance alongside the histograms in a ROOT file.
void embedCampaignMetadata(const char *rootFile, const char *metadataFile) {
  std::ifstream input(metadataFile);
  assert(input.good());
  const std::string contents((std::istreambuf_iterator<char>(input)),
                             std::istreambuf_iterator<char>());
  assert(!contents.empty());

  TFile output(rootFile,"UPDATE");
  assert(!output.IsZombie());
  output.cd();
  TObjString metadata(contents.c_str());
  metadata.Write("zjet_campaign_metadata",TObject::kOverwrite);
  output.Write();
  output.Close();
}
