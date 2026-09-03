// Infer flavor-tagging transitions and flavor response residuals from the
// compact FlavorMatrix output written by zjet.C.
//
// Important: data provide reconstructed-tag marginals, not true labels.  The
// inferred data transition matrix is therefore the KL projection closest to
// the MC joint matrix with fixed MC truth marginals and observed data reco-tag
// marginals.  Individual transition scale factors remain model dependent.

#include "FlavorMatrixTools.h"

#include "TDecompSVD.h"
#include "TDirectory.h"
#include "TFile.h"
#include "TGraphErrors.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TH3D.h"
#include "TMatrixD.h"
#include "TMatrixDSym.h"
#include "TObjString.h"
#include "TParameter.h"
#include "TProfile.h"
#include "TProfile3D.h"
#include "TSystem.h"
#include "TVectorD.h"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace FlavorMatrixAnalysis {

struct ProfileSummary {
  bool valid = false;
  double mean = 0.;
  double error = 0.;
  double sumWeights = 0.;
};

struct FitResult {
  TVectorD residual;
  TMatrixDSym covariance;
  TVectorD singularValues;
  int rank = 0;
  int rows = 0;
  double nonzeroCondition = std::numeric_limits<double>::infinity();
  double fullCondition = std::numeric_limits<double>::infinity();
  double chi2 = 0.;
  int zeroErrorFallbacks = 0;
  bool solved = false;

  explicit FitResult(int nFlavor)
      : residual(nFlavor), covariance(nFlavor), singularValues(0) {
    residual = 1.;
    covariance.Zero();
  }

  FitResult(const FitResult &other)
      : residual(other.residual), covariance(other.covariance),
        singularValues(other.singularValues), rank(other.rank),
        rows(other.rows), nonzeroCondition(other.nonzeroCondition),
        fullCondition(other.fullCondition), chi2(other.chi2),
        zeroErrorFallbacks(other.zeroErrorFallbacks), solved(other.solved) {}

  FitResult &operator=(const FitResult &other) {
    if (this==&other) return *this;
    residual.ResizeTo(other.residual);
    residual = other.residual;
    covariance.ResizeTo(other.covariance);
    covariance = other.covariance;
    singularValues.ResizeTo(other.singularValues);
    singularValues = other.singularValues;
    rank = other.rank;
    rows = other.rows;
    nonzeroCondition = other.nonzeroCondition;
    fullCondition = other.fullCondition;
    chi2 = other.chi2;
    zeroErrorFallbacks = other.zeroErrorFallbacks;
    solved = other.solved;
    return *this;
  }
};

struct TruthGroup {
  int outputId;
  const char *name;
  std::vector<int> sourceIds;
};

struct PtRangeResult {
  FitResult responseFit;
  FitResult responseFitRuSlope;
  double ruSlope = std::numeric_limits<double>::quiet_NaN();
  FitResult fsrScaleFit;
  FitResult ueScaleFit;
  FitResult mnuScaleFit;
  FitResult fnuScaleFit;
  FitResult mnuFsrScaleFit;
  std::vector<ProfileSummary> dataTagResponse;
  std::vector<ProfileSummary> mcRawTagResponse;
  std::vector<ProfileSummary> mcCompositionCorrectedTagResponse;
  std::vector<ProfileSummary> rawTagRatio;
  std::vector<ProfileSummary> compositionCorrectedTagRatio;
  std::vector<ProfileSummary> mcFsrFraction;
  std::vector<ProfileSummary> dataFsrFraction;
  std::vector<ProfileSummary> fsrRatio;
  std::vector<ProfileSummary> mcUeFraction;
  std::vector<ProfileSummary> dataUeFraction;
  std::vector<ProfileSummary> ueRatio;
  std::vector<ProfileSummary> mcMnuFraction;
  std::vector<ProfileSummary> dataMnuFraction;
  std::vector<ProfileSummary> mnuRatio;
  std::vector<ProfileSummary> mcFnuFraction;
  std::vector<ProfileSummary> dataFnuFraction;
  std::vector<ProfileSummary> fnuRatio;
  std::vector<ProfileSummary> mcMnuFsrFraction;
  std::vector<ProfileSummary> dataMnuFsrFraction;
  std::vector<ProfileSummary> mnuFsrRatio;

  PtRangeResult(int nReco, int nTruth)
      : responseFit(nTruth), responseFitRuSlope(nTruth),
        fsrScaleFit(nTruth), ueScaleFit(nTruth),
        mnuScaleFit(nTruth), fnuScaleFit(nTruth), mnuFsrScaleFit(nTruth),
        dataTagResponse(nReco), mcRawTagResponse(nReco),
        mcCompositionCorrectedTagResponse(nReco), rawTagRatio(nReco),
        compositionCorrectedTagRatio(nReco), mcFsrFraction(nTruth),
        dataFsrFraction(nTruth), fsrRatio(nTruth), mcUeFraction(nTruth),
        dataUeFraction(nTruth), ueRatio(nTruth), mcMnuFraction(nTruth),
        dataMnuFraction(nTruth), mnuRatio(nTruth),
        mcFnuFraction(nTruth), dataFnuFraction(nTruth), fnuRatio(nTruth),
        mcMnuFsrFraction(nTruth), dataMnuFsrFraction(nTruth),
        mnuFsrRatio(nTruth) {}
};

// The undefined reconstructed tag is intentionally omitted: it is empty for
// valid NanoAOD v15 hybrid scores.  Combining d/u+s and undefined+g gives the
// four physically useful truth groups constrained by the four measured tag
// rows.  The undefined parton category is gluon-like in the tagger controls.
const std::vector<int> recoIds = {1,4,5,6};
const std::vector<TruthGroup> truthGroups = {
  {1,"uds",{1,3}},
  {4,"c",{4}},
  {5,"b",{5}},
  {6,"g",{0,6}},
};
const std::vector<int> allTruthSourceIds = {0,1,3,4,5,6};

const char *recoName(int id) {
  if (id==0) return "undefined";
  if (id==1) return "uds";
  if (id==4) return "c";
  if (id==5) return "b";
  if (id==6) return "g";
  return "unused";
}

const char *truthName(int id) {
  if (id==0) return "unused";
  if (id==1) return "uds";
  if (id==3) return "unused";
  if (id==4) return "c";
  if (id==5) return "b";
  if (id==6) return "g";
  return "unused";
}

void labelRecoAxis(TAxis *axis) {
  if (!axis) return;
  axis->SetBinLabel(1,"undefined");
  axis->SetBinLabel(2,"uds");
  axis->SetBinLabel(3,"unused");
  axis->SetBinLabel(4,"unused");
  axis->SetBinLabel(5,"c");
  axis->SetBinLabel(6,"b");
  axis->SetBinLabel(7,"g");
}

void labelTruthAxis(TAxis *axis) {
  if (!axis) return;
  axis->SetBinLabel(1,"unused");
  axis->SetBinLabel(2,"uds");
  axis->SetBinLabel(3,"unused");
  axis->SetBinLabel(4,"unused");
  axis->SetBinLabel(5,"c");
  axis->SetBinLabel(6,"b");
  axis->SetBinLabel(7,"g");
}

bool selectedPtBin(const TAxis *axis, int bin, double minimumPt,
                   double maximumPt) {
  return axis && bin>=1 && bin<=axis->GetNbins() &&
         axis->GetBinUpEdge(bin)>minimumPt+1.e-9 &&
         axis->GetBinLowEdge(bin)<maximumPt-1.e-9;
}

double integratedCount(const TH3D *histogram, int recoId,
                       const std::vector<int> &truthIds,
                       double minimumPt, double maximumPt) {
  if (!histogram) return 0.;
  double result = 0.;
  const int iy = histogram->GetYaxis()->FindFixBin(recoId);
  for (int ix=1; ix<=histogram->GetNbinsX(); ++ix) {
    if (!selectedPtBin(
          histogram->GetXaxis(),ix,minimumPt,maximumPt)) continue;
    for (int truthId : truthIds) {
      const int iz = histogram->GetZaxis()->FindFixBin(truthId);
      result += histogram->GetBinContent(ix,iy,iz);
    }
  }
  return result;
}

double integratedDataCount(const TH3D *histogram, int recoId,
                           double minimumPt, double maximumPt) {
  if (!histogram) return 0.;
  double result = 0.;
  const int iy = histogram->GetYaxis()->FindFixBin(recoId);
  for (int ix=1; ix<=histogram->GetNbinsX(); ++ix) {
    if (!selectedPtBin(
          histogram->GetXaxis(),ix,minimumPt,maximumPt)) continue;
    for (int iz=1; iz<=histogram->GetNbinsZ(); ++iz)
      result += histogram->GetBinContent(ix,iy,iz);
  }
  return result;
}

ProfileSummary summarizeProfile(const TProfile3D *profile, int recoId,
                                const std::vector<int> &truthIds,
                                double minimumPt, double maximumPt) {
  ProfileSummary result;
  if (!profile) return result;
  const int iy = profile->GetYaxis()->FindFixBin(recoId);
  double weightedSum = 0.;
  double numeratorVariance = 0.;
  for (int ix=1; ix<=profile->GetNbinsX(); ++ix) {
    if (!selectedPtBin(
          profile->GetXaxis(),ix,minimumPt,maximumPt)) continue;
    for (int truthId : truthIds) {
      const int iz = profile->GetZaxis()->FindFixBin(truthId);
      const int globalBin = profile->GetBin(ix,iy,iz);
      const double entries = profile->GetBinEntries(globalBin);
      const double value = profile->GetBinContent(globalBin);
      const double error = profile->GetBinError(globalBin);
      if (!std::isfinite(entries) || !std::isfinite(value) ||
          std::fabs(entries)<1.e-15)
        continue;
      result.sumWeights += entries;
      weightedSum += entries*value;
      if (std::isfinite(error) && error>=0.)
        numeratorVariance += entries*entries*error*error;
    }
  }
  if (std::fabs(result.sumWeights)<1.e-12) return result;
  result.mean = weightedSum/result.sumWeights;
  result.error = std::sqrt(std::max(0.,numeratorVariance)) /
                 std::fabs(result.sumWeights);
  result.valid = std::isfinite(result.mean);
  return result;
}

ProfileSummary summarizeProfile(const TProfile *profile, double minimumPt,
                                double maximumPt) {
  ProfileSummary result;
  if (!profile) return result;
  double weightedSum = 0.;
  double numeratorVariance = 0.;
  for (int bin=1; bin<=profile->GetNbinsX(); ++bin) {
    if (!selectedPtBin(profile->GetXaxis(),bin,minimumPt,maximumPt)) continue;
    const double entries = profile->GetBinEntries(bin);
    const double value = profile->GetBinContent(bin);
    const double error = profile->GetBinError(bin);
    if (!std::isfinite(entries) || !std::isfinite(value) ||
        std::fabs(entries)<1.e-15)
      continue;
    result.sumWeights += entries;
    weightedSum += entries*value;
    if (std::isfinite(error) && error>=0.)
      numeratorVariance += entries*entries*error*error;
  }
  if (std::fabs(result.sumWeights)<1.e-12) return result;
  result.mean = weightedSum/result.sumWeights;
  result.error = std::sqrt(std::max(0.,numeratorVariance))/
                 std::fabs(result.sumWeights);
  result.valid = std::isfinite(result.mean) && std::isfinite(result.error);
  return result;
}

double inclusiveRuSlope(TFile *file, const char *variant,
                        double minimumPt, double maximumPt) {
  if (!file) return std::numeric_limits<double>::quiet_NaN();
  const std::string axis = std::string(variant)=="tc" ? "zmmjet" :
                           (std::string(variant)=="pf" ? "jetpt" : "ptave");
  const std::string directory = "truth_hdm/parallel/"+axis+"/";
  const ProfileSummary product = summarizeProfile(dynamic_cast<TProfile*>(
    file->Get((directory+"mpfu_reco_gen_product").c_str())),
    minimumPt,maximumPt);
  const ProfileSummary square = summarizeProfile(dynamic_cast<TProfile*>(
    file->Get((directory+"mpfu_gen_squared").c_str())),
    minimumPt,maximumPt);
  if (!product.valid || !square.valid || std::fabs(square.mean)<1.e-12)
    return std::numeric_limits<double>::quiet_NaN();
  const double slope = product.mean/square.mean;
  return std::isfinite(slope) && slope>0.05 && slope<2.0
    ? slope : std::numeric_limits<double>::quiet_NaN();
}

ProfileSummary scaledSummary(ProfileSummary input, double scale) {
  if (!input.valid || !std::isfinite(scale)) return ProfileSummary();
  input.mean *= scale;
  input.error *= std::fabs(scale);
  input.valid = std::isfinite(input.mean) && std::isfinite(input.error);
  return input;
}

ProfileSummary sumSummaries(const ProfileSummary &first,
                            const ProfileSummary &second) {
  ProfileSummary result;
  if (!first.valid || !second.valid) return result;
  result.mean = first.mean+second.mean;
  // Old productions do not contain an event-level fnu profile. Its central
  // value is exactly linear in the two component means; use a conservative
  // covariance-free uncertainty until the next production provides p3fnu.
  result.error = std::hypot(first.error,second.error);
  result.sumWeights = first.sumWeights;
  result.valid = std::isfinite(result.mean) && std::isfinite(result.error);
  return result;
}

ProfileSummary ratioSummary(const ProfileSummary &numerator,
                            const ProfileSummary &denominator) {
  ProfileSummary result;
  if (!numerator.valid || !denominator.valid ||
      !std::isfinite(denominator.mean) ||
      std::fabs(denominator.mean)<1.e-12)
    return result;
  result.mean = numerator.mean/denominator.mean;
  result.error = std::hypot(
    numerator.error/denominator.mean,
    numerator.mean*denominator.error/
      (denominator.mean*denominator.mean));
  result.sumWeights = numerator.sumWeights;
  result.valid = std::isfinite(result.mean) && std::isfinite(result.error);
  return result;
}

ProfileSummary readComponent(
  TFile *file, const char *component, const char *variant, int recoId,
  const std::vector<int> &truthIds, double minimumPt, double maximumPt,
  double scale=1.) {
  if (!file) return ProfileSummary();
  const std::string path = std::string("FlavorMatrix/p3")+component+
    variant+"_parallel_flavormatrix";
  TProfile3D *profile = dynamic_cast<TProfile3D*>(file->Get(path.c_str()));
  return scaledSummary(summarizeProfile(
    profile,recoId,truthIds,minimumPt,maximumPt),scale);
}

ProfileSummary readComponentAcrossReco(
  TFile *file, const char *component, const char *variant,
  const std::vector<int> &truthIds, double minimumPt, double maximumPt) {
  ProfileSummary result;
  double weightedSum = 0.;
  double numeratorVariance = 0.;
  for (int recoId : recoIds) {
    const ProfileSummary cell = readComponent(
      file,component,variant,recoId,truthIds,minimumPt,maximumPt);
    if (!cell.valid || std::fabs(cell.sumWeights)<1.e-15) continue;
    result.sumWeights += cell.sumWeights;
    weightedSum += cell.sumWeights*cell.mean;
    numeratorVariance += std::pow(cell.sumWeights*cell.error,2);
  }
  if (std::fabs(result.sumWeights)<1.e-12) return result;
  result.mean = weightedSum/result.sumWeights;
  result.error = std::sqrt(std::max(0.,numeratorVariance))/
                 std::fabs(result.sumWeights);
  result.valid = std::isfinite(result.mean) && std::isfinite(result.error);
  return result;
}

ProfileSummary componentHDMFallback(
  TFile *file, int recoId, const std::vector<int> &truthIds,
  double minimumPt, double maximumPt, const char *variant="tc",
  double responseU=ZJetFlavorMatrix::responseU) {
  ProfileSummary result;
  if (!file) return result;
  ProfileSummary component[3];
  component[0] = readComponent(
    file,"m0",variant,recoId,truthIds,minimumPt,maximumPt);
  component[1] = readComponent(
    file,"mn",variant,recoId,truthIds,minimumPt,maximumPt);
  component[2] = readComponent(
    file,"mu",variant,recoId,truthIds,minimumPt,maximumPt);
  for (const ProfileSummary &entry : component)
    if (!entry.valid) return result;
  const double m0 = component[0].mean;
  const double mn = component[1].mean;
  const double mu = component[2].mean;
  const double denominator =
    1.-mn/ZJetFlavorMatrix::responseN-mu/responseU;
  if (std::fabs(denominator)<1.e-9) return result;
  result.mean = (m0-mn-mu)/denominator;
  // This propagation omits component covariances. The merged component means
  // are nevertheless mandatory because an event-wise nonlinear HDM average
  // is both singular and non-mergeable.
  const double numerator = m0-mn-mu;
  const double dM0 = 1./denominator;
  const double dMn =
    (-denominator+numerator/ZJetFlavorMatrix::responseN)/
    (denominator*denominator);
  const double dMu =
    (-denominator+numerator/responseU)/
                     (denominator*denominator);
  result.error = std::sqrt(
    dM0*dM0*component[0].error*component[0].error+
    dMn*dMn*component[1].error*component[1].error+
    dMu*dMu*component[2].error*component[2].error);
  result.sumWeights = component[0].sumWeights;
  result.valid = std::isfinite(result.mean) && std::isfinite(result.error);
  return result;
}

ProfileSummary readHDM(TFile *file, int recoId,
                       const std::vector<int> &truthIds,
                       double minimumPt, double maximumPt,
                       const char *variant="tc",
                       double responseU=ZJetFlavorMatrix::responseU) {
  // HDM is a nonlinear ratio.  Averaging an event-wise HDM profile is not the
  // same as constructing HDM from the merged component means, and individual
  // events can have a nearly singular denominator.  Always use the mergeable
  // m0/mn/mu component profiles here.
  return componentHDMFallback(
    file,recoId,truthIds,minimumPt,maximumPt,variant,responseU);
}

std::string jsonEscape(const std::string &input) {
  std::ostringstream output;
  for (char character : input) {
    if (character=='\\' || character=='\"') output << '\\';
    if (character=='\n') output << "\\n";
    else output << character;
  }
  return output.str();
}

FitResult fitResponse(const TMatrixD &design, const TVectorD &data,
                      const TVectorD &errors, double priorSigma) {
  const int rows = design.GetNrows();
  const int flavors = design.GetNcols();
  FitResult result(flavors);
  result.rows = rows;
  if (flavors==0) return result;

  TMatrixD whitened(rows,flavors);
  TVectorD whitenedData(rows);
  for (int row=0; row!=rows; ++row) {
    double sigma = errors[row];
    if (!std::isfinite(sigma) || sigma<=0.) {
      sigma = 0.05;
      ++result.zeroErrorFallbacks;
    }
    sigma = std::max(sigma,1.e-8);
    whitenedData[row] = data[row]/sigma;
    for (int flavor=0; flavor!=flavors; ++flavor)
      whitened(row,flavor) = design(row,flavor)/sigma;
  }

  if (rows>0) {
    // ROOT's TDecompSVD requires rows >= columns.  Transposition preserves
    // the non-zero singular values for the common underdetermined J<F case.
    TMatrixD diagnosticMatrix = whitened;
    if (rows<flavors) {
      TMatrixD transpose(TMatrixD::kTransposed,whitened);
      diagnosticMatrix.ResizeTo(transpose);
      diagnosticMatrix = transpose;
    }
    TDecompSVD diagnosticSVD(diagnosticMatrix);
    const TVectorD singularValues = diagnosticSVD.GetSig();
    result.singularValues.ResizeTo(singularValues);
    result.singularValues = singularValues;
    double largest = 0.;
    for (int index=0; index!=result.singularValues.GetNrows(); ++index)
      largest = std::max(largest,std::fabs(result.singularValues[index]));
    const double threshold = std::max(1.e-12,largest*1.e-10);
    double smallest = std::numeric_limits<double>::infinity();
    for (int index=0; index!=result.singularValues.GetNrows(); ++index) {
      const double singular = std::fabs(result.singularValues[index]);
      if (singular>threshold) {
        ++result.rank;
        smallest = std::min(smallest,singular);
      }
    }
    if (largest>0. && std::isfinite(smallest))
      result.nonzeroCondition = largest/smallest;
    if (result.rank>=flavors)
      result.fullCondition = result.nonzeroCondition;
  }

  const bool usePrior = std::isfinite(priorSigma) && priorSigma>0.;
  if (!usePrior && rows<flavors) {
    std::cerr << "Response system has fewer measured tag rows than true "
              << "flavors; an unregularized solution is not identifiable. "
              << "Returning unity residuals.\n";
    return result;
  }
  const int augmentedRows = rows+(usePrior ? flavors : 0);
  TMatrixD augmented(augmentedRows,flavors);
  TVectorD target(augmentedRows);
  for (int row=0; row!=rows; ++row) {
    target[row] = whitenedData[row];
    for (int flavor=0; flavor!=flavors; ++flavor)
      augmented(row,flavor) = whitened(row,flavor);
  }
  if (usePrior) {
    for (int flavor=0; flavor!=flavors; ++flavor) {
      augmented(rows+flavor,flavor) = 1./priorSigma;
      target[rows+flavor] = 1./priorSigma;
    }
  }
  if (augmentedRows==0) return result;

  TDecompSVD fitSVD(augmented);
  Bool_t solved = false;
  result.residual = fitSVD.Solve(target,solved);
  result.solved = solved;
  if (!solved) {
    result.residual = 1.;
    return result;
  }

  TMatrixD hessian(flavors,flavors);
  hessian.Zero();
  for (int row=0; row!=augmentedRows; ++row)
    for (int first=0; first!=flavors; ++first)
      for (int second=0; second!=flavors; ++second)
        hessian(first,second) +=
          augmented(row,first)*augmented(row,second);
  TDecompSVD hessianSVD(hessian);
  Bool_t inverted = false;
  TMatrixD inverse = hessianSVD.Invert(inverted);
  if (inverted) {
    for (int first=0; first!=flavors; ++first)
      for (int second=0; second!=flavors; ++second)
        result.covariance(first,second) =
          0.5*(inverse(first,second)+inverse(second,first));
  }

  for (int row=0; row!=rows; ++row) {
    double prediction = 0.;
    for (int flavor=0; flavor!=flavors; ++flavor)
      prediction += design(row,flavor)*result.residual[flavor];
    double sigma = errors[row];
    if (!std::isfinite(sigma) || sigma<=0.) sigma = 0.05;
    result.chi2 += std::pow((data[row]-prediction)/sigma,2);
  }
  return result;
}

} // namespace FlavorMatrixAnalysis

void analyzeFlavorMatrix(
  const char *dataFileName="/tmp/zjet_flavormatrix_DATA.root",
  const char *mcFileName="/tmp/zjet_flavormatrix_MC.root",
  const char *outputDirectory="output/flavorMatrix",
  double minimumPt=30., double responsePriorSigma=0.20) {
  using namespace FlavorMatrixAnalysis;

  std::unique_ptr<TFile> dataFile(TFile::Open(dataFileName,"READ"));
  std::unique_ptr<TFile> mcFile(TFile::Open(mcFileName,"READ"));
  if (!dataFile || dataFile->IsZombie())
    throw std::runtime_error(std::string("Failed to open data file ")+
                             dataFileName);
  if (!mcFile || mcFile->IsZombie())
    throw std::runtime_error(std::string("Failed to open MC file ")+
                             mcFileName);

  TH3D *dataCounts = dynamic_cast<TH3D*>(dataFile->Get(
    "FlavorMatrix/h3counts_parallel_flavormatrix"));
  TH3D *mcCounts = dynamic_cast<TH3D*>(mcFile->Get(
    "FlavorMatrix/h3counts_parallel_flavormatrix"));
  if (!dataCounts || !mcCounts)
    throw std::runtime_error(
      "Missing FlavorMatrix/h3counts_parallel_flavormatrix");

  const int nReco = recoIds.size();
  const int nTruth = truthGroups.size();
  const double maximumPt = 600.;
  TMatrixD mcRaw(nReco,nTruth);
  TVectorD dataRaw(nReco);
  mcRaw.Zero();
  dataRaw.Zero();
  double totalMC = 0.;
  double totalData = 0.;
  for (int reco=0; reco!=nReco; ++reco) {
    const double dataCount =
      integratedDataCount(
        dataCounts,recoIds[reco],minimumPt,maximumPt);
    if (dataCount < -1.e-9)
      throw std::runtime_error(
        "Negative integrated data parallel count: the IPF likelihood "
        "requires a non-negative population");
    dataRaw[reco] = std::max(0.,dataCount);
    totalData += dataRaw[reco];
    for (int truth=0; truth!=nTruth; ++truth) {
      const double mcCount = integratedCount(
        mcCounts,recoIds[reco],truthGroups[truth].sourceIds,
        minimumPt,maximumPt);
      if (mcCount < -1.e-9)
        throw std::runtime_error(
          "Negative integrated MC parallel count: do not clip signed NLO "
          "weights in a likelihood; use an unweighted population matrix or "
          "a signed-weight treatment");
      mcRaw(reco,truth) = std::max(0.,mcCount);
      totalMC += mcRaw(reco,truth);
    }
  }
  if (!(totalMC>0.) || !(totalData>0.))
    throw std::runtime_error(
      "No positive FlavorMatrix parallel counts above the requested pT");

  for (int truth=0; truth!=nTruth; ++truth) {
    double column = 0.;
    for (int reco=0; reco!=nReco; ++reco) column += mcRaw(reco,truth);
    if (!(column>1.e-12*totalMC))
      throw std::runtime_error(
        "The four-flavor response model has an empty MC truth group " +
        std::string(truthGroups[truth].name));
  }

  TVectorD truthPrior(nTruth);
  TVectorD dataRecoFraction(nReco);
  TMatrixD efficiencyMC(nReco,nTruth);
  TMatrixD jointMC(nReco,nTruth);
  for (int reco=0; reco!=nReco; ++reco)
    dataRecoFraction[reco] = dataRaw[reco]/totalData;
  const double activeTotalMC = totalMC;
  for (int flavor=0; flavor!=nTruth; ++flavor) {
    double column = 0.;
    for (int reco=0; reco!=nReco; ++reco) column += mcRaw(reco,flavor);
    truthPrior[flavor] = column/activeTotalMC;
    for (int reco=0; reco!=nReco; ++reco) {
      efficiencyMC(reco,flavor) = column>0.
        ? mcRaw(reco,flavor)/column : 0.;
      jointMC(reco,flavor) = mcRaw(reco,flavor)/activeTotalMC;
    }
  }

  // Add a very small pseudocount only to avoid structural zeros in a sparse
  // smoke sample.  IPF restores the original MC truth marginals exactly.
  const double pseudocount = std::max(1.e-12,activeTotalMC*1.e-12);
  TMatrixD inferredJoint(nReco,nTruth);
  double initialTotal = 0.;
  for (int reco=0; reco!=nReco; ++reco)
    for (int flavor=0; flavor!=nTruth; ++flavor) {
      inferredJoint(reco,flavor) =
        std::max(mcRaw(reco,flavor),pseudocount);
      initialTotal += inferredJoint(reco,flavor);
    }
  inferredJoint *= 1./initialTotal;

  bool ipfConverged = false;
  int ipfIterations = 0;
  double ipfMaximumError = std::numeric_limits<double>::infinity();
  for (ipfIterations=1; ipfIterations<=10000; ++ipfIterations) {
    for (int reco=0; reco!=nReco; ++reco) {
      double row = 0.;
      for (int flavor=0; flavor!=nTruth; ++flavor)
        row += inferredJoint(reco,flavor);
      if (row>0.) {
        const double scale = dataRecoFraction[reco]/row;
        for (int flavor=0; flavor!=nTruth; ++flavor)
          inferredJoint(reco,flavor) *= scale;
      }
    }
    for (int flavor=0; flavor!=nTruth; ++flavor) {
      double column = 0.;
      for (int reco=0; reco!=nReco; ++reco)
        column += inferredJoint(reco,flavor);
      if (column>0.) {
        const double scale = truthPrior[flavor]/column;
        for (int reco=0; reco!=nReco; ++reco)
          inferredJoint(reco,flavor) *= scale;
      }
    }
    ipfMaximumError = 0.;
    for (int reco=0; reco!=nReco; ++reco) {
      double row = 0.;
      for (int flavor=0; flavor!=nTruth; ++flavor)
        row += inferredJoint(reco,flavor);
      ipfMaximumError = std::max(
        ipfMaximumError,std::fabs(row-dataRecoFraction[reco]));
    }
    for (int flavor=0; flavor!=nTruth; ++flavor) {
      double column = 0.;
      for (int reco=0; reco!=nReco; ++reco)
        column += inferredJoint(reco,flavor);
      ipfMaximumError = std::max(
        ipfMaximumError,std::fabs(column-truthPrior[flavor]));
    }
    if (ipfMaximumError<1.e-10) {
      ipfConverged = true;
      break;
    }
  }
  if (!ipfConverged)
    throw std::runtime_error(
      "Flavor tagging IPF did not converge; refusing to write a seemingly "
      "valid transition matrix");

  TMatrixD efficiencyData(nReco,nTruth);
  TMatrixD transitionSF(nReco,nTruth);
  TMatrixD compositionMC(nReco,nTruth);
  TMatrixD compositionData(nReco,nTruth);
  TMatrixD compositionRatio(nReco,nTruth);
  for (int reco=0; reco!=nReco; ++reco)
    for (int flavor=0; flavor!=nTruth; ++flavor) {
      efficiencyData(reco,flavor) = truthPrior[flavor]>0.
        ? inferredJoint(reco,flavor)/truthPrior[flavor] : 0.;
      transitionSF(reco,flavor) = efficiencyMC(reco,flavor)>0.
        ? efficiencyData(reco,flavor)/efficiencyMC(reco,flavor) : 0.;
      double mcRecoFraction = 0.;
      for (int other=0; other!=nTruth; ++other)
        mcRecoFraction += jointMC(reco,other);
      compositionMC(reco,flavor) = mcRecoFraction>0.
        ? jointMC(reco,flavor)/mcRecoFraction : 0.;
      compositionData(reco,flavor) = dataRecoFraction[reco]>0.
        ? inferredJoint(reco,flavor)/dataRecoFraction[reco] : 0.;
      compositionRatio(reco,flavor) = compositionMC(reco,flavor)>0.
        ? compositionData(reco,flavor)/compositionMC(reco,flavor) : 0.;
    }

  std::vector<ProfileSummary> dataResponse(nReco);
  std::vector<std::vector<ProfileSummary> > mcResponse(
    nReco,std::vector<ProfileSummary>(nTruth));
  int cellFallbacks = 0;
  for (int reco=0; reco!=nReco; ++reco) {
    dataResponse[reco] = readHDM(
      dataFile.get(),recoIds[reco],std::vector<int>{0},
      minimumPt,maximumPt);
    for (int flavor=0; flavor!=nTruth; ++flavor) {
      mcResponse[reco][flavor] = readHDM(
        mcFile.get(),recoIds[reco],truthGroups[flavor].sourceIds,
        minimumPt,maximumPt);
    }
  }

  // A sparse matrix may have an allowed transition with no response-profile
  // entries.  Use the count-weighted true-flavor response as a documented
  // fallback, then the global response only if the whole flavor is empty.
  std::vector<double> truthResponseFallback(nTruth,
    std::numeric_limits<double>::quiet_NaN());
  double globalNumerator = 0.;
  double globalDenominator = 0.;
  for (int flavor=0; flavor!=nTruth; ++flavor) {
    double numerator = 0.;
    double denominator = 0.;
    for (int reco=0; reco!=nReco; ++reco) {
      if (!mcResponse[reco][flavor].valid) continue;
      const double count = mcRaw(reco,flavor);
      numerator += count*mcResponse[reco][flavor].mean;
      denominator += count;
    }
    if (denominator>0.) truthResponseFallback[flavor] = numerator/denominator;
    if (denominator>0.) {
      globalNumerator += numerator;
      globalDenominator += denominator;
    }
  }
  const double globalResponseFallback = globalDenominator>0.
    ? globalNumerator/globalDenominator : 1.;
  for (int reco=0; reco!=nReco; ++reco)
    for (int flavor=0; flavor!=nTruth; ++flavor)
      if (!mcResponse[reco][flavor].valid) {
        mcResponse[reco][flavor].valid = true;
        mcResponse[reco][flavor].mean =
          std::isfinite(truthResponseFallback[flavor])
            ? truthResponseFallback[flavor] : globalResponseFallback;
        mcResponse[reco][flavor].error = 0.;
        ++cellFallbacks;
      }

  std::vector<int> fitReco;
  for (int reco=0; reco!=nReco; ++reco)
    if (dataRecoFraction[reco]>0. && dataResponse[reco].valid)
      fitReco.push_back(reco);
  TMatrixD design(fitReco.size(),nTruth);
  TVectorD observed(fitReco.size());
  TVectorD observedError(fitReco.size());
  for (int row=0; row!=static_cast<int>(fitReco.size()); ++row) {
    const int reco = fitReco[row];
    observed[row] = dataResponse[reco].mean;
    observedError[row] = dataResponse[reco].error;
    for (int flavor=0; flavor!=nTruth; ++flavor)
      design(row,flavor) = compositionData(reco,flavor)*
                           mcResponse[reco][flavor].mean;
  }
  FitResult fit = fitResponse(
    design,observed,observedError,responsePriorSigma);
  if (!fit.solved)
    throw std::runtime_error(
      "Flavor-response SVD solve failed; refusing to write unity residuals "
      "with zero uncertainty");

  TVectorD responseTemplateUncertainty(nTruth);
  TVectorD purityStatUncertainty(nTruth);
  TVectorD taggingEfficiency5PctUncertainty(nTruth);
  TVectorD flavorFraction10PctSensitivity(nTruth);
  responseTemplateUncertainty.Zero();
  purityStatUncertainty.Zero();
  taggingEfficiency5PctUncertainty.Zero();
  flavorFraction10PctSensitivity.Zero();

  for (int row=0; row!=static_cast<int>(fitReco.size()); ++row) {
    const int reco = fitReco[row];
    for (int flavor=0; flavor!=nTruth; ++flavor) {
      const double sigma = mcResponse[reco][flavor].error;
      if (!std::isfinite(sigma) || sigma<=0.) continue;
      TMatrixD plus = design;
      TMatrixD minus = design;
      plus(row,flavor) += compositionData(reco,flavor)*sigma;
      minus(row,flavor) -= compositionData(reco,flavor)*sigma;
      FitResult plusFit = fitResponse(
        plus,observed,observedError,responsePriorSigma);
      FitResult minusFit = fitResponse(
        minus,observed,observedError,responsePriorSigma);
      if (!plusFit.solved || !minusFit.solved) continue;
      for (int resultFlavor=0; resultFlavor!=nTruth; ++resultFlavor) {
        const double shift = 0.5*(plusFit.residual[resultFlavor]-
                                  minusFit.residual[resultFlavor]);
        responseTemplateUncertainty[resultFlavor] += shift*shift;
      }
    }
  }
  for (int flavor=0; flavor!=nTruth; ++flavor)
    responseTemplateUncertainty[flavor] =
      std::sqrt(responseTemplateUncertainty[flavor]);

  for (int row=0; row!=static_cast<int>(fitReco.size()); ++row) {
    const int reco = fitReco[row];
    double rowCount = 0.;
    for (int flavor=0; flavor!=nTruth; ++flavor)
      rowCount += mcRaw(reco,flavor);
    if (!(rowCount>0.)) continue;
    for (int variedFlavor=0; variedFlavor!=nTruth; ++variedFlavor) {
      const double probability = compositionData(reco,variedFlavor);
      const double sigma = std::sqrt(
        std::max(0.,probability*(1.-probability)/rowCount));
      if (!(sigma>0.)) continue;
      TMatrixD variedDesign[2] = {design,design};
      for (int direction=0; direction!=2; ++direction) {
        std::vector<double> varied(nTruth);
        double normalization = 0.;
        for (int flavor=0; flavor!=nTruth; ++flavor) {
          varied[flavor] = compositionData(reco,flavor);
          if (flavor==variedFlavor)
            varied[flavor] = std::max(
              0.,varied[flavor]+(direction==0 ? sigma : -sigma));
          normalization += varied[flavor];
        }
        if (normalization>0.)
          for (int flavor=0; flavor!=nTruth; ++flavor)
            variedDesign[direction](row,flavor) =
              varied[flavor]/normalization*mcResponse[reco][flavor].mean;
      }
      FitResult plusFit = fitResponse(
        variedDesign[0],observed,observedError,responsePriorSigma);
      FitResult minusFit = fitResponse(
        variedDesign[1],observed,observedError,responsePriorSigma);
      if (!plusFit.solved || !minusFit.solved) continue;
      for (int resultFlavor=0; resultFlavor!=nTruth; ++resultFlavor) {
        const double shift = 0.5*(plusFit.residual[resultFlavor]-
                                  minusFit.residual[resultFlavor]);
        purityStatUncertainty[resultFlavor] += shift*shift;
      }
    }
  }
  for (int flavor=0; flavor!=nTruth; ++flavor)
    purityStatUncertainty[flavor] =
      std::sqrt(purityStatUncertainty[flavor]);

  // Preliminary independent 5% relative shifts of each allowed tagging
  // transition.  Renormalizing the affected reconstructed-tag row keeps the
  // varied coefficients interpretable as purities.  This is a sensitivity
  // benchmark, not yet a covariance model for calibrated tagger SFs.
  for (int row=0; row!=static_cast<int>(fitReco.size()); ++row) {
    const int reco = fitReco[row];
    for (int variedFlavor=0; variedFlavor!=nTruth; ++variedFlavor) {
      TMatrixD variedDesign[2] = {design,design};
      for (int direction=0; direction!=2; ++direction) {
        std::vector<double> varied(nTruth);
        double normalization = 0.;
        for (int flavor=0; flavor!=nTruth; ++flavor) {
          varied[flavor] = compositionData(reco,flavor);
          if (flavor==variedFlavor)
            varied[flavor] *= direction==0 ? 1.05 : 0.95;
          normalization += varied[flavor];
        }
        if (!(normalization>0.)) continue;
        for (int flavor=0; flavor!=nTruth; ++flavor)
          variedDesign[direction](row,flavor) = varied[flavor]/normalization*
            mcResponse[reco][flavor].mean;
      }
      FitResult plusFit = fitResponse(
        variedDesign[0],observed,observedError,responsePriorSigma);
      FitResult minusFit = fitResponse(
        variedDesign[1],observed,observedError,responsePriorSigma);
      if (!plusFit.solved || !minusFit.solved) continue;
      for (int resultFlavor=0; resultFlavor!=nTruth; ++resultFlavor) {
        const double shift = 0.5*(plusFit.residual[resultFlavor]-
                                  minusFit.residual[resultFlavor]);
        taggingEfficiency5PctUncertainty[resultFlavor] += shift*shift;
      }
    }
  }
  for (int flavor=0; flavor!=nTruth; ++flavor)
    taggingEfficiency5PctUncertainty[flavor] =
      std::sqrt(taggingEfficiency5PctUncertainty[flavor]);

  for (int variedFlavor=0; variedFlavor!=nTruth; ++variedFlavor) {
    TMatrixD variedDesign[2] = {design,design};
    for (int direction=0; direction!=2; ++direction)
      for (int row=0; row!=static_cast<int>(fitReco.size()); ++row) {
        const int reco = fitReco[row];
        std::vector<double> varied(nTruth);
        double normalization = 0.;
        for (int flavor=0; flavor!=nTruth; ++flavor) {
          varied[flavor] = compositionData(reco,flavor)*
            (flavor==variedFlavor
              ? (direction==0 ? 1.10 : 0.90) : 1.00);
          normalization += varied[flavor];
        }
        for (int flavor=0; flavor!=nTruth; ++flavor)
          variedDesign[direction](row,flavor) = varied[flavor]/normalization*
            mcResponse[reco][flavor].mean;
      }
    FitResult plusFit = fitResponse(
      variedDesign[0],observed,observedError,responsePriorSigma);
    FitResult minusFit = fitResponse(
      variedDesign[1],observed,observedError,responsePriorSigma);
    if (!plusFit.solved || !minusFit.solved) continue;
    for (int resultFlavor=0; resultFlavor!=nTruth; ++resultFlavor) {
      const double shift = 0.5*(plusFit.residual[resultFlavor]-
                                minusFit.residual[resultFlavor]);
      flavorFraction10PctSensitivity[resultFlavor] += shift*shift;
    }
  }
  for (int flavor=0; flavor!=nTruth; ++flavor)
    flavorFraction10PctSensitivity[flavor] =
      std::sqrt(flavorFraction10PctSensitivity[flavor]);

  auto fitPtRange = [&](double lowPt, double highPt,
                        const char *variant) {
    PtRangeResult rangeResult(nReco,nTruth);
    TMatrixD rangeMC(nReco,nTruth);
    TVectorD rangeData(nReco);
    rangeMC.Zero();
    rangeData.Zero();
    double rangeTotalMC = 0.;
    double rangeTotalData = 0.;
    for (int reco=0; reco!=nReco; ++reco) {
      rangeData[reco] = std::max(0.,readComponent(
        dataFile.get(),"m0",variant,recoIds[reco],
        std::vector<int>{0},lowPt,highPt).sumWeights);
      rangeTotalData += rangeData[reco];
      for (int flavor=0; flavor!=nTruth; ++flavor) {
        rangeMC(reco,flavor) = std::max(0.,readComponent(
          mcFile.get(),"m0",variant,recoIds[reco],
          truthGroups[flavor].sourceIds,lowPt,highPt).sumWeights);
        rangeTotalMC += rangeMC(reco,flavor);
      }
    }
    if (!(rangeTotalMC>0.) || !(rangeTotalData>0.)) return rangeResult;
    rangeResult.ruSlope = inclusiveRuSlope(
      mcFile.get(),variant,lowPt,highPt);

    TVectorD rangeTruthPrior(nTruth);
    TVectorD rangeRecoFraction(nReco);
    for (int reco=0; reco!=nReco; ++reco)
      rangeRecoFraction[reco] = rangeData[reco]/rangeTotalData;
    for (int flavor=0; flavor!=nTruth; ++flavor) {
      double column = 0.;
      for (int reco=0; reco!=nReco; ++reco)
        column += rangeMC(reco,flavor);
      if (!(column>0.)) return rangeResult;
      rangeTruthPrior[flavor] = column/rangeTotalMC;
    }

    TMatrixD rangeJoint(nReco,nTruth);
    const double rangePseudocount = std::max(1.e-12,rangeTotalMC*1.e-12);
    double rangeInitialTotal = 0.;
    for (int reco=0; reco!=nReco; ++reco)
      for (int flavor=0; flavor!=nTruth; ++flavor) {
        rangeJoint(reco,flavor) =
          std::max(rangeMC(reco,flavor),rangePseudocount);
        rangeInitialTotal += rangeJoint(reco,flavor);
      }
    rangeJoint *= 1./rangeInitialTotal;
    bool rangeConverged = false;
    for (int iteration=0; iteration!=10000; ++iteration) {
      for (int reco=0; reco!=nReco; ++reco) {
        double row = 0.;
        for (int flavor=0; flavor!=nTruth; ++flavor)
          row += rangeJoint(reco,flavor);
        if (row>0.) {
          const double scale = rangeRecoFraction[reco]/row;
          for (int flavor=0; flavor!=nTruth; ++flavor)
            rangeJoint(reco,flavor) *= scale;
        }
      }
      for (int flavor=0; flavor!=nTruth; ++flavor) {
        double column = 0.;
        for (int reco=0; reco!=nReco; ++reco)
          column += rangeJoint(reco,flavor);
        if (column>0.) {
          const double scale = rangeTruthPrior[flavor]/column;
          for (int reco=0; reco!=nReco; ++reco)
            rangeJoint(reco,flavor) *= scale;
        }
      }
      double maximumError = 0.;
      for (int reco=0; reco!=nReco; ++reco) {
        double row = 0.;
        for (int flavor=0; flavor!=nTruth; ++flavor)
          row += rangeJoint(reco,flavor);
        maximumError = std::max(
          maximumError,std::fabs(row-rangeRecoFraction[reco]));
      }
      if (maximumError<1.e-10) {
        rangeConverged = true;
        break;
      }
    }
    if (!rangeConverged) return rangeResult;

    TMatrixD rangeComposition(nReco,nTruth);
    TMatrixD rangeCompositionMC(nReco,nTruth);
    for (int reco=0; reco!=nReco; ++reco) {
      double mcRow = 0.;
      for (int flavor=0; flavor!=nTruth; ++flavor)
        mcRow += rangeMC(reco,flavor);
      for (int flavor=0; flavor!=nTruth; ++flavor) {
        rangeComposition(reco,flavor) = rangeRecoFraction[reco]>0.
          ? rangeJoint(reco,flavor)/rangeRecoFraction[reco] : 0.;
        rangeCompositionMC(reco,flavor) = mcRow>0.
          ? rangeMC(reco,flavor)/mcRow : 0.;
      }
    }

    std::vector<ProfileSummary> rangeDataResponse(nReco);
    std::vector<std::vector<ProfileSummary> > rangeMCResponse(
      nReco,std::vector<ProfileSummary>(nTruth));
    for (int reco=0; reco!=nReco; ++reco) {
      rangeDataResponse[reco] = readHDM(
        dataFile.get(),recoIds[reco],std::vector<int>{0},lowPt,highPt,
        variant);
      rangeResult.dataTagResponse[reco] = rangeDataResponse[reco];
      rangeResult.mcRawTagResponse[reco] = readHDM(
        mcFile.get(),recoIds[reco],allTruthSourceIds,lowPt,highPt,variant);
      for (int flavor=0; flavor!=nTruth; ++flavor)
        rangeMCResponse[reco][flavor] = readHDM(
          mcFile.get(),recoIds[reco],truthGroups[flavor].sourceIds,
          lowPt,highPt,variant);
    }
    for (int flavor=0; flavor!=nTruth; ++flavor) {
      double numerator = 0.;
      double denominator = 0.;
      for (int reco=0; reco!=nReco; ++reco) {
        if (!rangeMCResponse[reco][flavor].valid) continue;
        numerator += rangeMC(reco,flavor)*
                     rangeMCResponse[reco][flavor].mean;
        denominator += rangeMC(reco,flavor);
      }
      if (!(denominator>0.)) return rangeResult;
      const double fallback = numerator/denominator;
      for (int reco=0; reco!=nReco; ++reco)
        if (!rangeMCResponse[reco][flavor].valid) {
          rangeMCResponse[reco][flavor].valid = true;
          rangeMCResponse[reco][flavor].mean = fallback;
          rangeMCResponse[reco][flavor].error = 0.;
        }
    }

    for (int reco=0; reco!=nReco; ++reco) {
      ProfileSummary corrected;
      double variance = 0.;
      double sumPurity = 0.;
      for (int flavor=0; flavor!=nTruth; ++flavor) {
        if (!rangeMCResponse[reco][flavor].valid) continue;
        const double purity = rangeComposition(reco,flavor);
        corrected.mean += purity*rangeMCResponse[reco][flavor].mean;
        variance += std::pow(
          purity*rangeMCResponse[reco][flavor].error,2);
        sumPurity += purity;
      }
      corrected.error = std::sqrt(std::max(0.,variance));
      corrected.sumWeights = rangeData[reco];
      corrected.valid = sumPurity>1.-1.e-6 &&
        std::isfinite(corrected.mean) && std::isfinite(corrected.error);
      rangeResult.mcCompositionCorrectedTagResponse[reco] = corrected;
      rangeResult.rawTagRatio[reco] = ratioSummary(
        rangeDataResponse[reco],rangeResult.mcRawTagResponse[reco]);
      rangeResult.compositionCorrectedTagRatio[reco] = ratioSummary(
        rangeDataResponse[reco],corrected);
    }

    std::vector<int> rangeFitReco;
    for (int reco=0; reco!=nReco; ++reco)
      if (rangeRecoFraction[reco]>0. && rangeDataResponse[reco].valid)
        rangeFitReco.push_back(reco);
    TMatrixD rangeDesign(rangeFitReco.size(),nTruth);
    TVectorD rangeObserved(rangeFitReco.size());
    TVectorD rangeErrors(rangeFitReco.size());
    for (int row=0; row!=static_cast<int>(rangeFitReco.size()); ++row) {
      const int reco = rangeFitReco[row];
      rangeObserved[row] = rangeDataResponse[reco].mean;
      rangeErrors[row] = rangeDataResponse[reco].error;
      for (int flavor=0; flavor!=nTruth; ++flavor)
        rangeDesign(row,flavor) = rangeComposition(reco,flavor)*
          rangeMCResponse[reco][flavor].mean;
    }
    rangeResult.responseFit = fitResponse(
      rangeDesign,rangeObserved,rangeErrors,responsePriorSigma);

    // Diagnostic alternate HDM fit: derive one merge-safe zero-intercept
    // effective R_u slope from truth-matched MC in this pT interval, then use
    // that same response in the data and MC HDM definitions.  This is not a
    // data-only calibration and is intentionally kept alongside, rather than
    // replacing, the nominal R_u=0.92 result.
    if (std::isfinite(rangeResult.ruSlope)) {
      std::vector<ProfileSummary> slopeDataResponse(nReco);
      std::vector<std::vector<ProfileSummary> > slopeMCResponse(
        nReco,std::vector<ProfileSummary>(nTruth));
      for (int reco=0; reco!=nReco; ++reco) {
        slopeDataResponse[reco] = readHDM(
          dataFile.get(),recoIds[reco],std::vector<int>{0},lowPt,highPt,
          variant,rangeResult.ruSlope);
        for (int flavor=0; flavor!=nTruth; ++flavor)
          slopeMCResponse[reco][flavor] = readHDM(
            mcFile.get(),recoIds[reco],truthGroups[flavor].sourceIds,
            lowPt,highPt,variant,rangeResult.ruSlope);
      }
      for (int flavor=0; flavor!=nTruth; ++flavor) {
        double numerator = 0.;
        double denominator = 0.;
        for (int reco=0; reco!=nReco; ++reco) {
          if (!slopeMCResponse[reco][flavor].valid) continue;
          numerator += rangeMC(reco,flavor)*
                       slopeMCResponse[reco][flavor].mean;
          denominator += rangeMC(reco,flavor);
        }
        if (!(denominator>0.)) continue;
        const double fallback = numerator/denominator;
        for (int reco=0; reco!=nReco; ++reco)
          if (!slopeMCResponse[reco][flavor].valid) {
            slopeMCResponse[reco][flavor].valid = true;
            slopeMCResponse[reco][flavor].mean = fallback;
            slopeMCResponse[reco][flavor].error = 0.;
          }
      }
      std::vector<int> slopeFitReco;
      for (int reco=0; reco!=nReco; ++reco)
        if (rangeRecoFraction[reco]>0. && slopeDataResponse[reco].valid)
          slopeFitReco.push_back(reco);
      TMatrixD slopeDesign(slopeFitReco.size(),nTruth);
      TVectorD slopeObserved(slopeFitReco.size());
      TVectorD slopeErrors(slopeFitReco.size());
      for (int row=0; row!=static_cast<int>(slopeFitReco.size()); ++row) {
        const int reco = slopeFitReco[row];
        slopeObserved[row] = slopeDataResponse[reco].mean;
        slopeErrors[row] = slopeDataResponse[reco].error;
        for (int flavor=0; flavor!=nTruth; ++flavor)
          slopeDesign(row,flavor) = rangeComposition(reco,flavor)*
            slopeMCResponse[reco][flavor].mean;
      }
      rangeResult.responseFitRuSlope = fitResponse(
        slopeDesign,slopeObserved,slopeErrors,responsePriorSigma);
    }

    auto fitRecoilFraction = [&](const char *component, double scale,
                                 FitResult &componentFit,
                                 std::vector<ProfileSummary> &mcTruth,
                                 std::vector<ProfileSummary> &dataTruth,
                                 std::vector<ProfileSummary> &ratio) {
      std::vector<ProfileSummary> dataComponent(nReco);
      std::vector<std::vector<ProfileSummary> > mcComponent(
        nReco,std::vector<ProfileSummary>(nTruth));
      auto readFraction = [&](TFile *file, int reco,
                              const std::vector<int> &truthIds) {
        ProfileSummary direct = readComponent(
          file,component,variant,reco,truthIds,lowPt,highPt,scale);
        if (direct.valid || std::string(component)!="fnu") return direct;
        return sumSummaries(
          readComponent(file,"mn",variant,reco,truthIds,lowPt,highPt,
                        1./ZJetFlavorMatrix::responseN),
          readComponent(file,"mu",variant,reco,truthIds,lowPt,highPt,
                        1./ZJetFlavorMatrix::responseU));
      };
      for (int reco=0; reco!=nReco; ++reco) {
        dataComponent[reco] = readFraction(
          dataFile.get(),recoIds[reco],std::vector<int>{0});
        for (int flavor=0; flavor!=nTruth; ++flavor)
          mcComponent[reco][flavor] = readFraction(
            mcFile.get(),recoIds[reco],truthGroups[flavor].sourceIds);
      }

      for (int flavor=0; flavor!=nTruth; ++flavor) {
        double weightedMean = 0.;
        double weightedVariance = 0.;
        double weights = 0.;
        for (int reco=0; reco!=nReco; ++reco) {
          if (!mcComponent[reco][flavor].valid) continue;
          const double count = rangeMC(reco,flavor);
          weightedMean += count*mcComponent[reco][flavor].mean;
          weightedVariance += std::pow(
            count*mcComponent[reco][flavor].error,2);
          weights += count;
        }
        if (!(weights>0.)) return;
        mcTruth[flavor].mean = weightedMean/weights;
        mcTruth[flavor].error = std::sqrt(weightedVariance)/weights;
        mcTruth[flavor].sumWeights = weights;
        mcTruth[flavor].valid = std::isfinite(mcTruth[flavor].mean) &&
                                std::isfinite(mcTruth[flavor].error);
        for (int reco=0; reco!=nReco; ++reco) {
          if (mcComponent[reco][flavor].valid) continue;
          mcComponent[reco][flavor] = mcTruth[flavor];
        }
      }

      std::vector<int> componentReco;
      for (int reco=0; reco!=nReco; ++reco)
        if (rangeRecoFraction[reco]>0. && dataComponent[reco].valid)
          componentReco.push_back(reco);
      TMatrixD componentDesign(componentReco.size(),nTruth);
      TVectorD componentObserved(componentReco.size());
      TVectorD componentErrors(componentReco.size());
      for (int row=0; row!=static_cast<int>(componentReco.size()); ++row) {
        const int reco = componentReco[row];
        componentObserved[row] = dataComponent[reco].mean;
        componentErrors[row] = dataComponent[reco].error;
        for (int flavor=0; flavor!=nTruth; ++flavor)
          componentDesign(row,flavor) = rangeComposition(reco,flavor)*
            mcComponent[reco][flavor].mean;
      }
      componentFit = fitResponse(
        componentDesign,componentObserved,componentErrors,responsePriorSigma);
      if (!componentFit.solved) return;
      for (int flavor=0; flavor!=nTruth; ++flavor) {
        const double fitError = componentFit.covariance(flavor,flavor)>0.
          ? std::sqrt(componentFit.covariance(flavor,flavor)) : 0.;
        ratio[flavor].mean = componentFit.residual[flavor];
        ratio[flavor].error = fitError;
        ratio[flavor].valid = std::isfinite(ratio[flavor].mean) &&
                              std::isfinite(ratio[flavor].error);
        dataTruth[flavor].mean = mcTruth[flavor].mean*ratio[flavor].mean;
        dataTruth[flavor].error = std::hypot(
          ratio[flavor].mean*mcTruth[flavor].error,
          mcTruth[flavor].mean*ratio[flavor].error);
        dataTruth[flavor].sumWeights = rangeTotalData;
        dataTruth[flavor].valid = mcTruth[flavor].valid &&
          std::isfinite(dataTruth[flavor].mean) &&
          std::isfinite(dataTruth[flavor].error);
      }
    };

    fitRecoilFraction(
      "mn",1./ZJetFlavorMatrix::responseN,rangeResult.fsrScaleFit,
      rangeResult.mcFsrFraction,rangeResult.dataFsrFraction,
      rangeResult.fsrRatio);
    fitRecoilFraction(
      "mu",1./ZJetFlavorMatrix::responseU,rangeResult.ueScaleFit,
      rangeResult.mcUeFraction,rangeResult.dataUeFraction,
      rangeResult.ueRatio);
    // m_n+m_u is independent of the arbitrary 15 GeV Type-I boundary.  Keep
    // it unscaled so a later UE-hole subtraction and effective R_u model can
    // be applied transparently.
    fitRecoilFraction(
      "mnu",1.,rangeResult.mnuScaleFit,
      rangeResult.mcMnuFraction,rangeResult.dataMnuFraction,
      rangeResult.mnuRatio);
    fitRecoilFraction(
      "fnu",1.,rangeResult.fnuScaleFit,
      rangeResult.mcFnuFraction,rangeResult.dataFnuFraction,
      rangeResult.fnuRatio);
    // Available from productions made after the projected-area control was
    // introduced. Older files simply leave this optional fit unsolved.
    fitRecoilFraction(
      "mnufsr",1.,rangeResult.mnuFsrScaleFit,
      rangeResult.mcMnuFsrFraction,rangeResult.dataMnuFsrFraction,
      rangeResult.mnuFsrRatio);
    return rangeResult;
  };

  const std::vector<double> ptFitEdges = {
    30.,40.,60.,85.,125.,180.,250.,400.,600.,
  };
  std::map<std::string,std::vector<PtRangeResult> > variantPtFits;
  for (const std::string variant : {"ab","ad","tc","pf"})
    for (size_t bin=0; bin+1<ptFitEdges.size(); ++bin)
      variantPtFits[variant].push_back(
        fitPtRange(ptFitEdges[bin],ptFitEdges[bin+1],variant.c_str()));
  const std::vector<PtRangeResult> &ptFits = variantPtFits.at("tc");

  std::vector<double> responsePrediction(nReco,
    std::numeric_limits<double>::quiet_NaN());
  for (int reco=0; reco!=nReco; ++reco) {
    if (dataRecoFraction[reco]<=0.) continue;
    double prediction = 0.;
    for (int flavor=0; flavor!=nTruth; ++flavor)
      prediction += compositionData(reco,flavor)*
                    mcResponse[reco][flavor].mean*fit.residual[flavor];
    responsePrediction[reco] = prediction;
  }

  if (gSystem->mkdir(outputDirectory,true)!=0 &&
      gSystem->AccessPathName(outputDirectory))
    throw std::runtime_error(std::string("Failed to create ")+
                             outputDirectory);
  const std::string rootOutput =
    std::string(outputDirectory)+"/flavorMatrixAnalysis.root";
  std::unique_ptr<TFile> output(TFile::Open(rootOutput.c_str(),"RECREATE"));
  if (!output || output->IsZombie())
    throw std::runtime_error("Failed to create "+rootOutput);

  TDirectory *tagging = output->mkdir("tagging");
  TDirectory *response = output->mkdir("response");
  TDirectory *diagnostics = output->mkdir("diagnostics");
  const double idBins[] = {-0.5,0.5,1.5,2.5,3.5,4.5,5.5,6.5};
  auto makeMatrix = [&](const char *name, const char *title) {
    TH2D *histogram = new TH2D(name,title,7,idBins,7,idBins);
    histogram->SetDirectory(tagging);
    labelRecoAxis(histogram->GetXaxis());
    labelTruthAxis(histogram->GetYaxis());
    return histogram;
  };
  TH2D *hEfficiencyMC = makeMatrix(
    "h2_efficiency_mc",";Reco hybrid flavor;True flavor;#epsilon^{MC}_{t|f}");
  TH2D *hEfficiencyData = makeMatrix(
    "h2_efficiency_data_inferred",
    ";Reco hybrid flavor;True flavor;Inferred #epsilon^{data}_{t|f}");
  TH2D *hTransitionSF = makeMatrix(
    "h2_transition_sf",";Reco hybrid flavor;True flavor;Transition SF");
  TH2D *hJointMC = makeMatrix(
    "h2_joint_mc",";Reco hybrid flavor;True flavor;P_{MC}(t,f)");
  TH2D *hJointData = makeMatrix(
    "h2_joint_data_inferred",
    ";Reco hybrid flavor;True flavor;Inferred P_{data}(t,f)");
  TH2D *hCompositionData = makeMatrix(
    "h2_composition_data",
    ";Reco hybrid flavor;True flavor;Inferred P_{data}(f|t)");
  TH2D *hCompositionMC = makeMatrix(
    "h2_composition_mc",
    ";Reco hybrid flavor;True flavor;P_{MC}(f|t)");
  TH2D *hCompositionRatio = makeMatrix(
    "h2_composition_ratio_data_over_mc",
    ";Reco hybrid flavor;True flavor;P_{data}(f|t)/P_{MC}(f|t)");
  TH1D *hTruthPrior = new TH1D(
    "h1_truth_prior_mc",";True flavor;MC truth prior",7,idBins);
  TH1D *hRecoFractionData = new TH1D(
    "h1_reco_fraction_data",";Reco hybrid flavor;Data tag fraction",7,idBins);
  hTruthPrior->SetDirectory(tagging);
  hRecoFractionData->SetDirectory(tagging);
  labelTruthAxis(hTruthPrior->GetXaxis());
  labelRecoAxis(hRecoFractionData->GetXaxis());
  for (int reco=0; reco!=nReco; ++reco) {
    const int bin = hRecoFractionData->GetXaxis()->FindFixBin(recoIds[reco]);
    hRecoFractionData->SetBinContent(bin,dataRecoFraction[reco]);
  }
  for (int flavor=0; flavor!=nTruth; ++flavor) {
    const int trueId = truthGroups[flavor].outputId;
    const int bin = hTruthPrior->GetXaxis()->FindFixBin(trueId);
    hTruthPrior->SetBinContent(bin,truthPrior[flavor]);
  }
  for (int reco=0; reco!=nReco; ++reco)
    for (int flavor=0; flavor!=nTruth; ++flavor) {
      const int xbin = hEfficiencyMC->GetXaxis()->FindFixBin(recoIds[reco]);
      const int trueId = truthGroups[flavor].outputId;
      const int ybin = hEfficiencyMC->GetYaxis()->FindFixBin(trueId);
      hEfficiencyMC->SetBinContent(xbin,ybin,efficiencyMC(reco,flavor));
      hEfficiencyData->SetBinContent(xbin,ybin,efficiencyData(reco,flavor));
      hTransitionSF->SetBinContent(xbin,ybin,transitionSF(reco,flavor));
      hJointMC->SetBinContent(xbin,ybin,jointMC(reco,flavor));
      hJointData->SetBinContent(xbin,ybin,inferredJoint(reco,flavor));
      hCompositionMC->SetBinContent(
        xbin,ybin,compositionMC(reco,flavor));
      hCompositionData->SetBinContent(xbin,ybin,compositionData(reco,flavor));
      hCompositionRatio->SetBinContent(
        xbin,ybin,compositionRatio(reco,flavor));
    }

  response->cd();
  TH2D *hResponseMC = new TH2D(
    "h2_response_mc_by_transition",
    ";Reco hybrid flavor;True flavor;MC HDM response",7,idBins,7,idBins);
  TH2D *hResponseData = new TH2D(
    "h2_response_data_estimated_by_transition",
    ";Reco hybrid flavor;True flavor;Estimated data HDM response",
    7,idBins,7,idBins);
  TH2D *hResponseRatio = new TH2D(
    "h2_response_ratio_data_over_mc_by_transition",
    ";Reco hybrid flavor;True flavor;Estimated data/MC HDM response",
    7,idBins,7,idBins);
  labelRecoAxis(hResponseMC->GetXaxis());
  labelTruthAxis(hResponseMC->GetYaxis());
  labelRecoAxis(hResponseData->GetXaxis());
  labelTruthAxis(hResponseData->GetYaxis());
  labelRecoAxis(hResponseRatio->GetXaxis());
  labelTruthAxis(hResponseRatio->GetYaxis());
  for (int reco=0; reco!=nReco; ++reco)
    for (int flavor=0; flavor!=nTruth; ++flavor) {
      const int xbin = hResponseMC->GetXaxis()->FindFixBin(recoIds[reco]);
      const int trueId = truthGroups[flavor].outputId;
      const int ybin = hResponseMC->GetYaxis()->FindFixBin(trueId);
      hResponseMC->SetBinContent(xbin,ybin,mcResponse[reco][flavor].mean);
      hResponseMC->SetBinError(xbin,ybin,mcResponse[reco][flavor].error);
      hResponseData->SetBinContent(
        xbin,ybin,mcResponse[reco][flavor].mean*fit.residual[flavor]);
      hResponseData->SetBinError(
        xbin,ybin,std::hypot(
          fit.residual[flavor]*mcResponse[reco][flavor].error,
          mcResponse[reco][flavor].mean*
            std::sqrt(std::max(0.,fit.covariance(flavor,flavor)))));
      hResponseRatio->SetBinContent(xbin,ybin,fit.residual[flavor]);
      hResponseRatio->SetBinError(
        xbin,ybin,std::sqrt(std::max(0.,fit.covariance(flavor,flavor))));
    }
  TH1D *hDataReco = new TH1D(
    "h1_response_data_reco",";Reco hybrid flavor;Data HDM response",7,idBins);
  TH1D *hPredictionReco = new TH1D(
    "h1_response_prediction_reco",
    ";Reco hybrid flavor;Fitted data-response prediction",7,idBins);
  labelRecoAxis(hDataReco->GetXaxis());
  labelRecoAxis(hPredictionReco->GetXaxis());
  for (int reco=0; reco!=nReco; ++reco) {
    const int bin = hDataReco->GetXaxis()->FindFixBin(recoIds[reco]);
    if (dataResponse[reco].valid) {
      hDataReco->SetBinContent(bin,dataResponse[reco].mean);
      hDataReco->SetBinError(bin,dataResponse[reco].error);
    }
    if (std::isfinite(responsePrediction[reco]))
      hPredictionReco->SetBinContent(bin,responsePrediction[reco]);
  }
  TH1D *hResidual = new TH1D(
    "h1_response_residual_data_over_mc",
    ";True flavor;R^{data}_{f}/R^{MC}_{f}",7,idBins);
  TH1D *hCorrection = new TH1D(
    "h1_response_correction_mc_over_data",
    ";True flavor;R^{MC}_{f}/R^{data}_{f}",7,idBins);
  labelTruthAxis(hResidual->GetXaxis());
  labelTruthAxis(hCorrection->GetXaxis());
  TGraphErrors *gResidual = new TGraphErrors(nTruth);
  gResidual->SetName("g_response_residual_data_over_mc");
  gResidual->SetTitle(";True-flavor ID;R^{data}_{f}/R^{MC}_{f}");
  for (int flavor=0; flavor!=nTruth; ++flavor) {
    const int trueId = truthGroups[flavor].outputId;
    const int bin = hResidual->GetXaxis()->FindFixBin(trueId);
    const double value = fit.residual[flavor];
    const double error = fit.covariance(flavor,flavor)>0.
      ? std::sqrt(fit.covariance(flavor,flavor)) : 0.;
    hResidual->SetBinContent(bin,value);
    hResidual->SetBinError(bin,error);
    if (value!=0.) {
      hCorrection->SetBinContent(bin,1./value);
      hCorrection->SetBinError(bin,error/(value*value));
    }
    gResidual->SetPoint(flavor,trueId,value);
    gResidual->SetPointError(flavor,0.,error);
  }
  gResidual->Write();
  fit.covariance.Write("response_residual_covariance");
  for (const auto &variantEntry : variantPtFits) {
    const std::string &variant = variantEntry.first;
    const std::vector<PtRangeResult> &fits = variantEntry.second;
    TGraphErrors ruSlopeGraph;
    ruSlopeGraph.SetName(Form("g_ru_slope_used_vs_pt_%s",variant.c_str()));
    ruSlopeGraph.SetTitle(
      ";reference p_{T} (GeV);MC-derived effective R_{u} slope");
    for (size_t bin=0; bin<fits.size(); ++bin) {
      if (!std::isfinite(fits[bin].ruSlope)) continue;
      const double low = ptFitEdges[bin];
      const double high = ptFitEdges[bin+1];
      const int point = ruSlopeGraph.GetN();
      ruSlopeGraph.SetPoint(point,std::sqrt(low*high),fits[bin].ruSlope);
      ruSlopeGraph.SetPointError(point,0.5*(high-low),0.);
    }
    ruSlopeGraph.Write();
    for (int flavor=0; flavor!=nTruth; ++flavor) {
      TGraphErrors *graph = new TGraphErrors();
      TGraphErrors *slopeGraph = new TGraphErrors();
      graph->SetName(Form("g_response_residual_vs_pt_%s_%s",
                          variant.c_str(),truthGroups[flavor].name));
      slopeGraph->SetName(Form(
        "g_response_residual_ru_slope_vs_pt_%s_%s",
        variant.c_str(),truthGroups[flavor].name));
      graph->SetTitle(Form(
        ";reference p_{T} (GeV);R^{data}_{%s}/R^{MC}_{%s}",
        truthGroups[flavor].name,truthGroups[flavor].name));
      slopeGraph->SetTitle(graph->GetTitle());
      for (size_t bin=0; bin<fits.size(); ++bin) {
        // Keep the graph quantitative: rank-deficient or nearly singular
        // bins remain in the TSV diagnostics but are not drawn.
        const FitResult &ptFit = fits[bin].responseFit;
        if (!ptFit.solved || ptFit.rank<nTruth ||
            !std::isfinite(ptFit.nonzeroCondition) ||
            ptFit.nonzeroCondition>100.)
          continue;
        const double low = ptFitEdges[bin];
        const double high = ptFitEdges[bin+1];
        const double error = ptFit.covariance(flavor,flavor)>0.
          ? std::sqrt(ptFit.covariance(flavor,flavor)) : 0.;
        const int point = graph->GetN();
        graph->SetPoint(point,std::sqrt(low*high),ptFit.residual[flavor]);
        graph->SetPointError(point,0.5*(high-low),error);
        const FitResult &slopeFit = fits[bin].responseFitRuSlope;
        if (slopeFit.solved && slopeFit.rank>=nTruth &&
            std::isfinite(slopeFit.nonzeroCondition) &&
            slopeFit.nonzeroCondition<=100.) {
          const double slopeError = slopeFit.covariance(flavor,flavor)>0.
            ? std::sqrt(slopeFit.covariance(flavor,flavor)) : 0.;
          const int slopePoint = slopeGraph->GetN();
          slopeGraph->SetPoint(
            slopePoint,std::sqrt(low*high),slopeFit.residual[flavor]);
          slopeGraph->SetPointError(
            slopePoint,0.5*(high-low),slopeError);
        }
      }
      graph->Write();
      slopeGraph->Write();
      if (variant=="tc")
        graph->Write(Form("g_response_residual_vs_pt_%s",
                          truthGroups[flavor].name));
    }
  }

  auto addSummaryPoint = [&](TGraphErrors *graph, size_t ptBin,
                             const ProfileSummary &summary) {
    if (!graph || !summary.valid) return;
    const double low = ptFitEdges[ptBin];
    const double high = ptFitEdges[ptBin+1];
    const int point = graph->GetN();
    graph->SetPoint(point,std::sqrt(low*high),summary.mean);
    graph->SetPointError(point,0.5*(high-low),summary.error);
  };
  for (int reco=0; reco!=nReco; ++reco) {
    const std::string tag = recoName(recoIds[reco]);
    TGraphErrors *dataGraph = new TGraphErrors();
    TGraphErrors *rawMcGraph = new TGraphErrors();
    TGraphErrors *correctedMcGraph = new TGraphErrors();
    TGraphErrors *rawRatioGraph = new TGraphErrors();
    TGraphErrors *correctedRatioGraph = new TGraphErrors();
    dataGraph->SetName(Form("g_tag_response_data_vs_pt_%s",tag.c_str()));
    rawMcGraph->SetName(Form(
      "g_tag_response_mc_raw_vs_pt_%s",tag.c_str()));
    correctedMcGraph->SetName(Form(
      "g_tag_response_mc_composition_corrected_vs_pt_%s",tag.c_str()));
    rawRatioGraph->SetName(Form(
      "g_tag_response_raw_data_over_mc_vs_pt_%s",tag.c_str()));
    correctedRatioGraph->SetName(Form(
      "g_tag_response_composition_corrected_data_over_mc_vs_pt_%s",
      tag.c_str()));
    for (size_t bin=0; bin<ptFits.size(); ++bin) {
      addSummaryPoint(dataGraph,bin,ptFits[bin].dataTagResponse[reco]);
      addSummaryPoint(rawMcGraph,bin,ptFits[bin].mcRawTagResponse[reco]);
      addSummaryPoint(correctedMcGraph,bin,
        ptFits[bin].mcCompositionCorrectedTagResponse[reco]);
      addSummaryPoint(rawRatioGraph,bin,ptFits[bin].rawTagRatio[reco]);
      addSummaryPoint(correctedRatioGraph,bin,
        ptFits[bin].compositionCorrectedTagRatio[reco]);
    }
    dataGraph->Write();
    rawMcGraph->Write();
    correctedMcGraph->Write();
    rawRatioGraph->Write();
    correctedRatioGraph->Write();
  }

  auto componentFit = [](const PtRangeResult &range,
                         const std::string &component) -> const FitResult& {
    if (component=="fsr") return range.fsrScaleFit;
    if (component=="ue") return range.ueScaleFit;
    if (component=="mnu") return range.mnuScaleFit;
    if (component=="fnu") return range.fnuScaleFit;
    return range.mnuFsrScaleFit;
  };
  auto componentValues = [](const PtRangeResult &range,
                            const std::string &component, int kind)
      -> const std::vector<ProfileSummary>& {
    if (component=="fsr")
      return kind==0 ? range.mcFsrFraction
                     : (kind==1 ? range.dataFsrFraction : range.fsrRatio);
    if (component=="ue")
      return kind==0 ? range.mcUeFraction
                     : (kind==1 ? range.dataUeFraction : range.ueRatio);
    if (component=="mnu")
      return kind==0 ? range.mcMnuFraction
                     : (kind==1 ? range.dataMnuFraction : range.mnuRatio);
    if (component=="fnu")
      return kind==0 ? range.mcFnuFraction
                     : (kind==1 ? range.dataFnuFraction : range.fnuRatio);
    return kind==0 ? range.mcMnuFsrFraction
                   : (kind==1 ? range.dataMnuFsrFraction
                              : range.mnuFsrRatio);
  };
  for (const auto &variantEntry : variantPtFits) {
    const std::string &variant = variantEntry.first;
    const std::vector<PtRangeResult> &fits = variantEntry.second;
    for (int flavor=0; flavor!=nTruth; ++flavor) {
      const std::string name = truthGroups[flavor].name;
      for (const std::string component : {"fsr","ue","mnu","fnu","mnufsr"}) {
        TGraphErrors *mcGraph = new TGraphErrors();
        TGraphErrors *dataGraph = new TGraphErrors();
        TGraphErrors *ratioGraph = new TGraphErrors();
        const std::string variantSuffix = "_"+variant+"_"+name;
        mcGraph->SetName(Form(
          "g_%s_fraction_mc_vs_pt%s",component.c_str(),
          variantSuffix.c_str()));
        dataGraph->SetName(Form(
          "g_%s_fraction_data_vs_pt%s",component.c_str(),
          variantSuffix.c_str()));
        ratioGraph->SetName(Form(
          "g_%s_fraction_ratio_data_over_mc_vs_pt%s",
          component.c_str(),variantSuffix.c_str()));
        for (size_t bin=0; bin<fits.size(); ++bin) {
          const FitResult &fitComponent = componentFit(fits[bin],component);
          if (!fitComponent.solved || fitComponent.rank<nTruth ||
              !std::isfinite(fitComponent.nonzeroCondition) ||
              fitComponent.nonzeroCondition>100.)
            continue;
          addSummaryPoint(mcGraph,bin,
                          componentValues(fits[bin],component,0)[flavor]);
          addSummaryPoint(dataGraph,bin,
                          componentValues(fits[bin],component,1)[flavor]);
          addSummaryPoint(ratioGraph,bin,
                          componentValues(fits[bin],component,2)[flavor]);
        }
        mcGraph->Write();
        dataGraph->Write();
        ratioGraph->Write();
        if (variant=="tc") {
          mcGraph->Write(Form("g_%s_fraction_mc_vs_pt_%s",
                              component.c_str(),name.c_str()));
          dataGraph->Write(Form("g_%s_fraction_data_vs_pt_%s",
                                component.c_str(),name.c_str()));
          ratioGraph->Write(Form(
            "g_%s_fraction_ratio_data_over_mc_vs_pt_%s",
            component.c_str(),name.c_str()));
        }
      }
    }
  }

  // Flavor-resolved effective recoil responses for the next production.
  // The ratio of means becomes unstable when the generator component is
  // close to zero, so write both it and the zero-intercept regression slope.
  for (const std::string variant : {"ab","ad","tc","pf"})
    for (int flavor=0; flavor!=nTruth; ++flavor)
      for (const std::string component : {"n","u"}) {
        TGraphErrors meanRatio;
        TGraphErrors slope;
        meanRatio.SetName(Form(
          "g_r%s_mean_ratio_mc_vs_pt_%s_%s",component.c_str(),
          variant.c_str(),truthGroups[flavor].name));
        slope.SetName(Form(
          "g_r%s_slope_mc_vs_pt_%s_%s",component.c_str(),
          variant.c_str(),truthGroups[flavor].name));
        for (size_t bin=0; bin+1<ptFitEdges.size(); ++bin) {
          const double low = ptFitEdges[bin];
          const double high = ptFitEdges[bin+1];
          const std::vector<int> &truthIds = truthGroups[flavor].sourceIds;
          const ProfileSummary reco = readComponentAcrossReco(
            mcFile.get(),component=="n" ? "recomnmatched" : "recoumatched",
            variant.c_str(),truthIds,low,high);
          const ProfileSummary gen = readComponentAcrossReco(
            mcFile.get(),component=="n" ? "genmn" : "genmu",
            variant.c_str(),truthIds,low,high);
          const ProfileSummary product = readComponentAcrossReco(
            mcFile.get(),component=="n" ? "recogenmn" : "recogenmu",
            variant.c_str(),truthIds,low,high);
          const ProfileSummary square = readComponentAcrossReco(
            mcFile.get(),component=="n" ? "genmn2" : "genmu2",
            variant.c_str(),truthIds,low,high);
          const ProfileSummary mean = ratioSummary(reco,gen);
          const ProfileSummary fitSlope = ratioSummary(product,square);
          const double x = std::sqrt(low*high);
          if (mean.valid) {
            const int point = meanRatio.GetN();
            meanRatio.SetPoint(point,x,mean.mean);
            meanRatio.SetPointError(point,0.5*(high-low),mean.error);
          }
          if (fitSlope.valid) {
            const int point = slope.GetN();
            slope.SetPoint(point,x,fitSlope.mean);
            slope.SetPointError(point,0.5*(high-low),fitSlope.error);
          }
        }
        meanRatio.Write();
        slope.Write();
        if (variant=="tc") {
          meanRatio.Write(Form("g_r%s_mean_ratio_mc_vs_pt_%s",
                               component.c_str(),truthGroups[flavor].name));
          slope.Write(Form("g_r%s_slope_mc_vs_pt_%s",
                           component.c_str(),truthGroups[flavor].name));
        }
      }

  // Data-accessible R_u closure slope.  This uses
  // f_u^closure=1-m_2/R_2-m_n/R_n with R_2=R_n=1 and therefore remains a
  // conditional estimator.  The MC inclusive and flavor-resolved graphs
  // expose the leading-jet-resolution and composition bias directly.
  for (const std::string variant : {"ab","ad","tc","pf"}) {
    TGraphErrors dataClosure;
    TGraphErrors mcClosure;
    dataClosure.SetName(Form(
      "g_ru_closure_slope_data_vs_pt_%s",variant.c_str()));
    mcClosure.SetName(Form(
      "g_ru_closure_slope_mc_vs_pt_%s",variant.c_str()));
    std::vector<TGraphErrors> mcFlavorClosure(nTruth);
    for (int flavor=0; flavor!=nTruth; ++flavor)
      mcFlavorClosure[flavor].SetName(Form(
        "g_ru_closure_slope_mc_vs_pt_%s_%s",variant.c_str(),
        truthGroups[flavor].name));
    for (size_t bin=0; bin+1<ptFitEdges.size(); ++bin) {
      const double low = ptFitEdges[bin];
      const double high = ptFitEdges[bin+1];
      const double x = std::sqrt(low*high);
      auto appendClosure = [&](TGraphErrors &graph, TFile *file,
                               const std::vector<int> &truthIds) {
        const ProfileSummary product = readComponentAcrossReco(
          file,"mufuclosure",variant.c_str(),truthIds,low,high);
        const ProfileSummary square = readComponentAcrossReco(
          file,"fuclosure2",variant.c_str(),truthIds,low,high);
        const ProfileSummary slope = ratioSummary(product,square);
        if (!slope.valid) return;
        const int point = graph.GetN();
        graph.SetPoint(point,x,slope.mean);
        graph.SetPointError(point,0.5*(high-low),slope.error);
      };
      appendClosure(dataClosure,dataFile.get(),std::vector<int>{0});
      appendClosure(mcClosure,mcFile.get(),allTruthSourceIds);
      for (int flavor=0; flavor!=nTruth; ++flavor)
        appendClosure(mcFlavorClosure[flavor],mcFile.get(),
                      truthGroups[flavor].sourceIds);
    }
    dataClosure.Write();
    mcClosure.Write();
    for (TGraphErrors &graph : mcFlavorClosure) graph.Write();
  }

  diagnostics->cd();
  TH1D *hSingular = new TH1D(
    "h1_singular_values",";SVD mode;Singular value",
    std::max(1,fit.singularValues.GetNrows()),0.5,
    std::max(1,fit.singularValues.GetNrows())+0.5);
  for (int index=0; index!=fit.singularValues.GetNrows(); ++index)
    hSingular->SetBinContent(index+1,fit.singularValues[index]);
  TParameter<int>("response_rank",fit.rank).Write();
  TParameter<int>("response_fit_solved",fit.solved ? 1 : 0).Write();
  TParameter<int>("response_rows",fit.rows).Write();
  TParameter<int>("active_truth_flavors",nTruth).Write();
  TParameter<int>("ipf_converged",ipfConverged ? 1 : 0).Write();
  TParameter<int>("ipf_iterations",ipfIterations).Write();
  TParameter<double>("ipf_maximum_marginal_error",ipfMaximumError).Write();
  TParameter<double>("response_nonzero_condition_number",
                     fit.nonzeroCondition).Write();
  TParameter<double>("response_full_condition_number",
                     fit.fullCondition).Write();
  TParameter<double>("response_chi2",fit.chi2).Write();
  TParameter<double>("response_prior_sigma",responsePriorSigma).Write();
  TParameter<double>("minimum_pt",minimumPt).Write();
  TParameter<double>("maximum_pt",maximumPt).Write();
  TParameter<int>("hdm_constructed_from_component_means",1).Write();
  TParameter<int>("sparse_cell_response_fallbacks",cellFallbacks).Write();
  TParameter<int>("zero_profile_error_fallbacks",fit.zeroErrorFallbacks).Write();
  TH2D *hUncertainty = new TH2D(
    "h2_response_uncertainty_components",
    ";True flavor;Source;Absolute uncertainty",
    nTruth,-0.5,nTruth-0.5,5,-0.5,4.5);
  const char *uncertaintySources[] = {
    "data stat.", "purity stat.",
    "response stat.", "efficiency 5%", "fractions 10%",
  };
  for (int flavor=0; flavor!=nTruth; ++flavor) {
    hUncertainty->GetXaxis()->SetBinLabel(
      flavor+1,truthGroups[flavor].name);
    hUncertainty->SetBinContent(
      flavor+1,1,std::sqrt(std::max(0.,fit.covariance(flavor,flavor))));
    hUncertainty->SetBinContent(
      flavor+1,2,purityStatUncertainty[flavor]);
    hUncertainty->SetBinContent(
      flavor+1,3,responseTemplateUncertainty[flavor]);
    hUncertainty->SetBinContent(
      flavor+1,4,taggingEfficiency5PctUncertainty[flavor]);
    hUncertainty->SetBinContent(
      flavor+1,5,flavorFraction10PctSensitivity[flavor]);
  }
  for (int source=0; source!=5; ++source)
    hUncertainty->GetYaxis()->SetBinLabel(
      source+1,uncertaintySources[source]);
  TObjString(
    "Data transition efficiencies are the KL/IPF projection closest to the "
    "MC joint matrix with fixed MC truth marginals and observed data reco "
    "marginals. Tagging marginals use the un-subtracted parallel population "
    "to keep a non-negative likelihood. The undefined reco tag is omitted; "
    "truth d/u+s and undefined+g are combined into four physical groups. "
    "Individual transition SFs "
    "are model dependent and are not independently identified by data tag "
    "fractions.")
    .Write("tagging_inference_model",TObject::kOverwrite);
  TObjString(
    "Response fit assumes one multiplicative data/MC residual per true "
    "flavor, common to all reconstructed tags. HDM is constructed after "
    "merging from the m0, mn and mu profile means with Rn=1 and Ru=0.92; "
    "the event-wise p3hdm profile is deliberately not used. A separate "
    "diagnostic rebuild uses the zero-intercept truth-matched MC Ru slope in "
    "each pT interval for both data and MC; it is not a data-only "
    "calibration. New inputs may also contain a conditional data closure "
    "proxy with R2=Rn=1 and its MC truth-bias control. Flavor response "
    "uses the pure parallel barrel population; the transverse sideband is "
    "retained elsewhere as a pileup control. The fit uses "
    "cell-specific MC responses and a Gaussian prior centered at one. "
    "MC/template and "
    "truth-fraction uncertainties are not included in the reported "
    "conditional covariance. Separate diagnostics estimate MC template "
    "statistics, an approximate multinomial purity term, and the response "
    "sensitivity to independent 5% relative tagging-transition changes and "
    "independent 10% relative truth-flavor-fraction changes. These "
    "terms are not yet a complete systematic covariance.")
    .Write("response_inference_model",TObject::kOverwrite);
  output->Write("",TObject::kOverwrite);
  output->Purge();
  output->Close();

  const std::string taggingTable =
    std::string(outputDirectory)+"/tagging_matrices.tsv";
  std::ofstream taggingStream(taggingTable.c_str());
  taggingStream << std::setprecision(10)
    << "pt_min\ttrue_id\ttrue_flavor\treco_id\treco_flavor"
       "\ttruth_prior_mc\treco_fraction_data\tjoint_mc"
       "\tjoint_data_inferred\tefficiency_mc\tefficiency_data_inferred"
       "\ttransition_sf\tpurity_mc\tpurity_data_inferred\tpurity_ratio\n";
  for (int flavor=0; flavor!=nTruth; ++flavor) {
    const int trueId = truthGroups[flavor].outputId;
    for (int reco=0; reco!=nReco; ++reco)
      taggingStream
        << minimumPt << '\t' << trueId << '\t' << truthName(trueId) << '\t'
        << recoIds[reco] << '\t' << recoName(recoIds[reco]) << '\t'
        << truthPrior[flavor] << '\t' << dataRecoFraction[reco] << '\t'
        << jointMC(reco,flavor) << '\t' << inferredJoint(reco,flavor) << '\t'
        << efficiencyMC(reco,flavor) << '\t'
        << efficiencyData(reco,flavor) << '\t'
        << transitionSF(reco,flavor) << '\t'
        << compositionMC(reco,flavor) << '\t'
        << compositionData(reco,flavor) << '\t'
        << compositionRatio(reco,flavor) << '\n';
  }

  const std::string responseTable =
    std::string(outputDirectory)+"/response_residuals.tsv";
  std::ofstream responseStream(responseTable.c_str());
  responseStream << std::setprecision(10)
    << "pt_min\ttrue_id\ttrue_flavor\tdata_over_mc\tuncertainty"
       "\tmc_over_data\tcorrection_uncertainty\n";
  for (int flavor=0; flavor!=nTruth; ++flavor) {
    const int trueId = truthGroups[flavor].outputId;
    const double value = fit.residual[flavor];
    const double error = fit.covariance(flavor,flavor)>0.
      ? std::sqrt(fit.covariance(flavor,flavor)) : 0.;
    responseStream << minimumPt << '\t' << trueId << '\t'
      << truthName(trueId) << '\t' << value << '\t' << error << '\t'
      << (value!=0. ? 1./value : 0.) << '\t'
      << (value!=0. ? error/(value*value) : 0.) << '\n';
  }

  const std::string responsePtTable =
    std::string(outputDirectory)+"/response_residuals_vs_pt.tsv";
  std::ofstream responsePtStream(responsePtTable.c_str());
  responsePtStream << std::setprecision(10)
    << "pt_low\tpt_high\ttrue_id\ttrue_flavor\tdata_over_mc"
       "\tuncertainty\tfit_solved\trank\tcondition\n";
  for (size_t ptBin=0; ptBin<ptFits.size(); ++ptBin) {
    for (int flavor=0; flavor!=nTruth; ++flavor) {
      const int trueId = truthGroups[flavor].outputId;
      const FitResult &ptFit = ptFits[ptBin].responseFit;
      const double error = ptFit.covariance(flavor,flavor)>0.
        ? std::sqrt(ptFit.covariance(flavor,flavor)) : 0.;
      responsePtStream
        << ptFitEdges[ptBin] << '\t' << ptFitEdges[ptBin+1] << '\t'
        << trueId << '\t' << truthGroups[flavor].name << '\t'
        << ptFit.residual[flavor] << '\t' << error << '\t'
        << (ptFit.solved ? 1 : 0) << '\t'
        << ptFit.rank << '\t'
        << ptFit.nonzeroCondition << '\n';
    }
  }

  const std::string responseBinningTable =
    std::string(outputDirectory)+
      "/response_residuals_vs_binning_and_pt.tsv";
  std::ofstream responseBinningStream(responseBinningTable.c_str());
  responseBinningStream << std::setprecision(10)
    << "pt_binning\tpt_low\tpt_high\ttrue_id\ttrue_flavor\tdata_over_mc"
       "\tuncertainty\tfit_solved\trank\tcondition\n";
  for (const auto &variantEntry : variantPtFits)
    for (size_t ptBin=0; ptBin<variantEntry.second.size(); ++ptBin)
      for (int flavor=0; flavor!=nTruth; ++flavor) {
        const FitResult &ptFit = variantEntry.second[ptBin].responseFit;
        const double error = ptFit.covariance(flavor,flavor)>0.
          ? std::sqrt(ptFit.covariance(flavor,flavor)) : 0.;
        responseBinningStream
          << variantEntry.first << '\t'
          << ptFitEdges[ptBin] << '\t' << ptFitEdges[ptBin+1] << '\t'
          << truthGroups[flavor].outputId << '\t'
          << truthGroups[flavor].name << '\t'
          << ptFit.residual[flavor] << '\t' << error << '\t'
          << (ptFit.solved ? 1 : 0) << '\t' << ptFit.rank << '\t'
          << ptFit.nonzeroCondition << '\n';
      }

  const std::string ruImpactTable =
    std::string(outputDirectory)+"/response_ru_slope_impact.tsv";
  std::ofstream ruImpactStream(ruImpactTable.c_str());
  ruImpactStream << std::setprecision(10)
    << "pt_binning\tpt_low\tpt_high\tru_slope_mc\ttrue_flavor"
       "\tnominal_data_over_mc\tslope_data_over_mc\tshift_per_mille"
       "\tslope_fit_solved\trank\tcondition\n";
  for (const auto &variantEntry : variantPtFits)
    for (size_t ptBin=0; ptBin<variantEntry.second.size(); ++ptBin)
      for (int flavor=0; flavor!=nTruth; ++flavor) {
        const PtRangeResult &range = variantEntry.second[ptBin];
        const FitResult &alternate = range.responseFitRuSlope;
        ruImpactStream
          << variantEntry.first << '\t'
          << ptFitEdges[ptBin] << '\t' << ptFitEdges[ptBin+1] << '\t'
          << range.ruSlope << '\t' << truthGroups[flavor].name << '\t'
          << range.responseFit.residual[flavor] << '\t'
          << alternate.residual[flavor] << '\t'
          << 1000.*(alternate.residual[flavor]-
                    range.responseFit.residual[flavor]) << '\t'
          << (alternate.solved ? 1 : 0) << '\t' << alternate.rank << '\t'
          << alternate.nonzeroCondition << '\n';
      }

  const std::string tagResponseTable =
    std::string(outputDirectory)+"/tag_response_vs_pt.tsv";
  std::ofstream tagResponseStream(tagResponseTable.c_str());
  tagResponseStream << std::setprecision(10)
    << "pt_low\tpt_high\treco_id\treco_flavor"
       "\tdata_hdm\tdata_hdm_error\tmc_raw_hdm\tmc_raw_hdm_error"
       "\tdata_over_mc_raw\tdata_over_mc_raw_error"
       "\tmc_composition_corrected_hdm"
       "\tmc_composition_corrected_hdm_error"
       "\tdata_over_mc_composition_corrected"
       "\tdata_over_mc_composition_corrected_error\n";
  for (size_t ptBin=0; ptBin<ptFits.size(); ++ptBin)
    for (int reco=0; reco!=nReco; ++reco) {
      const PtRangeResult &range = ptFits[ptBin];
      tagResponseStream
        << ptFitEdges[ptBin] << '\t' << ptFitEdges[ptBin+1] << '\t'
        << recoIds[reco] << '\t' << recoName(recoIds[reco]) << '\t'
        << range.dataTagResponse[reco].mean << '\t'
        << range.dataTagResponse[reco].error << '\t'
        << range.mcRawTagResponse[reco].mean << '\t'
        << range.mcRawTagResponse[reco].error << '\t'
        << range.rawTagRatio[reco].mean << '\t'
        << range.rawTagRatio[reco].error << '\t'
        << range.mcCompositionCorrectedTagResponse[reco].mean << '\t'
        << range.mcCompositionCorrectedTagResponse[reco].error << '\t'
        << range.compositionCorrectedTagRatio[reco].mean << '\t'
        << range.compositionCorrectedTagRatio[reco].error << '\n';
    }

  const std::string recoilFractionTable =
    std::string(outputDirectory)+"/recoil_fractions_vs_pt.tsv";
  std::ofstream recoilFractionStream(recoilFractionTable.c_str());
  recoilFractionStream << std::setprecision(10)
    << "pt_binning\tpt_low\tpt_high\ttrue_id\ttrue_flavor\tcomponent"
       "\tmc_fraction\tmc_fraction_error\tdata_fraction_inferred"
       "\tdata_fraction_error\tdata_over_mc\tdata_over_mc_error"
       "\tfit_rank\tfit_condition\n";
  for (const auto &variantEntry : variantPtFits)
    for (size_t ptBin=0; ptBin<variantEntry.second.size(); ++ptBin)
      for (int flavor=0; flavor!=nTruth; ++flavor)
        for (const std::string component : {"fsr","ue","mnu","fnu","mnufsr"}) {
          const PtRangeResult &range = variantEntry.second[ptBin];
          const std::vector<ProfileSummary> &mcValues =
            componentValues(range,component,0);
          const std::vector<ProfileSummary> &dataValues =
            componentValues(range,component,1);
          const std::vector<ProfileSummary> &ratios =
            componentValues(range,component,2);
          const FitResult &fitComponent = componentFit(range,component);
          recoilFractionStream
            << variantEntry.first << '\t'
            << ptFitEdges[ptBin] << '\t' << ptFitEdges[ptBin+1] << '\t'
            << truthGroups[flavor].outputId << '\t'
            << truthGroups[flavor].name << '\t' << component << '\t'
            << mcValues[flavor].mean << '\t'
            << mcValues[flavor].error << '\t'
            << dataValues[flavor].mean << '\t'
            << dataValues[flavor].error << '\t'
            << ratios[flavor].mean << '\t' << ratios[flavor].error << '\t'
            << fitComponent.rank << '\t'
            << fitComponent.nonzeroCondition << '\n';
        }

  const std::string uncertaintyTable =
    std::string(outputDirectory)+"/response_uncertainties.tsv";
  std::ofstream uncertaintyStream(uncertaintyTable.c_str());
  uncertaintyStream << std::setprecision(10)
    << "true_id\ttrue_flavor\tdata_response_stat_conditional"
       "\tmc_purity_stat_approx\tmc_response_template_stat"
       "\tresponse_shift_for_independent_5pct_tagging_efficiencies"
       "\tresponse_shift_for_independent_10pct_flavor_fractions\n";
  for (int flavor=0; flavor!=nTruth; ++flavor)
    uncertaintyStream
      << truthGroups[flavor].outputId << '\t'
      << truthGroups[flavor].name << '\t'
      << std::sqrt(std::max(0.,fit.covariance(flavor,flavor))) << '\t'
      << purityStatUncertainty[flavor] << '\t'
      << responseTemplateUncertainty[flavor] << '\t'
      << taggingEfficiency5PctUncertainty[flavor] << '\t'
      << flavorFraction10PctSensitivity[flavor] << '\n';

  const std::string summaryName =
    std::string(outputDirectory)+"/fit_summary.json";
  std::ofstream summary(summaryName.c_str());
  summary << std::setprecision(12)
    << "{\n"
    << "  \"data_file\": \"" << jsonEscape(dataFileName) << "\",\n"
    << "  \"mc_file\": \"" << jsonEscape(mcFileName) << "\",\n"
    << "  \"minimum_pt\": " << minimumPt << ",\n"
    << "  \"tagging_model\": \"KL/IPF projection with fixed MC truth and "
       "observed data reco marginals\",\n"
    << "  \"tagging_count_source\": \"unsubtracted parallel window\",\n"
    << "  \"response_source\": \"pure parallel barrel component "
       "profiles\",\n"
    << "  \"transition_scale_factors_are_model_dependent\": true,\n"
    << "  \"ipf_converged\": " << (ipfConverged ? "true" : "false")
    << ",\n"
    << "  \"ipf_iterations\": " << ipfIterations << ",\n"
    << "  \"ipf_maximum_marginal_error\": " << ipfMaximumError << ",\n"
    << "  \"response_model\": \"one true-flavor data/MC multiplier common "
       "to all reconstructed tags\",\n"
    << "  \"response_prior_sigma\": " << responsePriorSigma << ",\n"
    << "  \"response_rows\": " << fit.rows << ",\n"
    << "  \"active_truth_flavors\": " << nTruth << ",\n"
    << "  \"response_rank\": " << fit.rank << ",\n"
    << "  \"response_fit_solved\": "
    << (fit.solved ? "true" : "false") << ",\n"
    << "  \"response_rank_deficient\": "
    << (fit.rank<nTruth ? "true" : "false") << ",\n"
    << "  \"response_nonzero_condition_number\": "
    << (std::isfinite(fit.nonzeroCondition)
          ? std::to_string(fit.nonzeroCondition) : "null") << ",\n"
    << "  \"response_full_condition_number\": "
    << (std::isfinite(fit.fullCondition)
          ? std::to_string(fit.fullCondition) : "null") << ",\n"
    << "  \"response_chi2\": " << fit.chi2 << ",\n"
    << "  \"hdm_constructed_from_component_means\": true,\n"
    << "  \"sparse_cell_response_fallbacks\": " << cellFallbacks << ",\n"
    << "  \"conditional_covariance_only\": true\n"
    << "}\n";

  std::cout << "Flavor-matrix analysis written to " << rootOutput << "\n"
            << "IPF: " << (ipfConverged ? "converged" : "NOT converged")
            << " after " << ipfIterations << " iterations; max marginal error "
            << ipfMaximumError << "\n"
            << "Response system: " << fit.rows << " rows, " << nTruth
            << " active true flavors, rank " << fit.rank
            << ", nonzero-mode condition number " << fit.nonzeroCondition
            << (fit.rank<nTruth
                  ? "; rank deficient, so the prior controls null modes" : "")
            << "\n"
            << "Model warning: per-transition data SFs are the minimum-KL "
               "inference, not independent measurements.\n";
}
