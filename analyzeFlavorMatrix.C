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
#include "TProfile3D.h"
#include "TSystem.h"
#include "TVectorD.h"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
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
};

struct TruthGroup {
  int outputId;
  const char *name;
  std::vector<int> sourceIds;
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

ProfileSummary componentHDMFallback(
  TFile *file, int recoId, const std::vector<int> &truthIds,
  double minimumPt, double maximumPt) {
  ProfileSummary result;
  if (!file) return result;
  const char *names[] = {
    "FlavorMatrix/p3m0tc_parallel_flavormatrix",
    "FlavorMatrix/p3mntc_parallel_flavormatrix",
    "FlavorMatrix/p3mutc_parallel_flavormatrix",
  };
  ProfileSummary component[3];
  for (int index=0; index!=3; ++index) {
    TProfile3D *profile = dynamic_cast<TProfile3D*>(file->Get(names[index]));
    component[index] = summarizeProfile(
      profile,recoId,truthIds,minimumPt,maximumPt);
    if (!component[index].valid) return result;
  }
  const double m0 = component[0].mean;
  const double mn = component[1].mean;
  const double mu = component[2].mean;
  const double denominator =
    1.-mn/ZJetFlavorMatrix::responseN-mu/ZJetFlavorMatrix::responseU;
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
    (-denominator+numerator/ZJetFlavorMatrix::responseU)/
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
                       double minimumPt, double maximumPt) {
  // HDM is a nonlinear ratio.  Averaging an event-wise HDM profile is not the
  // same as constructing HDM from the merged component means, and individual
  // events can have a nearly singular denominator.  Always use the mergeable
  // m0/mn/mu component profiles here.
  return componentHDMFallback(
    file,recoId,truthIds,minimumPt,maximumPt);
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
  const double maximumPt = std::numeric_limits<double>::infinity();
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
  TVectorD flavorFractionSensitivity(nTruth);
  responseTemplateUncertainty.Zero();
  purityStatUncertainty.Zero();
  flavorFractionSensitivity.Zero();

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

  for (int variedFlavor=0; variedFlavor!=nTruth; ++variedFlavor) {
    TMatrixD variedDesign = design;
    for (int row=0; row!=static_cast<int>(fitReco.size()); ++row) {
      const int reco = fitReco[row];
      std::vector<double> varied(nTruth);
      double normalization = 0.;
      for (int flavor=0; flavor!=nTruth; ++flavor) {
        varied[flavor] = compositionData(reco,flavor)*
          (flavor==variedFlavor ? 1.01 : 1.00);
        normalization += varied[flavor];
      }
      for (int flavor=0; flavor!=nTruth; ++flavor)
        variedDesign(row,flavor) = varied[flavor]/normalization*
          mcResponse[reco][flavor].mean;
    }
    FitResult variedFit = fitResponse(
      variedDesign,observed,observedError,responsePriorSigma);
    if (!variedFit.solved) continue;
    for (int resultFlavor=0; resultFlavor!=nTruth; ++resultFlavor) {
      const double shift = variedFit.residual[resultFlavor]-
                           fit.residual[resultFlavor];
      flavorFractionSensitivity[resultFlavor] += shift*shift;
    }
  }
  for (int flavor=0; flavor!=nTruth; ++flavor)
    flavorFractionSensitivity[flavor] =
      std::sqrt(flavorFractionSensitivity[flavor]);

  auto fitPtRange = [&](double lowPt, double highPt) {
    FitResult rangeResult(nTruth);
    TMatrixD rangeMC(nReco,nTruth);
    TVectorD rangeData(nReco);
    rangeMC.Zero();
    rangeData.Zero();
    double rangeTotalMC = 0.;
    double rangeTotalData = 0.;
    for (int reco=0; reco!=nReco; ++reco) {
      rangeData[reco] = std::max(0.,integratedDataCount(
        dataCounts,recoIds[reco],lowPt,highPt));
      rangeTotalData += rangeData[reco];
      for (int flavor=0; flavor!=nTruth; ++flavor) {
        rangeMC(reco,flavor) = std::max(0.,integratedCount(
          mcCounts,recoIds[reco],truthGroups[flavor].sourceIds,
          lowPt,highPt));
        rangeTotalMC += rangeMC(reco,flavor);
      }
    }
    if (!(rangeTotalMC>0.) || !(rangeTotalData>0.)) return rangeResult;

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
    for (int reco=0; reco!=nReco; ++reco)
      for (int flavor=0; flavor!=nTruth; ++flavor)
        rangeComposition(reco,flavor) = rangeRecoFraction[reco]>0.
          ? rangeJoint(reco,flavor)/rangeRecoFraction[reco] : 0.;

    std::vector<ProfileSummary> rangeDataResponse(nReco);
    std::vector<std::vector<ProfileSummary> > rangeMCResponse(
      nReco,std::vector<ProfileSummary>(nTruth));
    for (int reco=0; reco!=nReco; ++reco) {
      rangeDataResponse[reco] = readHDM(
        dataFile.get(),recoIds[reco],std::vector<int>{0},lowPt,highPt);
      for (int flavor=0; flavor!=nTruth; ++flavor)
        rangeMCResponse[reco][flavor] = readHDM(
          mcFile.get(),recoIds[reco],truthGroups[flavor].sourceIds,
          lowPt,highPt);
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
    return fitResponse(
      rangeDesign,rangeObserved,rangeErrors,responsePriorSigma);
  };

  const std::vector<double> ptFitEdges = {
    30.,40.,60.,85.,125.,180.,250.,400.,
  };
  std::vector<FitResult> ptFits;
  for (size_t bin=0; bin+1<ptFitEdges.size(); ++bin)
    ptFits.push_back(fitPtRange(ptFitEdges[bin],ptFitEdges[bin+1]));

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
  for (int flavor=0; flavor!=nTruth; ++flavor) {
    TGraphErrors *graph = new TGraphErrors();
    graph->SetName(Form("g_response_residual_vs_pt_%s",
                        truthGroups[flavor].name));
    graph->SetTitle(Form(
      ";p_{T,Z} (GeV);R^{data}_{%s}/R^{MC}_{%s}",
      truthGroups[flavor].name,truthGroups[flavor].name));
    for (size_t bin=0; bin<ptFits.size(); ++bin) {
      // Keep the graph quantitative: rank-deficient or nearly singular bins
      // remain in the TSV diagnostics but are not drawn as measurements.
      if (!ptFits[bin].solved || ptFits[bin].rank<nTruth ||
          !std::isfinite(ptFits[bin].nonzeroCondition) ||
          ptFits[bin].nonzeroCondition>100.)
        continue;
      const double low = ptFitEdges[bin];
      const double high = ptFitEdges[bin+1];
      const double x = std::sqrt(low*high);
      const double error = ptFits[bin].covariance(flavor,flavor)>0.
        ? std::sqrt(ptFits[bin].covariance(flavor,flavor)) : 0.;
      const int point = graph->GetN();
      graph->SetPoint(point,x,ptFits[bin].residual[flavor]);
      graph->SetPointError(point,0.5*(high-low),error);
    }
    graph->Write();
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
  TParameter<int>("hdm_constructed_from_component_means",1).Write();
  TParameter<int>("sparse_cell_response_fallbacks",cellFallbacks).Write();
  TParameter<int>("zero_profile_error_fallbacks",fit.zeroErrorFallbacks).Write();
  TH2D *hUncertainty = new TH2D(
    "h2_response_uncertainty_components",
    ";True flavor;Source;Absolute uncertainty",
    nTruth,-0.5,nTruth-0.5,4,-0.5,3.5);
  const char *uncertaintySources[] = {
    "data stat.", "purity stat.",
    "response stat.", "fractions 1%",
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
      flavor+1,4,flavorFractionSensitivity[flavor]);
  }
  for (int source=0; source!=4; ++source)
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
    "the event-wise p3hdm profile is deliberately not used. Flavor response "
    "uses the pure parallel barrel population; the transverse sideband is "
    "retained elsewhere as a pileup control. The fit uses "
    "cell-specific MC responses and a Gaussian prior centered at one. "
    "MC/template and "
    "truth-fraction uncertainties are not included in the reported "
    "conditional covariance. Separate diagnostics estimate MC template "
    "statistics, an approximate multinomial purity term, and the response "
    "sensitivity to independent 1% relative flavor-fraction changes. These "
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
      const double error = ptFits[ptBin].covariance(flavor,flavor)>0.
        ? std::sqrt(ptFits[ptBin].covariance(flavor,flavor)) : 0.;
      responsePtStream
        << ptFitEdges[ptBin] << '\t' << ptFitEdges[ptBin+1] << '\t'
        << trueId << '\t' << truthGroups[flavor].name << '\t'
        << ptFits[ptBin].residual[flavor] << '\t' << error << '\t'
        << (ptFits[ptBin].solved ? 1 : 0) << '\t'
        << ptFits[ptBin].rank << '\t'
        << ptFits[ptBin].nonzeroCondition << '\n';
    }
  }

  const std::string uncertaintyTable =
    std::string(outputDirectory)+"/response_uncertainties.tsv";
  std::ofstream uncertaintyStream(uncertaintyTable.c_str());
  uncertaintyStream << std::setprecision(10)
    << "true_id\ttrue_flavor\tdata_response_stat_conditional"
       "\tmc_purity_stat_approx\tmc_response_template_stat"
       "\tresponse_shift_for_independent_1pct_flavor_fractions\n";
  for (int flavor=0; flavor!=nTruth; ++flavor)
    uncertaintyStream
      << truthGroups[flavor].outputId << '\t'
      << truthGroups[flavor].name << '\t'
      << std::sqrt(std::max(0.,fit.covariance(flavor,flavor))) << '\t'
      << purityStatUncertainty[flavor] << '\t'
      << responseTemplateUncertainty[flavor] << '\t'
      << flavorFractionSensitivity[flavor] << '\n';

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
