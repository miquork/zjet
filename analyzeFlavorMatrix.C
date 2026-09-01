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

const std::vector<int> recoIds = {0,1,4,5,6};
const std::vector<int> truthIds = {0,1,3,4,5,6};

const char *recoName(int id) {
  if (id==0) return "undefined";
  if (id==1) return "uds";
  if (id==4) return "c";
  if (id==5) return "b";
  if (id==6) return "g";
  return "unused";
}

const char *truthName(int id) {
  if (id==0) return "undefined";
  if (id==1) return "d+u";
  if (id==3) return "s";
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
  axis->SetBinLabel(1,"undefined");
  axis->SetBinLabel(2,"d+u");
  axis->SetBinLabel(3,"unused");
  axis->SetBinLabel(4,"s");
  axis->SetBinLabel(5,"c");
  axis->SetBinLabel(6,"b");
  axis->SetBinLabel(7,"g");
}

bool selectedPtBin(const TAxis *axis, int bin, double minimumPt) {
  return axis && bin>=1 && bin<=axis->GetNbins() &&
         axis->GetBinUpEdge(bin)>minimumPt+1.e-9;
}

double integratedCount(const TH3D *histogram, int recoId, int truthId,
                       double minimumPt) {
  if (!histogram) return 0.;
  double result = 0.;
  const int iy = histogram->GetYaxis()->FindFixBin(recoId);
  const int iz = histogram->GetZaxis()->FindFixBin(truthId);
  for (int ix=1; ix<=histogram->GetNbinsX(); ++ix) {
    if (!selectedPtBin(histogram->GetXaxis(),ix,minimumPt)) continue;
    result += histogram->GetBinContent(ix,iy,iz);
  }
  return result;
}

double integratedDataCount(const TH3D *histogram, int recoId,
                           double minimumPt) {
  if (!histogram) return 0.;
  double result = 0.;
  const int iy = histogram->GetYaxis()->FindFixBin(recoId);
  for (int ix=1; ix<=histogram->GetNbinsX(); ++ix) {
    if (!selectedPtBin(histogram->GetXaxis(),ix,minimumPt)) continue;
    for (int iz=1; iz<=histogram->GetNbinsZ(); ++iz)
      result += histogram->GetBinContent(ix,iy,iz);
  }
  return result;
}

ProfileSummary summarizeProfile(const TProfile3D *profile, int recoId,
                                int truthId, bool sumTruth,
                                double minimumPt) {
  ProfileSummary result;
  if (!profile) return result;
  const int iy = profile->GetYaxis()->FindFixBin(recoId);
  const int firstTruth = sumTruth ? 1 : profile->GetZaxis()->FindFixBin(truthId);
  const int lastTruth = sumTruth ? profile->GetNbinsZ() : firstTruth;
  double weightedSum = 0.;
  double numeratorVariance = 0.;
  for (int ix=1; ix<=profile->GetNbinsX(); ++ix) {
    if (!selectedPtBin(profile->GetXaxis(),ix,minimumPt)) continue;
    for (int iz=firstTruth; iz<=lastTruth; ++iz) {
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

ProfileSummary componentHDMFallback(TFile *file, int recoId, int truthId,
                                    bool sumTruth, double minimumPt) {
  ProfileSummary result;
  if (!file) return result;
  const char *names[] = {
    "FlavorMatrix/p3m0tc_flavormatrix",
    "FlavorMatrix/p3mntc_flavormatrix",
    "FlavorMatrix/p3mutc_flavormatrix",
  };
  ProfileSummary component[3];
  for (int index=0; index!=3; ++index) {
    TProfile3D *profile = dynamic_cast<TProfile3D*>(file->Get(names[index]));
    component[index] = summarizeProfile(
      profile,recoId,truthId,sumTruth,minimumPt);
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

ProfileSummary readHDM(TFile *file, int recoId, int truthId, bool sumTruth,
                       double minimumPt) {
  // HDM is a nonlinear ratio.  Averaging an event-wise HDM profile is not the
  // same as constructing HDM from the merged component means, and individual
  // events can have a nearly singular denominator.  Always use the mergeable
  // m0/mn/mu component profiles here.
  return componentHDMFallback(file,recoId,truthId,sumTruth,minimumPt);
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
  const int nTruthAll = truthIds.size();
  TMatrixD mcRaw(nReco,nTruthAll);
  TVectorD dataRaw(nReco);
  mcRaw.Zero();
  dataRaw.Zero();
  double totalMC = 0.;
  double totalData = 0.;
  for (int reco=0; reco!=nReco; ++reco) {
    const double dataCount =
      integratedDataCount(dataCounts,recoIds[reco],minimumPt);
    if (dataCount < -1.e-9)
      throw std::runtime_error(
        "Negative integrated data parallel count: the IPF likelihood "
        "requires a non-negative population");
    dataRaw[reco] = std::max(0.,dataCount);
    totalData += dataRaw[reco];
    for (int truth=0; truth!=nTruthAll; ++truth) {
      const double mcCount = integratedCount(
        mcCounts,recoIds[reco],truthIds[truth],minimumPt);
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

  std::vector<int> activeTruth;
  for (int truth=0; truth!=nTruthAll; ++truth) {
    double column = 0.;
    for (int reco=0; reco!=nReco; ++reco) column += mcRaw(reco,truth);
    if (column>1.e-12*totalMC) activeTruth.push_back(truth);
  }
  const int nTruth = activeTruth.size();
  if (nTruth==0) throw std::runtime_error("No active MC truth flavors");

  TVectorD truthPrior(nTruth);
  TVectorD dataRecoFraction(nReco);
  TMatrixD efficiencyMC(nReco,nTruth);
  TMatrixD jointMC(nReco,nTruth);
  for (int reco=0; reco!=nReco; ++reco)
    dataRecoFraction[reco] = dataRaw[reco]/totalData;
  double activeTotalMC = 0.;
  for (int index : activeTruth)
    for (int reco=0; reco!=nReco; ++reco)
      activeTotalMC += mcRaw(reco,index);
  for (int flavor=0; flavor!=nTruth; ++flavor) {
    const int truth = activeTruth[flavor];
    double column = 0.;
    for (int reco=0; reco!=nReco; ++reco) column += mcRaw(reco,truth);
    truthPrior[flavor] = column/activeTotalMC;
    for (int reco=0; reco!=nReco; ++reco) {
      efficiencyMC(reco,flavor) = column>0. ? mcRaw(reco,truth)/column : 0.;
      jointMC(reco,flavor) = mcRaw(reco,truth)/activeTotalMC;
    }
  }

  // Add a very small pseudocount only to avoid structural zeros in a sparse
  // smoke sample.  IPF restores the original MC truth marginals exactly.
  const double pseudocount = std::max(1.e-12,activeTotalMC*1.e-12);
  TMatrixD inferredJoint(nReco,nTruth);
  double initialTotal = 0.;
  for (int reco=0; reco!=nReco; ++reco)
    for (int flavor=0; flavor!=nTruth; ++flavor) {
      const int truth = activeTruth[flavor];
      inferredJoint(reco,flavor) =
        std::max(mcRaw(reco,truth),pseudocount);
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
  TMatrixD compositionData(nReco,nTruth);
  for (int reco=0; reco!=nReco; ++reco)
    for (int flavor=0; flavor!=nTruth; ++flavor) {
      efficiencyData(reco,flavor) = truthPrior[flavor]>0.
        ? inferredJoint(reco,flavor)/truthPrior[flavor] : 0.;
      transitionSF(reco,flavor) = efficiencyMC(reco,flavor)>0.
        ? efficiencyData(reco,flavor)/efficiencyMC(reco,flavor) : 0.;
      compositionData(reco,flavor) = dataRecoFraction[reco]>0.
        ? inferredJoint(reco,flavor)/dataRecoFraction[reco] : 0.;
    }

  std::vector<ProfileSummary> dataResponse(nReco);
  std::vector<std::vector<ProfileSummary> > mcResponse(
    nReco,std::vector<ProfileSummary>(nTruth));
  int cellFallbacks = 0;
  for (int reco=0; reco!=nReco; ++reco) {
    dataResponse[reco] = readHDM(
      dataFile.get(),recoIds[reco],0,true,minimumPt);
    for (int flavor=0; flavor!=nTruth; ++flavor) {
      mcResponse[reco][flavor] = readHDM(
        mcFile.get(),recoIds[reco],truthIds[activeTruth[flavor]],false,
        minimumPt);
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
    const int truth = activeTruth[flavor];
    for (int reco=0; reco!=nReco; ++reco) {
      if (!mcResponse[reco][flavor].valid) continue;
      const double count = mcRaw(reco,truth);
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
    const int trueId = truthIds[activeTruth[flavor]];
    const int bin = hTruthPrior->GetXaxis()->FindFixBin(trueId);
    hTruthPrior->SetBinContent(bin,truthPrior[flavor]);
  }
  for (int reco=0; reco!=nReco; ++reco)
    for (int flavor=0; flavor!=nTruth; ++flavor) {
      const int xbin = hEfficiencyMC->GetXaxis()->FindFixBin(recoIds[reco]);
      const int trueId = truthIds[activeTruth[flavor]];
      const int ybin = hEfficiencyMC->GetYaxis()->FindFixBin(trueId);
      hEfficiencyMC->SetBinContent(xbin,ybin,efficiencyMC(reco,flavor));
      hEfficiencyData->SetBinContent(xbin,ybin,efficiencyData(reco,flavor));
      hTransitionSF->SetBinContent(xbin,ybin,transitionSF(reco,flavor));
      hJointMC->SetBinContent(xbin,ybin,jointMC(reco,flavor));
      hJointData->SetBinContent(xbin,ybin,inferredJoint(reco,flavor));
      hCompositionData->SetBinContent(xbin,ybin,compositionData(reco,flavor));
    }

  response->cd();
  TH2D *hResponseMC = new TH2D(
    "h2_response_mc_by_transition",
    ";Reco hybrid flavor;True flavor;MC HDM response",7,idBins,7,idBins);
  labelRecoAxis(hResponseMC->GetXaxis());
  labelTruthAxis(hResponseMC->GetYaxis());
  for (int reco=0; reco!=nReco; ++reco)
    for (int flavor=0; flavor!=nTruth; ++flavor) {
      const int xbin = hResponseMC->GetXaxis()->FindFixBin(recoIds[reco]);
      const int trueId = truthIds[activeTruth[flavor]];
      const int ybin = hResponseMC->GetYaxis()->FindFixBin(trueId);
      hResponseMC->SetBinContent(xbin,ybin,mcResponse[reco][flavor].mean);
      hResponseMC->SetBinError(xbin,ybin,mcResponse[reco][flavor].error);
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
    const int trueId = truthIds[activeTruth[flavor]];
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
  TObjString(
    "Data transition efficiencies are the KL/IPF projection closest to the "
    "MC joint matrix with fixed MC truth marginals and observed data reco "
    "marginals. Tagging marginals use the un-subtracted parallel population "
    "to keep a non-negative likelihood, while response profiles use the "
    "signed parallel-minus-transverse estimator. Individual transition SFs "
    "are model dependent and are not independently identified by data tag "
    "fractions.")
    .Write("tagging_inference_model",TObject::kOverwrite);
  TObjString(
    "Response fit assumes one multiplicative data/MC residual per true "
    "flavor, common to all reconstructed tags. HDM is constructed after "
    "merging from the m0, mn and mu profile means with Rn=1 and Ru=0.92; "
    "the event-wise p3hdm profile is deliberately not used. The fit uses "
    "cell-specific MC responses and a Gaussian prior centered at one. "
    "MC/template and "
    "truth-fraction uncertainties are not yet included in the reported "
    "conditional covariance.")
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
       "\ttransition_sf\n";
  for (int flavor=0; flavor!=nTruth; ++flavor) {
    const int trueId = truthIds[activeTruth[flavor]];
    for (int reco=0; reco!=nReco; ++reco)
      taggingStream
        << minimumPt << '\t' << trueId << '\t' << truthName(trueId) << '\t'
        << recoIds[reco] << '\t' << recoName(recoIds[reco]) << '\t'
        << truthPrior[flavor] << '\t' << dataRecoFraction[reco] << '\t'
        << jointMC(reco,flavor) << '\t' << inferredJoint(reco,flavor) << '\t'
        << efficiencyMC(reco,flavor) << '\t'
        << efficiencyData(reco,flavor) << '\t'
        << transitionSF(reco,flavor) << '\n';
  }

  const std::string responseTable =
    std::string(outputDirectory)+"/response_residuals.tsv";
  std::ofstream responseStream(responseTable.c_str());
  responseStream << std::setprecision(10)
    << "pt_min\ttrue_id\ttrue_flavor\tdata_over_mc\tuncertainty"
       "\tmc_over_data\tcorrection_uncertainty\n";
  for (int flavor=0; flavor!=nTruth; ++flavor) {
    const int trueId = truthIds[activeTruth[flavor]];
    const double value = fit.residual[flavor];
    const double error = fit.covariance(flavor,flavor)>0.
      ? std::sqrt(fit.covariance(flavor,flavor)) : 0.;
    responseStream << minimumPt << '\t' << trueId << '\t'
      << truthName(trueId) << '\t' << value << '\t' << error << '\t'
      << (value!=0. ? 1./value : 0.) << '\t'
      << (value!=0. ? error/(value*value) : 0.) << '\n';
  }

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
    << "  \"response_source\": \"signed parallel-minus-transverse "
       "component profiles\",\n"
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
