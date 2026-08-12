#ifndef ZJET_MUON_CORRECTIONS_H
#define ZJET_MUON_CORRECTIONS_H

#include "data/MuonCorrections/2024_Summer24_generated.h"

#include "TRandom3.h"
#include "Math/QuantFuncMathCore.h"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <initializer_list>
#include <iterator>
#include <vector>

namespace zjetcorrections {

class SeedSequence {
 public:
  explicit SeedSequence(std::initializer_list<uint32_t> seeds)
    : seeds_(seeds) {}

  template <typename Iterator>
  void generate(Iterator begin, Iterator end) const {
    const size_t size = std::distance(begin,end);
    if (size==0) return;
    const uint32_t multiplier = 0x9e3779b9;
    const uint32_t mix = 0x85ebca6b;
    std::vector<uint32_t> buffer(size,0x8b8b8b8b);
    size_t index = 0;
    for (; index<std::min(size,seeds_.size()); ++index)
      buffer[index] ^= seeds_[index]+multiplier*index;
    for (; index<size; ++index) buffer[index] ^= multiplier*index;
    for (size_t item=0; item<size; ++item) {
      const uint32_t value =
        buffer[(item+size-1)%size]^(buffer[item]>>27);
      buffer[item] = (value*mix)^(buffer[item]<<13);
    }
    std::copy(buffer.begin(),buffer.end(),begin);
  }

 private:
  std::vector<uint32_t> seeds_;
};

class CrystalBall {
 public:
  CrystalBall(double mean, double sigma, double alpha, double exponent)
    : mean_(mean), sigma_(sigma), alpha_(alpha), exponent_(exponent) {
    initialize();
  }

  double inverseCdf(double probability) const {
    if (probability<cdfMinusAlpha_)
      return mean_+g_*(f_-std::pow(nc_/probability,k_));
    if (probability>cdfPlusAlpha_)
      return mean_-g_*(f_-std::pow(c_-probability/nc_,-k_));
    const double erfArgument =
      (d_-probability/ns_)/sqrtPiOver2_;
    // erf^-1(x) = Phi^-1((x+1)/2)/sqrt(2). ROOT's MathCore normal
    // quantile keeps the worker build independent of a system Boost install.
    return mean_-sigma_*ROOT::Math::normal_quantile(
      0.5*(erfArgument+1.),1.);
  }

 private:
  double cdf(double value) const {
    const double normalized = (value-mean_)/sigma_;
    if (normalized<-alpha_)
      return nc_/std::pow(f_-sigma_*normalized/g_,exponent_-1.);
    if (normalized>alpha_)
      return nc_*(c_-std::pow(f_+sigma_*normalized/g_,1.-exponent_));
    return ns_*(d_-sqrtPiOver2_*std::erf(-normalized/sqrt2_));
  }

  void initialize() {
    const double absoluteAlpha = std::fabs(alpha_);
    const double exponential = std::exp(-absoluteAlpha*absoluteAlpha/2.);
    const double tail = exponent_/absoluteAlpha/(exponent_-1.)*exponential;
    const double gaussian =
      2.*sqrtPiOver2_*std::erf(absoluteAlpha/sqrt2_);
    c_ = (gaussian+2.*tail)/tail;
    d_ = (gaussian+2.*tail)/2.;
    const double normalization = 1./sigma_/(gaussian+2.*tail);
    k_ = 1./(exponent_-1.);
    ns_ = normalization*sigma_;
    nc_ = ns_*tail;
    f_ = 1.-absoluteAlpha*absoluteAlpha/exponent_;
    g_ = sigma_*exponent_/absoluteAlpha;
    cdfMinusAlpha_ = cdf(mean_-alpha_*sigma_);
    cdfPlusAlpha_ = cdf(mean_+alpha_*sigma_);
  }

  static constexpr double pi_ = 3.14159;
  const double sqrtPiOver2_ = std::sqrt(pi_/2.);
  const double sqrt2_ = std::sqrt(2.);
  double mean_;
  double sigma_;
  double alpha_;
  double exponent_;
  double c_ = 0.;
  double d_ = 0.;
  double k_ = 0.;
  double ns_ = 0.;
  double nc_ = 0.;
  double f_ = 0.;
  double g_ = 0.;
  double cdfMinusAlpha_ = 0.;
  double cdfPlusAlpha_ = 0.;
};

template <size_t Size>
size_t clampedBin(double value, const double (&edges)[Size]) {
  if (value<edges[0]) return 0;
  if (value>=edges[Size-1]) return Size-2;
  return std::upper_bound(edges,edges+Size,value)-edges-1;
}

class Summer24MuonCorrections {
 public:
  double scaleFactor(bool isData, double pt, double eta, double phi,
                     int charge) const {
    if (pt==0. || pt<26.) return 1.;
    const size_t etaBin = clampedBin(
      eta,ZJetMuonCorrectionData::etaEdges);
    const size_t phiBin = clampedBin(
      phi,ZJetMuonCorrectionData::phiEdges);
    const double m = isData
      ? ZJetMuonCorrectionData::mData[etaBin][phiBin]
      : ZJetMuonCorrectionData::mMc[etaBin][phiBin];
    const double a = isData
      ? ZJetMuonCorrectionData::aData[etaBin][phiBin]
      : ZJetMuonCorrectionData::aMc[etaBin][phiBin];
    return 1./(m+charge*a*pt);
  }

  double resolutionFactor(double pt, double eta, double phi,
                          double trackerLayers, int eventNumber,
                          int luminosityBlock) const {
    if (pt==0. || pt<26. || pt>200.) return 1.;
    const size_t etaBin = clampedBin(
      std::fabs(eta),ZJetMuonCorrectionData::absoluteEtaEdges);
    const size_t layerBin = clampedBin(
      trackerLayers,ZJetMuonCorrectionData::trackerLayerEdges);
    const double *cb = ZJetMuonCorrectionData::crystalBall[etaBin][layerBin];
    CrystalBall distribution(cb[0],cb[1],cb[3],cb[2]);
    const int64_t phiSeed = static_cast<int64_t>(
      (phi/std::acos(-1.))*((1LL<<31)-1))&0xFFF;
    SeedSequence sequence({static_cast<uint32_t>(eventNumber),
                           static_cast<uint32_t>(luminosityBlock),
                           static_cast<uint32_t>(phiSeed)});
    uint32_t seed = 0;
    sequence.generate(&seed,&seed+1);
    TRandom3 random(seed);
    const double randomValue = distribution.inverseCdf(random.Rndm());

    const double *polynomial =
      ZJetMuonCorrectionData::resolutionPolynomial[etaBin][layerBin];
    const double standardDeviation = std::max(
      0.,polynomial[0]+polynomial[1]*pt+polynomial[2]*pt*pt);
    const double kData = ZJetMuonCorrectionData::kData[etaBin];
    const double kMc = ZJetMuonCorrectionData::kMc[etaBin];
    const double residual = kMc<kData
      ? std::sqrt(kData*kData-kMc*kMc) : 0.;
    const double correctedPt =
      pt*(1.+residual*standardDeviation*randomValue);
    const double factor = correctedPt/pt;
    if (!std::isfinite(factor) || factor>2. || factor<0.1 || correctedPt<0.)
      return 1.;
    return factor;
  }
};

} // namespace zjetcorrections

#endif
