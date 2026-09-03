#ifndef ZJET_FLAVOR_MATRIX_TOOLS_H
#define ZJET_FLAVOR_MATRIX_TOOLS_H

#include "TProfile3D.h"

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace ZJetFlavorMatrix {

constexpr double response2 = 1.00;
constexpr double responseN = 1.00;
constexpr double responseU = 0.92;

// HDM is nonlinear, so <HDM(event)> is not HDM(<m0>,<mn>,<mu>) and cannot be
// merged correctly across worker files. Rebuild the derived profile from the
// already-merged raw component means. The uncertainty below propagates the
// three profile errors without their (important) covariance; the raw profiles
// remain authoritative for precision uncertainty work.
inline void finalizeHDMProfile(TProfile3D *hdm, const TProfile3D *m0,
                               const TProfile3D *mn,
                               const TProfile3D *mu) {
  if (!hdm || !m0 || !mn || !mu)
    throw std::runtime_error("Missing FlavorMatrix HDM component profile");
  if (hdm->GetNcells()!=m0->GetNcells() ||
      hdm->GetNcells()!=mn->GetNcells() ||
      hdm->GetNcells()!=mu->GetNcells())
    throw std::runtime_error("Incompatible FlavorMatrix HDM profile shapes");

  hdm->Reset("ICES");
  if (hdm->GetSumw2N()==0 || hdm->GetBinSumw2()->GetSize()==0)
    hdm->Sumw2();
  for (int bin=0; bin<hdm->GetNcells(); ++bin) {
    const double sumw = m0->GetBinEntries(bin);
    if (!(sumw>0.)) continue;
    const double value0 = m0->GetBinContent(bin);
    const double valueN = mn->GetBinContent(bin);
    const double valueU = mu->GetBinContent(bin);
    const double numerator = value0-valueN-valueU;
    const double denominator =
      1.-valueN/responseN-valueU/responseU;
    if (!std::isfinite(numerator) || !std::isfinite(denominator) ||
        std::fabs(denominator)<1.e-8)
      continue;
    const double value = numerator/denominator;
    if (!std::isfinite(value)) continue;

    const double d0 = 1./denominator;
    const double dN =
      (-denominator+numerator/responseN)/(denominator*denominator);
    const double dU =
      (-denominator+numerator/responseU)/(denominator*denominator);
    const double error2 =
      std::pow(d0*m0->GetBinError(bin),2) +
      std::pow(dN*mn->GetBinError(bin),2) +
      std::pow(dU*mu->GetBinError(bin),2);

    double sumw2 = 0.;
    if (const TArrayD *weights2 = m0->GetBinSumw2())
      if (bin<weights2->GetSize()) sumw2 = weights2->At(bin);
    if (!(sumw2>0.)) sumw2 = std::fabs(sumw);
    const double effectiveEntries = sumw*sumw/sumw2;
    if (!(effectiveEntries>0.)) continue;

    hdm->SetBinEntries(bin,sumw);
    hdm->GetBinSumw2()->SetAt(sumw2,bin);
    hdm->SetBinContent(bin,sumw*value);
    // TProfile stores sum(w*y^2), not the displayed error itself.
    hdm->GetSumw2()->SetAt(
      sumw*(value*value+std::max(0.,error2)*effectiveEntries),bin);
  }
  hdm->SetEntries(m0->GetEntries());
}

} // namespace ZJetFlavorMatrix

#endif
