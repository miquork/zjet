#ifndef ZJET_JER_RESOLUTION_H
#define ZJET_JER_RESOLUTION_H

#include <algorithm>
#include <cmath>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace zjetcorrections {

// Minimal reader for the public JME PtResolution payload used by this
// analysis. The payload has JetEta and Rho binning, one clamped JetPt
// variable, and the formula
//   sqrt(p0*abs(p0)/pt^2 + p1^2*pt^p3 + p2^2).
// Evaluating it directly avoids constructing an additional globally named
// ROOT TFormula alongside the JEC correctors on bare HTCondor workers.
class JetPtResolution {
 public:
  explicit JetPtResolution(const std::string &fileName) {
    std::ifstream input(fileName);
    if (!input)
      throw std::runtime_error("Cannot open JER resolution file " + fileName);

    std::string line;
    while (std::getline(input,line)) {
      const std::string::size_type first = line.find_first_not_of(" \t\r");
      if (first==std::string::npos || line[first]=='#' || line[first]=='{')
        continue;
      std::istringstream values(line);
      Record record;
      unsigned parameterCount = 0;
      if (!(values >> record.etaMin >> record.etaMax
                   >> record.rhoMin >> record.rhoMax
                   >> parameterCount >> record.ptMin >> record.ptMax
                   >> record.p0 >> record.p1 >> record.p2 >> record.p3) ||
          parameterCount!=6) {
        throw std::runtime_error(
          "Invalid JER resolution record in " + fileName + ": " + line);
      }
      std::string trailing;
      if (values >> trailing)
        throw std::runtime_error(
          "Unexpected trailing JER resolution field in " + fileName);
      records_.push_back(record);
    }
    if (records_.empty())
      throw std::runtime_error("No JER resolution records in " + fileName);
  }

  double resolution(double eta, double rho, double pt) const {
    for (const Record &record : records_) {
      // Match the JME reference reader: ranges include both endpoints and the
      // first matching record wins at a shared boundary.
      if (eta<record.etaMin || eta>record.etaMax ||
          rho<record.rhoMin || rho>record.rhoMax)
        continue;
      const double x = std::max(record.ptMin,std::min(pt,record.ptMax));
      const double variance =
        record.p0*std::fabs(record.p0)/(x*x) +
        record.p1*record.p1*std::pow(x,record.p3) +
        record.p2*record.p2;
      return std::sqrt(std::max(0.,variance));
    }
    return 1.;
  }

  size_t size() const { return records_.size(); }

 private:
  struct Record {
    double etaMin = 0.;
    double etaMax = 0.;
    double rhoMin = 0.;
    double rhoMax = 0.;
    double ptMin = 0.;
    double ptMax = 0.;
    double p0 = 0.;
    double p1 = 0.;
    double p2 = 0.;
    double p3 = 0.;
  };

  std::vector<Record> records_;
};

// The companion scale-factor payload is also a fixed JetEta-binned formula.
// Keep it independent of ROOT formula state for the same worker robustness.
class JetResolutionScaleFactor {
 public:
  explicit JetResolutionScaleFactor(const std::string &fileName) {
    std::ifstream input(fileName);
    if (!input)
      throw std::runtime_error("Cannot open JER scale-factor file " + fileName);

    std::string line;
    while (std::getline(input,line)) {
      const std::string::size_type first = line.find_first_not_of(" \t\r");
      if (first==std::string::npos || line[first]=='#' || line[first]=='{')
        continue;
      std::istringstream values(line);
      Record record;
      unsigned parameterCount = 0;
      if (!(values >> record.etaMin >> record.etaMax >> parameterCount
                   >> record.ptMin >> record.ptMax
                   >> record.p0 >> record.p1 >> record.p2
                   >> record.p3 >> record.p4 >> record.p5) ||
          parameterCount!=8) {
        throw std::runtime_error(
          "Invalid JER scale-factor record in " + fileName + ": " + line);
      }
      std::string trailing;
      if (values >> trailing)
        throw std::runtime_error(
          "Unexpected trailing JER scale-factor field in " + fileName);
      records_.push_back(record);
    }
    if (records_.empty())
      throw std::runtime_error("No JER scale-factor records in " + fileName);
  }

  double scaleFactor(double eta, double pt) const {
    for (const Record &record : records_) {
      if (eta<record.etaMin || eta>record.etaMax) continue;
      const double x = std::max(record.ptMin,std::min(pt,record.ptMax));
      const double numerator =
        record.p0*std::fabs(record.p0)/(x*x) +
        record.p1*record.p1/x + record.p2*record.p2;
      const double denominator =
        record.p3*std::fabs(record.p3)/(x*x) +
        record.p4*record.p4/x + record.p5*record.p5;
      return denominator>0. ?
        std::sqrt(std::max(0.,numerator)/denominator) : 1.;
    }
    return 1.;
  }

  size_t size() const { return records_.size(); }

 private:
  struct Record {
    double etaMin = 0.;
    double etaMax = 0.;
    double ptMin = 0.;
    double ptMax = 0.;
    double p0 = 0.;
    double p1 = 0.;
    double p2 = 0.;
    double p3 = 0.;
    double p4 = 0.;
    double p5 = 0.;
  };

  std::vector<Record> records_;
};

} // namespace zjetcorrections

#endif
