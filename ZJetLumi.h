#ifndef ZJET_LUMI_H
#define ZJET_LUMI_H

#include <cstdint>
#include <fstream>
#include <iostream>
#include <map>
#include <regex>
#include <sstream>
#include <string>
#include <vector>

class ZJetLumiData {
public:
  bool loadGoldenJson(const std::string &fileName) {
    std::ifstream input(fileName.c_str());
    if (!input.is_open()) return false;

    const std::string text((std::istreambuf_iterator<char>(input)),
                           std::istreambuf_iterator<char>());
    std::size_t position = 0;
    unsigned int runs = 0;
    unsigned long long lumisections = 0;
    const std::regex rangePattern("\\[\\s*([0-9]+)\\s*,\\s*([0-9]+)\\s*\\]");

    while ((position = text.find('"', position)) != std::string::npos) {
      const std::size_t quoteEnd = text.find('"', position + 1);
      if (quoteEnd == std::string::npos) break;
      const std::string runText = text.substr(position + 1, quoteEnd-position-1);
      if (runText.empty() || runText.find_first_not_of("0123456789") != std::string::npos) {
        position = quoteEnd + 1;
        continue;
      }

      const std::size_t arrayBegin = text.find('[', quoteEnd);
      if (arrayBegin == std::string::npos) break;
      int depth = 0;
      std::size_t arrayEnd = arrayBegin;
      for (; arrayEnd < text.size(); ++arrayEnd) {
        if (text[arrayEnd] == '[') ++depth;
        if (text[arrayEnd] == ']' && --depth == 0) break;
      }
      if (arrayEnd == text.size()) return false;

      const unsigned int run = static_cast<unsigned int>(std::stoul(runText));
      const std::string ranges = text.substr(arrayBegin, arrayEnd-arrayBegin+1);
      for (std::sregex_iterator it(ranges.begin(), ranges.end(), rangePattern), end;
           it != end; ++it) {
        const unsigned int first = static_cast<unsigned int>(std::stoul((*it)[1]));
        const unsigned int last = static_cast<unsigned int>(std::stoul((*it)[2]));
        if (last < first) return false;
        goodLumiRanges_[run].push_back(std::make_pair(first, last));
        lumisections += last-first+1;
      }
      ++runs;
      position = arrayEnd + 1;
    }

    std::cout << "Loaded " << runs << " runs and " << lumisections
              << " certified lumisections from " << fileName << std::endl;
    return !goodLumiRanges_.empty();
  }

  bool loadPileup(const std::string &fileName) {
    std::ifstream input(fileName.c_str());
    if (!input.is_open()) return false;

    std::string line;
    unsigned long long entries = 0;
    int avgpuColumn = -1;
    while (std::getline(input, line)) {
      if (line.empty()) continue;

      if (line[0] == '#') {
        std::stringstream header(line.substr(1));
        std::string column;
        int index = 0;
        while (std::getline(header,column,',')) {
          if (column == "avgpu") avgpuColumn = index;
          ++index;
        }
        continue;
      }

      unsigned int run = 0;
      unsigned int ls = 0;
      double mu = -1.;
      std::vector<std::string> fields;
      std::stringstream csv(line);
      std::string field;
      while (std::getline(csv, field, ',')) fields.push_back(field);

      // brilcalc CSV: #run:fill,ls,time,...,avgpu,source
      if (fields.size() >= 3) {
        try {
          run = static_cast<unsigned int>(std::stoul(fields[0]));
          ls = static_cast<unsigned int>(std::stoul(fields[1]));
          const int muIndex = (avgpuColumn >= 0 ? avgpuColumn
                                                : int(fields.size())-2);
          if (muIndex >= 0 && muIndex < int(fields.size()))
            mu = std::stod(fields[muIndex]);
        }
        catch (...) { mu = -1.; }
      }
      else {
        std::stringstream ascii(line);
        ascii >> run >> ls >> mu;
      }

      if (run > 0 && ls > 0 && mu >= 0.) {
        pileup_[key(run, ls)] = mu;
        ++entries;
      }
    }

    std::cout << "Loaded pileup mu for " << entries
              << " lumisections from " << fileName << std::endl;
    return !pileup_.empty();
  }

  bool accept(unsigned int run, unsigned int ls) const {
    if (goodLumiRanges_.empty()) return true;
    const auto runIt = goodLumiRanges_.find(run);
    if (runIt == goodLumiRanges_.end()) return false;
    for (const auto &range : runIt->second)
      if (ls >= range.first && ls <= range.second) return true;
    return false;
  }

  double pileup(unsigned int run, unsigned int ls) const {
    const auto it = pileup_.find(key(run, ls));
    return (it == pileup_.end() ? -1. : it->second);
  }

private:
  static std::uint64_t key(unsigned int run, unsigned int ls) {
    return (static_cast<std::uint64_t>(run) << 32) | ls;
  }

  std::map<unsigned int,
           std::vector<std::pair<unsigned int, unsigned int>>> goodLumiRanges_;
  std::map<std::uint64_t, double> pileup_;
};

#endif
