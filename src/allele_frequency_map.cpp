#include <fstream>
#include <iostream>
#include <memory>
#include <cstdlib>

#include "allele_frequency_map.hpp"
#include "utils.hpp"

AlleleFrequencyMap::AlleleFrequencyMap() {}

bool AlleleFrequencyMap::populateFromFile(const std::string& fname) {
  
  auto closeFile = [](std::ifstream* f) { f->close(); };
  // Open the file
  std::unique_ptr<std::ifstream, decltype(closeFile)> ifile(
      new std::ifstream(fname), closeFile);

  if (!ifile->is_open()) {
    std::cerr << "[Error] : Could not open allele frequency file " << fname << '\n';
    return false;
  }

  std::string line;
  while (std::getline(*ifile, line)) {
    auto toks = split(line, '\t');
    auto chr = toks[0];
    uint64_t pos = std::stoul(toks[1]);
    std::string ref = toks[2];
    std::string alt = toks[3];
    double   altFreq = std::stod(toks[4]);
    afMap_[std::make_tuple(chr, pos)] = {ref, alt, altFreq}; 
  }
  return true;
}

bool AlleleFrequencyMap::get(const std::string& chr, uint64_t pos, AlleleFrequencyValue& v) {
  auto afMapIt = afMap_.find(std::make_tuple(chr, pos));
  if (afMapIt != afMap_.end()) {
    v = afMapIt->second;
    return true;
  }
  return false;
}
