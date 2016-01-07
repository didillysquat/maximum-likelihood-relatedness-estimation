#ifndef ALLELE_FREQUENCY_MAP_HPP
#define ALLELE_FREQUENCY_MAP_HPP

#include <map>

struct AlleleFrequencyValue {
  //std::string chr;
  //uint64_t    pos;
  std::string ref;
  std::string     alt;
  double      altFreq;
};

class AlleleFrequencyMap {
  public:
    AlleleFrequencyMap();

    bool populateFromFile(const std::string& fname);
    bool get(const std::string& chr, uint64_t pos, AlleleFrequencyValue& v);
  private:
    std::map<std::tuple<std::string, uint64_t>, AlleleFrequencyValue> afMap_;
};

#endif // ALLELE_FREQUENCY_MAP_HPP
