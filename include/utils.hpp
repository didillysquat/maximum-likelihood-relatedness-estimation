#ifndef UTILS_H
#define UTILS_H

#include <string>
#include <vector>

#include <time.h>
#include <sys/resource.h>
#include <sys/time.h>

enum class InferenceType : uint8_t {
    BEST_GENOTYPE, ALL_GENOTYPES
};

enum class LikelihoodFormat : uint8_t {
    RAW, LOG, PHRED
};

struct StringItPair {
        StringItPair(std::string::iterator b, std::string::iterator e) :
            begin(b), end(e) {}
        StringItPair(StringItPair&& ) = default;
        StringItPair(const StringItPair& ) = default;
        std::string::iterator begin;
        std::string::iterator end;
};

std::vector<std::string> split(std::string &, char);
std::vector<StringItPair> split_it(std::string &, char);
std::vector<StringItPair> split_it(StringItPair&, char);
void print_time_elapsed(std::string, struct timeval*, struct timeval*);


#endif
