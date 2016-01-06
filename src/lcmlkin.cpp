#include "relatedness.hpp"
#include "allele_frequency_map.hpp"

#include "tclap/CmdLine.h"
#include "spdlog/spdlog.h"

#include <iostream>
#include <fstream>
#include <cstdlib>
#include <ctime>
#include <memory>

int main(int argc, char* argv[]){

    std::string versionStr = "v0.5.0";
    TCLAP::CmdLine cmd(
            "lcMLkin",
            ' ',
            versionStr);

    auto defaultNumThreads = std::thread::hardware_concurrency();
    TCLAP::ValueArg<std::string> input("i", "input", "The input VCF file",
                                       true, "", "path");
    TCLAP::ValueArg<std::string> output("o", "output", "Where output should be written",
                                       true, "", "path");
    TCLAP::ValueArg<uint32_t> numThreads("t", "threads", "Number of threads to use",
                                       false, defaultNumThreads, "integer >= 1");
    TCLAP::ValueArg<std::string> likelihoodFormat("l", "likelihoodFormat",
                                       "Type of genotype likelihood that should be used (raw | log | phred)",
                                       false, "raw", "string");
    TCLAP::ValueArg<std::string> type("g", "genotype",
                                      "Which inference algorithm to use; (all | best)",
                                      true, "all", "string");
    TCLAP::ValueArg<uint32_t> numBootstraps("b", "numBootstraps", "Number of bootstraps to perform",
                                       false, 0, "integer >= 1");
    TCLAP::SwitchArg testingFlag("s", "testing", "Only compute results for the first 16 individuals");
    TCLAP::ValueArg<std::string> alleleFreqFile("a", "allelefreq", "File containing pre-computed allele frequencies",
		    			       false, "", "path");
    TCLAP::ValueArg<std::string> unrelatedFile("u", "unrelated", "File containing list of unrelated individuals",
		    			       false, "", "path");
    cmd.add(input);
    cmd.add(output);
    cmd.add(numThreads);
    cmd.add(numBootstraps);
    cmd.add(likelihoodFormat);
    cmd.add(type);
    cmd.add(testingFlag);
    cmd.add(alleleFreqFile);
    cmd.add(unrelatedFile);
    try {
        cmd.parse(argc, argv);
        auto numWorkerThreads = numThreads.getValue();
        auto numBootstrapSamples = numBootstraps.getValue();
        std::cerr << "lcMLkin " << versionStr << '\n';
        std::cerr << "==============================" << '\n';
        std::cerr << "using " << numWorkerThreads << " threads\n";



        struct timeval start, end;
        struct timezone tzp;
        gettimeofday(&start, &tzp);

        relatedness r;

        auto inputFileName = input.getValue();
        auto outputFileName = output.getValue();
        r.set_infile(inputFileName.c_str());
        r.set_outfile(outputFileName.c_str());
        r.set_num_workers(numWorkerThreads);


        std::ofstream logStream(outputFileName + ".log");

        std::vector<spdlog::sink_ptr> sinks;
        sinks.push_back(std::make_shared<spdlog::sinks::stderr_sink_mt>());
        sinks.push_back(std::make_shared<spdlog::sinks::ostream_sink_mt>(logStream));
        auto combinedLogger = std::make_shared<spdlog::logger>("jointLog", std::begin(sinks), std::end(sinks));
        spdlog::register_logger(combinedLogger);

	combinedLogger->info("Will generate {} bootstrap samples", numBootstrapSamples);

        if (unrelatedFile.isSet()) {
            std::string& unrelatedFileName = unrelatedFile.getValue();
            r.parse_unrelated_individuals(unrelatedFileName);
        }

	std::unique_ptr<AlleleFrequencyMap> afMap;
	if (alleleFreqFile.isSet()) {
	    afMap.reset(new AlleleFrequencyMap);
	    std::string& alleleFreqFileName = alleleFreqFile.getValue();
	    bool success = afMap->populateFromFile(alleleFreqFileName);
	    if (!success) { 
	      combinedLogger->error("Could not load allele frequency map from [{}]",
			      	    alleleFreqFileName);
	      std::exit(1);
	    }
	}

        auto inferenceType = type.getValue();
        if (inferenceType == "all") {
            r.set_inference_type(InferenceType::ALL_GENOTYPES);
        } else if (inferenceType == "best") {
            r.set_inference_type(InferenceType::BEST_GENOTYPE);
        } else {
            std::cerr << "Error: don't understand inference type "
                << inferenceType << ", must be \"all\" or \"best\"\n";
            std::exit(1);
        }

        auto likelihoodFmt = likelihoodFormat.getValue();
        if (likelihoodFmt == "raw") {
            r.set_likelihood_format(LikelihoodFormat::RAW);
        } else if (likelihoodFmt == "log") {
            r.set_likelihood_format(LikelihoodFormat::LOG);
        } else if (likelihoodFmt == "phred") {
            r.set_likelihood_format(LikelihoodFormat::PHRED);
        } else {
            combinedLogger->error("{} is not a valid likelihood format!", likelihoodFmt);
            std::exit(1);
        }

        std::cout<<"Populating Data"<<std::endl;
        //r.populate_data(); //works correctly
        r.populate_data_new(afMap.get());

        ///std::cout<<"Calculating Allelle Frequencies"<<std::endl;
        //r.calculate_allele_frequencies(); //works correctly

        std::cout<<"Calculating Prestored IBS|IBD"<<std::endl;
        r.calculate_ibs(); //works correctly

        std::cout<<"Starting Pairwise IBD Computations"<<std::endl;
        bool testing = testingFlag.getValue();
        r.calculate_pairwise_ibd(testing, numBootstrapSamples);

        gettimeofday(&end, &tzp);
        print_time_elapsed("", &start, &end);

        combinedLogger->flush();
        logStream.close();

    } catch (TCLAP::ArgException& e) {
        std::cerr << "Exception [" << e.error() << "] when parsing argument "
            << e.argId() << '\n';
        return 1;
    }
    return 0;
}
