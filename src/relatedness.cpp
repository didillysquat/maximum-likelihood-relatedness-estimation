#include <numeric>
#include <iostream>
#include <fstream>
#include <sstream>

#include <cstdlib>
#include <ctime>
#include <iomanip>

#include <algorithm>
#include <math.h>
#include <limits>

#include "Variant.h"
// #include <omp.h>
#include "tclap/CmdLine.h"
#include "spdlog/spdlog.h"
#include "concurrentqueue.h"
#include "relatedness.hpp"
#include "utils.hpp"

#include "cppitertools/itertools.hpp"

void relatedness::populate_data_new() {
	auto logger = spdlog::get("jointLog");

    if (!vcfFile) {
        vcfFile.reset(new vcflib::VariantCallFile);
    }
	vcfFile->open(infile);
	if (!vcfFile->is_open()) {
	    std::cerr << "Error opening input VCF file :[" << infile << "]!\n";
	    std::exit(1);
	}

    auto& header = vcfFile->sampleNames;

	// Create the unrelated individual list
	if (haveUnrelatedList_) { // If we have a specific list of unrelated individuals, then use it
		logger->info("have unrelated list");
		//std::cerr << "index list = { ";
		for(int i=0; i<header.size();i++){
			if(std::find(unrelated_individuals.begin(), unrelated_individuals.end(), header[i]) != unrelated_individuals.end()) {
				unrelated_individual_index.push_back(i);
				//std::cerr << i << " ";
			}
		}
		//std::cerr << "}\n";
	} else { // Otherwise, assume everyone is unrelated
	    unrelated_individual_index.resize(header.size());
	    std::iota(unrelated_individual_index.begin(),
		      unrelated_individual_index.end(), 0);
	}


    std::vector<std::vector<std::vector<double>>> genotypeLikelihoods;
    std::vector<double> alleleFreqs;

    switch (likelihoodFormat) {
        case LikelihoodFormat::RAW:
            logger->info("Genotype likelihood format is RAW");
            break;
        case LikelihoodFormat::LOG:
            logger->info("Genotype likelihood format is LOG");
            break;
        case LikelihoodFormat::PHRED:
            logger->info("Genotype likelihood format is PHRED");
            break;
    }

    /**
     *  For *some* reason (WHY?) vcflib doesn't return the samples for a variant in the
     *  order in which the appear in the line (left-to-right).  Thus, we need to map from
     *  a sample's name to it's index so as to maintain a globally-consistent order.
     **/
    std::unordered_map<std::string, uint32_t> sampleNameToIndex;
    uint32_t idx{0};
    for (auto n : vcfFile->sampleNames) {
        sampleNameToIndex[n] = idx;
        ++idx;
    }
    // The number of samples we expect for each variant
    size_t numSamples = vcfFile->sampleNames.size();

    auto genotypeFieldName = (likelihoodFormat == LikelihoodFormat::PHRED) ? "PL" : "GL";
    bool excludeBroken{true};
    // Code modified from
    // https://github.com/ekg/vcflib/blob/b1e9b31d2d95f957f86cdfd1e7c9ec25b4950ee8/src/vcfglbound.cpp
    snp_count = 0;
    using std::find;
    using std::reverse;
	vcflib::Variant var(*vcfFile);
    // For each variant
	while (vcfFile->getNextVariant(var)) {
        ++snp_count;

        // If there is no GL field -- note this and continue
        if (find(var.format.begin(), var.format.end(), "GL") == var.format.end()) {
            logger->warn("No GL field in variant: {}", var.sequenceName);
            continue;
        }
        // If we didn't find the "GT" field (what does this mean?)
        if (find(var.format.begin(), var.format.end(), "GT") == var.format.end()) {
            logger->warn("No GT field in variant: {}", var.sequenceName);
            continue;
            //var.format.push_back("GT");
            //reverse(var.format.begin(), var.format.end());
        }

        // Reset for each variant
        bool isbroken = false;
        uint32_t numZeros{0};
        uint32_t totGT{0};

        // New variant entry
        std::vector<std::vector<double>> variantLikelihoods(numSamples, std::vector<double>());

        // For each sample in this variant
        for (auto s = var.samples.begin(); s != var.samples.end(); ++s) {

            auto& sample = s->second;

            // Get the genotype
            auto t = sample.find("GT");
            if (t == sample.end()) {
                logger->error("Couldn't find the required field GT in record {}", s->first);
                std::exit(1);
            }

            auto l = sample.find(genotypeFieldName);

            if (l == sample.end()) {
                logger->error("Couldn't find the required field {} in record {}",
                              genotypeFieldName, s->first);
                std::exit(1);
            }

            // At this point, we have both the GT and GL / PL fields.
            // Check that we have the same number of each

            // The "vector" of genotypes for this sample
            std::vector<std::string>& gtstrs = t->second;
            // The vector of genotype likelihoods for this sample
            std::vector<std::string>& glstrs = l->second;
            // The index of this sample in header order
            auto sampleIndex = sampleNameToIndex[s->first];

            auto& gtstr = gtstrs.front(); // Assume only one GT field per sample per variant.
            std::vector<double> gls;
            bool parseLikelihoods{true};
            // If GT is available, parse it.
            // If one or both of the genotypes are null ".", then
            // fill in the likelihoods as invalid and move on to the
            // next sample for this variant.
            if (gtstr != "./." and gtstr != ".|.") {
                auto gtVec = split(gtstr, "/|");
                // Check that there are only 2 alleles?
                if (gtVec[0] == "." or gtVec[1] == ".") {
                    gls = {-std::numeric_limits<double>::infinity(),
                           -std::numeric_limits<double>::infinity(),
                           -std::numeric_limits<double>::infinity()};
                    parseLikelihoods = false;
                }
                int ga = std::stoi(gtVec[0]);
                int gb = std::stoi(gtVec[1]);
                numZeros += (ga == 0) ? 1 : 0;
                numZeros += (gb == 0) ? 1 : 0;
                totGT += 2;
            } else {
                gls = {-std::numeric_limits<double>::infinity(),
                       -std::numeric_limits<double>::infinity(),
                       -std::numeric_limits<double>::infinity()};
                parseLikelihoods = false;
            }

            if (parseLikelihoods) {
                for (auto gl = glstrs.begin(); gl != glstrs.end(); ++gl) {
                    const auto& glstr = *gl;
                    double val = -std::numeric_limits<double>::infinity();
                    if (glstr != ".") {
                        const char* glcstr = glstr.c_str();
                        char* end;
                        val = std::strtod(glcstr, &end);
                    }
                    gls.push_back(val);
                }

                // Do checking of likelihoods and any conversion here
                // Here, we convert logspace -> rawspace --- the other way around
                // might be better
                switch (likelihoodFormat) {
                    case LikelihoodFormat::LOG :
                        for (auto& g : gls) {
                            if (std::isfinite(g)) {
                                g = std::pow(10.0, g);
                            }
                        }
                        break;
                    case LikelihoodFormat::PHRED :
                        for (auto& g : gls) {
                            if (std::isfinite(g)) {
                                g = std::pow(10.0, -g / 10.0);
                            }
                        }
                        break;
                    case LikelihoodFormat::RAW :
                        // Convert invalid values to our sentinel
                        // (-inf)
                        for (auto& g : gls) {
                            if (std::isfinite(g)) {
                                // Raw values should not be negative
                                // (check for > 1)?
                                if (g < 0.0) {
                                    g = -std::numeric_limits<double>::infinity();
                                }
                            }
                        }
                        break;
                }
            }
            variantLikelihoods[sampleIndex] = gls;
        }
        genotypeLikelihoods.push_back(std::move(variantLikelihoods));

        if (totGT == 0) {
            logger->error("No samples for variant {}!", var.sequenceName);
            std::exit(1);
        }

        // Push back the dominant allele frequency
        double af = static_cast<float>(numZeros) / totGT;
        alleleFreqs.push_back(af);

        if (excludeBroken && isbroken) {
            logger->warn("excluding VCF record @ {} : {} due to GLs > 0", var.sequenceName, var.position);
        }
    }
    logger->info("Parsed variants at {} total sites", snp_count);
    logger->info("allele freqs = {}", alleleFreqs.size());
    logger->info("gt likelihoods = {}", genotypeLikelihoods.size());

    std::swap(gtProbs, genotypeLikelihoods);
    // Probably a better way to do this, but deal with that later
    size_t zeroCount{0};
    size_t oneCount{0};
    allele_frequency.resize(alleleFreqs.size());
    for (size_t i = 0; i < alleleFreqs.size(); ++i) {
        allele_frequency(i) = alleleFreqs[i];
        if (allele_frequency(i) == 0.0) { zeroCount++; }
        if (allele_frequency(i) == 1.0) { oneCount++; }
    }
    std::cerr << "number of SNPs masked (freq 0) = " << zeroCount << "\n";
    std::cerr << "number of SNPs masked (freq 1) = " << oneCount << "\n";
}


void relatedness::populate_data(){

	//Read the VCF file line by line
	std::ifstream vcfFile (infile);
  	std::vector<std::string> data;
  	std::string line;
	std::string commentPrefix = "##";
  	while (std::getline(vcfFile, line)){
  		if(!line.empty()) {
		    // Ignore comment lines (i.e. lines starting with "##")
		    if (!std::equal(commentPrefix.begin(), commentPrefix.end(), line.begin())) {
 	 		data.push_back(line);
		    }
		}
	}
  	vcfFile.close();

	if (data[0].front() != '#') {
	    std::cerr << "I could not find a header line in VCF file. Please make sure it is properly formatted.\n";
	    std::exit(1);
	}

  	std::vector<std::string> head = split(data[0],'\t');
  	header = {head.begin()+9, head.end()};

	if (haveUnrelatedList_) { // If we have a specific list of unrelated individuals, then use it
		std::cerr << "have unrelated list\n";
		std::cerr << "index list = { ";
		for(int i=0; i<header.size();i++){
			if(std::find(unrelated_individuals.begin(), unrelated_individuals.end(), header[i]) != unrelated_individuals.end()) {
				unrelated_individual_index.push_back(i);
				std::cerr << i << " ";
			}
		}
		std::cerr << "}\n";
	} else { // Otherwise, assume everyone is unrelated
	    unrelated_individual_index.resize(header.size());
	    std::iota(unrelated_individual_index.begin(),
		      unrelated_individual_index.end(), 0);
	}


	std::string sampleLine;
	for(auto iter = data.begin() + 1; iter != data.end(); iter++){
		if (!(iter->find('#') == 0)) { // If this isn't a comment line
			// Copy over a sample line to analyze for GQ, GL and GT fields
			if (sampleLine.empty()) {
				sampleLine = *iter;
			}
			std::vector<std::string> elements = split(*iter,'\t');
			std::vector<std::string> elementsOut;
			std::copy(elements.begin()+9, elements.end(), std::back_inserter(elementsOut));
			if(elementsOut.size()>1){
				snp_data.push_back(elementsOut);
			}
		}
	}

	// Determine which indices in INFO are GT, GL and GQ ---
	// This code assumes it is the same for *all* records!
	size_t infoIndex = 8;
	auto info = split(sampleLine, '\t')[infoIndex];
	auto infoFields = split(info, ':');
	size_t fieldNum{0};
	for (auto& f : infoFields) {
	    if (f == "GT") {
	        GTIndex = fieldNum;
	    } else if ( f == "GQ") {
		GQIndex = fieldNum;
	    } else if ( f == "GL") {
		GLIndex = fieldNum;
	    }
	    ++fieldNum;
	}
	std::cerr << "GTIndex = " << GTIndex << ", GQIndex = " << GQIndex << ", GLIndex = " << GLIndex << "\n";

	snp_count = snp_data.size();
	std::cerr << "Parsed " << snp_count << " SNPs\n";
}

void relatedness::calculate_allele_frequencies(){

	allele_frequency = Eigen::VectorXd::Zero(snp_count);

	for(int i=0; i<snp_count; i++){

		std::vector<int> counts;

		for(int id: unrelated_individual_index){
			std::string geno = split(snp_data[i][id],':')[GTIndex];
			if(geno != "./."){
				counts.push_back(std::stoi(std::string(1,geno.front())));
				counts.push_back(std::stoi(std::string(1,geno.back())));
			}
		}
		int zeroes=0;
		for(auto iter=counts.begin(); iter!=counts.end(); iter++){
			if(*iter==0) zeroes++;
		}
		allele_frequency(i) = zeroes/(double)counts.size();
	}

}

void relatedness::calculate_ibs(){

	for(int i=0; i<GENOTYPE_COUNT; i++){
		ibs_all(i) = Eigen::MatrixXd::Zero(snp_count,IBD_COUNT);
	}

	//Populate matrices with actual probablities
	for(int i=0; i<snp_count; i++){

		double p = allele_frequency(i);
		double q = 1.0-p;

		ibs_all(0)(i,0) = pow(p,2)*pow(q,2); 	ibs_all(0)(i,1) = 0;			 		ibs_all(0)(i,2) = 0; //PP QQ
    	ibs_all(1)(i,0) = pow(q,2)*pow(p,2);	ibs_all(1)(i,1) = 0; 					ibs_all(1)(i,2) = 0; //QQ PP

		ibs_all(2)(i,0) = pow(p,2)*(2*p*q);		ibs_all(2)(i,1) = pow(p,2)*q; 			ibs_all(2)(i,2) = 0; //PP PQ
		ibs_all(3)(i,0) = (2*p*q)*pow(p,2);		ibs_all(3)(i,1) = (p*q)*p; 				ibs_all(3)(i,2) = 0; //PQ PP
		ibs_all(4)(i,0) = (2*p*q)*pow(q,2);		ibs_all(4)(i,1) = (p*q)*q; 				ibs_all(4)(i,2) = 0; //PQ QQ
		ibs_all(5)(i,0) = pow(q,2)*(2*p*q);		ibs_all(5)(i,1) = pow(q,2)*p;			ibs_all(5)(i,2) = 0; //QQ PQ

		ibs_all(6)(i,0) = pow(p,4); 			ibs_all(6)(i,1) = pow(p,3);				ibs_all(6)(i,2) = pow(p,2); //PP PP
		ibs_all(7)(i,0) = (2*p*q)*(2*p*q); 		ibs_all(7)(i,1) = p*q*(p+q);			ibs_all(7)(i,2) = 2*(p*q);  //PQ PQ
		ibs_all(8)(i,0) = pow(q,4);		 		ibs_all(8)(i,1) = pow(q,3);				ibs_all(8)(i,2) = pow(q,2); //QQ QQ

	}

}

inline bool isInvalidLikelihood(double l) {
    return std::isinf(l);
}

void relatedness::calculate_pairwise_likelihood(
		std::pair<int,int> pair,
		Eigen::MatrixXd& ibs_pairwise,
		Eigen::VectorXd& mask_snp,
		std::vector<double>& likelihood_1,
		std::vector<double>& likelihood_2){


    // NEW
    auto id1 = pair.first;
    auto id2 = pair.second;
    // For every variant (SNP)
	for(auto j : iter::range(snp_count)) {

        double p=allele_frequency(j);
		// double q=1.0-p;

		//Mask SNPs where allele frequency is fixed
		if (p == 1.0 or p == 0.0){
			mask_snp(j) = 1;
		}

        // Likelihoods for the first and second individuals
        auto& l1 = gtProbs[j][id1];
        auto& l2 = gtProbs[j][id2];
        for (size_t i = 0; i < l1.size(); ++i) {
            if (isInvalidLikelihood(l1[i]) or
                isInvalidLikelihood(l2[i])) {
                mask_snp(j) = 1;
            }
        }
		//Calculate the probability of observing all possible two genotype combinations by multiplying their likelihoods
		ibs_pairwise(j,0)=(l1[0] * l2[2]);
		ibs_pairwise(j,1)=(l1[2] * l2[0]);

		ibs_pairwise(j,2)=(l1[0] * l2[1]);
		ibs_pairwise(j,3)=(l1[1] * l2[0]);
		ibs_pairwise(j,4)=(l1[1] * l2[2]);
		ibs_pairwise(j,5)=(l1[2] * l2[1]);

		ibs_pairwise(j,6)=(l1[0] * l2[0]);
		ibs_pairwise(j,7)=(l1[1] * l2[1]);
		ibs_pairwise(j,8)=(l1[2] * l2[2]);
	}

    // OLD

    /*
	//Parallallize
	//Iterate through all SNPs
	for(int j=0; j<snp_count; j++){

		double p=allele_frequency(j);
		// double q=1.0-p;

		//Mask SNPs where allele frequency is fixed
		if (p == 1.0 or p == 0.0){
			mask_snp(j)=1;
		}

		//Pulls out info field for first individual from VCF, pulls out the three precomputed genotype likelihoods
		std::vector<StringItPair> l1 = split_it(split_it(snp_data[j][pair.first],':')[GLIndex],',');
		std::vector<StringItPair> l2 = split_it(split_it(snp_data[j][pair.second],':')[GLIndex],',');

		//Assert: l1==l2?
		//double* likelihood_1 = new double[l1.size()];
		//double* likelihood_2 = new double[l2.size()];
		if (likelihood_1.size() < l1.size()) {
			likelihood_1.resize(l1.size());
		}
		if (likelihood_2.size() < l2.size()) {
			likelihood_2.resize(l2.size());
		}

		//Convert likelihoods from strings to floating point numbers
		for(int k=0; k<l1.size(); k++){
			char* e1 = &*l1[k].end;
			char* e2 = &*l2[k].end;
			likelihood_1[k]=std::strtod(&(*l1[k].begin), &e1);
			likelihood_2[k]=std::strtod(&(*l2[k].begin), &e2);
			//If one of those likelihoods comes out as negative (anomaly), we mask those SNPs
			if(likelihood_1[k] == -9.0 or likelihood_2[k] == -9.0){
				mask_snp[j]=1;
			}
		}

		//Calculate the probability of observing all possible two genotype combinations by multiplying their likelihoods
		ibs_pairwise(j,0)=(likelihood_1[0]*likelihood_2[2]);
		ibs_pairwise(j,1)=(likelihood_1[2]*likelihood_2[0]);

		ibs_pairwise(j,2)=(likelihood_1[0]*likelihood_2[1]);
		ibs_pairwise(j,3)=(likelihood_1[1]*likelihood_2[0]);
		ibs_pairwise(j,4)=(likelihood_1[1]*likelihood_2[2]);
		ibs_pairwise(j,5)=(likelihood_1[2]*likelihood_2[1]);

		ibs_pairwise(j,6)=(likelihood_1[0]*likelihood_2[0]);
		ibs_pairwise(j,7)=(likelihood_1[1]*likelihood_2[1]);
		ibs_pairwise(j,8)=(likelihood_1[2]*likelihood_2[2]);

		//Clean up
		//delete[] likelihood_1;
		//delete[] likelihood_2;
	}
    */
}

std::vector<std::string>& relatedness::getHeader() { return vcfFile->sampleNames; }

void pairwiseLikelihoodWorker(
        Eigen::MatrixXd& ibs_pairwise,
        Eigen::MatrixXd& ibs_best,
        Eigen::Array<Eigen::MatrixXd, GENOTYPE_COUNT, 1>& ibs_all,
        Eigen::VectorXd& mask_snp,
        relatedness& relateObj,
        std::vector<double>& l1,
        std::vector<double>& l2,
        InferenceType inferenceType,
        std::atomic<bool>& queueFilled,
        std::atomic<uint32_t>& numWorking,
        moodycamel::ConcurrentQueue<std::pair<int, int>>& workQueue,
        moodycamel::ConcurrentQueue<std::string>& resultQueue) {


    auto snp_count = relateObj.getSNPCount();
    Eigen::Array<Eigen::MatrixXd, GENOTYPE_COUNT, 1> workingSpace;
	for(int i = 0; i < GENOTYPE_COUNT; i++){
		workingSpace(i) = Eigen::MatrixXd::Zero(snp_count,IBD_COUNT);
	}

    fmt::MemoryWriter resWriter;
    std::pair<int, int> indPair;
    bool gotNewPair{false};
    while (!queueFilled or (gotNewPair = workQueue.try_dequeue(indPair))) {
        if (!gotNewPair) { continue; }
        gotNewPair = false;
        mask_snp.setZero();
        ibs_pairwise.setZero();
        ibs_best.setZero();
        relateObj.calculate_pairwise_likelihood(indPair, ibs_pairwise,
                                                mask_snp, l1, l2);

        //Identify the most likely genotype combination: Index of genotype for that SNP
        //For each SNP, given the best genotype combination, pull out the appropriate P(IBS|IBD) for all three IBS possibilities
        //Eigen::MatrixXd ibs_best = Eigen::MatrixXd::Zero(snp_count,IBD_COUNT);
        //auto& ibs_best = ibs_best_mats[id];

        for(size_t j = 0; j < relateObj.getSNPCount(); j++){
            Eigen::MatrixXf::Index bestIndex;
            ibs_pairwise.row(j).maxCoeff(&bestIndex);
            for(int k=0; k<IBD_COUNT; k++){
                ibs_best(j,k) = ibs_all(bestIndex)(j,k);
            }
        }

        Eigen::Vector3d k_est;
        switch (inferenceType) {
            case InferenceType::BEST_GENOTYPE:
               k_est = relateObj.optimize_parameters(ibs_best, mask_snp);
               break;
            case InferenceType::ALL_GENOTYPES:
               k_est = relateObj.optimize_parameters(ibs_all,
                                                     mask_snp,
                                                     ibs_pairwise,
                                                     workingSpace);
               break;
        }

        auto& ind1 = relateObj.getHeader()[indPair.first];
        auto& ind2 = relateObj.getHeader()[indPair.second];
        if (indPair.first == indPair.second) {
            std::cerr << "i1 = "<< indPair.first << ", i2 = " << indPair.second << '\n';
            std::exit(1);
        }
        resWriter.write("{}\t{}\t{:.3f}\t{:.3f}\t{:.3f}\t{:.3f}\t{}\n",
                        ind1, ind2,
                        k_est(0), k_est(1), k_est(2),
                        0.5 * k_est(1) + k_est(2),
                        mask_snp.size() - mask_snp.sum());
        resultQueue.enqueue(resWriter.str());
        resWriter.clear();
        /*
		outfile << pairs[i].first+1 << "\t" << pairs[i].second+1 << "\t"
				<< std::setprecision(3)
				<< k_est(0) << "\t" << k_est(1) << "\t" << k_est(2) << "\t"
				<< 0.5*k_est(1)+k_est(2) << "\t"
				<< mask_snp.size()-mask_snp.sum() << "\n";
                */
    }
    --numWorking;

}

void relatedness::calculate_pairwise_ibd(bool testing){

	#ifdef DEBUG
	std::cout << "Ind1\tInd2\tk0_hat\tk1_hat\tk2_hat\tpi_HAT\tnbSNP\n";
	#endif

	std::string output_file (outfile);
	std::ofstream outfile (output_file);
	outfile << "Ind1\tInd2\tk0_hat\tk1_hat\tk2_hat\tpi_HAT\tnbSNP\n";

	std::srand(std::time(0));

    uint32_t numThreads = getNumWorkers();

    std::vector<std::vector<double>> l1s(numThreads);
    std::vector<std::vector<double>> l2s(numThreads);

    std::vector<Eigen::MatrixXd> ibs_pairwise_mats(numThreads,
                                  Eigen::MatrixXd::Zero(snp_count, GENOTYPE_COUNT));

    std::vector<Eigen::VectorXd> mask_snp_vecs(numThreads,
                                  Eigen::VectorXd::Zero(snp_count));

    std::vector<Eigen::MatrixXd> ibs_best_mats(numThreads,
                                  Eigen::MatrixXd::Zero(snp_count, IBD_COUNT));

    // Spawn off numThreads separate threads to do the work
    std::vector<std::thread> workers;
    // true when work queue is done being filled
    std::atomic<bool> queueFilled{false};
    // counts down and reaches 0 when all threads are done working
    std::atomic<uint32_t> numWorking{numThreads};
    moodycamel::ConcurrentQueue<std::pair<int, int>> workQueue;
    moodycamel::ConcurrentQueue<std::string> resultQueue;
    auto inferenceType = getInferenceType();
    for (size_t i = 0; i < numThreads; ++i) {
        workers.emplace_back(pairwiseLikelihoodWorker,
                              std::ref(ibs_pairwise_mats[i]),
                              std::ref(ibs_best_mats[i]),
                              std::ref(ibs_all),
                              std::ref(mask_snp_vecs[i]),
                              std::ref(*this),
                              std::ref(l1s[i]),
                              std::ref(l2s[i]),
                              inferenceType,
                              std::ref(queueFilled),
                              std::ref(numWorking),
                              std::ref(workQueue),
                              std::ref(resultQueue));
    }

    //Generate Pairs
    size_t maxI = testing ? 15 : vcfFile->sampleNames.size();//header.size();
    size_t maxJ = testing ? 16 : vcfFile->sampleNames.size();//header.size();
	for(int i = 0; i < maxI; i++){
		for(int j = i+1; j < maxJ; j++){
			//if(split(header[i],'_')[0]==split(header[j],'_')[0])
                	workQueue.enqueue(std::make_pair(i,j));
		}
	}
    queueFilled = true;

    std::string resLine;
    bool obtainedRes{false};
    size_t numRes{0};
    while ((obtainedRes = resultQueue.try_dequeue(resLine)) or numWorking > 0) {
        if (obtainedRes) {
            ++numRes;
            if (numRes % 100 == 0) {
                std::cerr << "\r\rComputed results for " << numRes << " pairs";
            }
            outfile << resLine;
        }
    }
    std::cerr << "\n";

    // join all the worker threads
    for (auto& t : workers) { t.join(); }
    /*
	//Iterate through all pairwise computations
	#pragma omp parallel for
	for(int i=0; i<pairs.size(); i++) {
        // The index of this thread
        size_t id = omp_get_thread_num();

		//Matrices for all possible pairs of genotypes for every SNP.
		//This will eventually store genotype likelihoods based on the calling likelihood (for example based on read depth)
		//Store these pairwise likelihoods in an array AATT,TTAA,AAAT,ATAA,ATTT,TTAT,AAAA,ATAT,TTTT
		//Eigen::MatrixXd ibs_pairwise = Eigen::MatrixXd::Zero(snp_count,GENOTYPE_COUNT);
        auto& ibs_pairwise = ibs_pairwise_mats[id];

		//A matrix to denote SNPs we may want to mask for two reasons(see below)
		//Eigen::VectorXd mask_snp = Eigen::VectorXd::Zero(snp_count);
        auto& mask_snp = mask_snp_vecs[id];

        auto& l1 = l1s[id];
        auto& l2 = l2s[id];
		calculate_pairwise_likelihood(pairs[i], ibs_pairwise, mask_snp, l1, l2);

		//Identify the most likely genotype combination: Index of genotype for that SNP
		//For each SNP, given the best genotype combination, pull out the appropriate P(IBS|IBD) for all three IBS possibilities
		//Eigen::MatrixXd ibs_best = Eigen::MatrixXd::Zero(snp_count,IBD_COUNT);
        auto& ibs_best = ibs_best_mats[id];

		for(int j=0; j<snp_count; j++){
			Eigen::MatrixXf::Index bestIndex;
			double max = ibs_pairwise.row(j).maxCoeff(&bestIndex);
			for(int k=0; k<IBD_COUNT; k++){
				ibs_best(j,k) = ibs_all(bestIndex)(j,k);
			}
		}

		Eigen::Vector3d k_est = optimize_parameters(ibs_best, mask_snp);
		//floor(x*10.0)/10.0
		#ifdef DEBUG
		std::cout << pairs[i].first+1 << "\t" << pairs[i].second+1 << "\t"
				<< std::setprecision(3)
				<< k_est(0) << "\t" << k_est(1) << "\t" << k_est(2) << "\t"
				<< 0.5*k_est(1)+k_est(2) << "\t"
				<< mask_snp.size()-mask_snp.sum() << "\n";
		#endif

		#pragma omp critical
		outfile << pairs[i].first+1 << "\t" << pairs[i].second+1 << "\t"
				<< std::setprecision(3)
				<< k_est(0) << "\t" << k_est(1) << "\t" << k_est(2) << "\t"
				<< 0.5*k_est(1)+k_est(2) << "\t"
				<< mask_snp.size()-mask_snp.sum() << "\n";

	}
    */

}


Eigen::Vector3d relatedness::optimize_parameters(
        Eigen::Array<Eigen::MatrixXd, GENOTYPE_COUNT, 1>& ibs_all,
        Eigen::VectorXd& mask_snp,
        Eigen::MatrixXd& pibs,
        Eigen::Array<Eigen::MatrixXd, GENOTYPE_COUNT, 1>& Xall){

	Eigen::Vector3d k_values = Eigen::Vector3d::Random();
	k_values = k_values.cwiseAbs();
	k_values /= k_values.sum();

	return em_optimization(k_values, ibs_all, mask_snp, pibs, Xall);
}



Eigen::Vector3d relatedness::optimize_parameters(
        Eigen::MatrixXd& ibs_best,
        Eigen::VectorXd& mask_snp){

	Eigen::Vector3d k_values = Eigen::Vector3d::Random();
	k_values = k_values.cwiseAbs();
	k_values /= k_values.sum();

	return em_optimization(k_values, ibs_best, mask_snp);


/*
	for(int j=0; j<3; j++){
		bool flag = false;
		while(flag==false){
			double k1 = rand()/(double)(RAND_MAX);
			double k2 = rand()/(double)(RAND_MAX);
			if(k1+k2<=1.0){
				double k0=1.0-(k1+k2);
				if(4*k2*k0<pow(k1,2)){
					std::pair<double,double> k12 = std::make_pair(k1,k2);
					if(kin(k12) != 100000){
						flag = true;
					}
				}
			}
		}
		//to-do: minimize function
		//to-do: opt_parameter <- parameter that minimizes the function
	}

	for(int j=0; j<3; j++){
		bool flag = false;
		while(flag==false){
			double k1 = rand()/(double)(RAND_MAX);
			double k2 = rand()/(double)(RAND_MAX);
			if(k1+k2<=1.0){
				double k0=1.0-(k1+k2);
				if(4*k2*k0<pow(k1,2)){
					std::pair<double,double> k12 = std::make_pair(k1,k2);
					if(gl_kin(k12) != 100000){
						flag = true;
					}
				}
			}
		}
		//to-do: minimize function
		//to-do: opt_parameter <- parameter that minimizes the function
	}
*/
}

double relatedness::kin(std::pair<double,double> k12){

	double k0 = 1.0-(k12.first+k12.second);
	double k1 = k12.first;
	double k2 = k12.second;

	Eigen::Vector3d k;
	k << k0,k1,k2;

	Eigen::VectorXd ibs_k_sum = Eigen::VectorXd::Zero(snp_count);
	double ibs_sum{0.0};
	for(int i=0; i<snp_count; i++){
		for(int j=0; j<IBD_COUNT;j++){
			ibs_k_sum(i)+=ibs_best(i,j)*k(j);
		}
		if(mask_snp(i)==1){
			ibs_k_sum(i)=0;
		}else{
			ibs_k_sum(i) = std::log(ibs_k_sum(i));
		}
		ibs_sum+=ibs_k_sum(i);
	}
	ibs_sum*=-1;

	bool flag=false;
	if(k0<0 || k0>1) {flag=true;}
	if(k1<0 || k1>1) {flag=true;}
	if(k2<0 || k2>1) {flag=true;}
    if(4*k2*k0 > pow(k1,2)) {flag=true;}
    if (ibs_sum==std::numeric_limits<double>::infinity()) {flag=true;}

    if(flag==true){
        ibs_sum=100000;
    }

    return ibs_sum;
}

double relatedness::gl_kin(std::pair<double,double> k12){

	double k0 = 1.0-(k12.first+k12.second);
	double k1 = k12.first;
	double k2 = k12.second;

	Eigen::Vector3d k_values;
	k_values << k0,k1,k2;

	Eigen::MatrixXd ibs_pairwise_t = Eigen::MatrixXd::Zero(GENOTYPE_COUNT,snp_count);
	ibs_pairwise_t = ibs_pairwise.transpose();

	Eigen::MatrixXd ibs_k = Eigen::MatrixXd::Zero(GENOTYPE_COUNT,snp_count);
	Eigen::VectorXd ibs_k_sum = Eigen::VectorXd::Zero(snp_count);
	double ibs_sum{0.0};

	for(int i=0; i<GENOTYPE_COUNT; i++){
		for(int j=0; j<snp_count; j++){
			for(int k=0; k<IBD_COUNT;k++){
				ibs_k(i,j)+=ibs_all(i)(j,k)*k_values(k);
			}
		}
		for(int j=0; j<snp_count; j++){
			ibs_k(i,j)*=ibs_pairwise_t(i,j);
		}
	}
	for(int i=0; i<snp_count; i++){
		for(int j=0; j<GENOTYPE_COUNT; j++){
			ibs_k_sum(i)+=ibs_k(j,i);
		}
		if(mask_snp(i)==1){
			ibs_k_sum(i)=0;
		}else{
			ibs_k_sum(i) = std::log(ibs_k_sum(i));
		}
		ibs_sum+=ibs_k_sum(i);
	}
	ibs_sum*=-1;

	bool flag=false;
	if(k0<0 || k0>1) {flag=true;}
	if(k1<0 || k1>1) {flag=true;}
	if(k2<0 || k2>1) {flag=true;}
    if(4*k2*k0 > pow(k1,2)) {flag=true;}
    if (ibs_sum==std::numeric_limits<double>::infinity()) {flag=true;}

    if(flag==true){
        ibs_sum=100000;
    }

    return ibs_sum;
}



// All-genotype variant
Eigen::Vector3d relatedness::em_optimization(
        Eigen::Vector3d k_values,
    	Eigen::Array<Eigen::MatrixXd, GENOTYPE_COUNT, 1>& ibs_all,
        Eigen::VectorXd& mask_snp,
        Eigen::MatrixXd& pibs,
        Eigen::Array<Eigen::MatrixXd, GENOTYPE_COUNT, 1>& Xall){

    //Probabilities of IBD for each SNP
	Eigen::MatrixXd ibd_probability = Eigen::MatrixXd::Zero(snp_count,IBD_COUNT);


    Eigen::Array<Eigen::MatrixXd, GENOTYPE_COUNT, 1> XPre;
    XPre = ibs_all;
    for(size_t gtIdx = 0; gtIdx < GENOTYPE_COUNT; ++gtIdx) {
        XPre(gtIdx).col(0).array() *= pibs.col(gtIdx).array();
        XPre(gtIdx).col(1).array() *= pibs.col(gtIdx).array();
        XPre(gtIdx).col(2).array() *= pibs.col(gtIdx).array();
    }



	//Difference bwtween subsequent parameter estimates
	double thresh = 100;
	//Iteration number
	int iter = 0;

	while(thresh > 1e-4){
		Eigen::MatrixXd X = Eigen::MatrixXd::Zero(snp_count,IBD_COUNT);
		Eigen::VectorXd XS = Eigen::VectorXd::Zero(snp_count); //Used to normalize X

        for(int j=0; j< IBD_COUNT; j++){
            for(size_t gtIdx = 0; gtIdx < GENOTYPE_COUNT; ++gtIdx) {
                Xall(gtIdx).col(j).noalias() = XPre(gtIdx).col(j)*k_values(j);
            }
        }
        /*
        for(int j=0; j< IBD_COUNT; j++){
            for(size_t gtIdx = 0; gtIdx < GENOTYPE_COUNT; ++gtIdx) {
                Xall(gtIdx).col(j).noalias() = ibs_all(gtIdx).col(j) * k_values(j);
            }
        }

        for(size_t gtIdx = 0; gtIdx < GENOTYPE_COUNT; ++gtIdx) {
            Xall(gtIdx).col(0).array() *= pibs.col(gtIdx).array();
            Xall(gtIdx).col(1).array() *= pibs.col(gtIdx).array();
            Xall(gtIdx).col(2).array() *= pibs.col(gtIdx).array();
        }
        */

        // Sum up all genotypes to get X
        for(size_t gtIdx = 0; gtIdx < GENOTYPE_COUNT; ++gtIdx) {
            X += Xall(gtIdx);
        }

		XS = X.rowwise().sum();
		for(int i=0;i <snp_count; i++){ //Mask
			if(mask_snp(i)==1){
				XS(i)=0;
			}
		}

		//Normalize X
		for(int i=0; i<snp_count; i++){
			if(XS(i)==0){
				for(int j=0; j<IBD_COUNT; j++){
					if(X(i,j)!=0){
						X(i,j)=1;
					}
				}
			}
			else{
				for(int j=0; j<IBD_COUNT; j++){
					X(i,j)/=XS(i);
				}
			}
		}

		//Copy X to PIBD
		ibd_probability=X;

		//New parameter estimates
		Eigen::Vector3d k_est = Eigen::Vector3d::Zero();

		//Sum of probabilities at each site
		for(int i=0; i<IBD_COUNT; i++){
			for(int j=0; j<snp_count; j++){
				if(mask_snp(j)!=1){ //Mask
					k_est(i) += X(j,i);
				}
			}
		}

		//Normalized estimates
		k_est /= k_est.sum();

		// Compute the difference between successive estimates to assess convergence
		thresh = (k_est-k_values).cwiseAbs().norm();
        k_values = k_est;
        iter++;
	}

	return k_values;
}


//Currently implements inference using the best genotypes (ibs_best)
Eigen::Vector3d relatedness::em_optimization(
        Eigen::Vector3d k_values,
        Eigen::MatrixXd& ibs_best,
        Eigen::VectorXd& mask_snp){

	//Probabilities of IBD for each SNP
	Eigen::MatrixXd ibd_probability = Eigen::MatrixXd::Zero(snp_count,IBD_COUNT);

	//Difference bwtween subsequent parameter estimates
	double thresh = 100;
	//Iteration number
	int iter = 0;

	while(thresh > 1e-4){

		Eigen::MatrixXd X = Eigen::MatrixXd::Zero(snp_count,IBD_COUNT);
		Eigen::VectorXd XS = Eigen::VectorXd::Zero(snp_count); //Used to normalize X

		for(int i=0; i<snp_count; i++){
			for(int j=0; j<IBD_COUNT; j++){
				X(i,j)=ibs_best(i,j)*k_values(j);
			}
		}

		XS = X.rowwise().sum();
		for(int i=0;i <snp_count; i++){ //Mask
			if(mask_snp(i)==1){
				XS(i)=0;
			}
		}

		//Normalize X
		for(int i=0; i<snp_count; i++){
			if(XS(i)==0){
				for(int j=0; j<IBD_COUNT; j++){
					if(X(i,j)!=0){
						X(i,j)=1;
					}
				}
			}
			else{
				for(int j=0; j<IBD_COUNT; j++){
					X(i,j)/=XS(i);
				}
			}
		}

		//Copy X to PIBD
		ibd_probability=X;

		//New parameter estimates
		Eigen::Vector3d k_est = Eigen::Vector3d::Zero();

		//Sum of probabilities at each site
		for(int i=0; i<IBD_COUNT; i++){
			for(int j=0; j<snp_count; j++){
				if(mask_snp(j)!=1){ //Mask
					k_est(i)+=ibd_probability(j,i);
				}
			}
		}

		//Normalized estimates
		k_est /= k_est.sum();

		// Compute the difference between successive estimates to assess convergence
		thresh = (k_est-k_values).cwiseAbs().norm();
        k_values = k_est;
        iter++;
	}

	return k_values;
}

void relatedness::set_infile(const char* filename){

	infile = std::string(filename);
}

void relatedness::set_outfile(const char* filename){

	outfile = std::string(filename);
}


void relatedness::parse_unrelated_individuals(const std::string& file) {
    std::ifstream inFile(file);
    if (!inFile.is_open()) {
        std::cerr << "Couldn't open unrelated individual file [" << file << "]!\n";
	std::exit(1);
    }
    for (std::string line; std::getline(inFile, line);) {
    	unrelated_individuals.emplace_back(line);
    }
    haveUnrelatedList_ = true;
}

/*
void usage(char *program) {
    std::cerr << "Usage: " << program << " <input_file_path> <output_file_path> <genotype>" << std::endl;
    exit(1);
}
*/

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
    TCLAP::SwitchArg testingFlag("s", "testing", "Only compute results for the first 16 individuals");
    TCLAP::ValueArg<std::string> unrelatedFile("u", "unrelated", "File containing list of unrelated individuals",
		    			       false, "", "path");
    cmd.add(input);
    cmd.add(output);
    cmd.add(numThreads);
    cmd.add(likelihoodFormat);
    cmd.add(type);
    cmd.add(testingFlag);
    cmd.add(unrelatedFile);
    try {
        cmd.parse(argc, argv);
        auto numWorkerThreads = numThreads.getValue();
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

	if (unrelatedFile.isSet()) {
	  std::string& unrelatedFileName = unrelatedFile.getValue();
	  r.parse_unrelated_individuals(unrelatedFileName);
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
        r.populate_data_new();

        ///std::cout<<"Calculating Allelle Frequencies"<<std::endl;
        //r.calculate_allele_frequencies(); //works correctly

        std::cout<<"Calculating Prestored IBS|IBD"<<std::endl;
        r.calculate_ibs(); //works correctly

        std::cout<<"Starting Pairwise IBD Computations"<<std::endl;
        bool testing = testingFlag.getValue();
        r.calculate_pairwise_ibd(testing);

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
