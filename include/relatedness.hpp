#ifndef RELATEDNESS_H
#define RELATEDNESS_H

#include <limits>
#include <string>
#include <vector>

#include <Eigen/Dense>

#include "utils.hpp"
#include <memory>

#define IBD_COUNT 3
#define GENOTYPE_COUNT 9

// forward declaration
namespace vcflib {
    class VariantCallFile;
}

class relatedness {

private:
    std::unique_ptr<vcflib::VariantCallFile> vcfFile{nullptr};
	std::string infile;

	std::string outfile;

	std::vector<std::string> header;

	std::vector<std::vector<std::string>> snp_data;
    std::vector<std::vector<std::vector<double>>> gtProbs;

	int snp_count;

	bool haveUnrelatedList_{false};

    InferenceType inferenceType;
    LikelihoodFormat likelihoodFormat;

    uint32_t num_worker_threads{1};

	std::vector<std::string> unrelated_individuals;
	std::vector<int> unrelated_individual_index;

	Eigen::VectorXd allele_frequency;

	Eigen::Array<Eigen::MatrixXd, GENOTYPE_COUNT, 1> ibs_all;	//Eigen::Array of GENOTYPE_COUNT Eigen::Matrix

	std::vector<std::pair<int,int>> pairs;

	Eigen::MatrixXd ibs_pairwise;

	Eigen::VectorXd mask_snp;

	Eigen::MatrixXd ibs_best;

	uint32_t GTIndex{std::numeric_limits<uint32_t>::max()};
	uint32_t GQIndex{std::numeric_limits<uint32_t>::max()};
	uint32_t GLIndex{std::numeric_limits<uint32_t>::max()};

public:

	void parse_unrelated_individuals(const std::string& file);

	void populate_data_new();

	void populate_data();

	void calculate_allele_frequencies();

	void calculate_ibs();

	void calculate_pairwise_likelihood(
            std::pair<int,int>,
            Eigen::MatrixXd&,
            Eigen::VectorXd&,
            std::vector<double>&,
            std::vector<double>&
            );

	void calculate_pairwise_ibd(bool testing);

	Eigen::Vector3d optimize_parameters(
            Eigen::Array<Eigen::MatrixXd, GENOTYPE_COUNT, 1>& ibs_all,
            Eigen::VectorXd& mask_snp,
            Eigen::MatrixXd& pibs,
            Eigen::Array<Eigen::MatrixXd, GENOTYPE_COUNT, 1>& workingSpace);

    Eigen::Vector3d optimize_parameters(
            Eigen::MatrixXd& ibs_best,
            Eigen::VectorXd& mask_snp);


	double kin(std::pair<double,double>);

	double gl_kin(std::pair<double,double>);


	Eigen::Vector3d em_optimization(
            Eigen::Vector3d k_values,
            Eigen::Array<Eigen::MatrixXd, GENOTYPE_COUNT, 1>& ibs_all,
            Eigen::VectorXd& mask_snp,
            Eigen::MatrixXd& pibs,
            Eigen::Array<Eigen::MatrixXd, GENOTYPE_COUNT, 1>& workingSpace);



	Eigen::Vector3d em_optimization(
            Eigen::Vector3d k_values,
             Eigen::MatrixXd& ibs_best,
            Eigen::VectorXd& mask_snp);

	void set_infile(const char*);

	void set_outfile(const char*);

    void set_num_workers(uint32_t num_workers) {
        num_worker_threads = num_workers;
    }

    void set_inference_type(InferenceType infType) {
        inferenceType = infType;
    }

    InferenceType getInferenceType() { return inferenceType; }

    void set_likelihood_format(LikelihoodFormat likeFmt) {
        likelihoodFormat = likeFmt;
    }

    LikelihoodFormat getLikelihoodFormat() { return likelihoodFormat; }

	int getSNPCount() { return snp_count; }

    uint32_t getNumWorkers() { return num_worker_threads; }

    std::vector<std::string>& getHeader();
};

#endif
