#include "cxxopts.hpp"
#include "func/Sequence.h"
#include "func/Entropy.h"
#include "openfhe.h"
#include "Task.h"
#include "MainAlgorithmSet.h"
#include "func/PQuantParams.h"
#include <algorithm>
#include <cstdlib>
#include <ctime>
#include <iostream>

using namespace std;
using namespace lbcrypto;

int main(int argc, char **argv) {
    cxxopts::Options options("pQuant_enc", "BFV implementation of secure pQuant algorithm");

    options.add_options()
        ("t,target", "target algorithm", cxxopts::value<std::string>()->default_value("test"))
        ("d,data", "name of dataset", cxxopts::value<std::string>()->default_value("data"))
        ("k,kmer", "Param kmer", cxxopts::value<int>()->default_value("15"))
        ("v,verbose", "enable verbose", cxxopts::value<bool>()->default_value("false"))
        ("b,bar", "print progress bar", cxxopts::value<bool>()->default_value("false"))
        ("g,gene", "gene path", cxxopts::value<std::string>()->default_value(""))
        ("r,read", "read path", cxxopts::value<std::string>()->default_value(""))
        ("e,thres", "entropy threshold", cxxopts::value<float>()->default_value("0.00001"))
        ("j,json", "use json format for kmer table", cxxopts::value<bool>()->default_value("false"))
        ("kmer_folder", "save/load kmerTable and kmerList from folder path; filename is automatically set by k/data/thres", cxxopts::value<std::string>()->default_value(""))
        ("kmer_table", "save/load kmerTable from file path", cxxopts::value<std::string>()->default_value(""))
        ("kmer_list", "save/load kmerList from file path", cxxopts::value<std::string>()->default_value(""))
        ("bfv_folder", "save/load BFV HE keys from folder path; filename is automatically set (pk, sk, params, etc)", cxxopts::value<std::string>()->default_value(""))
        ("ctxt_read_folder", "path to save encoded and encrypted reads", cxxopts::value<std::string>()->default_value(""))
        ("ctxt_out_folder", "path to save results of step 4", cxxopts::value<std::string>()->default_value(""))
        ("dng,debug_n_gene", "fix number of genes", cxxopts::value<int>()->default_value("-1"))
        ("gs,gene_start", "starting number index of gene", cxxopts::value<int>()->default_value("-1"))
        ("ge,gene_end", "ending number index of gene", cxxopts::value<int>()->default_value("-1"))
        ("ets,operate_then_serialize", "finish encrypt/decrypt and then serialize outputs", cxxopts::value<bool>()->default_value("false"))
        ("h,help", "Print usage")
    ;


    auto result = options.parse(argc, argv);

    if (result.count("help"))
    {
      std::cout << options.help() << std::endl;
      exit(0);
    }

    PQuantParams param(result);

    // test algorithms
    if (param.target.compare("bench") == 0) {
        Task::bfvBenchmark(param);
        return 0;
    }
    
    
    // main algorithms
    param.print();

    if (param.target.compare("STEP1") == 0) {
        MainAlgorithmSet::generateKmerTableFromReference(param);
    } else if (param.target.compare("STEP2") == 0) {
        MainAlgorithmSet::keyGenBFVandSerialize(param);
    } else if (param.target.compare("STEP3") == 0) {
        MainAlgorithmSet::encodeAndEncrypt(param);
    } else if (param.target.compare("STEP4") == 0) {
        MainAlgorithmSet::computeInnerProductBatch(param);
    } else if (param.target.compare("STEP5") == 0) {
        MainAlgorithmSet::decryptAndReturnGeneVector(param);
    } else if (param.target.compare("argument") == 0) {
    } else if (param.target.compare("fasta") == 0) {
        Task::readFastaFiles(param);
    } else if (param.target.compare("table") == 0) {
        Task::testKmerTable(param);
    } else if (param.target.compare("all") == 0) {
        Task::run_all(param);
    } else if (param.target.compare("main") == 0) {
        MainAlgorithmSet::generateKmerTableFromReference(param);
        MainAlgorithmSet::keyGenBFVandSerialize(param);
        MainAlgorithmSet::encodeAndEncrypt(param);
        MainAlgorithmSet::computeInnerProductBatch(param);
        MainAlgorithmSet::decryptAndReturnGeneVector(param);
    } else {
        std::cout << "Invalid target algorithm" << std::endl;
        exit(0);
    }
    return 0;
}
