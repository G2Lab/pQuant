#include "cxxopts.hpp"
#include "func/Sequence.h"
#include "func/Entropy.h"
#include "openfhe.h"
#include "Task.h"
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
        ("d,data", "dataset", cxxopts::value<std::string>()->default_value("toy"))
        ("k,kmer", "Param kmer", cxxopts::value<int>()->default_value("15"))
        ("v,verbose", "enable verbose", cxxopts::value<bool>()->default_value("false"))
        ("b,bar", "print progress bar", cxxopts::value<bool>()->default_value("false"))
        ("g,gene", "gene path", cxxopts::value<std::string>()->default_value(""))
        ("r,read", "read path", cxxopts::value<std::string>()->default_value(""))
        ("e,thres", "entropy threshold", cxxopts::value<float>()->default_value("0.00001"))
        ("m,kmer_matrix", "load kmer matrix from file path", cxxopts::value<std::string>()->default_value(""))
        ("dng,debug_n_gene", "fix number of genes", cxxopts::value<int>()->default_value("-1"))
        ("gs,gene_start", "starting number index of gene", cxxopts::value<int>()->default_value("-1"))
        ("ge,gene_end", "ending number index of gene", cxxopts::value<int>()->default_value("-1"))
        ("memory,divide_encode_mult", "divide encoding and multiplication steps", cxxopts::value<bool>()->default_value("false"))
        ("h,help", "Print usage")
    ;


    auto result = options.parse(argc, argv);

    if (result.count("help"))
    {
      std::cout << options.help() << std::endl;
      exit(0);
    }

    PQuantParams param(result);

    param.print();

    if (param.target.compare("argument") == 0) {
    } else if (param.target.compare("fasta") == 0) {
        Task::readFastaFiles(param);
    } else if (param.target.compare("table") == 0) {
        Task::testKmerTable(param);
    } else if (param.target.compare("bench") == 0) {
        Task::bfvBenchmark(param);
    } else if (param.target.compare("all") == 0) {
        Task::run_all(param);
    } else if (param.target.compare("json") == 0) {
        Task::testReadJson(param);
    } else {
        std::cout << "Invalid target algorithm" << std::endl;
        exit(0);
    }
    return 0;
}
