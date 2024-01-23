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
        ("s,serial", "serialize ciphertexts", cxxopts::value<bool>()->default_value("false"))
        ("g,gene", "gene path", cxxopts::value<std::string>()->default_value(""))
        ("r,read", "read path", cxxopts::value<std::string>()->default_value(""))
        ("o,out", "output path", cxxopts::value<std::string>()->default_value("../crypto/ctxt/"))
        ("h,help", "Print usage")
    ;


    auto result = options.parse(argc, argv);

    if (result.count("help"))
    {
      std::cout << options.help() << std::endl;
      exit(0);
    }

    string filename_read, filename_ref;
    if (result["gene"].as<std::string>().size() > 0 && result["read"].as<std::string>().size() > 0) {
        filename_read = result["read"].as<std::string>();
        filename_ref = result["gene"].as<std::string>();
    } else if (result["data"].as<std::string>().size() > 0) {
        string dataset = result["data"].as<std::string>();
        filename_read = datasetMap[dataset].first;
        filename_ref = datasetMap[dataset].second;
    } else {
        std::cout << "Please specify dataset or gene and read path" << std::endl;
        exit(0);
    }
    
    int k = result["kmer"].as<int>();
    string target = result["target"].as<std::string>();
    bool progress_bar = result["bar"].as<bool>();
    bool verbose = result["verbose"].as<bool>();
    bool serial = result["serial"].as<bool>();
    string output_path = result["out"].as<std::string>();

    PQuantParams param(target, filename_read, filename_ref, output_path, k, verbose, progress_bar, serial);

    param.print();


    if (target.compare("fasta") == 0) {
        Task::readFastaFiles(param);
    } else if (target.compare("kmer") == 0) {
        Task::kmerTables(param);
    } else if (target.compare("table2") == 0) {
        Task::kmerTable2(param);
    } else if (target.compare("bench") == 0) {
        Task::bfvBenchmark(param);
    } else if (target.compare("all") == 0) {
        Task::run_all(param);
    } else if (target.compare("json") == 0) {
        Task::testReadJson(param);
    } else {
        std::cout << "Invalid target algorithm" << std::endl;
        exit(0);
    }
    return 0;
}
