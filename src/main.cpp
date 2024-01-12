#include "Functions.h"
#include "PlainFunc.h"
#include "Sequence.h"
#include "cxxopts.hpp"
#include "entropy/Entropy.h"
#include "openfhe.h"
#include "test/Test.h"
#include "PQuantParams.h"
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
        ("h,help", "Print usage")
    ;


    auto result = options.parse(argc, argv);

    if (result.count("help"))
    {
      std::cout << options.help() << std::endl;
      exit(0);
    }

    string filename_read = datasetMap[result["data"].as<std::string>()].first;
    string filename_ref = datasetMap[result["data"].as<std::string>()].second;

    int k = result["kmer"].as<int>();
    string target = result["target"].as<std::string>();
    bool progress_bar = result["bar"].as<bool>();
    bool verbose = result["verbose"].as<bool>();

    PQuantParams param(target, filename_read, filename_ref, k, verbose, progress_bar);

    param.print();


    if (target.compare("fasta") == 0) {
        // filename_read = "../dataset/five_genes/five_gene_reads.fa";
        // string filename_ref = "../dataset/five_genes/five_gene_reference.fa";
        Test::testReadFastaFiles(filename_ref, filename_read);
    } else if (target.compare("bfv") == 0) {
        Test::previousAlgorithm();
    } else if (target.compare("kmer") == 0) {
        Test::kmerTables(param);
    } else if (target.compare("plain") == 0) {
        Test::plainExp();
    } else if (target.compare("bench") == 0) {
        Test::bfvBenchmark(param);
    } else if (target.compare("encread") == 0) {
        Test::encryptRead(param);
    }
    return 0;
}
