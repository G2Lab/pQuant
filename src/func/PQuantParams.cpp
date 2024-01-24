#include "func/PQuantParams.h"

PQuantParams::PQuantParams(cxxopts::ParseResult &result) {
    if (result["gene"].as<std::string>().size() > 0 && result["read"].as<std::string>().size() > 0) {
        this->path_filename_read = result["read"].as<std::string>();
        this->path_filename_ref = result["gene"].as<std::string>();
    } else if (result["data"].as<std::string>().size() > 0) {
        std::string dataset = result["data"].as<std::string>();
        this->path_filename_read = datasetMap[dataset].first;
        this->path_filename_ref = datasetMap[dataset].second;
    } else {
        std::cout << "Please specify dataset or gene and read path" << std::endl;
        exit(0);
    }

    this->path_kmer_matrix = result["kmer_matrix"].as<std::string>();
    
    this->k = result["kmer"].as<int>();
    this->target = result["target"].as<std::string>();
    this->progress_bar = result["bar"].as<bool>();
    this->verbose = result["verbose"].as<bool>();
    this->serial = result["serial"].as<bool>();
    this->path_output = result["out"].as<std::string>();
    this->debug_n_gene = result["debug_n_gene"].as<int>();
}

void PQuantParams::print() {
    std::cout << "path_filename_read = " << this->path_filename_read << std::endl;
    std::cout << "path_filename_ref = " << this->path_filename_ref << std::endl;
    std::cout << "path_output = " << this->path_output << std::endl;
    
    std::cout << "k = " << this->k << std::endl;
    std::cout << "target = " << this->target << std::endl;
    std::cout << "verbose = " << this->verbose << std::endl;
    std::cout << "serial = " << this->serial << std::endl;
    std::cout << "progress_bar = " << this->progress_bar << std::endl;
}
