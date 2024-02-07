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
    
    this->data = result["data"].as<std::string>();
    this->thres = result["thres"].as<float>();
    this->k = result["kmer"].as<int>();
    this->target = result["target"].as<std::string>();
    this->progress_bar = result["bar"].as<bool>();
    this->verbose = result["verbose"].as<bool>();
    this->debug_n_gene = result["debug_n_gene"].as<int>();
    this->gene_start = result["gene_start"].as<int>();
    this->gene_end = result["gene_end"].as<int>();
    this->memory = result["divide_encode_mult"].as<bool>();
}

void PQuantParams::print() {
    std::cout << " == file paths ==" << std::endl;
    std::cout << "data name = " << this->data << std::endl;
    std::cout << "path_filename_read = " << this->path_filename_read << std::endl;
    std::cout << "path_filename_ref = " << this->path_filename_ref << std::endl;
    std::cout << "path_kmer_matrix = " << this->path_kmer_matrix << std::endl;
    std::cout << std::endl;
    std::cout << " == default parameters ==" << std::endl;
    std::cout << "k = " << this->k << std::endl;
    std::cout << "thres = " << this->thres << std::endl;
    std::cout << "target = " << this->target << std::endl;
    std::cout << std::endl;
    std::cout << " == conditional parameters ==" << std::endl;
    std::cout << "debug_n_gene = " << this->debug_n_gene << std::endl;
    std::cout << "gene_start = " << this->gene_start << std::endl;
    std::cout << "gene_end = " << this->gene_end << std::endl;
    std::cout << std::endl;
    std::cout << " == flags ==" << std::endl;
    std::cout << "verbose = " << this->verbose << std::endl;
    std::cout << "divide_encode_mult (memory) = " << this->memory << std::endl;
    std::cout << "print_progress_bar = " << this->progress_bar << std::endl;
}
