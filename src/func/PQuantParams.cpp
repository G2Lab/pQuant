#include "func/PQuantParams.h"

PQuantParams::PQuantParams(cxxopts::ParseResult &result) {
    // exception for test functions
    if (result["target"].as<std::string>() == "bench") {
        this->target = "bench";
        return;
    }
    if (result["gene"].as<std::string>().size() == 0 && result["read"].as<std::string>().size() == 0 && result["data"].as<std::string>().size() == 0) {
        std::cerr << "Please specify dataset or gene and read path" << std::endl;
        exit(0);
    }
    if (result["gene"].as<std::string>().size() == 0 && result["read"].as<std::string>().size() == 0 && result["data"].as<std::string>().size() > 0) {
        // output error if result["data"] is not in datasetMap key list
        this->data = result["data"].as<std::string>();
        if (datasetMap.find(this->data) == datasetMap.end()) {
            std::cerr << "Invalid dataset: " << this->data << std::endl;
            std::cerr << "add gene & read paths, or try valid dataset name" << std::endl;
            exit(0);
        }
        this->filename_read = datasetMap[this->data].first;
        this->filename_ref = datasetMap[this->data].second;
    } else {
        if (result["gene"].as<std::string>().size() > 0) {
            this->filename_ref = result["gene"].as<std::string>();    
        }
        if (result["read"].as<std::string>().size() > 0) {
            this->filename_read = result["read"].as<std::string>();
        }
        if (result["data"].as<std::string>().size() > 0) {
            this->data = result["data"].as<std::string>();
        } else {
            this->data = "data";
        }
    }
    
    
    this->thres = result["thres"].as<float>();
    this->k = result["kmer"].as<int>();
    this->target = result["target"].as<std::string>();
    this->progress_bar = result["bar"].as<bool>();
    this->verbose = result["verbose"].as<bool>();
    this->debug_n_gene = result["debug_n_gene"].as<int>();
    this->gene_start = result["gene_start"].as<int>();
    this->gene_end = result["gene_end"].as<int>();
    this->batch_num = result["batch_num"].as<int>();
    this->batch_num_total = result["batch_num_total"].as<int>();
    this->json_format = result["json"].as<bool>();
    this->operate_then_serialize = result["operate_then_serialize"].as<bool>();

    this->foldername_kmer = result["kmer_folder"].as<std::string>();
    if (filename_kmerTable.size() == 0 && filename_kmerList.size() == 0 && foldername_kmer.size() > 0) {
        filename_kmerTable = foldername_kmer + "/kmertable_" + data + ".bin";
        filename_kmerList =  foldername_kmer + "/kmerlist_" + data + ".bin";
    }
    foldername_BFV = result["bfv_folder"].as<std::string>();
    foldername_ctxtread = result["ctxt_read_folder"].as<std::string>();
    foldername_ctxtout = result["ctxt_out_folder"].as<std::string>();
    if (foldername_ctxtread.size() == 0) {
        foldername_ctxtread = foldername_BFV + "/ctxtread";
    }
    if (foldername_ctxtout.size() == 0) {
        foldername_ctxtout = foldername_BFV + "/ctxtout";
    }
}

void PQuantParams::print() {
    std::cout << " == file paths ==" << std::endl;
    std::cout << "data name = " << this->data << std::endl;
    std::cout << "filename_read = " << this->filename_read << std::endl;
    std::cout << "filename_ref = " << this->filename_ref << std::endl;
    std::cout << "filename_kmerTable = " << this->filename_kmerTable << std::endl;
    std::cout << "filename_kmerList = " << this->filename_kmerList << std::endl;
    std::cout << "foldername_BFV = " << this->foldername_BFV << std::endl;
    std::cout << "foldername_ctxtread = " << this->foldername_ctxtread << std::endl;
    std::cout << "foldername_ctxtout = " << this->foldername_ctxtout << std::endl;
    
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
    std::cout << "batch_num = " << this->batch_num << std::endl;
    std::cout << "batch_num_total = " << this->batch_num_total << std::endl;
    std::cout << std::endl;
    std::cout << " == flags ==" << std::endl;
    std::cout << "verbose = " << this->verbose << std::endl;
    std::cout << "json_format = " << this->json_format << std::endl;
    std::cout << "print_progress_bar = " << this->progress_bar << std::endl;
    std::cout << "operate_then_serialize = " << this->operate_then_serialize << std::endl;
}
