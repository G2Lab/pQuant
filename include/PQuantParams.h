#ifndef __PQUANTPARAMS_H__
#define __PQUANTPARAMS_H__

#include <iostream>
#include <string>
#include <algorithm>


static  std::unordered_map<std::string, std::pair<const char*, const char*>> datasetMap = {
        {"five", { "../dataset/five_genes/five_gene_reads.fa", "../dataset/five_genes/five_gene_reference.fa"}},
        {"big", { "../dataset/big/reads_GM12878_20_genes.fa", "../dataset/big/exon_reference.fa"}},
        {"toy", { "../dataset/toy/reads_toy.fa", "../dataset/toy/refs_toy.fa"}},
        {"shorten", { "../dataset/five_genes/five_gene_reads_shorten.fa", "../dataset/five_genes/five_gene_reference.fa"}},
    };

class PQuantParams {
    public:
    long k;
    std::string filename_read;
    std::string filename_ref;
    bool verbose;
    bool progress_bar;
    std::string target;


    PQuantParams(std::string target, std::string filename_read, std::string filename_ref, long k, bool verbose, bool progress_bar) {
        this->k = k;
        this->filename_read = filename_read;
        this->filename_ref = filename_ref;
        this->verbose = verbose;
        this->progress_bar = progress_bar;
        this->target = target;
    }

    void print() {
        std::cout << "filename_read = " << filename_read << std::endl;
        std::cout << "filename_ref = " << filename_ref << std::endl;

        std::cout << "k = " << k << std::endl;
        std::cout << "target = " << target << std::endl;
        std::cout << "verbose = " << verbose << std::endl;
        std::cout << "progress_bar = " << progress_bar << std::endl;
    }
};


#endif // __PQUANTPARAMS_H__