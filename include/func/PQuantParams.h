#ifndef __PQUANTPARAMS_H__
#define __PQUANTPARAMS_H__

#include <iostream>
#include <string>
#include <algorithm>
#include <unordered_map>
#include "cxxopts.hpp"

static  std::unordered_map<std::string, std::pair<const char*, const char*>> datasetMap = {
        {"five", { "../dataset/five_genes/five_gene_reads.fa", "../dataset/five_genes/five_gene_reference.fa"}},
        {"big", { "../dataset/big/reads_GM12878_20_genes.fa", "../dataset/big/exon_reference.fa"}},
        {"toy", { "../dataset/toy/reads_toy.fa", "../dataset/toy/refs_toy.fa"}},
        {"shorten", { "../dataset/five_genes/five_gene_reads_shorten.fa", "../dataset/five_genes/five_gene_reference.fa"}},
        {"5k", { "/gpfs/commons/groups/gursoy_lab/cwalker/projects/pquant/workflow/data/test_fastqs/5k_random_protein_coding_genes.genes_only.fq", "/gpfs/commons/groups/gursoy_lab/cwalker/projects/pquant/workflow/data/reference/pquant/5k_random_protein_coding_genes.combined_exons.exons.fa"}}
    };

class PQuantParams {
    public:
    long k;
    float thres;
    std::string data;
    std::string filename_read;
    std::string filename_ref;
    std::string foldername_kmer;
    std::string filename_kmerTable;
    std::string filename_kmerList;
    std::string foldername_BFV;
    std::string foldername_ctxtread;
    std::string foldername_ctxtout;
    std::string target;

    //debug
    bool verbose;
    bool json_format;
    int debug_n_gene;
    int gene_start;
    int gene_end;
    bool progress_bar;
    bool operate_then_serialize;

    PQuantParams(cxxopts::ParseResult &result);

    void print();
};


#endif // __PQUANTPARAMS_H__