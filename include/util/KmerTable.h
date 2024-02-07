#ifndef KMERTABLE_H_
#define KMERTABLE_H_

#include <vector>
#include <map>
#include <iostream>
#include <fstream>

class KmerTable {
public:
    long K;
    size_t n_gene;
    size_t n_kmer_total;
    std::vector<std::string> geneNameIndex;
    std::map<size_t, std::map<size_t, size_t>> count; // gene, [kmer, num]
    std::map<size_t, size_t> countRead; // [kmer, num]
    std::map<size_t, float> entropy; // [kmer, entropy]
    bool isRef;

    KmerTable() {
        K = 0;
        n_gene = 0;
        n_kmer_total = 0;
        isRef = true;
    }

    void print();
    void save(std::string filename);
};

#endif