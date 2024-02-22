#ifndef KMERTABLE_H_
#define KMERTABLE_H_

#include <vector>
#include <map>
#include <iostream>
#include <fstream>
#include "util/util.h"

long convertKmerToNum(string kmer);

string sequenceReverseCompliment(string nt);

class KmerTable {
public:
    long K;
    float thres;
    size_t n_gene;
    size_t n_kmer_total;
    std::vector<std::string> geneNameIndex;
    std::map<size_t, std::map<size_t, size_t>> count; // gene, [kmer, num]
    std::map<size_t, std::map<size_t, size_t>> tableRef; // kmer, [(gene, num)] (entropy thresholded)
    std::map<size_t, size_t> countRead; // [kmer, num]
    std::map<size_t, float> entropy; // [kmer, entropy]
    bool isRef;

    KmerTable() {
        K = 0;
        n_gene = 0;
        n_kmer_total = 0;
        isRef = true;
    }

    KmerTable(vector<Sequence>& gene, PQuantParams &param, bool isRef_);

    bool operator==(const KmerTable& other) const;

    void print();
    void save(const std::string& filename);
    void saveKmerList(const std::string& filename);
    void load(const std::string& filename);
    
    // Function to save KmerTable data to a binary file
    void saveBinary(const std::string& filename);

    // Function to load KmerTable data from a binary file
    void loadBinary(const std::string& filename);

    void saveKmerListBinary(const std::string& filename);

    void filterGenes(size_t start, size_t end);
};

void loadKmerList(const std::string& filename, vector<size_t>& kmerList);
void loadKmerListBinary(const std::string& filename, vector<size_t>& kmerList);
#endif