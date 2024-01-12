#include "Entropy.h"

long convertKmerToNum(string kmer) {
    long num = 0;
    for (size_t i = 0; i < kmer.size(); i++) {
        if (kmer[i] == 'A') {
            num = num * 4 + 0;
        } else if (kmer[i] == 'C') {
            num = num * 4 + 1;
        } else if (kmer[i] == 'G') {
            num = num * 4 + 2;
        } else if (kmer[i] == 'T') {
            num = num * 4 + 3;
        } else {
            return -1;
        }
    }
    return num;
}

void computeKmerTable(vector<Sequence>& gene, long K, KmerTable &kmerTable) {
    kmerTable.K = K;
    if (kmerTable.geneNameIndex.size() == 0) {
        // save gene list into geneNameIndex
        for (size_t i = 0; i < gene.size(); i++) {
            kmerTable.geneNameIndex.push_back(gene[i].getGeneName());
        }
    }

    auto start_time = std::chrono::high_resolution_clock::now();
    // count kmers in gene
    for (size_t i = 0; i < gene.size(); i++) {
        for (int j = 0; j < gene[i].getNumSeq(); j++) {
            string seq = gene[i].getSeq(j);
            if (seq.size() < static_cast<size_t>(K)) {
                continue;
            }
            for (size_t k = 0; k < seq.size() - K + 1; k++) {
                string kmer = seq.substr(k, K);
                long num = convertKmerToNum(kmer);
                // cout << "kmer = " << kmer << " => num = " << num << endl;
                if (num == -1) {
                    continue;
                }
                if (kmerTable.count.find(num) == kmerTable.count.end()) {
                    kmerTable.count[num] = vector<long>(kmerTable.geneNameIndex.size(), 0);
                }
                kmerTable.count[num][i] += 1;
            }
        }
        print_progress_bar("countKmersPerGene", i, gene.size(), start_time);
    }

    // compute entropy
    start_time = std::chrono::high_resolution_clock::now();
    for (auto &kmerCount : kmerTable.count) {
        long sum = 0;
        for (size_t i = 0; i < kmerCount.second.size(); i++) {
            sum += kmerCount.second[i];
        }
        float entropy = 0;
        for (size_t i = 0; i < kmerCount.second.size(); i++) {
            float p = (float)kmerCount.second[i] / sum;
            if (p > 0) {
                entropy += p * log2(p);
            }
        }
        entropy = -entropy;
        kmerTable.entropy[kmerCount.first] = entropy;
        // print_progress_bar("ComputeEntropy", kmerCount.first, kmerTable.count.size(), start_time);
    }
}

void computeKmerTableForRead(vector<Sequence>& read, long K, KmerTable &kmerTable) {
    // for read, sum up all values and save them into countRead
    kmerTable.K = K;

    // count kmers in gene
    for (size_t i = 0; i < read.size(); i++) {
        for (int j = 0; j < read[i].getNumSeq(); j++) {
            string seq = read[i].getSeq(j);
            if (seq.size() < static_cast<size_t>(K)) {
                continue;
            }
            for (size_t k = 0; k < seq.size() - K + 1; k++) {
                string kmer = seq.substr(k, K);
                long num = convertKmerToNum(kmer);
                // cout << "kmer = " << kmer << " => num = " << num << endl;
                if (num == -1) {
                    continue;
                }
                if (kmerTable.countRead.find(num) == kmerTable.countRead.end()) {
                    kmerTable.countRead[num] = 0;
                }
                kmerTable.countRead[num] += 1;
            }
        }
    }

}