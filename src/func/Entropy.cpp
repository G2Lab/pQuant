#include "func/Entropy.h"

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

string sequenceReverseCompliment(string nt) {
    string ans = "";
    for (size_t i = 0; i < nt.size(); i++) {
        if (nt[i] == 'A') {
            ans = 'T' + ans;
        } else if (nt[i] == 'C') {
            ans = 'G' + ans;
        } else if (nt[i] == 'G') {
            ans = 'C' + ans;
        } else if (nt[i] == 'T') {
            ans = 'A' + ans;
        } else {
            return "";
        }
    }
    return ans;
}

void computeKmerTable(vector<Sequence>& gene, long K, KmerTable &kmerTable) {
    kmerTable.K = K;
    kmerTable.n_gene = gene.size();
    if (kmerTable.geneNameIndex.size() == 0) {
        // save gene list into geneNameIndex
        for (size_t i = 0; i < gene.size(); i++) {
            kmerTable.geneNameIndex.push_back(gene[i].getGeneName());
        }
    }
    
    auto start_time = std::chrono::high_resolution_clock::now();
    set<int> kmer_occurance;
    // map<size_t, map<size_t, size_t>> kmer_based_table;
    // count kmers in gene
    for (size_t i = 0; i < gene.size(); i++) {
        // gene i
        map<size_t, size_t> gene_kmer_count;
        for (size_t j = 0; j < static_cast<size_t>(gene[i].getNumSeq()); j++) {
            string seq = gene[i].getSeq(j);
            if (seq.size() < static_cast<size_t>(K)) {
                continue;
            }

            // original order
            for (size_t k = 0; k < seq.size() - K + 1; k++) {
                string kmer = seq.substr(k, K);
                long num = convertKmerToNum(kmer);
                if (num == -1) {
                    continue;
                }
                // if (gene_kmer_count.find(num) == gene_kmer_count.end()) {
                //     gene_kmer_count.insert(make_pair(num, 0));
                //     cout << "num = " << num << "is inserted" << endl;
                // }
                gene_kmer_count[num] += 1;
                // kmer_based_table[num][i] += 1;
                if (kmer_occurance.find(num) == kmer_occurance.end()) {
                    kmer_occurance.insert(num);
                }
            }

            // reverse compliment
            string seq_rc = sequenceReverseCompliment(seq);
            for (size_t k = 0; k < seq_rc.size() - K + 1; k++) {
                string kmer = seq_rc.substr(k, K);
                long num = convertKmerToNum(kmer);
                if (num == -1) {
                    continue;
                }
                gene_kmer_count[num] += 1;
                if (kmer_occurance.find(num) == kmer_occurance.end()) {
                    kmer_occurance.insert(num);
                }
            }
        }

        // concatenations
        for (size_t j = 0; j < static_cast<size_t>(gene[i].getNumSeq()); j++) {
            for (size_t k = j + 1; k < static_cast<size_t>(gene[i].getNumSeq()); k++) {
                string seq_j = gene[i].getSeq(j);
                string seq_k = gene[i].getSeq(k);
                string seq_combined;
                if (seq_j.size() >= static_cast<size_t>(K) && seq_k.size() >= static_cast<size_t>(K)) {
                    seq_combined = seq_j.substr(seq_j.size() - K + 1, K - 1) + seq_k.substr(0, K - 1);
                } else if (seq_j.size() >= static_cast<size_t>(K)) {
                    seq_combined = seq_j.substr(seq_j.size() - K + 1, K - 1) + seq_k;
                } else if (seq_k.size() >= static_cast<size_t>(K)) {
                    seq_combined = seq_j + seq_k.substr(0, K - 1);
                } else {
                    seq_combined = seq_j + seq_k;
                }
                if (seq_combined.size() < static_cast<size_t>(K)) {
                    continue;
                }
                // original order
                for (size_t l = 0; l < seq_combined.size() - K + 1; l++) {
                    string kmer = seq_combined.substr(l, K);
                    long num = convertKmerToNum(kmer);
                    if (num == -1) {
                        continue;
                    }
                    gene_kmer_count[num] += 1;
                    if (kmer_occurance.find(num) == kmer_occurance.end()) {
                        kmer_occurance.insert(num);
                    }
                }

                // reverse compliement
                string seq_combined_rc = sequenceReverseCompliment(seq_combined);
                for (size_t l = 0; l < seq_combined_rc.size() - K + 1; l++) {
                    string kmer = seq_combined_rc.substr(l, K);
                    long num = convertKmerToNum(kmer);
                    if (num == -1) {
                        continue;
                    }
                    gene_kmer_count[num] += 1;
                    if (kmer_occurance.find(num) == kmer_occurance.end()) {
                        kmer_occurance.insert(num);
                    }
                }
            }
        }

        // save gene_kmer_count into kmerTable.count
        kmerTable.count.insert(make_pair(i, gene_kmer_count));
        gene_kmer_count.clear();

        print_progress_bar("countKmersPerGene", i, gene.size(), start_time);
    }
    kmerTable.n_kmer_total = kmer_occurance.size();
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