#include "util/KmerTable.h"

void KmerTable::print() {
    if (isRef) {
        std::cout << "geneNames" << std::endl;
        for (auto &p : geneNameIndex) {
            std::cout << p << ", ";
        }
        std::cout << std::endl << std::endl;

        std::cout << " KmerTable for reference" << std::endl;
        std::cout << "n_gene = " << n_gene << std::endl;
        std::cout << "n_kmer_total = " << n_kmer_total << std::endl;

        for (const auto& geneEntry : count) {
            std::cout << "Gene: " << geneEntry.first << " (" << geneEntry.second.size() << " kmers)" << std::endl;
        
            for (const auto& kmerEntry : geneEntry.second) {
                 std::cout << "(" << kmerEntry.first << ", " << kmerEntry.second << ") ";
            }
             std::cout << std::endl;
        }
    } else {
        std::cout << " KmerTable for read" << std::endl;
        std::cout << "countRead.size() = " << countRead.size() << std::endl;
        for (auto &p : countRead) {
             std::cout << "(" << p.first << ", " << p.second << ") ";
        }
    }
    std::cout << std::endl;
}


void KmerTable::save(std::string filename) {
    // ofstream file(filename);
    // // if (!file.is_open()) {
    // //     throw runtime_error("Error: Failed to open file " + filename);
    // // }

    // file << "gene_name_vec: [";
    // for (size_t i = 0; i < geneNameIndex.size(); i++) {
    //     file << "\"" << geneNameIndex[i] << "\"";
    //     if (i < geneNameIndex.size() - 1) {
    //         file << ", ";
    //     }
    // }
    // file << "]" << std::endl;

    // file << "kmer_matrix: {" << std::endl;
    // for (size_t i = 0; i < n_gene; i++) {
    //     file << "  \"" << geneNameIndex[i] << "\": {" << std::endl;
    //     for (size_t j = 0; j < n_gene; j++) {
    //         file << "    \"" << geneNameIndex[j] << "\": " << count[i][j];
    //         if (j < n_gene - 1) {
    //             file << ",";
    //         }
    //         file << std::endl;
    //     }
    //     file << "  }";
    //     if (i < n_gene - 1) {
    //         file << ",";
    //     }
    //     file << std::endl;
    // }
    // file << "}" << std::endl;

    // file << "kmer_entropy: {" << std::endl;
    // for (size_t i = 0; i < n_kmer_total; i++) {
    //     file << "  \"" << i << "\": " << entropy.at(i);
    //     if (i < n_kmer_total - 1) {
    //         file << ",";
    //     }
    //     file << std::endl;
    // }
    // file << "}" << std::endl;

    // file.close();
}