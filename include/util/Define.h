#ifndef DEFINE_H_
#define DEFINE_H_

#include <vector>
#include "openfhe.h"

using namespace lbcrypto;

typedef vector<Plaintext> Plaintext_1d;
typedef vector<Plaintext_1d> Plaintext_2d;
typedef vector<Plaintext_2d> Plaintext_3d;
typedef vector<Plaintext_3d> Plaintext_4d;

typedef vector<Ciphertext<DCRTPoly>> Ciphertext_1d;
typedef vector<Ciphertext_1d> Ciphertext_2d;
typedef vector<Ciphertext_2d> Ciphertext_3d;
typedef vector<Ciphertext_3d> Ciphertext_4d;

typedef struct KmerTable {
    long K;
    size_t n_gene;
    size_t n_kmer_total;
    vector<string> geneNameIndex;
    map<size_t, map<size_t, size_t>> count; // gene, [kmer, num]
    map<size_t, size_t> countRead; // [kmer, num]
    map<size_t, float> entropy; // [kmer, entropy]
} KmerTable;

#endif