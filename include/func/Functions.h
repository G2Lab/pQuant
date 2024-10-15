#ifndef FUNCTIONS_H
#define FUNCTIONS_H

#include <vector>
#include <string>
#include <map>
#include <iostream>
#include "Sequence.h"
#include "util/FileReader.h"
#include "util/Define.h"
#include "util/util.h"
#include "openfhe.h"
#include <chrono>
#include <iomanip>
#include <cmath>
#include <numeric>

// header files needed for serialization
#include "ciphertext-ser.h"
#include "cryptocontext-ser.h"
#include "key/key-ser.h"
#include "scheme/bfvrns/bfvrns-ser.h"
#include "func/PQuantParams.h"

using namespace std;
using namespace lbcrypto;

void printKmerTable(KmerTable &kmerTable, bool isRef);

void encryptReadKmer(KmerTable &kmerTableRead, Ciphertext_1d &ct, CryptoContext<DCRTPoly> &cc, KeyPair<DCRTPoly> &keyPair, PQuantParams &param);

void encryptReadSparse(Ciphertext_1d &ct, KmerTable &kmerTableRead, KmerTable &kmerTableRef, CryptoContext<DCRTPoly> &cc, KeyPair<DCRTPoly> &keyPair, PQuantParams &param);

void encodeRefKmer(KmerTable &kmerTableRef, Plaintext_2d &pt_ref, CryptoContext<DCRTPoly> &cc, KeyPair<DCRTPoly> &keyPair, PQuantParams &param);

void multCtxtByRef(Ciphertext_2d &ct_out, Ciphertext_1d &ct, Plaintext_2d &pt_ref, CryptoContext<DCRTPoly> &cc);

void multCtxtByKmerTableRef(Ciphertext_1d &ct_out, Ciphertext_1d &ct, KmerTable kmerTableRef, CryptoContext<DCRTPoly> &cc, PQuantParams &param);

void multCtxtByEncodedRef(Ciphertext_1d &ct_out, Ciphertext_1d &ct, Plaintext_2d &pt_ref, CryptoContext<DCRTPoly> &cc, PQuantParams &param);

void decCtxtOut(Plaintext_1d &pt_out, Ciphertext_1d &ct_out, CryptoContext<DCRTPoly> &cc, KeyPair<DCRTPoly> &keyPair);

void sumUpCtxt(Ciphertext<DCRTPoly> &ct, Ciphertext_1d &ct_vec, CryptoContext<DCRTPoly> &cc);
void sumUpCtxtFromSerial(Ciphertext_1d &ct_out, size_t n_gene, size_t n_ctxt, CryptoContext<DCRTPoly> &cc, string path_output);

void generateSimulatedReads(vector<Sequence>& gene_seq, size_t read_len, size_t num_reads, vector<Sequence>& simulated_reads);
#endif