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

// header files needed for serialization
#include "ciphertext-ser.h"
#include "cryptocontext-ser.h"
#include "key/key-ser.h"
#include "scheme/bfvrns/bfvrns-ser.h"

#include <filesystem>
namespace fs = std::filesystem;

using namespace std;
using namespace lbcrypto;

void printKmerTable(KmerTable &kmerTable, bool isRef);
void printKmerTable2(KmerTable2 &kmerTable, bool isRef);

void encryptReadKmer(KmerTable &kmerTableRead, long K, Ciphertext_1d &ct, CryptoContext<DCRTPoly> &cc, KeyPair<DCRTPoly> &keyPair, bool=true);

void encodeRefKmer(KmerTable &kmerTableRef, long K, Plaintext_2d &pt_ref, CryptoContext<DCRTPoly> &cc, KeyPair<DCRTPoly> &keyPair, bool progress_bar);

void multCtxtByRef(Ciphertext_2d &ct_out, Ciphertext_1d &ct, Plaintext_2d &pt_ref, CryptoContext<DCRTPoly> &cc);

void multCtxtByKmerTableRef(Ciphertext_2d &ct_out, Ciphertext_1d &ct, KmerTable kmerTableRef, long K, CryptoContext<DCRTPoly> &cc);
void multCtxtByKmerTableRef2(Ciphertext_2d &ct_out, Ciphertext_1d &ct, KmerTable2 kmerTableRef, long K, CryptoContext<DCRTPoly> &cc);
void multCtxtByKmerTableRefFromSerial(Ciphertext_2d &ct_out, Ciphertext_1d &ct, KmerTable2 kmerTableRef, long K, CryptoContext<DCRTPoly> &cc, string out_path);
void decCtxtOut(Plaintext_1d &pt_out, Ciphertext_1d &ct_out, CryptoContext<DCRTPoly> &cc, KeyPair<DCRTPoly> &keyPair);

void sumUpCtxt(Ciphertext<DCRTPoly> &ct, Ciphertext_1d &ct_vec, CryptoContext<DCRTPoly> &cc);
void sumUpCtxtFromSerial(Ciphertext_1d &ct_out, size_t n_gene, size_t n_ctxt, CryptoContext<DCRTPoly> &cc, string out_path);
#endif