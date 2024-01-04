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

using namespace std;
using namespace lbcrypto;

void printKmerTable(KmerTable &kmerTable, bool isRef);

void encryptReadKmer(KmerTable &kmerTableRead, long K, vector<Ciphertext<DCRTPoly>> &ct, CryptoContext<DCRTPoly> &cc, KeyPair<DCRTPoly> &keyPair, bool=true);

void multAll(KmerTable &kmerTableRead, long K, vector<Ciphertext<DCRTPoly>> &ct, CryptoContext<DCRTPoly> &cc, KeyPair<DCRTPoly> &keyPair);

void multCtxtByRef(vector<Ciphertext<DCRTPoly>> &ct_out, vector<Ciphertext<DCRTPoly>> &ct, Plaintext_1d &pt, CryptoContext<DCRTPoly> &cc);

#endif