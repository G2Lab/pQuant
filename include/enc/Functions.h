#ifndef FUNCTIONS_H
#define FUNCTIONS_H

#include <vector>
#include <string>
#include <map>
#include <iostream>
#include "Sequence.h"
#include "util/FileReader.h"
#include "util/Define.h"
#include "openfhe.h"
using namespace std;
using namespace lbcrypto;

void encryptReadKmer(KmerTable &kmerTableRead, long K, vector<Ciphertext<DCRTPoly>> &ct, CryptoContext<DCRTPoly> &cc, KeyPair<DCRTPoly> &keyPair);

#endif