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


#endif