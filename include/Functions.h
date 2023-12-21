#ifndef Functions_h_
#define Functions_h_

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include "util/FileReader.h"
#include "util/Define.h"
#include "util/util.h"
#include "PlainFunc.h"
#include "Sequence.h"
#include "openfhe.h"

using namespace std;
using namespace lbcrypto;

void print_sequences(vector<Sequence>& reads_seq, vector<Sequence>& refs_seq);

void encode_seq(CryptoContext<DCRTPoly> &cryptoContext, Sequence &seq, Plaintext_2d &plaintexts, long K, long B, long slot_count);

void encode_refs(CryptoContext<DCRTPoly> &cryptoContext, vector<Sequence>& refs_seq, Plaintext_3d &ref_plain, long K, long B, long slot_count);

void encode_read(CryptoContext<DCRTPoly> &cryptoContext, vector<Sequence> &reads_seq, Plaintext_4d &read_plain, long K, long B, long slot_count);
void encode_read2(CryptoContext<DCRTPoly> &cryptoContext, vector<Sequence> &reads_seq, Plaintext_4d &read_plain, long K, long B, long slot_count);

void encrypt_reads(CryptoContext<DCRTPoly> &cryptoContext, KeyPair<DCRTPoly> &keyPair, Plaintext_4d &read_plain,
                   Ciphertext_4d &ctxt_kmers, long K, long B, long slot_count);

void eqaulity_test(CryptoContext<DCRTPoly> &cryptoContext, KeyPair<DCRTPoly> &keyPair, Ciphertext_4d &ctxt_kmers, Plaintext_3d &ref_plain, Ciphertext_3d &ctxt_out, long K, long B, long log_poly_modulus_degree, vector<Sequence> &refs_seq, vector<Sequence> &reads_seq);

void decrypt_output(CryptoContext<DCRTPoly> &cryptoContext, KeyPair<DCRTPoly> &keyPair, Ciphertext_3d &ctxt_out, vector<Sequence> &refs_seq, vector<Sequence> &reads_seq, long K, long B);

#endif