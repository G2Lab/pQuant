#ifndef Encode_h_
#define Encode_h_

#include "Sequence.h"
#include "openfhe.h"
using namespace lbcrypto;

int64_t encode_nt_to_num(string nt);

class EncodedSeq {
  private:
    string geneName;
    vector<vector<int64_t>> data;
    vector<vector<Plaintext>> encodedData;
    long n_seq;
    long n_ptxt;
    long n_block;
    long K;
    long B;

  public:
    EncodedSeq(Sequence seq, long K, long B, CryptoContext<DCRTPoly>& cc) {
        geneName = seq.getGeneName();
        n_block = K / B;
        n_seq = seq.getNumSeq();
        long n_slots = cc->GetCryptoParameters()->GetElementParams()->GetCyclotomicOrder() / 2;
        n_ptxt = (n_seq - 1) / n_slots + 1;

        encodedData.resize(n_ptxt);
        data.resize(n_seq, vector<int64_t>(n_block));

        for (int i = 0; i < n_seq; i++) {
            string kmer = seq.getSeq(i);
            for (int j = 0; j < n_block; j++) {
                string kmer_block = kmer.substr(j * B, B);
                long kmer_block_num = encode_nt_to_num(kmer_block);
                data[i][j] = kmer_block_num;
            }
        }

        for (int b = 0; b < n_block; b++) {
            for (int i = 0; i < n_ptxt; i++) {
                vector<int64_t> vec(n_slots);
                for (int j = 0; j < n_slots; j++) {
                    if (i * n_slots + j < n_seq) {
                        vec.push_back(data[i * n_slots + j][b]);
                    }
                }
                encodedData[i][b] = cc->MakePackedPlaintext(vec);
            }
        }
    }
};

#endif