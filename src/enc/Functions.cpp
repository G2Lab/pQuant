#include "enc/Functions.h"

void encryptReadKmer(KmerTable &kmerTableRead, long K, vector<Ciphertext<DCRTPoly>> &ct, CryptoContext<DCRTPoly> &cc, KeyPair<DCRTPoly> &keyPair) {
    // encrypt kmerTableRead.countRead into vector of ciphertexts (ct)
    // first, plain vector is created by vec[kmer(num)] = kmerTableRead.countRead[i]
    // then, plain vector is encrypted into ciphertexts

    // create plain vector
    long n_slots = cc->GetCryptoParameters()->GetElementParams()->GetCyclotomicOrder() / 2;
    long n_plain_vecs = pow(4, K) / n_slots;
    if (n_plain_vecs == 0) n_plain_vecs = 1;
    vector<vector<int64_t>> plain_vec(n_plain_vecs, vector<int64_t>(n_slots, 0));

    cout << "n_slots = " << n_slots << endl;
    cout << "n_plain_vecs = " << n_plain_vecs << endl;

    for (auto &p : kmerTableRead.countRead) {
        long num_vec = p.first / n_slots;
        long num_slot = p.first % n_slots;
        // cout << "num = " << p.first << " : num_vec = " << num_vec << ", num_slot = " << num_slot << endl;
        plain_vec[num_vec][num_slot] = p.second;
    }

    // // print to check plain_vec
    // for (int i = 0; i < n_plain_vecs; i++) {
    //     cout << "plain_vec[" << i << "] = " << endl;
    //     for (int j = 0; j < n_slots; j++) {
    //         cout << j << " : " << plain_vec[i][j] << endl;
    //     }
    // }
    
float progress = 0.0;
std::cout << std::endl;
    // encrypt plain vector into ciphertexts
    for (int i = 0; i < n_plain_vecs; i++) {
        Plaintext plain = cc->MakeCoefPackedPlaintext(plain_vec[i]);
        Ciphertext<DCRTPoly> ciphertext = cc->Encrypt(keyPair.secretKey, plain);
        ct.push_back(ciphertext);
        int barWidth = 70;

        std::cout << "[";
        int pos = barWidth * progress;
        for (int i = 0; i < barWidth; ++i) {
            if (i < pos) std::cout << "=";
            else if (i == pos) std::cout << ">";
            else std::cout << " ";
        }
        std::cout << "] " << int(progress * 100.0) << " %\r";
        std::cout.flush();

        progress += (double) 1. / n_plain_vecs; // for demonstration only      
    }
    cout << "done" << endl;
}