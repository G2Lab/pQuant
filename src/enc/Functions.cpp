#include "enc/Functions.h"

void printKmerTable(KmerTable &kmerTable, bool isRef) {
    if (isRef) {
        cout << " KmerTable for reference" << endl;
        cout << "geneNameIndex.size() = " << kmerTable.geneNameIndex.size() << endl;
        cout << "count.size() = " << kmerTable.count.size() << endl;
        cout << "entropy.size() = " << kmerTable.entropy.size() << endl;

        for (auto &p : kmerTable.count) {
            cout << "kmer = " << p.first << " => count : ";
            for (int i = 0; i < p.second.size(); i++) {
                cout << p.second[i] << " ";
            }
            // print entropy
            cout << " => entropy = " << kmerTable.entropy[p.first];
            cout << endl;
        }
    } else {
        cout << " KmerTable for read" << endl;
        cout << "countRead.size() = " << kmerTable.countRead.size() << endl;
        for (auto &p : kmerTable.countRead) {
            cout << "kmer = " << p.first << " => count = " << p.second << endl;
        }
    }
}

void encryptReadKmer(KmerTable &kmerTableRead, long K, Ciphertext_1d &ct, CryptoContext<DCRTPoly> &cc, KeyPair<DCRTPoly> &keyPair, bool progress_bar) {
    // encrypt kmerTableRead.countRead into vector of ciphertexts (ct)
    // first, plain vector is created by vec[kmer(num)] = kmerTableRead.countRead[i]
    // then, plain vector is encrypted into ciphertexts

    // create plain vector
    long n_slots = cc->GetCryptoParameters()->GetElementParams()->GetCyclotomicOrder() / 2;
    long n_plain_vecs = pow(4, K) / n_slots;
    if (n_plain_vecs == 0)
        n_plain_vecs = 1;
    vector<vector<int64_t>> plain_vec(n_plain_vecs, vector<int64_t>(n_slots, 0));

    cout << "n_slots = " << n_slots << endl;
    cout << "n_plain_vecs = " << n_plain_vecs << endl;

    auto start_time = std::chrono::high_resolution_clock::now();

    for (auto &p : kmerTableRead.countRead) {
        long num_vec = p.first / n_slots;
        long num_slot = p.first % n_slots;
        plain_vec[num_vec][num_slot] = p.second;
    }

    for (int i = 0; i < n_plain_vecs; i++) {
        Plaintext plain = cc->MakeCoefPackedPlaintext(plain_vec[i]);
        Ciphertext<DCRTPoly> ciphertext = cc->Encrypt(keyPair.secretKey, plain);
        ct.push_back(ciphertext);

        // update progress bar
        if (progress_bar)
            print_progress_bar(i, n_plain_vecs, start_time);
    }

    auto end_time = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::seconds>(end_time - start_time).count();
    cout << "encryptReadKmer duration = " << duration << " s" << endl;
}

void encodeRefKmer(KmerTable &kmerTableRef, long K, Plaintext_2d &pt_ref, CryptoContext<DCRTPoly> &cc, KeyPair<DCRTPoly> &keyPair, bool progress_bar) {
    // encode KmerTableRef into vector of plaintexts (pt_ref)
    // first, plain vector is created by vec[kmer(num)] += kmerTableRef.count[i] for all # of genes
    // then, plain vector is encoded to plaintext

    // create plain vector
    long n_slots = cc->GetCryptoParameters()->GetElementParams()->GetCyclotomicOrder() / 2;
    long n_genes = kmerTableRef.geneNameIndex.size();
    long n_vec_per_gene = ceil(pow(4, K) / (double)n_slots);

    vector<vector<vector<int64_t>>> plain_vec(n_genes, vector<vector<int64_t>>(n_vec_per_gene, vector<int64_t>(n_slots, 0)));

    cout << "n_slots = " << n_slots << endl;
    cout << "n_genes = " << n_genes << endl;
    cout << "n_vec_per_gene = " << n_vec_per_gene << endl;

    auto start_time = std::chrono::high_resolution_clock::now();

    for (auto &p : kmerTableRef.count) {
        long num_vec = p.first / n_slots;
        long num_slot = p.first % n_slots;

        // parameter checkpoint
        assert(n_genes == p.second.size());
        for (int i = 0; i < p.second.size(); i++) {
            // encode reverse polynomial version
            if (num_slot == 0) {
                plain_vec[i][num_vec][num_slot] += p.second[i];
            } else {
                plain_vec[i][num_vec][n_slots - num_slot] -= p.second[i];
            }
        }
    }

    pt_ref.resize(n_genes);
    for (int g = 0; g < n_genes; g++) {
        pt_ref[g].resize(n_vec_per_gene);
        for (int i = 0; i < n_vec_per_gene; i++) {
            Plaintext plain = cc->MakeCoefPackedPlaintext(plain_vec[g][i]);
            pt_ref[g][i] = plain;

            // update progress bar
            print_progress_bar(g * n_vec_per_gene + i, n_genes * n_vec_per_gene, start_time);
        }
    }

    auto end_time = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::seconds>(end_time - start_time).count();
    cout << "encodeRefKmer duration = " << duration << " s" << endl;
}

void multCtxtByRef(Ciphertext_2d &ct_out, Ciphertext_1d &ct, Plaintext_2d &pt_ref, CryptoContext<DCRTPoly> &cc) {
    int g_genes = pt_ref.size();
    ct_out.resize(g_genes);
    for (int g = 0; g < g_genes; g++) {
        assert(ct.size() == pt_ref[g].size());
        // multiply ciphertexts in ct with plaintexts in pt
        // result is stored in ct_out
        // ct_out[i] = ct[i] * pt[i]

        ct_out[g].resize(ct.size());
        for (int i = 0; i < ct.size(); i++) {
            ct_out[g][i] = cc->EvalMult(ct[i], pt_ref[g][i]);
        }
    }
}

void decCtxtOut(Plaintext_1d &pt_out, Ciphertext_1d &ct_out, CryptoContext<DCRTPoly> &cc, KeyPair<DCRTPoly> &keyPair) {
    auto start_time = std::chrono::high_resolution_clock::now();

    for (int i = 0; i < ct_out.size(); i++) {
        Plaintext plain_out;
        cc->Decrypt(keyPair.secretKey, ct_out[i], &plain_out);
        pt_out.push_back(plain_out);

        // update progress bar
        print_progress_bar(i, ct_out.size(), start_time);
    }
}

void sumUpCtxt(Ciphertext<DCRTPoly> &ct, Ciphertext_1d &ct_vec, CryptoContext<DCRTPoly> &cc) {
    auto start_time = std::chrono::high_resolution_clock::now();

    ct = ct_vec[0];
    for (int i = 1; i < ct_vec.size(); i++) {
        ct = cc->EvalAdd(ct, ct_vec[i]);

        // update progress bar
        print_progress_bar(i, ct_vec.size(), start_time);
    }
}