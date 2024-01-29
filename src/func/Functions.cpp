#include "func/Functions.h"

void printKmerTable(KmerTable &kmerTable, bool isRef) {
    if (isRef) {
        std::cout << "geneNames" << std::endl;
        for (auto &p : kmerTable.geneNameIndex) {
            std::cout << p << ", ";
        }
        std::cout << std::endl << std::endl;

        std::cout << " KmerTable for reference" << std::endl;
        std::cout << "n_gene = " << kmerTable.n_gene << std::endl;
        std::cout << "n_kmer_total = " << kmerTable.n_kmer_total << std::endl;

        for (const auto& geneEntry : kmerTable.count) {
            std::cout << "Gene: " << geneEntry.first << " (" << geneEntry.second.size() << " kmers)" << std::endl;
        
            for (const auto& kmerEntry : geneEntry.second) {
                 std::cout << "(" << kmerEntry.first << ", " << kmerEntry.second << ") ";
            }
             std::cout << std::endl;
        }

    } else {
        std::cout << " KmerTable for read" << std::endl;
        std::cout << "countRead.size() = " << kmerTable.countRead.size() << std::endl;
        for (auto &p : kmerTable.countRead) {
             std::cout << "(" << p.first << ", " << p.second << ") ";
        }
    }
}

void encryptReadKmer(KmerTable &kmerTableRead, Ciphertext_1d &ct, CryptoContext<DCRTPoly> &cc, KeyPair<DCRTPoly> &keyPair, PQuantParams &param) {
    int K = param.k;
    // encrypt kmerTableRead.countRead into vector of ciphertexts (ct)
    // first, plain vector is created by vec[kmer(num)] = kmerTableRead.countRead[i]
    // then, plain vector is encrypted into ciphertexts

    // create plain vector
    long n_slots = cc->GetCryptoParameters()->GetElementParams()->GetCyclotomicOrder() / 2;
    long n_plain_vecs = pow(4, K) / n_slots;
    if (n_plain_vecs == 0)
        n_plain_vecs = 1;
    vector<vector<int64_t>> plain_vec(n_plain_vecs, vector<int64_t>(n_slots, 0));

    std::cout << "n_slots = " << n_slots << std::endl;
    std::cout << "n_plain_vecs = " << n_plain_vecs << std::endl;

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
        if (param.progress_bar)
            print_progress_bar("EncryptReadKmer", i, n_plain_vecs, start_time);
    }

    auto end_time = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::seconds>(end_time - start_time).count();
    std::cout << "encryptReadKmer duration = " << duration << " s" << std::endl;
}

void encryptReadSparse(Ciphertext_1d &ct, KmerTable &kmerTableRead, KmerTable &kmerTableRef, CryptoContext<DCRTPoly> &cc, KeyPair<DCRTPoly> &keyPair, PQuantParams &param) {
    // encrypt kmerTableRead.countRead into vector of ciphertexts (ct)
    // first, plain vector is created by vec[kmer(num)] = kmerTableRead.countRead[i]
    // then, plain vector is encrypted into ciphertexts

    // create plain vector
    long n_slots = cc->GetCryptoParameters()->GetElementParams()->GetCyclotomicOrder() / 2;
    long n_plain_vecs = (kmerTableRef.n_kmer_total - 1) / n_slots + 1;
    vector<vector<int64_t>> plain_vec(n_plain_vecs, vector<int64_t>(n_slots, 0));

    std::cout << "n_slots = " << n_slots << std::endl;
    std::cout << "n_plain_vecs = " << n_plain_vecs << std::endl;

    auto start_time = std::chrono::high_resolution_clock::now();

    // set kmer index vector
    vector<size_t> kmer_index;
    for (auto &p : kmerTableRef.entropy) {
        kmer_index.push_back(p.first);
    }

    cout << "kmer_index.size() = " << kmer_index.size() << endl;

    // encode kmerTableRead.countRead based on kmer_index
    int count = 0;
    for (auto &p : kmerTableRead.countRead) {
        auto it = std::lower_bound(kmer_index.begin(), kmer_index.end(), p.first);
        if (it != kmer_index.end() && *it == p.first) {
            size_t index = std::distance(kmer_index.begin(), it);
            size_t num_vec = index / n_slots;
            size_t num_slot = index % n_slots;
            plain_vec[num_vec][num_slot] = p.second;
        }
        if (param.progress_bar) 
            print_progress_bar("encodeRead", count, kmerTableRead.countRead.size(), start_time);
        count += 1;
    }

    for (int i = 0; i < n_plain_vecs; i++) {
        Plaintext plain = cc->MakeCoefPackedPlaintext(plain_vec[i]);
        Ciphertext<DCRTPoly> ciphertext = cc->Encrypt(keyPair.secretKey, plain);
        ct.push_back(ciphertext);

        // update progress bar
        if (param.progress_bar)
            print_progress_bar("EncryptReadKmer", i, n_plain_vecs, start_time);
    }

    auto end_time = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::seconds>(end_time - start_time).count();
    std::cout << "encryptReadKmer duration = " << duration << " s" << std::endl;
}

// void encodeRefKmer(KmerTable &kmerTableRef, long K, Plaintext_2d &pt_ref, CryptoContext<DCRTPoly> &cc, KeyPair<DCRTPoly> &keyPair, bool progress_bar) {
//     // encode KmerTableRef into vector of plaintexts (pt_ref)
//     // first, plain vector is created by vec[kmer(num)] += kmerTableRef.count[i] for all # of genes
//     // then, plain vector is encoded to plaintext

//     // create plain vector
//     long n_slots = cc->GetCryptoParameters()->GetElementParams()->GetCyclotomicOrder() / 2;
//     long n_genes = kmerTableRef.geneNameIndex.size();
//     long n_vec_per_gene = ceil(pow(4, K) / (double)n_slots);

//     vector<vector<vector<int64_t>>> plain_vec(n_genes, vector<vector<int64_t>>(n_vec_per_gene, vector<int64_t>(n_slots, 0)));

//     std::cout << "n_slots = " << n_slots << std::endl;
//     std::cout << "n_genes = " << n_genes << std::endl;
//     std::cout << "n_vec_per_gene = " << n_vec_per_gene << std::endl;

//     auto start_time = std::chrono::high_resolution_clock::now();

//     for (auto &p : kmerTableRef.count) {
//         long num_vec = p.first / n_slots;
//         long num_slot = p.first % n_slots;

//         // parameter checkpoint
//         // assert(n_genes == p.second.size());
//         for (size_t i = 0; i < p.second.size(); i++) {
//             // encode reverse polynomial version
//             if (num_slot == 0) {
//                 plain_vec[i][num_vec][num_slot] += p.second[i];
//             } else {
//                 plain_vec[i][num_vec][n_slots - num_slot] -= p.second[i];
//             }
//         }
//     }
//     auto table_time = std::chrono::high_resolution_clock::now();
//     auto table_duration = std::chrono::duration_cast<std::chrono::seconds>(table_time - start_time).count();

//     std::cout << "Ref Kmer Table generated" << std::endl;
//     std::cout << "encodeRefKmer table duration = " << table_duration << " s" << std::endl;
//     std::cout << "Memory usage = " << getMemoryUsage() / 1024 / 1024 << " KB" << std::endl;

//     pt_ref.resize(n_genes);
//     for (int g = 0; g < n_genes; g++) {
//         pt_ref[g].resize(n_vec_per_gene);
//         for (int i = 0; i < n_vec_per_gene; i++) {
//             Plaintext plain = cc->MakeCoefPackedPlaintext(plain_vec[g][i]);
//             pt_ref[g][i] = plain;
//             // update progress bar
//             print_progress_bar("EncodeRefKmer", g * n_vec_per_gene + i, n_genes * n_vec_per_gene, start_time);
//         }
//     }

//     auto end_time = std::chrono::high_resolution_clock::now();
//     auto duration = std::chrono::duration_cast<std::chrono::seconds>(end_time - start_time).count();
//     std::cout << "encodeRefKmer duration = " << duration << " s" << std::endl;
// }

void multCtxtByRef(Ciphertext_2d &ct_out, Ciphertext_1d &ct, Plaintext_2d &pt_ref, CryptoContext<DCRTPoly> &cc) {
    int g_genes = pt_ref.size();
    ct_out.resize(g_genes);
    for (int g = 0; g < g_genes; g++) {
        assert(ct.size() == pt_ref[g].size());
        // multiply ciphertexts in ct with plaintexts in pt
        // result is stored in ct_out
        // ct_out[i] = ct[i] * pt[i]

        ct_out[g].resize(ct.size());
        for (size_t i = 0; i < ct.size(); i++) {
            ct_out[g][i] = cc->EvalMult(ct[i], pt_ref[g][i]);
        }
    }
}

void multCtxtByKmerTableRef2(Ciphertext_2d &ct_out, Ciphertext_1d &ct, KmerTable kmerTableRef, long K, CryptoContext<DCRTPoly> &cc) {
    long n_slots = cc->GetCryptoParameters()->GetElementParams()->GetCyclotomicOrder() / 2;
    long n_genes = kmerTableRef.n_gene;
    long n_vec_per_gene = ceil(pow(4, K) / (double)n_slots);
    ct_out.resize(n_genes);
    
    std::cout << "n_slots = " << n_slots << std::endl;
    std::cout << "n_genes = " << n_genes << std::endl;
    std::cout << "n_vec_per_gene = " << n_vec_per_gene << std::endl;

    auto start_time = std::chrono::high_resolution_clock::now();
    auto table_time = std::chrono::high_resolution_clock::now();
    auto table_duration = std::chrono::duration_cast<std::chrono::seconds>(table_time - start_time).count();

    std::cout << "Ref Kmer Table generated" << std::endl;
    std::cout << "encodeRefKmer table duration = " << table_duration << " s" << std::endl;
    std::cout << "Memory usage = " << getMemoryUsage() / 1024 / 1024 << " KB" << std::endl;


    for (int g = 0; g < n_genes; g++) {
        // assert(ct.size() == pt_ref[g].size());
        // multiply ciphertexts in ct with plaintexts in pt
        // result is stored in ct_out
        // ct_out[i] = ct[i] * pt[i]
        vector<vector<int64_t>> plain_vec(n_vec_per_gene, vector<int64_t>(n_slots, 0));
        Plaintext plain;
        
        map<size_t, size_t> geneEntry = kmerTableRef.count[g];        
        for (auto &kmerEntry: geneEntry) {
            long num_vec = kmerEntry.first / n_slots;
            long num_slot = kmerEntry.first % n_slots;
            if (num_slot == 0) {
                plain_vec[num_vec][num_slot] += kmerEntry.second;
            } else {
                plain_vec[num_vec][n_slots - num_slot] -= kmerEntry.second;
            }
        }
        
        ct_out[g].resize(ct.size());
        for (size_t i = 0; i < ct.size(); i++) {
            plain = cc->MakeCoefPackedPlaintext(plain_vec[i]);
            Ciphertext<DCRTPoly> ctxt = cc->EvalMult(ct[i], plain);
            ct_out[g][i] = cc->EvalMult(ct[i], plain);

            // update progress bar
            print_progress_bar("multCtxtByKmerTableRef", g * ct.size() + i, n_genes * ct.size(), start_time);
        }
        plain_vec.clear();
    }
}


void multCtxtByKmerTableRefFromSerial(Ciphertext_1d &ct_out, Ciphertext_1d &ct, KmerTable kmerTableRef, CryptoContext<DCRTPoly> &cc, PQuantParams &param) {
    auto start_time = std::chrono::high_resolution_clock::now();
    long n_slots = cc->GetCryptoParameters()->GetElementParams()->GetCyclotomicOrder() / 2;
    long n_genes = kmerTableRef.n_gene;
    long n_vec_per_gene = ceil(pow(4, param.k) / (double)n_slots);
    if (param.path_kmer_matrix.size() > 0) {
        n_vec_per_gene = (kmerTableRef.n_kmer_total - 1) / n_slots + 1;
    }
    ct_out.resize(n_genes);
    
    std::cout << "n_slots = " << n_slots << std::endl;
    std::cout << "n_genes = " << n_genes << std::endl;
    std::cout << "n_vec_per_gene = " << n_vec_per_gene << std::endl;

    
    auto table_time = std::chrono::high_resolution_clock::now();
    auto table_duration = std::chrono::duration_cast<std::chrono::seconds>(table_time - start_time).count();

    std::cout << "Ref Kmer Table generated" << std::endl;
    std::cout << "encodeRefKmer table duration = " << table_duration << " s" << std::endl;
    std::cout << "Memory usage = " << getMemoryUsage() / 1024 / 1024 << " KB" << std::endl;

    vector<size_t> kmer_index;
    if (param.path_kmer_matrix.size() > 0) {
        for (auto &p : kmerTableRef.entropy) {
            kmer_index.push_back(p.first);
        }
    }

    auto start_time_encode = std::chrono::high_resolution_clock::now();
    vector<vector<vector<int64_t>>> plain_vec(n_genes, vector<vector<int64_t>>(n_vec_per_gene, vector<int64_t>(n_slots, 0)));
    for (int g = 0; g < n_genes; g++) {
        // assert(ct.size() == pt_ref[g].size());
        // multiply ciphertexts in ct with plaintexts in pt
        // result is stored in ct_out
        // ct_out[i] = ct[i] * pt[i]
        
        size_t num_vec, num_slot;
        map<size_t, size_t> geneEntry = kmerTableRef.count[g];        
        for (auto &kmerEntry: geneEntry) {
            if (param.path_kmer_matrix.size() == 0) {
                num_vec = kmerEntry.first / n_slots;
                num_slot = kmerEntry.first % n_slots;
            } else {
                size_t index = std::distance(kmer_index.begin(), std::find(kmer_index.begin(), kmer_index.end(), kmerEntry.first));
                // pass if p.first is not in kmer_index
                if (index == kmer_index.size()) {
                    cout << "index = " << index << " : not found" << endl;
                    continue;
                }
                num_vec = index / n_slots;
                num_slot = index % n_slots;
            }
            
            if (num_slot == 0) {
                plain_vec[g][num_vec][num_slot] += kmerEntry.second;
            } else {
                plain_vec[g][num_vec][n_slots - num_slot] -= kmerEntry.second;
            }
        }
        // update progress bar
        print_progress_bar("encode reference", g, n_genes, start_time_encode);
    }
        // string path = DATAFOLDER + "/gene" + to_string(g);
        // fs::create_directory(DATAFOLDER);
        // fs::create_directory(path);
    auto start_time_mult = std::chrono::high_resolution_clock::now();
    for (int g = 0; g < n_genes; g++) {   
        for (size_t i = 0; i < ct.size(); i++) {
            Plaintext plain = cc->MakeCoefPackedPlaintext(plain_vec[g][i]);
            if (i == 0) {
                ct_out[g] = cc->EvalMult(ct[0], plain);
            } else {
                Ciphertext<DCRTPoly> ctxt = cc->EvalMult(ct[i], plain);
                cc->EvalAddInPlace(ct_out[g], ctxt);
            }

            // update progress bar
            print_progress_bar("multCtxtByKmerTableRefSerial", g * ct.size() + i, n_genes * ct.size(), start_time_mult);
        }
    }
}

void decCtxtOut(Plaintext_1d &pt_out, Ciphertext_1d &ct_out, CryptoContext<DCRTPoly> &cc, KeyPair<DCRTPoly> &keyPair) {
    // auto start_time = std::chrono::high_resolution_clock::now();

    for (size_t i = 0; i < ct_out.size(); i++) {
        Plaintext plain_out;
        cc->Decrypt(keyPair.secretKey, ct_out[i], &plain_out);
        pt_out.push_back(plain_out);

        // update progress bar
        // print_progress_bar(i, ct_out.size(), start_time);
    }
}

void sumUpCtxt(Ciphertext<DCRTPoly> &ct, Ciphertext_1d &ct_vec, CryptoContext<DCRTPoly> &cc) {
    // auto start_time = std::chrono::high_resolution_clock::now();
    std::string CTXT_FOLDER = "../crypto/ctxt";
    
    ct = ct_vec[0];
    for (size_t i = 1; i < ct_vec.size(); i++) {
        ct = cc->EvalAdd(ct, ct_vec[i]);

        // update progress bar
        // print_progress_bar(i, ct_vec.size(), start_time);
    }
}


void sumUpCtxtFromSerial(Ciphertext_1d &ct_out, size_t n_gene, size_t n_ctxt, CryptoContext<DCRTPoly> &cc, string path_output) {
    // auto start_time = std::chrono::high_resolution_clock::now();
    std::string CTXT_FOLDER = "../crypto/ctxt";
    if (path_output != "") {
        CTXT_FOLDER = path_output;
    }
    
    for (size_t g = 0; g < n_gene; g++) {
        string path = CTXT_FOLDER + "/gene" + to_string(g);
        Ciphertext<DCRTPoly> ct;
        if (Serial::DeserializeFromFile(path + "/" + "ctxt0.txt", ct, SerType::BINARY) == false) {
            std::cerr << "Could not read the ciphertext" << std::endl;
            return;
        }
        for (size_t i = 1; i < n_ctxt; i++) {
            Ciphertext<DCRTPoly> ct_load;
            if (Serial::DeserializeFromFile(path + "/" + "ctxt" + to_string(i) + ".txt", ct_load, SerType::BINARY) == false) {
                std::cerr << "Could not read the ciphertext" << std::endl;
                return;
            }
            ct = cc->EvalAdd(ct, ct_load);
            // update progress bar
            // print_progress_bar(i, ct_vec.size(), start_time);
        }
        ct_out.push_back(ct);
    }
}