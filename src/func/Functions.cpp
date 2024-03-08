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
    std::cout << std::endl;
}

void encryptReadKmer(KmerTable &kmerTableRead, Ciphertext_1d &ct, CryptoContext<DCRTPoly> &cc, KeyPair<DCRTPoly> &keyPair, PQuantParams &param) {
    int K = param.k;
    // encrypt kmerTableRead.countRead into vector of ciphertexts (ct)
    // first, plain vector is created by vec[kmer(num)] = kmerTableRead.countRead[i]
    // then, plain vector is encrypted into ciphertexts

    // create plain vector
    size_t n_slots = cc->GetCryptoParameters()->GetElementParams()->GetCyclotomicOrder() / 2;
    std::cout << "n_slots = " << n_slots << std::endl;
    size_t n_plain_vecs = (pow(4, K) - 1) / n_slots + 1;
    std::cout << "n_plain_vecs = " << n_plain_vecs << std::endl;
    vector<vector<int64_t>> plain_vec;

    size_t i = 0;
    auto start_time = std::chrono::high_resolution_clock::now();
    for (auto &p : kmerTableRead.countRead) {
        size_t num_vec = p.first / n_slots;
        size_t num_slot = p.first % n_slots;
        // Resize plain_vec if necessary
        if (num_vec >= plain_vec.size()) {
            plain_vec.resize(num_vec + 1, vector<int64_t>(n_slots, 0));
        }

        plain_vec[num_vec][num_slot] = p.second;
        print_progress_bar("crate plain vec", i, kmerTableRead.countRead.size(), start_time);
        i += 1;
    }

    start_time = std::chrono::high_resolution_clock::now();
    for (size_t i = 0; i < n_plain_vecs; i++) {
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
    size_t n_slots = cc->GetCryptoParameters()->GetElementParams()->GetCyclotomicOrder() / 2;
    size_t n_plain_vecs = (kmerTableRef.n_kmer_total - 1) / n_slots + 1;
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

    for (size_t i = 0; i < n_plain_vecs; i++) {
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

void encodeRefKmer(KmerTable &kmerTableRef, Plaintext_2d &pt_ref, CryptoContext<DCRTPoly> &cc, KeyPair<DCRTPoly> &keyPair, PQuantParams &param) {
    auto start_time = std::chrono::high_resolution_clock::now();
    size_t n_slots = cc->GetCryptoParameters()->GetElementParams()->GetCyclotomicOrder() / 2;
    size_t n_genes = kmerTableRef.n_gene;
    size_t n_vec_per_gene = ceil(pow(4, param.k) / (double)n_slots);
    if (param.filename_kmerTable.size() > 0) {
        n_vec_per_gene = (kmerTableRef.n_kmer_total - 1) / n_slots + 1;
    }
    
    
    std::cout << "n_slots = " << n_slots << std::endl;
    std::cout << "n_genes = " << n_genes << std::endl;
    std::cout << "n_vec_per_gene = " << n_vec_per_gene << std::endl;

    
    auto table_time = std::chrono::high_resolution_clock::now();
    auto table_duration = std::chrono::duration_cast<std::chrono::seconds>(table_time - start_time).count();

    std::cout << "Ref Kmer Table generated" << std::endl;
    std::cout << "encodeRefKmer table duration = " << table_duration << " s" << std::endl;
    std::cout << "Memory usage = " << getMemoryUsage() / 1024 / 1024 << " KB" << std::endl;

    vector<size_t> kmer_index;
    if (param.filename_kmerTable.size() > 0) {
        for (auto &p : kmerTableRef.entropy) {
            kmer_index.push_back(p.first);
        }
    }

    auto start_time_encode = std::chrono::high_resolution_clock::now();
    pt_ref.resize(n_genes);
    for (size_t g = 0; g < n_genes; g++) {
        pt_ref[g].resize(n_vec_per_gene);
        vector<vector<int64_t>> plain_vec(n_vec_per_gene, vector<int64_t>(n_slots, 0));

        // assert(ct.size() == pt_ref[g].size());
        // multiply ciphertexts in ct with plaintexts in pt
        // result is stored in ct_out
        // ct_out[i] = ct[i] * pt[i]
        
        size_t num_vec, num_slot;
        map<size_t, size_t> geneEntry = kmerTableRef.count[g];        
        for (auto &kmerEntry: geneEntry) {
            if (param.filename_kmerTable.size() == 0) {
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
                plain_vec[num_vec][num_slot] += kmerEntry.second;
            } else {
                plain_vec[num_vec][n_slots - num_slot] -= kmerEntry.second;
            }
        }

        for (size_t i = 0; i < n_vec_per_gene; i++) {
            pt_ref[g][i] = cc->MakeCoefPackedPlaintext(plain_vec[i]);
        }
        // update progress bar
        if (param.progress_bar)
            print_progress_bar("encode reference", g, n_genes, start_time_encode);
        plain_vec.clear();
    }
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
        for (size_t i = 0; i < ct.size(); i++) {
            ct_out[g][i] = cc->EvalMult(ct[i], pt_ref[g][i]);
        }
    }
}

void multCtxtByEncodedRef(Ciphertext_1d &ct_out, Ciphertext_1d &ct, Plaintext_2d &pt_ref, CryptoContext<DCRTPoly> &cc, PQuantParams &param) {
    // long n_slots = cc->GetCryptoParameters()->GetElementParams()->GetCyclotomicOrder() / 2;
    long n_genes = pt_ref.size();
    // long n_vec_per_gene = pt_ref[0].size();
    ct_out.resize(n_genes);
    
    cout << "before multCtxt" << endl;
    printMemoryUsage();
    cout << endl;
    
    auto start_time_mult = std::chrono::high_resolution_clock::now();
    for (int g = 0; g < n_genes; g++) {   
        for (size_t i = 0; i < ct.size(); i++) {
            if (i == 0) {
                ct_out[g] = cc->EvalMult(ct[0], pt_ref[g][i]);
            } else {
                Ciphertext<DCRTPoly> ctxt = cc->EvalMult(ct[i], pt_ref[g][i]);
                cc->EvalAddInPlace(ct_out[g], ctxt);
            }

            // update progress bar
            if (param.progress_bar)
                print_progress_bar("multCtxtByEncodedRef", g * ct.size() + i, n_genes * ct.size(), start_time_mult);
        }
    }
}

void multCtxtByKmerTableRef(Ciphertext_1d &ct_out, Ciphertext_1d &ct, KmerTable kmerTableRef, CryptoContext<DCRTPoly> &cc, PQuantParams &param) {
    auto start_time = std::chrono::high_resolution_clock::now();
    size_t encode_total_time = 0;
    size_t mult_total_time = 0;
    long n_slots = cc->GetCryptoParameters()->GetElementParams()->GetCyclotomicOrder() / 2;
    long n_genes = kmerTableRef.n_gene;
    long n_vec_per_gene = ceil(pow(4, param.k) / (double)n_slots);
    if (param.filename_kmerTable.size() > 0) {
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
    if (param.filename_kmerTable.size() > 0) {
        for (auto &p : kmerTableRef.entropy) {
            kmer_index.push_back(p.first);
        }
    }

    for (int g = 0; g < n_genes; g++) {
        // assert(ct.size() == pt_ref[g].size());
        // multiply ciphertexts in ct with plaintexts in pt
        // result is stored in ct_out
        // ct_out[i] = ct[i] * pt[i]
        vector<vector<int64_t>> plain_vec = vector<vector<int64_t>>(n_vec_per_gene, vector<int64_t>(n_slots, 0));
        
        size_t num_vec, num_slot;
        map<size_t, size_t> geneEntry = kmerTableRef.count[g];        
        auto start_time_encode = std::chrono::high_resolution_clock::now();
        for (auto &kmerEntry: geneEntry) {
            if (param.filename_kmerTable.size() == 0) {
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
                plain_vec[num_vec][num_slot] += kmerEntry.second;
            } else {
                plain_vec[num_vec][n_slots - num_slot] -= kmerEntry.second;
            }
        }
        auto end_time_encode = std::chrono::high_resolution_clock::now();
        auto duration_encode = std::chrono::duration_cast<std::chrono::milliseconds>(end_time_encode - start_time_encode).count();
        encode_total_time += duration_encode;

        auto start_time_mult = std::chrono::high_resolution_clock::now();
        for (size_t i = 0; i < ct.size(); i++) {
            Plaintext plain = cc->MakeCoefPackedPlaintext(plain_vec[i]);
            if (i == 0) {
                ct_out[g] = cc->EvalMult(ct[0], plain);
            } else {
                Ciphertext<DCRTPoly> ctxt = cc->EvalMult(ct[i], plain);
                cc->EvalAddInPlace(ct_out[g], ctxt);
            }

            // update progress bar
            print_progress_bar("multCtxtByKmerTableRefSerial", g * ct.size() + i, n_genes * ct.size(), start_time_mult);
        }
        plain_vec.clear();
        geneEntry.clear();
        auto end_time_mult = std::chrono::high_resolution_clock::now();
        auto duration_mult = std::chrono::duration_cast<std::chrono::milliseconds>(end_time_mult - start_time_mult).count();
        mult_total_time += duration_mult;
    }
    
    cout << "end multCtxt" << endl;
    cout << "encode_total_time = " << encode_total_time << " ms" << endl;
    cout << "mult_total_time = " << mult_total_time << " ms" << endl;
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