#include "MainAlgorithmSet.h"

void MainAlgorithmSet::generateKmerTableFromReference(PQuantParams &param) {
    // read fasta file (reference)
    vector<Sequence> seq_vec;
    auto start_time = std::chrono::high_resolution_clock::now();
    cout << "=== read fasta file ===" << endl;
    readFastaFile(param.filename_ref, seq_vec);
    auto end_time = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time).count();
    std::cout << "readFastaFile duration = " << duration << " ms" << std::endl;
    std::cout << std::endl;

    // generate kmerTable

    cout << "=== compute kmerTable ===" << endl;
    auto start_time_kmerTable = std::chrono::high_resolution_clock::now();
    KmerTable kmerTable(seq_vec, param, true);
    auto end_time_kmerTable = std::chrono::high_resolution_clock::now();
    auto duration_kmerTable = std::chrono::duration_cast<std::chrono::milliseconds>(end_time_kmerTable - start_time_kmerTable).count();
    std::cout << "computeKmerTable duration = " << duration_kmerTable << " ms" << std::endl;
    std::cout << std::endl;
    if (param.verbose) {
        kmerTable.print();
    }
    // save kmerTable & kmerList to file
    cout << "=== save kmerTable & kmerList to file ===" << endl;
    std::filesystem::path dir_path = param.foldername_kmer;
    if (!std::filesystem::exists(dir_path)) {
        if(std::filesystem::create_directories(dir_path)) {
            std::cout << "Directory created: " << dir_path << std::endl;
        } else {
            std::cerr << "Failed to create directory: " << dir_path << std::endl;
        }
    }
    auto start_time_save = std::chrono::high_resolution_clock::now();
    if (param.json_format) {
        kmerTable.save(param.filename_kmerTable);
        kmerTable.saveKmerList(param.filename_kmerList);
    } else {
        kmerTable.saveBinary(param.filename_kmerTable);
        kmerTable.saveKmerListBinary(param.filename_kmerList);
    }
    auto end_time_save = std::chrono::high_resolution_clock::now();
    auto duration_save = std::chrono::duration_cast<std::chrono::milliseconds>(end_time_save - start_time_save).count();
    std::cout << "save duration = " << duration_save << " ms" << std::endl;
    std::cout << std::endl;

    cout << " === Duration summaries ===" << endl;
    cout << "readFastaFile duration = " << duration << " ms" << endl;
    cout << "computeKmerTable duration = " << duration_kmerTable << " ms" << endl;
    cout << "save duration = " << duration_save << " ms" << endl;
    cout << "total duration = " << duration + duration_kmerTable + duration_save << " ms" << endl;
    printFileSize(param.filename_kmerTable);
    printFileSize(param.filename_kmerList);
    printMemoryUsage();
}

void MainAlgorithmSet::keyGenBFVandSerialize(PQuantParams &param) {
    CCParams<CryptoContextBFVRNS> parameters;
    parameters.SetPlaintextModulus(536903681);
    parameters.SetMultiplicativeDepth(1);
    parameters.SetMaxRelinSkDeg(2);

    std::cout << "=== GenCryptoContext ===" << std::endl;
    auto start_time = std::chrono::high_resolution_clock::now();
    CryptoContext<DCRTPoly> cryptoContext = GenCryptoContext(parameters);
    // enable features that you wish to use
    cryptoContext->Enable(PKE);
    cryptoContext->Enable(KEYSWITCH);
    cryptoContext->Enable(LEVELEDSHE);
    cryptoContext->Enable(ADVANCEDSHE);
    auto end_time = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time).count();
    std::cout << "GenCryptoContext duration = " << duration << " ms" << std::endl;
    std::cout << std::endl;

    std::cout << "\np = " << cryptoContext->GetCryptoParameters()->GetPlaintextModulus() << std::endl;
    auto n = cryptoContext->GetCryptoParameters()->GetElementParams()->GetCyclotomicOrder() / 2;
    std::cout << "n = " << n << endl;
    std::cout << "log2 q = "
              << log2(cryptoContext->GetCryptoParameters()->GetElementParams()->GetModulus().ConvertToDouble())
              << std::endl;
    std::cout << std::endl;

    // Initialize Public Key Containers
    std::cout << "=== KeyGen ===" << std::endl;
    auto start_time_keygen = std::chrono::high_resolution_clock::now();
    KeyPair<DCRTPoly> keyPair = cryptoContext->KeyGen();
    cryptoContext->EvalMultKeysGen(keyPair.secretKey);
    auto end_time_keygen = std::chrono::high_resolution_clock::now();
    auto duration_keygen = std::chrono::duration_cast<std::chrono::milliseconds>(end_time_keygen - start_time_keygen).count();
    std::cout << "KeyGen duration = " << duration_keygen << " ms" << std::endl;
    std::cout << std::endl;

    // create fildername_BFV
    std::filesystem::path dir_path = param.foldername_BFV;
    if (!std::filesystem::exists(dir_path)) {
        if(std::filesystem::create_directories(dir_path)) {
            std::cout << "Directory created: " << dir_path << std::endl;
        } else {
            std::cerr << "Failed to create directory: " << dir_path << std::endl;
        }
    }

    std::string filename_context = param.foldername_BFV + "/context.txt";
    std::string filename_public = param.foldername_BFV + "/key-public.txt";
    std::string filename_private = param.foldername_BFV + "/key-private.txt";
    std::string filename_eval_mult = param.foldername_BFV + "/key-eval-mult.txt";

    std::cout << "=== Serialize ===" << std::endl;
    auto start_time_serialize = std::chrono::high_resolution_clock::now();
    // Serialize cryptocontext
    if (!Serial::SerializeToFile(filename_context, cryptoContext, SerType::BINARY)) {
        std::cerr << "Error writing serialization of the crypto context to "
                     "cryptocontext.txt"
                  << std::endl;
        return;
    }
    std::cout << "The cryptocontext has been serialized in " << filename_context << std::endl;
    
    // Serialize the public key
    if (!Serial::SerializeToFile(filename_public, keyPair.publicKey, SerType::BINARY)) {
        std::cerr << "Error writing serialization of public key to key-public.txt" << std::endl;
        return;
    }
    std::cout << "The public key has been serialized in " << filename_private << std::endl;
    
    // Serialize the secret key
    if (!Serial::SerializeToFile(filename_private, keyPair.secretKey, SerType::BINARY)) {
        std::cerr << "Error writing serialization of private key to key-private.txt" << std::endl;
        return;
    }
    std::cout << "The secret key has been serialized in " << filename_private << std::endl;

    // multiplication
    std::ofstream emkeyfile(filename_eval_mult, std::ios::out | std::ios::binary);
    if (emkeyfile.is_open()) {
        if (cryptoContext->SerializeEvalMultKey(emkeyfile, SerType::BINARY) == false) {
            std::cerr << "Error writing serialization of the eval mult keys to "
                         "key-eval-mult.txt"
                      << std::endl;
            return;
        }
        std::cout << "The eval mult keys have been serialized in " << filename_eval_mult << std::endl;
        emkeyfile.close();
    }
    auto end_time_serialize = std::chrono::high_resolution_clock::now();
    auto duration_serialize = std::chrono::duration_cast<std::chrono::milliseconds>(end_time_serialize - start_time_serialize).count();
    std::cout << "Serialize duration = " << duration_serialize << " ms" << std::endl;

    std::cout << " === Duration summaries ===" << endl;
    std::cout << "GenCryptoContext duration = " << duration << " ms" << std::endl;
    std::cout << "KeyGen duration = " << duration_keygen << " ms" << std::endl;
    std::cout << "Serialize duration = " << duration_serialize << " ms" << std::endl;
    std::cout << "total duration = " << duration + duration_keygen + duration_serialize << " ms" << std::endl;
    printFileSize(filename_context);
    printFileSize(filename_public);
    printFileSize(filename_private);
    printFileSize(filename_eval_mult);
    printMemoryUsage();
}

void MainAlgorithmSet::encodeAndEncrypt(PQuantParams &param) {
    std::string filename_context = param.foldername_BFV + "/context.txt";
    std::string filename_private = param.foldername_BFV + "/key-private.txt";

    std::cout << "=== Deserialize ===" << std::endl;
    auto start_time = std::chrono::high_resolution_clock::now();
    // Deserialize the crypto context
    CryptoContext<DCRTPoly> cc;
    if (!Serial::DeserializeFromFile(filename_context, cc, SerType::BINARY)) {
        std::cerr << "I cannot read serialization from " << filename_context << std::endl;
        return;
    }
    std::cout << "The cryptocontext has been deserialized from " << filename_context << std::endl;

    PrivateKey<DCRTPoly> sk;
    if (!Serial::DeserializeFromFile(filename_private, sk, SerType::BINARY)) {
        std::cerr << "Could not read secret key" << std::endl;
        return;
    }
    std::cout << "The secret key has been deserialized from " << filename_private << std::endl;
    auto end_time = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time).count();
    std::cout << "Deserialize duration = " << duration << " ms" << std::endl;
    std::cout << std::endl;

    // read kmerList
    std::cout << "=== load kmerList ===" << std::endl;
    auto start_time_kmerList = std::chrono::high_resolution_clock::now();
    vector<size_t> kmer_list;
    if (param.json_format) {
        loadKmerList(param.filename_kmerList, kmer_list);
    } else {
        loadKmerListBinary(param.filename_kmerList, kmer_list);
    }
    auto end_time_kmerList = std::chrono::high_resolution_clock::now();
    auto duration_kmerList = std::chrono::duration_cast<std::chrono::milliseconds>(end_time_kmerList - start_time_kmerList).count();
    std::cout << "load kmerList duration = " << duration_kmerList << " ms" << std::endl;
    std::cout << std::endl;

    cout << "=== read reads and compute kmerTableRead ===" << endl;
    auto start_time_read = std::chrono::high_resolution_clock::now();
    vector<Sequence> reads_seq;
    // check filename_read if it is .fa or .fq
    if (param.filename_read.find(".fa") != std::string::npos && param.filename_read.find(".fastq") == std::string::npos) {
        cout << "read file is .fa format" << endl;
        readFastaFile(param.filename_read, reads_seq);
    } else if (param.filename_read.find(".fq") != std::string::npos || param.filename_read.find(".fastq") != std::string::npos) {
        cout << "read file is .fq or .fastq format" << endl;
        readFastQFile(param.filename_read, reads_seq);
    } else {
        std::cout << "Invalid file format" << std::endl;
        std::cout << "Assume it is .fastq format" << std::endl;
        readFastQFile(param.filename_read, reads_seq);
    }
    // readFastQFile(param.filename_read, reads_seq);
    // KmerTable kmerTableRead(reads_seq, param, false);
    KmerTable kmerTableRead(reads_seq, param, kmer_list);
    reads_seq.clear();
    auto end_time_read = std::chrono::high_resolution_clock::now();
    auto duration_read = std::chrono::duration_cast<std::chrono::milliseconds>(end_time_read - start_time_read).count();
    std::cout << "readFastaFile duration = " << duration_read << " ms" << std::endl;
    if (param.verbose) {
        kmerTableRead.print();
    }
    std::cout << std::endl;

    if (param.verbose) {
        for (size_t i = 0; i < kmer_list.size(); i++) {
            std::cout << kmer_list[i] << " ";
        }
        cout << "kmer_list.size() = " << kmer_list.size() << endl;
    }

    size_t n_kmer_total = kmer_list.size();
    size_t n_slots = cc->GetCryptoParameters()->GetElementParams()->GetCyclotomicOrder() / 2;
    size_t n_plain_vecs = (n_kmer_total - 1) / n_slots + 1;
    vector<vector<int64_t>> plain_vec(n_plain_vecs, vector<int64_t>(n_slots, 0));
    
    std::cout << "n_slots = " << n_slots << std::endl;
    std::cout << "n_plain_vecs = " << n_plain_vecs << std::endl;
    cout << "n_kmer_total = " << n_kmer_total << endl;

    // encode kmerTableRead.countRead based on kmer_index
    std::cout << "=== encode kmerTableRead.countRead ===" << std::endl;
    auto start_time_encode = std::chrono::high_resolution_clock::now();
    int count = 0;
    size_t index;
    for (auto &p : kmerTableRead.countRead) {
        bool search = binary_search_util(kmer_list, p.first, index);
        // auto it = std::lower_bound(kmer_list.begin(), kmer_list.end(), p.first);
        if (search) {
        // if (it != kmer_list.end() && *it == p.first) {
            // size_t index = std::distance(kmer_list.begin(), it);
            size_t num_vec = index / n_slots;
            size_t num_slot = index % n_slots;
            plain_vec[num_vec][num_slot] = p.second;
        }
        if (param.progress_bar) 
            print_progress_bar("encodeRead", count, kmerTableRead.countRead.size(), start_time_encode);
        count += 1;
    }
    auto end_time_encode = std::chrono::high_resolution_clock::now();
    auto duration_encode = std::chrono::duration_cast<std::chrono::milliseconds>(end_time_encode - start_time_encode).count();
    std::cout << "encodeRead duration = " << duration_encode << " ms" << std::endl;
    std::cout << std::endl;

    // encrypt plain_vec
    std::cout << "=== encrypt plain_vec ===" << std::endl;
    
    std::string filename_ctxtRead = param.foldername_ctxtread;
    // create folder if doesn't exist
    if (!std::filesystem::exists(filename_ctxtRead)) {
        std::filesystem::create_directory(filename_ctxtRead);
    }

    size_t duration_encrypt = 0;
    size_t duration_serialize = 0;
    size_t duration_encrypt_then_serialize = 0;
    
    cout << "before enc, check memory" << endl;
    printMemoryUsage();
    if (param.operate_then_serialize) {
        auto start_time_encrypt = std::chrono::high_resolution_clock::now();
        Ciphertext_1d ct;
        for (size_t i = 0; i < n_plain_vecs; i++) {
            Plaintext plain = cc->MakeCoefPackedPlaintext(plain_vec[i]);
            Ciphertext<DCRTPoly> ciphertext = cc->Encrypt(sk, plain);
            ct.push_back(ciphertext);

            // update progress bar
            if (param.progress_bar)
                print_progress_bar("EncryptReadKmer", i, n_plain_vecs, start_time_encrypt);
        }
        auto end_time_encrypt = std::chrono::high_resolution_clock::now();
        duration_encrypt = std::chrono::duration_cast<std::chrono::milliseconds>(end_time_encrypt - start_time_encrypt).count();
        std::cout << "EncryptReadKmer duration = " << duration_encrypt << " ms" << std::endl;
        std::cout << std::endl;
        
        // serialize ct
        std::cout << "=== serialize ct ===" << std::endl;
        auto start_time_serialize = std::chrono::high_resolution_clock::now();
        // create folder if doesn't exist
        if (!std::filesystem::exists(filename_ctxtRead)) {
            std::filesystem::create_directory(filename_ctxtRead);
        }
        for (size_t i = 0; i < ct.size(); i++) {
            std::string filename_ctxtRead_i = filename_ctxtRead + "/ct_" + std::to_string(i) + ".txt";
            if (!Serial::SerializeToFile(filename_ctxtRead_i, ct[i], SerType::BINARY)) {
                std::cerr << "Error writing serialization of the ctxtRead to " << filename_ctxtRead_i << std::endl;
                return;
            }
            if (param.progress_bar) {
                print_progress_bar("serializeCtxtRead", i, ct.size(), start_time_serialize);
            }
        }
        auto end_time_serialize = std::chrono::high_resolution_clock::now();
        duration_serialize = std::chrono::duration_cast<std::chrono::milliseconds>(end_time_serialize - start_time_serialize).count();
        std::cout << "serializeCtxtRead duration = " << duration_serialize << " ms" << std::endl;
        std::cout << std::endl;
    } else {
        auto start_time_encrypt_and_serialize = std::chrono::high_resolution_clock::now();
        Plaintext plain;
        Ciphertext<DCRTPoly> ciphertext;
        std::string filename_ctxtRead_i;
        for (size_t i = 0; i < n_plain_vecs; i++) {
            // encode & encrypt
            plain = cc->MakeCoefPackedPlaintext(plain_vec[i]);
            ciphertext = cc->Encrypt(sk, plain);

            // serialize ctxt
            filename_ctxtRead_i = filename_ctxtRead + "/ct_" + std::to_string(i) + ".txt";
            if (!Serial::SerializeToFile(filename_ctxtRead_i, ciphertext, SerType::BINARY)) {
                std::cerr << "Error writing serialization of the ctxtRead to " << filename_ctxtRead_i << std::endl;
                return;
            }

            // update progress bar
            if (param.progress_bar)
                print_progress_bar("EncryptAndSerializeReadKmer", i, n_plain_vecs, start_time_encrypt_and_serialize);
        }
        auto end_time_encrypt_and_serialize = std::chrono::high_resolution_clock::now();
        duration_encrypt_then_serialize = std::chrono::duration_cast<std::chrono::milliseconds>(end_time_encrypt_and_serialize - start_time_encrypt_and_serialize).count();
        std::cout << "EncryptAndSerializeReadKmer duration = " << duration_encrypt_then_serialize << " ms" << std::endl;
        std::cout << std::endl;
    }
    
    auto total_duration = duration_read + duration_kmerList + duration_encode;
    std::cout << " === Duration summaries ===" << endl;
    std::cout << "readFastaFile duration = " << duration_read << " ms" << endl;
    std::cout << "load kmerList duration = " << duration_kmerList << " ms" << endl;
    std::cout << "encodeRead duration = " << duration_encode << " ms" << endl;
    if (param.operate_then_serialize) {
        std::cout << "EncryptReadKmer duration = " << duration_encrypt << " ms" << endl;
        std::cout << "serializeCtxtRead duration = " << duration_serialize << " ms" << std::endl;
        total_duration = total_duration + duration_encrypt + duration_serialize;
    } else {
        std::cout << "EncryptAndSerializeReadKmer duration = " << duration_encrypt_then_serialize << " ms" << endl;
        total_duration = total_duration + duration_encrypt_then_serialize;
    }
    
    std::cout << "total duration = " << total_duration << " ms" << std::endl;
    std::cout << std::endl;
    std::cout << "single ctxt size" << endl;
    printFileSize(filename_ctxtRead + "/ct_0.txt");
    printFolderSize(filename_ctxtRead);
    printMemoryUsage();
}
    
void MainAlgorithmSet::computeInnerProductBatch(PQuantParams &param) {
    std::string filename_context = param.foldername_BFV + "/context.txt";
    std::string filename_eval_mult = param.foldername_BFV + "/key-eval-mult.txt";

    // Deserialize the crypto context
    std::cout << "=== Deserialize ===" << std::endl;
    auto start_time = std::chrono::high_resolution_clock::now();
    CryptoContext<DCRTPoly> cc;
    if (!Serial::DeserializeFromFile(filename_context, cc, SerType::BINARY)) {
        std::cerr << "I cannot read serialization from " << filename_context << std::endl;
        return;
    }
    std::cout << "The cryptocontext has been deserialized from " << filename_context << std::endl;

    // Deserialize the eval mult keys
    std::ifstream emkeyfile(filename_eval_mult, std::ios::in | std::ios::binary);
    if (!emkeyfile.is_open()) {
        std::cerr << "Could not read eval mult keys" << std::endl;
        return;
    }
    cc->DeserializeEvalMultKey(emkeyfile, SerType::BINARY);
    std::cout << "The eval mult keys have been deserialized from " << filename_eval_mult << std::endl;
    auto end_time = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time).count();
    std::cout << "Deserialize duration = " << duration << " ms" << std::endl;
    std::cout << std::endl;

    // read kmerTable
    std::cout << "=== read kmerTable ===" << std::endl;
    auto start_time_kmerTable = std::chrono::high_resolution_clock::now();
    KmerTable kmerTableRef;
    vector<size_t> kmer_list;
    if (param.json_format) {
        kmerTableRef.load(param.filename_kmerTable);
    } else {
        kmerTableRef.loadBinary(param.filename_kmerTable);
    }
    if (param.verbose) {
        kmerTableRef.print();
    }
    auto end_time_kmerTable = std::chrono::high_resolution_clock::now();
    auto duration_kmerTable = std::chrono::duration_cast<std::chrono::milliseconds>(end_time_kmerTable - start_time_kmerTable).count();
    std::cout << "load kmerTable duration = " << duration_kmerTable << " ms" << std::endl;
    std::cout << std::endl;

    size_t start = 0;
    size_t end = kmerTableRef.n_gene;
    if (param.batch_num_total >= 0) {
        long n_genes_per_batch = (kmerTableRef.n_gene - 1) / param.batch_num_total + 1;
        if (param.batch_num >= 0) {
            start = param.batch_num * n_genes_per_batch;
            end = (param.batch_num + 1) * n_genes_per_batch;
            end = end < kmerTableRef.n_gene ? end : kmerTableRef.n_gene;
        } else {
            start = 0;
            end = kmerTableRef.n_gene;
        }
    } else if (param.gene_start >= 0 && param.gene_end >= 0) {
        start = param.gene_start;
        end = param.gene_end;
    } 
    cout << "Run genes from " << start << " to " << end - 1 << endl;

    // read kmerList
    std::cout << "=== load kmerList ===" << std::endl;
    auto start_time_kmerList = std::chrono::high_resolution_clock::now();
    if (param.json_format) {
        loadKmerList(param.filename_kmerList, kmer_list);
    } else {
        loadKmerListBinary(param.filename_kmerList, kmer_list);
    }
    auto end_time_kmerList = std::chrono::high_resolution_clock::now();
    auto duration_kmerList = std::chrono::duration_cast<std::chrono::milliseconds>(end_time_kmerList - start_time_kmerList).count();
    std::cout << "load kmerList duration = " << duration_kmerList << " ms" << std::endl;
    std::cout << std::endl;

    // set parameters
    size_t n_slots = cc->GetCryptoParameters()->GetElementParams()->GetCyclotomicOrder() / 2;
    size_t n_genes = end - start;
    size_t n_vec_per_gene = (kmerTableRef.n_kmer_total - 1) / n_slots + 1;

    std::cout << "check parameters" << std::endl;
    std::cout << "n_slots = " << n_slots << std::endl;
    std::cout << "n_genes = " << n_genes << std::endl;
    std::cout << "n_vec_per_gene = " << n_vec_per_gene << std::endl;
    std::cout << "n_kmer_total = " << kmerTableRef.n_kmer_total << std::endl;
    std::cout << std::endl;
    
    Ciphertext_1d ct(n_vec_per_gene);

    // debug, don't know why this is needed
    vector<int64_t> dummy_vec(n_slots, 0);
    Plaintext dummy = cc->MakeCoefPackedPlaintext(dummy_vec);

    // read filename_ctxtRead
    std::string filename_ctxtRead = param.foldername_ctxtread;
    std::cout << "=== read ctxtRead ===" << std::endl;
    auto start_time_ctxtRead = std::chrono::high_resolution_clock::now();
    for (size_t r = 0; r < n_vec_per_gene; r++) {
        if (!Serial::DeserializeFromFile(filename_ctxtRead + "/ct_" + std::to_string(r) + ".txt", ct[r], SerType::BINARY)) {
            std::cerr << "Error reading serialization of the ctxtRead" << std::endl;
            return;
        }
        if (param.progress_bar)
            print_progress_bar("read ctxtRead", r, n_vec_per_gene, start_time_ctxtRead);
    }
    auto end_time_ctxtRead = std::chrono::high_resolution_clock::now();
    auto duration_ctxtRead = std::chrono::duration_cast<std::chrono::milliseconds>(end_time_ctxtRead - start_time_ctxtRead).count();
    std::cout << "read ctxtRead duration = " << duration_ctxtRead << " ms" << std::endl;
    std::cout << std::endl;
    
    // encode kmerTable
    std::cout << "=== encode kmerTable and mult ===" << std::endl;
    auto start_time_encode_and_mult = std::chrono::high_resolution_clock::now();
    Ciphertext_1d ct_out(n_genes);
    for (size_t g = start; g < end; g++) {
        vector<vector<int64_t>> plain_vec(n_vec_per_gene, vector<int64_t>(n_slots, 0));
        Plaintext_1d pt_ref(n_vec_per_gene);
        size_t num_vec, num_slot;
        for (auto &kmerEntry: kmerTableRef.count[g]) {
            size_t index;
            if (!binary_search_util(kmer_list, kmerEntry.first, index)) {
                continue;
            }
            num_vec = index / n_slots;
            num_slot = index % n_slots;
            
            if (num_slot == 0) {
                plain_vec[num_vec][num_slot] += kmerEntry.second;
            } else {
                plain_vec[num_vec][n_slots - num_slot] -= kmerEntry.second;
            }
        }
        for (size_t i = 0; i < ct.size(); i++) {
            Plaintext pt_ref = cc->MakeCoefPackedPlaintext(plain_vec[i]);
            if (i == 0) {
                ct_out[g - start] = cc->EvalMult(ct[0], pt_ref);
            } else {
                Ciphertext<DCRTPoly> ctxt = cc->EvalMult(ct[i], pt_ref);
                cc->EvalAddInPlace(ct_out[g - start], ctxt);
            }
            // update progress bar
            if (param.progress_bar)
                print_progress_bar("multCtxtByEncodedRef", g * ct.size() + i, n_genes * ct.size(), start_time_encode_and_mult);
        }
    }
    auto end_time_encode_and_mult = std::chrono::high_resolution_clock::now();
    auto duration_encode_end_mult = std::chrono::duration_cast<std::chrono::milliseconds>(end_time_encode_and_mult - start_time_encode_and_mult).count();
    std::cout << "multCtxtByEncodedRef duration = " << duration_encode_end_mult << " ms" << std::endl;
    std::cout << std::endl;

    // serialize ct_out
    std::cout << "=== serialize ct_out ===" << std::endl;
    auto start_time_serialize = std::chrono::high_resolution_clock::now();
    std::string filename_ctxtOut = param.foldername_ctxtout;
    // create folder if doesn't exist
    if (!std::filesystem::exists(filename_ctxtOut)) {
        std::filesystem::create_directory(filename_ctxtOut);
    }
    for (size_t i = start; i < end; i++) {
        std::string filename_ctxtOut_i = filename_ctxtOut + "/ct_" + std::to_string(i) + ".txt";
        if (!Serial::SerializeToFile(filename_ctxtOut_i, ct_out[i - start], SerType::BINARY)) {
            std::cerr << "Error writing serialization of the ctxt_out to " << filename_ctxtOut_i << std::endl;
            return;
        }
        if (param.progress_bar)
            print_progress_bar("serializeCtxtOut", i, n_genes, start_time_serialize);
    }
    auto end_time_serialize = std::chrono::high_resolution_clock::now();
    auto duration_serialize = std::chrono::duration_cast<std::chrono::milliseconds>(end_time_serialize - start_time_serialize).count();
    std::cout << "serializeCtxtOut duration = " << duration_serialize << " ms" << std::endl;
    std::cout << std::endl;

    std::cout << " === Duration summaries ===" << endl;
    std::cout << "load kmerTable duration = " << duration_kmerTable << " ms" << endl;
    std::cout << "load kmerList duration = " << duration_kmerList << " ms" << endl;
    std::cout << "read ctxtRead duration = " << duration_ctxtRead << " ms" << endl;
    std::cout << "encode_and_mult duration = " << duration_encode_end_mult << " ms" << std::endl;
    std::cout << "serializeCtxtOut duration = " << duration_serialize << " ms" << std::endl;
    std::cout << "total duration = " << duration_kmerTable + duration_kmerList + duration_encode_end_mult + duration_ctxtRead + duration_serialize << " ms" << std::endl;

    std::cout << "single ctxt size" << endl;
    printFileSize(filename_ctxtOut + "/ct_" + to_string(start) + ".txt");
    printFolderSize(filename_ctxtOut);
    printMemoryUsage();
}

void MainAlgorithmSet::decryptAndReturnGeneVector(PQuantParams &param) {
    std::string filename_context = param.foldername_BFV + "/context.txt";
    std::string filename_private = param.foldername_BFV + "/key-private.txt";

    // Deserialize the crypto context
    std::cout << "=== Deserialize ===" << std::endl;
    auto start_time = std::chrono::high_resolution_clock::now();
    CryptoContext<DCRTPoly> cc;
    if (!Serial::DeserializeFromFile(filename_context, cc, SerType::BINARY)) {
        std::cerr << "I cannot read serialization from " << filename_context << std::endl;
        return;
    }
    std::cout << "The cryptocontext has been deserialized from " << filename_context << std::endl;

    PrivateKey<DCRTPoly> sk;
    if (Serial::DeserializeFromFile(filename_private, sk, SerType::BINARY) == false) {
        std::cerr << "Could not read secret key" << std::endl;
        return;
    }
    std::cout << "The secret key has been deserialized from " << filename_private << std::endl;
    auto end_time = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time).count();
    std::cout << "Deserialize duration = " << duration << " ms" << std::endl;
    std::cout << std::endl;

    // read kmerTable
    // std::cout << "=== read kmerTable ===" << std::endl;
    // auto start_time_kmerTable = std::chrono::high_resolution_clock::now();
    // KmerTable kmerTableRef;
    //     if (param.json_format) {
    //     kmerTableRef.load(param.filename_kmerTable);
    // } else {
    //     kmerTableRef.loadBinary(param.filename_kmerTable);
    // }
    // if (param.verbose) {
    //     kmerTableRef.print();
    // }
    // auto end_time_kmerTable = std::chrono::high_resolution_clock::now();
    // auto duration_kmerTable = std::chrono::duration_cast<std::chrono::milliseconds>(end_time_kmerTable - start_time_kmerTable).count();
    // std::cout << "load kmerTable duration = " << duration_kmerTable << " ms" << std::endl;
    // std::cout << std::endl;

    // size_t n_genes = kmerTableRef.n_gene;
    // param.gene_start = 0;
    // param.gene_end = n_genes;

    size_t n_genes = 0;

    try {
        for (const auto& entry : std::filesystem::directory_iterator(param.foldername_ctxtout)) {
            if (entry.is_regular_file()) {
                ++n_genes;
            }
        }
    } catch (std::filesystem::filesystem_error& e) {
        std::cerr << "Filesystem error: " << e.what() << std::endl;
    } catch (std::exception& e) {
        std::cerr << "General error: " << e.what() << std::endl;
    }
    
    size_t duration_out = 0;
    size_t duration_dec = 0;
    size_t duration_out_then_dec = 0;
    vector<long> final(n_genes, 0);
    if (param.operate_then_serialize) {
        // read ctxt_out
        std::cout << "=== read ctxt_out ===" << std::endl;
        auto start_time_ctxtOut = std::chrono::high_resolution_clock::now();
        Ciphertext_1d ct_out;
        //std::string filename_ctxtOut = param.foldername_BFV + "/ctxt_out";
        std::string filename_ctxtOut = param.foldername_ctxtout;

        for (size_t i = 0; i < n_genes; i++) {
            std::string filename_ctxtOut_i = filename_ctxtOut + "/ct_" + std::to_string(i) + ".txt";
            Ciphertext<DCRTPoly> ctxt;
            if (!Serial::DeserializeFromFile(filename_ctxtOut_i, ctxt, SerType::BINARY)) {
                std::cerr << "Error reading serialization of the ctxt_out from " << filename_ctxtOut_i << std::endl;
                return;
            }
            ct_out.push_back(ctxt);
            if (param.progress_bar)
                print_progress_bar("read ctxtOut", i, n_genes, start_time_ctxtOut);
        }
        cout << "ctxt_out loaded: num = " << ct_out.size() << endl;
        auto end_time_ctxtOut = std::chrono::high_resolution_clock::now();
        duration_out = std::chrono::duration_cast<std::chrono::milliseconds>(end_time_ctxtOut - start_time_ctxtOut).count();
        std::cout << "read ctxtOut duration = " << duration_out << " ms" << std::endl;
        std::cout << std::endl;

        cout << " === run decrypt_output === " << endl;
        auto start_time_dec = std::chrono::high_resolution_clock::now();
        Plaintext_1d pt_sum(ct_out.size());
        for (size_t i = 0; i < ct_out.size(); i++) {
            cc->Decrypt(sk, ct_out[i], &pt_sum[i]);
            if (param.progress_bar)
                print_progress_bar("decryptOutput", i, ct_out.size(), start_time_dec);
        }
        auto end_time_dec = std::chrono::high_resolution_clock::now();
        duration_dec = std::chrono::duration_cast<std::chrono::milliseconds>(end_time_dec - start_time_dec).count();
        std::cout << "decryptOutput duration = " << duration_dec << " ms" << std::endl;
        std::cout << std::endl;

        for (size_t i = 0; i < pt_sum.size(); i++) {
            auto plain = pt_sum[i]->GetCoefPackedValue();
            final[i] = plain[0];
        }
    } else {
        // read ctxt_out
        std::cout << "=== read ctxt_out ===" << std::endl;
        auto start_time_ctxt = std::chrono::high_resolution_clock::now();
        //std::string filename_ctxtOut = param.foldername_BFV + "/ctxt_out";
        std::string filename_ctxtOut = param.foldername_ctxtout;


        Plaintext plain;
        Ciphertext<DCRTPoly> ctxt;
        for (size_t i = 0; i < n_genes; i++) {
            std::string filename_ctxtOut_i = filename_ctxtOut + "/ct_" + std::to_string(i) + ".txt";
            if (!Serial::DeserializeFromFile(filename_ctxtOut_i, ctxt, SerType::BINARY)) {
                std::cerr << "Error reading serialization of the ctxt_out from " << filename_ctxtOut_i << std::endl;
                return;
            }
            if (param.progress_bar)
                print_progress_bar("read ctxtOut", i, n_genes, start_time_ctxt);

            cc->Decrypt(sk, ctxt, &plain);
            final[i] = plain->GetCoefPackedValue()[0];
            if (param.progress_bar)
                print_progress_bar("decryptOutput", i, n_genes, start_time_ctxt);
        }
        auto end_time_dec = std::chrono::high_resolution_clock::now();
        duration_out_then_dec = std::chrono::duration_cast<std::chrono::milliseconds>(end_time_dec - start_time_ctxt).count();
        std::cout << "serialAndDec duration = " << duration_out_then_dec << " ms" << std::endl;
        std::cout << std::endl;
    }
    
    // read kmerTable
    std::cout << "=== read kmerTable ===" << std::endl;
    auto start_time_kmerTable = std::chrono::high_resolution_clock::now();
    KmerTable kmerTableRef;
        if (param.json_format) {
        kmerTableRef.load(param.filename_kmerTable);
    } else {
        kmerTableRef.loadBinary(param.filename_kmerTable);
    }
    if (param.verbose) {
        kmerTableRef.print();
    }
    auto end_time_kmerTable = std::chrono::high_resolution_clock::now();
    auto duration_kmerTable = std::chrono::duration_cast<std::chrono::milliseconds>(end_time_kmerTable - start_time_kmerTable).count();
    std::cout << "load kmerTable duration = " << duration_kmerTable << " ms" << std::endl;
    std::cout << std::endl;

    cout << " === print gene vector === " << endl;
    for (size_t i = 0; i < final.size(); i++) {
        cout << kmerTableRef.geneNameIndex[i] << " : " << final[i] << endl;
    }
    cout << endl;
    
    auto total_duration = duration_kmerTable + duration_dec;
    std::cout << " === Duration summaries ===" << endl;
    std::cout << "load kmerTable duration = " << duration_kmerTable << " ms" << endl;
    if (param.operate_then_serialize) {
        std::cout << "read ctxtOut duration = " << duration_out << " ms" << endl;
        std::cout << "decryptOutput duration = " << duration_dec << " ms" << endl;
        total_duration = total_duration + duration_out + duration_dec;
    } else {
        std::cout << "serialAndDec duration = " << duration_out_then_dec << " ms" << endl;
        total_duration = total_duration + duration_out_then_dec;
    }
    std::cout << "total duration = " << total_duration << " ms" << std::endl;
    printMemoryUsage();
}
