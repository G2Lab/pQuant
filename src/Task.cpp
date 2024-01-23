#include "Task.h"

void Task::readFastaFiles(PQuantParams &param) {
    // Read reference sequences from the file
    vector<Sequence> refs_seq;
    readFastaFile(param.filename_read, refs_seq);

    // Read read sequences from the file
    vector<Sequence> reads_seq;
    readFastaFile(param.filename_ref, reads_seq);

    long total_n_ref = 0;
    long n_gene_ref = 0;
    long total_n_read = 0;
    long n_gene_read = 0;

    long len_max = 0;
    long len_min = 999;
    long avg = 0;
    // Print some information about the reference sequences
    cout << "Reference sequences:" << endl;
    cout << param.filename_ref << endl;
    for (auto &ref_seq : refs_seq) {
        cout << "Gene name: " << ref_seq.getGeneName() << ", Num seq: " << ref_seq.getNumSeq() << endl;
        n_gene_ref += 1;
        for (size_t i = 0; i < static_cast<size_t>(ref_seq.getNumSeq()); i++) {
            // cout << "  -  " << ref_seq.getSeq(i) << endl;
            total_n_ref += 1;
            long len = ref_seq.getSeq(i).size();
            if (len_max < len)
                len_max = len;
            if (len_min > len)
                len_min = len;
            avg += len;
        }
    }

    // Print some information about the read sequences
    cout << "Read sequences:" << endl;
    cout << param.filename_read << endl;
    for (auto &read_seq : reads_seq) {
        // cout << "Gene name: " << read_seq.getGeneName() << ", Num seq: " << read_seq.getNumSeq() << endl;
        n_gene_read += 1;
        for (size_t i = 0; i < static_cast<size_t>(read_seq.getNumSeq()); i++) {
            // cout << "   " << read_seq.getSeq(i) << endl;
            total_n_read += 1;
        }
    }
    cout << "n_ref = " << total_n_ref << endl;
    cout << "n_read = " << total_n_read << endl;
    cout << "n_gene_ref = " << n_gene_ref << endl;
    cout << "n_gene_read = " << n_gene_read << endl;

    cout << "max = " << len_max << endl;
    cout << "min = " << len_min << endl;
    cout << "avg = " << avg / total_n_ref << endl;
}

void Task::testKmerTable(PQuantParams &param) {
    string filename_ref = param.filename_ref;
    string filename_read = param.filename_read;
    long K = param.k;

    vector<Sequence> refs_seq;
    vector<Sequence> reads_seq;
    readFastaFile(filename_ref, refs_seq);
    readFastaFile(filename_read, reads_seq);

    cout << "check reads_Seq" << endl;
    for (size_t i = 0; i < reads_seq.size(); i++) {
        cout << "reads_seq[" << i << "].getGeneName() = " << reads_seq[i].getGeneName() << endl;
        cout << "reads_seq[" << i << "].getNumSeq() = " << reads_seq[i].getNumSeq() << endl;
        for (size_t j = 0; j < static_cast<size_t>(reads_seq[i].getNumSeq()); j++) {
            cout << "reads_seq[" << i << "].getSeq(" << j << ") = " << reads_seq[i].getSeq(j) << endl;
        }
    }

    KmerTable kmerTableRef;
    KmerTable kmerTableRead;
    computeKmerTable(refs_seq, K, kmerTableRef);
    computeKmerTableForRead(reads_seq, K, kmerTableRead);

    cout << "geneNameIndex.size() = " << kmerTableRef.n_gene << endl;
    cout << "count.size() = " << kmerTableRef.count.size() << endl;
    cout << "entropy.size() = " << kmerTableRef.entropy.size() << endl;

    cout << "geneNameIndex.size() = " << kmerTableRead.geneNameIndex.size() << endl;
    cout << "count.size() = " << kmerTableRead.count.size() << endl;
    cout << "entropy.size() = " << kmerTableRead.entropy.size() << endl;

    cout << "=== kmerTableRef ===" << endl;
    printKmerTable(kmerTableRef, true);

    cout << "=== kmerTableRead ===" << endl;
    printKmerTable(kmerTableRead, false);
}

void Task::bfvBenchmark(PQuantParams &param) {
    TimeVar t;
    double processingTime(0.0);

    size_t memory_usage = getMemoryUsage();

    // Print memory usage in a human-readable format
    double memory_usage_gb = static_cast<double>(memory_usage) / (1024.0 * 1024.0 * 1024.0); // Convert to GB

    std::cout << "Memory Usage: " << memory_usage_gb << " GB" << std::endl;

    CCParams<CryptoContextBFVRNS> parameters;
    parameters.SetSecurityLevel(HEStd_128_classic);
    parameters.SetPlaintextModulus(65537);
    parameters.SetMultiplicativeDepth(2);
    parameters.SetMaxRelinSkDeg(2);
    parameters.SetScalingModSize(20);

    CryptoContext<DCRTPoly> cryptoContext = GenCryptoContext(parameters);
    // enable features that you wish to use
    cryptoContext->Enable(PKE);
    cryptoContext->Enable(KEYSWITCH);
    cryptoContext->Enable(LEVELEDSHE);
    cryptoContext->Enable(ADVANCEDSHE);


    std::cout << "\np = " << cryptoContext->GetCryptoParameters()->GetPlaintextModulus() << std::endl;
    std::cout << "n = " << cryptoContext->GetCryptoParameters()->GetElementParams()->GetCyclotomicOrder() / 2
              << std::endl;
    std::cout << "log2 q = "
              << log2(cryptoContext->GetCryptoParameters()->GetElementParams()->GetModulus().ConvertToDouble())
              << std::endl;

    // Initialize Public Key Containers
    KeyPair<DCRTPoly> keyPair = cryptoContext->KeyGen();
    cryptoContext->EvalMultKeysGen(keyPair.secretKey);
    Plaintext_3d ref_plain;
    Plaintext_4d read_plain;
    long slot_count = cryptoContext->GetCryptoParameters()->GetElementParams()->GetCyclotomicOrder() / 2;
    long log_poly_modulus_size = (long)log2(slot_count);
    cout << "log_slot_count = " << log_poly_modulus_size << endl;
    vector<int32_t> rotList(log_poly_modulus_size);


    // benchmark each operations
    long iteration = 1000;
    cout << "<<< every algorithm is iterated " << iteration << " times >>> " << endl << endl;

    // generate random plaintexts
    vector<int64_t> plain_vec(slot_count, 0);
    for (int i = 0; i < slot_count; i++) {
        plain_vec[i] = rand() % 2;
    }

    // encode
    Plaintext plain;
    TIC(t);
    for (int i = 0; i < iteration; i++) {
        plain = cryptoContext->MakeCoefPackedPlaintext(plain_vec);
    }
    processingTime = TOC(t);
    std::cout << "encode time: " << processingTime / iteration << "ms" << std::endl;

    //encrypt
    Ciphertext<DCRTPoly> ciphertext, ciphertextMult;
    TIC(t);
    for (int i = 0; i < iteration; i++) {
        ciphertext = cryptoContext->Encrypt(keyPair.secretKey, plain);
        ciphertextMult = cryptoContext->Encrypt(keyPair.secretKey, plain);
    }
    processingTime = TOC(t);
    std::cout << "Encrypt time: " << processingTime / iteration << "ms" << std::endl;

    // decrypt
    Plaintext plain2;
    TIC(t);
    for (int i = 0; i < iteration; i++) {
        cryptoContext->Decrypt(keyPair.secretKey, ciphertext, &plain2);
    }
    processingTime = TOC(t);
    std::cout << "Decrypt time: " << processingTime / iteration << "ms" << std::endl;

    // mult ctxt * plain
    Ciphertext<DCRTPoly> ciphertext2;
    TIC(t);
    for (int i = 0; i < iteration; i++) {
        ciphertext2 = cryptoContext->EvalMult(ciphertext, plain);
    }
    processingTime = TOC(t);
    std::cout << "Mult ctxt * plain time: " << processingTime / iteration << "ms" << std::endl;

    // mult ctxt * vector
    Ciphertext<DCRTPoly> ciphertext4;
    TIC(t);
    for (int i = 0; i < iteration; i++) {
        Plaintext plain = cryptoContext->MakeCoefPackedPlaintext(plain_vec);
        ciphertext4 = cryptoContext->EvalMult(ciphertext, plain);
    }
    processingTime = TOC(t);
    std::cout << "Mult ctxt * vector (including encoding) time: " << processingTime / iteration << "ms" << std::endl;

    // mult ctxt * ctxt
    TIC(t);
    for (int i = 0; i < iteration; i++) {
        ciphertext2 = cryptoContext->EvalMult(ciphertext, ciphertextMult);
    }
    processingTime = TOC(t);
    std::cout << "Mult ctxt * ctxt time: " << processingTime / iteration << "ms" << std::endl;
    // add
    Ciphertext<DCRTPoly> ciphertext3, ciphertextAdd;
    ciphertextAdd = cryptoContext->Encrypt(keyPair.secretKey, plain);
    TIC(t);
    for (int i = 0; i < iteration; i++) {
        ciphertext3 = cryptoContext->EvalAdd(ciphertext, ciphertextAdd);
    }
    processingTime = TOC(t);
    std::cout << "Add time: " << processingTime / iteration << "ms" << std::endl;

    memory_usage = getMemoryUsage();

    // Print memory usage in a human-readable format
    memory_usage_gb = static_cast<double>(memory_usage) / (1024.0 * 1024.0 * 1024.0); // Convert to GB

    std::cout << "Memory Usage: " << memory_usage_gb << " GB" << std::endl;
}

void Task::run_all(PQuantParams &param) {
    // Read read sequences from the file
    vector<Sequence> reads_seq, refs_seq;
    readFastaFile(param.filename_read, reads_seq);
    readFastaFile(param.filename_ref, refs_seq);

    CCParams<CryptoContextBFVRNS> parameters;
    parameters.SetPlaintextModulus(65537);
    parameters.SetMultiplicativeDepth(1);
    parameters.SetMaxRelinSkDeg(2);

    CryptoContext<DCRTPoly> cryptoContext = GenCryptoContext(parameters);
    // enable features that you wish to use
    cryptoContext->Enable(PKE);
    cryptoContext->Enable(KEYSWITCH);
    cryptoContext->Enable(LEVELEDSHE);
    cryptoContext->Enable(ADVANCEDSHE);

    std::cout << "\np = " << cryptoContext->GetCryptoParameters()->GetPlaintextModulus() << std::endl;
    std::cout << "n = " << cryptoContext->GetCryptoParameters()->GetElementParams()->GetCyclotomicOrder() / 2
              << std::endl;
    auto n = cryptoContext->GetCryptoParameters()->GetElementParams()->GetCyclotomicOrder() / 2;
    std::cout << "log2 q = "
              << log2(cryptoContext->GetCryptoParameters()->GetElementParams()->GetModulus().ConvertToDouble())
              << std::endl;

    // Initialize Public Key Containers
    KeyPair<DCRTPoly> keyPair = cryptoContext->KeyGen();
    cryptoContext->EvalMultKeysGen(keyPair.secretKey);

    KmerTable kmerTableRef;
    computeKmerTable(refs_seq, param.k, kmerTableRef);
    KmerTable kmerTableRead;
    computeKmerTableForRead(reads_seq, param.k, kmerTableRead);

    if (param.verbose) {
        printKmerTable(kmerTableRef, true);
        printKmerTable(kmerTableRead, false);
    }

    // remove all but x gene (debugging)
    if (param.debug_n_gene > 0) {
        kmerTableRef.n_gene = param.debug_n_gene;
        // erase all but debug_n_gene
        kmerTableRef.geneNameIndex.erase(std::next(kmerTableRef.geneNameIndex.begin(), param.debug_n_gene),
                                         kmerTableRef.geneNameIndex.end());
    }
    

    
    // check kmerTable lengths
    cout << endl;
    cout << endl;
    cout << "kmerTableRef.count.size() = " << kmerTableRef.count.size() << endl;
    cout << "kmerTableRead.countRead.size() = " << kmerTableRead.countRead.size() << endl;
    cout << "kmerTableRef.entropy.size() = " << kmerTableRef.entropy.size() << endl;
    cout << "kmerRef.n_gene = " << kmerTableRef.n_gene << endl;
    cout << endl;

    auto start_time = std::chrono::high_resolution_clock::now();
    auto end_time = std::chrono::high_resolution_clock::now();
    
    start_time = std::chrono::high_resolution_clock::now();
    cout << endl;
    cout << " === run encryptReadKmer === " << endl;
    Ciphertext_1d ct;
    encryptReadKmer(kmerTableRead, param.k, ct, cryptoContext, keyPair);
    cout << endl;
    cout << "ct.size() = " << ct.size() << endl;
    cout << endl;
    end_time = std::chrono::high_resolution_clock::now();
    auto duration_encread = std::chrono::duration_cast<std::chrono::seconds>(end_time - start_time).count();
    cout << "encryptReadKmer duration = " << duration_encread << " s" << endl;

    start_time = std::chrono::high_resolution_clock::now();
    cout << endl;
    cout << " === run multCtxtByRef === " << endl;
    Ciphertext_2d ct_out;
    if (param.serial) {
        multCtxtByKmerTableRefFromSerial(ct_out, ct, kmerTableRef, param.k, cryptoContext, param.out_path);
    } else {
        multCtxtByKmerTableRef2(ct_out, ct, kmerTableRef, param.k, cryptoContext);
    }
    cout << endl;
    cout << "ct_out.size() = " << ct_out.size() << endl;
    cout << endl;
    end_time = std::chrono::high_resolution_clock::now();
    auto duration_mult = std::chrono::duration_cast<std::chrono::seconds>(end_time - start_time).count();
    cout << "multCtxtByRef duration = " << duration_mult << " s" << endl;

    start_time = std::chrono::high_resolution_clock::now();
    cout << endl;
    cout << " === run sumUpCtxt === " << endl;
    Ciphertext_1d ct_sum;
    if (param.serial) {
        long n_gene = kmerTableRef.n_gene;
        long n_ctxt = pow(4, param.k) / n;
        sumUpCtxtFromSerial(ct_sum, n_gene, n_ctxt, cryptoContext, param.out_path);
    } else {
        ct_sum.resize(ct_out.size());
        for (size_t i = 0; i < ct_out.size(); i++) {
            sumUpCtxt(ct_sum[i], ct_out[i], cryptoContext);
        }
    }
    
    cout << endl;
    cout << endl;
    end_time = std::chrono::high_resolution_clock::now();
    auto duration_sum = std::chrono::duration_cast<std::chrono::seconds>(end_time - start_time).count();
    cout << "sumUpCtxt duration = " << duration_sum << " s" << endl;

    start_time = std::chrono::high_resolution_clock::now();
    cout << endl;
    cout << " === run decrypt_output === " << endl;
    Plaintext_1d pt_sum(ct_sum.size());
    for (size_t i = 0; i < ct_sum.size(); i++) {
        cryptoContext->Decrypt(keyPair.secretKey, ct_sum[i], &pt_sum[i]);
    }
    cout << endl;
    cout << "decrypt done" << endl;
    cout << endl;
    for (size_t i = 0; i < pt_sum.size(); i++) {
        auto plain = pt_sum[i]->GetCoefPackedValue();
        cout << plain[0] << endl;
        // cout << plain[0] << endl;
        // cout << "pt_sum[" << i << "] = " << plain[0] << endl;
        // cout << plain << endl;
    }
    cout << endl;
    end_time = std::chrono::high_resolution_clock::now();
    auto duration_dec = std::chrono::duration_cast<std::chrono::seconds>(end_time - start_time).count();
    cout << "dec duration = " << duration_dec << " s" << endl;

    // vector<Plaintext> pt_out;
    // decCtxtOut(pt_out, ct_out, cryptoContext, keyPair);
    // cout << "pt_out.size() = " << pt_out.size() << endl;
    cout << "== duration summary ==" << endl;
    cout << "encryptReadKmer duration = " << duration_encread << " s" << endl;
    cout << "multCtxtByRef duration   = " << duration_mult << " s" << endl;
    cout << "sumUpCtxt duration       = " << duration_sum << " s" << endl;
    cout << "dec duration             = " << duration_dec << " s" << endl;
}

void Task::testReadJson(PQuantParams &param) {
    // Specify the path to your JSON file
    std::string filename = "../kmer.json";
    std::ifstream file(filename);
    if (!file.is_open()) {
        throw runtime_error("Error: Failed to open file " + filename);
    } else {
        std::stringstream buffer;
        buffer << file.rdbuf();
        std::string jsonStr = buffer.str();

        KmerTable kmerTable;
        parseJson(filename, kmerTable);
        printKmerTable(kmerTable, true);
    }
}
