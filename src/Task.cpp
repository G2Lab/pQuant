#include "Task.h"

void Task::readFastaFiles(PQuantParams &param) {
    // Read reference sequences from the file
    vector<Sequence> refs_seq;
    readFastaFile(param.path_filename_read, refs_seq);

    // Read read sequences from the file
    vector<Sequence> reads_seq;
    readFastaFile(param.path_filename_ref, reads_seq);

    long total_n_ref = 0;
    long n_gene_ref = 0;
    long total_n_read = 0;
    long n_gene_read = 0;

    long len_max = 0;
    long len_min = 999;
    long avg = 0;
    // Print some information about the reference sequences
    cout << "Reference sequences:" << endl;
    cout << param.path_filename_ref << endl;
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
    cout << param.path_filename_read << endl;
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
    string path_filename_ref = param.path_filename_ref;
    string path_filename_read = param.path_filename_read;

    vector<Sequence> refs_seq;
    vector<Sequence> reads_seq;
    readFastaFile(path_filename_ref, refs_seq);
    readFastaFile(path_filename_read, reads_seq);

    // cout << "check reads_Seq" << endl;
    // for (size_t i = 0; i < reads_seq.size(); i++) {
    //     cout << "reads_seq[" << i << "].getGeneName() = " << reads_seq[i].getGeneName() << endl;
    //     cout << "reads_seq[" << i << "].getNumSeq() = " << reads_seq[i].getNumSeq() << endl;
    //     for (size_t j = 0; j < static_cast<size_t>(reads_seq[i].getNumSeq()); j++) {
    //         cout << "reads_seq[" << i << "].getSeq(" << j << ") = " << reads_seq[i].getSeq(j) << endl;
    //     }
    // }

    KmerTable kmerTableRef(refs_seq, param, true);
    KmerTable kmerTableRead(reads_seq, param, false);

    cout << "geneNameIndex.size() = " << kmerTableRef.n_gene << endl;
    cout << "count.size() = " << kmerTableRef.count.size() << endl;
    cout << "entropy.size() = " << kmerTableRef.entropy.size() << endl;

    cout << "geneNameIndex.size() = " << kmerTableRead.geneNameIndex.size() << endl;
    cout << "count.size() = " << kmerTableRead.count.size() << endl;
    cout << "entropy.size() = " << kmerTableRead.entropy.size() << endl;

    cout << "=== kmerTableRef ===" << endl;
    // printKmerTable(kmerTableRef, true);
    if (param.verbose)
        kmerTableRef.print();

    cout << "=== kmerTableRead ===" << endl;
    // printKmerTable(kmerTableRead, false);
    // kmerTableRead.print();

    cout << "== save kmerTableRef to txt ==" << endl;
    kmerTableRef.save("kmerTableRef");

    cout << "print file size" << endl;
    ifstream file("kmerTableRef.txt");
    ifstream file2("kmerTableRef_kmer_list.txt");
    if (!file.is_open()) {
        throw runtime_error("Error: Failed to open file kmerTableRef");
    } else {
        file.seekg(0, ios::end);
        std::streampos filesize = file.tellg();
        double filesizeMB = static_cast<double>(filesize) / 1024 / 1024;
        // print filesize. If the file size is large, print as GB
        if (filesizeMB > 1024) {
            cout << "filesize = " << filesizeMB / 1024 << " GB" << endl;
        } else {
            cout << "filesize = " << filesizeMB << " MB" << endl;
        }
    }
    if (!file2.is_open()) {
        throw runtime_error("Error: Failed to open file kmerTableRef_kmer_list");
    } else {
        file2.seekg(0, ios::end);
        std::streampos filesize = file2.tellg();
        double filesizeMB = static_cast<double>(filesize) / 1024 / 1024;
        // print filesize. If the file size is large, print as GB
        if (filesizeMB > 1024) {
            cout << "filesize = " << filesizeMB / 1024 << " GB" << endl;
        } else {
            cout << "filesize = " << filesizeMB << " MB" << endl;
        }
    }

    cout << "== save kmerTableRef to bin ==" << endl;
    kmerTableRef.save_binary("kmerTableRef");

    cout << "print file size" << endl;
    ifstream file3("kmerTableRef.bin");
    if (!file3.is_open()) {
        throw runtime_error("Error: Failed to open file kmerTableRef.bin");
    } else {
        file3.seekg(0, ios::end);
        std::streampos filesize = file3.tellg();
        double filesizeMB = static_cast<double>(filesize) / 1024 / 1024;
        // print filesize. If the file size is large, print as GB
        if (filesizeMB > 1024) {
            cout << "filesize = " << filesizeMB / 1024 << " GB" << endl;
        } else {
            cout << "filesize = " << filesizeMB << " MB" << endl;
        }
    }
    ifstream file4("kmerTableRef_kmer_list.bin");
    if (!file4.is_open()) {
        throw runtime_error("Error: Failed to open file kmerTableRef_kmer_list.bin");
    } else {
        file4.seekg(0, ios::end);
        std::streampos filesize = file4.tellg();
        double filesizeMB = static_cast<double>(filesize) / 1024 / 1024;
        // print filesize. If the file size is large, print as GB
        if (filesizeMB > 1024) {
            cout << "filesize = " << filesizeMB / 1024 << " GB" << endl;
        } else {
            cout << "filesize = " << filesizeMB << " MB" << endl;
        }
    }

}

void Task::bfvBenchmark(PQuantParams &param) {
    TimeVar t;
    double processingTime(0.0);

    printMemoryUsage();

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

    printMemoryUsage();
}

void Task::run_all(PQuantParams &param) {
    auto start_time = std::chrono::high_resolution_clock::now();
    auto end_time = std::chrono::high_resolution_clock::now();
    
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
    auto n = cryptoContext->GetCryptoParameters()->GetElementParams()->GetCyclotomicOrder() / 2;
    std::cout << "n = " << n << endl;
    std::cout << "log2 q = "
              << log2(cryptoContext->GetCryptoParameters()->GetElementParams()->GetModulus().ConvertToDouble())
              << std::endl;

    // Initialize Public Key Containers
    KeyPair<DCRTPoly> keyPair = cryptoContext->KeyGen();
    cryptoContext->EvalMultKeysGen(keyPair.secretKey);

    KmerTable kmerTableRef;
    start_time = std::chrono::high_resolution_clock::now();
    if (param.path_kmer_matrix.size() > 0) {
        cout << "=== read kmerTableRef from json ===" << endl;
        // parseJson(param.path_kmer_matrix, kmerTableRef, param);
    } else {
        cout << "=== read refs_seq ===" << endl;
        vector<Sequence> refs_seq;
        readFastaFile(param.path_filename_ref, refs_seq);
        kmerTableRef = KmerTable(refs_seq, param, true);
        // computeKmerTable(refs_seq, param.k, kmerTableRef);
        refs_seq.clear();
    }
    if (param.verbose) {
        printKmerTable(kmerTableRef, true);
    }
    end_time = std::chrono::high_resolution_clock::now();
    auto duration_ref = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time).count();
    cout << "computeKmerTable duration = " << duration_ref << " ms" << endl;
    
    
    vector<Sequence> reads_seq;
    cout << "=== read reads_seq ===" << endl;
    readFastaFile(param.path_filename_read, reads_seq);
    // KmerTable kmerTableRead;
    cout << "=== compute kmerTableRead ===" << endl;
    // computeKmerTableForRead(reads_seq, param.k, kmerTableRead);
    KmerTable kmerTableRead(reads_seq, param, false);
    reads_seq.clear();
    if (param.verbose) {
        printKmerTable(kmerTableRead, false);
    }

    // remove all but x gene (debugging)
    if (param.debug_n_gene > 0) {
        kmerTableRef.n_gene = param.debug_n_gene;
        // erase all but debug_n_gene
        kmerTableRef.geneNameIndex.erase(std::next(kmerTableRef.geneNameIndex.begin(), param.debug_n_gene), kmerTableRef.geneNameIndex.end());

        cout << "Used gene list" << endl;
        for (size_t i = 0; i < kmerTableRef.geneNameIndex.size(); i++) {
            cout << kmerTableRef.geneNameIndex[i] << endl;
        }
    } else if (param.gene_start >= 0 && param.gene_end >= 0) {
        // erase all but gene_start to gene_end
        kmerTableRef.n_gene = param.gene_end - param.gene_start + 1;
        kmerTableRef.geneNameIndex.erase(std::next(kmerTableRef.geneNameIndex.begin(), param.gene_end + 1), kmerTableRef.geneNameIndex.end());
        kmerTableRef.geneNameIndex.erase(kmerTableRef.geneNameIndex.begin(), std::next(kmerTableRef.geneNameIndex.begin(), param.gene_start));
        cout << "Used gene list" << endl;
        for (size_t i = 0; i < kmerTableRef.geneNameIndex.size(); i++) {
            cout << kmerTableRef.geneNameIndex[i] << endl;
        }
    }
    
    // check kmerTable lengths
    cout << endl;
    cout << endl;
    cout << "kmerTableRef.count.size() = " << kmerTableRef.count.size() << endl;
    cout << "kmerTableRead.countRead.size() = " << kmerTableRead.countRead.size() << endl;
    cout << "kmerTableRef.entropy.size() = " << kmerTableRef.entropy.size() << endl;
    cout << "kmerRef.n_gene = " << kmerTableRef.n_gene << endl;
    cout << endl;
    
    start_time = std::chrono::high_resolution_clock::now();
    cout << endl;
    cout << " === run encryptReadKmer === " << endl;
    Ciphertext_1d ct;
    if (param.path_kmer_matrix.size() == 0) {
        cout << "run full table version" << endl;
        encryptReadKmer(kmerTableRead, ct, cryptoContext, keyPair, param);
    } else {
        cout << "run sparse version" << endl;
        encryptReadSparse(ct, kmerTableRead, kmerTableRef, cryptoContext, keyPair, param);
    }
    cout << endl;
    cout << "ct.size() = " << ct.size() << endl;
    cout << endl;
    end_time = std::chrono::high_resolution_clock::now();
    auto duration_encread = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time).count();
    cout << "encryptReadKmer duration = " << duration_encread << " ms" << endl;

    size_t duration_mult, duration_encode, duration_encode_and_mult;
    Ciphertext_1d ct_out;
    if (param.memory) {
        start_time = std::chrono::high_resolution_clock::now();
        cout << endl;
        cout << " === run encodeRefKmer === " << endl;
        Plaintext_2d pt_ref;
        encodeRefKmer(kmerTableRef, pt_ref, cryptoContext, keyPair, param);
        cout << endl;
        cout << "pt_ref.size() = " << pt_ref.size() << endl;
        cout << endl;
        end_time = std::chrono::high_resolution_clock::now();
        duration_encode = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time).count();
        cout << "encodeRefKmer duration = " << duration_encode << " ms" << endl;
        start_time = std::chrono::high_resolution_clock::now();
        cout << endl;
        cout << " === run multCtxtByRef === " << endl;
        multCtxtByEncodedRef(ct_out, ct, pt_ref, cryptoContext, param);
        cout << endl;
        cout << "ct_out.size() = " << ct_out.size() << endl;
        cout << endl;
        end_time = std::chrono::high_resolution_clock::now();
        duration_mult = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time).count();
        cout << "multCtxtByRef duration = " << duration_mult << " ms" << endl;
    } else {
        start_time = std::chrono::high_resolution_clock::now();
        cout << endl;
        cout << " === run multCtxtByKmerTableRefFromSerial === " << endl;
        multCtxtByKmerTableRef(ct_out, ct, kmerTableRef, cryptoContext, param);
        cout << endl;
        cout << "ct_out.size = " << ct_out.size() << endl;
        cout << endl;
        end_time = std::chrono::high_resolution_clock::now();
        duration_encode_and_mult = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time).count();
        cout << "multCtxtByRef duration = " << duration_mult << " ms" << endl;
    }
    
    

    start_time = std::chrono::high_resolution_clock::now();
    cout << endl;
    cout << " === run decrypt_output === " << endl;
    Plaintext_1d pt_sum(ct_out.size());
    for (size_t i = 0; i < ct_out.size(); i++) {
        cryptoContext->Decrypt(keyPair.secretKey, ct_out[i], &pt_sum[i]);
    }
    end_time = std::chrono::high_resolution_clock::now();
    auto duration_dec = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time).count();
    cout << endl;
    cout << "decrypt done" << endl;
    cout << endl;
    for (size_t i = 0; i < pt_sum.size(); i++) {
        auto plain = pt_sum[i]->GetCoefPackedValue();
        cout << kmerTableRef.geneNameIndex[i] << " : " << plain[0] << endl;
    }
    cout << endl;
    
    cout << "dec duration = " << duration_dec << " ms" << endl;

    // vector<Plaintext> pt_out;
    // decCtxtOut(pt_out, ct_out, cryptoContext, keyPair);
    // cout << "pt_out.size() = " << pt_out.size() << endl;

    cout << "== Outputs summary == " << endl;
    cout << "kmerTableRef.entropy.size() = " << kmerTableRef.entropy.size() << endl;
    printMemoryUsage();
    cout << endl;
    cout << "== duration summary ==" << endl;
    cout << "computeKmerTable duration = " << duration_ref << " ms" << endl;
    cout << "encryptReadKmer duration = " << duration_encread << " ms" << endl;
    if (param.memory) {
        cout << "encodeRefKmer duration   = " << duration_encode << " ms" << endl;
        cout << "multCtxtByRef duration   = " << duration_mult << " ms" << endl;
    } else {
        cout << "multCtxtByRef duration   = " << duration_encode_and_mult << " ms" << endl;
    }
    cout << "dec duration             = " << duration_dec << " ms" << endl;

}

// void Task::testReadJson(PQuantParams &param) {
//     // Specify the path to your JSON file
//     std::string filename = "../../pQuant_rust/kmer.json";
//     std::ifstream file(filename);
//     if (!file.is_open()) {
//         throw runtime_error("Error: Failed to open file " + filename);
//     } else {
//         KmerTable kmerTable;
//         parseJson(filename, kmerTable, param);
//         printKmerTable(kmerTable, true);
//     }
// }
