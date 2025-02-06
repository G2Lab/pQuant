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

    vector<Sequence> refs_seq;
    vector<Sequence> reads_seq;
    std::cout << "read reference sequences" << endl;
    readFastaFile(filename_ref, refs_seq);
    std::cout << "read read sequences" << endl;
    if (param.filename_read.find(".fa") != std::string::npos && param.filename_read.find(".fastq") == std::string::npos) {
        std::cout << "confirmed filename ends with .fa" << endl;
        readFastaFile(param.filename_read, reads_seq);
    } else if (param.filename_read.find(".fq") != std::string::npos || param.filename_read.find(".fastq") != std::string::npos) {
        std::cout << "confirmed filename ends with .fq or .fastq" << endl;
        readFastQFile(param.filename_read, reads_seq);
    } else {
        std::cout << "Invalid file format" << std::endl;
        std::cout << "Assume it is .fa format" << std::endl;
        readFastaFile(param.filename_read, reads_seq);
    }

    KmerTable kmerTableRef(refs_seq, param, true);
    KmerTable kmerTableRead(reads_seq, param, false);

    cout << "geneNameIndex.size() = " << kmerTableRef.n_gene << endl;
    cout << "count.size() = " << kmerTableRef.count.size() << endl;
    cout << "entropy.size() = " << kmerTableRef.entropy.size() << endl;

    cout << "geneNameIndex.size() = " << kmerTableRead.geneNameIndex.size() << endl;
    cout << "count.size() = " << kmerTableRead.count.size() << endl;
    cout << "entropy.size() = " << kmerTableRead.entropy.size() << endl;

    cout << "=== kmerTableRef ===" << endl;
    // if (param.verbose)
    //     kmerTableRef.print();

    cout << "=== kmerTableRead ===" << endl;
    if (param.verbose)
        kmerTableRead.print();
    return;

    cout << "== save kmerTableRef to txt ==" << endl;
    kmerTableRef.save("kmerTableRef.txt");
    kmerTableRef.saveKmerList("kmerTableRef_kmer_list.txt");

    cout << "print file size" << endl;
    printFileSize("kmerTableRef.txt");
    printFileSize("kmerTableRef_kmer_list.txt");

    cout << "print file size" << endl;

    // load and check ==
    KmerTable kmerTableRef2;
    kmerTableRef2.load("kmerTableRef.txt");
    vector<size_t> kmer_list;
    loadKmerList("kmerTableRef_kmer_list.txt", kmer_list);

    cout << "check kmer list ==" << endl;
    int i = 0;
    for (auto p: kmerTableRef.entropy) {
        if (p.first != kmer_list[i]) {
            cout << "error at " << i << endl;
            cout << p.first << " != " << kmer_list[i] << endl;
            return;
        }
        i += 1;
    }
    cout << "done" << endl;
    cout << "== save kmerTableRef to bin ==" << endl;
    kmerTableRef.saveBinary("kmerTableRef.bin");
    kmerTableRef.saveKmerListBinary("kmerTableRef_kmer_list.bin");

    cout << "== load kmerTableRef from bin ==" << endl;
    KmerTable kmerTableRef3;
    kmerTableRef3.loadBinary("kmerTableRef.bin");
    vector<size_t> kmer_list3;
    loadKmerListBinary("kmerTableRef_kmer_list.bin", kmer_list3);
    cout << "check kmer list ==" << endl;
    i = 0;
    for (auto p: kmerTableRef.entropy) {
        if (param.verbose) {
            cout << "p.first vs kmre_list3[" << i << "] = " << p.first << " - " << kmer_list3[i] << " (length = " << kmerTableRef.entropy.size() << ", " << kmer_list3.size() << endl;     
        }
        if (p.first != kmer_list3[i]) {
            cout << "error at " << i << endl;
            cout << p.first << " != " << kmer_list3[i] << endl;
            return;
        }
        i += 1;
    }
}

void Task::bfvBenchmark(PQuantParams &param) {
    TimeVar t;
    double processingTime(0.0);

    printMemoryUsage();

    CCParams<CryptoContextBFVRNS> parameters;
    parameters.SetSecurityLevel(HEStd_128_classic);
    parameters.SetPlaintextModulus(536903681);
    parameters.SetMultiplicativeDepth(2);
    parameters.SetMaxRelinSkDeg(2);
    // parameters.SetScalingModSize(20);

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
    long iteration = 10;
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
        // plain = cryptoContext->MakePackedPlaintext(plain_vec);
    }
    processingTime = TOC(t);
    std::cout << "encode time: " << processingTime / iteration << "ms" << std::endl;
    // DCRTPoly poly = plain->GetElement<DCRTPoly>();
    // std::stringstream s;
    // DCRTPoly deser;
    // Serial::Serialize(poly, s, SerType::JSON);
    // cout << "serializetion confirmed" << endl;
    // Serial::Deserialize(deser, s, SerType::JSON);
    // cout << "deserializetion confirmed" << endl;
    // std::ofstream f;
    // f.open("test.txt");
    // f << s.str();
    // f.close();
    // cout << "save confirmed" << endl;
    // cout << "plain" << endl;
    // cout << poly << endl;
    // cout << "deserialized" << endl;
    // cout << deser << endl;
    // Plaintext plain_deser = CryptoContext::MakePlaintext(PACKED_ENCODING, cryptoContext, deser);

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

    KmerTable *kmerTableGen;
    start_time = std::chrono::high_resolution_clock::now();
    if (param.filename_kmerTable.size() > 0) {
        cout << "=== read kmerTableRef from json ===" << endl;
        kmerTableGen = new KmerTable();
        kmerTableGen->load(param.filename_kmerTable);
    } else {
        cout << "=== read refs_seq ===" << endl;
        vector<Sequence> refs_seq;
        readFastaFile(param.filename_ref, refs_seq);
        kmerTableGen = new KmerTable(refs_seq, param, true);
        refs_seq.clear();
    }
    if (param.verbose) {
        kmerTableGen->print();
    }
    end_time = std::chrono::high_resolution_clock::now();
    auto duration_ref = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time).count();
    cout << "computeKmerTable duration = " << duration_ref << " ms" << endl;
    
    cout << " === Save & Load === " << endl;
    cout << "save kmerTableGen to txt" << endl;
    if (param.json_format) {
        kmerTableGen->save("kmerTableGen.txt");
        kmerTableGen->saveKmerList("kmerTableGen_kmer_list.txt");
    } else {
        kmerTableGen->saveBinary("kmerTableGen.txt");
        kmerTableGen->saveKmerListBinary("kmerTableGen_kmer_list.txt");
    }
    

    cout << "print file size" << endl;
    printFileSize("kmerTableGen.txt");
    printFileSize("kmerTableGen_kmer_list.txt");

    cout << "load kmerTableRef from txt" << endl;
    KmerTable *kmerTableRef;
    kmerTableRef = new KmerTable();
    vector<size_t> kmer_list;
    if (param.json_format) {
        kmerTableRef->load("kmerTableGen.txt");
        loadKmerList("kmerTableGen_kmer_list.txt", kmer_list);
    } else {
        kmerTableRef->loadBinary("kmerTableGen.txt");
        loadKmerListBinary("kmerTableGen_kmer_list.txt", kmer_list);
    }
    

    cout << "check kmer list ==" << endl;
    int i = 0;
    for (auto p: kmerTableRef->entropy) {
        if (p.first != kmer_list[i]) {
            cout << "error at " << i << endl;
            cout << p.first << " != " << kmer_list[i] << endl;
            return;
        }
        i += 1;
    }

    
    vector<Sequence> reads_seq;
    cout << "=== read reads_seq ===" << endl;
    readFastaFile(param.filename_read, reads_seq);
    // KmerTable kmerTableRead;
    cout << "=== compute kmerTableRead ===" << endl;
    // computeKmerTableForRead(reads_seq, param.k, kmerTableRead);
    KmerTable kmerTableRead(reads_seq, param, false);
    reads_seq.clear();
    if (param.verbose) {
        kmerTableRead.print();
    }

    // remove all but x gene (debugging)
    if (param.debug_n_gene > 0) {
        kmerTableRef->n_gene = param.debug_n_gene;
        // erase all but debug_n_gene
        kmerTableRef->geneNameIndex.erase(std::next(kmerTableRef->geneNameIndex.begin(), param.debug_n_gene), kmerTableRef->geneNameIndex.end());

        cout << "Used gene list" << endl;
        for (size_t i = 0; i < kmerTableRef->geneNameIndex.size(); i++) {
            cout << kmerTableRef->geneNameIndex[i] << endl;
        }
    } else if (param.gene_start >= 0 && param.gene_end >= 0) {
        // erase all but gene_start to gene_end
        kmerTableRef->n_gene = param.gene_end - param.gene_start + 1;
        kmerTableRef->geneNameIndex.erase(std::next(kmerTableRef->geneNameIndex.begin(), param.gene_end + 1), kmerTableRef->geneNameIndex.end());
        kmerTableRef->geneNameIndex.erase(kmerTableRef->geneNameIndex.begin(), std::next(kmerTableRef->geneNameIndex.begin(), param.gene_start));
        cout << "Used gene list" << endl;
        for (size_t i = 0; i < kmerTableRef->geneNameIndex.size(); i++) {
            cout << kmerTableRef->geneNameIndex[i] << endl;
        }
    }
    
    // check kmerTable lengths
    cout << endl;
    cout << endl;
    cout << "kmerTableRef->count.size() = " << kmerTableRef->count.size() << endl;
    cout << "kmerTableRead.countRead.size() = " << kmerTableRead.countRead.size() << endl;
    cout << "kmerTableRef->entropy.size() = " << kmerTableRef->entropy.size() << endl;
    cout << "kmerRef.n_gene = " << kmerTableRef->n_gene << endl;
    cout << endl;
    
    start_time = std::chrono::high_resolution_clock::now();
    cout << endl;
    cout << " === run encryptReadKmer === " << endl;
    Ciphertext_1d ct;
    if (param.filename_kmerTable.size() == 0) {
        cout << "run full table version" << endl;
        encryptReadKmer(kmerTableRead, ct, cryptoContext, keyPair, param);
    } else {
        cout << "run sparse version" << endl;
        encryptReadSparse(ct, kmerTableRead, *kmerTableRef, cryptoContext, keyPair, param);
    }
    cout << endl;
    cout << "ct.size() = " << ct.size() << endl;
    cout << endl;
    end_time = std::chrono::high_resolution_clock::now();
    auto duration_encread = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time).count();
    cout << "encryptReadKmer duration = " << duration_encread << " ms" << endl;

    size_t duration_mult, duration_encode, duration_encode_and_mult;
    Ciphertext_1d ct_out;
    if (param.operate_then_serialize) {
        start_time = std::chrono::high_resolution_clock::now();
        cout << endl;
        cout << " === run encodeRefKmer === " << endl;
        Plaintext_2d pt_ref;
        encodeRefKmer(*kmerTableRef, pt_ref, cryptoContext, keyPair, param);
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
        multCtxtByKmerTableRef(ct_out, ct, *kmerTableRef, cryptoContext, param);
        cout << endl;
        cout << "ct_out.size = " << ct_out.size() << endl;
        cout << endl;
        end_time = std::chrono::high_resolution_clock::now();
        duration_encode_and_mult = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time).count();
        cout << "multCtxtByRef duration = " << duration_encode_and_mult << " ms" << endl;
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
        cout << kmerTableRef->geneNameIndex[i] << " : " << plain[0] << endl;
    }
    cout << endl;
    
    cout << "dec duration = " << duration_dec << " ms" << endl;

    // vector<Plaintext> pt_out;
    // decCtxtOut(pt_out, ct_out, cryptoContext, keyPair);
    // cout << "pt_out.size() = " << pt_out.size() << endl;

    cout << "== Outputs summary == " << endl;
    cout << "kmerTableRef->entropy.size() = " << kmerTableRef->entropy.size() << endl;
    printMemoryUsage();
    cout << endl;
    cout << "== duration summary ==" << endl;
    cout << "computeKmerTable duration = " << duration_ref << " ms" << endl;
    cout << "encryptReadKmer duration = " << duration_encread << " ms" << endl;
    if (param.operate_then_serialize) {
        cout << "encodeRefKmer duration   = " << duration_encode << " ms" << endl;
        cout << "multCtxtByRef duration   = " << duration_mult << " ms" << endl;
    } else {
        cout << "multCtxtByRef duration   = " << duration_encode_and_mult << " ms" << endl;
    }
    cout << "dec duration             = " << duration_dec << " ms" << endl;

}

void Task::testSimulatedReadEncoding(PQuantParams &param) {
    vector<Sequence> seq_vec;
    auto start_time_readFastaFile = std::chrono::high_resolution_clock::now();
    cout << "=== read fasta file ===" << endl;
    readFastaFile(param.filename_ref, seq_vec);
    auto end_time_readFastaFile = std::chrono::high_resolution_clock::now();
    auto duration_readFastaFile = std::chrono::duration_cast<std::chrono::milliseconds>(end_time_readFastaFile - start_time_readFastaFile).count();
    std::cout << "readFastaFile duration = " << duration_readFastaFile << " ms" << std::endl;
    std::cout << std::endl;

    // simulate read
    auto start_time_simulateRead = std::chrono::high_resolution_clock::now();
    cout << "=== simulate read ===" << endl;
    vector<Sequence> reads_seq;

    generateSimulatedReads(seq_vec, param.sim_len, param.sim_num, param.sim_err, reads_seq);
    auto end_time_simulateRead = std::chrono::high_resolution_clock::now();
    auto duration_simulateRead = std::chrono::duration_cast<std::chrono::milliseconds>(end_time_simulateRead - start_time_simulateRead).count();
    std::cout << "simulateRead duration = " << duration_simulateRead << " ms" << std::endl;
    std::cout << std::endl;

    // // print simulatedRead to check
    // cout << "=== print simulated read ===" << endl;
    for (auto &read_seq : reads_seq) {
        if (read_seq.getNumSeq() > 0) {
            cout << "Gene name: " << read_seq.getGeneName() << ", Num seq: " << read_seq.getNumSeq() << endl;
            for (size_t i = 0; i < static_cast<size_t>(read_seq.getNumSeq()); i++) {
                cout << "   " << read_seq.getSeq(i) << endl;
            }
        }
    }
    
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

    auto start_time_read = std::chrono::high_resolution_clock::now();
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
        if (search) {
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