#include "test/Test.h"


size_t getMemoryUsage() {
    size_t memory_usage = 0;

    struct rusage r_usage;
    getrusage(RUSAGE_SELF, &r_usage);
    memory_usage = static_cast<size_t>(r_usage.ru_maxrss) * 1024; // Convert to bytes

    return memory_usage;
}

void Test::testReadFastaFiles(const string &ref_filename, const string &read_filename) {
    // Read reference sequences from the file
    vector<Sequence> refs_seq;
    readFastaFile(ref_filename, refs_seq);

    // Read read sequences from the file
    vector<Sequence> reads_seq;
    readFastaFile(read_filename, reads_seq);

    long total_n_ref = 0;
    long n_gene_ref = 0;
    long total_n_read = 0;
    long n_gene_read = 0;

    long len_max = 0;
    long len_min = 999;
    long avg = 0;
    // Print some information about the reference sequences
    cout << "Reference sequences:" << endl;
    cout << ref_filename << endl;
    for (auto &ref_seq : refs_seq) {
        cout << "Gene name: " << ref_seq.getGeneName() << ", Num seq: " << ref_seq.getNumSeq() << endl;
        n_gene_ref += 1;
        for (int i = 0; i < ref_seq.getNumSeq(); i++) {
            // cout << "  -  " << ref_seq.getSeq(i) << endl;
            total_n_ref += 1;
            long len = ref_seq.getSeq(i).size();
            if (len_max < len) len_max = len;
            if (len_min > len) len_min = len;
            avg += len;
        }
    }

    // Print some information about the read sequences
    cout << "Read sequences:" << endl;
    cout << read_filename << endl;
    for (auto &read_seq : reads_seq) {
        // cout << "Gene name: " << read_seq.getGeneName() << ", Num seq: " << read_seq.getNumSeq() << endl;
        n_gene_read += 1;
        for (int i = 0; i < read_seq.getNumSeq(); i++) {
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


void Test::previousAlgorithm() {
    TimeVar t;
    double processingTime(0.0);

    string filename_read = "../dataset/five_genes/five_gene_reads_shorten.fa";
    string filename_ref = "../dataset/five_genes/five_gene_reference.fa";
    long K = 32;
    long B = 8;

    vector<Sequence> refs_seq;
    readFastaFile(filename_ref, refs_seq);

    // Read read sequences from the file
    vector<Sequence> reads_seq;
    readFastaFile(filename_read, reads_seq);

    // benchmarking variables
    // TimeVar t;
    // double processingTime(0.0);

    CCParams<CryptoContextBFVRNS> parameters;
    parameters.SetPlaintextModulus(65537);
    parameters.SetMultiplicativeDepth(18);
    parameters.SetMaxRelinSkDeg(2);
    parameters.SetScalingModSize(40);

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
    for (int i = 0; i < rotList.size(); i++) {
        rotList[i] = (1 << i);
    }
    cryptoContext->EvalRotateKeyGen(keyPair.secretKey, rotList);

    cout << "=== Encode Refs ===" << endl;
    TIC(t);
    encode_refs(cryptoContext, refs_seq, ref_plain, K, B, slot_count);
    processingTime = TOC(t);
    std::cout << "Encoding ref time: " << processingTime << "ms" << std::endl;
    cout << "=== Encode reads ===" << endl;
    TIC(t);
    encode_read(cryptoContext, reads_seq, read_plain, K, B, slot_count);
    processingTime = TOC(t);
    std::cout << "Encoding read time: " << processingTime << "ms" << std::endl;
    TIC(t);
    Ciphertext_4d ctxt_kmers;
    cout << "=== Encrypt reads ===" << endl;
    encrypt_reads(cryptoContext, keyPair, read_plain, ctxt_kmers, K, B, slot_count);
    processingTime = TOC(t);
    std::cout << "Encrypt Reads time: " << processingTime << "ms" << std::endl;
    Ciphertext_3d ctxt_out;
    cout << "=== Equality Test ===" << endl;
    TIC(t);
    eqaulity_test(cryptoContext, keyPair, ctxt_kmers, ref_plain, ctxt_out, K, B, log_poly_modulus_size, refs_seq, reads_seq);
    processingTime = TOC(t);
    std::cout << "EQtest time: " << processingTime << "ms" << std::endl;
    cout << "=== Decrypt ===" << endl;
    TIC(t);
    decrypt_output(cryptoContext, keyPair, ctxt_out, refs_seq, reads_seq, K, B);
    processingTime = TOC(t);
    std::cout << "Decrypt time: " << processingTime << "ms" << std::endl;
}

void Test::bfvBasicAlgorithmBenchmarks() {
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
    
    // generate random plaintexts
    vector<int64_t> plain_vec(slot_count, 0);
    for (int i = 0; i < slot_count; i++) {
        plain_vec[i] = rand() % 2;
    }

    // encode
    Plaintext plain;
    TIC(t);
    for (int i = 0; i < 1000; i++) {
        plain = cryptoContext->MakeCoefPackedPlaintext(plain_vec);
    }
    processingTime = TOC(t);
    std::cout << "encode time: " << processingTime / 1000 << "ms" << std::endl;

    //encrypt
    Ciphertext<DCRTPoly> ciphertext, ciphertextMult;
    TIC(t);
    for (int i = 0; i < 1000; i++) {
        ciphertext = cryptoContext->Encrypt(keyPair.secretKey, plain);
        ciphertextMult = cryptoContext->Encrypt(keyPair.secretKey, plain);
    }
    processingTime = TOC(t);
    std::cout << "Encrypt time: " << processingTime / 1000 << "ms" << std::endl;

    // decrypt
    Plaintext plain2;
    TIC(t);
    for (int i = 0; i < 1000; i++) {
        cryptoContext->Decrypt(keyPair.secretKey, ciphertext, &plain2);
    }
    processingTime = TOC(t);
    std::cout << "Decrypt time: " << processingTime / 1000 << "ms" << std::endl;

    // mult ctxt * plain
    Ciphertext<DCRTPoly> ciphertext2;
    TIC(t);
    for (int i = 0; i < 1000; i++) {
        ciphertext2 = cryptoContext->EvalMult(ciphertext, plain);
    }
    processingTime = TOC(t);
    std::cout << "Mult ctxt * plain time: " << processingTime / 1000 << "ms" << std::endl;

    // mult ctxt * ctxt
    TIC(t);
    for (int i = 0; i < 1000; i++) {
        ciphertext2 = cryptoContext->EvalMult(ciphertext, ciphertextMult);
    }
    processingTime = TOC(t);
    std::cout << "Mult ctxt * ctxt time: " << processingTime / 1000 << "ms" << std::endl;
    // add
    Ciphertext<DCRTPoly> ciphertext3, ciphertextAdd;
    ciphertextAdd = cryptoContext->Encrypt(keyPair.secretKey, plain);
    TIC(t);
    for (int i = 0; i < 1000; i++) {
        ciphertext3 = cryptoContext->EvalAdd(ciphertext, ciphertextAdd);
    }
    processingTime = TOC(t);
    std::cout << "Add time: " << processingTime / 1000 << "ms" << std::endl;

    memory_usage = getMemoryUsage();

    // Print memory usage in a human-readable format
    memory_usage_gb = static_cast<double>(memory_usage) / (1024.0 * 1024.0 * 1024.0); // Convert to GB

    std::cout << "Memory Usage: " << memory_usage_gb << " GB" << std::endl;
    
}


void Test::kmerTables() {
    // string filename_ref = "../dataset/five_genes/five_gene_reference.fa";
    // string filename_read = "../dataset/five_genes/five_gene_reads.fa";
    string filename_ref = "../dataset/refs_toy.fa";
    string filename_read = "../dataset/reads_toy.fa";

    vector<Sequence> refs_seq;
    vector<Sequence> reads_seq;
    readFastaFile(filename_ref, refs_seq);
    readFastaFile(filename_read, reads_seq);


    KmerTable kmerTableRef, kmerTableRead;
    long K = 7;
    computeKmerTable(refs_seq, K, kmerTableRef);
    computeKmerTableForRead(reads_seq, K, kmerTableRead);

    cout << "geneNameIndex.size() = " << kmerTableRef.geneNameIndex.size() << endl;
    cout << "count.size() = " << kmerTableRef.count.size() << endl;
    cout << "entropy.size() = " << kmerTableRef.entropy.size() << endl;

    cout << "geneNameIndex.size() = " << kmerTableRead.geneNameIndex.size() << endl;
    cout << "count.size() = " << kmerTableRead.count.size() << endl;
    cout << "entropy.size() = " << kmerTableRead.entropy.size() << endl;

    cout << "=== kmerTableRef ===" << endl;
    for (auto &p : kmerTableRef.count) {
        cout << "kmer = " << p.first << " => count = ";
        for (int i = 0; i < p.second.size(); i++) {
            cout << p.second[i] << " ";
        }
        cout << endl;
    }
    for (auto &p : kmerTableRef.entropy) {
        cout << "kmer = " << p.first << " => entropy = " << p.second << endl;
    }

    cout << "=== kmerTableRead ===" << endl;
    for (auto &p : kmerTableRead.countRead) {
        cout << "kmer = " << p.first << " => count = " << p.second << endl;
    }
}

void Test::plainExp() {
    string filename_read = "../dataset/reads_toy.fa";
    string filename_ref = "../dataset/refs_toy.fa";
    long K = 8;
    long B = 4;

    // string filename_read = "../dataset/five_genes/five_gene_reads.fa";
    // string filename_ref = "../dataset/five_genes/five_gene_reference.fa";
    // string filename_read = "../dataset/big/reads_GM12878_20_genes.fa";
    // string filename_ref = "../dataset/big/exon_reference.fa";
    // long K = 32;
    // long B = 8;

    std::srand(std::time(nullptr));

    vector<Sequence> refs_seq;
    readFastaFile(filename_ref, refs_seq);

    // Read read sequences from the file
    vector<Sequence> reads_seq_vec;
    readFastaFile(filename_read, reads_seq_vec);

    // print_sequences(refs_seq, reads_seq_vec);

    unordered_map<int, int> check;
    long total_n_read = 0;
    for (auto &read_seq : reads_seq_vec) {
        total_n_read += read_seq.getNumSeq();
    }
    long current_n = 0;
    while (current_n < 5000) {
        int random_index = std::rand() % reads_seq_vec.size();
        Sequence read_seq = reads_seq_vec[random_index];
        long numSeq = read_seq.getNumSeq();
        string readGeneName = read_seq.getGeneName();
        int random2 = std::rand() % numSeq;
        string read = read_seq.getSeq(random2);
        // for (int i = 0; i < numSeq; i++) {
        vector<long> count_vec;
        // string read = read_seq.getSeq(i);
        cout << "(" << current_n << " / " << total_n_read << ") ";
        current_n++;
        cout << "read from " << readGeneName << " : " << read << endl;
        cout << random2 << "th read from " << random_index << "th gene" << endl;
        long match_sort = 0;
        long max_count = 0;
        long match_count = 0;
        for (auto &refs : refs_seq) {
            long count = count_matching(read, refs, K, B);
            count_vec.push_back(count);
            // cout << "   vs ref gene " << refs.getGeneName();
            if (refs.getGeneName() == readGeneName) {
                // cout << "(v)";
                match_count = count;
            }
            // cout << endl;
            // cout << "   match = " << count << endl;
        }
        sort(count_vec.rbegin(), count_vec.rend());
        max_count = count_vec[0];
        auto it = find(count_vec.begin(), count_vec.end(), match_count);
        match_sort = distance(count_vec.begin(), it);

        cout << "max, match, order = " << max_count << ", " << match_count << ", " << match_sort << endl;
        if (match_count == 0) {
            check[-1] += 1;
        } else {
            check[match_sort] += 1;
        }
        // }
    }
    cout << "check order" << endl;
    for (auto &p : check) {
        cout << "order " << p.first << " : " << p.second << endl;
    }
}

void Test::encryptRead() {
    TimeVar t;
    double processingTime(0.0);
    string filename_read = "../dataset/five_genes/five_gene_reads_shorten.fa";
    // string filename_read = "../dataset/reads_toy.fa";
    long K = 15;

    cout << "read filename = " << filename_read << endl;
    cout << "K = " << K << endl;

    // Read read sequences from the file
    vector<Sequence> reads_seq;
    readFastaFile(filename_read, reads_seq);

    // benchmarking variables
    // TimeVar t;
    // double processingTime(0.0);

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
    std::cout << "log2 q = "
              << log2(cryptoContext->GetCryptoParameters()->GetElementParams()->GetModulus().ConvertToDouble())
              << std::endl;

    // Initialize Public Key Containers
    KeyPair<DCRTPoly> keyPair = cryptoContext->KeyGen();
    cryptoContext->EvalMultKeysGen(keyPair.secretKey);


    KmerTable kmerTableRead;
    computeKmerTableForRead(reads_seq, K, kmerTableRead);


    vector<Ciphertext<DCRTPoly>> ct;
    TIC(t);
    encryptReadKmer(kmerTableRead, K, ct, cryptoContext, keyPair);
    processingTime = TOC(t);
    std::cout << "Encrypt Reads time: " << processingTime << "ms" << std::endl;
}