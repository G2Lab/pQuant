#include "Functions.h"

void print_sequences(vector<Sequence> &refs_seq, vector<Sequence> &reads_seq) {
    cout << "== References ==" << endl;
    for (auto &ref : refs_seq) {
        cout << "Gene Name: " << ref.getGeneName() << endl;
        for (long i = 0; i < ref.getNumSeq(); i++) {
            cout << "Exon " << i + 1 << ": " << ref.getSeq(i) << endl;
        }
        cout << endl;
    }

    cout << "== Reads ==" << endl;
    for (auto &read : reads_seq) {
        cout << "Gene Name: " << read.getGeneName() << endl;
        for (size_t i = 0; i < static_cast<size_t>(read.getNumSeq()); i++) {
            cout << "Exon " << i + 1 << ": " << read.getSeq(i) << endl;
        }
        cout << endl;
    }
}
void encode_seq(CryptoContext<DCRTPoly> &cryptoContext, Sequence &seq, Plaintext_2d &plaintexts, long K, long B, long slot_count) {
    // pre-compute total number of substrings
    // (then compute pack_num)
    long total_num_substr = 0;
    for (size_t i = 0; i < static_cast<size_t>(seq.getNumSeq()); i++) {
        total_num_substr += seq.getSeq(i).size() - K + 1;
    }

    long pack_num = (total_num_substr - 1) / slot_count + 1;

    // vector exon_encoded is 3-dim
    // pack_num * block_num (=K/B) * slot_count
    vector<vector<vector<int64_t>>> exon_encoded;
    exon_encoded.resize(
        pack_num,
        vector<vector<int64_t>>(K / B, vector<int64_t>(slot_count, 0)));

    // encode substrings
    long pack = 0;
    long slot = 0;
    for (size_t i = 0; i < static_cast<size_t>(seq.getNumSeq()); i++) {
        const auto &ref = seq.getSeq(i);
        long N = ref.size();
        for (int j = 0; j < N - K + 1; j++) {
            string subref = ref.substr(j, K);
            for (int k = 0; k < K / B; k++) {
                string sub_subref = subref.substr(k * B, B);
                // remark: encode with negation, then addtion will be comparison
                exon_encoded[pack][k][slot] = -encode_nt_to_num(sub_subref);
            }
            slot += 1;
            if (slot == slot_count) {
                slot = 0;
                pack += 1;
            }
        }
    }

    plaintexts.resize(exon_encoded.size());
    for (size_t i = 0; i < exon_encoded.size(); i++) {
        plaintexts[i].resize(K / B);
        for (int j = 0; j < K / B; j++) {
            plaintexts[i][j] = cryptoContext->MakePackedPlaintext(exon_encoded[i][j]);
        }
    }
}

void encode_refs(CryptoContext<DCRTPoly> &cryptoContext, vector<Sequence> &refs_seq, Plaintext_3d &ref_plain, long K, long B, long slot_count) {
    for (auto &seq : refs_seq) {
        Plaintext_2d exon_plain;
        encode_seq(cryptoContext, seq, exon_plain, K, B, slot_count);
        ref_plain.push_back(exon_plain);
    }

    bool debug = true;
    if (debug) {
        cout << endl;
        cout << ref_plain[0].size() << " * " << ref_plain[0][0].size() << endl;
        for (size_t e = 0; e < ref_plain.size(); e++) {
            cout << " ==" << e << "th exon ==" << endl;
            for (size_t i = 0; i < ref_plain[e].size(); i++) {
                cout << i << "th block" << endl;
                for (size_t j = 0; j < ref_plain[e][i].size(); j++) {
                    Plaintext dummy = ref_plain[e][i][j];
                    dummy->SetLength(16);
                    cout << dummy << endl;
                }
            }
        }
    }
}

void encode_read(CryptoContext<DCRTPoly> &cryptoContext, vector<Sequence> &reads_seq, Plaintext_4d &read_plain, long K, long B, long slot_count) {
    // read_plain is
    // -> 4-dim : exon_num * #reads * #kmers * blocks(K/B)
    // #kmers = read size / K (generally)

    for (auto &seq : reads_seq) {

        // define 3-dim vector for the current sequence
        Plaintext_3d read_encoded;
        // for each exon, define block vectors
        for (size_t i = 0; i < static_cast<size_t>(seq.getNumSeq()); i++) {
            const auto &read = seq.getSeq(i);

            // # of kmers = size / K
            // each kmer is divided into B-size blocks
            // total # of blocks = (size / K) * (K / B) = size / B
            long kmer_num = read.size() / K;
            Plaintext_2d encoded(kmer_num, Plaintext_1d(K / B));

            // cout << "Read = " << read << endl;
            for (int j = 0; j < kmer_num; j++) {
                string kmer = read.substr(K * j, K);
                for (int k = 0; k < K / B; k++) {
                    string kmer_block = kmer.substr(k * B, B);
                    long kmer_block_num = encode_nt_to_num(kmer_block);
                    // cout << "block = " << kmer_block << " -> " << kmer_block_num << endl;
                    vector<int64_t> dummy(slot_count, kmer_block_num);
                    encoded[j][k] = cryptoContext->MakePackedPlaintext(dummy);
                }
            }
            read_encoded.push_back(encoded);
        }
        read_plain.push_back(read_encoded);
    }

    bool debug = true;
    if (debug) {
        cout << endl;
        cout << "K = " << K << ", B = " << B << endl;
        cout << read_plain[0].size() << " * " << read_plain[0][0].size() << read_plain[0][0][0].size() << endl;
        for (size_t e = 0; e < read_plain.size(); e++) {
            cout << " ==" << e << "th exon ==" << endl;
            for (size_t i = 0; i < read_plain[e].size(); i++) {
                cout << i << "th read" << endl;
                for (size_t j = 0; j < read_plain[e][i].size(); j++) {
                    cout << j << "th block / " << read_plain[e][i][j].size() << endl;
                    for (size_t k = 0; k < read_plain[e][i][j].size(); k++) {
                        Plaintext dummy = read_plain[e][i][j][k];
                        dummy->SetLength(16);
                        cout << dummy << endl;
                    }
                }
            }
        }
    }
}

void encode_read2(CryptoContext<DCRTPoly> &cryptoContext, vector<Sequence> &reads_seq, Plaintext_4d &read_plain, long K, long B, long slot_count) {
    // read_plain is
    // -> 4-dim : exon_num * #reads * #kmers * blocks(K/B)
    // #kmers = read size / K (generally)

    for (auto &seq : reads_seq) {

        // define 3-dim vector for the current sequence
        Plaintext_3d read_encoded;
        // for each exon, define block vectors
        vector<int64_t> dummy(slot_count);
        long encoding_loc = 0;
        for (size_t i = 0; i < static_cast<size_t>(seq.getNumSeq()); i++) {
            const auto &read = seq.getSeq(i);

            // # of kmers = size / K
            // each kmer is divided into B-size blocks
            // total # of blocks = (size / K) * (K / B) = size / B
            long kmer_num = read.size() / K;
            Plaintext_2d encoded(kmer_num, Plaintext_1d(K / B));

            cout << "Read = " << read << endl;
            for (int j = 0; j < kmer_num; j++) {
                string kmer = read.substr(K * j, K);
                for (int k = 0; k < K / B; k++) {
                    string kmer_block = kmer.substr(k * B, B);
                    long kmer_block_num = encode_nt_to_num(kmer_block);
                    cout << "block = " << kmer_block << " -> " << kmer_block_num << endl;
                    dummy[encoding_loc] = kmer_block_num;
                    encoding_loc += 1;
                    if (encoding_loc == slot_count) {
                        Plaintext plain = cryptoContext->MakePackedPlaintext(dummy);
                        encoded[j].push_back(plain);
                        dummy.resize(slot_count);
                        encoding_loc = 0;
                    }
                    // vector<int64_t> dummy(slot_count, kmer_block_num);
                    // encoded[j][k] = cryptoContext->MakePackedPlaintext(dummy);
                }
            }
            read_encoded.push_back(encoded);
        }
        read_plain.push_back(read_encoded);
    }

    bool debug = true;
    if (debug) {
        cout << endl;
        cout << "K = " << K << ", B = " << B << endl;
        cout << read_plain[0].size() << " * " << read_plain[0][0].size() << read_plain[0][0][0].size() << endl;
        for (size_t e = 0; e < read_plain.size(); e++) {
            cout << " ==" << e << "th exon ==" << endl;
            for (size_t i = 0; i < read_plain[e].size(); i++) {
                cout << i << "th read" << endl;
                for (size_t j = 0; j < read_plain[e][i].size(); j++) {
                    cout << j << "th block / " << read_plain[e][i][j].size() << endl;
                    for (size_t k = 0; k < read_plain[e][i][j].size(); k++) {
                        Plaintext dummy = read_plain[e][i][j][k];
                        dummy->SetLength(16);
                        cout << dummy << endl;
                    }
                }
            }
        }
    }
}

void encrypt_reads(CryptoContext<DCRTPoly> &cryptoContext, KeyPair<DCRTPoly> &keyPair, Plaintext_4d &read_plain,
                   Ciphertext_4d &ctxt_kmers, long K, long B, long slot_count) {
    // kmer_encoded is
    // -> 4-dim : exon_num * #read * #kmer * blocks
    for (auto &kmer_exon : read_plain) {
        Ciphertext_3d ctxt_exon;
        ctxt_exon.resize(kmer_exon.size());
        for (size_t i = 0; i < kmer_exon.size(); i++) {
            ctxt_exon[i].resize(kmer_exon[i].size());
            for (size_t j = 0; j < kmer_exon[i].size(); j++) {
                ctxt_exon[i][j].resize(kmer_exon[i][j].size());
                for (size_t k = 0; k < kmer_exon[i][j].size(); k++) {
                    // cout << "i,j,k = " << i << ", " << j << ", " << k << endl;
                    ctxt_exon[i][j][k] = cryptoContext->Encrypt(keyPair.publicKey, kmer_exon[i][j][k]);
                }
            }
        }
        ctxt_kmers.push_back(ctxt_exon);
    }
    // output kmer has same structure as kmer_encoded
    // -> 4-dim : exon_num * #read * #kmer * blocks
}

void eqaulity_test(CryptoContext<DCRTPoly> &cryptoContext, KeyPair<DCRTPoly> &keyPair, Ciphertext_4d &ctxt_kmers, Plaintext_3d &ref_plain, Ciphertext_3d &ctxt_out, long K, long B, long log_poly_modulus_degree, vector<Sequence> &refs_seq, vector<Sequence> &reads_seq) {
    // ref_plain(Encoding) structure:
    // num_exon_ref * pack_num * block_num (=K/B) (plaintext)

    // ctxt_kmers structure:
    // num_exon_read * #reads * #kmer * block_num (=k/B) (ciphertext)
    // (read_size / B) kmers come from the same read

    // output structure:
    // num_exon_read * #reads * num_exon_ref
    // compute eqtest for all kmers & packings and aggregate all to save in one ciphertext

    ctxt_out.resize(ctxt_kmers.size());
    for (size_t n_exon = 0; n_exon < ctxt_kmers.size(); n_exon++) {
        ctxt_out[n_exon].resize(ctxt_kmers[n_exon].size());

        for (size_t n_read = 0; n_read < ctxt_kmers[n_exon].size(); n_read++) {
            // cout << "read = " << reads_seq[n_exon].getSeq(n_read) << endl;
            ctxt_out[n_exon][n_read].resize(ref_plain.size());

            for (size_t n_exon_ref = 0; n_exon_ref < ref_plain.size(); n_exon_ref++) {
                Ciphertext_1d ctxt_pack_list;

                for (size_t n_pack = 0; n_pack < ref_plain[n_exon_ref].size(); n_pack++) {
                    // cout << "ref = " << refs_seq[n_exon_ref].getSeq(n_pack) << endl;
                    Ciphertext_1d ctxt_kmers_list;
                    for (size_t n_kmer = 0; n_kmer < ctxt_kmers[n_exon][n_read].size(); n_kmer++) {
                        Ciphertext_1d ctxt_blocks = ctxt_kmers[n_exon][n_read][n_kmer];
                        Ciphertext_1d ctxt_sub(ctxt_blocks.size());
                        vector<Plaintext> plain_blocks = ref_plain[n_exon_ref][n_pack];
                        // ctxt_blocks: K/B-size

                        // substitution
                        long n_blocks = K / B;
                        // cout << "substitute test (" << n_pack << ", " << n_kmer << ")" << endl;
                        for (int i = 0; i < n_blocks; i++) {
                            ctxt_sub[i] = cryptoContext->EvalAdd(ctxt_blocks[i], plain_blocks[i]);

                            // Plaintext ctxt_dec;
                            // cryptoContext->Decrypt(keyPair.secretKey, ctxt_blocks[i], &ctxt_dec);
                            // Plaintext plain_dec;
                            // cryptoContext->Decrypt(keyPair.secretKey, ctxt_sub[i], &plain_dec);
                            // ctxt_dec->SetLength(16);
                            // plain_dec->SetLength(16);
                            // plain_blocks[i]->SetLength(16);
                            // cout << ctxt_dec << endl;
                            // cout << " + " << endl;
                            // cout << plain_blocks[i] << endl;
                            // cout << " = " << endl;
                            // cout << plain_dec << endl;
                            // cout << endl;
                        }

                        // eqtest: sub-square-sum
                        // ctxt_eqtest = (ctxt0 - "0")^2 + (ctxt0 - ctxt1)^2 + ... + (ctxtn-1 - ctxtn)^2 +
                        Ciphertext<DCRTPoly> ctxt_eqtest = cryptoContext->EvalSquare(ctxt_sub[0]);

                        for (int i = 0; i < K / B - 1; i++) {
                            Ciphertext<DCRTPoly> ctxt_dummy = cryptoContext->EvalSub(ctxt_sub[i], ctxt_sub[i + 1]);
                            cryptoContext->EvalSquareInPlace(ctxt_dummy);
                            cryptoContext->EvalAddInPlace(ctxt_eqtest, ctxt_dummy);
                        }

                        // // // instead of square sum, just compute direct sums
                        // for (auto& ctxt : ctxt_sub) {
                        //     cryptoContext->EvalSquareInPlace(ctxt);
                        // }
                        // // Ciphertext<DCRTPoly> ctxt_eqtest = cryptoContext->EvalAddMany(ctxt_sub);
                        // Ciphertext<DCRTPoly> ctxt_eqtest = ctxt_sub[0];
                        // for (auto i = 1; i < ctxt_sub.size(); i++) {
                        //     cryptoContext->EvalAddInPlace(ctxt_eqtest, ctxt_sub[i]);
                        // }

                        // cout << "after square sum computation.." << endl;
                        // Plaintext plain_dec;
                        // cryptoContext->Decrypt(keyPair.secretKey, ctxt_eqtest, &plain_dec);
                        // plain_dec->SetLength(16);
                        // cout << plain_dec << endl;

                        // ctxt <- 1 - ctxt^{p-1}
                        for (size_t i = 0; i < 16; i++) {
                            // cout << i << "th squaring.." << endl;
                            cryptoContext->EvalSquareInPlace(ctxt_eqtest);
                            // Plaintext plain_dec;
                            // cryptoContext->Decrypt(keyPair.secretKey, ctxt_eqtest, &plain_dec);
                            // plain_dec->SetLength(16);
                            // cout << plain_dec << endl;
                        }
                        vector<int64_t> one((1 << log_poly_modulus_degree), -1);
                        Plaintext plain_one = cryptoContext->MakePackedPlaintext(one);
                        cryptoContext->EvalAddInPlace(ctxt_eqtest, plain_one);
                        cryptoContext->EvalNegateInPlace(ctxt_eqtest);

                        // cout << "1 - ^p-1 test" << endl;
                        // cryptoContext->Decrypt(keyPair.secretKey, ctxt_eqtest, &plain_dec);
                        // plain_dec->SetLength(32);
                        // cout << plain_dec << endl;
                        // rotsum
                        for (long i = 0; i < log_poly_modulus_degree; i++) {
                            // cout << "rotate " << i << endl;
                            Ciphertext<DCRTPoly> rot_dummy = cryptoContext->EvalRotate(ctxt_eqtest, (1 << i));
                            cryptoContext->EvalAddInPlace(ctxt_eqtest, rot_dummy);
                            // Plaintext plain_dec;
                            // cryptoContext->Decrypt(keyPair.secretKey, ctxt_eqtest, &plain_dec);
                            // plain_dec->SetLength(16);
                            // cout << plain_dec << endl;
                        }
                        ctxt_kmers_list.push_back(ctxt_eqtest);
                    }

                    Ciphertext<DCRTPoly> ctxt_kmer_aggregate = cryptoContext->EvalAddMany(ctxt_kmers_list);
                    // cout << "check kmers_list" << endl;
                    // for (auto i = 0; i < ctxt_kmers_list.size(); i++) {
                    //     Plaintext plain_kmers;
                    //     cryptoContext->Decrypt(keyPair.secretKey, ctxt_kmers_list[i], &plain_kmers);
                    //     plain_kmers->SetLength(16);
                    //     cout << i << " : " << plain_kmers << endl;
                    // }
                    // cout << "sum" << endl;
                    // Plaintext plain_agg;
                    // cryptoContext->Decrypt(keyPair.secretKey, ctxt_kmer_aggregate, &plain_agg);
                    // plain_agg->SetLength(16);
                    // cout << plain_agg << endl;
                    ctxt_kmers_list.clear();
                    ctxt_pack_list.push_back(ctxt_kmer_aggregate);
                }
                ctxt_out[n_exon][n_read][n_exon_ref] = cryptoContext->EvalAddMany(ctxt_pack_list);
                ctxt_pack_list.clear();
                // cout << "push back ref" << n_exon << ", " << n_read << ", " << endl;
                // cout << "final test in ctxt_out : " << n_exon << ", " << n_read << ", " << n_exon_ref << endl;
                // cout << "checkc aggregate" << endl;
                // for (size_t i = 0; i < ctxt_pack_list.size(); i++) {
                //     Plaintext plain_pack;
                //     cryptoContext->Decrypt(keyPair.secretKey, ctxt_pack_list[i], &plain_pack);
                //     plain_pack->SetLength(16);
                //     cout << i << " : " << plain_pack << endl;
                // }
                // Plaintext plain_test;
                // cryptoContext->Decrypt(keyPair.secretKey, ctxt_out[n_exon][n_read][n_exon_ref], &plain_test);
                // plain_test->SetLength(16);
                // cout << plain_test << endl;
            }
        }
    }
}

void decrypt_output(CryptoContext<DCRTPoly> &cryptoContext, KeyPair<DCRTPoly> &keyPair, Ciphertext_3d &ctxt_out, vector<Sequence> &refs_seq, vector<Sequence> &reads_seq, long K, long B) {
    cout << "ctxt_out size check" << endl;
    cout << "exon_reads = " << ctxt_out.size() << ", reads " << ctxt_out[0].size() << ", exon_num_refs = " << ctxt_out[0][0].size() << endl;

    for (size_t n_exon_read = 0; n_exon_read < ctxt_out.size(); n_exon_read++) {
        for (size_t n_read = 0; n_read < ctxt_out[n_exon_read].size(); n_read++) {
            cout << "read gene " << reads_seq[n_exon_read].getGeneName() << endl;
            cout << " (" << reads_seq[n_exon_read].getSeq(n_read) << ")" << endl;
            for (size_t n_exon_ref = 0; n_exon_ref < ctxt_out[n_exon_read][n_read].size();
                 n_exon_ref++) {
                cout << "   vs ref gene " << refs_seq[n_exon_ref].getGeneName();
                if (refs_seq[n_exon_ref].getGeneName() == reads_seq[n_exon_read].getGeneName()) cout << "(v)";
                cout << endl;
                // for (size_t i = 0; i < refs_seq[n_exon_ref].getNumSeq(); i++) {
                //     cout << "   " << refs_seq[n_exon_ref].getSeq(i) << endl;
                // }
                Plaintext plain_out;
                vector<int64_t> out;
                cryptoContext->Decrypt(keyPair.secretKey, ctxt_out[n_exon_read][n_read][n_exon_ref], &plain_out);
                plain_out->SetLength(1);
                string read = reads_seq[n_exon_read].getSeq(n_read);
                long match_real = count_matching(read, refs_seq[n_exon_ref], K, B);
                cout << "   match(real) = " << plain_out << " (" << match_real << ")" << endl;
            }
        }
    }
}
