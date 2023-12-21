#include "PlainFunc.h"

long count_matching(string &read, Sequence &ref_seq, long K, long B) {
    long count = 0;
    long kmer_num = read.size() / K;
    for (long j = 0; j < kmer_num; j++) {
        string kmer = read.substr(K * j, K);
        for (int i = 0; i < ref_seq.getNumSeq(); i++) {
            string ref = ref_seq.getSeq(i);
            long ref_size = ref.size();
            for (int s = 0; s < ref_size - K + 1; s++) {
                if (kmer == ref.substr(s, K)) {
                    count += 1;
                }
            }
        }
    }
    return count;
}