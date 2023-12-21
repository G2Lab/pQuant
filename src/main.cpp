#include "Functions.h"
#include "Sequence.h"
#include "PlainFunc.h"
#include <algorithm>
#include <cstdlib>
#include <ctime>
#include <iostream>
#include "entropy/Entropy.h"
#include "openfhe.h"
#include "test/Test.h"


using namespace std;
using namespace lbcrypto;


int main(int argc, char **argv) {
    if (argc < 2) {
        cout << "No arguments found" << endl;
        return 0;
    }
    // string filename_read = "../dataset/big/reads_GM12878_20_genes.fa";
    // string filename_ref = "../dataset/big/exon_reference.fa";
    // string filename_read = "../dataset/five_genes/five_gene_reads_shorten.fa";
    if (strcmp(argv[1], "test") == 0) {
        if (strcmp(argv[2], "fasta") == 0) {
            string filename_read = "../dataset/five_genes/five_gene_reads.fa";
            string filename_ref = "../dataset/five_genes/five_gene_reference.fa";
            Test::testReadFastaFiles(filename_ref,filename_read);
        } else if (strcmp(argv[2], "bfv") == 0) {
            Test::previousAlgorithm();
        } else if (strcmp(argv[2], "bfvbench") == 0) {
            Test::bfvBasicAlgorithmBenchmarks();
        } else if (strcmp(argv[2], "kmer") == 0) {
            Test::kmerTables();
        } else if (strcmp(argv[2], "plain") == 0) {
            Test::plainExp();
        } else if (strcmp(argv[2], "encread") == 0) {
            Test::encryptRead();
        }
    } else {
        // printout error: no matching arguments found
        cout << "No matching arguments found" << endl;
    }
    return 0;
}
