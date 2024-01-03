#include "util/FileReader.h"

#include <iostream>
#include <gtest/gtest.h>

// Example function to be tested
int add(int a, int b) {
    return a + b;
}

// Test case to check the add function
TEST(AddFunctionTest, PositiveNumbers) {
    EXPECT_EQ(add(2, 3), 5);
}

TEST(FileReaderTest, ReadFastaFile) {
    std::string filename_read = "../dataset/test/reads_toy.fa";
    std::string filename_ref = "../dataset/test/refs_toy.fa";
    // Read reference sequences from the file
    std::vector<Sequence> refs_seq;
    readFastaFile(filename_ref, refs_seq);

    // Read read sequences from the file
    vector<Sequence> reads_seq;
    readFastaFile(filename_read, reads_seq);
    // Read reference sequences from the file

    EXPECT_EQ(refs_seq.size(), 3);
    EXPECT_EQ(refs_seq[0].getGeneName(), "ENSG00000085984");
    EXPECT_EQ(refs_seq[1].getGeneName(), "ENSG00000085983");
    EXPECT_EQ(refs_seq[2].getGeneName(), "ENSG00000085982");
    
    EXPECT_EQ(refs_seq[0].getNumSeq(), 2);
    EXPECT_EQ(refs_seq[1].getNumSeq(), 2);
    EXPECT_EQ(refs_seq[2].getNumSeq(), 4);

    long len_max = 0;
    long len_min = 999;    
    for (auto &ref_seq : refs_seq) {
        for (int i = 0; i < ref_seq.getNumSeq(); i++) {
            long len = ref_seq.getSeq(i).size();
            if (len_max < len)
                len_max = len;
            if (len_min > len)
                len_min = len;
        }
    }

    EXPECT_EQ(len_max, 16);
    EXPECT_EQ(len_min, 16);

    long n_gene_read = 0;
    long total_n_read = 0;
    for (auto &read_seq : reads_seq) {
        n_gene_read += 1;
        for (int i = 0; i < read_seq.getNumSeq(); i++) {
            total_n_read += 1;
        }
    }

    EXPECT_EQ(n_gene_read, 1);
    EXPECT_EQ(total_n_read, 4);
}

int main(int argc, char **argv) {
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}