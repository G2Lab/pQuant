#ifndef TEST_H
#define TEST_H

#include <iostream>
#include "util/FileReader.h"
#include "util/Define.h"
#include "Functions.h"
#include "entropy/Entropy.h"
#include "enc/Functions.h"
#include "openfhe.h"


using namespace std;

size_t getMemoryUsage();

class Test {
public:
    static void testReadFastaFiles(const string &ref_filename, const string &read_filename);

    static void previousAlgorithm();

    static void bfvBasicAlgorithmBenchmarks();

    static void kmerTables();

    static void plainExp();

    static void encryptRead();
};

#endif