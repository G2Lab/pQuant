#ifndef TASK_H
#define TASK_H

#include <iostream>
#include "util/FileReader.h"
#include "util/Define.h"
#include "util/KmerTable.h"
#include "func/Functions.h"
#include "func/Entropy.h"
#include "func/PQuantParams.h"
#include "openfhe.h"


using namespace std;


class Task {
public:
    static void readFastaFiles(PQuantParams &param);

    static void testKmerTable(PQuantParams &param);
    
    static void bfvBenchmark(PQuantParams &param);

    static void run_all(PQuantParams &param);

    // static void testReadJson(PQuantParams &param);
};

#endif