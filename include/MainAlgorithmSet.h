#ifndef MAINALGORITHMSET_H
#define MAINALGORITHMSET_H

#include <iostream>
#include "util/FileReader.h"
#include "util/Define.h"
#include "util/KmerTable.h"
#include "func/Functions.h"
#include "func/Entropy.h"
#include "func/PQuantParams.h"
#include "openfhe.h"
#include <filesystem>

class MainAlgorithmSet {
public:
    /*
    step 1
     - compute KmerTable from reference(Gene) dataset
     - save kmerTable to file (will be used in computation)
     - save kmerList to file (will be used in encryption)
    */
    static void generateKmerTableFromReference(PQuantParams &param);

    /*
    step 2
     - generate BFV keys and serialize
    */
    static void keyGenBFVandSerialize(PQuantParams &param);
    /*
    step 3
     - load kmerList
     - encode & encrypt read
    */
    static void encodeAndEncrypt(PQuantParams &param);
    
    /*
    step 4
     - load kmerTable (batch)
     - compute inner product b/w ctxts and table
    */
    static void computeInnerProductBatch(PQuantParams &param);

    /*
    step 5
     - load final ctxt
     - decrypt and return final gene vector
    */
    static void decryptAndReturnGeneVector(PQuantParams &param);
};

#endif // MAINALGORITHMSET_H