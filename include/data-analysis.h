#ifndef DATA_ANALYSIS_H
#define DATA_ANALYSIS_H

// Include necessary headers
#include <iostream>
#include <vector>
#include "func/PQuantParams.h"
#include "util/FileReader.h"
#include "util/Define.h"
#include "util/KmerTable.h"
#include "func/Functions.h"
#include "func/PQuantParams.h"
#include <chrono>
#include <cstdlib>
#include <fstream>
#include <sstream>

using namespace std;

void analyze_dataset(PQuantParams param);

void compute_tpm(PQuantParams param);

#endif // DATA_ANALYSIS_H