#ifndef FileReader_h_
#define FileReader_h_

#include <iomanip>
#include <fstream>
#include <iostream>
#include <string>
#include <unordered_map>
#include <map>
#include <vector>
#include "func/Sequence.h"
#include "func/PQuantParams.h"
#include "util/Define.h"
#include "util/util.h"
#include "util/KmerTable.h"

using namespace std;

void readFastaFile(const string &filename, vector<Sequence>& refs_seq);

void readFastQFile(const string& filename, vector<Sequence>& seq_vec);

// void saveKmerTable(const string &filename, const KmerTable& kmerTable);

// void parseJson(const std::string& jsonFile, KmerTable& kmerTable, PQuantParams& param);

#endif