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
#include "util/Define.h"
#include "json.hpp"

using namespace std;

void readFastaFile(const string &filename, vector<Sequence>& refs_seq);


void parseJson(const std::string& jsonStr, KmerTable& kmerTable);

#endif