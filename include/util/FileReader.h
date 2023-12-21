#ifndef FileReader_h_
#define FileReader_h_

#include <iomanip>
#include <fstream>
#include <iostream>
#include <string>
#include <unordered_map>
#include <vector>
#include "Sequence.h"

using namespace std;


void readFastaFile(const string &filename, vector<Sequence>& refs_seq);

#endif