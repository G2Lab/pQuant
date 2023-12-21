#ifndef PlainFunc_h_
#define PlainFunc_h_

#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include "util/FileReader.h"
#include "util/util.h"
#include "Sequence.h"
using namespace std;

long count_matching(string &read, Sequence &ref_seq, long K, long B);


#endif