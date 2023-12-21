#ifndef ENTROPY_H
#define ENTROPY_H

#include <vector>
#include <string>
#include <map>
#include <iostream>
#include "Sequence.h"
#include "util/Define.h"
#include "util/FileReader.h"
using namespace std;


void computeKmerTable(vector<Sequence>& gene, long K, KmerTable &kmerTable);

void computeKmerTableForRead(vector<Sequence>& read, long K, KmerTable &kmerTable);

#endif