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

using namespace std;

void readFastaFile(const string &filename, vector<Sequence>& refs_seq);

// typedef struct KmerTable {
//     long K;
//     vector<string> geneNameIndex;
//     map<long, vector<long>> count;
//     map<long, long> countRead;
//     map<long, float> entropy;
// } KmerTable;


// void parseJson(const std::string& jsonStr, KmerTable& kmerTable) {
//     // Find kmer_matrix
//     size_t kmerMatrixStart = jsonStr.find("\"kmer_matrix\":{");
//     size_t kmerMatrixEnd = jsonStr.find("}", kmerMatrixStart);
//     std::string kmerMatrixStr = jsonStr.substr(kmerMatrixStart + 14, kmerMatrixEnd - kmerMatrixStart - 14);

//     std::istringstream kmerMatrixStream(kmerMatrixStr);
//     int key, innerKey, innerValue;
//     char comma;
//     while (kmerMatrixStream >> key >> comma >> innerKey >> comma >> innerValue) {
//         kmerTable.kmer_matrix[key][innerKey] = innerValue;
//         kmerMatrixStream >> comma; // consume the comma
//     }

//     // Find kmer_entropy
//     size_t kmerEntropyStart = jsonStr.find("\"kmer_entropy\":{");
//     size_t kmerEntropyEnd = jsonStr.find("}", kmerEntropyStart);
//     std::string kmerEntropyStr = jsonStr.substr(kmerEntropyStart + 16, kmerEntropyEnd - kmerEntropyStart - 16);

//     std::istringstream kmerEntropyStream(kmerEntropyStr);
//     while (kmerEntropyStream >> key >> comma >> innerValue) {
//         kmerTable.kmer_entropy[key] = innerValue;
//         if (kmerEntropyStream.peek() == ',') {
//             kmerEntropyStream.ignore(); // consume the comma
//         }
//     }

//     // Find gene_name_vec
//     size_t geneNameVecStart = jsonStr.find("\"gene_name_vec\":[");
//     size_t geneNameVecEnd = jsonStr.find("]", geneNameVecStart);
//     std::string geneNameVecStr = jsonStr.substr(geneNameVecStart + 16, geneNameVecEnd - geneNameVecStart - 16);

//     std::istringstream geneNameVecStream(geneNameVecStr);
//     std::string geneName;
//     while (std::getline(geneNameVecStream, geneName, ',')) {
//         // geneName.erase(std::remove_if(geneName.begin(), geneName.end(), ::isspace), geneName.end());
//         kmerTable.gene_name_vec.push_back(geneName);
//     }
// }

#endif