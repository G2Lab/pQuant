#include "util/FileReader.h"

void readFastaFile(const string& filename, vector<Sequence>& seq_vec) {
    ifstream file(filename);
    if (!file.is_open()) {
        throw runtime_error("Error: Failed to open file " + filename);
    }

    string line;
    string gene_name;
    // long exon_num = -1;
    vector<string> exon_seq;

    // Keep track of Sequence objects for each gene name
    unordered_map<string, Sequence> gene_to_seq;
    while (getline(file, line)) {
        if (line.empty()) {
            continue;
        }

        if (line[0] == '>') {
            // Add the previous sequence to the appropriate Sequence object
            if (!exon_seq.empty()) {
                auto& seq = gene_to_seq[gene_name];
                if (seq.getGeneName() == "NA") {
                    seq.setGeneName(gene_name);
                }
                for (size_t i = 0; i < exon_seq.size(); i++) {
                    seq.addSeq(exon_seq[i]);
                }
            }

            // Parse gene name and exon number from the header line
            size_t gene_pos = 1;
            size_t exon_pos = line.find("exon");
            size_t chr_pos = line.find("_");
            if (exon_pos == string::npos || chr_pos == string::npos) {
                throw runtime_error("Error: Invalid header line in file " + filename);
            }
            gene_name = line.substr(gene_pos, exon_pos - 1 - gene_pos);
            
            // Reset exon_seq for the new sequence
            exon_seq.clear();
        } else {
            // Concatenate the current line to the current exon sequence
            if (exon_seq.empty()) {
                exon_seq.push_back(line);
            } else {
                exon_seq.back() += line;
            }
        }
    }

    // Add the final sequence to the appropriate Sequence object
    if (!exon_seq.empty()) {
        auto& seq = gene_to_seq[gene_name];
        if (seq.getGeneName() == "NA") {
            seq.setGeneName(gene_name);
        }
        for (size_t i = 0; i < exon_seq.size(); i++) {
            seq.addSeq(exon_seq[i]);
        }
    }

    // Add all Sequence objects to the seq_vec vector
    for (auto& entry : gene_to_seq) {
        seq_vec.push_back(entry.second);
    }

    file.close();
}

void readFastQFile(const string& filename, vector<Sequence>& seq_vec) {
    ifstream file(filename);
    if (!file.is_open()) {
        throw runtime_error("Error: Failed to open file " + filename);
    }

    string line;
    int line_number = 0;
    vector<string> exon_seq;

    // Keep track of Sequence objects for each gene name
    unordered_map<string, Sequence> gene_to_seq;

    while (getline(file, line)) {
        if (line.empty()) {
            continue;
        }

        line_number++;

        if (line_number % 4 == 2) {
            auto& seq = gene_to_seq["gene"];
            seq.setGeneName("gene");
            seq.addSeq(line);

        }
    }

    file.close();
}


// using json = nlohmann::json;

// void parseJson(const std::string& jsonFile, KmerTable& kmerTable, PQuantParams& param) {
//     std::ifstream file(jsonFile);
//     if (!file.is_open()) {
//         std::cerr << "Error opening file: " << jsonFile << std::endl;
//         return;
//     }

//     json jsonData;
//     auto start_time = std::chrono::high_resolution_clock::now();
//     try {
//         cout << "read json file" << endl;
//         file >> jsonData;
//     } catch (const std::exception& e) {
//         std::cerr << "Error parsing JSON: " << e.what() << std::endl;
//         return;
//     }
//     cout << "json file read" << endl;
//     auto end_time = std::chrono::high_resolution_clock::now();
//     auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time).count();
//     cout << "Time to read json file: " << duration << " ms" << endl;
    
//     // Extract gene_name_vec
//     if (jsonData.contains("gene_name_vec")) {
//         kmerTable.geneNameIndex = jsonData["gene_name_vec"].get<std::vector<std::string>>();
//     }

//     kmerTable.n_gene = kmerTable.geneNameIndex.size();

//     // Extract kmer_matrix
//     cout << " extracting kmer matrix table" << endl;
//     start_time = std::chrono::high_resolution_clock::now();
//     if (jsonData.contains("kmer_matrix")) {
//         for (json::iterator it = jsonData["kmer_matrix"].begin(); it != jsonData["kmer_matrix"].end(); ++it) {
//             for (json::iterator it2 = it.value().begin(); it2 != it.value().end(); ++it2) {
//                 stringstream s2(it2.key());
//                 stringstream s1(it.key());
//                 size_t it2_key, it_key;
//                 s2 >> it2_key;
//                 s1 >> it_key;
//                 kmerTable.count[it2_key][it_key] = it2.value();
//             }
//         }
//     }
//     end_time = std::chrono::high_resolution_clock::now();
//     duration = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time).count();
//     cout << "Time to extract kmer matrix table: " << duration << " ms" << endl;

//     // Extract kmer_entropy
//     kmerTable.n_kmer_total = 0;
//     start_time = std::chrono::high_resolution_clock::now();
//     if (jsonData.contains("kmer_entropy")) {
//         int i = 0;
//         for (json::iterator it = jsonData["kmer_entropy"].begin(); it != jsonData["kmer_entropy"].end(); ++it) {
//             stringstream s1(it.key());
//             size_t it_key;
//             s1 >> it_key;
//             kmerTable.entropy.insert(make_pair(it_key, it.value()));
//             kmerTable.n_kmer_total += 1;
//             if (param.progress_bar)
//                 print_progress_bar("load kmer entropy table", i, jsonData["kmer_entropy"].size(), start_time);
//             i += 1;
//         }
//     }
//     end_time = std::chrono::high_resolution_clock::now();
//     duration = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time).count();
//     cout << "Time to extract kmer entropy table: " << duration << " ms" << endl;
// }