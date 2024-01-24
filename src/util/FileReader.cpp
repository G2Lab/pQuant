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
            auto& seq = gene_to_seq[0];
            seq.setGeneName(0);
            seq.addSeq(line);

        }
    }

    file.close();
}

using json = nlohmann::json;

void parseJson(const std::string& jsonFile, KmerTable& kmerTable) {
    std::ifstream file(jsonFile);
    if (!file.is_open()) {
        std::cerr << "Error opening file: " << jsonFile << std::endl;
        return;
    }

    json jsonData;
    try {
        file >> jsonData;
    } catch (const std::exception& e) {
        std::cerr << "Error parsing JSON: " << e.what() << std::endl;
        return;
    }
    
    // Extract gene_name_vec
    if (jsonData.contains("gene_name_vec")) {
        kmerTable.geneNameIndex = jsonData["gene_name_vec"].get<std::vector<std::string>>();
    }

    kmerTable.n_gene = kmerTable.geneNameIndex.size();
    // for (size_t i = 0; i < kmerTable.n_gene; i++) {
    //     kmerTable.count[i] = {};
    // }
    // Extract kmer_matrix
    if (jsonData.contains("kmer_matrix")) {
        for (json::iterator it = jsonData["kmer_matrix"].begin(); it != jsonData["kmer_matrix"].end(); ++it) {
            for (json::iterator it2 = it.value().begin(); it2 != it.value().end(); ++it2) {
                kmerTable.count[stoi(it2.key())][stoi(it.key())] = it2.value();
            }
        }
    }

    // Extract kmer_entropy
    kmerTable.n_kmer_total = 0;
    if (jsonData.contains("kmer_entropy")) {
        for (json::iterator it = jsonData["kmer_entropy"].begin(); it != jsonData["kmer_entropy"].end(); ++it) {
            kmerTable.entropy.insert(make_pair(stoi(it.key()), it.value()));
            kmerTable.n_kmer_total += 1;
        }
    }
}