#include "util/KmerTable.h"


long convertKmerToNum(string kmer) {
    long num = 0;
    for (size_t i = 0; i < kmer.size(); i++) {
        if (kmer[i] == 'A') {
            num = num * 4 + 0;
        } else if (kmer[i] == 'C') {
            num = num * 4 + 1;
        } else if (kmer[i] == 'G') {
            num = num * 4 + 2;
        } else if (kmer[i] == 'T') {
            num = num * 4 + 3;
        } else {
            return -1;
        }
    }
    return num;
}

string sequenceReverseCompliment(string nt) {
    string ans = "";
    for (size_t i = 0; i < nt.size(); i++) {
        if (nt[i] == 'A') {
            ans = 'T' + ans;
        } else if (nt[i] == 'C') {
            ans = 'G' + ans;
        } else if (nt[i] == 'G') {
            ans = 'C' + ans;
        } else if (nt[i] == 'T') {
            ans = 'A' + ans;
        } else {
            return "";
        }
    }
    return ans;
}

KmerTable::KmerTable(vector<Sequence>& gene, PQuantParams &param, bool isRef_) {
    K = param.k;
    thres = param.thres;
    isRef = isRef_;
    n_gene = gene.size();

    // case 1: from Reference
    if (isRef) {
        if (geneNameIndex.size() == 0) {
            // save gene list into geneNameIndex
            for (size_t i = 0; i < gene.size(); i++) {
                geneNameIndex.push_back(gene[i].getGeneName());
            }
        }

        auto start_time = std::chrono::high_resolution_clock::now();
        map<size_t, map<size_t, size_t>> kmer_occurance;
        // map<size_t, map<size_t, size_t>> kmer_based_table;
        // count kmers in gene
        for (size_t i = 0; i < gene.size(); i++) {
            // gene i
            map<size_t, size_t> gene_kmer_count;
            for (size_t j = 0; j < static_cast<size_t>(gene[i].getNumSeq()); j++) {
                string seq = gene[i].getSeq(j);
                if (seq.size() < static_cast<size_t>(K)) {
                    continue;
                }

                // original order
                for (size_t k = 0; k < seq.size() - K + 1; k++) {
                    string kmer = seq.substr(k, K);
                    long num = convertKmerToNum(kmer);
                    if (num == -1) {
                        continue;
                    }
                    gene_kmer_count[num] += 1;
                    kmer_occurance[num][i] += 1;
                }

                // reverse compliment
                string seq_rc = sequenceReverseCompliment(seq);
                for (size_t k = 0; k < seq_rc.size() - K + 1; k++) {
                    string kmer = seq_rc.substr(k, K);
                    long num = convertKmerToNum(kmer);
                    if (num == -1) {
                        continue;
                    }
                    gene_kmer_count[num] += 1;
                    kmer_occurance[num][i] += 1;
                }
            }

            // concatenations
            for (size_t j = 0; j < static_cast<size_t>(gene[i].getNumSeq()); j++) {
                for (size_t k = j + 1; k < static_cast<size_t>(gene[i].getNumSeq()); k++) {
                    string seq_j = gene[i].getSeq(j);
                    string seq_k = gene[i].getSeq(k);
                    string seq_combined;
                    if (seq_j.size() >= static_cast<size_t>(K) && seq_k.size() >= static_cast<size_t>(K)) {
                        seq_combined = seq_j.substr(seq_j.size() - K + 1, K - 1) + seq_k.substr(0, K - 1);
                    } else if (seq_j.size() >= static_cast<size_t>(K)) {
                        seq_combined = seq_j.substr(seq_j.size() - K + 1, K - 1) + seq_k;
                    } else if (seq_k.size() >= static_cast<size_t>(K)) {
                        seq_combined = seq_j + seq_k.substr(0, K - 1);
                    } else {
                        seq_combined = seq_j + seq_k;
                    }
                    if (seq_combined.size() < static_cast<size_t>(K)) {
                        continue;
                    }
                    // original order
                    for (size_t l = 0; l < seq_combined.size() - K + 1; l++) {
                        string kmer = seq_combined.substr(l, K);
                        long num = convertKmerToNum(kmer);
                        if (num == -1) {
                            continue;
                        }
                        gene_kmer_count[num] += 1;
                        kmer_occurance[num][i] += 1;
                    }

                    // reverse compliement
                    string seq_combined_rc = sequenceReverseCompliment(seq_combined);
                    for (size_t l = 0; l < seq_combined_rc.size() - K + 1; l++) {
                        string kmer = seq_combined_rc.substr(l, K);
                        long num = convertKmerToNum(kmer);
                        if (num == -1) {
                            continue;
                        }
                        gene_kmer_count[num] += 1;
                        kmer_occurance[num][i] += 1;
                    }
                }
            }

            // save gene_kmer_count into count
            count.insert(make_pair(i, gene_kmer_count));
            gene_kmer_count.clear();

            print_progress_bar("countKmersPerGene", i, gene.size(), start_time);
        }

        std::map<size_t, float> entropy_total;
        // compute entropy
        for (auto it = kmer_occurance.begin(); it != kmer_occurance.end(); ++it) {
            size_t kmer = it->first;
            float entropy_ = 0;
            map<size_t, size_t> kmer_vector = kmer_occurance[kmer];
            float sum = 0;
            for (auto it2 = kmer_vector.begin(); it2 != kmer_vector.end(); ++it2) {
                sum += (float) it2->second;
            }
            for (auto it2 = kmer_vector.begin(); it2 != kmer_vector.end(); ++it2) {
                float p = static_cast<float>(it2->second) / sum;
                entropy_ += -p * log2(p);
            }
            entropy_total.insert(make_pair(kmer, entropy_));
        }

        // threshold entropy
        for (auto it = entropy_total.begin(); it != entropy_total.end(); ++it) {
            if (it->second < thres) {
                entropy.insert(make_pair(it->first, it->second));
            }
        }
        entropy_total.clear();

        // threshold kmer_occurance
        for (auto it = kmer_occurance.begin(); it != kmer_occurance.end(); ++it) {
            if (entropy.find(it->first) != entropy.end()) {
                tableRef.insert(make_pair(it->first, it->second));
            }
        }
        n_kmer_total = entropy.size();
        kmer_occurance.clear();
    } else {
        // case 2: from read
        // count kmers in gene
        for (size_t i = 0; i < gene.size(); i++) {
            for (int j = 0; j < gene[i].getNumSeq(); j++) {
                string seq = gene[i].getSeq(j);
                if (seq.size() < static_cast<size_t>(K)) {
                    continue;
                }
                for (size_t k = 0; k < seq.size() - K + 1; k++) {
                    string kmer = seq.substr(k, K);
                    long num = convertKmerToNum(kmer);
                    if (num == -1) {
                        continue;
                    }
                    if (countRead.find(num) == countRead.end()) {
                        countRead[num] = 0;
                    }
                    countRead[num] += 1;
                }
            }
        }
    }
}

bool KmerTable::operator==(const KmerTable& other) const {
    if (K != other.K) {
        return false;
    }
    if (thres != other.thres) {
        return false;
    }
    if (n_gene != other.n_gene) {
        return false;
    }
    if (n_kmer_total != other.n_kmer_total) {
        return false;
    }
    if (geneNameIndex != other.geneNameIndex) {
        return false;
    }
    if (count != other.count) {
        return false;
    }
    if (tableRef != other.tableRef) {
        return false;
    }
    if (countRead != other.countRead) {
        return false;
    }
    if (entropy != other.entropy) {
        return false;
    }
    if (isRef != other.isRef) {
        return false;
    }
    return true;

}

void KmerTable::print() {
    if (isRef) {
        std::cout << "K = " << K << std::endl;
        std::cout << "thres = " << thres << std::endl;
        std::cout << "n_gene = " << n_gene << std::endl;
        std::cout << "n_kmer_total = " << n_kmer_total << std::endl;
        std::cout << std::endl;
        std::cout << "geneNames" << std::endl;
        for (auto &p : geneNameIndex) {
            std::cout << p << ", ";
        }
        std::cout << std::endl
                  << std::endl;

        std::cout << " KmerTable for reference" << std::endl;
        std::cout << "n_gene = " << n_gene << std::endl;
        std::cout << "n_kmer_total = " << n_kmer_total << std::endl;

        for (const auto &geneEntry : count) {
            std::cout << "Gene: " << geneEntry.first << " (" << geneEntry.second.size() << " kmers)" << std::endl;

            for (const auto &kmerEntry : geneEntry.second) {
                std::cout << "(" << kmerEntry.first << ", " << kmerEntry.second << ") ";
            }
            std::cout << std::endl;
        }

        std::cout << "kmer entropy" << std::endl;
        for (auto &p : entropy) {
            std::cout << p.first << ", " << p.second << std::endl;
        }

        std::cout << "kmer table" << std::endl;
        for (auto &p : tableRef) {
            std::cout << p.first << ": ";
            for (auto &q : p.second) {
                std::cout << "(" << q.first << ", " << q.second << ") ";
            }
            std::cout << std::endl;
        }
    } else {
        std::cout << " KmerTable for read" << std::endl;
        std::cout << "countRead.size() = " << countRead.size() << std::endl;
        for (auto &p : countRead) {
            std::cout << "(" << p.first << ", " << p.second << ") ";
        }
    }
    std::cout << std::endl;
}

void KmerTable::save(const std::string& filename) {
    std::ofstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Error opening file: " << filename << std::endl;
        return;
    }
    // save K, thres, n_gene, and n_kmer_total
    file << "K: " << K << std::endl;
    file << "thres: " << thres << std::endl;
    file << "n_gene: " << n_gene << std::endl;
    file << "n_kmer_total: " << n_kmer_total << std::endl;
    file << "isRef = " << isRef << std::endl;
    file << "geneNameIndex: [";
    for (size_t i = 0; i < geneNameIndex.size(); i++) {
        file << geneNameIndex[i];
        if (i < geneNameIndex.size() - 1) {
            file << ", ";
        }
    }
    file << "]" << std::endl;

    file << "count: {" << std::endl;
    auto start_time = std::chrono::high_resolution_clock::now();
    int i = 0;
    for (auto p : count) {
        file << p.first << ": {";
        for (auto q : p.second) {
            file << q.first << ": " << q.second;
            if (q.first != p.second.rbegin()->first) {
                file << ",";
            }
            file << " ";
        }
        file << "}";
        if (p.first != count.rbegin()->first) {
            file << ",";
        }
        file << std::endl;
        print_progress_bar("count saving..", i, count.size(), start_time);
        i += 1;
    }
    file << "}" << std::endl;

    file << "tableRef: {" << std::endl;
    i = 0;
    start_time = std::chrono::high_resolution_clock::now();
    for (auto p : tableRef) {
        file << p.first << ": {";
        for (auto q : p.second) {
            file << q.first << ": " << q.second;
            if (q.first != p.second.rbegin()->first) {
                file << ",";
            }
            file << " ";
        }
        file << "}";
        if (p.first != tableRef.rbegin()->first) {
            file << ",";
        }
        file << std::endl;
        print_progress_bar("tableRef saving..", i, tableRef.size(), start_time);
        i += 1;
    }
    file << "}" << std::endl;

    file << "entropy: {";
    std::cout << entropy.size() << std::endl;
    i = 0;
    start_time = std::chrono::high_resolution_clock::now();
    for (auto p : entropy) {
        file << p.first << ": " << p.second;
        if (p.first != entropy.rbegin()->first) {
            file << ",";
        }
        file << " ";
        print_progress_bar("entropy saving..", i, entropy.size(), start_time);
        i += 1;
    }
    file << "}";

    file.close();
}

void KmerTable::saveKmerList(const std::string& filename) {
    std::ofstream file_kmer(filename);
    if (!file_kmer.is_open()) {
        std::cerr << "Error opening file: " << filename + "_kmer_list.txt" << std::endl;
        return;
    }
    for (auto p: entropy) {
        file_kmer << p.first << std::endl;
    }
    file_kmer.close();
}

void KmerTable::load(const std::string& filename) {
    std::ifstream file(filename);

    if (!file.is_open()) {
        std::cerr << "Error opening file: " << filename << ".txt" << std::endl;
        return;
    }

    // Load K, thres, n_gene, and n_kmer_total
    // Load K from the file
    std::string line;
    std::getline(file, line); // Read the first line containing K data
    std::istringstream iss(line);
    std::string prefix;
    std::getline(iss, prefix, ' '); // Remove "K: "
    std::string temp;
    std::getline(iss, temp); // Read the remaining value
    K = std::stol(temp); // Convert string to long

    // Load thres from the file
    std::getline(file, line); // Read the line containing thres data
    std::istringstream iss_thres(line);
    std::getline(iss_thres, prefix, ' '); // Remove "thres: "
    std::getline(iss_thres, temp); // Read the remaining value
    thres = std::stof(temp); // Convert string to float containing n_kmer_total data

    // Load n_gene from the file
    std::getline(file, line); // Read the line containing n_gene data
    std::istringstream iss_n_gene(line);
    std::getline(iss_n_gene, prefix, ' '); // Remove "n_gene: "
    std::getline(iss_n_gene, temp); // Read the remaining value
    n_gene = std::stol(temp); // Convert string to long

    // Load n_kmer_total from the file
    std::getline(file, line); // Read the line containing n_kmer_total data
    std::istringstream iss_n_kmer_total(line);
    std::getline(iss_n_kmer_total, prefix, ' '); // Remove "n_kmer_total: "
    std::getline(iss_n_kmer_total, temp); // Read the remaining value
    n_kmer_total = std::stol(temp); // Convert string to long

    // Load isRef from the file
    std::getline(file, line); // Read the line containing isRef data
    std::istringstream iss_isRef(line);
    std::getline(iss_isRef, prefix, ' '); // Remove "isRef = "
    std::getline(iss_isRef, temp); // Read the remaining value
    isRef = (temp == "true") ? true : false; // Convert string to boolean


    // Load geneNameIndex
    std::getline(file, line); // Read the next line containing geneNameIndex data
    std::istringstream iss5(line);
    std::getline(iss5, temp, '['); // Remove "geneNameIndex: ["
    while (std::getline(iss5, temp, ',')) {
        if (!temp.empty()) {
            // remove "]" from the last gene name
            if (temp.find("]") != std::string::npos) {
                temp = temp.substr(0, temp.size() - 1);
            }
            geneNameIndex.push_back(temp);
        }
    }

    // Load count
    while (std::getline(file, line)) {
        if (line.find("count: {") != std::string::npos) {
            break;
        }
    }
    count.clear(); // Clear existing data
    while (std::getline(file, line)) {
        // break if line ONLY contains "}"
        if (line == "}") {
            break;
        }
        size_t gene;
        std::map<size_t, size_t> kmerData;
        std::istringstream iss6(line);
        iss6 >> gene;
        char seperator;
        cout << "gene = " << gene << endl;
        while (iss6 >> seperator) {
            if (seperator == ',' || seperator == '{') {
                size_t key, val;
                iss6 >> key >> seperator >> val;
                kmerData[key] = val;
                cout << "key = " << key << ", val = " << val << endl;
            }
        }
        count[gene] = kmerData;
    }

    // Load tableRef
    while (std::getline(file, line)) {
        if (line.find("tableRef: {") != std::string::npos) {
            break;
        }
    }
    tableRef.clear(); // Clear existing data
    while (std::getline(file, line)) {
        // break if line ONLY contains "}"
        if (line == "}") {
            break;
        }
        size_t kmer;
        std::map<size_t, size_t> kmerData;
        std::istringstream iss6(line);
        iss6 >> kmer;
        char seperator;
        while (iss6 >> seperator) {
            if (seperator == ',' || seperator == '{') {
                size_t key, val;
                iss6 >> key >> seperator >> val;
                kmerData[key] = val;
            }
        }
        tableRef[kmer] = kmerData;
    }

    // Load entropy
    while (std::getline(file, line)) {
        if (line.find("entropy:") != std::string::npos) {
            break;
        }
    }
    entropy.clear(); // Clear existing data
    std::istringstream iss6(line);
    char seperator;
    while (iss6 >> seperator) {
        if (seperator == ',' || seperator == '{') {
            size_t key;
            float val;
            iss6 >> key >> seperator >> val;
            entropy[key] = val;
        }
    }

    file.close();
}

// Function to save KmerTable data to a binary file
void KmerTable::save_binary(const std::string &filename) {
    std::ofstream file(filename + ".bin", std::ios::binary);
    std::ofstream file_kmer(filename + "_kmer_list.bin", std::ios::binary);
    if (!file.is_open()) {
        std::cerr << "Error opening file: " << filename << std::endl;
        return;
    }
    if (!file_kmer.is_open()) {
        std::cerr << "Error opening file: " << filename + "_kmer_list.bin" << std::endl;
        return;
    }

    // Write K, thres, n_gene, and n_kmer_total
    file.write(reinterpret_cast<char *>(&K), sizeof(K));
    file.write(reinterpret_cast<char *>(&thres), sizeof(thres));
    file.write(reinterpret_cast<char *>(&n_gene), sizeof(n_gene));
    file.write(reinterpret_cast<char *>(&n_kmer_total), sizeof(n_kmer_total));

    // Write geneNameIndex
    size_t geneNameIndexSize = geneNameIndex.size();
    file.write(reinterpret_cast<char *>(&geneNameIndexSize), sizeof(geneNameIndexSize));
    for (const auto &geneName : geneNameIndex) {
        file.write(geneName.c_str(), geneName.size() + 1); // Include null terminator
    }

    // Write count
    size_t countSize = count.size();
    file.write(reinterpret_cast<char *>(&countSize), sizeof(countSize));
    for (const auto &geneEntry : count) {
        size_t gene = geneEntry.first;
        file.write(reinterpret_cast<char *>(&gene), sizeof(gene));
        size_t innerMapSize = geneEntry.second.size();
        file.write(reinterpret_cast<char *>(&innerMapSize), sizeof(innerMapSize));
        for (const auto &kmerEntry : geneEntry.second) {
            size_t kmer = kmerEntry.first;
            size_t num = kmerEntry.second;
            file.write(reinterpret_cast<char *>(&kmer), sizeof(kmer));
            file.write(reinterpret_cast<char *>(&num), sizeof(num));
        }
    }

    // Write countRead
    size_t countReadSize = countRead.size();
    file.write(reinterpret_cast<char *>(&countReadSize), sizeof(countReadSize));
    for (const auto &entry : countRead) {
        size_t kmer = entry.first;
        size_t num = entry.second;
        file.write(reinterpret_cast<char *>(&kmer), sizeof(kmer));
        file.write(reinterpret_cast<char *>(&num), sizeof(num));
    }

    // Write entropy
    size_t entropySize = entropy.size();
    file.write(reinterpret_cast<char *>(&entropySize), sizeof(entropySize));
    for (const auto &entry : entropy) {
        size_t kmer = entry.first;
        float entropyValue = entry.second;
        file.write(reinterpret_cast<char *>(&kmer), sizeof(kmer));
        file_kmer.write(reinterpret_cast<char *>(&kmer), sizeof(kmer));
        file.write(reinterpret_cast<char *>(&entropyValue), sizeof(entropyValue));
    }

    // Write isRef
    file.write(reinterpret_cast<char *>(&isRef), sizeof(isRef));

    // Close the file
    file.close();
}

// Function to load KmerTable data from a binary file
void KmerTable::load_binary(const std::string &filename) {
    std::ifstream file(filename, std::ios::binary);
    if (!file.is_open()) {
        std::cerr << "Error opening file: " << filename << std::endl;
        return;
    }

    // Read K, n_gene, and n_kmer_total
    file.read(reinterpret_cast<char *>(&K), sizeof(K));
    file.read(reinterpret_cast<char *>(&n_gene), sizeof(n_gene));
    file.read(reinterpret_cast<char *>(&n_kmer_total), sizeof(n_kmer_total));

    // Read geneNameIndex
    size_t geneNameIndexSize;
    file.read(reinterpret_cast<char *>(&geneNameIndexSize), sizeof(geneNameIndexSize));
    geneNameIndex.resize(geneNameIndexSize);
    for (size_t i = 0; i < geneNameIndexSize; ++i) {
        std::getline(file, geneNameIndex[i], '\0'); // Read null-terminated strings
    }

    // Read count
    size_t countSize;
    file.read(reinterpret_cast<char *>(&countSize), sizeof(countSize));
    count.clear(); // Clear existing data
    for (size_t i = 0; i < countSize; ++i) {
        size_t gene;
        file.read(reinterpret_cast<char *>(&gene), sizeof(gene));
        size_t innerMapSize;
        file.read(reinterpret_cast<char *>(&innerMapSize), sizeof(innerMapSize));
        count[gene] = std::map<size_t, size_t>(); // Initialize inner map
        for (size_t j = 0; j < innerMapSize; ++j) {
            size_t kmer, num;
            file.read(reinterpret_cast<char *>(&kmer), sizeof(kmer));
            file.read(reinterpret_cast<char *>(&num), sizeof(num));
            count[gene][kmer] = num;
        }
    }

    // Read countRead
    size_t countReadSize;
    file.read(reinterpret_cast<char *>(&countReadSize), sizeof(countReadSize));
    countRead.clear(); // Clear existing data
    for (size_t i = 0; i < countReadSize; ++i) {
        size_t kmer, num;
        file.read(reinterpret_cast<char *>(&kmer), sizeof(kmer));
        file.read(reinterpret_cast<char *>(&num), sizeof(num));
        countRead[kmer] = num;
    }

    // Read entropy
    size_t entropySize;
    file.read(reinterpret_cast<char *>(&entropySize), sizeof(entropySize));
    entropy.clear(); // Clear existing data
    for (size_t i = 0; i < entropySize; ++i) {
        size_t kmer;
        float entropyValue;
        file.read(reinterpret_cast<char *>(&kmer), sizeof(kmer));
        file.read(reinterpret_cast<char *>(&entropyValue), sizeof(entropyValue));
        entropy[kmer] = entropyValue;
    }

    // Read isRef
    file.read(reinterpret_cast<char *>(&isRef), sizeof(isRef));

    // Close the file
    file.close();
}

void loadKmerList(const std::string& filename, vector<size_t>& kmerList) {
    std::ifstream file_kmer(filename);
    if (!file_kmer.is_open()) {
        std::cerr << "Error opening file: " << filename << std::endl;
        return;
    }
    std::string line;
    while (std::getline(file_kmer, line)) {
        std::istringstream iss(line);
        size_t kmer;
        iss >> kmer;
        kmerList.push_back(kmer);
    } 
    file_kmer.close();
}