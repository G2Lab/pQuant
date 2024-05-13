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
            // if not AGCT, then continue with the same character
            ans = nt[i] + ans;
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
        map<size_t, map<size_t, size_t>> count_total;
        if (geneNameIndex.size() == 0) {
            // save gene list into geneNameIndex
            for (size_t i = 0; i < gene.size(); i++) {
                geneNameIndex.push_back(gene[i].getGeneName());
            }
        }

        auto start_time = std::chrono::high_resolution_clock::now();
        map<size_t, map<size_t, size_t>> kmer_occurance;
        // count kmers in gene
        size_t start = 0;
        size_t end = gene.size();
        if (param.gene_start >= 0 && param.gene_end >= 0) {
            start = param.gene_start;
            end = param.gene_end;
        }
        for (size_t i = start; i < end; i++) {
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
                    if (seq_j.size() == 0 || seq_k.size() == 0) {
                        continue;
                    }
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
            count_total.insert(make_pair(i, gene_kmer_count));
            gene_kmer_count.clear();

            if (param.progress_bar) {
                print_progress_bar("countKmersPerGene", i, end - start, start_time);
            }
        }

        std::map<size_t, float> entropy_total;
        // compute entropy
        int progress_count = 0;
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
            if (param.progress_bar) {
                print_progress_bar("computeEntropy", progress_count, kmer_occurance.size(), start_time);
                progress_count += 1;
            }
        }

        // filter entropy
        progress_count = 0;
        for (auto it = entropy_total.begin(); it != entropy_total.end(); ++it) {
            if (it->second < thres) {
                entropy.insert(make_pair(it->first, it->second));
            }
            if (param.progress_bar) {
                print_progress_bar("filterEntropy", progress_count, entropy_total.size(), start_time);
                progress_count += 1;
            }
        }
        entropy_total.clear();

        // filter kmer_occurance
        progress_count = 0;
        for (auto it = kmer_occurance.begin(); it != kmer_occurance.end(); ++it) {
            if (entropy.find(it->first) != entropy.end()) {
                tableRef.insert(make_pair(it->first, it->second));
            }
            if (param.progress_bar) {
                print_progress_bar("filterKmerOccurance", progress_count, kmer_occurance.size(), start_time);
                progress_count += 1;
            }
        }
        n_kmer_total = entropy.size();
        kmer_occurance.clear();

        // filter count_total
        progress_count = 0;
        for (auto it = count_total.begin(); it != count_total.end(); ++it) {
            map<size_t, size_t> gene_kmer_count = it->second;
            map<size_t, size_t> gene_kmer_count_filtered;
            for (auto it2 = gene_kmer_count.begin(); it2 != gene_kmer_count.end(); ++it2) {
                if (entropy.find(it2->first) != entropy.end()) {
                    gene_kmer_count_filtered.insert(make_pair(it2->first, it2->second));
                }
            }
            count.insert(make_pair(it->first, gene_kmer_count_filtered));
            if (param.progress_bar) {
                print_progress_bar("filterCountTotal", progress_count, count_total.size(), start_time);
                progress_count += 1;
            }
        }
        count_total.clear();
    } else {
        auto start_time = std::chrono::high_resolution_clock::now();
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
                    countRead[num] += 1;
                }
                if (param.progress_bar) {
                    print_progress_bar("Encode strings", j, gene[i].getNumSeq(), start_time);
                }
            }
        }
    }
}

bool KmerTable::operator==(const KmerTable& other) const {
    if (K != other.K || thres != other.thres || n_gene != other.n_gene || n_kmer_total != other.n_kmer_total || isRef != other.isRef) {
        cout << "K, thres, n_gene, n_kmer_total, or isRef is different\n"
        << "K: " << K << ", " << other.K << "\n"
        << "thres: " << thres << ", " << other.thres << "\n"
        << "n_gene: " << n_gene << ", " << other.n_gene << "\n"
        << "n_kmer_total: " << n_kmer_total << ", " << other.n_kmer_total << "\n"
        << "isRef: " << isRef << ", " << other.isRef << endl;
        return false;
    }
    if (geneNameIndex.size() != other.geneNameIndex.size()) {
        cout << "geneNameIndex size is different\n"
        << "geneNameIndex.size(): " << geneNameIndex.size() << ", " << other.geneNameIndex.size() << endl;
        return false;
    }
    for (size_t i = 0; i < geneNameIndex.size(); i++) {
        if (geneNameIndex[i] != other.geneNameIndex[i]) {
            cout << "geneNameIndex[" << i << "] is different\n"
            << "geneNameIndex[" << i << "]: " << geneNameIndex[i] << ", " << other.geneNameIndex[i] << endl;
            return false;
        }
    }
    if (count.size() != other.count.size()) {
        cout << "count size is different\n"
        << "count.size(): " << count.size() << ", " << other.count.size() << endl;
        return false;
    }
    for (auto &p : count) {
        if (other.count.find(p.first) == other.count.end()) {
            cout << "gene " << p.first << " is not found in other.count\n";
            return false;
        }
        if (p.second.size() != other.count.at(p.first).size()) {
            cout << "p.second.size() is different\n"
            << "p.second.size(): " << p.second.size() << ", " << other.count.at(p.first).size() << endl;
            return false;
        }
        for (auto &q : p.second) {
            if (other.count.at(p.first).find(q.first) == other.count.at(p.first).end()) {
                cout << "kmer " << q.first << " is not found in other.count[" << p.first << "]\n";
                return false;
            }
            if (q.second != other.count.at(p.first).at(q.first)) {
                cout << "count is different\n"
                << "count[" << p.first << "][" << q.first << "]: " << q.second << ", " << other.count.at(p.first).at(q.first) << endl;
                return false;
            }
        }
    }

    if (tableRef.size() != other.tableRef.size()) {
        cout << "tableRef size is different\n"
        << "tableRef.size(): " << tableRef.size() << ", " << other.tableRef.size() << endl;
        return false;
    }
    for (auto &p : tableRef) {
        if (other.tableRef.find(p.first) == other.tableRef.end()) {
            cout << "kmer " << p.first << " is not found in other.tableRef\n";
            return false;
        }
        if (p.second.size() != other.tableRef.at(p.first).size()) {
            cout << "p.second.size() is different\n"
            << "p.second.size(): " << p.second.size() << ", " << other.tableRef.at(p.first).size() << endl;
            return false;
        }
        for (auto &q : p.second) {
            if (other.tableRef.at(p.first).find(q.first) == other.tableRef.at(p.first).end()) {
                cout << "gene " << q.first << " is not found in other.tableRef[" << p.first << "]\n";
                return false;
            }
            if (q.second != other.tableRef.at(p.first).at(q.first)) {
                cout << "count is different\n"
                << "count[" << p.first << "][" << q.first << "]: " << q.second << ", " << other.tableRef.at(p.first).at(q.first) << endl;
                return false;
            }
        }
    }

    if (countRead.size() != other.countRead.size()) {
        cout << "countRead size is different\n"
        << "countRead.size(): " << countRead.size() << ", " << other.countRead.size() << endl;
        return false;
    }
    for (auto &p : countRead) {
        if (other.countRead.find(p.first) == other.countRead.end()) {
            cout << "kmer " << p.first << " is not found in other.countRead\n";
            return false;
        }
        if (p.second != other.countRead.at(p.first)) {
            cout << "counRead is different\n"
            << "counRead[" << p.first << "]: " << p.second << ", " << other.countRead.at(p.first) << endl;
            return false;
        }
    }
    if (entropy.size() != other.entropy.size()) {
        cout << "entropy size is different\n"
        << "entropy.size(): " << entropy.size() << ", " << other.entropy.size() << endl;
        return false;
    }
    for (auto &p : entropy) {
        if (other.entropy.find(p.first) == other.entropy.end()) {
            cout << "kmer " << p.first << " is not found in other.entropy\n";
            return false;
        }
        if (p.second - other.entropy.at(p.first) >= 0.0001) {
            cout << "entropy is different\n"
            << "entropy[" << p.first << "]: " << p.second << ", " << other.entropy.at(p.first) << endl;
            return false;
        }
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
    file << "isRef: " << isRef << std::endl;
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
        file_kmer << p.first << " ";
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
    isRef = (temp == "1");

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
            // remove ' ' from prefix of temp
            if (temp[0] == ' ') {
                temp = temp.substr(1, temp.size());
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
        while (iss6 >> seperator) {
            if (seperator == ',' || seperator == '{') {
                size_t key, val;
                iss6 >> key >> seperator >> val;
                kmerData[key] = val;
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
void KmerTable::saveBinary(const std::string& filename) {
    std::ofstream file(filename, std::ios::out | std::ios::binary);
    if (!file.is_open()) {
        std::cerr << "Error opening file: " << filename << std::endl;
        return;
    }
    // Save K, thres, n_gene, and n_kmer_total as binary
    file.write(reinterpret_cast<const char*>(&K), sizeof(K));
    file.write(reinterpret_cast<const char*>(&thres), sizeof(thres));
    file.write(reinterpret_cast<const char*>(&n_gene), sizeof(n_gene));
    file.write(reinterpret_cast<const char*>(&n_kmer_total), sizeof(n_kmer_total));
    file.write(reinterpret_cast<const char*>(&isRef), sizeof(isRef));

    // Save geneNameIndex as binary
    size_t geneNameIndexSize = geneNameIndex.size();
    file.write(reinterpret_cast<const char*>(&geneNameIndexSize), sizeof(geneNameIndexSize));
    for (const std::string& geneName : geneNameIndex) {
        size_t geneNameSize = geneName.size();
        file.write(reinterpret_cast<const char*>(&geneNameSize), sizeof(geneNameSize));
        file.write(geneName.c_str(), geneNameSize);
    }
    
    // Save count as binary
    size_t countSize = count.size();
    file.write(reinterpret_cast<const char*>(&countSize), sizeof(countSize));
    for (auto& p : count) {
        file.write(reinterpret_cast<const char*>(&p.first), sizeof(p.first));
        size_t sizeP = p.second.size();
        file.write(reinterpret_cast<const char*>(&sizeP), sizeof(sizeP));
        for (auto& q : p.second) {
            file.write(reinterpret_cast<const char*>(&q.first), sizeof(q.first));
            file.write(reinterpret_cast<const char*>(&q.second), sizeof(q.second));
        }
    }

    // Save tableRef as binary
    size_t tableRefSize = tableRef.size();
    file.write(reinterpret_cast<const char*>(&tableRefSize), sizeof(tableRefSize));
    for (auto& p : tableRef) {
        file.write(reinterpret_cast<const char*>(&p.first), sizeof(p.first));
        size_t sizeP = p.second.size();
        file.write(reinterpret_cast<const char*>(&sizeP), sizeof(sizeP));
        for (auto& q : p.second) {
            file.write(reinterpret_cast<const char*>(&q.first), sizeof(q.first));
            file.write(reinterpret_cast<const char*>(&q.second), sizeof(q.second));
        }
    }

    // Save entropy as binary
    size_t entropySize = entropy.size();
    file.write(reinterpret_cast<const char*>(&entropySize), sizeof(entropySize));
    for (auto& p : entropy) {
        file.write(reinterpret_cast<const char*>(&p.first), sizeof(p.first));
        file.write(reinterpret_cast<const char*>(&p.second), sizeof(p.second));
    }

    file.close();
}

void KmerTable::loadBinary(const std::string& filename) {
    std::ifstream file(filename, std::ios::in | std::ios::binary);
    if (!file.is_open()) {
        std::cerr << "Error opening file: " << filename << std::endl;
        return;
    }

    // Load K, thres, n_gene, and n_kmer_total from binary
    file.read(reinterpret_cast<char*>(&K), sizeof(K));
    file.read(reinterpret_cast<char*>(&thres), sizeof(thres));
    file.read(reinterpret_cast<char*>(&n_gene), sizeof(n_gene));
    file.read(reinterpret_cast<char*>(&n_kmer_total), sizeof(n_kmer_total));
    file.read(reinterpret_cast<char*>(&isRef), sizeof(isRef));

    // Load geneNameIndex from binary
    size_t geneNameIndexSize;
    file.read(reinterpret_cast<char*>(&geneNameIndexSize), sizeof(geneNameIndexSize));
    geneNameIndex.resize(geneNameIndexSize);
    for (size_t i = 0; i < geneNameIndexSize; ++i) {
        size_t geneNameSize;
        file.read(reinterpret_cast<char*>(&geneNameSize), sizeof(geneNameSize));
        std::string geneName(geneNameSize, '\0');
        file.read(&geneName[0], geneNameSize);
        geneNameIndex[i] = geneName;
    }

    // Load count from binary
    size_t countSize;
    file.read(reinterpret_cast<char*>(&countSize), sizeof(countSize));
    for (size_t i = 0; i < countSize; ++i) {
        size_t key;
        file.read(reinterpret_cast<char*>(&key), sizeof(key));
        size_t mapSize;
        file.read(reinterpret_cast<char*>(&mapSize), sizeof(mapSize));
        std::map<size_t, size_t> map_dummy;
        for (size_t j = 0; j < mapSize; ++j) {
            size_t innerKey, innerValue;
            file.read(reinterpret_cast<char*>(&innerKey), sizeof(innerKey));
            file.read(reinterpret_cast<char*>(&innerValue), sizeof(innerValue));
            map_dummy[innerKey] = innerValue;
        }
        count[key] = map_dummy;
    }

    // Load tableRef from binary
    size_t tableRefSize;
    file.read(reinterpret_cast<char*>(&tableRefSize), sizeof(tableRefSize));
    for (size_t i = 0; i < tableRefSize; ++i) {
        size_t key;
        file.read(reinterpret_cast<char*>(&key), sizeof(key));
        size_t mapSize;
        file.read(reinterpret_cast<char*>(&mapSize), sizeof(mapSize));
        std::map<size_t, size_t> map_;
        for (size_t j = 0; j < mapSize; ++j) {
            size_t innerKey, innerValue;
            file.read(reinterpret_cast<char*>(&innerKey), sizeof(innerKey));
            file.read(reinterpret_cast<char*>(&innerValue), sizeof(innerValue));
            map_[innerKey] = innerValue;
        }
        tableRef[key] = map_;
    }

    // Load entropy from binary
    size_t entropySize;
    file.read(reinterpret_cast<char*>(&entropySize), sizeof(entropySize));
    for (size_t i = 0; i < entropySize; ++i) {
        size_t key;
        double value;
        file.read(reinterpret_cast<char*>(&key), sizeof(key));
        file.read(reinterpret_cast<char*>(&value), sizeof(value));
        entropy[key] = value;
    }

    file.close();
}

void KmerTable::filterGenes(size_t start, size_t end) {
    n_gene = end - start + 1;
    geneNameIndex.erase(std::next(geneNameIndex.begin(), end + 1), geneNameIndex.end());
    geneNameIndex.erase(geneNameIndex.begin(), std::next(geneNameIndex.begin(), start));
    if (n_gene != geneNameIndex.size()) {
        std::cerr << "KmerTable::filterGenes::Error: n_gene is not equal to geneNameIndex.size()" << std::endl;
        exit(0);
    }
    cout << "Used gene list" << endl;
    for (size_t i = 0; i < n_gene; i++) {
        cout << geneNameIndex[i] << endl;
    }
}

void loadKmerList(const std::string& filename, vector<size_t>& kmerList) {
    std::ifstream file_kmer(filename);
    if (!file_kmer.is_open()) {
        std::cerr << "Error opening file: " << filename << std::endl;
        return;
    }
    std::string line;
    std::getline(file_kmer, line);
    std::istringstream iss(line);
    size_t kmer;
    while(iss >> kmer) {
        kmerList.push_back(kmer);
    } 
    file_kmer.close();
}

void KmerTable::saveKmerListBinary(const std::string& filename) {
    std::ofstream file(filename, std::ios::out | std::ios::binary);
    if (!file.is_open()) {
        std::cerr << "Error opening file: " << filename << std::endl;
        return;
    }

    // Write the size of the entropy map
    size_t entropySize = entropy.size();
    file.write(reinterpret_cast<const char*>(&entropySize), sizeof(entropySize));

    // Write each key-value pair in the entropy map
    for (const auto& p : entropy) {
        file.write(reinterpret_cast<const char*>(&p.first), sizeof(p.first));
    }

    file.close();
}

// Binary load function for entropy data
void loadKmerListBinary(const std::string& filename, vector<size_t>& kmerList) {
    std::ifstream file(filename, std::ios::in | std::ios::binary);
    if (!file.is_open()) {
        std::cerr << "Error opening file: " << filename << std::endl;
        return;
    }

    // Read the size of the entropy map
    size_t entropySize;
    file.read(reinterpret_cast<char*>(&entropySize), sizeof(entropySize));

    // Read each key-value pair and populate the entropy map
    for (size_t i = 0; i < entropySize; ++i) {
        size_t key;
        file.read(reinterpret_cast<char*>(&key), sizeof(key));
        kmerList.push_back(key);
    }

    file.close();
}