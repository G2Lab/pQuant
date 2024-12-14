#include "data-analysis.h"

void analyze_dataset(PQuantParams param) {

    // Convert reads to FASTA format with quality filtering
    std::cout << "=== Running seqtk to filter reads by quality ===" << std::endl;
    std::string seqtk_cmd = "seqtk seq -q " + std::to_string(param.quality) + " -l 0 " + param.filename_read + " | wc -l | awk -v N=2 '{print $1 / N}'";
    FILE* pipe = popen(seqtk_cmd.c_str(), "r");
    if (!pipe) {
        std::cerr << "Failed to run command." << std::endl;
        return;
    }
    char buffer[128];
    std::string result = "";
    while (fgets(buffer, sizeof(buffer), pipe) != nullptr) {
        result += buffer;
    }
    pclose(pipe);
    int count_reads = std::stoi(result);

    // Load kmer list
    std::cout << "=== load kmerList ===" << std::endl;
    auto start_time_kmerList = std::chrono::high_resolution_clock::now();
    vector<size_t> kmer_list;
    if (param.json_format) {
        loadKmerList(param.filename_kmerList, kmer_list);
    } else {
        loadKmerListBinary(param.filename_kmerList, kmer_list);
    }
    auto end_time_kmerList = std::chrono::high_resolution_clock::now();
    auto duration_kmerList = std::chrono::duration_cast<std::chrono::milliseconds>(end_time_kmerList - start_time_kmerList).count();
    std::cout << "load kmerList duration = " << duration_kmerList << " ms" << std::endl;
    std::cout << std::endl;

    // // Read filtered FASTA file
    // std::cout << "=== Counting reads in filtered FASTA file ===" << std::endl;
    // size_t count_reads = 0;
    // std::ifstream infile("out.fa");
    // std::string line;
    // while (std::getline(infile, line)) {
    //     if (!line.empty() && line[0] == '>') {
    //         ++count_reads;  // Each header line corresponds to a new read
    //     }
    // }
    // infile.close();

    // std::cout << "Number of reads in out.fa with quality > " << param.quality << ": " << count_reads << std::endl;

    // Read original FASTQ file
    std::cout << "=== read reads ===" << std::endl;
    auto start_time_read = std::chrono::high_resolution_clock::now();
    vector<Sequence> reads_seq;
    if (param.filename_read.find(".fa") != std::string::npos && param.filename_read.find(".fastq") == std::string::npos) {
        std::cout << "read file is .fa format" << std::endl;
        readFastaFile(param.filename_read, reads_seq);
    } else if (param.filename_read.find(".fq") != std::string::npos || param.filename_read.find(".fastq") != std::string::npos) {
        std::cout << "read file is .fq or .fastq format" << std::endl;
        readFastQFile(param.filename_read, reads_seq);
    } else {
        std::cout << "Invalid file format" << std::endl;
        std::cout << "Assume it is .fastq format" << std::endl;
        readFastQFile(param.filename_read, reads_seq);
    }
    auto end_time_read = std::chrono::high_resolution_clock::now();
    auto duration_read = std::chrono::duration_cast<std::chrono::milliseconds>(end_time_read - start_time_read).count();
    std::cout << "readFastaFile duration = " << duration_read << " ms" << std::endl;

    // Report statistics
    size_t n_kmer_total = kmer_list.size();
    size_t n_reads = reads_seq.size();

    std::cout << "Total number of kmers: " << n_kmer_total << std::endl;
    std::cout << "Total number of reads: " << n_reads << std::endl;
    std::cout << "Number of reads in out.fa with quality > " << param.quality << ": " << count_reads << std::endl;
}

void compute_tpm(PQuantParams param) {
    // Load gene vector from file
    std::string filename_geneVector = param.foldername_kmer + "/geneVector.txt";
    std::ifstream geneVectorFile(filename_geneVector);
    if (!geneVectorFile.is_open()) {
        std::cerr << "Could not open file " << filename_geneVector << std::endl;
        return;
    }

    std::vector<std::pair<std::string, double>> geneVector;
    std::string geneName;
    double value;
    while (geneVectorFile >> geneName >> value) {
        geneVector.emplace_back(geneName, value);
    }
    geneVectorFile.close();

    // load gene data
    // read fasta file (reference)
    vector<Sequence> seq_vec;
    auto start_time_readFastaFile = std::chrono::high_resolution_clock::now();
    cout << "=== read fasta file ===" << endl;
    readFastaFile(param.filename_ref, seq_vec);
    auto end_time_readFastaFile = std::chrono::high_resolution_clock::now();
    auto duration_readFastaFile = std::chrono::duration_cast<std::chrono::milliseconds>(end_time_readFastaFile - start_time_readFastaFile).count();
    std::cout << "readFastaFile duration = " << duration_readFastaFile << " ms" << std::endl;
    std::cout << std::endl;

    // compute gene length
    std::unordered_map<std::string, long> gene_lengths;
    for (auto& seq : seq_vec) {
        gene_lengths[seq.getGeneName()] = 0;
        for (long i = 0; i < seq.getNumSeq(); i++) {
            gene_lengths[seq.getGeneName()] += seq.getSeq(i).size();
        }
    }

    // Check that each gene in geneVector has a corresponding length in gene_lengths
    double totalRPK = 0.0;
    std::unordered_map<std::string, double> rpks;

    // Step 1: Calculate RPK for each gene
    for (const auto& gene : geneVector) {
        const std::string& gene_name = gene.first;
        double read_count = gene.second;

        auto it = gene_lengths.find(gene_name);
        if (it != gene_lengths.end()) {
            double gene_length_kb = it->second / 1000.0;  // Convert gene length to kilobases
            double rpk = read_count / gene_length_kb;
            rpks[gene_name] = rpk;
            totalRPK += rpk;
        } else {
            std::cerr << "Error: Gene length not found for " << gene_name << std::endl;
            return;
        }
    }

    // Step 2: Calculate TPM for each gene
    std::unordered_map<std::string, double> tpm_values;
    for (const auto& gene : rpks) {
        double tpm = (gene.second / totalRPK) * 1e6;  // Normalize to TPM
        tpm_values[gene.first] = tpm;
    }

    // Output the TPM values
    for (const auto& gene : tpm_values) {
        std::cout << "Gene: " << gene.first << ", TPM: " << gene.second << std::endl;
    }

    double thres = param.tpm_thres;

    // count # of genes with TPM > thres
    size_t count = 0;
    for (const auto& gene : tpm_values) {
        if (gene.second > thres) {
            count++;
        }
    }

    std::cout << "Number of genes with TPM > " << thres << ": " << count << std::endl;
}