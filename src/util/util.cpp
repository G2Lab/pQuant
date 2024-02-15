#include "util/util.h"

bool binary_search_util(std::vector<size_t>& v, size_t target, size_t& index) {
    size_t start = 0;
    size_t end = v.size();
    while(start + 1 < end) {
        size_t pivot = (start + end) / 2;
        if (target <= v[pivot]) {
            end = pivot;
        } else {
            start = pivot;
        }
    }
    if (v[start] == target) {
        index = start;
        return true;
    } else if (v[end] == target) {
        index = end;
        return true;
    } else {
        index = 0;
        return false;
    }
}

std::string ElapsedTime(std::chrono::seconds secs) {
    using namespace std;
    using namespace std::chrono;

    bool neg = secs < 0s;
    if (neg)
        secs = -secs;

    auto h = duration_cast<hours>(secs);
    secs -= h;
    auto m = duration_cast<minutes>(secs);
    secs -= m;

    stringstream result;

    if (neg)
        result << '-';
    
    result << setw(2) << setfill('0') << h.count() << ':';
    result << setw(2) << setfill('0') << m.count() << ':';
    result << setw(2) << setfill('0') << secs.count();

    return result.str();
}

void print_progress_bar(std::string instruction, int iter, int total, std::chrono::time_point<std::chrono::high_resolution_clock> start_time, int bar_width, bool print) {
    long divide = total < 100 ? 1 : total / 100;
    if (iter == 0) {
        std::cout << std::endl;
    }
    if (print && (iter + 1) % divide == 0) {
        // Calculate progress and time
        double progress = (double)(iter + 1) / (double)total;
        auto current_time = std::chrono::high_resolution_clock::now();
        auto elapsed_time = std::chrono::duration_cast<std::chrono::seconds>(current_time - start_time).count();

        // Calculate estimated time remaining
        auto average_time_per_iteration = elapsed_time / (double)(iter + 1);
        auto estimated_time_remaining = average_time_per_iteration * (total - iter - 1);

        // Print progress bar
        std::cout << "\r " << instruction << " : [";
        int pos = bar_width * progress;
        for (int j = 0; j < bar_width; ++j) {
            if (j < pos)
                std::cout << "=";
            else if (j == pos)
                std::cout << ">";
            else
                std::cout << " ";
        }
        std::cout << "] ";
        std::cout << " (" << iter + 1 << "/" << total << ") ";
        std::cout << std::setw(3) << int(progress * 100.0) << "% ";

        // Print time information
        std::cout << "| Elapsed Time: " << ElapsedTime(std::chrono::seconds(elapsed_time)) << "s ";
        std::cout << "| ETA: " << ElapsedTime(std::chrono::seconds(static_cast<long long>(estimated_time_remaining))) << "s ";
        cout << "Memory : " << (double) getMemoryUsage() / (1024.0 * 1024.0 * 1024.0) << " GB";

        std::cout.flush();
    }
    if (iter == total - 1) {
        std::cout << std::endl;
    }
}

size_t getMemoryUsage() {
    size_t memory_usage = 0;

    struct rusage r_usage;
    getrusage(RUSAGE_SELF, &r_usage);
    memory_usage = static_cast<size_t>(r_usage.ru_maxrss) * 1024; // Convert to bytes

    return memory_usage;
}

size_t getPeakMemoryUsage() {
    size_t peak_memory_usage = 0;

    // Open the status file for the current process
    std::ifstream status_file("/proc/self/status");
    std::string line;

    // Search for the line containing peak memory usage
    while (std::getline(status_file, line)) {
        if (line.find("VmPeak") != std::string::npos) {
            std::istringstream iss(line);
            std::string label;
            size_t value;
            iss >> label >> value;
            peak_memory_usage = value * 1024; // Convert to bytes
            break;
        }
    }

    return peak_memory_usage;
}

void printMemoryUsage() {
    size_t memory_usage_mb = getMemoryUsage() / (1024. * 1024.);
    size_t pick_memory_usage = getPeakMemoryUsage() / (1024. * 1024.);
    std::string pick_memory_surfix = "MB";
    if (pick_memory_usage > 1024) {
        pick_memory_usage = pick_memory_usage / 1024;
        pick_memory_surfix = "GB";
    }

    if (memory_usage_mb < 1024) {
        std::cout << "Memory usage = " << memory_usage_mb << " MB (pick = " << pick_memory_usage << " " << pick_memory_surfix << ")" << std::endl;
    } else {
        std::cout << "Memory usage = " << memory_usage_mb / 1024. << " GB (pick = " << pick_memory_usage << " " << pick_memory_surfix << ")" << std::endl;
    }
}

long encode_nt_to_num(string nt) {
    long ans = 0;
    for (size_t i = 0; i < nt.size(); i++) {
        char a = nt[i];
        long num;
        switch (a) {
            case 'A':
                num = 0;
                break;
            case 'C':
                num = 1;
                break;
            case 'G':
                num = 2;
                break;
            case 'T':
                num = 3;
                break;
            default:
                num = -1;
                break;
        }
        ans = ans * 4 + num;
    }
    return ans;
}

template <typename T>
inline void print_vector(std::vector<T> vec, std::size_t print_size, int prec) {
    /*
    Save the formatting information for std::cout.
    */
    std::ios old_fmt(nullptr);
    old_fmt.copyfmt(std::cout);

    std::size_t slot_count = vec.size();

    std::cout << std::fixed << std::setprecision(prec);
    std::cout << std::endl;
    if (slot_count <= 2 * print_size) {
        std::cout << "    [";
        for (std::size_t i = 0; i < slot_count; i++) {
            std::cout << " " << vec[i] << ((i != slot_count - 1) ? "," : " ]\n");
        }
    } else {
        vec.resize(std::max(vec.size(), 2 * print_size));
        std::cout << "    [";
        for (std::size_t i = 0; i < print_size; i++) {
            std::cout << " " << vec[i] << ",";
        }
        if (vec.size() > 2 * print_size) {
            std::cout << " ...,";
        }
        for (std::size_t i = slot_count - print_size; i < slot_count; i++) {
            std::cout << " " << vec[i] << ((i != slot_count - 1) ? "," : " ]\n");
        }
    }
    std::cout << std::endl;

    /*
    Restore the old std::cout formatting.
    */
    std::cout.copyfmt(old_fmt);
}

template <typename T>
inline void print_matrix(std::vector<T> matrix, std::size_t row_size) {
    /*
    We're not going to print every column of the matrix (there are 2048).
    Instead print this many slots from beginning and end of the matrix.
    */
    std::size_t print_size = 5;

    std::cout << std::endl;
    std::cout << "    [";
    for (std::size_t i = 0; i < print_size; i++) {
        std::cout << std::setw(3) << std::right << matrix[i] << ",";
    }
    std::cout << std::setw(3) << " ...,";
    for (std::size_t i = row_size - print_size; i < row_size; i++) {
        std::cout << std::setw(3) << matrix[i] << ((i != row_size - 1) ? "," : " ]\n");
    }
    std::cout << "    [";
    for (std::size_t i = row_size; i < row_size + print_size; i++) {
        std::cout << std::setw(3) << matrix[i] << ",";
    }
    std::cout << std::setw(3) << " ...,";
    for (std::size_t i = 2 * row_size - print_size; i < 2 * row_size; i++) {
        std::cout << std::setw(3) << matrix[i]
                  << ((i != 2 * row_size - 1) ? "," : " ]\n");
    }
    std::cout << std::endl;
}

void printFileSize(const std::string& filename) {
    ifstream file(filename);
    if (!file.is_open()) {
        throw runtime_error("Error: Failed to open file " + filename);
    } else {
        file.seekg(0, ios::end);
        std::streampos filesize = file.tellg();
        double filesizeKB = static_cast<double>(filesize) / 1024;
        // print filesize. If the file size is large, print as GB
        cout << "filesize of " << filename << " = ";
        if (filesizeKB < 1024) {
            cout << filesizeKB << " KB" << endl;
        } else {
            double filesizeMB = filesizeKB / 1024;
            if (filesizeMB < 1024) {
                cout << filesizeMB << " MB" << endl;
            } else {
                double filesizeGB = filesizeMB / 1024;
                cout << filesizeGB << " GB" << endl;
            }
            
        }
    }
}

void printFolderSize(const std::string& foldername) {
    std::string command = "du -sh " + foldername;
    int status = system(command.c_str());
    if (status != 0) {
        std::cerr << "Error executing command: " << command << std::endl;
    }
}