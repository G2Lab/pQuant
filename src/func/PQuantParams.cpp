#include "func/PQuantParams.h"


PQuantParams::PQuantParams(std::string target, std::string filename_read, std::string filename_ref, std::string out_path, long k, bool verbose, bool progress_bar, bool serial) {
    this->k = k;
    this->filename_read = filename_read;
    this->filename_ref = filename_ref;
    this->out_path = out_path;
    this->verbose = verbose;
    this->progress_bar = progress_bar;
    this->target = target;
    this->serial = serial;
}

void PQuantParams::print() {
    std::cout << "filename_read = " << this->filename_read << std::endl;
    std::cout << "filename_ref = " << this->filename_ref << std::endl;
    std::cout << "out_path = " << this->out_path << std::endl;
    
    std::cout << "k = " << this->k << std::endl;
    std::cout << "target = " << this->target << std::endl;
    std::cout << "verbose = " << this->verbose << std::endl;
    std::cout << "serial = " << this->serial << std::endl;
    std::cout << "progress_bar = " << this->progress_bar << std::endl;
}
