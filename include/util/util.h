#ifndef util_h_
#define util_h_

#include <sys/wait.h>
#include <sys/resource.h>
#include <unistd.h>

#include <chrono>
#include <cmath>
#include <functional>
#include <iostream>
#include <string>
#include <unordered_map>
#include <vector>
#include <iomanip>
#include <chrono>
#include <ctime>

#include "util/date.h"
#include "util/FileReader.h"
using namespace std;

size_t getMemoryUsage();

void printMemoryUsage();

long encode_nt_to_num(string nt);

template <typename T>
inline void print_vector(std::vector<T> vec, std::size_t print_size, int prec);

template <typename T>
void print_matrix(std::vector<T> matrix, std::size_t row_size);

void print_progress_bar(std::string instruction, int iter, int total, std::chrono::time_point<std::chrono::high_resolution_clock> start_time, int=70, bool=true);

std::string ElapsedTime(std::chrono::seconds secs);

#endif