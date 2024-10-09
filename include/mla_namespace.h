#ifndef MLA_NAMESPACE_H
#define MLA_NAMESPACE_H

#include <cassert>
#include <iostream>
#include <fstream>
#include <vector>
#include <random>

#include "mla_types.h"

namespace awgn {

std::vector<double> addNoise(std::vector<int> encodedMsg, double SNR);

} // namespace awgn

namespace crc {

// binary sum, used in crc_check
int bin_sum(int i, int j);

// converts decimal input to binary output, with a given number of bits
void dec_to_binary(int input, std::vector<int>& output, int bit_number);

// converts decimal output to n-bit BPSK
std::vector<int> get_point(int output, int n);

// checks the decoded message against the crc
bool crc_check(std::vector<int> input_data, int crc_bits_num, int crc_dec);

void crc_calculation(std::vector<int>& input_data, int crc_bits_num, int crc_dec);

} // namespace crc

namespace utils {

// prints a vector of doubles, with commas seperating elements
void print_double_vector(std::vector<double> vector);

// prints a vector of ints, with commas seperating elements
void print_int_vector(std::vector<int> vector);

// outputs a vector of ints to a file
void output_int_vector(std::vector<int> vector, std::ofstream& file);

} // namespace turbo_elf_utils


int make_file_interleaver(char interleaver_file[],
                          unsigned short int interleaver[], int n);

void elf_turbo_simulation(CodeInformation code);

#endif