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
std::vector<double> generateStandardNormalNoise(size_t l);
std::vector<double> scaleNoise(std::vector<double> input, double target_noise_power);

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

// compute vector energy (sum of squares)
double compute_vector_energy(std::vector<double> v);

// compute hamming distance between two integer vectors
int compute_hamming_distance(std::vector<int> v1, std::vector<int> v2);

int compute_hamming_distance_with_puncturing(std::vector<int> v1, std::vector<int> v2, std::vector<int> punctured_indices);


} // namespace turbo_elf_utils

#endif