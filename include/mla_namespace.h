#ifndef MLA_NAMESPACE_H
#define MLA_NAMESPACE_H

#include <cassert>
#include <iostream>
#include <fstream>
#include <vector>
#include <random>
#include <stdexcept>
#include <algorithm>

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

// computes vector energy, aka sum of squares
double compute_vector_energy(std::vector<double> vector);

// Euclidean distance metric
template <typename T1, typename T2>
double euclidean_distance(
    const std::vector<T1>& v1, 
    const std::vector<T2>& v2, 
    const std::vector<int>& punctured_indices) 
{
    if (v1.size() != v2.size()) {
        throw std::invalid_argument("Vectors must be of the same size");
    }

    double sum = 0.0;
    for (size_t i = 0; i < v1.size(); i++) {
        // Skip if index is in punctured_indices
        if (std::find(punctured_indices.begin(), punctured_indices.end(), i) != punctured_indices.end()) {
            continue;
        }
        sum += std::pow(static_cast<double>(v1[i]) - static_cast<double>(v2[i]), 2);
    }
    return std::sqrt(sum);
}

// Euclidean distance metric
template <typename T1, typename T2>
double sum_of_squares(
    const std::vector<T1>& v1, 
    const std::vector<T2>& v2, 
    const std::vector<int>& punctured_indices) 
{
    if (v1.size() != v2.size()) {
        throw std::invalid_argument("Vectors must be of the same size");
    }

    double sum = 0.0;
    for (size_t i = 0; i < v1.size(); i++) {
        // Skip if index is in punctured_indices
        if (std::find(punctured_indices.begin(), punctured_indices.end(), i) != punctured_indices.end()) {
            continue;
        }
        sum += std::pow(static_cast<double>(v1[i]) - static_cast<double>(v2[i]), 2);
    }
    return sum;
}


// Element-wise Squared Distance
template <typename T1, typename T2>
std::vector<double> elementwise_squared_distance(
    const std::vector<T1>& v1, 
    const std::vector<T2>& v2, 
    const std::vector<int>& punctured_indices) 
{
    if (v1.size() != v2.size()) {
        throw std::invalid_argument("Vectors must be of the same size");
    }

    std::vector<double> distances;
    for (size_t i = 0; i < v1.size(); i++) {
        // Skip if index is in punctured_indices
        if (std::find(punctured_indices.begin(), punctured_indices.end(), i) != punctured_indices.end()) {
            continue;
        }
        distances.push_back(std::pow(static_cast<double>(v1[i]) - static_cast<double>(v2[i]), 2));
    }
    return distances;
}

} // namespace utils

#endif