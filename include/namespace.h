#ifndef MLA_NAMESPACE_H
#define MLA_NAMESPACE_H

#include <cassert>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <random>
#include <stdexcept>
#include <algorithm>

#include "types.h"

namespace awgn {

extern std::default_random_engine generator;

std::vector<float> addNoise(std::vector<int> encodedMsg, float esno_dB);

float normpdf(float x, float mu = 0.0, float sigma = 1.0);

float log_normpdf(float x, float mu, float sigma);

std::vector<float> addAWNGNoise(std::vector<int> transmittedMessage, std::vector<int> puncturedIndices, float esno_dB, bool noiseless);

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

std::vector<int> generateRandomCRCMessage(CodeInformation code);

} // namespace crc

namespace utils {

// prints a vector of floats, with commas seperating elements
void print_float_vector(std::vector<float> vector);

// prints a vector of ints, with commas seperating elements
void print_int_vector(std::vector<int> vector);

// outputs a vector of ints to a file
void output_int_vector(std::vector<int> vector, std::ofstream& file);

// computes vector energy, aka sum of squares
template <typename T1>
float compute_vector_energy(
    std::vector<T1> vector) 
{
    if (vector.size() == 0) std::cerr << "EMPTY VECTOR!" << std::endl;
	float sum_of_squares = 0.0;
	for (size_t i = 0; i < vector.size(); i++) {
		sum_of_squares += float(vector[i] * vector[i]);
	}
	return sum_of_squares;
}

// computes angle (radians) between two vectors
template <typename T1, typename T2>
float compute_angle_between_vectors_rad(
    const std::vector<T1>& vec1,
    const std::vector<T2>& vec2,
    const std::vector<int>& punctured_indices) 
{
    // computes the angle (in radians) between a float vector and an integer vector
	// assumes the energy of the integer vector is 128.
	if (vec1.size() != vec2.size()) {std::cerr << "INVALID INNER PRODUCT DUE TO UNCOMPATIBLE SHAPE! ABORT!" << std::endl; exit(1);}
	float inner_product = 0.0f;
    float vec1_energy_sqrt = 0.0f;
    float vec2_energy_sqrt = 0.0f;
    for (size_t i = 0; i < vec1.size(); i++) {
        if (std::find(punctured_indices.begin(), punctured_indices.end(), i) != punctured_indices.end()) {
            continue;
        } else {
            inner_product += vec1[i] * vec2[i];
            vec1_energy_sqrt += vec1[i] * vec1[i];
            vec2_energy_sqrt += vec2[i] * vec2[i];
        }
    }
	vec1_energy_sqrt = std::sqrt(vec1_energy_sqrt);
	vec2_energy_sqrt = std::sqrt(vec2_energy_sqrt);
	float angle_rad = std::acos( inner_product/(vec1_energy_sqrt * vec2_energy_sqrt) );
	return angle_rad;
}

// normalize vector to unit energy
std::vector<float> normalize_to_unit_energy(std::vector<float> vec);

// Euclidean distance metric
template <typename T1, typename T2>
float euclidean_distance(
    const std::vector<T1>& v1, 
    const std::vector<T2>& v2, 
    const std::vector<int>& punctured_indices) 
{
    if (v1.size() != v2.size()) {
        throw std::invalid_argument("Vectors must be of the same size");
        exit(1);
    }

    float sum = 0.0;
    for (size_t i = 0; i < v1.size(); i++) {
        // Skip if index is in punctured_indices
        if (std::find(punctured_indices.begin(), punctured_indices.end(), i) != punctured_indices.end()) {
            continue;
        }
        sum += std::pow(static_cast<float>(v1[i]) - static_cast<float>(v2[i]), 2);
    }
    return std::sqrt(sum);
}

// Euclidean distance metric
template <typename T1, typename T2>
float sum_of_squares(
    const std::vector<T1>& v1, 
    const std::vector<T2>& v2, 
    const std::vector<int>& punctured_indices) 
{
    if (v1.size() != v2.size()) {
        throw std::invalid_argument("Vectors must be of the same size");
    }

    float sum = 0.0;
    for (size_t i = 0; i < v1.size(); i++) {
        // Skip if index is in punctured_indices
        if (std::find(punctured_indices.begin(), punctured_indices.end(), i) != punctured_indices.end()) {
            continue;
        }
        sum += std::pow(static_cast<float>(v1[i]) - static_cast<float>(v2[i]), 2);
    }
    return sum;
}


// Element-wise Squared Distance
template <typename T1, typename T2>
std::vector<float> elementwise_squared_distance(
    const std::vector<T1>& v1, 
    const std::vector<T2>& v2, 
    const std::vector<int>& punctured_indices) 
{
    if (v1.size() != v2.size()) {
        throw std::invalid_argument("Vectors must be of the same size");
    }

    std::vector<float> distances;
    for (size_t i = 0; i < v1.size(); i++) {
        // Skip if index is in punctured_indices
        if (std::find(punctured_indices.begin(), punctured_indices.end(), i) != punctured_indices.end()) {
            continue;
        }
        distances.push_back(std::pow(static_cast<float>(v1[i]) - static_cast<float>(v2[i]), 2));
    }
    return distances;
}

} // namespace utils

#endif