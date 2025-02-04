#include "../include/mla_namespace.h"


namespace awgn {

std::default_random_engine generator(41);


std::vector<double> addNoise(std::vector<int> encodedMsg, double SNR) {
  std::vector<double> noisyMsg;

  double variance = pow(10.0, -SNR / 10.0);
  double sigma = sqrt(variance);
  std::normal_distribution<double> distribution(0.0, sigma);

  for (int i = 0; i < encodedMsg.size(); i++) {
    noisyMsg.push_back(encodedMsg[i] + distribution(generator));
  }
  return noisyMsg;
}

std::vector<double> generateStandardNormalNoise(size_t l) {
	/* Generate I.I.D Standard Normal Noise of vector length l */
	std::vector<double> out(l, 0.0);
	std::normal_distribution<double> standard_normal_distribution(0.0, 1);
	for (size_t i = 0; i < l; i++) {
		out[i] = standard_normal_distribution(generator);
	}
	return out;
}

std::vector<double> scaleNoise(std::vector<double> input, double target_noise_power) {
	/* Scales input noise vector to have target_noise_power */
	int l = input.size();
	double input_sum_of_squares = 0.0;
	for (size_t i = 0; i < l; i++) {
		input_sum_of_squares += input[i] * input[i];
	}
	
	double scale = std::sqrt(target_noise_power / input_sum_of_squares);
	
	
	std::vector<double> out(l, 0.0);
	for (size_t i = 0; i < input.size(); i++) {
		out[i] = input[i] * scale;
	}

	return out;
}

} // namespace awgn

namespace crc {

// converts decimal input to binary output, with a given number of bits
// since we need to keep track of leading zeros
void dec_to_binary(int input, std::vector<int>& output, int bit_number) {
	output.assign(bit_number, -1);
	for (int i = bit_number - 1; i >= 0; i--) {
		int k = input >> i;
		if (k & 1)
			output[bit_number - 1 - i] = 1;
		else
			output[bit_number - 1 - i] = 0;
	}
}

// converts decimal output to n-bit BPSK
std::vector<int> get_point(int output, int n) {
	std::vector<int> bin_output;
	dec_to_binary(output, bin_output, n);
	for (int i=0; i<n; i++){
		bin_output[i] = -2 * bin_output[i] + 1;
	}
	return bin_output;
}


// binary sum, used in crc_check
int bin_sum(int i, int j) {
	return (i + j) % 2;
}

// checks the decoded message against the crc
bool crc_check(std::vector<int> input_data, int crc_bits_num, int crc_dec) {
	std::vector<int> CRC;
	dec_to_binary(crc_dec, CRC, crc_bits_num);

	for (int ii = 0; ii <= (int)input_data.size() - crc_bits_num; ii++) {
		if (input_data[ii] == 1) {
			// Note: transform doesn't include .end
			std::transform(input_data.begin() + ii, input_data.begin() + (ii + crc_bits_num), CRC.begin(), input_data.begin() + ii, bin_sum);
		}
	}
	bool zeros = std::all_of(input_data.begin(), input_data.end(), [](int i) { return i == 0; });
	return zeros;
}

void crc_calculation(std::vector<int>& input_data, int crc_bits_num, int crc_dec){
	// crc_bits_num: the number of CRC bits, redundancy bits number is 1 less.
	int length = (int)input_data.size();
	std::vector<int> CRC;
	dec_to_binary(crc_dec, CRC, crc_bits_num);
	input_data.resize(length + crc_bits_num - 1, 0);

	std::vector<int> output_data = input_data;
	for (int ii = 0; ii <= length - 1; ii++)
	{
		if (output_data[ii] == 1)
		{
			std::transform(output_data.begin() + ii, output_data.begin() + (ii + crc_bits_num), CRC.begin(), output_data.begin() + ii, bin_sum);
		}
	}

	for (int ii = length; ii < (int)output_data.size(); ii++){ input_data[ii] = output_data[ii];}
}

} // namespace crc

namespace utils {

// prints a vector of doubles, with commas seperating elements
void print_double_vector(std::vector<double> vector){
	if(vector.size() == 0)
		return;
	for(int i = 0; i < vector.size() - 1; i++){
		std::cout << vector[i] << ", ";
	}
	std::cout << vector[vector.size() - 1] << std::endl;
}

// prints a vector of ints, with commas seperating elements
void print_int_vector(std::vector<int> vector){
	if(vector.size() == 0)
		return;
	for(int i = 0; i < vector.size() - 1; i++){
		std::cout << vector[i] << ", ";
	}
	std::cout << vector[vector.size() - 1] << std::endl;
}

// outputs a vector of ints to a file
void output_int_vector(std::vector<int> vector, std::ofstream& file){
	if(vector.size() == 0)
		return;
	for(int i = 0; i < vector.size() - 1; i++){
		file << vector[i] << ", ";
	}
	file << vector[vector.size() - 1] << std::endl;
}

// computes vector energy
double compute_vector_energy ( std::vector<double> v) {
	double sum_of_squares = 0.0;
	for ( int i = 0; i < v.size(); i++ ) {
		sum_of_squares += v[i] * v[i];
	}
	
	return std::sqrt(sum_of_squares);
}

} // namespace turbo_elf_utils