#include "../include/namespace.h"



namespace awgn {

std::default_random_engine generator;

std::vector<float> addNoise(std::vector<int> encodedMsg, float SNR) {
  std::vector<float> noisyMsg;

  float variance = pow(10.0, -SNR / 10.0) / 2.0;
	// std::cout << "variance" << std::fixed << std::setprecision(4) << variance << std::endl;

  float sigma = sqrt(variance);
  std::normal_distribution<float> distribution(0.0, sigma);

  for (int i = 0; i < encodedMsg.size(); i++) {
    noisyMsg.push_back(encodedMsg[i] + distribution(generator));
  }
  return noisyMsg;
}

float normpdf(float x, float mu, float sigma) {
	/* returns the probability density function of the standard normal distribution, evaluated at the values in x 
	*/
	float pdf = (float)1.0/(sqrt(2*M_PI)*sigma) * std::exp(pow(((x - mu) / sigma),2) * -0.5);
	return pdf;
}

float log_normpdf(float x, float mu, float sigma) {
	/* returns the log of the probability density function of the standard normal distribution, evaluated at the values in x 
	*/
	return std::log(1.0/(sqrt(2*M_PI)*sigma)) - 0.5 * pow(((x - mu) / sigma), 2);
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

// prints a vector of floats, with commas seperating elements
void print_float_vector(std::vector<float> vector){
	if(vector.size() == 0)
		return;
	for(int i = 0; i < vector.size() - 1; i++){
		std::cout << std::setprecision(5) << vector[i] << ", ";
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

float compute_vector_energy(std::vector<float> vector){
	if (vector.size() == 0) std::cerr << "EMPTY VECTOR!" << std::endl;
	float sum_of_squares = 0.0;
	for (int i = 0; i < vector.size(); i++) {
		sum_of_squares += vector[i] * vector[i];
	}
	return sum_of_squares;
}

float compute_angle_between_vectors_rad(std::vector<float> vec1, std::vector<int> vec2) {
	// computes the angle (in radians) between a float vector and an integer vector
	// assumes the energy of the integer vector is 128.
	if (vec1.size() != vec2.size()) {std::cerr << "INVALID INNER PRODUCT DUE TO UNCOMPATIBLE SHAPE! ABORT!" << std::endl; exit(1);}
	float inner_product = std::inner_product(vec1.begin(), vec1.end(), vec2.begin(), 0.0);
	float vec1_energy_sqrt = std::sqrt(compute_vector_energy(vec1));
	float vec2_energy_sqrt = std::sqrt(128);
	float angle_rad = std::acos( inner_product/(vec1_energy_sqrt * vec2_energy_sqrt) );
	return angle_rad;
}

} // namespace utils