#include "../include/namespace.h"



namespace awgn {

std::default_random_engine generator;


std::vector<float> addNoise(std::vector<int> encodedMsg, float esno_dB) {
  std::vector<float> noisyMsg;

  float variance = pow(10.0, -esno_dB / 10.0) / 2.0;
	// std::cout << "variance" << std::fixed << std::setprecision(4) << variance << std::endl;

  float sigma = sqrt(variance);
  std::normal_distribution<float> distribution(0.0, sigma);

	float esno = pow(10.0, esno_dB / 10.0);
	std::normal_distribution<float> llr_distribution(4*esno, std::sqrt(8*esno));

  for (size_t i = 0; i < encodedMsg.size(); i++) {
    // noisyMsg.push_back(encodedMsg[i] + distribution(generator));
		noisyMsg.push_back(encodedMsg[i] * llr_distribution(generator));
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

// this takes the transmitted message and adds AWGN noise to it
// it also punctures the bits that are not used in the trellis
std::vector<float> addAWNGNoise(std::vector<int> transmittedMessage, std::vector<int> puncturedIndices, float esno_dB, bool noiseless){
	std::vector<float> receivedMessage;
	if(noiseless){
		for(size_t i = 0; i < transmittedMessage.size(); i++)
			receivedMessage.push_back((float)transmittedMessage[i]);
	} else {
		receivedMessage = awgn::addNoise(transmittedMessage, esno_dB);
	}

	// puncture the bits. it is more convenient to puncture on this side than on the 
	// decoder, so we insert zeros which provide no information to the decoder
	for(size_t index = 0; index < puncturedIndices.size(); index++) {
		if (puncturedIndices[index] > (int)receivedMessage.size()) {
			std::cout << "out of bounds index: " << puncturedIndices[index] << std::endl;
			std::cerr << "Puncturing index out of bounds" << std::endl;
			exit(1);
		}
		receivedMessage[puncturedIndices[index]] = 0;
	}

	return receivedMessage;
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

// this generates a random binary string of length code.numInfoBits, and appends the appropriate CRC bits
std::vector<int> generateRandomCRCMessage(CodeInformation code){
	std::vector<int> info_crc;
	for(int i = 0; i < code.numInfoBits; i++)
		info_crc.push_back(rand()%2);

	// compute the CRC
	crc_calculation(info_crc, code.crcDeg, code.crc);
	return info_crc;
}

} // namespace crc

namespace utils {

// prints a vector of floats, with commas seperating elements
void print_float_vector(std::vector<float> vector){
	if(vector.size() == 0)
		return;
	for(size_t i = 0; i < vector.size() - 1; i++){
		std::cout << std::setprecision(5) << vector[i] << ", ";
	}
	std::cout << vector[vector.size() - 1] << std::endl;
}

// prints a vector of ints, with commas seperating elements
void print_int_vector(std::vector<int> vector){
	if(vector.size() == 0)
		return;
	for(size_t i = 0; i < vector.size() - 1; i++){
		std::cout << vector[i] << ", ";
	}
	std::cout << vector[vector.size() - 1] << std::endl;
}

// outputs a vector of ints to a file
void output_int_vector(std::vector<int> vector, std::ofstream& file){
	if(vector.size() == 0)
		return;
	for(size_t i = 0; i < vector.size() - 1; i++){
		file << vector[i] << ", ";
	}
	file << vector[vector.size() - 1] << std::endl;
}


std::vector<float> normalize_to_unit_energy(std::vector<float> vec) {
	if (vec.empty()) return vec;
	// Compute energy (sum of squares)
	double energy = 0.0;
	for (float v : vec) {
			energy += static_cast<double>(v) * v;
	}
	if (energy == 0.0) {
			throw std::runtime_error("Cannot normalize: zero energy vector.");
	}
	// Scale factor is 1/sqrt(energy)
	double scale = 1.0 / std::sqrt(energy);
	for (float& v : vec) {
			v = static_cast<float>(v * scale);
	}
	return vec;
}

} // namespace utils