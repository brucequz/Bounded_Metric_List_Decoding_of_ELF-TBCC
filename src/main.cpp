#include <iostream>
#include <fstream>
#include <vector>

#include "../include/mla_consts.h"
#include "../include/mla_types.h"
#include "../include/mla_namespace.h"
#include "../include/feedForwardTrellis.h"
#include "../include/lowRateListDecoder.h"

void Noise_injection_sim(CodeInformation code, double injected_power);
std::vector<int> generateRandomCRCMessage(CodeInformation code, bool noiseless);
std::vector<int> generateTransmittedMessage(std::vector<int> originalMessage, FeedForwardTrellis encodingTrellis);
std::vector<double> addAWNGNoiseAndPuncture(std::vector<int> transmittedMessage, std::vector<int> puncturedIndices, double snr, bool noiseless);

int main() {

  std::cout << "K: " << K << std::endl;
  std::cout << "N: " << N << std::endl;
  std::cout << "V: " << V << std::endl;
  std::cout << "M: " << M << std::endl;
    
  CodeInformation code;
  code.k = K;         // numerator of the rate
  code.n = N;         // denominator of the rate
  code.v = V;         // number of memory elements
  code.crcDeg = M+1;  // m+1, degree of CRC, # bits of CRC polynomial
  code.crc = CRC;     // CRC polynomial
  code.numInfoBits = NUM_INFO_BITS; // number of information bits
  code.numerators = {POLY1, POLY2};

  Noise_injection_sim(code, TARGET_NOISE_ENERGY);

  return 0;
}

void Noise_injection_sim(CodeInformation code, double injected_power) {
	/** Injecting noise with certain power to a fixed codeword and observe list size
	 * 
	 * Input:
	 * 	- code: code information such as K, N, V, etc
	 * 	- injected_noise_power: 
	 */


	/// starting with all zero message bits 
	std::vector<int> all_zero_m_bits(K, 0);

	FeedForwardTrellis encodingTrellis(code.k, code.n, code.v, code.numerators);

	std::vector<int> message_with_crc = generateRandomCRCMessage( code, NOISELESS );
	std::vector<int> transmitted_message = generateTransmittedMessage( message_with_crc, encodingTrellis );
	
	std::cout << "transmitted message size: " << transmitted_message.size() << std::endl;

	/// inject random noise with fixed energy (sum of squares)
	for ( int i = 0; i < 10; i++ ) {
		std::vector<double> standard_noise = awgn::generateStandardNormalNoise( (N/K) * NUM_INFO_BITS );
		std::vector<double> scaled_noise 	 = awgn::scaleNoise(standard_noise, injected_power);
	}

	// std::vector<double> noisy_message = 
	

	// hard-decoding

	// comparing with soft-decoding
	
}

// this generates a random binary string of length code.numInfoBits, and appends the appropriate CRC bits
std::vector<int> generateRandomCRCMessage(CodeInformation code, bool noiseless){
	
	std::vector<int> message ( code.numInfoBits );

	if (noiseless == 0) {
		for(int i = 0; i < code.numInfoBits; i++)
			message[i] = (rand()%2);
	} else {
		std::fill( message.begin(), message.end(), 0 ); // noiseless
	}
	// compute the CRC
	crc::crc_calculation(message, code.crcDeg, code.crc);
	return message;
}

// this takes the message bits, including the CRC, and encodes them using the trellis
std::vector<int> generateTransmittedMessage(std::vector<int> originalMessage, FeedForwardTrellis encodingTrellis){
	// encode the message
	std::vector<int> encodedMessage = encodingTrellis.encode(originalMessage);

	return encodedMessage;
}

// this takes the transmitted message and adds AWGN noise to it
// it also punctures the bits that are not used in the trellis
std::vector<double> addAWNGNoiseAndPuncture(std::vector<int> transmittedMessage, std::vector<int> puncturedIndices, double snr, bool noiseless){
	std::vector<double> receivedMessage;
	if(noiseless){
		for(int i = 0; i < transmittedMessage.size(); i++)
			receivedMessage.push_back((double)transmittedMessage[i]);
	} else {
		receivedMessage = awgn::addNoise(transmittedMessage, snr);
	}

	// puncture the bits. it is more convenient to puncture on this side than on the 
	// decoder, so we insert zeros which provide no information to the decoder
	for(int index = 0; index < puncturedIndices.size(); index++) {
		if (puncturedIndices[index] > receivedMessage.size()) {
			std::cout << "out of bounds index: " << puncturedIndices[index] << std::endl;
			std::cerr << "Puncturing index out of bounds" << std::endl;
			exit(1);
		}
		receivedMessage[puncturedIndices[index]] = 0;
	}

	return receivedMessage;
}