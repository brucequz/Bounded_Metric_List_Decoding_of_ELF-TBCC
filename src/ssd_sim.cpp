#include "../include/consts_K16N32.h"

#include <cassert>
#include <cstddef>
#include <cstdlib>
#include <iostream>
#include <vector>

#include "../consts.h"
#include "../include/types.h"
#include "../include/namespace.h"
#include "../include/feedForwardTrellis.h"
#include "../include/lowRateListDecoder.h"

std::vector<int> generateTransmittedMessage(std::vector<int> info_crc, const FeedForwardTrellis& encodingTrellis);

int main()
{

  /* Code config */
  CodeInformation code;
  code.k = k;         // numerator of the rate
  code.n = n;         // denominator of the rate
  code.v = V;         // number of memory elements
  code.crcDeg = M+1;  // m+1, degree of CRC, # bits of CRC polynomial
  code.crc = CRC;     // CRC polynomial
  code.numInfoBits = K; // number of information bits
  code.numerators = NUMERATORS;

  /* - Simulation esno_dB setup - */
  std::vector<int> puncturedIndices = PUNCTURING_INDICES;
  float esno_dB = 0.0;
  float offset = 10 * log10((float)K/N);
  esno_dB = EBN0 + offset;

  /* - Trellis setup - */
  FeedForwardTrellis encodingTrellis(code.k, code.n, code.v, code.numerators);

  /* - Decoder setup - */
  LowRateListDecoder listDecoder(encodingTrellis, MAX_LISTSIZE, code.crcDeg, code.crc, STOPPING_RULE);
  
  /* ==== SIMULATION begins ==== */
  std::cout << std::endl << "**- Simulation Started for EbN0 = " << std::fixed << std::setprecision(2) << EBN0 << " -**" << std::endl;
  int num_trials	 	= 0;
  
  /* Generate info with crc */
  std::vector<int> originalMessage = crc::generateRandomCRCMessage(code);
  utils::print_int_vector(originalMessage);

  /* Convolutional encoder */
  std::vector<int> transmittedMessage = generateTransmittedMessage(originalMessage, encodingTrellis);

  /* LLR generation */
  std::vector<float> receivedMessage = awgn::addAWNGNoise(transmittedMessage, puncturedIndices, esno_dB, NOISELESS);
  utils::print_float_vector(receivedMessage);

  // /* Channel */
  // float esno = pow(10.0, esno_dB / 10.0);
  // for (size_t i = 0; i < receivedMessage.size(); i++) {
  //   receivedMessage[i] = receivedMessage[i] / (4 * esno);
  // }
  

}

// this takes the message bits, including the CRC, and encodes them using the trellis
std::vector<int> generateTransmittedMessage(std::vector<int> info_crc, const FeedForwardTrellis& encodingTrellis){
	/*
	encodes to get the transmitted message bits (info + zero termination + crc) before modulation.
	*/ 
	if (ENCODING_RULE != 'T' && ENCODING_RULE != 'Z') {std::cerr << "ISSUE: INCORRECT ENCODING_RULE" << std::endl;}
	// encode the message
	std::vector<int> encodedMessage;
	if (ENCODING_RULE == 'T') {
		encodedMessage = encodingTrellis.encode(info_crc);
		assert(encodedMessage.size() == (K+M) / k * n);
	} else if (ENCODING_RULE == 'Z') {
		for (int i=0; i<V; i++){
			info_crc.push_back(0);
		}
		// std::cout << "info crc with termination: ";
		// utils::print_int_vector(info_crc);
		// std::cout << std::endl;
		encodedMessage = encodingTrellis.encode_zt(info_crc);
		assert(encodedMessage.size() == (K+M+V) / k * n);
		
	}
	return encodedMessage;
}