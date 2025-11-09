#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <cstdlib>
#include <cassert>
#include <numeric>
#include <string>
#include <sstream>
#include "mpi.h"

#include "../consts.h"
#include "../include/types.h"
#include "../include/namespace.h"
#include "../include/feedForwardTrellis.h"
#include "../include/lowRateListDecoder.h"


std::vector<int> generateTransmittedMessage(std::vector<int> info_crc, FeedForwardTrellis encodingTrellis);

int main() {

  float EbN0 = EBN0[0];
  std::ostringstream ebn0_str;
  ebn0_str.precision(2);
  ebn0_str << std::fixed << EbN0;


  CodeInformation code;
  code.k = k;         // numerator of the rate
  code.n = n;         // denominator of the rate
  code.v = V;         // number of memory elements
  code.crcDeg = M+1;  // m+1, degree of CRC, # bits of CRC polynomial
  code.crc = CRC;     // CRC polynomial
  code.numInfoBits = K; // number of information bits
  code.numerators = NUMERATORS;

  float esno_dB = 0.0;
  float offset = 10 * log10((float)K/N);
  esno_dB = EbN0 + offset;

  /* - Trellis setup - */
  FeedForwardTrellis encodingTrellis(code.k, code.n, code.v, code.numerators);

  /* - Decoder setup - */
  LowRateListDecoder listDecoder(encodingTrellis, MAX_LISTSIZE, code.crcDeg, code.crc, STOPPING_RULE);



  std::vector<int> originalMessage = crc::generateRandomCRCMessage(code);
  // std::cout << "original message: ";
  // utils::print_int_vector(originalMessage);
  // std::cout << std::endl;
  std::cout << "printing original message size: " << originalMessage.size() << std::endl; 
  std::vector<int> transmittedMessage = generateTransmittedMessage(originalMessage, encodingTrellis);
  // std::cout << "transmitted message: ";
  // utils::print_int_vector(transmittedMessage);
  std::cout << ", length = " << transmittedMessage.size() << std::endl;
  std::vector<float> receivedMessage = awgn::addAWNGNoise(transmittedMessage, PUNCTURING_INDICES, esno_dB, NOISELESS);
  // std::cout << "received message size: " << receivedMessage.size() << std::endl;
  // std::cout << "received message: ";
  // utils::print_float_vector(receivedMessage);
  
  // Produce Raw Channel values
  // float esno = pow(10.0, esno_dB / 10.0);
  // for (size_t i = 0; i < receivedMessage.size(); i++) {
  //   receivedMessage[i] = receivedMessage[i] / (4 * esno);
  // }

  std::vector<float> test_y = {0.589739,-0.284312,-1.08152,-0.0696566,-0.134552,-2.526,1.72583,1.17851,0.952543,-1.21649,-1.29977,2.2999,-1.43011,1.11148,-1.28908,-0.689244,2.09608,1.19325,1.30767,-0.570555,-0.210941,1.81212,-0.463859,-1.13916,-1.53398,-0.905185,1.06429,1.34352,0.946723,1.00888,1.06857,1.53815,-1.69318,-0.606677,0.327464,-0.806208,0.26177,1.08627,1.11096,-0.195175,0.535952,1.1712,-2.68258,2.68117,0.739682,-1.95271,0.909489,0.867118,2.08665,-1.7794,0.615295,-2.01334,1.56331,1.76181,1.75014,1.96397,-0.678962,1.21492,-1.92767,0.19834,0.71163,1.43103,-0.167239,-2.01844,0.486671,-0.845129,2.02314,-1.52479,-2.65461,-0.251462,1.37641,0.756147,1.10593,-2.13908,2.20357,0.294537,1.15499,-2.24204,-1.36013,2.05836,-1.58403,-0.800992,-0.0780751,2.03898,1.34308,-1.44849,-0.186794,1.4361,1.61541,1.16107,-0.143325,1.35607,1.16623,-1.88832,0.272814,0.164786,1.87965,1.13101,0.973139,-0.10615,-1.28584,-1.50928,-0.145866,-1.17508,-0.464871,-1.23001,1.78707,-0.0238718,1.52092,-0.68787,-1.00671,1.7719,-0.328828,-1.24241,-2.25049,-1.28118,1.68136,-1.11321,0.7377,1.58502,0.0245299,-0.102473,-0.386636,-0.0467234,-1.2811,-0.877867,-2.03576,2.01943,-0.70135,1.41664,0.727573,-0.711986,0.0230203,1.77984,-0.407218,-2.81605,-1.12916,-1.50594,-1.61196,0.457685,2.01931,-0.967515,-0.597745,1.29109,0.887739,-1.08186,-1.13637,0.975162,1.02225,0.649267,-3.00996,1.98873,-0.570733,-1.25564,1.4082,-1.35827,-1.73009,-1.91991,0.681697,-1.73102,1.28938,2.5335,2.66804,-1.15876,0.564328,0.652644,-0.260803,1.3354,0.888268,-2.16116,-1.83161,-0.655841};

  listDecoder._bcjr_log_gamma(test_y, 0.6735);
}



// this takes the message bits, including the CRC, and encodes them using the trellis
std::vector<int> generateTransmittedMessage(std::vector<int> info_crc, FeedForwardTrellis encodingTrellis){
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