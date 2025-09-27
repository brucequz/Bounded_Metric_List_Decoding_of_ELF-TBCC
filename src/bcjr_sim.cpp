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


int main() {
  CodeInformation code;
  code.k = k;         // numerator of the rate
  code.n = n;         // denominator of the rate
  code.v = V;         // number of memory elements
  code.crcDeg = M+1;  // m+1, degree of CRC, # bits of CRC polynomial
  code.crc = CRC;     // CRC polynomial
  code.numInfoBits = K; // number of information bits
  code.numerators = NUMERATORS;



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
  std::vector<float> receivedMessage = awgn::addAWNGNoise(transmittedMessage, puncturedIndices, esno_dB, NOISELESS);
  // std::cout << "received message size: " << receivedMessage.size() << std::endl;
  // std::cout << "received message: ";
  // utils::print_float_vector(receivedMessage);
  
  // Produce Raw Channel values
  float esno = pow(10.0, esno_dB / 10.0);
  for (size_t i = 0; i < receivedMessage.size(); i++) {
    receivedMessage[i] = receivedMessage[i] / (4 * esno);
  }
}


