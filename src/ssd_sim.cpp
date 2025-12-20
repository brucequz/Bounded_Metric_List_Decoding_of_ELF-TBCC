#include "../consts.h"
#include "../include/feedForwardTrellis.h"
#include "../include/lowRateListDecoder.h"
#include "../include/namespace.h"
#include "../include/types.h"

#include <cassert>
#include <cstddef>
#include <cstdlib>
#include <iostream>
#include <vector>

#include "../include/namespace.h"

std::vector<int> generateTransmittedMessage(std::vector<int> info_crc,
                                            const FeedForwardTrellis& encodingTrellis);

int main()
{

  /* Code config */
  // CodeInformation code;
  CodeInformation code{.k = k,
                       .n = n,
                       .v = V,
                       .crcDeg = M + 1,
                       .crc = CRC,
                       .numInfoBits = K,
                       .numerators = NUMERATORS};

  /* - esno_dB setup - */
  float offset = 10 * log10((float)K / N);
  float esno_dB = EBN0 + offset;

  /* - Trellis setup - */
  FeedForwardTrellis encodingTrellis(code.k, code.n, code.v, code.numerators);

  /* - Decoder setup - */
  LowRateListDecoder listDecoder(encodingTrellis, MAX_LISTSIZE, code.crcDeg, code.crc,
                                 STOPPING_RULE);

  /* ==== SIMULATION begins ==== */
  std::cout << std::endl
            << "**- Simulation Started for EbN0 = " << std::fixed << std::setprecision(2) << EBN0
            << " -**" << std::endl;
  int num_trials = 0;

  /* Generate info */
  std::vector<int> info(K, 0);

  /* CRC encoder */
  std::vector<int> info_crc(Ncrc, 0);
  unsigned long remainder = crc::remdr(info_crc, CRC, Ncrc, M);
  std::cout << "remainder at encoding: " << remainder << std::endl;
  
  /* Convolutional encoder */
  std::vector<int> transmittedMessage = generateTransmittedMessage(info_crc, encodingTrellis);
  std::cout << "printing transmitted message: " << std::endl;
  utils::print_int_vector(transmittedMessage);
  std::cout << std::endl;

  /* LLR generation */
  std::vector<float> receivedMessage =
      awgn::addAWNGNoise(transmittedMessage, PUNCTURING_INDICES, esno_dB, NOISELESS);
  utils::print_float_vector(receivedMessage);

  // /* Channel */
  // float esno = pow(10.0, esno_dB / 10.0);
  // for (size_t i = 0; i < receivedMessage.size(); i++) {
  //   receivedMessage[i] = receivedMessage[i] / (4 * esno);
  // }

  /* List decoding */
  MessageInformation decodingResult;
  decodingResult = listDecoder.forceDecoding_MaxListsize(receivedMessage, PUNCTURING_INDICES, transmittedMessage, 6);
  utils::print_int_vector(decodingResult.message);

  std::cout << "decoding result metric: " << decodingResult.metric << std::endl;

  
}

// this takes the message bits, including the CRC, and encodes them using the trellis
std::vector<int> generateTransmittedMessage(std::vector<int> info_crc,
                                            const FeedForwardTrellis& encodingTrellis)
{
  /*
  encodes to get the transmitted message bits (info + zero termination + crc) before modulation.
  */
  if (ENCODING_RULE != 'T' && ENCODING_RULE != 'Z')
  {
    std::cerr << "ISSUE: INCORRECT ENCODING_RULE" << std::endl;
  }
  // encode the message
  std::vector<int> encodedMessage;
  if (ENCODING_RULE == 'T')
  {
    encodedMessage = encodingTrellis.encode(info_crc);
    assert(encodedMessage.size() == (K + M) / k * n);
  }
  else if (ENCODING_RULE == 'Z')
  {
    for (int i = 0; i < V; i++)
    {
      info_crc.push_back(0);
    }
    // std::cout << "info crc with termination: ";
    // utils::print_int_vector(info_crc);
    // std::cout << std::endl;
    encodedMessage = encodingTrellis.encode_zt(info_crc);
    assert(encodedMessage.size() == (K + M + V) / k * n);
  }
  return encodedMessage;
}