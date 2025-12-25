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

std::vector<int> generateTransmittedMessage(std::vector<int> info_crc,
                                            const FeedForwardTrellis& encodingTrellis);
void find_gabriel_neighbors(const CodeInformation& code);

int main() {

  /* Code config */
  // CodeInformation code;
  CodeInformation code{.k = k,
                       .n = n,
                       .v = V,
                       .crcLen = M + 1,
                       .crc = CRC,
                       .numInfoBits = K,
                       .numerators = NUMERATORS};

  /* - esno_dB setup - */
  float offset = 10 * log10((float) K / N);
  float esno_dB = EBN0 + offset;

  /* - Trellis setup - */
  FeedForwardTrellis encodingTrellis(code.k, code.n, code.v, code.numerators);

  /* - Decoder setup - */
  LowRateListDecoder listDecoder(encodingTrellis, MAX_LISTSIZE, code.crcLen, code.crc,
                                 STOPPING_RULE);

  /* ==== SIMULATION begins ==== */
  std::cout << std::endl
            << "**- Simulation Started for EbN0 = " << std::fixed << std::setprecision(2) << EBN0
            << " -**" << std::endl;

  /* Generate info */
  std::vector<int> info(K, 0);

  /* CRC encoder */
  std::vector<int> info_crc(Ncrc, 0);
  for (size_t i = 0; i < info.size(); i++) {
    info_crc[i] = info[i];
  }
  unsigned long remainder = crc::remdr(info_crc, CRC, Ncrc, M);
  for (size_t i = 0; i < M; i++) {
    std::cout << "remainder " << ((remainder >> (M - 1 - i)) & 1) << std::endl;
    info_crc[Kcrc + i] = (remainder >> (M - 1 - i)) & 1;
  }

  /* Convolutional encoder */
  std::vector<int> transmittedMessage = generateTransmittedMessage(info_crc, encodingTrellis);

  /* Convolutional code generator matrix */
  encodingTrellis.computeGeneratorMatrix();

  /* LLR generation */
  std::vector<float> receivedMessage =
      awgn::addAWNGNoise(transmittedMessage, PUNCTURING_INDICES, esno_dB, NOISELESS);

  std::cout << "info: ";
  utils::print_int_vector(info);
  std::cout << "info_crc: ";
  utils::print_int_vector(info_crc);
  std::cout << "remainder at encoding: " << remainder << std::endl;
  std::cout << "transmitted codeword: ";
  utils::print_int_vector(transmittedMessage);
  std::cout << "received sequence: ";
  utils::print_float_vector(receivedMessage);

  // /* Channel */
  // float esno = pow(10.0, esno_dB / 10.0);
  // for (size_t i = 0; i < receivedMessage.size(); i++) {
  //   receivedMessage[i] = receivedMessage[i] / (4 * esno);
  // }

  /* List decoding */
  MessageInformation decodingResult;
  decodingResult = listDecoder.forceDecoding_MaxListsize(receivedMessage, PUNCTURING_INDICES,
                                                         transmittedMessage, 1000);
  utils::print_int_vector(decodingResult.message);

  /* Find Gabriel Neighbors */
  find_gabriel_neighbors(code);
}

void find_gabriel_neighbors(const CodeInformation& code) {
  srand(0);

  if ((code.numInfoBits + code.crcLen - 1) % code.k != 0) {
    std::cout << "invalid msg + crc length" << std::endl;
    return;
  }

  FeedForwardTrellis encodingTrellis(code.k, code.n, code.v, code.numerators);

  int listSize = 1e4; // normal Viterbi
  LowRateListDecoder listDecoder(encodingTrellis, listSize, code.crcLen, code.crc, STOPPING_RULE);
  std::cout << "beginning simulations" << std::endl;
  // simulate the comms system
  std::string filename = "params/K24N48_TBCC_GenMatrix.txt";
  std::vector<std::vector<int>> G_mat = io::read2DVectorFromFile(filename);

  std::vector<float> zeroscwdb(Nconv, 1.0f);
  MessageInformation decodedMessage = listDecoder.lowrateDecoding_neighbors(zeroscwdb, G_mat);
}

// this takes the message bits, including the CRC, and encodes them using the trellis
std::vector<int> generateTransmittedMessage(std::vector<int> info_crc,
                                            const FeedForwardTrellis& encodingTrellis) {
  /*
  encodes to get the transmitted message bits (info + zero termination + crc) before modulation.
  */
  if (ENCODING_RULE != 'T' && ENCODING_RULE != 'Z') {
    std::cerr << "ISSUE: INCORRECT ENCODING_RULE" << std::endl;
  }
  // encode the message
  std::vector<int> encodedMessage;
  if (ENCODING_RULE == 'T') {
    encodedMessage = encodingTrellis.encode(info_crc);
    assert(encodedMessage.size() == (K + M) / k * n);
  } else if (ENCODING_RULE == 'Z') {
    for (int i = 0; i < V; i++) {
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