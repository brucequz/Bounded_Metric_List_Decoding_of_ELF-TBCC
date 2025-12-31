#include "../consts.h"
#include "../include/feedForwardTrellis.h"
#include "../include/lowRateListDecoder.h"
#include "../include/namespace.h"
#include "../include/types.h"

#include <cassert>
#include <cstddef>
#include <cstdlib>
#include <exception>
#include <iostream>
#include <unordered_map>
#include <utility>
#include <vector>

std::vector<int> generateTransmittedMessage(std::vector<int> info_crc,
                                            const FeedForwardTrellis& encodingTrellis);
std::vector<int> find_positive_divisor(int num);
void find_gabriel_neighbors(const CodeInformation& code);
std::vector<std::vector<int>> find_lyndon_words(int s, size_t n);
void search_message_gabriel_neighbors(const FeedForwardTrellis& trellis, const LowRateListDecoder& decoder);

int main() {

  /* Code config */
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
  // MessageInformation decodingResult;
  // decodingResult = listDecoder.forceDecoding_MaxListsize(receivedMessage, PUNCTURING_INDICES,
  //                                                        transmittedMessage, 866);
  // utils::print_int_vector(decodingResult.message);

  /* Find Gabriel Neighbors */
  // find_gabriel_neighbors(code);

  /* Search for Gabriel neighbors by using Lyndon messages */
  search_message_gabriel_neighbors(encodingTrellis, listDecoder);

}

std::vector<int> find_positive_divisor(int num) {
  if (num <= 0) {
    std::cerr << "num needs to be positive !" << std::endl;
  }
  std::vector<int> positive_divisors;
  for (int i = 1; i <= num; ++i){
    if (num % i == 0)
        positive_divisors.push_back(i);
  }
  return positive_divisors;
}

std::vector<std::vector<int>> find_lyndon_words(int s, size_t n) {
  /** Generate nonempty Lyndon words of length <= n over an s-symbol alphabet.
  The words are generated in lexicographic order, using an algorithm from
  J.-P. Duval, Theor. Comput. Sci. 1988, doi:10.1016/0304-3975(88)90113-2.
  As shown by Berstel and Pocchiola, it takes constant average time
  per generated word.
  */
  std::vector<std::vector<int>> result;
  std::vector<int> w = {-1};
  while (!w.empty()) {
    w.back() += 1;
    if (w.size() == n) {
      result.push_back(w);
    }
    size_t m = w.size();
    while (w.size() < n) {
      w.push_back(w[w.size()-m]);
    }
    while (!w.empty() && w.back() == s-1) {
      w.pop_back();
    }
  }
  return result;
}

void search_message_gabriel_neighbors(const FeedForwardTrellis& trellis, const LowRateListDecoder& decoder) {

  /* - Compute Lyndon words */
  std::vector<std::vector<int>> lyndon_messages;

  /* Find positive divisors */
  std::vector<int> positive_divisors = find_positive_divisor(Kconv);
  for (int div : positive_divisors) {
    std::vector<std::vector<int>> truncated_lyndon_messages = find_lyndon_words(2, div);
    
    for (size_t i = 0; i < truncated_lyndon_messages.size(); i++) {
      std::vector<int> lyndon_message(Kconv);
      for (size_t j = 0; j < lyndon_message.size(); j++) {
        lyndon_message[j] = truncated_lyndon_messages[i][j % div];
      }
      lyndon_messages.push_back(lyndon_message);
    }
  }

  std::cout << "number of lyndon messages: " << lyndon_messages.size() << std::endl;

  /* Records messages that lead to weight 12 codewords */
  std::string LyndonMFile = "output/lyndon_message_cwd_12.txt";
  std::vector<std::vector<int>> weight_12_messages;

  /* Read in generator matrix */
  std::string filename = "params/K24N48_TBCC_GenMatrix.txt";
  std::vector<std::vector<int>> G_mat = io::read2DVectorFromFile(filename);

  /* Initialize neighbor spectra */
  //!< key = hamming weight; pair.first = num_nei; pair.second = num_nonnei
  std::unordered_map<int, std::pair<int, int>> neighbor_spectra;

  for (size_t i_lyndon_w = 0; i_lyndon_w < lyndon_messages.size(); i_lyndon_w++) {
    /* Encode Lyndon messages */
    std::vector<int> m = lyndon_messages[i_lyndon_w];
    if (m.size() == 1) {
      int val = m[0];
      m = std::vector<int>(Kconv, val);
    } 
    std::vector<int> modulated_codeword = generateTransmittedMessage(m, trellis);
    // std::cout << "printing modulated_codeword: "; utils::print_int_vector(modulated_codeword);

    /* Check if they are Gabriel neighbor of the all-zero codeword */
    bool is_neighbor = decoder.isGabrielNeighbor(G_mat, modulated_codeword);

    int hamming_distance = 0;
    for (size_t i = 0; i < modulated_codeword.size(); i++) {
      if (modulated_codeword[i] < 0) {
        hamming_distance++;
      }
    }

    if (hamming_distance == 12) {
      weight_12_messages.push_back(m);
    }

    if (neighbor_spectra.find(hamming_distance) != neighbor_spectra.end()) {
      //!< if hamming distance already discovered
      if (is_neighbor) {
        neighbor_spectra[hamming_distance].first++;
      } else {
        neighbor_spectra[hamming_distance].second++;
      }
    } else {
      //!< if hamming distance entry not discovered
      if (is_neighbor) {
        neighbor_spectra[hamming_distance].first = 1;
        neighbor_spectra[hamming_distance].second = 0;
      } else {
        neighbor_spectra[hamming_distance].first = 0;
        neighbor_spectra[hamming_distance].second = 1;
      }
    }
  }

  /* Output lyndon messages with codeword weight 12 to file */
  try {
    OutputFile LyndonMessageCwd12_File(LyndonMFile);
    LyndonMessageCwd12_File.write2DVector(weight_12_messages);
    weight_12_messages.clear();
  } catch (const std::exception& e) {
    std::cerr << "Error: " << e.what() << std::endl;
  }

  /* Output neighbor spectra to console */
  for (const auto& kvpair : neighbor_spectra) {
    std::cout << "Hamming weight: " << kvpair.first << "; " << "nei: " << kvpair.second.first << "; non-nei: " << kvpair.second.second << std::endl;
  }
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