#include "../consts.h"
#include "../include/feedForwardTrellis.h"
#include "../include/lowRateListDecoder.h"
#include "../include/namespace.h"
#include "../include/types.h"
#include "../include/json.hpp"
#include "../include/fileHandler.h"

#include <cassert>
#include <cmath>
#include <cstddef>
#include <cstdlib>
#include <exception>
#include <fmt/core.h>
#include <fmt/os.h>
#include <fmt/ranges.h>
#include <iostream>
#include <unordered_map>
#include <map>
#include <random>
#include <utility>
#include <vector>

using json = nlohmann::json;

std::vector<int> generateTransmittedMessage(std::vector<int> info_crc,
                                            const FeedForwardTrellis& encodingTrellis);
std::vector<int> find_positive_divisor(int num);
void find_gabriel_neighbors(const CodeInformation& code);
std::vector<std::vector<int>> find_lyndon_words(int s, size_t n);
void search_message_gabriel_neighbors(const FeedForwardTrellis& trellis,
                                      const LowRateListDecoder& decoder);
void compute_serial_coded_spectra(const FeedForwardTrellis& trellis);
void FER_simulation(const FeedForwardTrellis& trellis,
                    LowRateListDecoder decoder);

int main() {

  /* - Generator setup - */
  std::mt19937 gen(BASE_SEED);
  std::bernoulli_distribution d(0.5);

  /* - Code config - */
  CodeInformation code{.kconv = kconv,
                       .nconv = nconv,
                       .v = V,
                       .crcLen = M + 1,
                       .crc = CRC,
                       .numInfoBits = K,
                       .numerators = NUMERATORS};

  /* - Trellis setup - */
  FeedForwardTrellis encodingTrellis(code.kconv, code.nconv, code.v, code.numerators);

  /* - Decoder setup - */
  LowRateListDecoder listDecoder(&encodingTrellis, MAX_LISTSIZE, code.crcLen, code.crc,
                                 STOPPING_RULE);


  std::map<std::pair<int, int>, std::vector<std::vector<int>>> cosetLeadersMsgs;

  std::ifstream Msg_inFile(OUTPUTFILEPATH + "nonTB_coset_leaders_msgs.json");
  if (!Msg_inFile) {
      fmt::print(stderr, "Error: Could not open file for reading!\n");
      return 1;
  }

  // 2. Parse JSON
  json j_msg_read;
  Msg_inFile >> j_msg_read;

  // 3. Convert JSON back to the complex map
  for (auto const& [strKey, value] : j_msg_read.items()) {
      // Find the comma we used as a delimiter
      size_t commaPos = strKey.find(',');
      if (commaPos != std::string::npos) {
          // Extract the two integers from the string "10,20"
          int x = std::stoi(strKey.substr(0, commaPos));
          int y = std::stoi(strKey.substr(commaPos + 1));

          // Convert the JSON array of arrays back to vector<vector<int>>
          cosetLeadersMsgs[{x, y}] = value.get<std::vector<std::vector<int>>>();
      }
  }

  std::map<std::pair<int, int>, std::vector<std::vector<int>>> cosetLeadersCwds;

  std::ifstream Cwd_inFile(OUTPUTFILEPATH + "nonTB_coset_leaders_cwds.json");
  if (!Cwd_inFile) {
      fmt::print(stderr, "Error: Could not open file for reading!\n");
      return 1;
  }

  // 2. Parse JSON
  json j_cwd_read;
  Cwd_inFile >> j_cwd_read;

  // 3. Convert JSON back to the complex map
  for (auto const& [strKey, value] : j_cwd_read.items()) {
      // Find the comma we used as a delimiter
      size_t commaPos = strKey.find(',');
      if (commaPos != std::string::npos) {
          // Extract the two integers from the string "10,20"
          int x = std::stoi(strKey.substr(0, commaPos));
          int y = std::stoi(strKey.substr(commaPos + 1));

          // Convert the JSON array of arrays back to vector<vector<int>>
          cosetLeadersCwds[{x, y}] = value.get<std::vector<std::vector<int>>>();
      }
  }

  std::vector<std::vector<int>> gabrielNeighbors;
  try {
    std::vector<int> neightborList = {0, 8};
    for (int hamming_dist : neightborList) {
      std::string gabrielNeighborFileName = OUTPUTFILEPATH + fmt::format("codeword_{}.txt", hamming_dist);
      FileHandler gabrielNeighborFile(gabrielNeighborFileName);
      std::vector<std::vector<int>> temp = gabrielNeighborFile.read2DVector<int>();
      gabrielNeighbors.insert(gabrielNeighbors.end(), temp.begin(), temp.end());
    }

  } catch (const std::exception& e) {
    std::cerr << "Error: " << e.what() << std::endl;
  }

  /* Simulation parameters */
  for (size_t i_eb = 0; i_eb < EBN0.size(); i_eb++) {
    int num_sim = 0;
    int num_ssd_err = 0;
    // int num_slvd_err = 0;

    while (num_ssd_err < MAX_ERRORS) {
      num_sim++;

      /* Generate info */
      std::vector<int> info(K, 0);
      for (size_t i = 0; i < info.size(); ++i) {
        int bit = d(gen);
        info[i] = bit;
      }

      /* CRC encoder */
      std::vector<int> info_crc(Ncrc, 0);
      for (size_t i = 0; i < info.size(); i++) {
        info_crc[i] = info[i];
      }
      std::vector<int> remainder = crc::remdr_slidingWindow(info_crc, CRC_VEC, false);
      for (size_t i = 0; i < remainder.size(); i++) {
        info_crc[Kcrc + i] = remainder[i];
      }

      /* Convolutional encoder */
      std::vector<int> modulatedCodeword = generateTransmittedMessage(info_crc, encodingTrellis);

      /* Channel */
      float offset = 10 * log10((float) K / N);
      float esno_dB = EBN0[i_eb] + offset;
      float esno = pow(10.0, esno_dB / 10.0);
      std::vector<float> receivedWord =
          awgn::addAWNGNoise(modulatedCodeword, PUNCTURING_INDICES, esno_dB, NOISELESS);
      for (size_t i = 0; i < receivedWord.size(); i++) {
        receivedWord[i] = receivedWord[i] / (4 * esno);
      }

      // float transmittedCorr = 0.0f;
      // for (size_t i = 0; i < modulatedCodeword.size(); i++) {
      //   transmittedCorr += modulatedCodeword[i] * receivedWord[i];
      // }
      // fmt::print("transmittedCorr: {}\n", transmittedCorr);
      
      // std::cout << "info: " << std::endl;
      // utils::print_int_vector(info);
      // std::cout << "info_crc: " << std::endl;
      // utils::print_int_vector(info_crc);
      // std::cout << "remainder: ";
      // utils::print_int_vector(remainder);
      // std::cout << "transmitted codeword: " << std::endl;
      // utils::print_int_vector(modulatedCodeword);
      // std::cout << "received sequence: " << std::endl;
      // utils::print_float_vector(receivedWord);

      MessageInformation decodingResult;
      decodingResult = listDecoder.ssdSLVDDecoding(receivedWord, modulatedCodeword, PUNCTURING_INDICES, cosetLeadersMsgs, cosetLeadersCwds, gabrielNeighbors);
      if (decodingResult.codeword != modulatedCodeword) {
        num_ssd_err++;
      }

      /* SLVD decoding */
      // MessageInformation SLVDResult;
      // SLVDResult = listDecoder.lowRateDecoding_MaxListsize(receivedWord, PUNCTURING_INDICES);
      // if (SLVDResult.message != info_crc) {
      //   num_slvd_err++;
      // }
    }
    fmt::print("*-----------------\n");
    fmt::print("ebno_dB: {}\n", EBN0[i_eb]);
    fmt::print("num_sims: {}\n", num_sim);
    fmt::print("num_ssd_err: {}\n", num_ssd_err);
    // fmt::print("num_slvd_err: {}\n", num_slvd_err);
  }

  /* Find Coset Leaders */
  // listDecoder.findNonTBCosetLeaders(receivedWord, PUNCTURING_INDICES, std::pow(2, 15));

  /* FER simulation */
  // FER_simulation(encodingTrellis, listDecoder);

  // MessageInformation decodingResult;
  // decodingResult = listDecoder.forceDecoding_MaxListsize(receivedWord, PUNCTURING_INDICES,
  //                                                        modulatedCodeword, std::pow(2,15));
  // utils::print_int_vector(decodingResult.message);

  /* Compute distance spectra for serially-concatenated code */
  // compute_serial_coded_spectra(encodingTrellis);

  /* Find Gabriel Neighbors */
  // find_gabriel_neighbors(code);

  /* Search for Gabriel neighbors by using Lyndon messages */
  // search_message_gabriel_neighbors(encodingTrellis, listDecoder);

  return 0;
}

std::vector<int> find_positive_divisor(int num) {
  if (num <= 0) {
    std::cerr << "num needs to be positive !" << std::endl;
  }
  std::vector<int> positive_divisors;
  for (int i = 1; i <= num; ++i) {
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
      w.push_back(w[w.size() - m]);
    }
    while (!w.empty() && w.back() == s - 1) {
      w.pop_back();
    }
  }
  return result;
}

void search_message_gabriel_neighbors(const FeedForwardTrellis& trellis,
                                      const LowRateListDecoder& decoder) {

  /* - Compute Lyndon words */
  std::vector<std::vector<int>> lyndon_messages;

  /* Find positive divisors */
  std::vector<int> positive_divisors = find_positive_divisor(Kconv);
  for (int div : positive_divisors) {
    std::vector<std::vector<int>> truncated_lyndon_messages = find_lyndon_words(2, div);
    std::cout << "Lyndon word length: " << div
              << "; number of truncated lyndon messages: " << truncated_lyndon_messages.size()
              << std::endl;
    for (size_t i = 0; i < truncated_lyndon_messages.size(); i++) {
      std::vector<int> lyndon_message(Kconv);
      for (size_t j = 0; j < lyndon_message.size(); j++) {
        lyndon_message[j] = truncated_lyndon_messages[i][j % div];
      }
      lyndon_messages.push_back(lyndon_message);
    }
  }

  /* - Records messages that lead to weight 12 codewords */
  int hamming_weight_of_interest = 8;
  std::string LyndonMFile =
      OUTPUTFILEPATH + fmt::format("lyndon_message_cwd_{}.txt", hamming_weight_of_interest);
  std::vector<std::vector<int>> weight_12_messages;

  /* - Read in generator matrix */
  std::vector<std::vector<int>> G_mat = io::read2DVectorFromFile(GENMATRIXFILEPATH);

  /* Initialize neighbor spectra */
  //!< key = hamming weight; pair.first = num_nei; pair.second = num_nonnei
  std::unordered_map<int, std::pair<int, int>> neighbor_spectra;

  /* Initialize codeword distance spectra */
  std::map<int, int> cwd_dist_spectra;

  for (size_t i_lyndon_w = 0; i_lyndon_w < lyndon_messages.size(); i_lyndon_w++) {
    /* Encode Lyndon messages */
    std::vector<int> m = lyndon_messages[i_lyndon_w];

    std::vector<int> modulated_codeword = generateTransmittedMessage(m, trellis);

    /* Check if they are Gabriel neighbor of the all-zero codeword */
    bool is_neighbor = decoder.isGabrielNeighbor(G_mat, modulated_codeword);

    int hamming_distance = 0;
    for (size_t i = 0; i < modulated_codeword.size(); i++) {
      if (modulated_codeword[i] < 0) {
        hamming_distance++;
      }
    }

    if (hamming_distance == hamming_weight_of_interest) {
      weight_12_messages.push_back(m);
      // std::cout << "printing modulated_codeword: "; utils::print_int_vector(modulated_codeword);
    }

    /* cwd dist spectra */
    if (cwd_dist_spectra.find(hamming_distance) != cwd_dist_spectra.end()) {
      //!< if hamming distance already discovered
      cwd_dist_spectra[hamming_distance]++;
    } else {
      //!< if hamming distance entry not discovered
      cwd_dist_spectra[hamming_distance] = 1;
    }

    /* neighbor spectra */
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
    FileHandler LyndonMessageCwd12_File(LyndonMFile);
    LyndonMessageCwd12_File.write2DVector(weight_12_messages);
    weight_12_messages.clear();
  } catch (const std::exception& e) {
    std::cerr << "Error: " << e.what() << std::endl;
  }

  /* Output neighbor spectra to console */
  for (const auto& kvpair : neighbor_spectra) {
    std::cout << "Hamming weight: " << kvpair.first << "; " << "nei: " << kvpair.second.first
              << "; non-nei: " << kvpair.second.second << std::endl;
  }

  std::cout << "printing cwd dist spectra: " << std::endl;
  /* Output cwd dist spectra to console */
  for (const auto& kvpair : cwd_dist_spectra) {
    std::cout << "Hamming weight: " << kvpair.first << "; " << "number of cwds: " << kvpair.second << std::endl;
  }
}

void find_gabriel_neighbors(const CodeInformation& code) {
  srand(0);

  if ((code.numInfoBits + code.crcLen - 1) % code.kconv != 0) {
    std::cout << "invalid msg + crc length" << std::endl;
    return;
  }

  FeedForwardTrellis encodingTrellis(code.kconv, code.nconv, code.v, code.numerators);

  int listSize = 1e4; // normal Viterbi
  LowRateListDecoder listDecoder(&encodingTrellis, listSize, code.crcLen, code.crc, STOPPING_RULE);
  std::cout << "beginning simulations" << std::endl;
  // simulate the comms system
  std::vector<std::vector<int>> G_mat = io::read2DVectorFromFile(GENMATRIXFILEPATH);

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
    assert(encodedMessage.size() == (K + M) / kconv * nconv);
  } else if (ENCODING_RULE == 'Z') {
    for (int i = 0; i < V; i++) {
      info_crc.push_back(0);
    }
    // std::cout << "info crc with termination: ";
    // utils::print_int_vector(info_crc);
    // std::cout << std::endl;
    encodedMessage = encodingTrellis.encode_zt(info_crc);
    assert(encodedMessage.size() == (K + M + V) / kconv * nconv);
  }
  return encodedMessage;
}

void compute_serial_coded_spectra(const FeedForwardTrellis& trellis) {
  /* Initialize neighbor spectra */
  //!< key = hamming weight; pair.first = num_nei; pair.second = num_nonnei
  std::map<int, int> dist_spectra;

  std::vector<int> info;

  for (int i = 0; i < std::pow(2, K); i++) {
    info.clear();

    /* Generate info */
    for (int j = 0; j < K; j++) {
      info.push_back(((i >> j) & 1));
    }

    /* CRC encoder */
    std::vector<int> info_crc(Ncrc, 0);
    for (size_t i = 0; i < info.size(); i++) {
      info_crc[i] = info[i];
    }
    std::vector<int> remainder = crc::remdr_slidingWindow(info_crc, CRC_VEC, false);
    for (size_t i = 0; i < remainder.size(); i++) {
      info_crc[Kcrc + i] = remainder[i];
    }

    /* Convolutional encoder */
    std::vector<int> modulatedCodeword = generateTransmittedMessage(info_crc, trellis);

    /* Compute hamming weight of codeword */
    int hamming_weight = 0;
    for (size_t i_codeword = 0; i_codeword < modulatedCodeword.size(); i_codeword++) {
      if (modulatedCodeword[i_codeword] < 0) {
        hamming_weight++;
      }
    }

    /*record to distance spectrum */
    if (dist_spectra.find(hamming_weight) != dist_spectra.end()) {
      //!< if hamming distance already discovered
      dist_spectra[hamming_weight]++;
    } else {
      //!< if hamming distance entry not discovered
      dist_spectra[hamming_weight] = 1;
    }
  }

  /* Output neighbor spectra to console */
  std::cout << "CRC: ";
  utils::print_int_vector(CRC_VEC);
  for (const auto& kvpair : dist_spectra) {
    std::cout << "Hamming weight: " << kvpair.first << "; " << "num of codewords: " << kvpair.second
              << std::endl;
  }
}

void FER_simulation(const FeedForwardTrellis& trellis, LowRateListDecoder decoder) {

  /* - Generator setup - */
  std::mt19937 gen(BASE_SEED);
  std::bernoulli_distribution d(0.5);

  /* Reading in coset leaders */
  std::vector<std::vector<int>> cosetLeaderMsgs;
  std::vector<std::vector<int>> cosetLeaderCwds;
  std::vector<std::vector<int>> gabrielNeighbors;
  try {
    std::string cosetLeaderMsgName = OUTPUTFILEPATH + "coset_leaders_msgs.txt";
    FileHandler cosetLeaderMsgFile(cosetLeaderMsgName);
    cosetLeaderMsgs = cosetLeaderMsgFile.read2DVector<int>();

    std::string cosetLeaderCwdFileName = OUTPUTFILEPATH + "coset_leaders_cwds.txt";
    FileHandler cosetLeaderCwdFile(cosetLeaderCwdFileName);
    cosetLeaderCwds = cosetLeaderCwdFile.read2DVector<int>(' ');

    std::vector<int> neightborList = {0, 8, 10, 12, 14, 16};
    for (int hamming_dist : neightborList) {
      std::string gabrielNeighborFileName = OUTPUTFILEPATH + fmt::format("gabriel_codeword_{}.txt", hamming_dist);
      FileHandler gabrielNeighborFile(gabrielNeighborFileName);
      std::vector<std::vector<int>> temp = gabrielNeighborFile.read2DVector<int>();
      gabrielNeighbors.insert(gabrielNeighbors.end(), temp.begin(), temp.end());
    }

  } catch (const std::exception& e) {
    std::cerr << "Error: " << e.what() << std::endl;
  }

  // std::cout << "cosetLeaderCwds[10]: ";
  // utils::print_int_vector(cosetLeaderCwds[10]);
  std::cout << "gabrielNeighbor size: " << gabrielNeighbors.size() << std::endl;
  // std::cout << "gabrielNeighbor[0]: ";
  // utils::print_int_vector(gabrielNeighbors[0]);

  /* - esno_dB setup - */
  float offset = 10 * log10((float) K / N);
  
  for (size_t i_ebno = 0; i_ebno < EBN0.size(); i_ebno++) {
    float esno_dB = EBN0[i_ebno] + offset;
    std::cout << "esno_dB: " << std::fixed << std::setprecision(4) << esno_dB << std::endl;
  
  
    int num_ssd_slvd_error = 0;
    int num_slvd_error = 0;
    int num_sims = 0;

    while (num_slvd_error < MAX_ERRORS) {
      num_sims++;
      if (num_sims % LOGGING_ITERS == 0) {
        std::cout << "num_sims: " << num_sims << "; num_slvd_error: " << num_slvd_error << std::endl;
      }

      /* Generate info */
      std::vector<int> info(K, 0);
      for (size_t i = 0; i < info.size(); ++i) {
        int bit = d(gen);
        info[i] = bit;
      }

      /* CRC encoder */
      std::vector<int> info_crc(Ncrc, 0);
      for (size_t i = 0; i < info.size(); i++) {
        info_crc[i] = info[i];
      }
      std::vector<int> remainder = crc::remdr_slidingWindow(info_crc, CRC_VEC, false);
      for (size_t i = 0; i < remainder.size(); i++) {
        info_crc[Kcrc + i] = remainder[i];
      }

      /* Convolutional encoder */
      std::vector<int> modulatedCodeword = generateTransmittedMessage(info_crc, trellis);

      /* Channel */
      std::vector<float> receivedWord =
          awgn::addAWNGNoise(modulatedCodeword, PUNCTURING_INDICES, esno_dB, NOISELESS);
      float esno = pow(10.0, esno_dB / 10.0);
      for (size_t i = 0; i < receivedWord.size(); i++) {
        receivedWord[i] = receivedWord[i] / (4 * esno);
      }

      /* Correlation metric */
      // float corr = 0.0f;
      // for (size_t i = 0; i < receivedWord.size(); i++) {
      //   corr += modulatedCodeword[i] * receivedWord[i];
      // }
      // std::cout << "transmitted codeword corr: " << corr << std::endl;

      

      /* SSD-SLVD decoding*/
      // MessageInformation decodingResult;
      // decodingResult = decoder.ssdSLVDDecoding(receivedWord, PUNCTURING_INDICES, modulatedCodeword, cosetLeaderMsgs, cosetLeaderCwds, gabrielNeighbors);
      // if (decodingResult.codeword != modulatedCodeword) {
      //   num_ssd_slvd_error++;
      // }

      /* SLVD decoding */
      MessageInformation SLVDResult;
      SLVDResult = decoder.lowRateDecoding_MaxListsize(receivedWord, PUNCTURING_INDICES);
      if (SLVDResult.message != info_crc) {
        // std::cout << "SLVD failed! " << std::endl;
        num_slvd_error++;
      }
    } // while (num_ssd_slvd_error < MAX_ERRORS)

    std::cout << "num_sim: " << num_sims << std::endl;
    std::cout << "num_slvd_error: " << num_slvd_error << std::endl;
    std::cout << "num_ssd_slvd_error: " << num_ssd_slvd_error << std::endl;
  } // for (size_t i_ebno = 0; i_ebno < EBN0.size(); i_ebno++)
}