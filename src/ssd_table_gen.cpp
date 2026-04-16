#include <cmath>
#include <vector>
#include <iostream>
#include <random>
#include <time.h>
#include <fmt/os.h>
#include <fmt/ranges.h>

#include "../consts.h"
#include "../include/json.hpp"
#include "../include/lowRateListDecoder.h"
#include "../include/minHeap.h"
#include "../include/namespace.h"
#include "../include/feedForwardTrellis.h"

using json = nlohmann::json;
static std::random_device rd{};
static std::mt19937 generator{rd()};

struct codeInformation {
  int kconv;
  int nconv;
  int v;
  int crcLen;
  int crc;
  int numInfoBits;
  double dmin;
  std::vector<int> numerators;
  int denominator;
  std::vector<std::vector<int>> hMatrix;
  std::vector<int> punc_idx;
};

void lowrate_TBferPerformance(codeInformation code);

void LE_SLVD_updated(codeInformation code);
void generateNeighborList_ELF(codeInformation code);
void generateNeighborList_ESD0(codeInformation code);
void generateNeighborList_ESD_ED_unsorted(codeInformation code);

std::vector<int> generateRandomCRCMessage(codeInformation code);
std::vector<double> generateTransmittedMessage(std::vector<int> originalMessage,
                                               FeedForwardTrellis encodingTrellis,
                                               double snr,
                                               std::vector<int> puncturedIndices,
                                               bool noiseless = false);

// int main() {
//   std::cout << __cplusplus << std::endl;
//   codeInformation code;

//   // code.kconv = 1;						// numerator of the rate
//   // code.nconv = 2;						// denominator of the rate
//   // code.v = 8;						// number of memory elements
//   // code.crcLen = 13;				// degree of the crc (m + 1)
//   // code.crc = 0x1565;			// crc in decimal, convert from hex if needed
//   // code.numInfoBits = 64;

//   // // optimal code numerators and denominator are known, and are given in octal, we use octal to calculate nextStates
//   // code.numerators = {561, 753}; // dmin = 16 (hamming) *4
//   // PUNCTURING_INDICES = {4, 10, 21, 24, 31, 37, 42, 48, 59, 62, 69, 75, 80, 86, 97, 100, 107, 113, 118, 124, 135, 138, 145, 151};

//   /*
// 	code.kconv = 1;						// numerator of the rate
// 	code.nconv = 2;						// denominator of the rate
// 	code.v = 5;						// number of memory elements
//     code.crcLen = 6;				// degree of the crc (m + 1)
// 	code.crc = 0x3F;				// crc in decimal, convert from hex if needed
// 	code.numInfoBits = 64;

// 	// optimal code numerators and denominator are known, and are given in octal, we use octal to calculate nextStates
// 	code.numerators = {53, 75}; //
// 	PUNCTURING_INDICES = {};

// 	*/

//   code.kconv       = 1;     // numerator of the rate
//   code.nconv       = 2;     // denominator of the rate
//   code.v           = 6;     // number of memory elements
//   code.crcLen      = 9;     // degree of the crc (m + 1)
//   code.crc         = 0x123; // crc in decimal, convert from hex if needed
//   code.numInfoBits = 16;
//   // dmin = 8
//   // optimal code numerators and denominator are known, and are given in octal, we use octal to calculate nextStates
//   code.numerators = {133, 171}; //
//   // PUNCTURING_INDICES = {3, 4, 9, 10, 15, 16, 21, 22, 27, 28, 33, 34, 39, 40, 45, 46};  // std::vector<int> puncpat = {1, 1, 1, 0, 0 ,1};
//   // PUNCTURING_INDICES = {};  // std::vector<int> puncpat = {1, 1, 1, 0, 0 ,1};

//   /*
// 	code.kconv = 1;						// numerator of the rate
// 	code.nconv = 4;						// denominator of the rate
// 	code.v = 8;						// number of memory elements
//     code.crcLen = 13;				// degree of the crc (m + 1)
// 	code.crc = 0x1005;				// crc in decimal, convert from hex if needed
// 	code.numInfoBits = 64;
// 	// dmin = 35

// 	// optimal code numerators and denominator are known, and are given in octal, we use octal to calculate nextStates
// 	code.numerators = {533, 575, 647, 711}; //
// 	PUNCTURING_INDICES = {4, 10, 21, 24, 31, 37, 42, 48, 59, 62, 69, 75, 80, 86, 97, 100, 107, 113, 118, 124, 135, 138, 145, 151, 156, 162, 173, 176, 183, 189, 194, 200, 211, 214, 221, 227, 232, 238, 249, 252, 259, 265, 270, 276, 287, 290, 297, 303};
// 	*/

//   // code.kconv = 1;						// numerator of the rate
//   // code.nconv = 1;						// denominator of the rate
//   // code.v = 8;						// number of memory elements
//   // code.crcLen = 9;				// degree of the crc (m + 1)
//   // code.crc = 0x18B;				// crc in decimal, convert from hex if needed
//   // code.numInfoBits = 64;
//   // // dmin = 5

//   // // optimal code numerators and denominator are known, and are given in octal, we use octal to calculate nextStates
//   // code.numerators = {205}; //
//   // PUNCTURING_INDICES = {};

//   // generateNeighborList_ELF(code);
//   // generateNeighborList_ESD0(code);
//   // LE_SLVD_updated(code);
//   // generateNeighborList_ESD_ED_unsorted(code);
//   // generateNeighborList_sequential(code);
//   // generateNeighborList_sequential_ZT(code);

//   return 0;
// }

// void generateNeighborList_sequential_ZT(codeInformation code){
// 	std::cout << "generating the list of neighbors sequentially" << std::endl;

// 	// set random seed for message generation
// 	srand(time(NULL));

// 	// check to make sure the code has valid values
// 	if ((code.numInfoBits + code.crcLen - 1) % code.kconv != 0) {
// 		std::cout << "invalid msg + crc length" << std::endl;
// 		return;
// 	}

// 	// this will decide how long we want each list to be
// 	int listSize = 1e8;
// 	double thre = 24; // find neighbors up to a certain distance

// 	// the below are the relevant initializations for low rate
// 	// decoding, this will be altered for high rate or ZTCC, for example
// 	FeedForwardTrellis encodingTrellis(code.kconv, code.nconv, code.v, code.numerators, PUNCTURING_INDICES);
// 	LowRateListDecoder listDecoder(encodingTrellis, listSize, code.crcLen, code.crc, PUNCTURING_INDICES);

// 	// initializing the zeros message to find neighbors
// 	std::vector<double> allZerosMessage;
// 	int messageLength = (code.numInfoBits + code.crcLen + code.v - 1) * code.nconv / code.kconv; //ZT
// 	std::cout << "message length: " << messageLength << std::endl;
// 	for(int i = 0; i < messageLength; i++){
// 		allZerosMessage.push_back(1);
// 	}
// 	// puncture the bits. we insert zeros which provide no information to the decoder
// 	for(int index = 0; index < PUNCTURING_INDICES.size(); index++)
// 		allZerosMessage[PUNCTURING_INDICES[index]] = 0;

// 	std::string dir = "kconv" + std::to_string(code.kconv) + "nconv" + std::to_string(code.nconv) + "v" + std::to_string(code.v) +
// 					 "m" + std::to_string(code.crcLen - 1) + "K" + std::to_string(code.numInfoBits) + "_ZT/";
// 	std::cout << dir << std::endl;
// 	listDecoder.generateNeighborList_sequential_ZT(allZerosMessage, dir, thre);

// }

// void generateNeighborList_sequential(codeInformation code) {
//   std::cout << "generating the list of neighbors sequentially" << std::endl;

//   // set random seed for message generation
//   srand(0);

//   // check to make sure the code has valid values
//   if ((code.numInfoBits + code.crcLen - 1) % code.kconv != 0) {
//     std::cout << "invalid msg + crc length" << std::endl;
//     return;
//   }

//   // this will decide how long we want each list to be
//   int listSize = 1e8;
//   double thre  = 20; // find neighbors up to a certain distance

//   // the below are the relevant initializations for low rate
//   // decoding, this will be altered for high rate or ZTCC, for example
//   FeedForwardTrellis encodingTrellis(code.kconv, code.nconv, code.v, code.numerators);
//   LowRateListDecoder listDecoder(encodingTrellis, listSize, code.crcLen, code.crc, STOPPING_RULE);

//   /* create all-zero, bpsk-modulated info_crc */
//   std::vector<float> allZerosMessage;
//   int messageLength = (code.numInfoBits + code.crcLen - 1) * code.nconv / code.kconv;
//   for (int i = 0; i < messageLength; i++) { allZerosMessage.push_back(1); }

//   /* puncturing */
//   for (size_t index = 0; index < PUNCTURING_INDICES.size(); index++) { allZerosMessage[PUNCTURING_INDICES[index]] = 0; }
//   std::cout << "printing allzero message: " << std::endl;
//   utils::print_float_vector(allZerosMessage);

//   std::string dir = "kconv" + std::to_string(code.kconv) + "nconv" + std::to_string(code.nconv) + "v" +
//                     std::to_string(code.v) + "m" + std::to_string(code.crcLen - 1) + "K" +
//                     std::to_string(code.numInfoBits) + "/";
//   std::cout << dir << std::endl;
//   // listDecoder.generateNeighborList_sequential(allZerosMessage, dir, thre);
//   listDecoder.generateNeighborList_sequential_TBonly(allZerosMessage, dir, thre);
//   // listDecoder.readNeighborList("test/");
// }

void LowRateListDecoder::findTBOffsets() {
  /* Output file name & output data structure */
  std::string TBOffsetsCWDFileName = OUTPUTFILEPATH + "TB_Offsets_CWD.json";
  std::vector<std::map<int, std::vector<std::vector<int>>>> TB_Offsets_CWD((int)std::pow(2, M));
  json j_TB_Offsets_CWD = json::array();

  std::string TBOffsetsMSGFileName = OUTPUTFILEPATH + "TB_Offsets_MSG.json";
  std::vector<std::map<int, std::vector<std::vector<int>>>> TB_Offsets_MSG((int)std::pow(2, M));
  json j_TB_Offsets_MSG = json::array();
  std::vector<std::map<int, int>> TB_Offsets_info((int)std::pow(2, M));

  std::vector<int> all_zeros_codeword(N, 0);
  std::vector<float> all_ones_modulated_codeword(N, 1.0f);
  std::vector<std::vector<cell>> trellisInfo; // trellisInfo is indexed [state][stage]
  trellisInfo = constructLowRateTrellis_Punctured(all_ones_modulated_codeword, PUNCTURING_INDICES);

  // start search
  MinHeap detourTree;
  std::vector<std::vector<int>> previousPaths;

  // create nodes for each valid ending state with no detours
  for (int i = 0; i < lowrate_numStates; i++) {
    DetourObject detour;
    detour.startingState = i;
    detour.pathMetric    = trellisInfo[i][lowrate_pathLength - 1].pathMetric;
    detourTree.insert(detour);
  }

  int numPathsSearched   = 0;
  int numTBPathsSearched = 0;

  while (numTBPathsSearched < std::pow(2, Kconv)) {
    DetourObject detour = detourTree.pop();
    std::vector<int> path(lowrate_pathLength);

    int newTracebackStage          = lowrate_pathLength - 1;
    float forwardPartialPathMetric = 0;
    int currentState               = detour.startingState;

    // if we are taking a detour from a previous path, we skip backwards to the point where we take
    // the detour from the previous path
    if (detour.originalPathIndex != -1) {
      forwardPartialPathMetric = detour.forwardPathMetric;
      newTracebackStage        = detour.detourStage;

      // while we only need to copy the path from the detour to the end, this simplifies things,
      // and we'll write over the earlier data in any case
      path         = previousPaths[detour.originalPathIndex];
      currentState = path[newTracebackStage];

      float suboptimalPathMetric = trellisInfo[currentState][newTracebackStage].suboptimalPathMetric;

      currentState = trellisInfo[currentState][newTracebackStage].suboptimalFatherState;
      newTracebackStage--;

      float prevPathMetric = trellisInfo[currentState][newTracebackStage].pathMetric;

      forwardPartialPathMetric += suboptimalPathMetric - prevPathMetric;
    }
    path[newTracebackStage] = currentState;

    // actually tracing back
    for (int stage = newTracebackStage; stage > 0; stage--) {
      float suboptimalPathMetric = trellisInfo[currentState][stage].suboptimalPathMetric;
      float currPathMetric       = trellisInfo[currentState][stage].pathMetric;

      // if there is a detour we add to the detourTree
      if (trellisInfo[currentState][stage].suboptimalFatherState != -1) {
        DetourObject localDetour;
        localDetour.detourStage       = stage;
        localDetour.originalPathIndex = numPathsSearched;
        localDetour.pathMetric        = suboptimalPathMetric + forwardPartialPathMetric;
        localDetour.forwardPathMetric = forwardPartialPathMetric;
        localDetour.startingState     = detour.startingState;
        detourTree.insert(localDetour);
      }
      currentState         = trellisInfo[currentState][stage].optimalFatherState;
      float prevPathMetric = trellisInfo[currentState][stage - 1].pathMetric;
      forwardPartialPathMetric += currPathMetric - prevPathMetric;
      path[stage - 1] = currentState;
    } // for(int stage = newTracebackStage; stage > 0; stage--)

    previousPaths.push_back(path);
    std::vector<int> message       = pathToMessage(path);
    std::vector<int> cand_codeword = pathToCodeword(path);
    std::vector<int> codeword      = bpsk::demodulate(cand_codeword);

    /* Computes ESD & ED */
    std::vector<int> ED = crc::remdr_slidingWindow(message, CRC_VEC);
    int ESD             = path[0] ^ path[this->lowrate_pathLength - 1];
    int ED_int          = utils::short_int_vector_to_int(ED);

    if (ESD == 0) {
      numTBPathsSearched++;
      /* Computes Hamming & Euclidean distances */
      int hamming_dist = 0;
      // float euclidean_dist = 0.0f;
      for (size_t i = 0; i < all_zeros_codeword.size(); i++) {
        if (std::find(PUNCTURING_INDICES.begin(), PUNCTURING_INDICES.end(), i) == PUNCTURING_INDICES.end()) {
          hamming_dist += (int)cand_codeword[i] != all_ones_modulated_codeword[i];
          // euclidean_dist += std::pow(std::abs(cand_codeword[i]-codeword[i]),2);
        }
      }

      if (TB_Offsets_info[ED_int].find(hamming_dist) == TB_Offsets_info[ED_int].end()) {
        // key does not exist
        TB_Offsets_info[ED_int][hamming_dist] = 1;
      } else {
        // key already exists
        TB_Offsets_info[ED_int][hamming_dist] += 1;
      }
      TB_Offsets_CWD[ED_int][hamming_dist].push_back(codeword);
      TB_Offsets_MSG[ED_int][hamming_dist].push_back(message);

      // std::cout << "candidate codeword: " << numTBPathsSearched << std::endl;
      // std::cout << "outputing syndrome: " << syndrome << std::endl;
      // std::cout << "path[0]: " << path[0] << "; path[Kconv]: " <<
      // path[this->lowrate_pathLength-1] << "; ESD: " << ESD << std::endl; std::cout << "message:
      // "; utils::print_int_vector(message); std::cout << "codeword: ";
      // utils::print_int_vector(cand_codeword); std::cout << "hamming dist: " << hamming_dist << ";
      // euclidean dist: " << euclidean_dist << "; forwardPathMetric: " << forwardPartialPathMetric
      // << std::endl; std::cout << std::endl;
    }

    numPathsSearched++;
  } // while(true)

  /* Write to json file */
  // Write down codewords
  for (size_t i_m = 0; i_m < TB_Offsets_CWD.size(); i_m++) {
    json current_map_json = json::object();
    for (auto const& [dist, vv_cwds] : TB_Offsets_CWD[i_m]) {
      assert((int)vv_cwds.size() == TB_Offsets_info[i_m][dist]);
      std::string strKey       = fmt::format("{}", dist);
      current_map_json[strKey] = vv_cwds;
    }
    j_TB_Offsets_CWD.push_back(current_map_json);
  }
  auto out_cwds = fmt::output_file(TBOffsetsCWDFileName);
  out_cwds.print("{}", j_TB_Offsets_CWD.dump(4));

  for (size_t i_m = 0; i_m < TB_Offsets_MSG.size(); i_m++) {
    json current_map_json = json::object();
    for (auto const& [dist, vv_msgs] : TB_Offsets_MSG[i_m]) {
      assert((int)vv_msgs.size() == TB_Offsets_info[i_m][dist]);
      std::string strKey       = fmt::format("{}", dist);
      current_map_json[strKey] = vv_msgs;
    }
    j_TB_Offsets_MSG.push_back(current_map_json);
  }
  auto out_msgs = fmt::output_file(TBOffsetsMSGFileName);
  out_msgs.print("{}", j_TB_Offsets_MSG.dump(4));
}

void LowRateListDecoder::genNonTBOffsets() {
  /* Generates lookup table for SSD, where HR is trellis paths and covers all ELF and TB syndromes. */

  /* Output file name & output data structure */
  std::string OffsetsCWDFileName = OUTPUTFILEPATH + "NonTB_Offsets_CWD.json";
  // The outer vector corresponds to all possible syndromes; the key to map is hamming distance.
  std::vector<std::map<int, std::vector<std::vector<int>>>> NonTB_Offsets_CWD((int)std::pow(2, M + V));
  json j_NonTB_Offsets_CWD = json::array();

  std::string OffsetsMSGFileName = OUTPUTFILEPATH + "NonTB_Offsets_MSG.json";
  // The outer vector corresponds to all possible syndromes; the key to map is hamming distance.
  std::vector<std::map<int, std::vector<std::vector<int>>>> NonTB_Offsets_MSG((int)std::pow(2, M + V));
  json j_NonTB_Offsets_MSG = json::array();
  std::vector<std::map<int, int>> NonTB_Offsets_info((int)std::pow(2, M + V));

  std::vector<int> all_zeros_codeword(N, 0);
  std::vector<float> all_ones_modulated_codeword(N, 1.0f);
  std::vector<std::vector<cell>> trellisInfo; // trellisInfo is indexed [state][stage]
  trellisInfo = constructLowRateTrellis_Punctured(all_ones_modulated_codeword, PUNCTURING_INDICES);

  // start search
  MinHeap detourTree;
  std::vector<std::vector<int>> previousPaths;

  // create nodes for each valid ending state with no detours
  for (int i = 0; i < lowrate_numStates; i++) {
    DetourObject detour;
    detour.startingState = i;
    detour.pathMetric    = trellisInfo[i][lowrate_pathLength - 1].pathMetric;
    detourTree.insert(detour);
  }

  int numPathsSearched   = 0;
  int numTBPathsSearched = 0;

  while (numTBPathsSearched < std::pow(2, Kconv)) {
    DetourObject detour = detourTree.pop();
    std::vector<int> path(lowrate_pathLength);

    int newTracebackStage          = lowrate_pathLength - 1;
    float forwardPartialPathMetric = 0;
    int currentState               = detour.startingState;

    // if we are taking a detour from a previous path, we skip backwards to the point where we take
    // the detour from the previous path
    if (detour.originalPathIndex != -1) {
      forwardPartialPathMetric = detour.forwardPathMetric;
      newTracebackStage        = detour.detourStage;

      // while we only need to copy the path from the detour to the end, this simplifies things,
      // and we'll write over the earlier data in any case
      path         = previousPaths[detour.originalPathIndex];
      currentState = path[newTracebackStage];

      float suboptimalPathMetric = trellisInfo[currentState][newTracebackStage].suboptimalPathMetric;

      currentState = trellisInfo[currentState][newTracebackStage].suboptimalFatherState;
      newTracebackStage--;

      float prevPathMetric = trellisInfo[currentState][newTracebackStage].pathMetric;

      forwardPartialPathMetric += suboptimalPathMetric - prevPathMetric;
    }
    path[newTracebackStage] = currentState;

    // actually tracing back
    for (int stage = newTracebackStage; stage > 0; stage--) {
      float suboptimalPathMetric = trellisInfo[currentState][stage].suboptimalPathMetric;
      float currPathMetric       = trellisInfo[currentState][stage].pathMetric;

      // if there is a detour we add to the detourTree
      if (trellisInfo[currentState][stage].suboptimalFatherState != -1) {
        DetourObject localDetour;
        localDetour.detourStage       = stage;
        localDetour.originalPathIndex = numPathsSearched;
        localDetour.pathMetric        = suboptimalPathMetric + forwardPartialPathMetric;
        localDetour.forwardPathMetric = forwardPartialPathMetric;
        localDetour.startingState     = detour.startingState;
        detourTree.insert(localDetour);
      }
      currentState         = trellisInfo[currentState][stage].optimalFatherState;
      float prevPathMetric = trellisInfo[currentState][stage - 1].pathMetric;
      forwardPartialPathMetric += currPathMetric - prevPathMetric;
      path[stage - 1] = currentState;
    } // for(int stage = newTracebackStage; stage > 0; stage--)

    previousPaths.push_back(path);
    std::vector<int> message       = pathToMessage(path);
    std::vector<int> cand_codeword = pathToCodeword(path);
    std::vector<int> codeword      = bpsk::demodulate(cand_codeword);

    /* Computes ESD & ED */
    std::vector<int> ED = crc::remdr_slidingWindow(message, CRC_VEC);
    int ESD             = path[0] ^ path[this->lowrate_pathLength - 1];
    int ED_int          = utils::short_int_vector_to_int(ED);
    int offset_map_idx  = ED_int * (std::pow(2, V)) + ESD;

    if (ESD == 0) { numTBPathsSearched++; }

    /* Computes Hamming & Euclidean distances */
    int hamming_dist = 0;
    // float euclidean_dist = 0.0f;
    for (size_t i = 0; i < all_zeros_codeword.size(); i++) {
      if (std::find(PUNCTURING_INDICES.begin(), PUNCTURING_INDICES.end(), i) == PUNCTURING_INDICES.end()) {
        hamming_dist += (int)cand_codeword[i] != all_ones_modulated_codeword[i];
        // euclidean_dist += std::pow(std::abs(cand_codeword[i]-codeword[i]),2);
      }
    }

    if (NonTB_Offsets_info[offset_map_idx].find(hamming_dist) == NonTB_Offsets_info[offset_map_idx].end()) {
      // key does not exist
      NonTB_Offsets_info[offset_map_idx][hamming_dist] = 1;
    } else {
      // key already exists
      NonTB_Offsets_info[offset_map_idx][hamming_dist] += 1;
    }
    NonTB_Offsets_MSG[offset_map_idx][hamming_dist].push_back(message);
    NonTB_Offsets_CWD[offset_map_idx][hamming_dist].push_back(codeword);

    // std::cout << "candidate codeword: " << numTBPathsSearched << std::endl;
    // std::cout << "outputing syndrome: " << syndrome << std::endl;
    // std::cout << "path[0]: " << path[0] << "; path[Kconv]: " <<
    // path[this->lowrate_pathLength-1] << "; ESD: " << ESD << std::endl; std::cout << "message:
    // "; utils::print_int_vector(message); std::cout << "codeword: ";
    // utils::print_int_vector(cand_codeword); std::cout << "hamming dist: " << hamming_dist << ";
    // euclidean dist: " << euclidean_dist << "; forwardPathMetric: " << forwardPartialPathMetric
    // << std::endl; std::cout << std::endl;

    numPathsSearched++;
  } // while (numTBPathsSearched < std::pow(2, Kconv))

  /* Write to json file */
  // Write down messages
  for (size_t i_m = 0; i_m < NonTB_Offsets_MSG.size(); i_m++) {
    json current_map_json = json::object();
    for (auto const& [dist, vv_msgs] : NonTB_Offsets_MSG[i_m]) {
      assert((int)vv_msgs.size() == NonTB_Offsets_info[i_m][dist]);
      std::string strKey       = fmt::format("{}", dist);
      current_map_json[strKey] = vv_msgs;
    }
    j_NonTB_Offsets_MSG.push_back(current_map_json);
  }
  auto out_msgs = fmt::output_file(OffsetsMSGFileName);
  out_msgs.print("{}", j_NonTB_Offsets_MSG.dump(4));
  // Write down codewords
  for (size_t i_m = 0; i_m < NonTB_Offsets_CWD.size(); i_m++) {
    json current_map_json = json::object();
    for (auto const& [dist, vv_cwds] : NonTB_Offsets_CWD[i_m]) {
      assert((int)vv_cwds.size() == NonTB_Offsets_info[i_m][dist]);
      std::string strKey       = fmt::format("{}", dist);
      current_map_json[strKey] = vv_cwds;
    }
    j_NonTB_Offsets_CWD.push_back(current_map_json);
  }
  auto out_cwds = fmt::output_file(OffsetsCWDFileName);
  out_cwds.print("{}", j_NonTB_Offsets_CWD.dump(4));
}

void LowRateListDecoder::findNonTBCosetLeaders(const std::vector<float>& receivedMessage,
                                               const std::vector<int>& punctured_indices,
                                               const int listsize) {
  /* Output file name & output data structure */
  std::string cosetLeaderMsgFileName = OUTPUTFILEPATH + "nonTB_coset_leaders_msgs.json";
  std::map<std::pair<int, int>, std::vector<std::vector<int>>> coset_leaders_msg;
  json j_coset_leaders_msgs;

  std::string cosetLeaderCwdFileName = OUTPUTFILEPATH + "nonTB_coset_leaders_cwds.json";
  std::map<std::pair<int, int>, std::vector<std::vector<int>>> coset_leaders_cwd;
  json j_coset_leaders_cwds;

  std::map<std::pair<int, int>, std::pair<int, int>> coset_info;

  // start search
  std::vector<std::vector<cell>> trellisInfo; // trellisInfo is indexed [state][stage]
  trellisInfo = constructLowRateTrellis_Punctured(receivedMessage, punctured_indices);
  MinHeap detourTree;
  std::vector<std::vector<int>> previousPaths;

  // create nodes for each valid ending state with no detours
  for (int i = 0; i < lowrate_numStates; i++) {
    DetourObject detour;
    detour.startingState = i;
    detour.pathMetric    = trellisInfo[i][lowrate_pathLength - 1].pathMetric;
    detourTree.insert(detour);
  }

  int numPathsSearched = 0;
  // int numTBPathsSearched = 0;

  while (numPathsSearched < listsize) {
    DetourObject detour = detourTree.pop();
    std::vector<int> path(lowrate_pathLength);

    int newTracebackStage          = lowrate_pathLength - 1;
    float forwardPartialPathMetric = 0;
    int currentState               = detour.startingState;

    // if we are taking a detour from a previous path, we skip backwards to the point where we take
    // the detour from the previous path
    if (detour.originalPathIndex != -1) {
      forwardPartialPathMetric = detour.forwardPathMetric;
      newTracebackStage        = detour.detourStage;

      // while we only need to copy the path from the detour to the end, this simplifies things,
      // and we'll write over the earlier data in any case
      path         = previousPaths[detour.originalPathIndex];
      currentState = path[newTracebackStage];

      float suboptimalPathMetric = trellisInfo[currentState][newTracebackStage].suboptimalPathMetric;

      currentState = trellisInfo[currentState][newTracebackStage].suboptimalFatherState;
      newTracebackStage--;

      float prevPathMetric = trellisInfo[currentState][newTracebackStage].pathMetric;

      forwardPartialPathMetric += suboptimalPathMetric - prevPathMetric;
    }
    path[newTracebackStage] = currentState;

    // actually tracing back
    for (int stage = newTracebackStage; stage > 0; stage--) {
      float suboptimalPathMetric = trellisInfo[currentState][stage].suboptimalPathMetric;
      float currPathMetric       = trellisInfo[currentState][stage].pathMetric;

      // if there is a detour we add to the detourTree
      if (trellisInfo[currentState][stage].suboptimalFatherState != -1) {
        DetourObject localDetour;
        localDetour.detourStage       = stage;
        localDetour.originalPathIndex = numPathsSearched;
        localDetour.pathMetric        = suboptimalPathMetric + forwardPartialPathMetric;
        localDetour.forwardPathMetric = forwardPartialPathMetric;
        localDetour.startingState     = detour.startingState;
        detourTree.insert(localDetour);
      }
      currentState         = trellisInfo[currentState][stage].optimalFatherState;
      float prevPathMetric = trellisInfo[currentState][stage - 1].pathMetric;
      forwardPartialPathMetric += currPathMetric - prevPathMetric;
      path[stage - 1] = currentState;
    } // for(int stage = newTracebackStage; stage > 0; stage--)

    previousPaths.push_back(path);
    std::vector<int> message       = pathToMessage(path);
    std::vector<int> cand_codeword = pathToCodeword(path);
    std::vector<int> codeword      = bpsk::demodulate(cand_codeword);

    /* Computes ESD & ED */
    std::vector<int> ED             = crc::remdr_slidingWindow(message, CRC_VEC);
    int ESD                         = path[0] ^ path[this->lowrate_pathLength - 1];
    int ED_int                      = utils::short_int_vector_to_int(ED);
    std::pair<int, int> ED_ESD_pair = std::make_pair(ED_int, ESD);

    /* Computes Hamming & Euclidean distances */
    int hamming_dist = 0;
    // float euclidean_dist = 0.0f;
    for (size_t i = 0; i < codeword.size(); i++) {
      if (std::find(punctured_indices.begin(), punctured_indices.end(), i) == punctured_indices.end()) {
        hamming_dist += (int)codeword[i] != 0;
        // euclidean_dist += std::pow(std::abs(cand_codeword[i]-codeword[i]),2);
      }
    }

    if (hamming_dist == 9) {
      fmt::print("non-TB coset leaders search ended.\n");
      break;
    }

    /* Record Coset Leaders */
    if (coset_info.find(ED_ESD_pair) == coset_info.end()) {
      coset_info[ED_ESD_pair] = std::make_pair(hamming_dist, 1);
      coset_leaders_msg[ED_ESD_pair].push_back(message);
      coset_leaders_cwd[ED_ESD_pair].push_back(codeword);
    } else {
      if (coset_info[ED_ESD_pair].first > hamming_dist) {
        coset_info[ED_ESD_pair] = std::make_pair(hamming_dist, 1);
        coset_leaders_msg[ED_ESD_pair].clear();
        coset_leaders_msg[ED_ESD_pair].push_back(message);
        coset_leaders_cwd[ED_ESD_pair].clear();
        coset_leaders_cwd[ED_ESD_pair].push_back(codeword);
      } else if (coset_info[ED_ESD_pair].first == hamming_dist) {
        coset_info[ED_ESD_pair].second++;
        coset_leaders_msg[ED_ESD_pair].push_back(message);
        coset_leaders_cwd[ED_ESD_pair].push_back(codeword);
      }
    }

    // numTBPathsSearched++;
    // std::cout << "candidate codeword: " << numTBPathsSearched << std::endl;
    // std::cout << "outputing syndrome: " << ED_int << std::endl;
    // std::cout << "path[0]: " << path[0] << "; path[Kconv]: " << path[this->lowrate_pathLength-1]
    // << "; ESD: " << ESD << std::endl; std::cout << "message: "; utils::print_int_vector(message);
    // std::cout << "codeword: "; utils::print_int_vector(cand_codeword);
    // std::cout << "hamming dist: " << hamming_dist << "; euclidean dist: " << euclidean_dist << ";
    // forwardPathMetric: " << forwardPartialPathMetric << std::endl; std::cout << std::endl;

    numPathsSearched++;
  } // while(numPathsSearched < this->listSize)

  try {
    std::string cosetInfoFileName = OUTPUTFILEPATH + "nonTB_coset_info.txt";
    auto CosetInfoFile            = fmt::output_file(cosetInfoFileName);
    for (auto kv_pair : coset_info) {
      CosetInfoFile.print("ED_ESD pair: {}, {}\n", kv_pair.first.first, kv_pair.first.second);
      CosetInfoFile.print("hamming_dist: {}\n", kv_pair.second.first);
      CosetInfoFile.print("count: {}\n\n", kv_pair.second.second);
    }
  } catch (const std::exception& e) { std::cerr << "Error: " << e.what() << std::endl; }

  /* Write to json file */
  // Write down messages
  for (auto const& [key, msgs] : coset_leaders_msg) {
    assert((int)msgs.size() == coset_info[key].second);
    std::string strKey           = fmt::format("{},{}", key.first, key.second);
    j_coset_leaders_msgs[strKey] = msgs;
  }
  auto out_msgs = fmt::output_file(cosetLeaderMsgFileName);
  out_msgs.print("{}", j_coset_leaders_msgs.dump(4));

  // Write down codewords
  for (auto const& [key, cwds] : coset_leaders_cwd) {
    assert((int)cwds.size() == coset_info[key].second);
    std::string strKey           = fmt::format("{},{}", key.first, key.second);
    j_coset_leaders_cwds[strKey] = cwds;
  }
  auto out_cwds = fmt::output_file(cosetLeaderCwdFileName);
  out_cwds.print("{}", j_coset_leaders_cwds.dump(4));
}
