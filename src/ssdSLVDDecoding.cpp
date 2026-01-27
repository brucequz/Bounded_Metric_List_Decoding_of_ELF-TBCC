#include "../consts.h"
#include "../include/fileHandler.h"
#include "../include/json.hpp"
#include "../include/lowRateListDecoder.h"
#include "../include/minHeap.h"
#include "../include/namespace.h"

#include <cassert>
#include <cstddef>
#include <fmt/os.h>
#include <fmt/ranges.h>
#include <map>
#include <utility>
#include <vector>
using json = nlohmann::json;

MessageInformation LowRateListDecoder::findTBCosetLeaders(const std::vector<float>& receivedMessage,
                                                          const std::vector<int>& punctured_indices,
                                                          const std::vector<int>& codeword,
                                                          const int listsize) {
  /* Output file name & output data structure */
  std::string cosetLeaderFileName = OUTPUTFILEPATH + "coset_leaders_msgs.txt";
  std::map<int, std::vector<std::vector<int>>> coset_leaders;

  std::vector<std::vector<cell>> trellisInfo; // trellisInfo is indexed [state][stage]
  trellisInfo = constructLowRateTrellis_Punctured(receivedMessage, punctured_indices);

  // start search
  MessageInformation output;
  MinHeap detourTree;
  std::vector<std::vector<int>> previousPaths;

  // create nodes for each valid ending state with no detours
  for (int i = 0; i < lowrate_numStates; i++) {
    DetourObject detour;
    detour.startingState = i;
    detour.pathMetric = trellisInfo[i][lowrate_pathLength - 1].pathMetric;
    detourTree.insert(detour);
  }

  int numPathsSearched = 0;
  int numTBPathsSearched = 0;

  while (numPathsSearched < listsize) {
    DetourObject detour = detourTree.pop();
    std::vector<int> path(lowrate_pathLength);

    int newTracebackStage = lowrate_pathLength - 1;
    float forwardPartialPathMetric = 0;
    int currentState = detour.startingState;

    // if we are taking a detour from a previous path, we skip backwards to the point where we take
    // the detour from the previous path
    if (detour.originalPathIndex != -1) {
      forwardPartialPathMetric = detour.forwardPathMetric;
      newTracebackStage = detour.detourStage;

      // while we only need to copy the path from the detour to the end, this simplifies things,
      // and we'll write over the earlier data in any case
      path = previousPaths[detour.originalPathIndex];
      currentState = path[newTracebackStage];

      float suboptimalPathMetric =
          trellisInfo[currentState][newTracebackStage].suboptimalPathMetric;

      currentState = trellisInfo[currentState][newTracebackStage].suboptimalFatherState;
      newTracebackStage--;

      float prevPathMetric = trellisInfo[currentState][newTracebackStage].pathMetric;

      forwardPartialPathMetric += suboptimalPathMetric - prevPathMetric;
    }
    path[newTracebackStage] = currentState;

    // actually tracing back
    for (int stage = newTracebackStage; stage > 0; stage--) {
      float suboptimalPathMetric = trellisInfo[currentState][stage].suboptimalPathMetric;
      float currPathMetric = trellisInfo[currentState][stage].pathMetric;

      // if there is a detour we add to the detourTree
      if (trellisInfo[currentState][stage].suboptimalFatherState != -1) {
        DetourObject localDetour;
        localDetour.detourStage = stage;
        localDetour.originalPathIndex = numPathsSearched;
        localDetour.pathMetric = suboptimalPathMetric + forwardPartialPathMetric;
        localDetour.forwardPathMetric = forwardPartialPathMetric;
        localDetour.startingState = detour.startingState;
        detourTree.insert(localDetour);
      }
      currentState = trellisInfo[currentState][stage].optimalFatherState;
      float prevPathMetric = trellisInfo[currentState][stage - 1].pathMetric;
      forwardPartialPathMetric += currPathMetric - prevPathMetric;
      path[stage - 1] = currentState;
    } // for(int stage = newTracebackStage; stage > 0; stage--)

    previousPaths.push_back(path);
    std::vector<int> message = pathToMessage(path);
    std::vector<int> cand_codeword = pathToCodeword(path);

    /* Computes ESD & ED */
    std::vector<int> ED = crc::remdr_slidingWindow(message, CRC_VEC);
    // int ESD = path[0] ^ path[this->lowrate_pathLength-1];
    int ED_int = utils::short_int_vector_to_int(ED);

    /* Computes Hamming & Euclidean distances */
    int hamming_dist = 0;
    // float euclidean_dist = 0.0f;
    for (size_t i = 0; i < codeword.size(); i++) {
      if (std::find(punctured_indices.begin(), punctured_indices.end(), i) ==
          punctured_indices.end()) {
        hamming_dist += (int) cand_codeword[i] != codeword[i];
        // euclidean_dist += std::pow(std::abs(cand_codeword[i]-codeword[i]),2);
      }
    }

    if (path[0] == path.back()) {

      coset_leaders[ED_int].push_back(message);
      numTBPathsSearched++;
      std::cout << "candidate codeword: " << numTBPathsSearched << std::endl;
      // std::cout << "outputing syndrome: " << syndrome << std::endl;
      // std::cout << "path[0]: " << path[0] << "; path[Kconv]: " <<
      // path[this->lowrate_pathLength-1] << "; ESD: " << ESD << std::endl; std::cout << "message:
      // "; utils::print_int_vector(message); std::cout << "codeword: ";
      // utils::print_int_vector(cand_codeword); std::cout << "hamming dist: " << hamming_dist << ";
      // euclidean dist: " << euclidean_dist << "; forwardPathMetric: " << forwardPartialPathMetric
      // << std::endl; std::cout << std::endl;
    }

    numPathsSearched++;
  } // while(numPathsSearched < this->listSize)

  for (auto kv_pair : coset_leaders) {
    try {
      FileHandler CosetLeadersFile(cosetLeaderFileName);
      CosetLeadersFile.write2DVector(kv_pair.second);
    } catch (const std::exception& e) {
      std::cerr << "Error: " << e.what() << std::endl;
    }
  }

  output.listSizeExceeded = true;
  return output;
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
    detour.pathMetric = trellisInfo[i][lowrate_pathLength - 1].pathMetric;
    detourTree.insert(detour);
  }

  int numPathsSearched = 0;
  // int numTBPathsSearched = 0;

  while (numPathsSearched < listsize) {
    DetourObject detour = detourTree.pop();
    std::vector<int> path(lowrate_pathLength);

    int newTracebackStage = lowrate_pathLength - 1;
    float forwardPartialPathMetric = 0;
    int currentState = detour.startingState;

    // if we are taking a detour from a previous path, we skip backwards to the point where we take
    // the detour from the previous path
    if (detour.originalPathIndex != -1) {
      forwardPartialPathMetric = detour.forwardPathMetric;
      newTracebackStage = detour.detourStage;

      // while we only need to copy the path from the detour to the end, this simplifies things,
      // and we'll write over the earlier data in any case
      path = previousPaths[detour.originalPathIndex];
      currentState = path[newTracebackStage];

      float suboptimalPathMetric =
          trellisInfo[currentState][newTracebackStage].suboptimalPathMetric;

      currentState = trellisInfo[currentState][newTracebackStage].suboptimalFatherState;
      newTracebackStage--;

      float prevPathMetric = trellisInfo[currentState][newTracebackStage].pathMetric;

      forwardPartialPathMetric += suboptimalPathMetric - prevPathMetric;
    }
    path[newTracebackStage] = currentState;

    // actually tracing back
    for (int stage = newTracebackStage; stage > 0; stage--) {
      float suboptimalPathMetric = trellisInfo[currentState][stage].suboptimalPathMetric;
      float currPathMetric = trellisInfo[currentState][stage].pathMetric;

      // if there is a detour we add to the detourTree
      if (trellisInfo[currentState][stage].suboptimalFatherState != -1) {
        DetourObject localDetour;
        localDetour.detourStage = stage;
        localDetour.originalPathIndex = numPathsSearched;
        localDetour.pathMetric = suboptimalPathMetric + forwardPartialPathMetric;
        localDetour.forwardPathMetric = forwardPartialPathMetric;
        localDetour.startingState = detour.startingState;
        detourTree.insert(localDetour);
      }
      currentState = trellisInfo[currentState][stage].optimalFatherState;
      float prevPathMetric = trellisInfo[currentState][stage - 1].pathMetric;
      forwardPartialPathMetric += currPathMetric - prevPathMetric;
      path[stage - 1] = currentState;
    } // for(int stage = newTracebackStage; stage > 0; stage--)

    previousPaths.push_back(path);
    std::vector<int> message = pathToMessage(path);
    std::vector<int> cand_codeword = pathToCodeword(path);
    std::vector<int> codeword = bpsk::demodulate(cand_codeword);

    /* Computes ESD & ED */
    std::vector<int> ED = crc::remdr_slidingWindow(message, CRC_VEC);
    int ESD = path[0] ^ path[this->lowrate_pathLength - 1];
    int ED_int = utils::short_int_vector_to_int(ED);
    std::pair<int, int> ED_ESD_pair = std::make_pair(ED_int, ESD);

    /* Computes Hamming & Euclidean distances */
    int hamming_dist = 0;
    // float euclidean_dist = 0.0f;
    for (size_t i = 0; i < codeword.size(); i++) {
      if (std::find(punctured_indices.begin(), punctured_indices.end(), i) ==
          punctured_indices.end()) {
        hamming_dist += (int) codeword[i] != 0;
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
    auto CosetInfoFile = fmt::output_file(cosetInfoFileName);
    for (auto kv_pair : coset_info) {
      CosetInfoFile.print("ED_ESD pair: {}, {}\n", kv_pair.first.first, kv_pair.first.second);
      CosetInfoFile.print("hamming_dist: {}\n", kv_pair.second.first);
      CosetInfoFile.print("count: {}\n\n", kv_pair.second.second);
    }
  } catch (const std::exception& e) {
    std::cerr << "Error: " << e.what() << std::endl;
  }

  /* Write to json file */
  // Write down messages
  for (auto const& [key, msgs] : coset_leaders_msg) {
    assert((int) msgs.size() == coset_info[key].second);
    std::string strKey = fmt::format("{},{}", key.first, key.second);
    j_coset_leaders_msgs[strKey] = msgs;
  }
  auto out_msgs = fmt::output_file(cosetLeaderMsgFileName);
  out_msgs.print("{}", j_coset_leaders_msgs.dump(4));

  // Write down codewords
  for (auto const& [key, cwds] : coset_leaders_cwd) {
    assert((int) cwds.size() == coset_info[key].second);
    std::string strKey = fmt::format("{},{}", key.first, key.second);
    j_coset_leaders_cwds[strKey] = cwds;
  }
  auto out_cwds = fmt::output_file(cosetLeaderCwdFileName);
  out_cwds.print("{}", j_coset_leaders_cwds.dump(4));
}

MessageInformation LowRateListDecoder::ssdSLVDDecoding(
    const std::vector<float>& receivedMessage, const std::vector<int>& transmittedWord,
    const std::vector<int>& punctured_indices,
    const std::map<std::pair<int, int>, std::vector<std::vector<int>>>& cosetLeadersMsgs,
    const std::map<std::pair<int, int>, std::vector<std::vector<int>>>& cosetLeadersCwds,
    const std::vector<std::vector<int>>& gabrielNeighbors) {

  // float transmitted_euclidean_dist = utils::euclidean_distance(receivedMessage, transmittedWord,
  // PUNCTURING_INDICES); fmt::print("transmitted_euclidean_dist: {}\n",
  // transmitted_euclidean_dist);

  // start search
  std::vector<std::vector<cell>> trellisInfo; // trellisInfo is indexed [state][stage]
  trellisInfo = constructLowRateTrellis_Punctured(receivedMessage, punctured_indices);
  MessageInformation output;
  MinHeap detourTree;
  std::vector<std::vector<int>> previousPaths;

  // create nodes for each valid ending state with no detours
  for (int i = 0; i < lowrate_numStates; i++) {
    DetourObject detour;
    detour.startingState = i;
    detour.pathMetric = trellisInfo[i][lowrate_pathLength - 1].pathMetric;
    detourTree.insert(detour);
  }

  int numPathsSearched = 0;
  // int numTBPathsSearched = 0;

  while (true) {

    DetourObject detour = detourTree.pop();
    std::vector<int> path(lowrate_pathLength);

    int newTracebackStage = lowrate_pathLength - 1;
    float forwardPartialPathMetric = 0;
    int currentState = detour.startingState;

    // if we are taking a detour from a previous path, we skip backwards to the point where we take
    // the detour from the previous path
    if (detour.originalPathIndex != -1) {
      forwardPartialPathMetric = detour.forwardPathMetric;
      newTracebackStage = detour.detourStage;

      // while we only need to copy the path from the detour to the end, this simplifies things,
      // and we'll write over the earlier data in any case
      path = previousPaths[detour.originalPathIndex];
      currentState = path[newTracebackStage];

      float suboptimalPathMetric =
          trellisInfo[currentState][newTracebackStage].suboptimalPathMetric;

      currentState = trellisInfo[currentState][newTracebackStage].suboptimalFatherState;
      newTracebackStage--;

      float prevPathMetric = trellisInfo[currentState][newTracebackStage].pathMetric;

      forwardPartialPathMetric += suboptimalPathMetric - prevPathMetric;
    }
    path[newTracebackStage] = currentState;

    // actually tracing back
    for (int stage = newTracebackStage; stage > 0; stage--) {
      float suboptimalPathMetric = trellisInfo[currentState][stage].suboptimalPathMetric;
      float currPathMetric = trellisInfo[currentState][stage].pathMetric;

      // if there is a detour we add to the detourTree
      if (trellisInfo[currentState][stage].suboptimalFatherState != -1) {
        DetourObject localDetour;
        localDetour.detourStage = stage;
        localDetour.originalPathIndex = numPathsSearched;
        localDetour.pathMetric = suboptimalPathMetric + forwardPartialPathMetric;
        localDetour.forwardPathMetric = forwardPartialPathMetric;
        localDetour.startingState = detour.startingState;
        detourTree.insert(localDetour);
      }
      currentState = trellisInfo[currentState][stage].optimalFatherState;
      float prevPathMetric = trellisInfo[currentState][stage - 1].pathMetric;
      forwardPartialPathMetric += currPathMetric - prevPathMetric;
      path[stage - 1] = currentState;
    } // for(int stage = newTracebackStage; stage > 0; stage--)

    previousPaths.push_back(path);
    // fmt::print("forwardPartialPathMetric: {}\n",forwardPartialPathMetric);

    /* Translates path to msg and cwd */
    std::vector<int> message = pathToMessage(path);
    std::vector<int> cand_codeword = pathToCodeword(path);
    std::vector<int> codeword = bpsk::demodulate(cand_codeword);

    /* Computes ESD & ED */
    std::vector<int> ED = crc::remdr_slidingWindow(message, CRC_VEC);
    int ESD = path[0] ^ path[this->lowrate_pathLength - 1];
    int ED_int = utils::short_int_vector_to_int(ED);
    std::pair<int, int> ED_ESD_pair = std::make_pair(ED_int, ESD);

    // std::cout << "candidate codeword: " << numPathsSearched << std::endl;
    // std::cout << "outputing syndrome: " << ED_int << std::endl;
    // std::cout << "path[0]: " << path[0] << "; path[Kconv]: " << path[this->lowrate_pathLength -
    // 1]
    //           << "; ESD: " << ESD << std::endl;
    // std::cout << "cand message: ";
    // utils::print_int_vector(message);
    // std::cout << "cand codeword: ";
    // utils::print_int_vector(cand_codeword);

    /* Return if ED == 0 and ESD == 0 */
    if (ED_int == 0 && ESD == 0) {
      output.message = message;
      output.codeword = cand_codeword;
      break;
    }

    if (ESD == 0) {
      /* ED correction with coset leaders */
      std::vector<int> best_coset_corrected_cwd;
      float best_coset_corrected_metric = 100000.0f; //<! some big number
      for (int i_cl = 0; i_cl < cosetLeadersMsgs.at(ED_ESD_pair).size(); i_cl++) {
        // fmt::print("There are {} coset leaders. \n", cosetLeadersMsgs.at(ED_ESD_pair).size());

        /* - Message ED correction */
        std::vector<int> ED_corrected_msg;
        for (size_t i_m = 0; i_m < message.size(); i_m++) {
          ED_corrected_msg.push_back(message[i_m] ^ cosetLeadersMsgs.at(ED_ESD_pair)[i_cl][i_m]);
        }
        std::vector<int> new_ED = crc::remdr_slidingWindow(ED_corrected_msg, CRC_VEC);
        // std::cout << "new_ED: " << utils::short_int_vector_to_int(new_ED) << std::endl;
        // std::cout << "ED_corrected_msg: ";
        // utils::print_int_vector(ED_corrected_msg);

        /* - Codeword ED correction */
        std::vector<int> ED_corrected_cwd;
        for (size_t i_c = 0; i_c < cand_codeword.size(); i_c++) {
          ED_corrected_cwd.push_back(codeword[i_c] ^ cosetLeadersCwds.at(ED_ESD_pair)[i_cl][i_c]);
        }
        // std::cout << "ED_corrected_cwd: ";
        // utils::print_int_vector(ED_corrected_cwd);

        /* BPSK modulated the corrected cwd */
        std::vector<int> modulated_ED_corrected_cwd;
        for (size_t i_c = 0; i_c < ED_corrected_cwd.size(); i_c++) {
          modulated_ED_corrected_cwd.push_back(-2 * ED_corrected_cwd[i_c] + 1);
        }
        // std::cout << "modulated_ED_corrected_cwd: ";
        // utils::print_int_vector(modulated_ED_corrected_cwd);

        /* Compute distance with received word */
        float euclidean_dist = utils::euclidean_distance(modulated_ED_corrected_cwd,
                                                         receivedMessage, PUNCTURING_INDICES);
        // fmt::print("euclidean distance to received word: ", euclidean_dist);

        /* Update the best coset corrected metric & cwd */
        if (euclidean_dist < best_coset_corrected_metric) {
          best_coset_corrected_metric = euclidean_dist;
          best_coset_corrected_cwd = modulated_ED_corrected_cwd;
        }
      }

      // fmt::print("best_coset_corrected_metric: {}\n", best_coset_corrected_metric);

      /* XOR the best coset corrected cwd with Gabriel Neighbors */
      std::vector<int> best_demod_coset_corrected_cwd(best_coset_corrected_cwd.size(), 0);
      for (size_t i = 0; i < best_coset_corrected_cwd.size(); i++) {
        best_demod_coset_corrected_cwd[i] = (1 - best_coset_corrected_cwd[i]) / 2;
      }
      std::vector<int> best_gabriel_neighbor;
      float best_correlation_metric = -1.0f;
      for (size_t i_g = 0; i_g < gabrielNeighbors.size(); i_g++) {
        // for each neighbor, compute metric
        float correlation = 0.0f;
        std::vector<int> bch_tbcc_candidate(N, 0);
        for (size_t i = 0; i < best_demod_coset_corrected_cwd.size(); i++) {
          bch_tbcc_candidate[i] =
              -2.0f * (best_demod_coset_corrected_cwd[i] ^ gabrielNeighbors[i_g][i]) + 1;
          correlation += receivedMessage[i] * bch_tbcc_candidate[i];
        }
        if (correlation > best_correlation_metric) {
          best_correlation_metric = correlation;
          best_gabriel_neighbor = bch_tbcc_candidate;
        }
      }

      // fmt::print("best correlation metric: {}\n", best_correlation_metric);
      // fmt::print("best_gabriel_neighbor : {}\n", best_gabriel_neighbor);

      output.codeword = best_gabriel_neighbor;
      break;
    }
    numPathsSearched++; //<! this has to stay here at the end
  } // while(true)

  return output;
}