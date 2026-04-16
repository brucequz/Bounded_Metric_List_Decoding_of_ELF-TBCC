#include "../consts.h"
#include "../include/json.hpp"
#include "../include/lowRateListDecoder.h"
#include "../include/minHeap.h"
#include "../include/namespace.h"

#include <cassert>
#include <cmath>
#include <cstddef>
#include <fmt/os.h>
#include <fmt/ranges.h>
#include <map>
#include <utility>
#include <vector>
using json = nlohmann::json;

MessageInformation LowRateListDecoder::ssdSLVDDecoding(
    const std::vector<float>& receivedMessage,
    const std::vector<int>& transmittedWord,
    const std::vector<int>& punctured_indices,
    const std::map<std::pair<int, int>, std::vector<std::vector<int>>>& cosetLeadersMsgs,
    const std::map<std::pair<int, int>, std::vector<std::vector<int>>>& cosetLeadersCwds,
    const std::vector<std::vector<int>>& gabrielNeighbors,
    SSD_TYPE type) {
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
    detour.pathMetric    = trellisInfo[i][lowrate_pathLength - 1].pathMetric;
    detourTree.insert(detour);
  }

  int numPathsSearched = 0;
  // int numTBPathsSearched = 0;

  while (true) {
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
    // fmt::print("forwardPartialPathMetric: {}\n",forwardPartialPathMetric);

    /* Translates path to msg and cwd */
    std::vector<int> message       = pathToMessage(path);
    std::vector<int> cand_codeword = pathToCodeword(path);
    std::vector<int> codeword      = bpsk::demodulate(cand_codeword);

    /* Computes ESD & ED */
    std::vector<int> ED             = crc::remdr_slidingWindow(message, CRC_VEC);
    int ESD                         = path[0] ^ path[this->lowrate_pathLength - 1];
    int ED_int                      = utils::short_int_vector_to_int(ED);
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
      output.message  = message;
      output.codeword = cand_codeword;
      break;
    }

    if (type == SSD_TYPE::TRELLIS_TYPE || (type == SSD_TYPE::TBCC_TYPE && ESD == 0)) {
      /* ED correction with coset leaders */
      std::vector<int> best_coset_corrected_cwd;
      float best_coset_corrected_metric = 100000.0f; //<! some big number
      for (size_t i_cl = 0; i_cl < cosetLeadersMsgs.at(ED_ESD_pair).size(); i_cl++) {
        // fmt::print("There are {} coset leaders. \n", cosetLeadersMsgs.at(ED_ESD_pair).size());

        /* - Message ED correction */
        std::vector<int> ED_corrected_msg;
        for (size_t i_m = 0; i_m < message.size(); i_m++) {
          ED_corrected_msg.push_back(message[i_m] ^ cosetLeadersMsgs.at(ED_ESD_pair)[i_cl][i_m]);
        }
        // std::vector<int> new_ED = crc::remdr_slidingWindow(ED_corrected_msg, CRC_VEC);
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
        float euclidean_dist =
            distance::compute_euclidean_distance(modulated_ED_corrected_cwd, receivedMessage, PUNCTURING_INDICES);
        // fmt::print("euclidean distance to received word: ", euclidean_dist);

        /* Update the best coset corrected metric & cwd */
        if (euclidean_dist < best_coset_corrected_metric) {
          best_coset_corrected_metric = euclidean_dist;
          best_coset_corrected_cwd    = modulated_ED_corrected_cwd;
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
          bch_tbcc_candidate[i] = -2.0f * (best_demod_coset_corrected_cwd[i] ^ gabrielNeighbors[i_g][i]) + 1;
          correlation += receivedMessage[i] * bch_tbcc_candidate[i];
        }
        if (correlation > best_correlation_metric) {
          best_correlation_metric = correlation;
          best_gabriel_neighbor   = bch_tbcc_candidate;
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

MessageInformation LowRateListDecoder::ssd_TrellisPath(
    const std::vector<float>& receivedMessage,
    const std::vector<std::map<int, std::vector<std::vector<int>>>>& TB_Offsets_CWD,
    const std::vector<std::map<int, std::vector<std::vector<int>>>>& TB_Offsets_MSG) {
  // start search
  std::vector<std::vector<cell>> trellisInfo; // trellisInfo is indexed [state][stage]
  trellisInfo = constructLowRateTrellis_Punctured(receivedMessage, PUNCTURING_INDICES);
  MessageInformation output;
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
    path                       = previousPaths[detour.originalPathIndex];
    currentState               = path[newTracebackStage];
    float suboptimalPathMetric = trellisInfo[currentState][newTracebackStage].suboptimalPathMetric;
    currentState               = trellisInfo[currentState][newTracebackStage].suboptimalFatherState;
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

  /* Translates path to msg and cwd */
  std::vector<int> message       = pathToMessage(path);
  std::vector<int> cand_codeword = pathToCodeword(path);
  std::vector<int> codeword      = bpsk::demodulate(cand_codeword);

  /* Computes ESD & ED */
  std::vector<int> ED = crc::remdr_slidingWindow(message, CRC_VEC);
  int ESD             = path[0] ^ path[this->lowrate_pathLength - 1];
  int ED_int          = utils::short_int_vector_to_int(ED);
  int offset_map_idx  = ED_int * (std::pow(2, V)) + ESD;

  // if ESD == 0 && ED_int == 0, return the current answer
  if (ESD == 0 && ED_int == 0) {
    output.message  = message;
    output.codeword = codeword;
    return output;
  }
  // if syndrome != 0, use offset_map_idx to find the correct map to correct syndromes
  float corr_metric = 0.0f;
  std::vector<int> ED_corrected_cwd;
  for (auto const& [h_dist, vv_cwds] : TB_Offsets_CWD[offset_map_idx]) {
    if (h_dist <= OFFSET_SPHERE_RADIUS) {
      for (auto const& v_cwd : vv_cwds) {
        /* - Codeword ED correction */
        std::vector<int> cand_correction;
        for (size_t i_c = 0; i_c < codeword.size(); i_c++) { cand_correction.push_back(codeword[i_c] ^ v_cwd[i_c]); }
        // fmt::print("cand_correction: {}\n", cand_correction);

        float temp_corr = 0.0f;
        for (size_t i = 0; i < cand_correction.size(); i++) {
          temp_corr += receivedMessage[i] * (-2 * cand_correction[i] + 1);
        }
        // fmt::print("temp_corr: {}\n", temp_corr);

        ED_corrected_cwd = (temp_corr > corr_metric) ? cand_correction : ED_corrected_cwd;
        corr_metric      = (temp_corr > corr_metric) ? temp_corr : corr_metric;
      } // for (auto const& v_cwd : vv_cwds)
    }
  } // for (auto const& [h_dist, vv_cwds] : TB_Offsets_CWD[offset_map_idx])

  output.codeword = ED_corrected_cwd;
  output.message  = message;
  return output;

  numPathsSearched++; //<! this has to stay here at the end
}