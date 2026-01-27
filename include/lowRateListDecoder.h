#ifndef LOWRATELISTDECODER_H
#define LOWRATELISTDECODER_H

class FeedForwardTrellis;
class MinHeap;
class FileHandler;
#include "types.h"

#include <algorithm>
#include <climits>
#include <cmath>
#include <iostream>
#include <limits>
#include <map>
#include <vector>

enum class METRIC_TYPE { PRODUCT_METRIC = 1, EUCLIDEAN_METRIC = 2, LOG_EUCLIDEAN_METRIC = 3 };

class LowRateListDecoder {
public:
  LowRateListDecoder(FeedForwardTrellis* FT, int listSize, int crcLength, int crc,
                     char stopping_rule);

  MessageInformation decode(std::vector<float> receivedMessage, std::vector<int> punctured_indices);

  // TB
  MessageInformation lowRateDecoding_MaxListsize(const std::vector<float>& receivedMessage,
                                                 const std::vector<int>& punctured_indices);
  MessageInformation findTBCosetLeaders(const std::vector<float>& receivedMessage,
                                        const std::vector<int>& punctured_indices,
                                        const std::vector<int>& codeword, const int listsize);
  void findNonTBCosetLeaders(const std::vector<float>& receivedMessage,
                             const std::vector<int>& punctured_indices, const int listsize);
  MessageInformation lowRateDecoding_MaxMetric(std::vector<float> receivedMessage,
                                               std::vector<int> punctured_indices);
  MessageInformation lowRateDecoding_MaxAngle(std::vector<float> receivedMessage,
                                              std::vector<int> punctured_indices);
  MessageInformation lowRateDecoding_MaxAngle_ProductMetric(std::vector<float> receivedMessage,
                                                            std::vector<int> punctured_indices);
  MessageInformation ssdSLVDDecoding(
      const std::vector<float>& receivedMessage, const std::vector<int>& transmittedWord,
      const std::vector<int>& punctured_indices,
      const std::map<std::pair<int, int>, std::vector<std::vector<int>>>& cosetLeadersMsgs,
      const std::map<std::pair<int, int>, std::vector<std::vector<int>>>& cosetLeadersCwds,
      const std::vector<std::vector<int>>& gabrielNeighbors);

  // ZT
  MessageInformation lowRateDecoding_MaxListsize_ZT(std::vector<float>& receivedMessage);
  MessageInformation lowRateDecoding_MaxAngle_ProductMetric_ZT(std::vector<float> receivedMessage);
  MessageInformation
  lowRateDecoding_MaxMetric_EuclideanMetric_ZT(std::vector<float> receivedMessage);

  // BCJR
  MessageInformation lowRateDecoding_BCJR(std::vector<float> receivedMessage, float sigma_sqrd);

  // ROVA
  MessageInformation
  lowRateDecoding_SquaredDistanceMetric_ROVA_ZT(std::vector<float> receivedMessage,
                                                float sigma_sqrd, float rova_t);

  // genie-Aided
  MessageInformation genieAided_LowRateDecoding_MaxAngle_ProductMetric(
      std::vector<float> receivedMessage, std::vector<int> transmittedCodeword,
      std::vector<int> punctured_indices, bool push_to_boundary);
  MessageInformation genieAided_LowRateDecoding_MaxListsize_Collect(
      std::vector<float> receivedMessage, std::vector<int> punctured_indices,
      std::vector<int> transmittedCodeword, std::vector<float> sampling_points,
      std::vector<float> listsize_collect_sample_points);
  MessageInformation genieAided_LowRateDecoding_MaxAngle_ProductMetric_Collector(
      std::vector<float> receivedMessage, std::vector<int> transmittedCodeword,
      std::vector<int> punctured_indices, std::vector<float> sampling_points,
      std::vector<float> listsize_collect_sample_points);
  static std::vector<float> push_to_angle_boundary(std::vector<float> receivedMessage,
                                                   std::vector<int> transmittedCodeword,
                                                   std::vector<int> punctured_indices, float angle);
  std::vector<int> max_listsize_upto_metric;
  std::vector<float> average_listsize_upto_metric;
  std::vector<float> max_lower_envelop;

  // Linearity
  void generateNeighborList_sequential(const std::vector<float>& allZerosMessage, std::string dir,
                                       double thre);
  void generateNeighborList_sequential_TBonly(const std::vector<float>& allZerosMessage,
                                              double thre);
  MessageInformation lowrateDecoding_neighbors(const std::vector<float>& receivedMessage,
                                               std::vector<std::vector<int>> G_mat);
  bool isGabrielNeighbor(
      std::vector<std::vector<int>> G_mat,
      const std::vector<int>& modulated_codeword) const; //!< pass in by value is intentional
  void writeDataToFile(const std::vector<std::vector<int>>& data, const std::string& filename);

private:
  int numForwardPaths;
  int listSize;
  int crcLength;
  int crc;
  int nconv;
  char stopping_rule;

  std::vector<std::vector<int>> lowrate_nextStates;
  std::vector<std::vector<int>> lowrate_outputs;
  int lowrate_numStates;
  int lowrate_symbolLength;
  int lowrate_pathLength;

  struct cell {
    int optimalFatherState = -1;
    int suboptimalFatherState = -1;
    float pathMetric = std::numeric_limits<float>::max();
    float suboptimalPathMetric = std::numeric_limits<float>::max();
    bool init = false;
  };

  std::vector<int> pathToMessage(const std::vector<int>& path) const;
  std::vector<int> pathToCodeword(const std::vector<int>& path) const;
  std::vector<int> pathToMessage_ZT(std::vector<int> path);

  /* ZT */
  // construct the trellis using either product/euclidean distance metric
  std::vector<std::vector<cell>> constructLowRateTrellis_ZT(std::vector<float> receivedMessage,
                                                            METRIC_TYPE metric_type);

  /* TB */
  std::vector<std::vector<cell>> constructLowRateTrellis(const std::vector<float>& receivedMessage);
  std::vector<std::vector<cell>>
  constructLowRateTrellis_Punctured(const std::vector<float>& receivedMessage,
                                    const std::vector<int>& punctured_indices);
  std::vector<std::vector<cell>>
  constructLowRateTrellis_Punctured_ProductMetric(const std::vector<float>& receivedMessage,
                                                  const std::vector<int>& punctured_indices);
};

#endif
