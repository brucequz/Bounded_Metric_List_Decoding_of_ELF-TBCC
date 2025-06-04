
#ifndef LOWRATELISTDECODER_H
#define LOWRATELISTDECODER_H

#include <climits>
#include <iostream>
#include <algorithm>

#include "feedForwardTrellis.h"
#include "minHeap.h"
#include "types.h"
#include "namespace.h"
#include "../consts.h"

class LowRateListDecoder{
public:
	LowRateListDecoder(FeedForwardTrellis FT, int listSize, int crcDegree, int crc, char stopping_rule);

	/* - Floating Point - */
	MessageInformation decode(std::vector<float> receivedMessage, std::vector<int> punctured_indices, float sigma_sqrd, float rova_t);
	MessageInformation lowRateDecoding_MaxListsize(std::vector<float> receivedMessage, std::vector<int> punctured_indices);
	MessageInformation lowRateDecoding_MaxMetric(std::vector<float> receivedMessage, std::vector<int> punctured_indices);
	MessageInformation lowRateDecoding_MaxAngle(std::vector<float> receivedMessage, std::vector<int> punctured_indices);
	MessageInformation lowRateDecoding_MaxAngle_ProductMetric(std::vector<float> receivedMessage, std::vector<int> punctured_indices);

	// ZT
	MessageInformation lowRateDecoding_MaxAngle_ProductMetric_ZT(std::vector<float> receivedMessage);
	MessageInformation lowRateDecoding_SquaredDistanceMetric_ROVA_ZT(std::vector<float> receivedMessage, float sigma_sqrd, float rova_t);

private:
	int numForwardPaths;
	int listSize;
	int crcDegree;
	int crc;
	int n;
	char stopping_rule;

	std::vector<std::vector<int>> lowrate_nextStates;
	std::vector<std::vector<int>> lowrate_outputs;
	int lowrate_numStates;
	int lowrate_symbolLength;
	int lowrate_pathLength;

	// ROVA
	std::vector<std::vector<std::vector<float>>> log_gammas_;

	struct cell {
		int optimalFatherState = -1;
		int suboptimalFatherState = -1;
		float pathMetric = INT_MAX;
		float suboptimalPathMetric = INT_MAX;
		bool init = false;
	};

	struct rova_cell {
		// ROVA
		float log_Gamma = -INFINITY;
		float log_Z = -INFINITY;
		// regular cell
		int optimalFatherState = -1;
		int suboptimalFatherState = -1;
		float pathMetric = INT_MAX;
		float suboptimalPathMetric = INT_MAX;
		bool init = false;
	};

	std::vector<int> pathToMessage(std::vector<int>); 
  std::vector<int> pathToCodeword(std::vector<int>); 
	std::vector<int> pathToMessage_ZT(std::vector<int> path);

	/* - Floating Point - */
	std::vector<std::vector<cell>> constructLowRateTrellis(std::vector<float> receivedMessage);

	// ZT
	std::vector<std::vector<cell>> constructLowRateTrellis_ZT(std::vector<float> receivedMessage);
	// ZT ROVA
	std::vector<std::vector<rova_cell>> constructLowRateTrellis_ROVA_Alg2_ZT(std::vector<float> receivedMessage);
	std::vector<std::vector<rova_cell>> constructLowRateTrellis_ROVA_Alg4_ZT(std::vector<float> receivedMessage, float sigma_sqrd);
	float compute_logGamma(std::vector<float> receivedMessage, std::vector<int> codeword, float sigma_sqrd);
	// computes max star approximation in BCJR decoding
	float max_star(float lnx, float lny);

	// TB Punctured
  std::vector<std::vector<cell>> constructLowRateTrellis_Punctured(std::vector<float> receivedMessage, std::vector<int> punctured_indices);
	std::vector<std::vector<cell>> constructLowRateTrellis_Punctured_ProductMetric(std::vector<float> receivedMessage, std::vector<int> punctured_indices);
};


#endif
