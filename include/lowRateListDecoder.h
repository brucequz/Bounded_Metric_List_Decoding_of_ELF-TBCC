#ifndef LOWRATELISTDECODER_H
#define LOWRATELISTDECODER_H

#include <climits>
#include <cmath>
#include <iostream>
#include <algorithm>

#include "feedForwardTrellis.h"
#include "minHeap.h"
#include "types.h"

enum class METRIC_TYPE {
	PRODUCT_METRIC = 1,
	EUCLIDEAN_METRIC = 2,
	LOG_EUCLIDEAN_METRIC = 3
};

class LowRateListDecoder{
public:
	LowRateListDecoder(FeedForwardTrellis FT, int listSize, int crcDegree, int crc, char stopping_rule);

	MessageInformation decode(std::vector<float> receivedMessage, std::vector<int> punctured_indices);

	// TB
	MessageInformation lowRateDecoding_MaxListsize(std::vector<float> receivedMessage, std::vector<int> punctured_indices);
	MessageInformation lowRateDecoding_MaxMetric(std::vector<float> receivedMessage, std::vector<int> punctured_indices);
	MessageInformation lowRateDecoding_MaxAngle(std::vector<float> receivedMessage, std::vector<int> punctured_indices);
	MessageInformation lowRateDecoding_MaxAngle_ProductMetric(std::vector<float> receivedMessage, std::vector<int> punctured_indices);

	// ZT
	MessageInformation lowRateDecoding_MaxAngle_ProductMetric_ZT(std::vector<float> receivedMessage);
	MessageInformation lowRateDecoding_MaxMetric_EuclideanMetric_ZT(std::vector<float> receivedMessage);

	// BCJR
	MessageInformation lowRateDecoding_BCJR(std::vector<float> receivedMessage, float sigma_sqrd);

	// ROVA
	MessageInformation decode_ROVA(std::vector<float> receivedMessage, std::vector<int> punctured_indices, float sigma_sqrd, float rova_t);
	MessageInformation lowRateDecoding_SquaredDistanceMetric_ROVA_ZT(std::vector<float> receivedMessage, float sigma_sqrd, float rova_t);

	// genie-Aided
	MessageInformation genieAided_LowRateDecoding_MaxAngle_ProductMetric(std::vector<float> receivedMessage, std::vector<int> transmittedCodeword, std::vector<int> punctured_indices, bool push_to_boundary);
	MessageInformation genieAided_LowRateDecoding_MaxListsize_Collect(std::vector<float> receivedMessage, std::vector<int> punctured_indices, std::vector<int> transmittedCodeword, std::vector<float> sampling_points, std::vector<float> listsize_collect_sample_points);
	MessageInformation genieAided_LowRateDecoding_MaxAngle_ProductMetric_Collector(std::vector<float> receivedMessage, std::vector<int> transmittedCodeword, std::vector<int> punctured_indices, std::vector<float> sampling_points, std::vector<float> listsize_collect_sample_points);
	static std::vector<float> push_to_angle_boundary(std::vector<float> receivedMessage, std::vector<int> transmittedCodeword, std::vector<int> punctured_indices, float angle);
	std::vector<int> max_listsize_upto_metric;
	std::vector<float> average_listsize_upto_metric;
	std::vector<float> max_lower_envelop;

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

	struct rova_cell : public cell {
		// ROVA
		float log_Gamma = -INFINITY;
		float log_Z = -INFINITY;
	};

	struct bcjr_cell : public cell {
		// BCJR
		float _log_alpha 	= -INFINITY;
		float _log_beta 	= -INFINITY;
	};

	std::vector<int> pathToMessage(std::vector<int>); 
  std::vector<int> pathToCodeword(std::vector<int>); 
	std::vector<int> pathToMessage_ZT(std::vector<int> path);

	/* - Floating Point - */
	std::vector<std::vector<cell>> constructLowRateTrellis(std::vector<float> receivedMessage);

	// ZT
	// constructs the trellis with rova information
	std::vector<std::vector<rova_cell>> constructLowRateTrellis_ROVA_Alg4_ZT(std::vector<float> receivedMessage, float sigma_sqrd);
	// construct the trellis using either product/euclidean distance metric
	std::vector<std::vector<cell>> constructLowRateTrellis_ZT(std::vector<float> receivedMessage, METRIC_TYPE metric_type);
	// computes the path likelihood Gamma in logarithm
	float compute_logGamma(std::vector<float> receivedMessage, std::vector<int> codeword, float sigma_sqrd);
	// computes max star approximation in BCJR decoding
	float max_star(float lnx, float lny);

	// TB Punctured
  std::vector<std::vector<cell>> constructLowRateTrellis_Punctured(std::vector<float> receivedMessage, std::vector<int> punctured_indices);
	std::vector<std::vector<cell>> constructLowRateTrellis_Punctured_ProductMetric(std::vector<float> receivedMessage, std::vector<int> punctured_indices);

	// BCJR
	void _bcjr_forward_pass(std::vector<std::vector<bcjr_cell>> &trellis);
	void _bcjr_backward_pass(std::vector<std::vector<bcjr_cell>> &trellis);
	void _bcjr_log_gamma(std::vector<float> receivedMessage, float sigma_sqrd);
};


#endif
