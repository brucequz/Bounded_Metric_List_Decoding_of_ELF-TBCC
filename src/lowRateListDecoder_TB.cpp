#include "../include/lowRateListDecoder.h"
#include <algorithm>
#include <cstdlib>
#include <exception>
#include <map>
#include <vector>

LowRateListDecoder::LowRateListDecoder(FeedForwardTrellis feedforwardTrellis, int listSize, int crcLength, int crc, char stopping_rule) {
  this->lowrate_nextStates    = feedforwardTrellis.getNextStates();
	this->lowrate_outputs       = feedforwardTrellis.getOutputs();
	this->lowrate_numStates     = feedforwardTrellis.getNumStates();
	this->lowrate_symbolLength  = feedforwardTrellis.getN();
	this->numForwardPaths       = lowrate_nextStates[0].size();
  this->listSize              = listSize;
  this->crcLength             = crcLength;
  this->crc                   = crc;
	this->stopping_rule					= stopping_rule;

	this->max_listsize_upto_metric = std::vector<int>(100, -1);
	this->average_listsize_upto_metric = std::vector<float>(100, 0.0f);
	this->max_lower_envelop = std::vector<float>(100, -1.0f);


	if (this->stopping_rule != 'M' && this->stopping_rule != 'L' && this->stopping_rule != 'A' && this->stopping_rule != 'R') {
		std::cerr << "[ERROR] INVALID STOPPING RULE" << std::endl;
	}
}

MessageInformation LowRateListDecoder::decode(std::vector<float> receivedMessage, std::vector<int> punctured_indices) {
	/** Decode according to a policy passed into the constructor
	 * 
	 */

	if (ENCODING_RULE == 'Z' && STOPPING_RULE == 'M') {
		return lowRateDecoding_MaxMetric_EuclideanMetric_ZT(receivedMessage);
	} else if (ENCODING_RULE == 'Z' && STOPPING_RULE == 'A') {
		return lowRateDecoding_MaxAngle_ProductMetric_ZT(receivedMessage);
	} else if (ENCODING_RULE == 'Z' && STOPPING_RULE == 'L') {
		return lowRateDecoding_MaxListsize_ZT(receivedMessage);
	} 
	if (this->stopping_rule == 'L') {
		// max listsize restriction
		return lowRateDecoding_MaxListsize(receivedMessage, punctured_indices);
	} else if (this->stopping_rule == 'M') {
		// max metric restriction
		return lowRateDecoding_MaxMetric(receivedMessage, punctured_indices);
	} else if (ENCODING_RULE == 'T' && STOPPING_RULE == 'A') {
		// max angle restriction
		return lowRateDecoding_MaxAngle_ProductMetric(receivedMessage, punctured_indices);
	}
	throw std::invalid_argument("INVALID DECODING CHOICE!");
}

MessageInformation LowRateListDecoder::forceDecoding_MaxListsize(const std::vector<float>& receivedMessage, const std::vector<int>& punctured_indices, const std::vector<int>& codeword, const int listsize) {
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
	for(int i = 0; i < lowrate_numStates; i++){
		DetourObject detour;
		detour.startingState = i;
		detour.pathMetric = trellisInfo[i][lowrate_pathLength - 1].pathMetric;
		detourTree.insert(detour);
	}

	int numPathsSearched = 0;
	int numTBPathsSearched = 0;
  
	while(numPathsSearched < listsize){
		DetourObject detour = detourTree.pop();
		std::vector<int> path(lowrate_pathLength);

		int newTracebackStage = lowrate_pathLength - 1;
		float forwardPartialPathMetric = 0;
		int currentState = detour.startingState;

		// if we are taking a detour from a previous path, we skip backwards to the point where we take the
		// detour from the previous path
		if(detour.originalPathIndex != -1){
			forwardPartialPathMetric = detour.forwardPathMetric;
			newTracebackStage = detour.detourStage;

			// while we only need to copy the path from the detour to the end, this simplifies things,
			// and we'll write over the earlier data in any case
			path = previousPaths[detour.originalPathIndex];
			currentState = path[newTracebackStage];

			float suboptimalPathMetric = trellisInfo[currentState][newTracebackStage].suboptimalPathMetric;

			currentState = trellisInfo[currentState][newTracebackStage].suboptimalFatherState;
			newTracebackStage--;
			
			float prevPathMetric = trellisInfo[currentState][newTracebackStage].pathMetric;

			forwardPartialPathMetric += suboptimalPathMetric - prevPathMetric;
			
		}
		path[newTracebackStage] = currentState;

		// actually tracing back
		for(int stage = newTracebackStage; stage > 0; stage--) {
			float suboptimalPathMetric = trellisInfo[currentState][stage].suboptimalPathMetric;
			float currPathMetric = trellisInfo[currentState][stage].pathMetric;

			// if there is a detour we add to the detourTree
			if(trellisInfo[currentState][stage].suboptimalFatherState != -1){
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
			if (std::find(punctured_indices.begin(), punctured_indices.end(), i) == punctured_indices.end()) {
				hamming_dist += (int)cand_codeword[i]!=codeword[i];
				// euclidean_dist += std::pow(std::abs(cand_codeword[i]-codeword[i]),2);
			}
		}

		if (hamming_dist == 8) {
			break;
		}

		if (path[0] == path.back()) {
			
			coset_leaders[ED_int].push_back(message);
			numTBPathsSearched++;
			std::cout << "candidate codeword: " << numTBPathsSearched << std::endl;
			// std::cout << "outputing syndrome: " << syndrome << std::endl;
			// std::cout << "path[0]: " << path[0] << "; path[Kconv]: " << path[this->lowrate_pathLength-1] << "; ESD: " << ESD << std::endl;
			// std::cout << "message: "; utils::print_int_vector(message);
			// std::cout << "codeword: "; utils::print_int_vector(cand_codeword);
			// std::cout << "hamming dist: " << hamming_dist << "; euclidean dist: " << euclidean_dist << "; forwardPathMetric: " << forwardPartialPathMetric << std::endl;
			// std::cout << std::endl;
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

MessageInformation LowRateListDecoder::lowRateDecoding_MaxListsize(const std::vector<float>& receivedMessage, const std::vector<int>& punctured_indices) {
	// trellisInfo is indexed [state][stage]
	std::vector<std::vector<cell>> trellisInfo;
	trellisInfo = constructLowRateTrellis_Punctured(receivedMessage, punctured_indices);

	// start search
	MessageInformation output;
	//RBTree detourTree;
	MinHeap detourTree;
	std::vector<std::vector<int>> previousPaths;
	

	// create nodes for each valid ending state with no detours
	// std::cout<< "end path metrics:" <<std::endl;
	for(int i = 0; i < lowrate_numStates; i++){
		DetourObject detour;
		detour.startingState = i;
		detour.pathMetric = trellisInfo[i][lowrate_pathLength - 1].pathMetric;
		detourTree.insert(detour);
	}

	int numPathsSearched = 0;
	int TBPathsSearched = 0;
  
	while(numPathsSearched < this->listSize){
		DetourObject detour = detourTree.pop();
		std::vector<int> path(lowrate_pathLength);

		int newTracebackStage = lowrate_pathLength - 1;
		float forwardPartialPathMetric = 0;
		int currentState = detour.startingState;

		// if we are taking a detour from a previous path, we skip backwards to the point where we take the
		// detour from the previous path
		if(detour.originalPathIndex != -1){
			forwardPartialPathMetric = detour.forwardPathMetric;
			newTracebackStage = detour.detourStage;

			// while we only need to copy the path from the detour to the end, this simplifies things,
			// and we'll write over the earlier data in any case
			path = previousPaths[detour.originalPathIndex];
			currentState = path[newTracebackStage];

			float suboptimalPathMetric = trellisInfo[currentState][newTracebackStage].suboptimalPathMetric;

			currentState = trellisInfo[currentState][newTracebackStage].suboptimalFatherState;
			newTracebackStage--;
			
			float prevPathMetric = trellisInfo[currentState][newTracebackStage].pathMetric;

			forwardPartialPathMetric += suboptimalPathMetric - prevPathMetric;
			
		}
		path[newTracebackStage] = currentState;

		// actually tracing back
		for(int stage = newTracebackStage; stage > 0; stage--){
			float suboptimalPathMetric = trellisInfo[currentState][stage].suboptimalPathMetric;
			float currPathMetric = trellisInfo[currentState][stage].pathMetric;

			// if there is a detour we add to the detourTree
			if(trellisInfo[currentState][stage].suboptimalFatherState != -1){
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
		}
		
		previousPaths.push_back(path);

		std::vector<int> message = pathToMessage(path);
		std::vector<int> codeword = pathToCodeword(path);

		// std::cout << "forwardPartialPathMetric: " << forwardPartialPathMetric << std::endl;
		
		// one trellis decoding requires both a tb and crc check
		if(path[0] == path[lowrate_pathLength - 1] && crc::crc_check(message, crcLength, crc) && numPathsSearched <= this->listSize){
			output.message = message;
			output.codeword = codeword;
			output.path = path;
		 	output.listSize = numPathsSearched + 1;
			output.metric = forwardPartialPathMetric;
			output.TBListSize = TBPathsSearched + 1;
		 	return output;
		}

		numPathsSearched++;
		if(path[0] == path[lowrate_pathLength - 1])
			TBPathsSearched++;
	} // while(numPathsSearched < this->listSize)

	output.listSizeExceeded = true;
	std::cout << "listsize exceeded!!!" << std::endl;
	return output;
}



MessageInformation LowRateListDecoder::lowRateDecoding_MaxMetric(std::vector<float> receivedMessage, std::vector<int> punctured_indices){
	// trellisInfo is indexed [state][stage]
	std::vector<std::vector<cell>> trellisInfo;
	trellisInfo = constructLowRateTrellis_Punctured(receivedMessage, punctured_indices);

	// start search
	MessageInformation output;
	//RBTree detourTree;
	MinHeap detourTree;
	std::vector<std::vector<int>> previousPaths;
	

	// create nodes for each valid ending state with no detours
	// std::cout<< "end path metrics:" <<std::endl;
	for(int i = 0; i < lowrate_numStates; i++){
		DetourObject detour;
		detour.startingState = i;
		detour.pathMetric = trellisInfo[i][lowrate_pathLength - 1].pathMetric;
		detourTree.insert(detour);
	}

	int numPathsSearched = 0;
	int TBPathsSearched = 0;
	float currentMetricExplored = 0.0;
  
	while(currentMetricExplored < MAX_METRIC){
		DetourObject detour = detourTree.pop();
		std::vector<int> path(lowrate_pathLength);

		int newTracebackStage = lowrate_pathLength - 1;
		float forwardPartialPathMetric = 0;
		int currentState = detour.startingState;

		// if we are taking a detour from a previous path, we skip backwards to the point where we take the
		// detour from the previous path
		if(detour.originalPathIndex != -1){
			forwardPartialPathMetric = detour.forwardPathMetric;
			newTracebackStage = detour.detourStage;

			// while we only need to copy the path from the detour to the end, this simplifies things,
			// and we'll write over the earlier data in any case
			path = previousPaths[detour.originalPathIndex];
			currentState = path[newTracebackStage];

			float suboptimalPathMetric = trellisInfo[currentState][newTracebackStage].suboptimalPathMetric;

			currentState = trellisInfo[currentState][newTracebackStage].suboptimalFatherState;
			newTracebackStage--;
			
			float prevPathMetric = trellisInfo[currentState][newTracebackStage].pathMetric;

			forwardPartialPathMetric += suboptimalPathMetric - prevPathMetric;
			
		}
		path[newTracebackStage] = currentState;

		// actually tracing back
		for(int stage = newTracebackStage; stage > 0; stage--){
			float suboptimalPathMetric = trellisInfo[currentState][stage].suboptimalPathMetric;
			float currPathMetric = trellisInfo[currentState][stage].pathMetric;

			// if there is a detour we add to the detourTree
			if(trellisInfo[currentState][stage].suboptimalFatherState != -1){
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
		std::vector<int> codeword = pathToCodeword(path);
		currentMetricExplored = forwardPartialPathMetric;
		
		// one trellis decoding requires both a tb and crc check
		if(path[0] == path[lowrate_pathLength - 1] && crc::crc_check(message, crcLength, crc) && currentMetricExplored <= MAX_METRIC){
			output.message = message;
			output.path = path;
		 	output.listSize = numPathsSearched + 1;
			output.metric = forwardPartialPathMetric;
			output.TBListSize = TBPathsSearched + 1;
		 	return output;
		}

		numPathsSearched++;
		if(path[0] == path[lowrate_pathLength - 1])
			TBPathsSearched++;
	} // while(currentMetricExplored < MAX_METRIC)

	output.listSizeExceeded = true;
	// std::cerr << "[WARNING]: TC IS NOT FOUND!!! " << std::endl;
	return output;
}

MessageInformation LowRateListDecoder::lowRateDecoding_MaxAngle(std::vector<float> receivedMessage, std::vector<int> punctured_indices) {
	// trellisInfo is indexed [state][stage]
	std::vector<std::vector<cell>> trellisInfo;
	trellisInfo = constructLowRateTrellis_Punctured(receivedMessage, punctured_indices);

	// start search
	MessageInformation output;
	//RBTree detourTree;
	MinHeap detourTree;
	std::vector<std::vector<int>> previousPaths;
	

	// create nodes for each valid ending state with no detours
	// std::cout<< "end path metrics:" <<std::endl;
	for(int i = 0; i < lowrate_numStates; i++){
		DetourObject detour;
		detour.startingState = i;
		detour.pathMetric = trellisInfo[i][lowrate_pathLength - 1].pathMetric;
		detourTree.insert(detour);
	}

	int numPathsSearched = 0;
	int TBPathsSearched = 0;
	float currentAngleExplored = 0.0;
	
	while(currentAngleExplored < MAX_ANGLE){
		DetourObject detour = detourTree.pop();
		std::vector<int> path(lowrate_pathLength);

		int newTracebackStage = lowrate_pathLength - 1;
		float forwardPartialPathMetric = 0;
		int currentState = detour.startingState;

		// if we are taking a detour from a previous path, we skip backwards to the point where we take the
		// detour from the previous path
		if(detour.originalPathIndex != -1){
			forwardPartialPathMetric = detour.forwardPathMetric;
			newTracebackStage = detour.detourStage;

			// while we only need to copy the path from the detour to the end, this simplifies things,
			// and we'll write over the earlier data in any case
			path = previousPaths[detour.originalPathIndex];
			currentState = path[newTracebackStage];

			float suboptimalPathMetric = trellisInfo[currentState][newTracebackStage].suboptimalPathMetric;

			currentState = trellisInfo[currentState][newTracebackStage].suboptimalFatherState;
			newTracebackStage--;
			
			float prevPathMetric = trellisInfo[currentState][newTracebackStage].pathMetric;

			forwardPartialPathMetric += suboptimalPathMetric - prevPathMetric;
			
		}
		path[newTracebackStage] = currentState;

		// actually tracing back
		for(int stage = newTracebackStage; stage > 0; stage--){
			float suboptimalPathMetric = trellisInfo[currentState][stage].suboptimalPathMetric;
			float currPathMetric = trellisInfo[currentState][stage].pathMetric;

			// if there is a detour we add to the detourTree
			if(trellisInfo[currentState][stage].suboptimalFatherState != -1){
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
		std::vector<int> codeword = pathToCodeword(path);
		std::cout << "[DEBUG!] receivedMessage size: " << receivedMessage.size() << "codeword size: " << codeword.size() << std::endl;
		currentAngleExplored = utils::compute_angle_between_vectors_rad(receivedMessage, codeword, punctured_indices);
		
		// one trellis decoding requires both a tb and crc check
		if(path[0] == path[lowrate_pathLength - 1] && crc::crc_check(message, crcLength, crc) && currentAngleExplored < MAX_ANGLE){
			output.message = message;
			output.path = path;
			output.listSize = numPathsSearched + 1;
			output.metric = forwardPartialPathMetric;
			output.TBListSize = TBPathsSearched + 1;
			output.angle_received_decoded_rad = currentAngleExplored;
			return output;
		}

		numPathsSearched++;
		if(path[0] == path[lowrate_pathLength - 1])
			TBPathsSearched++;
	} // while(currentAngleExplored < MAX_ANGLE)

	output.listSizeExceeded = true;
	output.listSize = numPathsSearched;
	// std::cerr << "[WARNING]: TC IS NOT FOUND!!! " << std::endl;
	return output;
}

MessageInformation LowRateListDecoder::lowRateDecoding_MaxAngle_ProductMetric(std::vector<float> receivedMessage, std::vector<int> punctured_indices){
	std::vector<std::vector<cell>> trellisInfo;

	trellisInfo = constructLowRateTrellis_Punctured_ProductMetric(receivedMessage, punctured_indices);
	// std::cout << "printing the first row of fixedp pathMetric." << std::endl;
	// for (size_t i = 0; i < trellisInfo.size(); i++) {
	// 	std::cout << trellisInfo[i].back().pathMetric << ", ";
	// }
	// std::cout << std::endl;

	// start search
	MessageInformation output;
	//RBTree detourTree;
	MinHeap detourTree;
	std::vector<std::vector<int>> previousPaths;
	

	// create nodes for each valid ending state with no detours
	for(int i = 0; i < lowrate_numStates; i++){
		DetourObject detour;
		detour.startingState = i;
		detour.pathMetric = trellisInfo[i][lowrate_pathLength - 1].pathMetric;
		detourTree.insert(detour);
	}

	int numPathsSearched = 0;
	int TBPathsSearched = 0;
	float currentAngleExplored = 0.0;
	
	while(currentAngleExplored < MAX_ANGLE){
		DetourObject detour = detourTree.pop();
		// std::cout << "floatp detour tree item: " << detour.pathMetric << std::endl;
		std::vector<int> path(lowrate_pathLength);

		int newTracebackStage = lowrate_pathLength - 1;
		float forwardPartialPathMetric = 0;
		int currentState = detour.startingState;

		// if we are taking a detour from a previous path, we skip backwards to the point where we take the
		// detour from the previous path
		if(detour.originalPathIndex != -1){
			forwardPartialPathMetric = detour.forwardPathMetric;
			newTracebackStage = detour.detourStage;

			// while we only need to copy the path from the detour to the end, this simplifies things,
			// and we'll write over the earlier data in any case
			path = previousPaths[detour.originalPathIndex];
			currentState = path[newTracebackStage];

			float suboptimalPathMetric = trellisInfo[currentState][newTracebackStage].suboptimalPathMetric;

			currentState = trellisInfo[currentState][newTracebackStage].suboptimalFatherState;
			newTracebackStage--;
			
			float prevPathMetric = trellisInfo[currentState][newTracebackStage].pathMetric;

			forwardPartialPathMetric += suboptimalPathMetric - prevPathMetric;
			
		}
		path[newTracebackStage] = currentState;

		// actually tracing back
		for(int stage = newTracebackStage; stage > 0; stage--){
			float suboptimalPathMetric = trellisInfo[currentState][stage].suboptimalPathMetric;
			float currPathMetric = trellisInfo[currentState][stage].pathMetric;

			// if there is a detour we add to the detourTree
			if(trellisInfo[currentState][stage].suboptimalFatherState != -1){
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
		std::vector<int> codeword = pathToCodeword(path);

		// std::cout << "floatp forward partial path metric: " << forwardPartialPathMetric << std::endl;
		// std::cout << "printing codeword: " << std::endl;
		// utils::print_int_vector(codeword);

		

		// another way to compute the angle
		// currentAngleExplored = utils::compute_angle_between_vectors_rad(receivedMessage, codeword);
		// std::cout << "printing -forwardPartialPathMetric/N: " << -forwardPartialPathMetric/N << std::endl;
		currentAngleExplored = std::acos( std::max(-1.0f, std::min(1.0f, -forwardPartialPathMetric/N)) );
		// std::cout << "floatp current Angle: " << currentAngleExplored << std::endl;
		// std::cout << "max angle: " << MAX_ANGLE << std::endl;
		
		// one trellis decoding requires both a tb and crc check
		if(path[0] == path[lowrate_pathLength - 1] && crc::crc_check(message, crcLength, crc) && currentAngleExplored <= MAX_ANGLE){
			output.message = message;
			output.path = path;
			output.listSize = numPathsSearched + 1;
			output.metric = forwardPartialPathMetric;
			output.TBListSize = TBPathsSearched + 1;
			output.angle_received_decoded_rad = currentAngleExplored;
			// std::cout << "list size: " << numPathsSearched << ", float p returning angle decoded = " << currentAngleExplored << ", metric = " << forwardPartialPathMetric << std::endl;
			
			return output;
		}

		numPathsSearched++;
		if(path[0] == path[lowrate_pathLength - 1])
			TBPathsSearched++;
	} // while(currentAngleExplored < MAX_ANGLE)

	output.listSizeExceeded = true;
	output.listSize = numPathsSearched;
	// std::cerr << "[WARNING]: TC IS NOT FOUND!!! " << std::endl;
	return output;
}

std::vector<std::vector<LowRateListDecoder::cell>> LowRateListDecoder::constructLowRateTrellis(const std::vector<float>& receivedMessage){
	std::vector<std::vector<cell>> trellisInfo;
	lowrate_pathLength = (receivedMessage.size() / lowrate_symbolLength) + 1;

	trellisInfo = std::vector<std::vector<cell>>(lowrate_numStates, std::vector<cell>(lowrate_pathLength));

	// initializes all the valid starting states
	for(int i = 0; i < lowrate_numStates; i++){
		trellisInfo[i][0].pathMetric = 0;
		trellisInfo[i][0].init = true;
	}
	
	// building the trellis
	for(int stage = 0; stage < lowrate_pathLength - 1; stage++){
		for(int currentState = 0; currentState < lowrate_numStates; currentState++){
			// if the state / stage is invalid, we move on
			if(!trellisInfo[currentState][stage].init)
				continue;

			// otherwise, we compute the relevent information
			for(int forwardPathIndex = 0; forwardPathIndex < numForwardPaths; forwardPathIndex++){
				// since our transitions correspond to symbols, the forwardPathIndex has no correlation 
				// beyond indexing the forward path

				int nextState = lowrate_nextStates[currentState][forwardPathIndex];
				
				// if the nextState is invalid, we move on
				if(nextState < 0)
					continue;
				
				float branchMetric = 0;
				std::vector<int> output_point = crc::get_point(lowrate_outputs[currentState][forwardPathIndex], lowrate_symbolLength);
				
				for(int i = 0; i < lowrate_symbolLength; i++){
					branchMetric += std::pow(receivedMessage[lowrate_symbolLength * stage + i] - (float)output_point[i], 2);
					// branchMetric += std::abs(receivedMessage[lowrate_symbolLength * stage + i] - (float)output_point[i]);
				}
				float totalPathMetric = branchMetric + trellisInfo[currentState][stage].pathMetric;
				
				// dealing with cases of uninitialized states, when the transition becomes the optimal father state, and suboptimal father state, in order
				if(!trellisInfo[nextState][stage + 1].init){
					trellisInfo[nextState][stage + 1].pathMetric = totalPathMetric;
					trellisInfo[nextState][stage + 1].optimalFatherState = currentState;
					trellisInfo[nextState][stage + 1].init = true;
				}
				else if(trellisInfo[nextState][stage + 1].pathMetric > totalPathMetric){
					trellisInfo[nextState][stage + 1].suboptimalPathMetric = trellisInfo[nextState][stage + 1].pathMetric;
					trellisInfo[nextState][stage + 1].suboptimalFatherState = trellisInfo[nextState][stage + 1].optimalFatherState;
					trellisInfo[nextState][stage + 1].pathMetric = totalPathMetric;
					trellisInfo[nextState][stage + 1].optimalFatherState = currentState;
				}
				else{
					trellisInfo[nextState][stage + 1].suboptimalPathMetric = totalPathMetric;
					trellisInfo[nextState][stage + 1].suboptimalFatherState = currentState;
				}
			}

		}
	}
	return trellisInfo;
}

std::vector<std::vector<LowRateListDecoder::cell>> LowRateListDecoder::constructLowRateTrellis_Punctured(const std::vector<float>& receivedMessage, const std::vector<int>& punctured_indices) {
	/* Constructs a trellis for a low rate code, with puncturing
		Args:
			receivedMessage (std::vector<float>): the received message
			punctured_indices (std::vector<int>): the indices of the punctured bits

		Returns:
			std::vector<std::vector<cell>>: the trellis
	*/

	/* ---- Code Begins ---- */
	std::vector<std::vector<cell>> trellisInfo;
	lowrate_pathLength = (receivedMessage.size() / lowrate_symbolLength) + 1;

	trellisInfo = std::vector<std::vector<cell>>(lowrate_numStates, std::vector<cell>(lowrate_pathLength));

	// initializes all the valid starting states
	for(int i = 0; i < lowrate_numStates; i++){
		trellisInfo[i][0].pathMetric = 0;
		trellisInfo[i][0].init = true;
	}
	
	// building the trellis
	for(int stage = 0; stage < lowrate_pathLength - 1; stage++){
		for(int currentState = 0; currentState < lowrate_numStates; currentState++){
			// if the state / stage is invalid, we move on
			if(!trellisInfo[currentState][stage].init)
				continue;

			// otherwise, we compute the relevent information
			for(int forwardPathIndex = 0; forwardPathIndex < numForwardPaths; forwardPathIndex++){
				// since our transitions correspond to symbols, the forwardPathIndex has no correlation 
				// beyond indexing the forward path

				int nextState = lowrate_nextStates[currentState][forwardPathIndex];
				
				// if the nextState is invalid, we move on
				if(nextState < 0)
					continue;
				
				float branchMetric = 0;
				std::vector<int> output_point = crc::get_point(lowrate_outputs[currentState][forwardPathIndex], lowrate_symbolLength);
				
				for(int i = 0; i < lowrate_symbolLength; i++){
					if (std::find(punctured_indices.begin(), punctured_indices.end(), lowrate_symbolLength * stage + i) != punctured_indices.end()){
						branchMetric += 0;
					} else {
						branchMetric += std::pow(receivedMessage[lowrate_symbolLength * stage + i] - (float)output_point[i], 2);
					}
				}
				
				float totalPathMetric = branchMetric + trellisInfo[currentState][stage].pathMetric;
				
				// dealing with cases of uninitialized states, when the transition becomes the optimal father state, and suboptimal father state, in order
				if(!trellisInfo[nextState][stage + 1].init){
					trellisInfo[nextState][stage + 1].pathMetric = totalPathMetric;
					trellisInfo[nextState][stage + 1].optimalFatherState = currentState;
					trellisInfo[nextState][stage + 1].init = true;
				}
				else if(trellisInfo[nextState][stage + 1].pathMetric > totalPathMetric){
					trellisInfo[nextState][stage + 1].suboptimalPathMetric = trellisInfo[nextState][stage + 1].pathMetric;
					trellisInfo[nextState][stage + 1].suboptimalFatherState = trellisInfo[nextState][stage + 1].optimalFatherState;
					trellisInfo[nextState][stage + 1].pathMetric = totalPathMetric;
					trellisInfo[nextState][stage + 1].optimalFatherState = currentState;
				}
				else{
					trellisInfo[nextState][stage + 1].suboptimalPathMetric = totalPathMetric;
					trellisInfo[nextState][stage + 1].suboptimalFatherState = currentState;
				}
			}

		}
	}
	return trellisInfo;
}

std::vector<std::vector<LowRateListDecoder::cell>> LowRateListDecoder::constructLowRateTrellis_Punctured_ProductMetric(const std::vector<float>& receivedMessage, const std::vector<int>& punctured_indices){
	/* Constructs a trellis for a low rate code, with puncturing
		Args:
			receivedMessage (std::vector<float>): the received message
			punctured_indices (std::vector<int>): the indices of the punctured bits

		Returns:
			std::vector<std::vector<cell>>: the trellis
	*/

	/* ---- Code Begins ---- */
	std::vector<std::vector<cell>> trellisInfo;
	lowrate_pathLength = (receivedMessage.size() / lowrate_symbolLength) + 1;

	trellisInfo = std::vector<std::vector<cell>>(lowrate_numStates, std::vector<cell>(lowrate_pathLength));

	// initializes all the valid starting states
	for(int i = 0; i < lowrate_numStates; i++){
		// std::cout << "before setting to 0: " << trellisInfo[i][0].pathMetric << std::endl;
		trellisInfo[i][0].pathMetric = 0;
		trellisInfo[i][0].init = true;
	}
	
	// building the trellis
	for(int stage = 0; stage < lowrate_pathLength - 1; stage++){
		for(int currentState = 0; currentState < lowrate_numStates; currentState++){
			// if the state / stage is invalid, we move on
			if(!trellisInfo[currentState][stage].init)
				continue;

			// otherwise, we compute the relevent information
			for(int forwardPathIndex = 0; forwardPathIndex < numForwardPaths; forwardPathIndex++){
				// since our transitions correspond to symbols, the forwardPathIndex has no correlation 
				// beyond indexing the forward path

				int nextState = lowrate_nextStates[currentState][forwardPathIndex];
				
				// if the nextState is invalid, we move on
				if(nextState < 0)
					continue;
				
				float branchMetric = 0;
				std::vector<int> output_point = crc::get_point(lowrate_outputs[currentState][forwardPathIndex], lowrate_symbolLength);
				
				for(int i = 0; i < lowrate_symbolLength; i++){
					if (std::find(punctured_indices.begin(), punctured_indices.end(), lowrate_symbolLength * stage + i) != punctured_indices.end()){
						branchMetric += 0;
					} else {
						branchMetric += -receivedMessage[lowrate_symbolLength* stage + i] * (float)output_point[i];
					}
				}
				
				float totalPathMetric = branchMetric + trellisInfo[currentState][stage].pathMetric;
				
				// dealing with cases of uninitialized states, when the transition becomes the optimal father state, and suboptimal father state, in order
				if(!trellisInfo[nextState][stage + 1].init){
					trellisInfo[nextState][stage + 1].pathMetric = totalPathMetric;
					trellisInfo[nextState][stage + 1].optimalFatherState = currentState;
					trellisInfo[nextState][stage + 1].init = true;
				}
				else if(trellisInfo[nextState][stage + 1].pathMetric > totalPathMetric){
					trellisInfo[nextState][stage + 1].suboptimalPathMetric = trellisInfo[nextState][stage + 1].pathMetric;
					trellisInfo[nextState][stage + 1].suboptimalFatherState = trellisInfo[nextState][stage + 1].optimalFatherState;
					trellisInfo[nextState][stage + 1].pathMetric = totalPathMetric;
					trellisInfo[nextState][stage + 1].optimalFatherState = currentState;
				}
				else{
					trellisInfo[nextState][stage + 1].suboptimalPathMetric = totalPathMetric;
					trellisInfo[nextState][stage + 1].suboptimalFatherState = currentState;
				}
			}
		}
	}
	return trellisInfo;
}



// converts a path through the tb trellis to the binary message it corresponds with
std::vector<int> LowRateListDecoder::pathToMessage(const std::vector<int>& path) const {
	std::vector<int> message;
	for(size_t pathIndex = 0; pathIndex < path.size() - 1; pathIndex++){
		for(int forwardPath = 0; forwardPath < numForwardPaths; forwardPath++){
			if(lowrate_nextStates[path[pathIndex]][forwardPath] == path[pathIndex + 1])
				message.push_back(forwardPath);
		}
	}
	return message;
}

// converts a path through the tb trellis to the BPSK it corresponds with
// currently does NOT puncture the codeword
std::vector<int> LowRateListDecoder::pathToCodeword(const std::vector<int>& path) const {
	std::vector<int> nopunc_codeword;
	for(size_t pathIndex = 0; pathIndex < path.size() - 1; pathIndex++){
		for(int forwardPath = 0; forwardPath < numForwardPaths; forwardPath++){
			if(lowrate_nextStates[path[pathIndex]][forwardPath] == path[pathIndex + 1]){
				std::vector<int> output_bin;
				crc::dec_to_binary(lowrate_outputs[path[pathIndex]][forwardPath], output_bin, lowrate_symbolLength);
				for (int outbit = 0; outbit < lowrate_symbolLength; outbit ++){
					nopunc_codeword.push_back(-2 * output_bin[outbit] + 1);
				}
			}
		}
	}

	return nopunc_codeword;
}