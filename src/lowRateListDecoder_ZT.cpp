#include "lowRateListDecoder.h"

// construct ZT trellis
std::vector<std::vector<LowRateListDecoder::cell>> LowRateListDecoder::constructLowRateTrellis_ZT(std::vector<float> receivedMessage, METRIC_TYPE metric_type){
	std::vector<std::vector<cell>> trellisInfo;
	lowrate_pathLength = (receivedMessage.size() / lowrate_symbolLength) + 1;

	trellisInfo = std::vector<std::vector<cell>>(lowrate_numStates, std::vector<cell>(lowrate_pathLength));

	// initialize only 0 as the starting states
	trellisInfo[0][0].pathMetric = 0;
	trellisInfo[0][0].init = true;
	
	// building the trellis
	for(int stage = 0; stage < lowrate_pathLength - V - 1; stage++){
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
					if (metric_type == METRIC_TYPE::PRODUCT_METRIC) {
						branchMetric += -1 * (receivedMessage[lowrate_symbolLength * stage + i] * (float)output_point[i]);
					} else if (metric_type == METRIC_TYPE::EUCLIDEAN_METRIC) {
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
	// ZT stage
	for(int stage = lowrate_pathLength - V - 1; stage < lowrate_pathLength - 1; stage++){
		for(int currentState = 0; currentState < lowrate_numStates; currentState++){
			// if the state / stage is invalid, we move on
			if(!trellisInfo[currentState][stage].init)
				continue;

			// zero terminating
			int forwardPathIndex = 0;
			// since our transitions correspond to symbols, the forwardPathIndex has no correlation 
			// beyond indexing the forward path

			int nextState = lowrate_nextStates[currentState][forwardPathIndex];
			
			// if the nextState is invalid, we move on
			if(nextState < 0)
				continue;
			
			float branchMetric = 0;
			std::vector<int> output_point = crc::get_point(lowrate_outputs[currentState][forwardPathIndex], lowrate_symbolLength);
			
			for(int i = 0; i < lowrate_symbolLength; i++){
				if (metric_type == METRIC_TYPE::PRODUCT_METRIC) {
					branchMetric += -1 * (receivedMessage[lowrate_symbolLength * stage + i] * (float)output_point[i]);
				} else if (metric_type == METRIC_TYPE::EUCLIDEAN_METRIC) {
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
	return trellisInfo;
}

MessageInformation LowRateListDecoder::lowRateDecoding_MaxAngle_ProductMetric_ZT(std::vector<float> receivedMessage) {
	std::vector<std::vector<cell>> trellisInfo;

	trellisInfo = constructLowRateTrellis_ZT(receivedMessage, METRIC_TYPE::PRODUCT_METRIC);

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

		std::vector<int> message = pathToMessage_ZT(path);

		// another way to compute the angle
		// currentAngleExplored = utils::compute_angle_between_vectors_rad(receivedMessage, codeword);
		currentAngleExplored = std::acos( std::max(-1.0f, std::min(1.0f, -forwardPartialPathMetric/N)) );
		
		// one trellis decoding requires both a tb and crc check
		if(path[0] == path[lowrate_pathLength - 1] && path[0] == 0 && crc::crc_check(message, crcDegree, crc) && currentAngleExplored <= MAX_ANGLE){
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

MessageInformation LowRateListDecoder::lowRateDecoding_MaxMetric_EuclideanMetric_ZT(std::vector<float> receivedMessage) {
	std::vector<std::vector<cell>> trellisInfo;

	trellisInfo = constructLowRateTrellis_ZT(receivedMessage, METRIC_TYPE::EUCLIDEAN_METRIC);

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
	float currentMetricExplored = 0.0;
	
	while(currentMetricExplored < MAX_METRIC){
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
		std::vector<int> message = pathToMessage_ZT(path);
		currentMetricExplored = forwardPartialPathMetric;
		// std::cout << "current Metric Explored: " << currentMetricExplored << std::endl;
		
		// one trellis decoding requires both a tb and crc check
		if(path[0] == path[lowrate_pathLength - 1] && path[0] == 0 && crc::crc_check(message, crcDegree, crc) && currentMetricExplored <= MAX_METRIC){
			output.message = message;
			output.path = path;
			output.listSize = numPathsSearched + 1;
			output.metric = forwardPartialPathMetric;
			output.TBListSize = TBPathsSearched + 1;
			// std::cout << "returning metric = " << forwardPartialPathMetric << std::endl;
			// std::cout << "--" << std::endl;
			return output;
		}

		numPathsSearched++;
		if(path[0] == path[lowrate_pathLength - 1])
			TBPathsSearched++;
	} // while(currentMetricExplored < MAX_METRIC)

	output.listSizeExceeded = true;
	output.listSize = numPathsSearched;
	std::cerr << "[WARNING]: TC IS NOT FOUND!!! " << std::endl;
	return output;
}


// converts a path through the ztcc trellis to the binary message it corresponds with
std::vector<int> LowRateListDecoder::pathToMessage_ZT(std::vector<int> path){
	std::vector<int> message;
	for(int pathIndex = 0; pathIndex < path.size() - 1 - V; pathIndex++){
		for(int forwardPath = 0; forwardPath < numForwardPaths; forwardPath++){
			if(lowrate_nextStates[path[pathIndex]][forwardPath] == path[pathIndex + 1])
				message.push_back(forwardPath);
		}
	}
	return message;
}