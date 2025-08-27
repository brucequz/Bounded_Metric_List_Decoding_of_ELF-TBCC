#include "../include/lowRateListDecoder.h"


MessageInformation LowRateListDecoder::genieAided_LowRateDecoding_MaxListsize(std::vector<float> receivedMessage, std::vector<int> punctured_indices, std::vector<int> transmittedCodeword, std::vector<float> sampling_points){
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
  float running_min_highrate_to_tx = 1000.0;
  size_t last_end = 0;
  
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

    /// We compute and record the following,
    // 1. Euclidean distance from transmitted word to high-rate codeword.
    // 2. Euclidean distance from projected received word to high-rate codeword.

    float highrate_to_tx = utils::euclidean_distance(codeword, transmittedCodeword, punctured_indices);
    float highrate_to_projected_rx = utils::euclidean_distance(codeword, receivedMessage, punctured_indices);
    running_min_highrate_to_tx = std::min(running_min_highrate_to_tx, highrate_to_tx);

    
    for (size_t id = last_end; id < sampling_points.size(); id++) {
      if (highrate_to_projected_rx >= sampling_points[id]) {
        output.vec_running_minimum[id] = running_min_highrate_to_tx;
      } else {
        last_end = id;
        break;
      }
    }

    // std::cout << "printing just one vec_running_minimum: ";
    // utils::print_float_vector(output.vec_running_minimum);
    
    // std::cout << "listsize: " << numPathsSearched+1 << ", highrate_to_tx: " << std::setprecision(5) << highrate_to_tx << ", highrate_to_projected_rx: " << std::setprecision(5) << highrate_to_projected_rx <<
    // ", running min of the highrate_to_tx: " << std::setprecision(5) << running_min_highrate_to_tx << std::endl;
    // std::cout << std::endl;
		
		// one trellis decoding requires both a tb and crc check
		if(path[0] == path[lowrate_pathLength - 1] && crc::crc_check(message, crcDegree, crc) && numPathsSearched <= this->listSize){
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
	} // while(numPathsSearched < this->listSize)

	output.listSizeExceeded = true;
	return output;
}