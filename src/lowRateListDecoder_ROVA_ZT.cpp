#include "lowRateListDecoder.h"



// std::vector<std::vector<LowRateListDecoder::rova_cell>> LowRateListDecoder::constructLowRateTrellis_ROVA_Alg2_ZT(std::vector<float> receivedMessage){
// 	std::vector<std::vector<rova_cell>> trellisInfo;
// 	lowrate_pathLength = (receivedMessage.size() / lowrate_symbolLength) + 1;

// 	trellisInfo = std::vector<std::vector<rova_cell>>(lowrate_numStates, std::vector<rova_cell>(lowrate_pathLength));
//   int trellis_width = trellisInfo[0].size();
//   int trellis_height = trellisInfo.size();
//   std::vector<std::vector<std::vector<float>>> gammas(
//     trellis_width-1, std::vector<std::vector<float>>(
//       trellis_height, std::vector<float>(
//         numForwardPaths, 1.0
//       )
//     )
//   );

// 	// initialize only 0 as the starting states
// 	trellisInfo[0][0].pathMetric = 0;
// 	trellisInfo[0][0].init = true;
//   trellisInfo[0][0].Gamma = 1.0;
	
// 	// building the trellis
// 	for(int stage = 0; stage < lowrate_pathLength - V - 1; stage++){
// 		for(int currentState = 0; currentState < lowrate_numStates; currentState++){
// 			// if the state / stage is invalid, we move on
// 			if(!trellisInfo[currentState][stage].init)
// 				continue;

// 			// otherwise, we compute the relevent information
// 			for(int forwardPathIndex = 0; forwardPathIndex < numForwardPaths; forwardPathIndex++){
// 				// since our transitions correspond to symbols, the forwardPathIndex has no correlation 
// 				// beyond indexing the forward path

// 				int nextState = lowrate_nextStates[currentState][forwardPathIndex];
				
// 				// if the nextState is invalid, we move on
// 				if(nextState < 0)
// 					continue;
				
// 				float branchMetric = 0;
// 				std::vector<int> output_point = crc::get_point(lowrate_outputs[currentState][forwardPathIndex], lowrate_symbolLength);
				
// 				for(int i = 0; i < lowrate_symbolLength; i++){
// 					branchMetric += std::pow(receivedMessage[lowrate_symbolLength * stage + i] - (float)output_point[i], 2);
//           gammas[stage][currentState][forwardPathIndex] *= awgn::normpdf(receivedMessage[lowrate_symbolLength * stage + i], (float)output_point[i], 0.3749);
// 				}
        
// 				float totalPathMetric = branchMetric + trellisInfo[currentState][stage].pathMetric;
//         float input_Gamma = trellisInfo[currentState][stage].Gamma * gammas[stage][currentState][forwardPathIndex];
				
// 				// dealing with cases of uninitialized states, when the transition becomes the optimal father state, and suboptimal father state, in order
// 				if(!trellisInfo[nextState][stage + 1].init){
// 					trellisInfo[nextState][stage + 1].pathMetric = totalPathMetric;
// 					trellisInfo[nextState][stage + 1].optimalFatherState = currentState;
// 					trellisInfo[nextState][stage + 1].init = true;
// 				}
// 				else if(trellisInfo[nextState][stage + 1].pathMetric > totalPathMetric){
// 					trellisInfo[nextState][stage + 1].suboptimalPathMetric = trellisInfo[nextState][stage + 1].pathMetric;
// 					trellisInfo[nextState][stage + 1].suboptimalFatherState = trellisInfo[nextState][stage + 1].optimalFatherState;
// 					trellisInfo[nextState][stage + 1].pathMetric = totalPathMetric;
// 					trellisInfo[nextState][stage + 1].optimalFatherState = currentState;
// 				}
// 				else{
// 					trellisInfo[nextState][stage + 1].suboptimalPathMetric = totalPathMetric;
// 					trellisInfo[nextState][stage + 1].suboptimalFatherState = currentState;
// 				}

//         if (trellisInfo[nextState][stage + 1].Gamma != 0) {
//           trellisInfo[nextState][stage + 1].P = std::max(trellisInfo[nextState][stage + 1].Gamma, input_Gamma) / (trellisInfo[nextState][stage + 1].Gamma + input_Gamma);
//         }
//         trellisInfo[nextState][stage + 1].Gamma = std::max(trellisInfo[nextState][stage + 1].Gamma, input_Gamma); // update Gamma

// 			} // for(int forwardPathIndex = 0; forwardPathIndex < numForwardPaths; forwardPathIndex++)

// 		} // for(int currentState = 0; currentState < lowrate_numStates; currentState++)
// 	} // for(int stage = 0; stage < lowrate_pathLength - V - 1; stage++)

// 	// ZT stage
// 	for(int stage = lowrate_pathLength - V - 1; stage < lowrate_pathLength - 1; stage++){
// 		for(int currentState = 0; currentState < lowrate_numStates; currentState++){
// 			// if the state / stage is invalid, we move on
// 			if(!trellisInfo[currentState][stage].init)
// 				continue;

// 			// zero terminating
// 			int forwardPathIndex = 0;
// 			// since our transitions correspond to symbols, the forwardPathIndex has no correlation 
// 			// beyond indexing the forward path

// 			int nextState = lowrate_nextStates[currentState][forwardPathIndex];
			
// 			// if the nextState is invalid, we move on
// 			if(nextState < 0)
// 				continue;
			
// 			float branchMetric = 0;
// 			std::vector<int> output_point = crc::get_point(lowrate_outputs[currentState][forwardPathIndex], lowrate_symbolLength);
			
// 			for(int i = 0; i < lowrate_symbolLength; i++){
// 				branchMetric += std::pow(receivedMessage[lowrate_symbolLength * stage + i] - (float)output_point[i], 2);
//         gammas[stage][currentState][forwardPathIndex] *= awgn::normpdf(receivedMessage[lowrate_symbolLength * stage + i], (float)output_point[i], 0.3749);
// 			}
// 			float totalPathMetric = branchMetric + trellisInfo[currentState][stage].pathMetric;
//       float input_Gamma = trellisInfo[currentState][stage].Gamma * gammas[stage][currentState][forwardPathIndex];
			
// 			// dealing with cases of uninitialized states, when the transition becomes the optimal father state, and suboptimal father state, in order
// 			if(!trellisInfo[nextState][stage + 1].init){
// 				trellisInfo[nextState][stage + 1].pathMetric = totalPathMetric;
// 				trellisInfo[nextState][stage + 1].optimalFatherState = currentState;
// 				trellisInfo[nextState][stage + 1].init = true;
// 			}
// 			else if(trellisInfo[nextState][stage + 1].pathMetric > totalPathMetric){
// 				trellisInfo[nextState][stage + 1].suboptimalPathMetric = trellisInfo[nextState][stage + 1].pathMetric;
// 				trellisInfo[nextState][stage + 1].suboptimalFatherState = trellisInfo[nextState][stage + 1].optimalFatherState;
// 				trellisInfo[nextState][stage + 1].pathMetric = totalPathMetric;
// 				trellisInfo[nextState][stage + 1].optimalFatherState = currentState;
// 			}
// 			else{
// 				trellisInfo[nextState][stage + 1].suboptimalPathMetric = totalPathMetric;
// 				trellisInfo[nextState][stage + 1].suboptimalFatherState = currentState;
// 			}


//       trellisInfo[nextState][stage + 1].P = std::max(trellisInfo[nextState][stage + 1].Gamma, input_Gamma) / (trellisInfo[nextState][stage + 1].Gamma + input_Gamma);
      
//       trellisInfo[nextState][stage + 1].Gamma = std::max(trellisInfo[nextState][stage + 1].Gamma, input_Gamma); // update Gamma
//       // std::cout << "Zero terminating Gamma: nextState: " << nextState << ", stage+1: " << stage+1 << ", P: "
//       // << std::fixed << std::setprecision(50) << trellisInfo[nextState][stage + 1].P << std::endl;

// 		} // for(int currentState = 0; currentState < lowrate_numStates; currentState++
// 	} // for(int stage = lowrate_pathLength - V - 1; stage < lowrate_pathLength - 1; stage++)


//   // std::cout << "before returning trellis: " << std::endl;
//   // std::cout << trellisInfo[0][lowrate_pathLength - 1].P;


// 	return trellisInfo;
// }

std::vector<std::vector<LowRateListDecoder::rova_cell>> LowRateListDecoder::constructLowRateTrellis_ROVA_Alg4_ZT(std::vector<float> receivedMessage, float sigma_sqrd){
	std::vector<std::vector<rova_cell>> trellisInfo;
	lowrate_pathLength = (receivedMessage.size() / lowrate_symbolLength) + 1;

	// std::cout << "Sigma sqrd: " << std::setprecision(10) << sigma_sqrd << std::endl;

	trellisInfo = std::vector<std::vector<rova_cell>>(lowrate_numStates, std::vector<rova_cell>(lowrate_pathLength));
  int trellis_width = trellisInfo[0].size();
  int trellis_height = trellisInfo.size();
  
	log_gammas_ = std::vector<std::vector<std::vector<float>>>(
    trellis_width - 1,
    std::vector<std::vector<float>>(
        trellis_height,
        std::vector<float>(numForwardPaths, 0.0)
    )
	);

	// initialize only 0 as the starting states
	trellisInfo[0][0].pathMetric = 0;
	trellisInfo[0][0].init = true;
	trellisInfo[0][0].log_Gamma = 0; // log_Gamma is initialized to 0
	trellisInfo[0][0].log_Z = 0; // log_Z is initialized to 0
	
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
					branchMetric += std::pow(receivedMessage[lowrate_symbolLength * stage + i] - (float)output_point[i], 2);
					log_gammas_[stage][currentState][forwardPathIndex] += awgn::log_normpdf(receivedMessage[lowrate_symbolLength * stage + i], (float)output_point[i], sqrt(sigma_sqrd));
				}
        
				float totalPathMetric = branchMetric + trellisInfo[currentState][stage].pathMetric;
        float input_Gamma = trellisInfo[currentState][stage].log_Gamma + log_gammas_[stage][currentState][forwardPathIndex];
				
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

        if (trellisInfo[nextState][stage + 1].log_Gamma != 0) {
          trellisInfo[nextState][stage + 1].log_Z = max_star(trellisInfo[currentState][stage].log_Z + log_gammas_[stage][currentState][forwardPathIndex], trellisInfo[nextState][stage + 1].log_Z);
        } else {
					trellisInfo[nextState][stage + 1].log_Z = trellisInfo[currentState][stage].log_Z + log_gammas_[stage][currentState][forwardPathIndex];
				}
				// update Gamma
        trellisInfo[nextState][stage + 1].log_Gamma = std::max(trellisInfo[nextState][stage + 1].log_Gamma, input_Gamma);

			} // for(int forwardPathIndex = 0; forwardPathIndex < numForwardPaths; forwardPathIndex++)

		} // for(int currentState = 0; currentState < lowrate_numStates; currentState++)
	} // for(int stage = 0; stage < lowrate_pathLength - V - 1; stage++)

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
				branchMetric += std::pow(receivedMessage[lowrate_symbolLength * stage + i] - (float)output_point[i], 2);
        log_gammas_[stage][currentState][forwardPathIndex] += awgn::log_normpdf(receivedMessage[lowrate_symbolLength * stage + i], (float)output_point[i], sqrt(sigma_sqrd));
			}
			float totalPathMetric = branchMetric + trellisInfo[currentState][stage].pathMetric;
      float input_Gamma = trellisInfo[currentState][stage].log_Gamma + log_gammas_[stage][currentState][forwardPathIndex];
			
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

      if (trellisInfo[nextState][stage + 1].log_Gamma != 0) {
				trellisInfo[nextState][stage + 1].log_Z = max_star(trellisInfo[currentState][stage].log_Z + log_gammas_[stage][currentState][forwardPathIndex], trellisInfo[nextState][stage + 1].log_Z);
			} else {
				trellisInfo[nextState][stage + 1].log_Z = trellisInfo[currentState][stage].log_Z + log_gammas_[stage][currentState][forwardPathIndex];
			}
			// update Gamma
			trellisInfo[nextState][stage + 1].log_Gamma = std::max(trellisInfo[nextState][stage + 1].log_Gamma, input_Gamma);

		} // for(int currentState = 0; currentState < lowrate_numStates; currentState++
	} // for(int stage = lowrate_pathLength - V - 1; stage < lowrate_pathLength - 1; stage++)

  // std::cout << "before returning trellis: " << std::endl;
	// std::cout << "Gamma: " << trellisInfo[0][lowrate_pathLength - 1].log_Gamma << std::endl;
	// std::cout << "Z: " << trellisInfo[0][lowrate_pathLength - 1].log_Z << std::endl;
	// std::cout << "first row of trellisInfo: " << std::endl;
	// for (const auto& cell : trellisInfo[0]) {
	// 	std::cout << "log Gamma: " << std::fixed << std::setprecision(50) << cell.log_Gamma << ", Z: " << cell.log_Z << std::endl;
	// }
  // float p = std::exp(trellisInfo[0].back().log_Gamma - trellisInfo[0].back().log_Z);
	// std::cout << "initial P: " << std::fixed << std::setprecision(50) << p << std::endl;

	return trellisInfo;
}

float LowRateListDecoder::compute_logGamma(std::vector<float> receivedMessage, std::vector<int> codeword, float sigma_sqrd) {
	float logGamma = 0.0;
	for (size_t i = 0; i < codeword.size(); i++) {
		logGamma += awgn::log_normpdf(receivedMessage[i], (float)codeword[i], sqrt(sigma_sqrd));
	}
	return logGamma;
}


float LowRateListDecoder::max_star(float lnx, float lny) {
	/* computes max star approximation in BCJR decoding
	 * 
	 * max_star(x, y) = max(x, y) + log(1 + exp(-|x - y|))
	 * 
	 * This is a more numerically stable version of the max star approximation
	 * that avoids overflow and underflow issues.
	 * 
	 * Assumes x and y are log values. 
	 */
	if (lnx > lny) {
		return lnx + std::log(1.0 + std::exp(lny - lnx));
	} else {
		return lny + std::log(1.0 + std::exp(lnx - lny));
	}
}


MessageInformation LowRateListDecoder::lowRateDecoding_SquaredDistanceMetric_ROVA_ZT(std::vector<float> receivedMessage, float sigma_sqrd) {
	std::vector<std::vector<rova_cell>> trellisInfo;

	trellisInfo = constructLowRateTrellis_ROVA_Alg4_ZT(receivedMessage, sigma_sqrd);


	// start search
	MessageInformation output;
	std::vector<MessageInformation> out_queue;
	//RBTree detourTree;
	MinHeap detourTree;
	std::vector<std::vector<int>> previousPaths;
	

	// create nodes for zero terminating paths
	DetourObject detour;
	detour.startingState = 0;
	detour.pathMetric = trellisInfo[0][lowrate_pathLength - 1].pathMetric;
	detourTree.insert(detour);

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

			// if there is a detour we can add to the detourTree
			if(trellisInfo[currentState][stage].suboptimalFatherState != -1){
				DetourObject localDetour;
				localDetour.detourStage = stage;
				localDetour.originalPathIndex = numPathsSearched;
				localDetour.pathMetric = suboptimalPathMetric + forwardPartialPathMetric;
				localDetour.forwardPathMetric = forwardPartialPathMetric;
				localDetour.startingState = detour.startingState;
				detourTree.insert(localDetour);
				// std::cout << "Detour added: " << localDetour.pathMetric << ", stage: " << stage << ", originalPathIndex: " << localDetour.originalPathIndex << std::endl;
				// std::cout << "currentPathMetric: " << currPathMetric << ", suboptimalPathMetric: " << suboptimalPathMetric << std::endl;
				// std::cout << "forwardPartialPathMetric: " << forwardPartialPathMetric << std::endl;
				
			}
			currentState = trellisInfo[currentState][stage].optimalFatherState;
			float prevPathMetric = trellisInfo[currentState][stage - 1].pathMetric;
			forwardPartialPathMetric += currPathMetric - prevPathMetric;
			path[stage - 1] = currentState;
		} // for(int stage = newTracebackStage; stage > 0; stage--)
		
		previousPaths.push_back(path);
		std::vector<int> codeword = pathToCodeword(path);
		// std::cout << "Decoded codeword: ";
		// utils::print_int_vector(codeword);
		// std::cout << "codeword size: " << codeword.size() << std::endl;
		float logGamma = compute_logGamma(receivedMessage, codeword, sigma_sqrd);
		std::vector<int> message = pathToMessage_ZT(path);
		// float p = std::exp(logGamma - trellisInfo[0].back().log_Z);
		// std::cout << "Path Metric: " << std::setprecision(10) << forwardPartialPathMetric << ", P: " << std::fixed << std::setprecision(10) << p << ", logGamma: " << std::setprecision(10) << logGamma << std::endl;
		// if (p < 5e-5) {
		// 	std::cout << "P is too small, break!" << std::endl;
		// 	break;
		// }

		// another way to compute the angle
		// currentAngleExplored = utils::compute_angle_between_vectors_rad(receivedMessage, codeword);
		
		
		// one trellis decoding requires both a tb and crc check
		if(path[0] == path[lowrate_pathLength - 1] && path[0] == 0 && crc::crc_check(message, crcDegree, crc) && currentAngleExplored <= MAX_ANGLE){
 			output.message = message;
			output.path = path;
			output.listSize = numPathsSearched + 1;
			output.metric = forwardPartialPathMetric;
			output.TBListSize = TBPathsSearched + 1;
			output.angle_received_decoded_rad = currentAngleExplored;
			output.log_Gamma = logGamma;
			out_queue.push_back(output);
			if (out_queue.size() == 2) {
				float first_logGamma = out_queue[0].log_Gamma;
				float second_logGamma = out_queue[1].log_Gamma;
				float rova_prob = std::exp(-std::log(1+std::exp(-first_logGamma + second_logGamma)));
				out_queue[0].rova_probability = rova_prob;
				

				if (rova_prob < ROVA_THRESHOLD) {
					std::cout << "returning number 1 with rova prob = " << std::setprecision(10) << rova_prob 
					<< " is smaller than ROVA_THRESHOLD = " << std::setprecision(2) << ROVA_THRESHOLD << std::endl;
					break;
				}
				return out_queue[0];
			}
		}

		numPathsSearched++;
		if(path[0] == path[lowrate_pathLength - 1])
			TBPathsSearched++;
	} // while(currentAngleExplored < MAX_ANGLE)

	output.listSizeExceeded = true;
	output.listSize = numPathsSearched;
	std::cerr << "[WARNING]: TC IS NOT FOUND!!! " << std::endl;
	return output;
}