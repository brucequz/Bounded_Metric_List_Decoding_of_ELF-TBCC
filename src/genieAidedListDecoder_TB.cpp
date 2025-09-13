#include "../include/lowRateListDecoder.h"


std::vector<float> LowRateListDecoder::push_to_angle_boundary(std::vector<float> receivedMessage, std::vector<int> transmittedCodeword, std::vector<int> punctured_indices, float angle) {
  // push the received word (RX) to the verge of the decoding sphere around the transmitted word (TX)
  // This requires computing the vector indicating the direction of the push
  // push the received word to the location where angle is equal to the expected angle
  // and then normalize the vector so that it lands on the sphere with radius sqrt(N)

  size_t rec_length = receivedMessage.size();
  std::vector<float> direction_to_push(rec_length, 0.0);
	// std::cout << "received length: " << rec_length << std::endl;
  for (size_t i = 0; i < rec_length; i++) {
		if (std::find(punctured_indices.begin(), punctured_indices.end(), i) != punctured_indices.end()) {
			continue;
		} else {
			direction_to_push[i] = receivedMessage[i] - transmittedCodeword[i];
		}
  }
	// std::cout << "direction to push: ";
	// utils::print_float_vector(direction_to_push);
	// std::cout << std::endl;

  direction_to_push = utils::normalize_to_unit_energy(direction_to_push);

  float direction_to_push_energy = utils::compute_vector_energy(direction_to_push);
  // std::cout << "printing the direction to push energy: " << direction_to_push_energy << std::endl;
	// std::cout << "direction to push size: " << direction_to_push.size() << std::endl;

  float original_angle_rad = utils::compute_angle_between_vectors_rad(receivedMessage, transmittedCodeword, punctured_indices);
  float received_word_energy = utils::compute_vector_energy(receivedMessage);
  // std::cout << "printing the original angle (rad) between received message and transmitted: " << std::setprecision(4) << original_angle_rad << std::endl;
  // std::cout << "printing the original energy of received message: " << received_word_energy << std::endl;

  // two_alpha is the angle between RX and TX
  float two_alpha_rad = utils::compute_angle_between_vectors_rad(receivedMessage, transmittedCodeword, punctured_indices);
  float alpha = two_alpha_rad / 2.0f;
  float beta = angle - two_alpha_rad;
  // std::cout << "printing angle: " << std::setprecision(4) << angle << ", printing two_alphas: " << two_alpha_rad << std::endl;
  // std::cout << "printing alpha: " << alpha <<  ", printing beta: " << beta << std::endl;

  
  float ell =  0.5 * utils::euclidean_distance(receivedMessage, transmittedCodeword, punctured_indices);
  // std::cout << "printing ell: " << ell << std::endl;

  float r_prime = ell / std::tan(alpha);
  // std::cout << "printing r_prime: " << r_prime << std::endl;

  float tangent_ratio = std::tan(alpha) / std::tan(alpha+beta);
  // std::cout << "printing x: " << (ell / tangent_ratio - ell) << std::endl;
  // std::cout << "printing tangent ratio: " << tangent_ratio << std::endl;
  std::vector<float> push_vector(rec_length, 0.0);
  std::vector<float> received_on_verge(rec_length, 0.0);

  for (size_t i = 0; i < rec_length; i++) {
    push_vector[i] = direction_to_push[i] * (ell / tangent_ratio - ell);
    received_on_verge[i] = receivedMessage[i] + push_vector[i];
  }

  float length_push_vector = std::sqrt(utils::compute_vector_energy(push_vector));
  // std::cout << "printing the length of the push vector: " << length_push_vector << std::endl; 

  float test = utils::compute_angle_between_vectors_rad(received_on_verge, transmittedCodeword, punctured_indices);
  // std::cout << "printing angle (rad) between received on verge and transmitted: " << std::setprecision(4) << test << std::endl;

  std::vector<float> rx_verge_to_tx(rec_length, 0.0);
  std::vector<float> rx_to_tx(rec_length, 0.0);
  for (size_t i = 0; i < rx_verge_to_tx.size(); i++) {
    rx_verge_to_tx[i] = float(received_on_verge[i] - transmittedCodeword[i]);
    rx_to_tx[i] = float(receivedMessage[i] - transmittedCodeword[i]);
  }

  float rx_verge_energy = utils::compute_vector_energy(received_on_verge);
  float energy_normalize_factor = std::sqrt(N / rx_verge_energy);  // normalizing received message
  std::vector<float> projected_rx_verge(received_on_verge.size(), 0.0);
  for (size_t i = 0; i < received_on_verge.size(); i++) {
    projected_rx_verge[i] = received_on_verge[i] * energy_normalize_factor;
  }

  float energy_check = utils::compute_vector_energy(projected_rx_verge);
  // std::cout << "energy check: " << energy_check << std::endl;
  float test2 = utils::compute_angle_between_vectors_rad(projected_rx_verge, transmittedCodeword, punctured_indices);
  // std::cout << "printing angle (rad) between projected received on verge and transmitted: " << std::setprecision(4) << test2 << std::endl;

  return projected_rx_verge;

}

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


MessageInformation LowRateListDecoder::genieAided_LowRateDecoding_MaxAngle_ProductMetric(std::vector<float> receivedMessage, std::vector<int> transmittedCodeword, std::vector<int> punctured_indices, std::vector<float> sampling_points){
	
  // push to angle boundary
  receivedMessage = push_to_angle_boundary(receivedMessage, transmittedCodeword, punctured_indices, MAX_ANGLE);

  std::vector<std::vector<cell>> trellisInfo;
	trellisInfo = constructLowRateTrellis_Punctured_ProductMetric(receivedMessage, punctured_indices);

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
  float running_min_highrate_to_tx = 1000.0;
  size_t last_end = 0;
	
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


    /// We compute and record the following,
    // 1. Euclidean distance from transmitted word to high-rate codeword.
    // 2. Euclidean distance from projected received word to high-rate codeword.
    float highrate_to_tx = utils::euclidean_distance(codeword, transmittedCodeword, punctured_indices);
    float highrate_to_projected_rx = utils::euclidean_distance(codeword, receivedMessage, punctured_indices);
		// output.vec_highrate_to_tx.push_back(highrate_to_tx);
		// output.vec_highrate_to_projected_rx.push_back(highrate_to_projected_rx);
    
    running_min_highrate_to_tx = std::min(running_min_highrate_to_tx, highrate_to_tx);

    for (size_t id = last_end; id < sampling_points.size(); id++) {
      if (highrate_to_projected_rx >= sampling_points[id]) {
        output.vec_running_minimum[id] = running_min_highrate_to_tx;
      } else {
        last_end = id;
        break;
      }
    }

		// another way to compute the angle
		currentAngleExplored = std::acos( std::max(-1.0f, std::min(1.0f, -forwardPartialPathMetric/N)) );
		
		// one trellis decoding requires both a tb and crc check
		if(path[0] == path[lowrate_pathLength - 1] && crc::crc_check(message, crcDegree, crc) && currentAngleExplored <= MAX_ANGLE){
			output.message = message;
			output.path = path;
			output.listSize = numPathsSearched + 1;
			output.metric = forwardPartialPathMetric;
			output.TBListSize = TBPathsSearched + 1;
			output.angle_received_decoded_rad = currentAngleExplored;
			
			// std::string folder_name = "output/temp";
			// std::string highrate_to_tx_and_projected_rx_filename = folder_name + "/highrate_to_tx_and_projected_rx" + std::to_string(output.listSize) + ".txt";
			// std::ofstream highrate_to_tx_and_projected_rx_File(highrate_to_tx_and_projected_rx_filename);

			// if (highrate_to_tx_and_projected_rx_File.is_open()) {
			// 	for (size_t i = 0; i < output.vec_highrate_to_tx.size(); i++) {
			// 		highrate_to_tx_and_projected_rx_File << std::setprecision(3) << output.vec_highrate_to_tx[i] << ", " << output.vec_highrate_to_projected_rx[i] << std::endl;
			// 	}
			// } else {
			// 	std::cerr << "file not opened." << std::endl;
			// }

			// highrate_to_tx_and_projected_rx_File.close();
			
			return output;
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