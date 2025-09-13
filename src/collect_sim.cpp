#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <cstdlib>
#include <cassert>
#include <numeric>
#include <string>
#include <sstream>
#include "mpi.h"

#include "../consts.h"
#include "../include/types.h"
#include "../include/namespace.h"
#include "../include/feedForwardTrellis.h"
#include "../include/lowRateListDecoder.h"

void ISTC_sim(CodeInformation code, int rank);
std::vector<int> generateTransmittedMessage(std::vector<int> info_crc, FeedForwardTrellis encodingTrellis);
void logSimulationParams();

int main(int argc, char *argv[]) {
    
  CodeInformation code;
  code.k = k;         // numerator of the rate
  code.n = n;         // denominator of the rate
  code.v = V;         // number of memory elements
  code.crcDeg = M+1;  // m+1, degree of CRC, # bits of CRC polynomial
  code.crc = CRC;     // CRC polynomial
  code.numInfoBits = K; // number of information bits
  code.numerators = NUMERATORS;

	
	/* MPI Init */
	MPI_Init(&argc, &argv);
	int world_rank, world_size;
	MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
	MPI_Comm_size(MPI_COMM_WORLD, &world_size);

	awgn::generator.seed(BASE_SEED + world_rank);  // Seed with rank

	if (world_rank == 0) {
		logSimulationParams();
	}

	MPI_Barrier(MPI_COMM_WORLD);

	/* Check */
	if ((code.numInfoBits + code.crcDeg - 1) % code.k != 0) {
			std::cerr << "invalid msg + crc length" << std::endl;
			exit(1);
	}

	srand(BASE_SEED + world_rank);  // Reproducible random seed

	
	if (STOPPING_RULE == 'M' && DECODING_RULE == 'N') {
		if (MAX_METRIC_VEC.size() >= 1) {
			for (int i = MAX_METRIC_VEC.size()-1; i >= 0; i--) {
				MAX_METRIC = MAX_METRIC_VEC[i];
				std::cout << "MAX_METRIC = " << MAX_METRIC << std::endl;
				ISTC_sim(code, world_rank);  // Run simulation
			}
		}
	} else if ((STOPPING_RULE == 'A' || STOPPING_RULE == 'L') && DECODING_RULE == 'P') {
		if (MAX_ANGLE_VEC.size() >= 1) {
			for (int i = MAX_ANGLE_VEC.size()-1; i >= 0; i--) {
				MAX_ANGLE = MAX_ANGLE_VEC[i];
				std::cout << "MAX_ANGLE = " << MAX_ANGLE << std::endl;
				std::cout << "world rank: " << world_rank << std::endl;
				ISTC_sim(code, world_rank);  // Run simulation
			}
		}
	}

	MPI_Finalize();

  return 0;
}


void ISTC_sim(CodeInformation code, int rank){
 	for (size_t ebn0_id = 0; ebn0_id < EBN0.size(); ebn0_id++) {

		/* - Output files setup - */
		float EbN0 = EBN0[ebn0_id];
		std::ostringstream ebn0_str;
		ebn0_str.precision(2);
		ebn0_str << std::fixed << EbN0;

		std::ostringstream thetad_str;
		thetad_str.precision(4);
		thetad_str << std::fixed << MAX_ANGLE;

		std::ostringstream nonProjDist_str;
		nonProjDist_str.precision(4);
		nonProjDist_str << std::fixed << MAX_METRIC;

		std::ostringstream ude_error_cnt_str;
		ude_error_cnt_str.precision(1);
		ude_error_cnt_str << std::fixed << MAX_ERRORS;
		
		std::string folder_name;
		if (STOPPING_RULE == 'A' && DECODING_RULE == 'P') {
			// for projected & angle decoding
			folder_name = "output/BALD/Curve_Sim_thetad_" + thetad_str.str() + "/EbN0_" + ebn0_str.str() + "/Proc" + std::to_string(rank);
		} else if (DECODING_RULE == 'N' && STOPPING_RULE == 'M') {
			// for non-projected & metric decoding
			folder_name = "output/BDLD/Curve_Sim_dist_" + nonProjDist_str.str() + "/EbN0_" + ebn0_str.str() + "/Proc" + std::to_string(rank);
		} else {
			folder_name = "output/Proc" + std::to_string(rank) + "_EbN0_" + ebn0_str.str() + "_ude_" + ude_error_cnt_str.str();
		}
		system(("mkdir -p " + folder_name).c_str());
		
		std::string RtoT_Metric_filename = folder_name + "/transmitted_metric.txt";
		std::ofstream RRVtoTransmitted_MetricFile(RtoT_Metric_filename.c_str());

		std::string RtoD_Metric_filename = folder_name + "/decoded_metric.txt";
		std::ofstream RRVtoDecoded_MetricFile(RtoD_Metric_filename.c_str());

		std::string RtoD_LS_filename = folder_name + "/decoded_listsize.txt";
		std::ofstream RRVtoDecoded_ListSizeFile(RtoD_LS_filename.c_str());

		std::string RtoD_Angle_filename = folder_name + "/decoded_angle.txt";
		std::ofstream RRVtoDecoded_AngleFile(RtoD_Angle_filename);
		
		std::string RtoD_Type_filename = folder_name + "/decoded_type.txt";
		std::ofstream RRVtoDecoded_DecodeTypeFile(RtoD_Type_filename);

    std::string Running_minimum_highrate_to_tx_filename = folder_name + "/running_minimum_highrate_to_tx.txt";
    std::ofstream Running_minimum_highrate_to_tx_File(Running_minimum_highrate_to_tx_filename);

		
		/* - Simulation esno_dB setup - */
		std::vector<int> puncturedIndices = PUNCTURING_INDICES;
		float esno_dB = 0.0;
		float offset = 10 * log10((float)K/N);
		esno_dB = EbN0 + offset;
		
		/* - Trellis setup - */
		FeedForwardTrellis encodingTrellis(code.k, code.n, code.v, code.numerators);

		/* - Decoder setup - */
		LowRateListDecoder listDecoder(encodingTrellis, MAX_LISTSIZE, code.crcDeg, code.crc, STOPPING_RULE);

		/* - Output Temporary Holder setup - */
		std::vector<float> RRVtoTransmitted_Metric;
		std::vector<float> RRVtoDecoded_Metric;
		std::vector<int> 		RRVtoDecoded_ListSize;
		std::vector<float> RRVtoDecoded_Angle;
		std::vector<int>		RRV_DecodedType;

    std::vector<std::vector<float>> twoD_vec_running_minimum;

    std::vector<float> sampling_points;
    float sample_start = 5.3f;
    float sample_end   = 8.9f;
    int   num_samples  = 100;
    float sample_step  = (sample_end - sample_start) / (num_samples - 1);
    for (int id_sample = 0; id_sample < num_samples; id_sample++) {
      sampling_points.push_back(sample_start + id_sample * sample_step);
    }
    std::cout << "printing sampling points: ";
    utils::print_float_vector(sampling_points);


		/* ==== SIMULATION begins ==== */
		std::cout << std::endl << "**- Simulation Started for EbN0 = " << std::fixed << std::setprecision(2) << EbN0 << " -**" << std::endl;
		int num_mistakes 	= 0;
		int num_failures 	= 0;
		int num_errors 	 	= 0; // num_mistakes + num_failures
		int num_trials	 	= 0;

		// lambda to decide if we continue the loop or not
		auto should_continue = [&]() -> bool {
			if (ERROR_RUN_TYPE == 'U') {
					return num_mistakes < MAX_ERRORS;
			} else if (ERROR_RUN_TYPE == 'T') {
					return num_errors < MAX_ERRORS;
			} else {
					throw std::runtime_error("Unknown TER_TYPE");
			}
		};

		auto should_end_of_file_log = [&]() -> bool {
			if (ERROR_RUN_TYPE == 'U') {
				return num_mistakes == MAX_ERRORS;
			} else if (ERROR_RUN_TYPE == 'T') {
				// std::cout << "num_errors = " << num_errors << std::endl;
				return num_errors == MAX_ERRORS;
			} else {
				throw std::runtime_error("Unknown TER_TYPE");
			}
		};

		while (should_continue()) {
			

			std::vector<int> originalMessage = crc::generateRandomCRCMessage(code);
			// std::cout << "original message: ";
			// utils::print_int_vector(originalMessage);
			// std::cout << std::endl;
			std::vector<int> transmittedMessage = generateTransmittedMessage(originalMessage, encodingTrellis);
			// std::cout << "transmitted message: ";
			// utils::print_int_vector(transmittedMessage);
			// std::cout << ", length = " << transmittedMessage.size() << std::endl;
			std::vector<float> receivedMessage = awgn::addAWNGNoise(transmittedMessage, puncturedIndices, esno_dB, NOISELESS);
			// std::cout << "received message size: " << receivedMessage.size() << std::endl;
			// std::cout << "received message: ";
			// utils::print_float_vector(receivedMessage);
			
			// Produce Raw Channel values
			float esno = pow(10.0, esno_dB / 10.0);
			for (size_t i = 0; i < receivedMessage.size(); i++) {
				receivedMessage[i] = receivedMessage[i] / (4 * esno);
			}

			// Transmitted statistics
			RRVtoTransmitted_Metric.push_back(utils::sum_of_squares(receivedMessage, transmittedMessage, puncturedIndices));
			
			// Project Received Message onto the codeword sphere
			MessageInformation decodingResult;
			float sigma_sqrd = pow(10.0, -esno_dB / 10.0) / 2.0;
			if (DECODING_RULE == 'P') {
				float received_word_energy = utils::compute_vector_energy(receivedMessage);
				float energy_normalize_factor = std::sqrt(N / received_word_energy);  // normalizing received message
				std::vector<float> projected_received_word(receivedMessage.size(), 0.0);
				for (size_t i = 0; i < receivedMessage.size(); i++) {
					projected_received_word[i] = receivedMessage[i] * energy_normalize_factor;
				}
				// Decoding
				// decodingResult = listDecoder.genieAided_LowRateDecoding_MaxListsize(projected_received_word, puncturedIndices, transmittedMessage, sampling_points);
        decodingResult = listDecoder.genieAided_LowRateDecoding_MaxAngle_ProductMetric(projected_received_word, transmittedMessage, puncturedIndices, sampling_points);
			} else if (DECODING_RULE == 'N') {
				decodingResult = listDecoder.decode(receivedMessage, puncturedIndices, sigma_sqrd);
			}
			

			// RRV
			if (!decodingResult.listSizeExceeded && decodingResult.message == originalMessage) {
				// correct decoding
				RRV_DecodedType.push_back(0);
				RRVtoDecoded_ListSize.push_back(decodingResult.listSize);
				RRVtoDecoded_Metric.push_back(decodingResult.metric);
				RRVtoDecoded_Angle.push_back(decodingResult.angle_received_decoded_rad);
        twoD_vec_running_minimum.push_back(decodingResult.vec_running_minimum);
				std::cout << "Sim id " << num_trials << ": Correct decoding! " << std::endl;
			} else if(decodingResult.listSizeExceeded) {
				// list size exceeded
				RRV_DecodedType.push_back(1);
				RRVtoDecoded_ListSize.push_back(decodingResult.listSize);
				num_failures++;
				std::cout << "Sim id " << num_trials << ": List size exceeded! num_mistakes = " << num_mistakes << std::endl;
			} else { 
				// incorrect decoding
				RRV_DecodedType.push_back(2);
				RRVtoDecoded_ListSize.push_back(decodingResult.listSize);
				RRVtoDecoded_Metric.push_back(decodingResult.metric);
				RRVtoDecoded_Angle.push_back(decodingResult.angle_received_decoded_rad);
				// twoD_vec_running_minimum.push_back(decodingResult.vec_running_minimum);
				num_mistakes++;
				std::cout << "Sim id " << num_trials << ": Undetected error! num_mistakes = " << num_mistakes << std::endl;
			}

			// Increment errors and trials
			num_errors = num_mistakes + num_failures;
			num_trials += 1;

			if (num_trials % LOGGING_ITERS == 0 || should_end_of_file_log()) {
				if (ERROR_RUN_TYPE == 'U') {std::cout << "numTrials = " << num_trials << ", number of undetected errors = " << num_mistakes << ", number of total errors = " << num_errors << std::endl;}
				if (ERROR_RUN_TYPE == 'T') {std::cout << "numTrials = " << num_trials << ", number of total errors = " << num_errors << std::endl;}
				 
				// RRV Write to file
				if (RRVtoDecoded_AngleFile.is_open()) {
					for (size_t i = 0; i < RRVtoDecoded_Angle.size(); i++) {
						RRVtoDecoded_AngleFile << std::setprecision(5) << RRVtoDecoded_Angle[i] << std::endl;
					}
					RRVtoDecoded_Angle.clear();
				}
				if (RRVtoTransmitted_MetricFile.is_open()) {
					for (size_t i = 0; i < RRVtoTransmitted_Metric.size(); i++) {
						RRVtoTransmitted_MetricFile << RRVtoTransmitted_Metric[i] << std::endl;
					}
					RRVtoTransmitted_Metric.clear();
				}
				if (RRVtoDecoded_MetricFile.is_open()) {
					for (size_t i = 0; i < RRVtoDecoded_Metric.size(); i++) {
						RRVtoDecoded_MetricFile << RRVtoDecoded_Metric[i] << std::endl;
					}
					RRVtoDecoded_Metric.clear();
				}
				if (RRVtoDecoded_ListSizeFile.is_open()) {
					for (size_t i = 0; i < RRVtoDecoded_ListSize.size(); i++) {
						RRVtoDecoded_ListSizeFile << RRVtoDecoded_ListSize[i] << std::endl;
					}
					RRVtoDecoded_ListSize.clear();
				}
				if (RRVtoDecoded_DecodeTypeFile.is_open()) {
					for (size_t i = 0; i < RRV_DecodedType.size(); i++) {
						RRVtoDecoded_DecodeTypeFile << RRV_DecodedType[i] << std::endl;
					}
					RRV_DecodedType.clear();
				}
        if (Running_minimum_highrate_to_tx_File.is_open()) {
					if (!twoD_vec_running_minimum.empty()) {
						for (size_t i = 0; i < twoD_vec_running_minimum.size(); i++) {
							if (twoD_vec_running_minimum[i].empty()) {
								continue;
							} else {
								for (size_t j = 0; j < twoD_vec_running_minimum[0].size()-1; j++) {
									Running_minimum_highrate_to_tx_File << twoD_vec_running_minimum[i][j] << ", ";
								}
								Running_minimum_highrate_to_tx_File << twoD_vec_running_minimum[i].back() << std::endl;
							}
						}
						// clear the 2d vectors
						twoD_vec_running_minimum.clear();
					}
        }

				
			} // if (num_trials % LOGGING_ITERS == 0 || num_errors == MAX_ERRORS)
			if (num_trials == 100) {num_errors = MAX_ERRORS;}
		} // while (num_mistakes < MAX_ERRORS)

		std::cout << std::endl << "At Eb/N0 = " << std::fixed << std::setprecision(2) << EbN0 << std::endl;
		std::cout << "number of total errors: " << num_errors << std::endl;
		std::cout << "number of undetected errors: " << num_mistakes << std::endl;
		std::cout << "number of detected errors: " << num_failures << std::endl;
		std::cout << "Undetected Error Rate: " << std::scientific << (float)num_mistakes/num_trials << std::endl;
		std::cout << "Detected Error Rate: " << std::scientific << (float)num_failures/num_trials << std::endl;
		std::cout << "TFR: " << (float)num_errors/num_trials << std::endl;
		std::cout << "*- Simulation Concluded for EbN0 = " << std::fixed << std::setprecision(2) << EbN0 << " -*" << std::endl;

		

		// Close Files
		RRVtoTransmitted_MetricFile.close();
		RRVtoDecoded_MetricFile.close();
		RRVtoDecoded_ListSizeFile.close();
		RRVtoDecoded_DecodeTypeFile.close();
		RRVtoDecoded_AngleFile.close();
    Running_minimum_highrate_to_tx_File.close();
	} // for (size_t ebn0_id = 0; ebn0_id < EBN0.size(); ebn0_id++) 

	std::cout << "***--- Simulation Concluded ---***" << std::endl;
}




// this takes the message bits, including the CRC, and encodes them using the trellis
std::vector<int> generateTransmittedMessage(std::vector<int> info_crc, FeedForwardTrellis encodingTrellis){
	/*
	encodes to get the transmitted message bits (info + zero termination + crc) before modulation.
	*/ 
	if (ENCODING_RULE != 'T' && ENCODING_RULE != 'Z') {std::cerr << "ISSUE: INCORRECT ENCODING_RULE" << std::endl;}
	// encode the message
	std::vector<int> encodedMessage;
	if (ENCODING_RULE == 'T') {
		encodedMessage = encodingTrellis.encode(info_crc);
		assert(encodedMessage.size() == (K+M) / k * n);
	} else if (ENCODING_RULE == 'Z') {
		for (int i=0; i<V; i++){
			info_crc.push_back(0);
		}
		// std::cout << "info crc with termination: ";
		// utils::print_int_vector(info_crc);
		// std::cout << std::endl;
		encodedMessage = encodingTrellis.encode_zt(info_crc);
		assert(encodedMessage.size() == (K+M+V) / k * n);
		
	}
	return encodedMessage;
}

void logSimulationParams() {
	std::cout << "+----------------------+------------+\n";
	std::cout << "| Parameter           | Value      |\n";
	std::cout << "+----------------------+------------+\n";

	/// ---------------- CODE INFO ----------------
	std::cout << "| " << std::left << std::setw(20) << "K"
						<< "| " << std::setw(10) << K << "|\n";
	std::cout << "| " << std::left << std::setw(20) << "N"
						<< "| " << std::setw(10) << N << "|\n";
	// std::cout << "| " << std::left << std::setw(20) << "GEN POLY"
	// 					<< "| " << "{" << POLY1 << ", " << POLY2 << "}" << "|\n";
	std::cout << "| " << std::left << std::setw(20) << "ELF"
						<< "| " << "0x" << std::setw(10) << std::hex << CRC << std::dec << "|\n";
	std::cout << "| " << std::left << std::setw(20) << "STOPPING RULE"
						<< "| " << std::setw(10) << STOPPING_RULE << "|\n";
	/// ---------------- ERROR TYPE ----------------
	if (ERROR_RUN_TYPE == 'T' || ERROR_RUN_TYPE == 'U') {
		std::cout << "| " << std::left << std::setw(20) << "ERROR TYPE"
						<< "| " << std::setw(10) << ERROR_RUN_TYPE << "|\n";
		std::cout << "| " << std::left << std::setw(20) << "ACCUMULATE"
						<< "| " << std::setw(10) << MAX_ERRORS << "|\n";
	} else {std::cerr << "INCORRECT ERROR TYPE! ABORT!"; exit(1);}
	/// ---------------- STOPPING RULE ----------------
	if (STOPPING_RULE == 'M') {
		std::cout << "| " << std::left << std::setw(20) << "MAX METRIC"
						<< "| " << std::setw(10) << MAX_METRIC << "|\n";
	} else if (STOPPING_RULE == 'L') {
		std::cout << "| " << std::left << std::setw(20) << "MAX LISTSIZE"
						<< "| " << std::setw(10) << MAX_LISTSIZE << "|\n";
	} else if (STOPPING_RULE == 'A') {
		std::cout << "| " << std::left << std::setw(20) << "MAX ANGLE"
						<< "| " << std::setw(10) << MAX_ANGLE << "|\n";
	} else if (STOPPING_RULE == 'R') {
		std::cout << "| " << std::left << std::setw(20) << "ROVA THRESHOLD"
						<< "| " << std::setw(10) << 4 << "|\n";		
	} else {std::cerr << "INCORRECT STOPPING RULE! ABORT!"; exit(1);}
	if (DECODING_RULE == 'P' || DECODING_RULE == 'N') {
		std::cout << "| " << std::left << std::setw(20) << "DECODING RULE"
						<< "| " << std::setw(10) << DECODING_RULE << "|\n";
	} else {std::cerr << "INCORRECT DECODING RULE! ABORT!"; exit(1);}
	/// ---------------- SIMULATION PARAMS ----------------
	std::cout << "| " << std::left << std::setw(20) << "MAX ERRORS"
						<< "| " << std::setw(10) << MAX_ERRORS << "|\n";
	std::cout << "| " << std::left << std::setw(20) << "NOISELESS?"
						<< "| " << std::setw(10) << NOISELESS << "|\n";
	std::cout << "| " << std::left << std::setw(20) << "LOGGING ITERS"
						<< "| " << std::setw(10) << LOGGING_ITERS << "|\n";
	std::cout << "| " << std::left << std::setw(20) << "BASE SEED"
	<< "| " << std::setw(10) << BASE_SEED << "|\n";

	std::cout << "+----------------------+------------+\n";
}