#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <cstdlib>
#include <numeric>
#include <string>
#include <sstream>
#include "/opt/homebrew/Cellar/open-mpi/5.0.7/include/mpi.h"
// #include "mpi.h"

#include "../include/mla_consts.h"
#include "../include/mla_types.h"
#include "../include/mla_namespace.h"
#include "../include/feedForwardTrellis.h"
#include "../include/lowRateListDecoder.h"

void ISTC_sim(CodeInformation code, int rank);
std::vector<int> generateRandomCRCMessage(CodeInformation code);
std::vector<int> generateTransmittedMessage(std::vector<int> originalMessage, FeedForwardTrellis encodingTrellis, double snr, std::vector<int> puncturedIndices, bool noiseless);
std::vector<double> addAWNGNoise(std::vector<int> transmittedMessage, std::vector<int> puncturedIndices, double snr, bool noiseless);
void logSimulationParams();

int main(int argc, char *argv[]) {
    
  CodeInformation code;
  code.k = K;         // numerator of the rate
  code.n = N;         // denominator of the rate
  code.v = V;         // number of memory elements
  code.crcDeg = M+1;  // m+1, degree of CRC, # bits of CRC polynomial
  code.crc = CRC;     // CRC polynomial
  code.numInfoBits = NUM_INFO_BITS; // number of information bits
  code.numerators = {POLY1, POLY2};

	
	/* MPI Init */
	MPI_Init(&argc, &argv);
	int world_rank, world_size;
	MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
	MPI_Comm_size(MPI_COMM_WORLD, &world_size);

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

	ISTC_sim(code, world_rank);  // Run simulation

	MPI_Finalize();

  return 0;
}

void ISTC_sim(CodeInformation code, int rank){

	for (size_t ebn0_id = 0; ebn0_id < EBN0.size(); ebn0_id++) {

		/* - Output files setup - */
		double EbN0 = EBN0[ebn0_id];
		std::ostringstream ebn0_str;
		ebn0_str.precision(2);
		ebn0_str << std::fixed << EbN0;

		std::ostringstream ude_error_cnt_str;
		ude_error_cnt_str.precision(1);
		ude_error_cnt_str << std::fixed << MAX_ERRORS;
		
		std::string folder_name = "output/Proc" + std::to_string(rank) + "_EbN0_" + ebn0_str.str() + "_ude_" + ude_error_cnt_str.str();
		system(("mkdir -p " + folder_name).c_str());
		
		std::string RtoT_Metric_filename = folder_name + "/transmitted_metric.txt";
		std::ofstream RRVtoTransmitted_MetricFile(RtoT_Metric_filename.c_str());

		std::string RtoD_Metric_filename = folder_name + "/decoded_metric.txt";
		std::ofstream RRVtoDecoded_MetricFile(RtoD_Metric_filename.c_str());

		std::string RtoD_LS_filename = folder_name + "/decoded_listsize.txt";
		std::ofstream RRVtoDecoded_ListSizeFile(RtoD_LS_filename.c_str());
		
		std::string RtoD_Type_filename = folder_name + "/decoded_type.txt";
		std::ofstream RRVtoDecoded_DecodeTypeFile(RtoD_Type_filename);
		
		/* - Simulation SNR setup - */
		std::vector<int> puncturedIndices = PUNCTURING_INDICES;
		double snr = 0.0;
		double offset = 10 * log10((double)N/K *NUM_INFO_BITS / (NUM_CODED_SYMBOLS));
		snr = EbN0 + offset;
		
		/* - Trellis setup - */
		FeedForwardTrellis encodingTrellis(code.k, code.n, code.v, code.numerators);

		/* - Decoder setup - */
		LowRateListDecoder listDecoder(encodingTrellis, MAX_LISTSIZE, code.crcDeg, code.crc, STOPPING_RULE);

		/* - Output Temporary Holder setup - */
		std::vector<double> RRVtoTransmitted_Metric;
		std::vector<double> RRVtoDecoded_Metric;
		std::vector<int> 		RRVtoDecoded_ListSize;
		std::vector<int>		RRV_DecodedType;

		/* ==== SIMULATION begins ==== */
		std::cout << std::endl << "**- Simulation Started for EbN0 = " << std::fixed << std::setprecision(2) << EbN0 << " -**" << std::endl;
		int num_mistakes 	= 0;
		int num_failures 	= 0;
		int num_errors 	 	= 0; // num_mistakes + num_failures
		int num_trials	 	= 0;

		while (num_mistakes < MAX_ERRORS) {

			std::vector<int> originalMessage = generateRandomCRCMessage(code);
			std::vector<int> transmittedMessage = generateTransmittedMessage(originalMessage, encodingTrellis, snr, puncturedIndices, NOISELESS);
			std::vector<double> receivedMessage = addAWNGNoise(transmittedMessage, puncturedIndices, snr, NOISELESS);
		
			// Transmitted statistics
			RRVtoTransmitted_Metric.push_back(utils::sum_of_squares(receivedMessage, transmittedMessage, puncturedIndices));
			
			// Project Received Message onto the codeword sphere
			if (DECODING_RULE == 'P') {
				double received_word_energy = utils::compute_vector_energy(receivedMessage);
				double energy_normalize_factor = std::sqrt(128.0 / received_word_energy);  // normalizing received message
				std::vector<double> projected_received_word(receivedMessage.size(), 0.0);
				for (size_t i = 0; i < receivedMessage.size(); i++) {
					projected_received_word[i] = receivedMessage[i] * energy_normalize_factor;
				}
				// Decoding
				MessageInformation decodingResult = listDecoder.decode(projected_received_word, puncturedIndices);
			} else if (DECODING_RULE == 'N') {
				MessageInformation decodingResult = listDecoder.decode(receivedMessage, puncturedIndices);
			}
			

			// RRV
			if (decodingResult.message == originalMessage) {
				// correct decoding
				RRV_DecodedType.push_back(0);
				RRVtoDecoded_ListSize.push_back(decodingResult.listSize);
				RRVtoDecoded_Metric.push_back(decodingResult.metric);
			} else if(decodingResult.listSizeExceeded) {
				// list size exceeded
				RRV_DecodedType.push_back(1);
				num_failures++;
			} else { 
				// incorrect decoding
				RRV_DecodedType.push_back(2);
				RRVtoDecoded_ListSize.push_back(decodingResult.listSize);
				RRVtoDecoded_Metric.push_back(decodingResult.metric);
				num_mistakes++;
			}

			// Increment errors and trials
			num_errors = num_mistakes + num_failures;
			num_trials += 1;

			if (num_trials % LOGGING_ITERS == 0 || num_errors == MAX_ERRORS) {
				 std::cout << "numTrials = " << num_trials << ", numErrors = " << num_errors << std::endl; 

				// RRV Write to file
				if (RRVtoTransmitted_MetricFile.is_open()) {
					for (int i = 0; i < RRVtoTransmitted_Metric.size(); i++) {
						RRVtoTransmitted_MetricFile << RRVtoTransmitted_Metric[i] << std::endl;
					}
					RRVtoTransmitted_Metric.clear();
				}
				if (RRVtoDecoded_MetricFile.is_open()) {
					for (int i = 0; i < RRVtoDecoded_Metric.size(); i++) {
						RRVtoDecoded_MetricFile << RRVtoDecoded_Metric[i] << std::endl;
					}
					RRVtoDecoded_Metric.clear();
				}
				if (RRVtoDecoded_ListSizeFile.is_open()) {
					for (int i = 0; i < RRVtoDecoded_ListSize.size(); i++) {
						RRVtoDecoded_ListSizeFile << RRVtoDecoded_ListSize[i] << std::endl;
					}
					RRVtoDecoded_ListSize.clear();
				}
				if (RRVtoDecoded_DecodeTypeFile.is_open()) {
					for (int i = 0; i < RRV_DecodedType.size(); i++) {
						RRVtoDecoded_DecodeTypeFile << RRV_DecodedType[i] << std::endl;
					}
					RRV_DecodedType.clear();
				}
			} // if (num_trials % LOGGING_ITERS == 0 || num_errors == MAX_ERRORS)
		} // while (num_mistakes < MAX_ERRORS)

		std::cout << std::endl << "At Eb/N0 = " << std::fixed << std::setprecision(2) << EbN0 << std::endl;
		std::cout << "number of errors: " << num_errors << std::endl;
		std::cout << "number of mistakes: " << num_mistakes << std::endl;
		std::cout << "number of failures: " << num_failures << std::endl;
		std::cout << "Mistakes Error Rate: " << std::scientific << (double)num_mistakes/num_trials << std::endl;
		std::cout << "Failures Error Rate: " << std::scientific << (double)num_failures/num_trials << std::endl;
		std::cout << "TFR: " << (double)num_errors/num_trials << std::endl;
		std::cout << "*- Simulation Concluded for EbN0 = " << std::fixed << std::setprecision(2) << EbN0 << " -*" << std::endl;

		

		// RRV
		RRVtoTransmitted_MetricFile.close();
		RRVtoDecoded_MetricFile.close();
		RRVtoDecoded_ListSizeFile.close();
		RRVtoDecoded_DecodeTypeFile.close();
	} // for (size_t ebn0_id = 0; ebn0_id < EBN0.size(); ebn0_id++) 

	std::cout << "***--- Simulation Concluded ---***" << std::endl;
}


// this generates a random binary string of length code.numInfoBits, and appends the appropriate CRC bits
std::vector<int> generateRandomCRCMessage(CodeInformation code){
	std::vector<int> message;
	for(int i = 0; i < code.numInfoBits; i++)
		message.push_back(rand()%2);
	// compute the CRC
	crc::crc_calculation(message, code.crcDeg, code.crc);
	return message;
}

// this takes the message bits, including the CRC, and encodes them using the trellis
std::vector<int> generateTransmittedMessage(std::vector<int> originalMessage, FeedForwardTrellis encodingTrellis, double snr, std::vector<int> puncturedIndices, bool noiseless){
	// encode the message
	std::vector<int> encodedMessage = encodingTrellis.encode(originalMessage);

	return encodedMessage;
}

// this takes the transmitted message and adds AWGN noise to it
// it also punctures the bits that are not used in the trellis
std::vector<double> addAWNGNoise(std::vector<int> transmittedMessage, std::vector<int> puncturedIndices, double snr, bool noiseless){
	std::vector<double> receivedMessage;
	if(noiseless){
		for(int i = 0; i < transmittedMessage.size(); i++)
			receivedMessage.push_back((double)transmittedMessage[i]);
	} else {
		receivedMessage = awgn::addNoise(transmittedMessage, snr);
	}

	// puncture the bits. it is more convenient to puncture on this side than on the 
	// decoder, so we insert zeros which provide no information to the decoder
	for(int index = 0; index < puncturedIndices.size(); index++) {
		if (puncturedIndices[index] > receivedMessage.size()) {
			std::cout << "out of bounds index: " << puncturedIndices[index] << std::endl;
			std::cerr << "Puncturing index out of bounds" << std::endl;
			exit(1);
		}
		receivedMessage[puncturedIndices[index]] = 0;
	}

	return receivedMessage;
}

void logSimulationParams() {
	std::cout << "+----------------------+------------+\n";
	std::cout << "| Parameter           | Value      |\n";
	std::cout << "+----------------------+------------+\n";

	/// ---------------- CODE INFO ----------------
	std::cout << "| " << std::left << std::setw(20) << "K"
						<< "| " << std::setw(10) << NUM_INFO_BITS << "|\n";
	std::cout << "| " << std::left << std::setw(20) << "N"
						<< "| " << std::setw(10) << NUM_CODED_SYMBOLS << "|\n";
	std::cout << "| " << std::left << std::setw(20) << "GEN POLY"
						<< "| " << "{" << POLY1 << ", " << POLY2 << "}" << "|\n";
	std::cout << "| " << std::left << std::setw(20) << "ELF"
						<< "| " << "0x" << std::setw(10) << std::hex << CRC << std::dec << "|\n";
	std::cout << "| " << std::left << std::setw(20) << "STOPPING RULE"
						<< "| " << std::setw(10) << STOPPING_RULE << "|\n";
	/// ---------------- STOPPING RULE ----------------
	if (STOPPING_RULE == 'M') {
		std::cout << "| " << std::left << std::setw(20) << "MAX METRIC"
						<< "| " << std::setw(10) << MAX_METRIC << "|\n";
	} else if (STOPPING_RULE == 'L') {
		std::cout << "| " << std::left << std::setw(20) << "MAX LISTSIZE"
						<< "| " << std::setw(10) << MAX_LISTSIZE << "|\n";
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