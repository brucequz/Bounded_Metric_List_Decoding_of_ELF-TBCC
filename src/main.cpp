#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <cstdlib>
#include <string>
#include <sstream>

#include "../include/mla_consts.h"
#include "../include/mla_types.h"
#include "../include/mla_namespace.h"
#include "../include/feedForwardTrellis.h"
#include "../include/lowRateListDecoder.h"

void ISTC_sim(CodeInformation code);
std::vector<int> generateRandomCRCMessage(CodeInformation code);
std::vector<int> generateTransmittedMessage(std::vector<int> originalMessage, FeedForwardTrellis encodingTrellis, double snr, std::vector<int> puncturedIndices, bool noiseless);
std::vector<double> addAWNGNoise(std::vector<int> transmittedMessage, std::vector<int> puncturedIndices, double snr, bool noiseless);

int main() {

  std::cout << "K: " << K << std::endl;
  std::cout << "N: " << N << std::endl;
  std::cout << "V: " << V << std::endl;
  std::cout << "M: " << M << std::endl;
    
  CodeInformation code;
  code.k = K;         // numerator of the rate
  code.n = N;         // denominator of the rate
  code.v = V;         // number of memory elements
  code.crcDeg = M+1;  // m+1, degree of CRC, # bits of CRC polynomial
  code.crc = CRC;     // CRC polynomial
  code.numInfoBits = NUM_INFO_BITS; // number of information bits
  code.numerators = {POLY1, POLY2};

	srand(42);

	/* Check */
	assert( (N/K) * (NUM_INFO_BITS + M) - PUNCTURING_INDICES.size() == NUM_CODED_SYMBOLS);
	if ((code.numInfoBits + code.crcDeg - 1) % code.k != 0) {
		std::cerr << "invalid msg + crc length" << std::endl;
		exit(1);
	}
	
		
	ISTC_sim(code);
		
  return 0;
}

void ISTC_sim(CodeInformation code){

	for (size_t ebn0_id = 0; ebn0_id < EBN0.size(); ebn0_id++) {

		/* - Output files setup - */
		double EbN0 = EBN0[ebn0_id];
		std::ostringstream ebn0_str;
		ebn0_str.precision(2);
		ebn0_str << std::fixed << EbN0;
		
		std::string folder_name = "output/EbN0_" + ebn0_str.str();
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
		std::vector<int> 		RRVtoTransmitted_ListSize;
		std::vector<double> RRVtoDecoded_Metric;
		std::vector<int> 		RRVtoDecoded_ListSize;
		std::vector<int>		RRV_DecodedType;

		/* ==== SIMULATION begins ==== */
		std::cout << std::endl << "**- Simulation Started for EbN0 = " << std::fixed << std::setprecision(1) << EbN0 << " -**" << std::endl;
		int num_mistakes 	= 0;
		int num_failures 	= 0;
		int num_errors 	 	= 0; // num_mistakes + num_failures
		int num_trials	 	= 0;

		while (num_errors < MAX_ERRORS) {

			std::vector<int> originalMessage = generateRandomCRCMessage(code);
			std::vector<int> transmittedMessage = generateTransmittedMessage(originalMessage, encodingTrellis, snr, puncturedIndices, NOISELESS);
			std::vector<double> receivedMessage = addAWNGNoise(transmittedMessage, puncturedIndices, snr, NOISELESS);
			std::vector<int> zero_point(receivedMessage.size(), 0);
		

			// Transmitted statistics
			RRVtoTransmitted_Metric.push_back(utils::sum_of_squares(receivedMessage, transmittedMessage, puncturedIndices));
			
			// Decoding
			MessageInformation standardDecoding = listDecoder.decode(receivedMessage, puncturedIndices);
			

			// RRV
			if (standardDecoding.message == originalMessage) {
				// correct decoding
				RRV_DecodedType.push_back(0);
				RRVtoDecoded_ListSize.push_back(standardDecoding.listSize);
				RRVtoDecoded_Metric.push_back(standardDecoding.metric);
			} else if(standardDecoding.listSizeExceeded) {
				// list size exceeded
				RRV_DecodedType.push_back(1);
				num_failures++;
			} else { 
				// incorrect decoding
				RRV_DecodedType.push_back(2);
				RRVtoDecoded_ListSize.push_back(standardDecoding.listSize);
				RRVtoDecoded_Metric.push_back(standardDecoding.metric);
				num_mistakes++;
			}

			// Increment errors and trials
			num_errors = num_failures + num_mistakes;
			num_trials += 1;

			if (num_trials % 1000 == 0) { std::cout << "numTrials = " << num_trials << ", numErrors = " << num_errors << std::endl; }
		} // while (num_errors <= MAX_ERRORS) 

		std::cout << std::endl << "At Eb/N0 = " << std::fixed << std::setprecision(1) << EbN0 << std::endl;
		std::cout << "number of errors: " << num_errors << std::endl;
		std::cout << "number of mistakes: " << num_mistakes << std::endl;
		std::cout << "number of failures: " << num_failures << std::endl;
		std::cout << "Mistakes Error Rate: " << std::scientific << (double)num_mistakes/num_trials << std::endl;
		std::cout << "Failures Error Rate: " << std::scientific << (double)num_failures/num_trials << std::endl;
		std::cout << "TFR: " << (double)num_errors/num_trials << std::endl;
		std::cout << "*- Simulation Concluded for EbN0 = " << std::fixed << std::setprecision(1) << EbN0 << " -*" << std::endl;

		// RRV Write to file
		if (RRVtoTransmitted_MetricFile.is_open()) {
			for (int i = 0; i < RRVtoTransmitted_Metric.size(); i++) {
				RRVtoTransmitted_MetricFile << RRVtoTransmitted_Metric[i] << std::endl;
			}
		}
		if (RRVtoDecoded_MetricFile.is_open()) {
			for (int i = 0; i < RRVtoDecoded_Metric.size(); i++) {
				RRVtoDecoded_MetricFile << RRVtoDecoded_Metric[i] << std::endl;
			}
		}
		if (RRVtoDecoded_ListSizeFile.is_open()) {
			for (int i = 0; i < RRVtoDecoded_ListSize.size(); i++) {
				RRVtoDecoded_ListSizeFile << RRVtoDecoded_ListSize[i] << std::endl;
			}
		}
		if (RRVtoDecoded_DecodeTypeFile.is_open()) {
			for (int i = 0; i < RRV_DecodedType.size(); i++) {
				RRVtoDecoded_DecodeTypeFile << RRV_DecodedType[i] << std::endl;
			}
		}

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