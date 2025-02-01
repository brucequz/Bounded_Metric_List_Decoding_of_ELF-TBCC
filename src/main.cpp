#include <iostream>
#include <fstream>
#include <vector>

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

  ISTC_sim(code);

  return 0;
}

void ISTC_sim(CodeInformation code){
	
	/* - Output files setup - */
	std::ofstream RRVtoTransmitted_MetricFile;
  RRVtoTransmitted_MetricFile.open("output/RRV/transmitted_metric.txt");
  if (!RRVtoTransmitted_MetricFile.is_open()) {
      std::cerr << "Error: Could not open the file output/RRV/transmitted_metric.txt" << std::endl;
      return;
  }

	std::ofstream RRVtoTransmitted_ListSizeFile;
  RRVtoTransmitted_ListSizeFile.open("output/RRV/transmitted_listsize.txt");
  if (!RRVtoTransmitted_ListSizeFile.is_open()) {
      std::cerr << "Error: Could not open the file output/RRV/transmitted_listsize.txt" << std::endl;
      return;
  }

	std::ofstream RRVtoDecoded_MetricFile;
  RRVtoDecoded_MetricFile.open("output/RRV/decoded_metric.txt");
  if (!RRVtoDecoded_MetricFile.is_open()) {
      std::cerr << "Error: Could not open the file output/RRV/decoded_metric.txt" << std::endl;
      return;
  }

	std::ofstream RRVtoDecoded_ListSizeFile;
  RRVtoDecoded_ListSizeFile.open("output/RRV/decoded_listsize.txt");
  if (!RRVtoDecoded_ListSizeFile.is_open()) {
      std::cerr << "Error: Could not open the file output/RRV/decoded_listsize.txt" << std::endl;
      return;
  }

	std::ofstream RRVtoDecoded_DecodeTypeFile;
  RRVtoDecoded_DecodeTypeFile.open("output/RRV/decoded_type.txt");
  if (!RRVtoDecoded_DecodeTypeFile.is_open()) {
      std::cerr << "Error: Could not open the file output/RRV/decoded_type.txt" << std::endl;
      return;
  }

	// PRV
	std::ofstream PRVtoTransmitted_MetricFile;
  PRVtoTransmitted_MetricFile.open("output/PRV/transmitted_metric.txt");
  if (!PRVtoTransmitted_MetricFile.is_open()) {
      std::cerr << "Error: Could not open the file output/PRV/transmitted_metric.txt" << std::endl;
      return;
  }

	std::ofstream PRVtoTransmitted_ListSizeFile;
  PRVtoTransmitted_ListSizeFile.open("output/PRV/transmitted_listsize.txt");
  if (!PRVtoTransmitted_ListSizeFile.is_open()) {
      std::cerr << "Error: Could not open the file output/PRV/transmitted_listsize.txt" << std::endl;
      return;
  }

	std::ofstream PRVtoDecoded_MetricFile;
  PRVtoDecoded_MetricFile.open("output/PRV/decoded_metric.txt");
  if (!PRVtoDecoded_MetricFile.is_open()) {
      std::cerr << "Error: Could not open the file output/PRV/decoded_metric.txt" << std::endl;
      return;
  }

	std::ofstream PRVtoDecoded_ListSizeFile;
  PRVtoDecoded_ListSizeFile.open("output/PRV/decoded_listsize.txt");
  if (!PRVtoDecoded_ListSizeFile.is_open()) {
      std::cerr << "Error: Could not open the file output/PRV/decoded_listsize.txt" << std::endl;
      return;
  }

	std::ofstream PRVtoDecoded_DecodeTypeFile;
  PRVtoDecoded_DecodeTypeFile.open("output/PRV/decoded_type.txt");
  if (!PRVtoDecoded_DecodeTypeFile.is_open()) {
      std::cerr << "Error: Could not open the file output/PRV/decoded_type.txt" << std::endl;
      return;
  }


	// History
	std::ofstream PRVtoTransmittedCodeword_PathHistoryFile;
  PRVtoTransmittedCodeword_PathHistoryFile.open("output/history/PRVtoTC.txt");
  if (!PRVtoTransmittedCodeword_PathHistoryFile.is_open()) {
      std::cerr << "Error: Could not open the file output/history/PRVtoTC.txt" << std::endl;
      return;
  }

	std::ofstream RRVtoPRV_ScalingFactorFile;
  RRVtoPRV_ScalingFactorFile.open("output/history/RRVtoPRVscaling.txt");
  if (!RRVtoPRV_ScalingFactorFile.is_open()) {
      std::cerr << "Error: Could not open the file output/history/RRVtoPRVscaling.txt" << std::endl;
      return;
  }

	std::ofstream DecodedCodewordNoiseMagFile;
	DecodedCodewordNoiseMagFile.open("output/history/decodedCodewordNoiseMag.txt");
	if (!DecodedCodewordNoiseMagFile.is_open()) {
			std::cerr << "Error: Could not open the file output/history/decodedCodewordNoiseMag.txt" << std::endl;
			return;
	}

	std::cout << "running the ISTC sim" << std::endl;
	srand(42);
	if ((code.numInfoBits + code.crcDeg - 1) % code.k != 0) {
		std::cout << "invalid msg + crc length" << std::endl;
		return;
	}

	std::vector<double> EbN0 = EBN0;
  std::vector<int> puncturedIndices = PUNCTURING_INDICES;

	std::vector<double> SNR;
	double offset = 10 * log10((double)2*64 / (128));
	for (int i=0; i< EbN0.size(); i++)
		SNR.push_back(EbN0[i] + offset);

	FeedForwardTrellis encodingTrellis(code.k, code.n, code.v, code.numerators);

	// Raw Received Value (RRV) List Decoder
	LowRateListDecoder listDecoder(encodingTrellis, MAX_LISTSIZE, code.crcDeg, code.crc);
	std::vector<double> RRVtoTransmitted_Metric(MAX_ITERATIONS);
	std::vector<int> 		RRVtoTransmitted_ListSize(MAX_ITERATIONS);
	std::vector<double> RRVtoDecoded_Metric(MAX_ITERATIONS);
	std::vector<int> 		RRVtoDecoded_ListSize(MAX_ITERATIONS);
	std::vector<int>		RRV_DecodedType(MAX_ITERATIONS);


	// Projected Received Value (PRV) List Decoder
	LowRateListDecoder listDecoder_Normalized(encodingTrellis, MAX_LISTSIZE, code.crcDeg, code.crc);
	std::vector<double> PRVtoTransmitted_Metric(MAX_ITERATIONS);
	std::vector<int> 		PRVtoTransmitted_ListSize(MAX_ITERATIONS);
	std::vector<double> PRVtoDecoded_Metric(MAX_ITERATIONS);
	std::vector<int> 		PRVtoDecoded_ListSize(MAX_ITERATIONS);
	std::vector<int>		PRV_DecodedType(MAX_ITERATIONS);

	// Projected Received Value to Transmitted Codeword History
	std::vector<std::vector<double>> PRVtoTC_History(MAX_ITERATIONS, std::vector<double>());

	// RRV to PRV Scaling Factor History
	std::vector<double> RRVtoPRV_ScalingFactor(MAX_ITERATIONS);

	// Decoding Noise Magnitude
	std::vector<std::vector<double>> DecodedCodewordNoiseMag(MAX_ITERATIONS, std::vector<double>());

	/* ==== SIMULATION begins ==== */
	for(int snrIndex = 0; snrIndex < SNR.size(); snrIndex++) {
		double snr = SNR[snrIndex];
		double standardMeanListSize = 0;
		int standardNumErrors = 0;
		int standardListSizeExceeded = 0;



		for(int numTrials = 0; numTrials < MAX_ITERATIONS; numTrials++) {

			if (numTrials % 1000 == 0) { std::cout << "currently at " << numTrials << std::endl; }
			
			std::vector<int> originalMessage = generateRandomCRCMessage(code);
			std::vector<int> transmittedMessage = generateTransmittedMessage(originalMessage, encodingTrellis, snr, puncturedIndices, NOISELESS);
			std::vector<double> receivedMessage = addAWNGNoise(transmittedMessage, puncturedIndices, snr, NOISELESS);
			std::vector<int> zero_point(receivedMessage.size(), 0);
			
			
			int sum_abs = std::accumulate(transmittedMessage.begin(), transmittedMessage.end(), 0, [](int sum, int x) {
        return sum + std::abs(x);
			});
			
			// Normalize the received signals to 128 and decode
			double origial_mag_sqrt = utils::euclidean_distance(receivedMessage, zero_point, puncturedIndices);
			
			std::vector<double> normalized_receivedMessage(receivedMessage.size(), 0);
			for (int i = 0; i < receivedMessage.size(); i++) {
				if (std::find(puncturedIndices.begin(), puncturedIndices.end(), i) == puncturedIndices.end()) {
					normalized_receivedMessage[i] = sqrt(128) * receivedMessage[i] / origial_mag_sqrt;
				}
			}

			
			
			// RRV to PRV Scaling Factor
			RRVtoPRV_ScalingFactor[numTrials] = sqrt(128) / origial_mag_sqrt;
			

			// TRANSMITTED
			RRVtoTransmitted_Metric[numTrials] = (utils::sum_of_squares(receivedMessage, transmittedMessage, puncturedIndices));
			PRVtoTransmitted_Metric[numTrials] = (utils::sum_of_squares(normalized_receivedMessage, transmittedMessage, puncturedIndices));
			
			// DECODED
      MessageInformation standardDecoding = listDecoder.lowRateDecoding(receivedMessage, puncturedIndices);
			MessageInformation normalizedDecoding = listDecoder_Normalized.lowRateDecoding_mla(normalized_receivedMessage, puncturedIndices, transmittedMessage);


			// RRV
			if (standardDecoding.message == originalMessage) {
				// correct decoding
				RRV_DecodedType[numTrials] 				= 0;
				RRVtoDecoded_ListSize[numTrials] 	= standardDecoding.listSize;
				RRVtoDecoded_Metric[numTrials] 		= standardDecoding.metric;
			} else if(standardDecoding.listSizeExceeded) {
				// list size exceeded
				RRV_DecodedType[numTrials] = 1;
				standardListSizeExceeded++;
			} else { 
				// incorrect decoding
				RRV_DecodedType[numTrials] 				= 2;
				RRVtoDecoded_ListSize[numTrials] 	= standardDecoding.listSize;
				RRVtoDecoded_Metric[numTrials] 		= standardDecoding.metric;
				standardNumErrors++;
				standardMeanListSize += (double)standardDecoding.listSize;
				std::cout << "decoding error:  ----------"  << std::endl; 
				std::cout << "standardDecoding.listSize: " << standardDecoding.listSize << std::endl;
				std::cout << "standardDecoding.metric: " << standardDecoding.metric << std::endl;
			}

			// PRV
			if (normalizedDecoding.message == originalMessage) {
				// correct decoding
				PRV_DecodedType[numTrials] 				= 0;
				PRVtoDecoded_ListSize[numTrials] 	= normalizedDecoding.listSize;
				PRVtoDecoded_Metric[numTrials] 		= normalizedDecoding.metric;
				if (normalizedDecoding.metric >= 71 && normalizedDecoding.listSize < 1600) {
					std::cout << "Large metric but small list size correct decoding: ----------" << std::endl;
					std::cout << "metric: " << normalizedDecoding.metric << std::endl;
					std::cout << "list size: " << normalizedDecoding.listSize << std::endl;
					std::cout << "TB list size: " << normalizedDecoding.TBListSize << std::endl;
				}
			} else if (normalizedDecoding.listSizeExceeded) {
				// list size exceeded
				PRV_DecodedType[numTrials] = 1;
			} else {
				// incorrect decoding
				PRV_DecodedType[numTrials] 				= 2;
				PRVtoDecoded_ListSize[numTrials] 	= normalizedDecoding.listSize;
				PRVtoDecoded_Metric[numTrials]		= normalizedDecoding.metric;
			}

			// PRV to TC History
			PRVtoTC_History[numTrials] = normalizedDecoding.pathToTransmittedCodewordHistory;

			// Decoded Codeword Noise Magnitude
			DecodedCodewordNoiseMag[numTrials] = normalizedDecoding.decodedCodewordSquaredNoiseMag;
			
		} // for(int numTrials = 0; numTrials < MAX_ITERATIONS; numTrials++)


		std::cout << "at Eb/N0 = " << EbN0[snrIndex] << std::endl;
		std::cout << "mean list size: " << (double)standardMeanListSize/(MAX_ITERATIONS - standardListSizeExceeded) << std::endl;
		std::cout << "number of times list size exceeded: " << standardListSizeExceeded << std::endl;
		std::cout << "number of errors: " << standardNumErrors << std::endl;
		std::cout << "UER: " << (double)standardNumErrors/(MAX_ITERATIONS - standardListSizeExceeded) << std::endl;
		std::cout << "TFR: " << (double)(standardNumErrors + standardListSizeExceeded)/MAX_ITERATIONS << std::endl;
		
		double scaling_average = std::accumulate(RRVtoPRV_ScalingFactor.begin(), RRVtoPRV_ScalingFactor.end(), 0.0) / RRVtoPRV_ScalingFactor.size();
		std::cout << "scaling average: " << scaling_average << std::endl;

		// RRV Write to file
		for (int i = 0; i < RRVtoTransmitted_Metric.size(); i++) {
			RRVtoTransmitted_MetricFile << RRVtoTransmitted_Metric[i] << std::endl;
		}
		for (int i = 0; i < RRVtoDecoded_Metric.size(); i++) {
			RRVtoDecoded_MetricFile << RRVtoDecoded_Metric[i] << std::endl;
		}
		for (int i = 0; i < RRVtoDecoded_ListSize.size(); i++) {
			RRVtoDecoded_ListSizeFile << RRVtoDecoded_ListSize[i] << std::endl;
		}
		for (int i = 0; i < RRV_DecodedType.size(); i++) {
			RRVtoDecoded_DecodeTypeFile << RRV_DecodedType[i] << std::endl;
		}

		// PRV Write to file
		for (int i = 0; i < PRVtoTransmitted_Metric.size(); i++) {
			PRVtoTransmitted_MetricFile << PRVtoTransmitted_Metric[i] << std::endl;
		}
		for (int i = 0; i < PRVtoDecoded_Metric.size(); i++) {
			PRVtoDecoded_MetricFile << PRVtoDecoded_Metric[i] << std::endl;
		}
		for (int i = 0; i < PRVtoDecoded_ListSize.size(); i++) {
			PRVtoDecoded_ListSizeFile << PRVtoDecoded_ListSize[i] << std::endl;
		}
		for (int i = 0; i < PRV_DecodedType.size(); i++) {
			PRVtoDecoded_DecodeTypeFile << PRV_DecodedType[i] << std::endl;
		}

		// PRV to TC History Write to file
		for (int i_trial = 0; i_trial < MAX_ITERATIONS; i_trial++) {
			for (int i = 0; i < PRVtoTC_History[i_trial].size(); i++) {
				PRVtoTransmittedCodeword_PathHistoryFile << PRVtoTC_History[i_trial][i] << ", ";
			}
			PRVtoTransmittedCodeword_PathHistoryFile << std::endl;
		}

		// RRV to PRV Scaling Factor Write to file
		for (int i = 0; i < RRVtoPRV_ScalingFactor.size(); i++) {
			RRVtoPRV_ScalingFactorFile << RRVtoPRV_ScalingFactor[i] << std::endl;
		}

		// Decoded Codeword Noise Magnitude Write to file
		for (int i = 0; i < DecodedCodewordNoiseMag.size(); i++) {
			if (DecodedCodewordNoiseMag[i].size() == 0) {
				for (int j = 0; j < 2*NUM_INFO_BITS; j++) {
					DecodedCodewordNoiseMagFile << 0.0 << ", ";
				}
			} else {
				for (int j = 0; j < DecodedCodewordNoiseMag[i].size(); j++) {
					DecodedCodewordNoiseMagFile << DecodedCodewordNoiseMag[i][j] << ", ";
				}
			}
			DecodedCodewordNoiseMagFile << std::endl;
		}

	} // for(int snrIndex = 0; snrIndex < SNR.size(); snrIndex++)

	// RRV
  RRVtoTransmitted_MetricFile.close();
	RRVtoTransmitted_ListSizeFile.close();
	RRVtoDecoded_MetricFile.close();
	RRVtoDecoded_ListSizeFile.close();
	RRVtoDecoded_DecodeTypeFile.close();

	// PRV
	PRVtoTransmitted_MetricFile.close();
	PRVtoTransmitted_ListSizeFile.close();
	PRVtoDecoded_MetricFile.close();
	PRVtoDecoded_ListSizeFile.close();
	PRVtoDecoded_DecodeTypeFile.close();

	// Path History
	PRVtoTransmittedCodeword_PathHistoryFile.close();
	RRVtoPRV_ScalingFactorFile.close();
	DecodedCodewordNoiseMagFile.close();

	std::cout << "concluded simulation" << std::endl;
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