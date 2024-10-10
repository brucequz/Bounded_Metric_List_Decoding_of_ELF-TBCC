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
std::vector<double> generateTransmittedMessage(std::vector<int> originalMessage, FeedForwardTrellis encodingTrellis, double snr, std::vector<int> puncturedIndices, bool noiseless);

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
	std::ofstream outFile;
  outFile.open("output/ISTC_sim_mla_listsize.txt");

  // Check if the file was opened successfully
  if (!outFile.is_open()) {
      std::cerr << "Error: Could not open the file ../output/ISTC_sim.txt" << std::endl;
      return;
  }
	std::cout << "running the ISTC sim" << std::endl;

	// set random seed for message generation
	srand(42);

	// check to make sure the code has valid values
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
	std::vector<double> correct_decoding_metrics;
  std::vector<int> listsize_mla;
	

	// the below are the relevant initializations for low rate
	// decoding, this will be altered for high rate or ZTCC, for example
	FeedForwardTrellis encodingTrellis(code.k, code.n, code.v, code.numerators);
	LowRateListDecoder listDecoder(encodingTrellis, MAX_LISTSIZE, code.crcDeg, code.crc);

	/* ---- SIMULATION begins ---- */
	for(int snrIndex = 0; snrIndex < SNR.size(); snrIndex++) {
		double snr = SNR[snrIndex];
		double standardMeanListSize = 0;
		int standardNumErrors = 0;
		int standardListSizeExceeded = 0;

		for(int numTrials = 0; numTrials < MAX_ITERATIONS; numTrials++) {

			if (numTrials % 100 == 0) { std::cout << "currently at " << numTrials << std::endl; }
			
			std::vector<int> originalMessage = generateRandomCRCMessage(code);
			std::vector<double> transmittedMessage = generateTransmittedMessage(originalMessage, encodingTrellis, snr, puncturedIndices, NOISELESS);

			// MessageInformation standardDecoding = listDecoder.lowRateDecoding(transmittedMessage, puncturedIndices);
      MessageInformation standardDecoding = listDecoder.lowRateDecoding_mla(transmittedMessage, puncturedIndices, TARGET_METRIC);

      // mla append the list size
      listsize_mla.push_back(standardDecoding.listSize);
      std::cout << "decoding metric: " << standardDecoding.metric << std::endl;
			
			if(standardDecoding.listSizeExceeded){
				standardListSizeExceeded++;
			}
			else if (standardDecoding.message != originalMessage){
				standardNumErrors ++;
				standardMeanListSize += (double)standardDecoding.listSize;
			} else {
				// save the metrics for the correct decodings
				correct_decoding_metrics.push_back(standardDecoding.metric);
			}
			
		} // for(int numTrials = 0; numTrials < MAX_ITERATIONS; numTrials++)


		std::cout << "at Eb/N0 = " << EbN0[snrIndex] << std::endl;
		std::cout << "mean list size: " << (double)standardMeanListSize/(MAX_ITERATIONS - standardListSizeExceeded) << std::endl;
		std::cout << "number of times list size exceeded: " << standardListSizeExceeded << std::endl;
		std::cout << "number of errors: " << standardNumErrors << std::endl;
		std::cout << "UER: " << (double)standardNumErrors/(MAX_ITERATIONS - standardListSizeExceeded) << std::endl;
		std::cout << "TFR: " << (double)(standardNumErrors + standardListSizeExceeded)/MAX_ITERATIONS << std::endl;


		for (int listsize : listsize_mla){
			outFile << listsize << std::endl;
		}

	} // for(int snrIndex = 0; snrIndex < SNR.size(); snrIndex++)

  outFile.close();
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

// this takes the message bits, including the CRC, and encodes them using the trellis,
// then adds noise and punctures the result
std::vector<double> generateTransmittedMessage(std::vector<int> originalMessage, FeedForwardTrellis encodingTrellis, double snr, std::vector<int> puncturedIndices, bool noiseless){
	// encode the message
	std::vector<int> encodedMessage = encodingTrellis.encode(originalMessage);
	// add noise
	std::vector<double> finalMessage;
	if(noiseless){
		for(int i = 0; i < encodedMessage.size(); i++)
			finalMessage.push_back((double)encodedMessage[i]);
	}
	else
		finalMessage = awgn::addNoise(encodedMessage, snr);
	// puncture the bits. it is more convenient to puncture on this side than on the 
	// decoder, so we insert zeros which provide no information to the decoder
	for(int index = 0; index < puncturedIndices.size(); index++)
		finalMessage[puncturedIndices[index]] = 0;
	return finalMessage;
}