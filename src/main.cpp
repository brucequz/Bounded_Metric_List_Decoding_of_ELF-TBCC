#include <iostream>
#include <vector>

#include "../include/mla_consts.h"
#include "../include/mla_types.h"
#include "../include/feedForwardTrellis.h"
#include "../include/lowRateListDecoder.h"

void ISTC_sim(CodeInformation code);

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

  return 0;
}

void ISTC_sim(CodeInformation code){
  /* Runs simulation for the code specified in ISTC 2023 https://ieeexplore.ieee.org/document/10273535
      (128, 64) tail-biting convolutional code with ELF = 0x1565 (m=12), and puncturing pattern from the paper.
      The simulation is run for the EBN0 values specified in the MLA_CONSTS file.



  Args:
    code: CodeInformation struct containing the code parameters

  */

	std::vector<double> EbN0(EBN0.begin(), EBN0.end());
	std::vector<double> SNR;
	double offset = 10 * log10((double)2*NUM_INFO_BITS / (128)); 		// real rate of this code *2 in the log
	for (int i=0; i< EbN0.size(); i++)
		SNR.push_back(EbN0[i] + offset);

  FeedForwardTrellis encodingTrellis(code.k, code.n, code.v, code.numerators);
	LowRateListDecoder listDecoder(encodingTrellis, MAX_LISTSIZE, code.crcDeg, code.crc);

// 	// running the simulations. in this example, we are simulation to collect TFR data
// 	for(int snrIndex = 0; snrIndex < SNR.size(); snrIndex++){
// 		double snr = SNR[snrIndex];
// 		double standardMeanListSize = 0;
// 		int standardNumErrors = 0;
// 		int standardListSizeExceeded = 0;

// 		for(int numTrials = 0; numTrials < maxNumTrials; numTrials++){
// 			if (numTrials % 1000 == 0){
// 				std::cout << "currently at " << numTrials << std::endl;
// 			}
// 			// the following lines generate a message that we can decode
// 			std::vector<int> originalMessage = generateRandomCRCMessage(code);
// 			std::vector<double> transmittedMessage = generateTransmittedMessage(originalMessage, encodingTrellis, snr, puncturedIndices);

// 			LowRateListDecoder::messageInformation standardDecoding = listDecoder.lowRateDecoding(transmittedMessage);

// 			if(standardDecoding.listSizeExceeded){
// 				standardListSizeExceeded++;
// 			}
// 			else if (standardDecoding.message != originalMessage){
// 				standardNumErrors ++;
// 				standardMeanListSize += (double)standardDecoding.listSize;
// 			}
			
// 		}
// 		// outputting results at each SNR / EbN0 point. in this case, since we're working with TFR, the ratio
// 		// of incorrect decodings to total decodings is used.
// 		std::cout << "at Eb/N0 = " << EbN0[snrIndex] << std::endl;
// 		std::cout << "mean list size: " << (double)standardMeanListSize/(maxNumTrials - standardListSizeExceeded) << std::endl;
// 		std::cout << "number of times list size exceeded: " << standardListSizeExceeded << std::endl;
// 		std::cout << "UER: " << (double)standardNumErrors/(maxNumTrials - standardListSizeExceeded) << std::endl;
// 		std::cout << "TFR: " << (double)(standardNumErrors + standardListSizeExceeded)/maxNumTrials << std::endl;


// 		// it may be the case that we want to write the outputs to a file, whether it be for further processing
// 		// or convenience. the following code snipped illustrates how this would be done. typically, the opening
// 		// and closing of the file, and declaration of the file name would be done once each outside the loop,
// 		// but they are included here to keep things organized.
// 		/*
// 		std::ofstream outputFile;
// 		filename = "distance_spectra.txt";
// 		outputFile.open(filename, fstream::app);
// 		outputFile << "example information to be saved";
// 		outputFile.close();
// 		*/
// 	}
// 	std::cout << "concluded simulation" << std::endl;

}