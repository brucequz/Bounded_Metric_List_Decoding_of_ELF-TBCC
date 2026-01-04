#ifndef MLA_TYPES_H
#define MLA_TYPES_H

#include <vector>

struct CodeInformation {
  int kconv;              // numerator of the rate
  int nconv;              // denominator of the rate
  int v;              // number of memory elements
  int crcLen;         // m+1, length of CRC, degree of CRC + 1
  int crc;            // CRC polynomial
  int numInfoBits;    // number of information bits
  
  const std::vector<int> numerators; // optimal code numerators
};

struct MessageInformation{
	MessageInformation() {
		message 					= std::vector<int>();
		codeword 					= std::vector<int>();
		path 							= std::vector<int>();
		listSize 					= -1;
    TBListSize        = -1;
		listSizeExceeded 	= false;
		metric 						= -1.0;
	};
	std::vector<int> message;
	std::vector<int> codeword;
	std::vector<int> path;
	int listSize;
  int TBListSize;
	bool listSizeExceeded;
	double metric;
	float angle_received_decoded_rad;
};

#endif