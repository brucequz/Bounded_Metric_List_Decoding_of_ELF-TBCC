#ifndef MLA_TYPES_H
#define MLA_TYPES_H

#include <vector>

struct CodeInformation {
  int k;              // numerator of the rate
  int n;              // denominator of the rate
  int v;              // number of memory elements
  int crcDeg;         // m+1, degree of CRC, # bits of CRC polynomial
  int crc;            // CRC polynomial
  int numInfoBits;    // number of information bits
  
  std::vector<int> numerators; // optimal code numerators
};

struct MessageInformation{
	MessageInformation() {
		message 					= std::vector<int>();
		path 							= std::vector<int>();
		listSize 					= -1;
    TBListSize        = -1;
		listSizeExceeded 	= false;
		metric 						= -1.0;
		angle_received_decoded_rad = 0.0;
		log_Gamma  = 0.0;
		rova_probability = 0.0;

		vec_highrate_to_tx = std::vector<float>(100, 0);
		vec_highrate_to_projected_rx = std::vector<float>(100, 0);
		vec_running_minimum = std::vector<float>(100, 0);
	};
	std::vector<int> message;
	std::vector<int> path;
	int listSize;
  int TBListSize;
	bool listSizeExceeded;
	double metric;
	double angle_received_decoded_rad;
	double log_Gamma;
	double rova_probability;
	std::vector<float> vec_highrate_to_tx;
	std::vector<float> vec_highrate_to_projected_rx;
	std::vector<float> vec_running_minimum;
};

#endif