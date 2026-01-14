/* K16N48 PARAMETERS */
#ifndef K21N62_PARAM
#define K21N62_PARAM

#include <vector>
#include <string>


#define __INT_MAX__ 2147483647

/* --- Convolutional Code Parameters --- */
const int kconv = 1;  /* Number of input bits */
const int nconv = 2;  /* Number of output bits */
const int K = 21; /* Number of input bits */
const int N = 62; /* Number of output bits */
const int V = 6;  /* Number of memory elements */
const int M = 10;  /* Degree of CRC */
const int Kcrc = 21;
const int Ncrc = 31;
const int Kconv = 31;
const int Nconv = 62;

const unsigned int CRC = 0b11101101001;           /* CRC polynomial */
const std::vector<int> CRC_VEC = {1,1,1,0,1,1,0,1,0,0,1};
const std::vector<int> NUMERATORS = {133, 171}; /* in octal */
const std::string GENMATRIXFILEPATH = "params/K31N62_TBCC_GenMatrix.txt";
const std::string OUTPUTFILEPATH = "output/K31N62/";
const std::vector<int> PUNCTURING_INDICES = {};

/* --- List Decoder Parameters --- */
constexpr int MAX_LISTSIZE = 1e3;      /* Maximum list size */
constexpr double MAX_METRIC = 26.7669; /* Maximum decoding metric */
const std::vector<float> MAX_METRIC_VEC = {18.54, 20.05, 21.25, 22.35, 23.35, 24.16, 25.16, 25.96};
const std::vector<float> MAX_ANGLE_VEC = {0.515, 0.52, 0.525, 0.53, 0.535, 0.54};
constexpr double MAX_ANGLE = -1.0;  /* Maximum angle for the list decoder */
constexpr char ENCODING_RULE = 'T'; /* code type: {Z: zero-terminated CC, T: tail-biting CC} */
constexpr char DECODING_RULE = 'N'; /* Decoding rule: {P: projected, N: non-projected}*/
constexpr char STOPPING_RULE = 'L'; /* Stopping rule: {M: metric, L: listsize, A: angle, R: rova} */
constexpr char ERROR_RUN_TYPE = 'U'; /* Accumulate which type of error: {U: undetected, T: total}*/

/* --- Simulation Parameters --- */
constexpr int MAX_ERRORS = 200;     /* Maximum number of errors */
constexpr bool NOISELESS = false;    /* Noiseless simulation */
constexpr int LOGGING_ITERS = 10000; /* Logging Interval*/
constexpr int BASE_SEED = 42;       /* Fixed base seed for simulation */
const std::vector<float> EBN0 = {3,3.5, 4, 4.5, 5};
#endif
