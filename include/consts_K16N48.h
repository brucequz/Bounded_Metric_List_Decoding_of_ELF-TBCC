/* K16N48 PARAMETERS */
#ifndef K16N48_PARAM
#define K16N48_PARAM

#include <vector>

#define __INT_MAX__ 2147483647

/* --- Convolutional Code Parameters --- */
const int kconv = 1;  /* Number of input bits */
const int nconv = 2;  /* Number of output bits */
const int K = 16; /* Number of input bits */
const int N = 48; /* Number of output bits */
const int V = 6;  /* Number of memory elements */
const int M = 8;  /* Degree of CRC - 1 */
const int Kcrc = 16;
const int Ncrc = 24;
const int Kconv = 24;
const int Nconv = 48;

const unsigned int CRC = 0b110001001;           /* CRC polynomial */
const std::vector<int> CRC_VEC = {1, 1, 0, 0, 0, 1, 0, 0, 1};
const std::vector<int> NUMERATORS = {133, 171}; /* in octal */

// const std::vector<int> PUNCTURING_INDICES = {3,  4,  9,  10, 15, 16, 21, 22, 27,
//                                              28, 33, 34, 39, 40, 45, 46};
const std::vector<int> PUNCTURING_INDICES = {};

/* --- List Decoder Parameters --- */
constexpr int MAX_LISTSIZE = 1e7;      /* Maximum list size */
constexpr double MAX_METRIC = 26.7669; /* Maximum decoding metric */
const std::vector<float> MAX_METRIC_VEC = {18.54, 20.05, 21.25, 22.35, 23.35, 24.16, 25.16, 25.96};
const std::vector<float> MAX_ANGLE_VEC = {0.515, 0.52, 0.525, 0.53, 0.535, 0.54};
constexpr double MAX_ANGLE = -1.0;  /* Maximum angle for the list decoder */
constexpr char ENCODING_RULE = 'T'; /* code type: {Z: zero-terminated CC, T: tail-biting CC} */
constexpr char DECODING_RULE = 'N'; /* Decoding rule: {P: projected, N: non-projected}*/
constexpr char STOPPING_RULE = 'L'; /* Stopping rule: {M: metric, L: listsize, A: angle, R: rova} */
constexpr char ERROR_RUN_TYPE = 'U'; /* Accumulate which type of error: {U: undetected, T: total}*/

/* --- Simulation Parameters --- */
constexpr int MAX_ERRORS = 100;     /* Maximum number of errors */
constexpr bool NOISELESS = true;    /* Noiseless simulation */
constexpr int LOGGING_ITERS = 1000; /* Logging Interval*/
constexpr int BASE_SEED = 42;       /* Fixed base seed for simulation */
constexpr float EBN0 = 3.5;
#endif
