/* K16N48 PARAMETERS */
#ifndef K15N30_PARAM
#define K15N30_PARAM

#include <string>
#include <vector>

/* --- Convolutional Code Parameters --- */
const int kconv = 1;  /* Number of input bits */
const int nconv = 2;  /* Number of output bits */
const int K     = 11; /* Number of input bits */
const int N     = 30; /* Number of output bits */
const int V     = 6;  /* Number of memory elements */
const int M     = 4;  /* Degree of CRC - 1 */
const int Kcrc  = 11;
const int Ncrc  = 15;
const int Kconv = 15;
const int Nconv = 30;

const unsigned int CRC              = 0b10011; /* CRC polynomial */
const std::vector<int> CRC_VEC      = {1, 0, 0, 1, 1};
const std::vector<int> NUMERATORS   = {133, 171}; /* in octal */
const std::string GENMATRIXFILEPATH = "params/K15N30_TBCC_GenMatrix.txt";
const std::string OUTPUTFILEPATH    = "output/K11N30/";

// const std::vector<int> PUNCTURING_INDICES = {3,  4,  9,  10, 15, 16, 21, 22, 27,
//                                              28, 33, 34, 39, 40, 45, 46};
const std::vector<int> PUNCTURING_INDICES = {};

/* --- List Decoder Parameters --- */
constexpr int MAX_LISTSIZE              = 1e8;     /* Maximum list size */
constexpr double MAX_METRIC             = 26.7669; /* Maximum decoding metric */
const std::vector<float> MAX_METRIC_VEC = {18.54, 20.05, 21.25, 22.35, 23.35, 24.16, 25.16, 25.96};
const std::vector<float> MAX_ANGLE_VEC  = {0.515, 0.52, 0.525, 0.53, 0.535, 0.54};
constexpr double MAX_ANGLE              = -1.0; /* Maximum angle for the list decoder */
constexpr char ENCODING_RULE            = 'T';  /* code type: {Z: zero-terminated CC, T: tail-biting CC} */
constexpr char DECODING_RULE            = 'N';  /* Decoding rule: {P: projected, N: non-projected}*/
constexpr char STOPPING_RULE            = 'L';  /* Stopping rule: {M: metric, L: listsize, A: angle, R: rova} */
constexpr char ERROR_RUN_TYPE           = 'U';  /* Accumulate which type of error: {U: undetected, T: total}*/

/* --- SSD Decoder Parameters --- */
constexpr int OFFSET_SPHERE_RADIUS = 12;

/* --- Simulation Parameters --- */
constexpr int MAX_ERRORS      = 100;   /* Maximum number of errors */
constexpr bool NOISELESS      = false; /* Noiseless simulation */
constexpr int LOGGING_ITERS   = 1000;  /* Logging Interval*/
constexpr int BASE_SEED       = 42;    /* Fixed base seed for simulation */
const std::vector<float> EBN0 = {3, 3.5, 4, 4.5, 5, 5.5, 6};
#endif
