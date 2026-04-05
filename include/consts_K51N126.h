/* K16N48 PARAMETERS */
#ifndef K16N48_PARAM
#define K16N48_PARAM
#include <vector>
#include <string>

/* --- Convolutional Code Parameters --- */
const int kconv = 1;   /* Number of input bits */
const int nconv = 2;   /* Number of output bits */
const int K     = 51;  /* Number of input bits */
const int N     = 126; /* Number of output bits */
const int V     = 6;   /* Number of memory elements */
const int M     = 12;  /* Degree of CRC */
const int Kcrc  = 51;
const int Ncrc  = 63;
const int Kconv = 63;
const int Nconv = 126;

const unsigned int CRC                    = 0b1001110001101; /* CRC polynomial */
const std::vector<int> CRC_VEC            = {1,0,0,1,1,1,0,0,0,1,1,0,1};
const std::vector<int> NUMERATORS         = {133, 171}; /* in octal */
const std::string GENMATRIXFILEPATH       = "params/K51N126_TBCC_GenMatrix.txt";
const std::string OUTPUTFILEPATH          = "output/K51N126/";
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
constexpr int OFFSET_SPHERE_RADIUS = 6;

/* --- Simulation Parameters --- */
constexpr int MAX_ERRORS      = 10;   /* Maximum number of errors */
constexpr bool NOISELESS      = false; /* Noiseless simulation */
constexpr int LOGGING_ITERS   = 1000;  /* Logging Interval*/
constexpr int BASE_SEED       = 42;    /* Fixed base seed for simulation */
const std::vector<float> EBN0 = {3.9};
#endif
