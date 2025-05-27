/* K64N128 PARAMETERS */
#ifndef K64N128_PARAM
#define K64N128_PARAM

#include <vector>

/* --- Convolutional Code Parameters --- */
constexpr int k = 1;                    /* Number of input bits */
constexpr int n = 2;                    /* Number of output bits */
constexpr int K = 64;                   /* Number of input bits */
constexpr int N = 128;                  /* Number of output bits */
constexpr int V = 8;                    /* Number of memory elements */
constexpr int M = 12;                   /* Degree of CRC - 1 */
constexpr unsigned int CRC = 0x1565;    /* CRC polynomial */
inline const std::vector<int> NUMERATORS = {561, 753};  /* in octal */

const std::vector<int> PUNCTURING_INDICES 
    = {4, 10, 21, 24, 31, 37, 
       42, 48, 59, 62, 69, 75, 
       80, 86, 97, 100, 107, 113, 
       118, 124, 135, 138, 145, 151};   /* 2023 ISTC paper puncturing pattern */

/* --- List Decoder Parameters --- */
constexpr int MAX_LISTSIZE = 1e7;       /* Maximum list size */
constexpr double MAX_METRIC = 87.6;     /* Maximum decoding metric */
constexpr double MAX_ANGLE = 0.7498;    /* Maximum angle for the list decoder */
inline constexpr char STOPPING_RULE = 'A';     /* Stopping rule: {M: metric, L: listsize, A: angle} */
inline constexpr char ENCODING_RULE = 'T';     /* Encoding rule: {Z: zero-terminated, T: tail-biting} */
constexpr char DECODING_RULE = 'P';     /* Decoding rule: {P: projected, N: non-projected}*/
constexpr char ERROR_RUN_TYPE = 'T';    /* Accumulate to which type of error: {U: undetected, T: total}*/

/* --- Simulation Parameters --- */
constexpr int MAX_ERRORS = 20;           /* Maximum number of errors */
constexpr bool NOISELESS = false;        /* Noiseless simulation */
constexpr int LOGGING_ITERS = 1000;      /* Logging Interval*/
constexpr int BASE_SEED = 42;            /* Fixed base seed for simulation */
inline const std::vector<float> EBN0 = {2.5};

#endif