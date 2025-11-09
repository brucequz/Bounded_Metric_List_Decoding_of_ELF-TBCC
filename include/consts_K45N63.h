/* K64N128 PARAMETERS */
#ifndef K64N80_PARAM
#define K64N80_PARAM

#include <vector>

/* --- Convolutional Code Parameters --- */
inline constexpr int k = 1;                    /* Number of input bits */
inline constexpr int n = 1;                    /* Number of output bits */
inline constexpr int K = 45;                   /* Number of input bits */
inline constexpr int N = 63;                   /* Number of output bits */
inline constexpr int V = 6;                    /* Number of memory elements */
inline constexpr int M = 12;                    /* Degree of CRC - 1 */
inline constexpr unsigned int CRC = 0x1E85;     /* CRC polynomial */
inline const std::vector<int> NUMERATORS = {103};  /* in octal */
inline const std::vector<int> PUNCTURING_INDICES = {};

/* --- List Decoder Parameters --- */
inline constexpr int MAX_LISTSIZE = 1e7;       /* Maximum list size */
inline double MAX_METRIC = 26.7669;     /* Maximum decoding metric */
inline const std::vector<float> MAX_METRIC_VEC = {18.54,	20.05,	21.25,	22.35,	23.35,	24.16,	25.16, 25.96};
inline const std::vector<float> MAX_ANGLE_VEC = {0.515, 0.52, 0.525, 0.53, 0.535, 0.54};
inline double MAX_ANGLE = -1.0;              /* Maximum angle for the list decoder */
inline constexpr char ENCODING_RULE = 'Z';     /* Select code type: {Z: zero-terminated CC, T: tail-biting CC} */
inline constexpr char DECODING_RULE = 'N';     /* Decoding rule: {P: projected, N: non-projected}*/
inline constexpr char STOPPING_RULE = 'L';     /* Stopping rule: {M: metric, L: listsize, A: angle, R: rova} */
inline constexpr char ERROR_RUN_TYPE = 'U';    /* Accumulate to which type of error: {U: undetected, T: total}*/

/* --- Simulation Parameters --- */
inline constexpr int MAX_ERRORS = 100;           /* Maximum number of errors */
inline constexpr bool NOISELESS = false;        /* Noiseless simulation */
inline constexpr int LOGGING_ITERS = 1000;      /* Logging Interval*/
inline constexpr int BASE_SEED = 42;            /* Fixed base seed for simulation */
inline const std::vector<float> EBN0 = {3.5};

#endif