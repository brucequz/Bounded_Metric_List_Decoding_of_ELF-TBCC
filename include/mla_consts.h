/* Metric - Listsize Analysis (MLA) Constants */
#ifndef MLACONST_H
#define MLACONST_H

#include <vector>

/* --- Convolutional Code Parameters --- */
constexpr int K = 1;                    /* Number of input bits */
constexpr int N = 2;                    /* Number of output bits */
constexpr int V = 8;                    /* Number of memory elements */
constexpr int M = 12;                   /* Degree of CRC - 1 */
constexpr unsigned int CRC = 0x1565;    /* CRC polynomial */
constexpr int POLY1 = 561;              /* Polynomial 1 */
constexpr int POLY2 = 753;              /* Polynomial 2 */
constexpr int NUM_INFO_BITS = 64;       /* Number of information bits */
constexpr int NUM_CODED_SYMBOLS = 128;  /* Number of coded symbols */

const std::vector<int> PUNCTURING_INDICES 
    = {4, 10, 21, 24, 31, 37, 
       42, 48, 59, 62, 69, 75, 
       80, 86, 97, 100, 107, 113, 
       118, 124, 135, 138, 145, 151};   /* 2023 ISTC paper puncturing pattern */

/* --- List Decoder Parameters --- */
constexpr int MAX_LISTSIZE = 1e7;      /* Maximum list size */
constexpr int MAX_METRIC = 100;         /* Maximum decoding metric */
constexpr char STOPPING_RULE = 'L';     /* Stopping rule */

/* --- Simulation Parameters --- */
constexpr int MAX_ERRORS = 2;           /* Maximum number of errors */
constexpr bool NOISELESS = false;       /* Noiseless simulation */
const std::vector<double> EBN0 = {2.5}; /* Eb/N0 values */
constexpr int LOGGING_ITERS = 1000;     /* Logging Interval*/
constexpr int BASE_SEED = 42;           /* Fixed base seed for simulation */

#endif