
#ifndef FEEDFORWARDTRELLIS_H
#define FEEDFORWARDTRELLIS_H

#include <vector>
#include "mla_types.h"

struct FeedForwardTrellis{
  FeedForwardTrellis(CodeInformation code);

  std::vector<int> encode(const std::vector<int>& message);

  std::vector<std::vector<int>> getNextStates() {return nextStates_;};
  std::vector<std::vector<int>> getOutputs() {return output_;};
  int getNumStates() {return numStates_;};
  int getN() {return n_;};
  int getV() {return v_;};


  private:
    int k_;
    int n_;
    int v_;
    int numInputSymbols_;
    int numOutputSymbols_;
    int numStates_;
    int code_rate_;
    std::vector<int> polynomials_;  // generator polynomial in octal
    std::vector<std::vector<int>> nextStates_;
    std::vector<std::vector<int>> output_;

    void computeNextStates();
    void computeOutput();
};

#endif