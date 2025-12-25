#ifndef FEEDFORWARDTRELLIS_H
#define FEEDFORWARDTRELLIS_H
#include "../consts.h"
#include "namespace.h"
#include "types.h"

#include <string>
#include <vector>

class FeedForwardTrellis {
public:
  FeedForwardTrellis(int k, int n, int v, std::vector<int> numerators);

  std::vector<int> encode(const std::vector<int>& originalMessage) const;
  std::vector<int> encode_zt(const std::vector<int>& originalMessage) const;

  std::vector<std::vector<int>> getNextStates();
  std::vector<std::vector<int>> getOutputs();

  /* generator matrix for TBCC code */
  void computeGeneratorMatrix();

  int getNumInputSymbols();
  int getNumOutputSymbols();
  int getNumStates();
  int getV();
  int getN();

private:
  int k;
  int n;
  int v;
  int numInputSymbols;
  int numOutputSymbols;
  int numStates;
  std::vector<int> numerators;
  std::vector<std::vector<int>> nextStates;
  std::vector<std::vector<int>> outputs;
  std::vector<std::vector<int>> generatorMatrix;

  void computeNextStates();

  std::vector<int> dec2Bin(int decimal, int length);
  int bin2Dec(std::vector<int> binary);
};
#endif
