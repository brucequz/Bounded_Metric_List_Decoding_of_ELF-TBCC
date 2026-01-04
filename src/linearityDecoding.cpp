#include "../include/lowRateListDecoder.h"

#include <cstddef>
#include <exception>
#include <format>
#include <fstream>

void LowRateListDecoder::generateNeighborList_sequential(const std::vector<float>& allZerosMessage,
                                                         std::string dir, double thre) {
  // 	run standard SLVD starting from all zeros codeword

  // trellisInfo is indexed [state][stage]
  std::vector<std::vector<cell>> trellisInfo;
  trellisInfo = constructLowRateTrellis(allZerosMessage);
  std::cout << "threshold: " << thre << std::endl;

  // store neighbor messages and codewords of a CERTAIN DISTANCE
  // the outer array has (Numstates * 2^m) rows, we iterate through ESD first
  // ESD=0, ED=0; ESD=0, ED=1; ESD=0, ED=2 ...
  // inner array have Kconv/Nconv * x, x=0,1,2,3…
  // Initialize each row with placeholder 3, remove the space-holder once we find a neighbor with
  // the ELF-ED pair
  int numRows = lowrate_numStates * std::pow(2, crcLength - 1);
  std::vector<std::vector<int>> neighbor_msg(numRows, std::vector<int>{3});
  std::vector<std::vector<int>> neighbor_cwd(numRows, std::vector<int>{3});
  int cur_distance = 0;
  int cur_paths = 0;
  int temp = 0;

  // start search
  MessageInformation output;
  // RBTree detourTree;
  MinHeap detourTree;
  std::vector<std::vector<int>> previousPaths;

  // create nodes for each valid ending state with no detours
  // std::cout<< "end path metrics:" <<std::endl;
  for (int i = 0; i < lowrate_numStates; i++) {
    DetourObject detour;
    detour.startingState = i;
    detour.pathMetric = trellisInfo[i][lowrate_pathLength - 1].pathMetric;
    detourTree.insert(detour);
  }

  int numPathsSearched = 0;
  while (numPathsSearched < this->listSize) {
    DetourObject detour = detourTree.pop();
    std::vector<int> path(lowrate_pathLength);

    int newTracebackStage = lowrate_pathLength - 1;
    double forwardPartialPathMetric = 0;
    int currentState = detour.startingState;

    // if we are taking a detour from a previous path, we skip backwards to the point where we take
    // the detour from the previous path
    if (detour.originalPathIndex != -1) {
      forwardPartialPathMetric = detour.forwardPathMetric;
      newTracebackStage = detour.detourStage;
      // std::cout <<"newTracebackStage"<< newTracebackStage<<std::endl;

      // while we only need to copy the path from the detour to the end, this simplifies things,
      // and we'll write over the earlier data in any case
      path = previousPaths[detour.originalPathIndex];
      currentState = path[newTracebackStage];
      // std::cout <<"currentState"<< currentState<<std::endl;

      double suboptimalPathMetric =
          trellisInfo[currentState][newTracebackStage].suboptimalPathMetric;

      currentState = trellisInfo[currentState][newTracebackStage].suboptimalFatherState;
      newTracebackStage--;

      double prevPathMetric = trellisInfo[currentState][newTracebackStage].pathMetric;

      forwardPartialPathMetric += suboptimalPathMetric - prevPathMetric;
    }
    path[newTracebackStage] = currentState;

    // actually tracing back
    for (int stage = newTracebackStage; stage > 0; stage--) {
      double suboptimalPathMetric = trellisInfo[currentState][stage].suboptimalPathMetric;
      double currPathMetric = trellisInfo[currentState][stage].pathMetric;

      // if there is a detour we add to the detourTree
      if (trellisInfo[currentState][stage].suboptimalFatherState != -1) {
        DetourObject localDetour;
        localDetour.detourStage = stage;
        localDetour.originalPathIndex = numPathsSearched;
        localDetour.pathMetric = suboptimalPathMetric + forwardPartialPathMetric;
        localDetour.forwardPathMetric = forwardPartialPathMetric;
        localDetour.startingState = detour.startingState;
        detourTree.insert(localDetour);
        // std::cout << "detour inserted" <<std::endl;
      }
      currentState = trellisInfo[currentState][stage].optimalFatherState;
      double prevPathMetric = trellisInfo[currentState][stage - 1].pathMetric;
      forwardPartialPathMetric += currPathMetric - prevPathMetric;
      path[stage - 1] = currentState;
    }

    previousPaths.push_back(path);

    std::vector<int> message = pathToMessage(path);
    std::vector<int> codeword = pathToCodeword(path);
    // print_int_vector(message);
    // print_int_vector(codeword);
    // puncture and recalculate metric
    for (size_t index = 0; index < PUNCTURING_INDICES.size(); index++) {
      codeword[PUNCTURING_INDICES[index]] = 0;
    }

    // print_int_vector(codeword);
    int metric = 0;
    for (size_t i = 0; i < codeword.size(); i++) {
      if (codeword[i] != 0)
        metric += std::pow((codeword[i] - allZerosMessage[i]), 2);
    }
    // std::cout << metric << std::endl;

    // Once we arrive at a new squared euclidean distance, write the current stored neighbors to 2
    // text files, one for message and one for codeword
    if (metric > cur_distance) {
      // Keep track of how many paths we explored when we arrive at a new squared euclidean distance
      cur_paths = numPathsSearched - temp;
      temp = numPathsSearched;
      std::cout << "current distance: " << cur_distance << ", num path (L_d): " << cur_paths
                << std::endl;
      // write to text file
      std::string msgFilename = "output/m_" + dir + std::to_string(cur_distance) + ".txt";
      std::string cwdFilename = "output/c_" + dir + std::to_string(cur_distance) + ".txt";
      writeDataToFile(neighbor_msg, msgFilename);
      writeDataToFile(neighbor_cwd, cwdFilename);

      // re-initialize neighbor_msg & neighbor_cwd
      for (auto& row : neighbor_msg) {
        row.clear();
        row.push_back(3);
      }
      for (auto& row : neighbor_cwd) {
        row.clear();
        row.push_back(3);
      }
      // update current distance
      cur_distance = metric;
    }

    // terminate when we reach a certain distance
    if (metric > thre) {
      output.listSizeExceeded = false;
      std::cout << "found all neighbors up to distance " << (int) thre
                << ", total num paths: " << numPathsSearched << std::endl;
      return;
    }

    // calculate EDS, ED to store message and codeword neighbors
    int ESD = path[0] ^ path[lowrate_pathLength - 1];
    int ED = crc::remdr(message, CRC, Ncrc, M);
    int neighbor_idx = ESD * pow(2, crcLength - 1) + ED;
    // add a neighbor to a specific row
    if (neighbor_msg[neighbor_idx].size() == 1 && neighbor_msg[neighbor_idx][0] == 3) {
      neighbor_msg[neighbor_idx].clear(); // remove the placeholder
      neighbor_cwd[neighbor_idx].clear();
    }
    // append new message & codeword
    neighbor_msg[neighbor_idx].insert(neighbor_msg[neighbor_idx].end(), message.begin(),
                                      message.end());
    neighbor_cwd[neighbor_idx].insert(neighbor_cwd[neighbor_idx].end(), codeword.begin(),
                                      codeword.end());

    // // // one trellis decoding requires both a tb and crc check
    // // if(path[0] == path[lowrate_pathLength - 1] && crc_check(message, crcLength, crc)){
    // if(numPathsSearched > 1048574){
    // 	output.message = message;
    // 	output.codeword = codeword;
    // 	output.path = path;
    //  	output.listSize = numPathsSearched + 1;
    // 	output.metric = forwardPartialPathMetric;
    // 	// std::cout << output.listSize << std::endl;
    // 	// print_int_vector(message);
    // 	output.TBListSize = TBPathsSearched + 1;
    //  	// return output;
    // }

    numPathsSearched++;
    if (numPathsSearched % 100000 == 0) {
      std::cout << "cur path: " << numPathsSearched << std::endl;
    }
  }
  output.listSizeExceeded = true;
  return;
}

void LowRateListDecoder::generateNeighborList_sequential_TBonly(
    const std::vector<float>& allZerosMessage, std::string dir, double thre) {

  // trellisInfo is indexed [state][stage]
  std::vector<std::vector<cell>> trellisInfo;
  trellisInfo = constructLowRateTrellis_Punctured(allZerosMessage, PUNCTURING_INDICES);
  std::cout << "threshold: " << thre << std::endl;

  // store neighbor messages and codewords of a CERTAIN DISTANCE
  // all offsets are TB, so there are only 2^m rows, one for each ELF difference
  // inner array have Kconv/Nconv * x, x=0,1,2,3…
  // Initialize each row with placeholder 3, remove the space-holder once we find a neighbor with
  // the ELF-ED pair
  int numRows = std::pow(2, crcLength - 1);
  std::vector<std::vector<int>> neighbor_msg(numRows, std::vector<int>{3});
  std::vector<std::vector<int>> neighbor_cwd(numRows, std::vector<int>{3});
  int cur_distance = 0;
  int cur_paths = 0;
  int temp = 0;

  // start search
  MessageInformation output;
  MinHeap detourTree;
  std::vector<std::vector<int>> previousPaths;

  // create nodes for each valid ending state with no detours
  // std::cout<< "end path metrics:" <<std::endl;
  for (int i = 0; i < lowrate_numStates; i++) {
    DetourObject detour;
    detour.startingState = i;
    detour.pathMetric = trellisInfo[i][lowrate_pathLength - 1].pathMetric;
    detourTree.insert(detour);
  }

  int numPathsSearched = 0;
  int TBPathsSearched = 0;
  while (TBPathsSearched < this->listSize) {
    DetourObject detour = detourTree.pop();
    std::vector<int> path(lowrate_pathLength);

    int newTracebackStage = lowrate_pathLength - 1;
    double forwardPartialPathMetric = 0;
    int currentState = detour.startingState;

    // if we are taking a detour from a previous path, we skip backwards to the point where we take
    // the detour from the previous path
    if (detour.originalPathIndex != -1) {
      forwardPartialPathMetric = detour.forwardPathMetric;
      newTracebackStage = detour.detourStage;
      // std::cout <<"newTracebackStage"<< newTracebackStage<<std::endl;

      // while we only need to copy the path from the detour to the end, this simplifies things,
      // and we'll write over the earlier data in any case
      path = previousPaths[detour.originalPathIndex];
      currentState = path[newTracebackStage];
      // std::cout <<"currentState"<< currentState<<std::endl;

      double suboptimalPathMetric =
          trellisInfo[currentState][newTracebackStage].suboptimalPathMetric;

      currentState = trellisInfo[currentState][newTracebackStage].suboptimalFatherState;
      newTracebackStage--;

      double prevPathMetric = trellisInfo[currentState][newTracebackStage].pathMetric;

      forwardPartialPathMetric += suboptimalPathMetric - prevPathMetric;
    }
    path[newTracebackStage] = currentState;

    // actually tracing back
    for (int stage = newTracebackStage; stage > 0; stage--) {
      double suboptimalPathMetric = trellisInfo[currentState][stage].suboptimalPathMetric;
      double currPathMetric = trellisInfo[currentState][stage].pathMetric;

      // if there is a detour we add to the detourTree
      if (trellisInfo[currentState][stage].suboptimalFatherState != -1) {
        DetourObject localDetour;
        localDetour.detourStage = stage;
        localDetour.originalPathIndex = numPathsSearched;
        localDetour.pathMetric = suboptimalPathMetric + forwardPartialPathMetric;
        localDetour.forwardPathMetric = forwardPartialPathMetric;
        localDetour.startingState = detour.startingState;
        detourTree.insert(localDetour);
        // std::cout << "detour inserted" <<std::endl;
      }
      currentState = trellisInfo[currentState][stage].optimalFatherState;
      double prevPathMetric = trellisInfo[currentState][stage - 1].pathMetric;
      forwardPartialPathMetric += currPathMetric - prevPathMetric;
      path[stage - 1] = currentState;
    }

    previousPaths.push_back(path);
    std::vector<int> message = pathToMessage(path);
    std::vector<int> codeword = pathToCodeword(path);

    /* recording TB only */
    if (path[0] == path[lowrate_pathLength - 1]) {

      /* puncturing */
      for (size_t index = 0; index < PUNCTURING_INDICES.size(); index++) {
        codeword[PUNCTURING_INDICES[index]] = 0;
      }
      int metric = 0;
      for (size_t i = 0; i < codeword.size(); i++) {
        if (codeword[i] != 0)
          metric += std::pow((codeword[i] - allZerosMessage[i]), 2);
      }
      std::cout << "forwardPartialMetric: " << forwardPartialPathMetric << std::endl;
      std::cout << "TB codeword " << TBPathsSearched << ", path[0] = " << path[0]
                << ", metric = " << metric << std::endl;
      utils::print_int_vector(codeword);

      // Once we arrive at a new squared euclidean distance, write the current stored neighbors to 2
      // text files, one for message and one for codeword
      if (metric > cur_distance) {
        // Keep track of how many paths we explored when we arrive at a new squared euclidean
        // distance
        cur_paths = TBPathsSearched - temp;
        temp = TBPathsSearched;
        std::cout << "current distance: " << cur_distance << ", num TB offsets: " << cur_paths
                  << std::endl;

        /* record to text file */
        // std::string msgFilename = "m_" + dir + std::to_string(cur_distance) + "_TB.txt";
        // std::string cwdFilename = "c_" + dir + std::to_string(cur_distance) + "_TB.txt";
        // writeDataToFile(neighbor_msg, msgFilename);
        // writeDataToFile(neighbor_cwd, cwdFilename);

        // // re-initialize neighbor_msg & neighbor_cwd
        // for (auto& row : neighbor_msg) {
        // 	row.clear();
        // 	row.push_back(3);
        // }
        // for (auto& row : neighbor_cwd) {
        // 	row.clear();
        // 	row.push_back(3);
        // }

        /* update cur_distance */
        cur_distance = metric;
      } // if (metric > cur_distance)

      // terminate when we reach a certain distance
      if (metric > thre) {
        output.listSizeExceeded = false;
        std::cout << "found all neighbors up to distance " << (int) thre
                  << ", total TB paths: " << TBPathsSearched << std::endl;
        return;
      }

      /*
      // calculate ED to store message and codeword neighbors
      int ED = crc_remainder(message, crcLength, crc);
      int neighbor_idx = ED;
      // add a neighbor to a specific row
      if (neighbor_msg[neighbor_idx].size() == 1 && neighbor_msg[neighbor_idx][0] == 3) {
              neighbor_msg[neighbor_idx].clear();  // remove the placeholder
              neighbor_cwd[neighbor_idx].clear();
      }
      // append new message & codeword
      neighbor_msg[neighbor_idx].insert(neighbor_msg[neighbor_idx].end(), message.begin(),
      message.end()); neighbor_cwd[neighbor_idx].insert(neighbor_cwd[neighbor_idx].end(),
      codeword.begin(), codeword.end());
      */
      TBPathsSearched++;
    } // if (path[0] == path[lowrate_pathLength - 1])

    numPathsSearched++;
    if (numPathsSearched % 100000 == 0) {
      std::cout << "numPathsSearched: " << numPathsSearched << std::endl;
      std::cout << "TBPathsSearched: " << TBPathsSearched << std::endl;
    }

  } // while(TBPathsSearched < this->listSize)
  output.listSizeExceeded = true;
  return;
}

MessageInformation
LowRateListDecoder::lowrateDecoding_neighbors(const std::vector<float>& receivedMessage,
                                              std::vector<std::vector<int>> G_mat) {
  std::vector<std::vector<cell>> trellisInfo;
  trellisInfo = constructLowRateTrellis_Punctured(receivedMessage, PUNCTURING_INDICES);
  std::string neighborSpectraFileName = OUTPUTFILEPATH + "K24N48_Neighbor.txt";
  std::ofstream neighborSpectra(neighborSpectraFileName, std::ios::out);

  // start search
  MessageInformation output;
  MinHeap detourTree;
  std::vector<std::vector<int>> previousPaths;

  for (int i = 0; i < lowrate_numStates; i++) {
    DetourObject detour;
    detour.startingState = i;
    detour.pathMetric = trellisInfo[i][lowrate_pathLength - 1].pathMetric;
    detourTree.insert(detour);
  }
  int numPathsSearched = 0;
  int numTBPathsSearched = 0;
  // // for printing codewords
  // std::ofstream outputCodewords; // distance between the received word and the decoding center
  // std::string filename = "v14_codewords_ZT.txt";
  // outputCodewords.open(filename, std::fstream::app);
  // std::ofstream outputfile_m; // distance between the received word and the decoding center
  // std::string filename_m = "v14_messages_ZT.txt";
  // outputfile_m.open(filename_m, std::fstream::app);

  int currmetric = 0;
  int neighbors = 0;
  int non_neighbors = 0;
  std::vector<std::vector<int>> neighbor_codewords;
  std::vector<std::vector<int>> neighbor_messages;

  while (currmetric < 14) {
    int last_curr_metric = currmetric;
    DetourObject detour = detourTree.pop();
    std::vector<int> path(lowrate_pathLength);
    int newTracebackStage = lowrate_pathLength - 1;
    double forwardPartialPathMetric = 0;
    int currentState = detour.startingState;
    // if we are taking a detour from a previous path, we skip backwards to the point where we take
    // the detour from the previous path
    if (detour.originalPathIndex != -1) {
      forwardPartialPathMetric = detour.forwardPathMetric;
      newTracebackStage = detour.detourStage;
      // while we only need to copy the path from the detour to the end, this simplifies things,
      // and we'll write over the earlier data in any case
      path = previousPaths[detour.originalPathIndex];
      currentState = path[newTracebackStage];
      double suboptimalPathMetric =
          trellisInfo[currentState][newTracebackStage].suboptimalPathMetric;
      currentState = trellisInfo[currentState][newTracebackStage].suboptimalFatherState;
      newTracebackStage--;
      double prevPathMetric = trellisInfo[currentState][newTracebackStage].pathMetric;
      forwardPartialPathMetric += suboptimalPathMetric - prevPathMetric;
    }
    path[newTracebackStage] = currentState;
    // actually tracing back
    for (int stage = newTracebackStage; stage > 0; stage--) {
      double suboptimalPathMetric = trellisInfo[currentState][stage].suboptimalPathMetric;
      double currPathMetric = trellisInfo[currentState][stage].pathMetric;
      // if there is a detour we add to the detourTree
      if (trellisInfo[currentState][stage].suboptimalFatherState != -1) {
        DetourObject localDetour;
        localDetour.detourStage = stage;
        localDetour.originalPathIndex = numPathsSearched;
        localDetour.pathMetric = suboptimalPathMetric + forwardPartialPathMetric;
        localDetour.forwardPathMetric = forwardPartialPathMetric;
        localDetour.startingState = detour.startingState;
        detourTree.insert(localDetour);
      }
      currentState = trellisInfo[currentState][stage].optimalFatherState;
      double prevPathMetric = trellisInfo[currentState][stage - 1].pathMetric;
      forwardPartialPathMetric += currPathMetric - prevPathMetric;
      path[stage - 1] = currentState;
    }
    previousPaths.push_back(path);
    std::vector<int> message = pathToMessage(path);
    std::vector<int> codeword = pathToCodeword(path);

    /* - Recoding TB Gabriel Neighbors */
    if (path.front() == path.back()) {
      numTBPathsSearched++;
      // std::cout << "codeword " << numTBPathsSearched;
      // std::cout << ": path[0]: " << path[0] << "; path[Kconv]: " <<
      // path[this->lowrate_pathLength-1] << std::endl; std::cout << "printing candidate message: ";
      // utils::print_int_vector(message); std::cout << "printing candidate codeword " <<
      // numTBPathsSearched << ": "; utils::print_int_vector(codeword);

      /* Convert to on-off keying codewords + identify location of zeros */
      std::vector<int> ookcw;
      std::vector<int> positions_of_zeros;
      for (size_t i_iter = 0; i_iter < codeword.size(); i_iter++) {
        if (codeword[i_iter] < 0) {
          ookcw.push_back(1);
        } else {
          ookcw.push_back(0);
          positions_of_zeros.push_back(i_iter);
        }
      }

      /* Compute hamming weight*/
      int dist_from_all_zeros = 0;
      for (size_t i_iter = 0; i_iter < ookcw.size(); i_iter++) {
        dist_from_all_zeros = dist_from_all_zeros + ookcw[i_iter];
      }
      currmetric = dist_from_all_zeros;

      /* Output to file if the latest codeword has a larger hamming weight */
      if (currmetric != last_curr_metric) {
        std::string codewordFileName = OUTPUTFILEPATH + std::format("codeword_{}.txt", last_curr_metric);
        std::string messageFileName = OUTPUTFILEPATH + std::format("message_{}.txt", last_curr_metric);
        try {
          OutputFile CwdFile(codewordFileName);
          CwdFile.write2DVector(neighbor_codewords);
          neighbor_codewords.clear();

          OutputFile MsgFile(messageFileName);
          MsgFile.write2DVector(neighbor_messages);
          neighbor_messages.clear();
        } catch (const std::exception& e) {
          std::cerr << "Error: " << e.what() << std::endl;
        }

        neighborSpectra << "Neighbors of Hamming Weight " << last_curr_metric << ": " << neighbors
                        << std::endl;
        neighborSpectra << "Non-Neighbors of Hamming Weight " << last_curr_metric << ": "
                        << non_neighbors << std::endl;
        neighborSpectra << currmetric << " " << std::endl;
        neighbors = 0;
        non_neighbors = 0;

        std::cout << "current metric: " << currmetric << std::endl;
      }

      /* Agrell's algorithm to check if latest codeword is a Gabriel neighbor of all-zero */
      int app_i = 0;
      int pivot;
      std::vector<std::vector<int>> G_mat_temp = G_mat;
      int is_neighbor = 0;

      for (size_t p_iter = 0; p_iter < positions_of_zeros.size(); p_iter++) {
        int p = positions_of_zeros[p_iter];
        bool firstone = true;
        for (int j = app_i; j < Kconv; j++) {
          if (G_mat_temp[j][p] == 1) {
            if (firstone) {
              pivot = j;
              firstone = false;
            } else {
              std::vector<int> g_j = G_mat_temp[j];
              std::vector<int> g_pivot = G_mat_temp[pivot];
              for (size_t pivot_iter = 0; pivot_iter < g_j.size(); ++pivot_iter) {
                g_j[pivot_iter] ^= g_pivot[pivot_iter]; // a = a XOR b
              }
              G_mat_temp[j] = g_j;
            }
          }
        }
        if (firstone == false) {
          std::vector<int> temp = G_mat_temp[app_i];
          G_mat_temp[app_i] = G_mat_temp[pivot];
          G_mat_temp[pivot] = temp;
          app_i = app_i + 1;
          if (app_i == Kconv - 1) {
            neighbor_messages.push_back(message);
            neighbor_codewords.push_back(ookcw);
            neighbors++;
            is_neighbor = 1;
            break;
          }
        }
      }
      if (is_neighbor == 0) {
        non_neighbors++;
      }
    }

    numPathsSearched++;
  }
  neighborSpectra.close();
  return output;
}

bool LowRateListDecoder::isGabrielNeighbor(std::vector<std::vector<int>> G_mat,
                                           const std::vector<int>& modulated_codeword) const {

  /* demodulation */
  std::vector<int> ookcw;
  std::vector<int> positions_of_zeros;
  for (size_t i_iter = 0; i_iter < modulated_codeword.size(); i_iter++) {
    if (modulated_codeword[i_iter] < 0) {
      ookcw.push_back(1);
    } else {
      ookcw.push_back(0);
      positions_of_zeros.push_back(i_iter);
    }
  }

  int app_i = 0;
  int pivot;
  std::vector<std::vector<int>> G_mat_temp = G_mat;
  bool is_neighbor = false;

  for (size_t p_iter = 0; p_iter < positions_of_zeros.size(); p_iter++) {
    int p = positions_of_zeros[p_iter];
    bool firstone = true;
    for (int j = app_i; j < Kconv; j++) {
      if (G_mat_temp[j][p] == 1) {
        if (firstone) {
          pivot = j;
          firstone = false;
        } else {
          std::vector<int> g_j = G_mat_temp[j];
          std::vector<int> g_pivot = G_mat_temp[pivot];
          for (size_t pivot_iter = 0; pivot_iter < g_j.size(); ++pivot_iter) {
            g_j[pivot_iter] ^= g_pivot[pivot_iter]; // a = a XOR b
          }
          G_mat_temp[j] = g_j;
        }
      }
    }
    if (firstone == false) {
      std::vector<int> temp = G_mat_temp[app_i];
      G_mat_temp[app_i] = G_mat_temp[pivot];
      G_mat_temp[pivot] = temp;
      app_i = app_i + 1;
      if (app_i == Kconv - 1) {
        is_neighbor = true;
        break;
      }
    }
  }
  return is_neighbor;
}

void LowRateListDecoder::writeDataToFile(const std::vector<std::vector<int>>& data,
                                         const std::string& filename) {
  std::ofstream outFile(filename, std::ios::app);
  if (!outFile) {
    std::cerr << "Error opening file: " << filename << std::endl;
    return;
  }

  for (const auto& row : data) {
    if (row.empty())
      continue; // Skip empty rows
    for (size_t i = 0; i < row.size(); ++i) {
      outFile << row[i];
      if (i < row.size() - 1)
        outFile << " "; // separate values with space
    }
    outFile << "\n";
  }
  // std::cout << "Successfully written to text file" << std::endl;
  outFile.close();
}