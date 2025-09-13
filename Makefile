CXX = mpicxx
CXXFLAGS = -std=c++17 -Wall -Wextra -I include -O2
LDFLAGS = 

# Directories
SRC_DIR = src
INCLUDE_DIR = include
BUILD_DIR = build
CONFIG ?= K64N128

# Create a list of all source and object files
SRC_FILES = $(wildcard $(SRC_DIR)/*.cpp)
COMMON_SRC_FILES = $(SRC_DIR)/minHeap.cpp $(SRC_DIR)/namespace.cpp $(SRC_DIR)/feedForwardTrellis.cpp $(SRC_DIR)/lowRateListDecoder_TB.cpp $(SRC_DIR)/lowRateListDecoder_ZT.cpp 
# COMMON_SRC_FILES = $(filter-out $(SRC_DIR)/bald_sim.cpp $(SRC_DIR)/rova_sim.cpp $(SRC_DIR)/collect_sim.cpp $(SRC_DIR)/genieAidedListDecoder_TB.cpp $(SRC_DIR)/lowRateListDecoder_BCJR.cpp, $(SRC_FILES))

BALD_SRC_FILES = $(COMMON_SRC_FILES) $(SRC_DIR)/bald_sim.cpp
ROVA_SRC_FILES = $(COMMON_SRC_FILES) $(SRC_DIR)/rova_sim.cpp
COLLECT_SRC_FILES = $(COMMON_SRC_FILES) $(SRC_DIR)/collect_sim.cpp $(SRC_DIR)/genieAidedListDecoder_TB.cpp

BALD_OBJ_FILES = $(patsubst $(SRC_DIR)/%.cpp, $(BUILD_DIR)/%.o, $(BALD_SRC_FILES))
ROVA_OBJ_FILES = $(patsubst $(SRC_DIR)/%.cpp, $(BUILD_DIR)/%.o, $(ROVA_SRC_FILES))
COLLECT_OBJ_FILES = $(patsubst $(SRC_DIR)/%.cpp, $(BUILD_DIR)/%.o, $(COLLECT_SRC_FILES))

# Executable name
BALD_TARGET = bald
ROVA_TARGET = rova
COLLECT_TARGET = collect

# Default rule
$(BALD_TARGET): clean consts.h $(BALD_OBJ_FILES) 
	$(CXX) $(LDFLAGS) $(BALD_OBJ_FILES) -o $@

$(ROVA_TARGET): clean consts.h $(ROVA_OBJ_FILES) 
	$(CXX) $(LDFLAGS) $(ROVA_OBJ_FILES) -o $@

$(COLLECT_TARGET): clean consts.h $(COLLECT_OBJ_FILES)
	$(CXX) $(LDFLAGS) $(COLLECT_OBJ_FILES) -o $@

consts.h:
	cp $(INCLUDE_DIR)/consts_$(CONFIG).h consts.h

# Rule to compile source files into object files
$(BUILD_DIR)/%.o: $(SRC_DIR)/%.cpp | $(BUILD_DIR)
	$(CXX) $(CXXFLAGS) $(INCLUDES) -DCONFIG_$(CONFIG) -c $< -o $@

# Rule to create the build directory if it doesn't exist
$(BUILD_DIR):
	mkdir -p $(BUILD_DIR)

# Clean up build files
clean:
	rm -rf $(BUILD_DIR) $(BALD_TARGET) $(ROVA_TARGET) $(COLLECT_TARGET) consts.h

# Rule to clean up object files
clean_obj:
	rm -rf $(BUILD_DIR)/*.o