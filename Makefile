CXX 		 := mpicxx
CXXFLAGS := -O2 -std=c++20 -Wall -Wextra -Iinclude -I/opt/homebrew/include
LDFLAGS  := -L/opt/homebrew/lib -lfmt

# Directories
SRC_DIR = src
INCLUDE_DIR = include
BUILD_DIR = build
CONFIG ?= K64N128

# Create a list of all source and object files
COMMON_SRC_FILES = $(SRC_DIR)/minHeap.cpp $(SRC_DIR)/namespace.cpp $(SRC_DIR)/feedForwardTrellis.cpp $(SRC_DIR)/lowRateListDecoder_TB.cpp $(SRC_DIR)/lowRateListDecoder_ZT.cpp 

BALD_SRC_FILES = $(COMMON_SRC_FILES) $(SRC_DIR)/bald_sim.cpp
ROVA_SRC_FILES = $(COMMON_SRC_FILES) $(SRC_DIR)/rova_sim.cpp
COLLECT_SRC_FILES = $(COMMON_SRC_FILES) $(SRC_DIR)/collect_sim.cpp $(SRC_DIR)/genieAidedListDecoder_TB.cpp
SSD_SRC_FILES = $(COMMON_SRC_FILES) $(SRC_DIR)/ssd_sim.cpp $(SRC_DIR)/linearityDecoding.cpp $(SRC_DIR)/ssdSLVDDecoding.cpp
GEN_SRC_FILES = $(COMMON_SRC_FILES) $(SRC_DIR)/ssd_table_gen.cpp $(SRC_DIR)/linearityDecoding.cpp

BALD_OBJ_FILES = $(patsubst $(SRC_DIR)/%.cpp, $(BUILD_DIR)/%.o, $(BALD_SRC_FILES))
ROVA_OBJ_FILES = $(patsubst $(SRC_DIR)/%.cpp, $(BUILD_DIR)/%.o, $(ROVA_SRC_FILES))
COLLECT_OBJ_FILES = $(patsubst $(SRC_DIR)/%.cpp, $(BUILD_DIR)/%.o, $(COLLECT_SRC_FILES))
SSD_OBJ_FILES = $(patsubst $(SRC_DIR)/%.cpp, $(BUILD_DIR)/%.o, $(SSD_SRC_FILES))
GEN_OBJ_FILES = $(patsubst $(SRC_DIR)/%.cpp, $(BUILD_DIR)/%.o, $(GEN_SRC_FILES))

# Executable name
BALD_TARGET = bald
ROVA_TARGET = rova
COLLECT_TARGET = collect
SSD_TARGET = ssd
GEN_SSD_TABLE_TARGET = gen

# Default rule
$(BALD_TARGET): clean consts.h $(BALD_OBJ_FILES) 
	$(CXX) $(LDFLAGS) $(BALD_OBJ_FILES) -o $@

$(ROVA_TARGET): clean consts.h $(ROVA_OBJ_FILES) 
	$(CXX) $(LDFLAGS) $(ROVA_OBJ_FILES) -o $@

$(COLLECT_TARGET): clean consts.h $(COLLECT_OBJ_FILES)
	echo $(PWD)
	$(CXX) $(LDFLAGS) $(COLLECT_OBJ_FILES) -o $@

$(SSD_TARGET): $(INCLUDE_DIR)/consts_$(CONFIG).h $(SSD_OBJ_FILES) 
	$(CXX) $(LDFLAGS) $(SSD_OBJ_FILES) -o $@

$(GEN_SSD_TABLE_TARGET): $(INCLUDE_DIR)/consts_$(CONFIG).h $(GEN_OBJ_FILES) 
	$(CXX) $(LDFLAGS) $(GEN_OBJ_FILES) -o $@

consts.h: $(INCLUDE_DIR)/consts_$(CONFIG).h
	@echo "Updating consts.h from consts_$(CONFIG).h"
	cp $(INCLUDE_DIR)/consts_$(CONFIG).h consts.h

# Rule to compile source files into object files
$(BUILD_DIR)/%.o: $(SRC_DIR)/%.cpp consts.h | $(BUILD_DIR)
	$(CXX) $(CXXFLAGS) $(INCLUDES) -DCONFIG_$(CONFIG) -c $< -o $@

# Rule to create the build directory if it doesn't exist
$(BUILD_DIR):
	mkdir -p $(BUILD_DIR)

# Clean up build files
clean:
	rm -rf $(BUILD_DIR) $(BALD_TARGET) $(ROVA_TARGET) $(COLLECT_TARGET) $(SSD_TARGET) $(GEN_SSD_TABLE_TARGET) consts.h

# Rule to clean up object files
clean_obj:
	rm -rf $(BUILD_DIR)/*.o