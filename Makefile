OBJS = src/main.o src/feedForwardTrellis.o src/lowRateListDecoder.o
CXX = mpicxx
CXXFLAGS = -std=c++17 -Wall -I include -O2


# Directories
SRC_DIR = src
INCLUDE_DIR = include
BUILD_DIR = build
CONFIG ?= K64N128

# Create a list of all source and object files
SRC_FILES = $(wildcard $(SRC_DIR)/*.cpp)
OBJ_FILES = $(patsubst $(SRC_DIR)/%.cpp, $(BUILD_DIR)/%.o, $(SRC_FILES))

# Executable name
TARGET = main

# Default rule
all: clean consts.h $(TARGET)


consts.h:
	cp $(INCLUDE_DIR)/consts_$(CONFIG).h consts.h

# Rule to link object files and create the final
$(TARGET): $(OBJ_FILES)
		$(CXX) $(CXXFLAGS) -o $(TARGET) $(OBJ_FILES)

# Rule to compile source files into object files
$(BUILD_DIR)/%.o: $(SRC_DIR)/%.cpp | $(BUILD_DIR)
	$(CXX) $(CXXFLAGS) $(INCLUDES) -DCONFIG_$(CONFIG) -c $< -o $@

# Rule to create the build directory if it doesn't exist
$(BUILD_DIR):
	mkdir -p $(BUILD_DIR)

# Clean up build files
clean:
	rm -rf $(BUILD_DIR) $(TARGET) consts.h

# Rule to clean up object files
clean_obj:
	rm -rf $(BUILD_DIR)/*.o