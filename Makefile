# Compiler and flags
CXX = mpic++
CXXFLAGS = -std=c++17 -DNDEBUG

# Directories
PROJECT_DIR = $(shell pwd)
SRC_DIR = src
DOCS_DIR = docs
INCLUDE_DIR = include
LIB_DIR = lib
BUILD_DIR = build

# PACS directories
PACS_ROOT = /home/jammy/shared-folder/project/pacs-examples
PACS_INC_DIR = $(PACS_ROOT)/Examples/src/Utilities
PACS_EXTRAS_DIR = $(PACS_ROOT)/Extras
PACS_MUPARSER_DIR = $(PACS_EXTRAS_DIR)/muparser/include/
PACS_JSON_DIR = $(PACS_EXTRAS_DIR)/json/include/nlohmann/
PACS_LIB_DIR = $(PACS_ROOT)/Examples/lib

ifneq (,$(shell grep -rl '#include <Eigen/' $(SRC_DIR)))
CXXFLAGS += -I${mkEigenInc}
endif

# Source files
SRC_FILES = $(wildcard $(SRC_DIR)/*.cpp) \
            $(wildcard $(SRC_DIR)/*/*.cpp) \
            main_parallel.cpp

OBJ_FILES = $(SRC_FILES:%.cpp=$(BUILD_DIR)/%.o)

# Target executable
TARGET = main_executable

# Default rule
all: $(TARGET)

# Linking
$(TARGET): $(OBJ_FILES)
	$(CXX) $(OBJ_FILES) -o $(TARGET) $(LDFLAGS)

# Compilation
$(BUILD_DIR)/%.o: %.cpp
	@mkdir -p $(dir $@)
	$(CXX) $(CXXFLAGS) -c $< -o $@

# Generate documentation
docs:
	doxygen Doxyfile

# Clean rule
clean:
	rm -rf $(BUILD_DIR) $(TARGET) $(DOCS_DIR)
	find $(INCLUDE_DIR) -type f ! -name '.gitignore' -delete
	find $(INCLUDE_DIR) -type d -empty -delete
	rm -rf $(LIB_DIR)/*

# Phony targets
.PHONY: all clean docs install
