TARGET=Decompose

SRC_DIR=src
OBJ_DIR=obj
BIN_DIR=bin
HEADER_DIR=headers

UTILS_DIR = $(SRC_DIR)/utils
UTILS_OBJ_DIR = $(OBJ_DIR)/utils

CXX_VERSION = c++17
CFLAGS= -Wall -I $(HEADER_DIR)
LIB= -L/usr/lib
EXTRA= -lmpfr -lgmp -fopenmp

SRCS_UTILS = $(wildcard $(UTILS_DIR)/*.cpp)
OBJS_UTILS = $(SRCS_UTILS:$(UTILS_DIR)/%.cpp=$(UTILS_OBJ_DIR)/%.o)
SRCS = $(wildcard $(SRC_DIR)/*.cpp)
OBJS = $(SRCS:$(SRC_DIR)/%.cpp=$(OBJ_DIR)/%.o)

.PHONY: all clean clean_all

all: $(BIN_DIR)/$(TARGET)

$(BIN_DIR)/$(TARGET): $(OBJS_UTILS) $(OBJS) | $(BIN_DIR)
	g++ -std=$(CXX_VERSION) $(LIB) $^ -o $@ $(EXTRA)

$(BIN_DIR):
	mkdir -p $@

$(UTILS_OBJ_DIR)/%.o:: $(UTILS_DIR)/%.cpp | $(UTILS_OBJ_DIR)
	g++ -std=$(CXX_VERSION) $(CFLAGS) $(LIB) -c $< -o $@ $(EXTRA)

$(UTILS_OBJ_DIR):
	mkdir -p $@

$(OBJ_DIR)/%.o:: $(SRC_DIR)/%.cpp | $(OBJ_DIR)
	g++ -std=$(CXX_VERSION) $(CFLAGS) $(LIB) -c $< -o $@ $(EXTRA)

$(OBJ_DIR):
	mkdir -p $@

clean:
	@rm -rf $(OBJ_DIR)

clean_all:
	@rm -rf $(OBJ_DIR)
