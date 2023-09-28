FC=gfortran
CC=gcc-13
CXX=g++-13

BUILD_DIR := build
BIN_DIR := bin

EXE := $(BIN_DIR)/main

CPP_SRC_DIR := src
CPP_INCLUDE_DIR := include
CPP_SOURCES := $(wildcard $(CPP_SRC_DIR)/*.cc)
CPP_OBJ_FILES := $(patsubst $(CPP_SRC_DIR)/%.cc,$(BUILD_DIR)/%_c.o,$(CPP_SOURCES))

# CPP_FLAGS := -Iinclude -MMD -MP -Wpacked -malign-double -mpreferred-stack-boundary=8
CPP_FLAGS := -Iinclude -MMD -MP -Wpacked -std=c++11
FORTRAN_FLAGS := -g -fPIC -fno-automatic -fno-backslash -fno-second-underscore -falign-commons
CXXFLAGS := -Wall -g
LDFLAGS := -Llib -lstdc++
LDLIBS := -lm

all: $(EXE)
.PHONY: all clean

$(EXE): $(BUILD_DIR)/libFortran $(BUILD_DIR)/libCPP $(BUILD_DIR)/main | $(BIN_DIR)
	$(CXX)  $(LDFLAGS) $(LDLIBS) $(CPP_FLAGS) $(CXXFLAGS) -lgfortran -o $@ $^
$(BIN_DIR) $(BUILD_DIR):
	mkdir -p $@

$(BUILD_DIR)/main: main.cc | $(BUILD_DIR)
	$(CXX) $(CPP_FLAGS) $(CXXFLAGS) -c $< -o $@ $(LDFLAGS) $(LDLIBS)

$(BUILD_DIR)/libFortran: $(CPP_SRC_DIR)/radial.f | $(BUILD_DIR)
	$(FC) $(FORTRAN_FLAGS) -c $< -o $@ -J$(BUILD_DIR)

$(BUILD_DIR)/libCPP: $(CPP_OBJ_FILES) | $(BUILD_DIR)
	ld -l $^ -o $@

$(BUILD_DIR)/%_c.o: $(CPP_SRC_DIR)/%.cc | $(BUILD_DIR)
	$(CXX) $(CPP_FLAGS) $(CXXFLAGS) -c $< -o $@  $(LDFLAGS) $(LDLIBS)

clean:
	@$(RM) -rv $(BIN_DIR) $(BUILD_DIR)

-include $(OBJ:_c.o=.d)