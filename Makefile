ARCH := $(shell uname -m)

BUILD_DIR := build
BIN_DIR := bin

EXE := $(BIN_DIR)/main

CPP_SRC_DIR := src
CPP_INCLUDE_DIR := include
CPP_SOURCES := $(wildcard $(CPP_SRC_DIR)/*.cc)
CPP_OBJ_FILES := $(patsubst $(CPP_SRC_DIR)/%.cc,$(BUILD_DIR)/%_c.o,$(CPP_SOURCES))

FORTRAN_FLAGS := -g -fPIC -fno-automatic -fno-backslash -fno-second-underscore -falign-commons
CPP_FLAGS_X86 := -Iinclude -MMD -MP -Wpacked -malign-double -mpreferred-stack-boundary=8
CPP_FLAGS_ARM64 := -Iinclude -MMD -MP -Wpacked -std=c++11
CXXFLAGS := -Wall -g
LD_FLAGS_X86 := -Llib -lstdc++
LD_FLAGS_ARM64 := -L/usr/local/lib -lstdc++ -ld_classic
LDLIBS := -lm

ifeq (${ARCH},arm64)
	# arch should be "arm64"
    LDFLAGS := ${LD_FLAGS_ARM64}
	CPP_FLAGS := ${CPP_FLAGS_ARM64}
	FC=gfortran
	CC=gcc-13
	CXX=g++-13
else
	# arch should be "x86_64"
    LDFLAGS := ${LD_FLAGS_X86}
	CPP_FLAGS := ${CPP_FLAGS_X86}
	FC=gfortran
	CC=gcc
	CXX=g++
endif

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
	ld -r -arch ${ARCH} $^ -o $@

$(BUILD_DIR)/%_c.o: $(CPP_SRC_DIR)/%.cc | $(BUILD_DIR)
	$(CXX) $(CPP_FLAGS) $(CXXFLAGS) -c $< -o $@  $(LDFLAGS) $(LDLIBS)

clean:
	@$(RM) -rv $(BIN_DIR) $(BUILD_DIR)

arch_test:
	@echo "Running on architecture: ${ARCH}"
	@if [ "$(ARCH)" = "arm64" ]; then \
		echo "Hello ${ARCH}"; \
	fi

-include $(OBJ:_c.o=.d)