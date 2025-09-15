# Makefile to compile and link QSI simulator programs
#
# Created on 01/08/2018
# Copyright (C) 2018 Attila Szabó <as2372@cam.ac.uk>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the Creative Commons Attribution License (CC-BY),
# version 4.0 of the License.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
#
# Please find a copy of the Creative Commons CC-BY License on
# <https://creativecommons.org/licenses/by/4.0/>.

# ---------- Compiler and linker directives ----------
# CXX = g++

CXXFLAGS = -std=c++2a -pedantic -Wall
CXXOPTS = -O1 -g
# LDFLAGS = -Wl,--gc-sections

src_path := $(dir $(abspath $(lastword $(MAKEFILE_LIST))))


# Output directory after linking
SOURCE_DIR = .
BIN_DIR = ../bin
BUILD_DIR = ../build

# ------ Establishing objects and sources ----

# Search subdirectories as well
COMMON_DIR = common
COMMON_CC = $(wildcard $(SOURCE_DIR)/$(COMMON_DIR)/*.cc )
COMMON_O = $(patsubst $(SOURCE_DIR)/%.cc,$(BUILD_DIR)/%.o,$(COMMON_CC))

MANAGER_DIR = manager
MANAGER_CC = $(wildcard $(SOURCE_DIR)/$(MANAGER_DIR)/*.cc)
MANAGER_O = $(patsubst $(SOURCE_DIR)/%.cc,$(BUILD_DIR)/%.o,$(MANAGER_CC))

STRUCT_DIR = struct
STRUCT_CC = $(wildcard $(SOURCE_DIR)/$(STRUCT_DIR)/*.cc)
STRUCT_O = $(patsubst $(SOURCE_DIR)/%.cc,$(BUILD_DIR)/%.o,$(STRUCT_CC))


# Change ~/include as necessary for your architecture
INC_PATH = -I$(COMMON_DIR) -I$(STRUCT_DIR) -I$(MANAGER_DIR) $(INC)

# ---------- Libraries needed for linking ----------

# Change or remove ~/lib as necessary for your architecture
LIB_PATH = -L/opt/homebrew/lib $(LIB)

LIBS = -lstdc++ -lm -lstdc++fs -lfftw3 -lgsl -lgslcblas

# ---------- Compilation rules ----------

# Nontrivial header dependences
DEPENDS := $(patsubst %.cc,%.d,$(SOURCES))

# Generic .cc -> .o rule: note we generally don't have .hh files
$(BUILD_DIR)/%.o: $(SOURCE_DIR)/%.cc
	mkdir -p $(@D)
	$(CXX) $(CXXFLAGS) $(CXXOPTS) $(INC_PATH) -MMD -MP $< -c -o $@

# ---------- Linking rules ----------

# Generic .o -> exec rule: uses all prerequisites (meant to be .o files)
%: $(BUILD_DIR)/%.o $(COMMON_O) $(MANAGER_O) $(STRUCT_O)
	$(CXX) $(CXXFLAGS) $(INC_PATH) $(LDFLAGS) -MMD $^ \
	$(LIB_PATH) $(LIBS) \
	-o $(BIN_DIR)/$@

# Don't allow direct compilation 
%: %.cc
# Don't allow removal cof object files in chained implicit linking
.SECONDARY:

# Cleanup:
.PHONY: clean
clean:
	-rm -r $(BUILD_DIR)/**
