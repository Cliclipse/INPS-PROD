# ----------------------------------------------------------------------------
# Toolchain configuration
# ----------------------------------------------------------------------------
CXX      := g++
AR       := ar
CXXSTD   := -std=c++20
WARNINGS := -Wall -Wextra -Wpedantic
OPTIMIZE := -O2
CXXFLAGS := $(CXXSTD) $(WARNINGS) $(OPTIMIZE)

INCLUDES := -I. \
            -I"pr akito" \
            -I"armadillo-code-15.0.x/armadillo-code-15.0.x/include"

# ----------------------------------------------------------------------------
# Project layout
# ----------------------------------------------------------------------------
TARGET  := libinps.a
OBJ_DIR := build

Z_MINUS_SRC := z_minus.cpp
HERMIT_SRC  := pr\ akito/Hermit.cpp
HOS_SRC     := pr\ akito/OneDHOSolution.cpp

OBJECTS := $(OBJ_DIR)/z_minus.o \
           $(OBJ_DIR)/Hermit.o \
           $(OBJ_DIR)/OneDHOSolution.o

# ----------------------------------------------------------------------------
# Primary targets
# ----------------------------------------------------------------------------
.PHONY: all clean

all: $(TARGET)
	@echo "Built $(TARGET)"

$(TARGET): $(OBJECTS)
	$(AR) rcs $@ $^

$(OBJ_DIR):
	@mkdir -p "$@"

# ----------------------------------------------------------------------------
# Compilation rules
# ----------------------------------------------------------------------------
$(OBJ_DIR)/z_minus.o: $(Z_MINUS_SRC) z_minus.h | $(OBJ_DIR)
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c "$<" -o "$@"

$(OBJ_DIR)/Hermit.o: $(HERMIT_SRC) pr\ akito/Hermit.h | $(OBJ_DIR)
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c "$<" -o "$@"

$(OBJ_DIR)/OneDHOSolution.o: $(HOS_SRC) pr\ akito/OneDHOSolution.h | $(OBJ_DIR)
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c "$<" -o "$@"

# ----------------------------------------------------------------------------
# Clean-up helpers
# ----------------------------------------------------------------------------
clean:
	rm -rf "$(OBJ_DIR)" "$(TARGET)"
	@echo "Workspace cleaned"
