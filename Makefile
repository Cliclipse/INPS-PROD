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
# Toolchain configuration
# ----------------------------------------------------------------------------
CXX      := g++
CXXSTD   := -std=c++20
WARNINGS := -Wall -Wextra -Wpedantic
OPTIMIZE := -O2
CXXFLAGS := $(CXXSTD) $(WARNINGS) $(OPTIMIZE)

INCLUDES := -I. \
            -I"armadillo-code-15.0.x/armadillo-code-15.0.x/include"

LDLIBS   := -larmadillo

# ----------------------------------------------------------------------------
# Project layout
# ----------------------------------------------------------------------------
TARGET  := main
OBJ_DIR := build

SRCS := \
	Basis.cpp \
	Poly.cpp \
	Hermit.cpp \
	OneDHOSolution.cpp \
	rho.cpp \
	z_minus_r_minus.cpp \
	Main.cpp

OBJS := $(SRCS:%.cpp=$(OBJ_DIR)/%.o)

# ----------------------------------------------------------------------------
# Primary targets
# ----------------------------------------------------------------------------
.PHONY: all clean

all: $(TARGET)

$(TARGET): $(OBJS)
	$(CXX) $(CXXFLAGS) $^ $(LDLIBS) -o $@
	@echo "Built $(TARGET)"

$(OBJ_DIR):
	@mkdir -p "$@"

# ----------------------------------------------------------------------------
# Compilation rules
# ----------------------------------------------------------------------------
$(OBJ_DIR)/%.o: %.cpp | $(OBJ_DIR)
	@mkdir -p "$(@D)"
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c "$<" -o "$@"

# ----------------------------------------------------------------------------
# Clean-up helpers
# ----------------------------------------------------------------------------
clean:
	rm -rf "$(OBJ_DIR)" "$(TARGET)"
	@echo "Workspace cleaned"
