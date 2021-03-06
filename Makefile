# Makefile to compile HMx

# Standard HMx flags
HMX_FFLAGS = \
	-std=gnu \
	-fcheck=all,no-array-temps \
	-fmax-errors=4 \
	-ffpe-trap=invalid,zero,overflow \
	-fimplicit-none \
	-O3 \
	-fdefault-real-8 \
	-fdefault-double-8 \
	-fopenmp #\
	-Warray-bounds #\
	-lgfortran \
	-lm

# Extra debugging flags
DEBUG_FLAGS = \
	-Wall \
	-fcheck=all \
	-fbounds-check \
	-fbacktrace \
	-Og

# Extra profiling flags
#PROFILE_FLAGS = \
#	-lprofiler

ifeq ($(COSMOSIS_SRC_DIR),)
# No cosmosis
FC = gfortran
FFLAGS = $(HMX_FFLAGS) -ffree-line-length-none 
all: bin lib
else
# With cosmosis
include $(COSMOSIS_SRC_DIR)/config/compilers.mk
COSMOSIS_FFLAGS := $(FFLAGS)
USER_LDFLAGS = -lcosmosis_fortran
FFLAGS = $(HMX_FFLAGS) $(COSMOSIS_FFLAGS)
all: bin lib cosmosis
endif

# Source-code directory
SRC_DIR = src

# Build directory
BUILD_DIR = build

# Module directory
MOD_DIR = /Users/Mead/Physics/library/src

# Debug build directory
DEBUG_BUILD_DIR = debug_build

# Profile build directory
PROFILE_BUILD_DIR = profile_build

# Library directory
LIB_DIR = lib

# Executable directory
BIN_DIR = bin

# Tests directory
TEST_DIR = tests

# Objects
_OBJ = \
	precision.o \
	constants.o \
	physics.o \
	sorting.o \
	special_functions.o \
	basic_operations.o \
	array_operations.o \
	file_info.o \
	io.o \
	random_numbers.o \
	table_integer.o \
	interpolate.o \
	solve_equations.o \
	string_operations.o \
	calculus_table.o \
	statistics.o \
	calculus.o \
	minimization.o \
	camb_stuff.o \
	cosmology_functions.o \
	hmx.o \
	limber.o \
	cosmic_emu_stuff.o \
	owls_stuff.o \
	owls_extras.o \
	multidark_stuff.o

# Add prefixes of build directory to objects
OBJ = $(addprefix $(BUILD_DIR)/,$(_OBJ))
DEBUG_OBJ = $(addprefix $(DEBUG_BUILD_DIR)/,$(_OBJ))
PROFILE_OBJ = $(addprefix $(PROFILE_BUILD_DIR)/,$(_OBJ))

# Make directories if they do not exist
make_dirs = @mkdir -p $(@D)

# Standard HMx
bin: $(BIN_DIR)/HMx

# TILMAN: Library
lib: $(LIB_DIR)/libhmx.a

# TILMAN: Cosmosis interface
cosmosis: lib $(LIB_DIR)/HMx_cosmosis_interface.so

# TILMAN: Test
test: $(TEST_DIR)/test_gas_gas

# HMx debugging rules
debug: FFLAGS += $(DEBUG_FLAGS)
debug: $(BIN_DIR)/HMx_debug

# HMx profile rules
#profile: FFLAGS += $(PROFILE_FLAGS)
profile: $(BIN_DIR)/HMx_profile

# Fitting
fitting: $(BIN_DIR)/HMx_fitting

# Fitting debugging
fitting_debug: FFLAGS += $(DEBUG_FLAGS)
fitting_debug: $(BIN_DIR)/HMx_fitting_debug

# Rule to make object files
$(BUILD_DIR)/%.o: $(MOD_DIR)/%.f90
	$(make_dirs)
	$(FC) -c -o $@ $< -J$(BUILD_DIR) $(LDFLAGS) $(FFLAGS)

# Rule to make HMx executable
$(BIN_DIR)/HMx: $(OBJ) $(SRC_DIR)/HMx_driver.f90
	@echo "\nBuilding executable.\n"
	$(make_dirs)
	$(FC) -o $@ $^ -J$(BUILD_DIR) $(LDFLAGS) $(FFLAGS)

# Rules to make test executables
$(TEST_DIR)/test_gas_gas: $(OBJ) $(TEST_DIR)/test_gas_gas.f90
	@echo "\nBuilding tests.\n"
	$(FC) -o $@ $^ -J$(BUILD_DIR) $(LDFLAGS) $(FFLAGS)

# Rule to make debugging objects
$(DEBUG_BUILD_DIR)/%.o: $(MOD_DIR)/%.f90
	$(make_dirs)
	$(FC) -c -o $@ $< -J$(DEBUG_BUILD_DIR) $(LDFLAGS) $(FFLAGS)

# Rule to make debugging executable
$(BIN_DIR)/HMx_debug: $(DEBUG_OBJ) $(SRC_DIR)/HMx_driver.f90
	@echo "\nBuilding debugging executable.\n"
	$(FC) -o $@ $^ -J$(DEBUG_BUILD_DIR) $(LDFLAGS) $(FFLAGS)

# Rule to make profiling objects
$(PROFILE_BUILD_DIR)/%.o: $(MOD_DIR)/%.f90
	$(make_dirs)
	$(FC) -c -o $@ $< -J$(PROFILE_BUILD_DIR) $(LDFLAGS) -lprofiler $(FFLAGS)

# Rule to make profile executable
$(BIN_DIR)/HMx_profile: $(PROFILE_OBJ) $(SRC_DIR)/HMx_driver.f90
	@echo "\nBuilding profiling executable.\n"
	$(FC) -o $@ $^ -J$(PROFILE_BUILD_DIR) $(LDFLAGS) -lprofiler $(FFLAGS)

# Rules to make fitting executables
$(BIN_DIR)/HMx_fitting: $(OBJ) $(SRC_DIR)/HMx_fitting.f90
	@echo "\nBuilding fitting.\n"
	$(make_dirs)
	$(FC) -o $@ $^ -J$(BUILD_DIR) $(LDFLAGS) $(FFLAGS)

# Rule to make fitting debugging executable
$(BIN_DIR)/HMx_fitting_debug: $(DEBUG_OBJ) $(SRC_DIR)/HMx_fitting.f90
	@echo "\nBuilding fitting debugging executable.\n"
	$(FC) -o $@ $^ -J$(DEBUG_BUILD_DIR) $(LDFLAGS) $(FFLAGS)

# Rule to make HMx static library
$(LIB_DIR)/libhmx.a: $(OBJ)
	@echo "\nBuilding static library.\n"
	$(make_dirs)
	$(AR) rc $@ $^

# Rule to make cosmosis interface
$(LIB_DIR)/HMx_cosmosis_interface.so: $(SRC_DIR)/cosmosis_interface.f90
	@echo "\nBuilding cosmosis interface.\n"
	$(FC) $(FFLAGS) -shared -o $@ $^ -L$(LIB_DIR) -lhmx $(LDFLAGS) -lcosmosis -I$(BUILD_DIR) -J$(BUILD_DIR)

# Clean up
.PHONY: clean
clean:
	rm -f $(BIN_DIR)/HMx
	rm -f $(BIN_DIR)/HMx_debug
	rm -f $(BIN_DIR)/HMx_profile
	rm -f $(BIN_DIR)/HMx_fitting
	rm -f $(BIN_DIR)/HMx_fitting_debug
	rm -f $(LIB_DIR)/libhmx.a
	rm -f $(LIB_DIR)/HMx_cosmosis_interface.so
	rm -f $(BUILD_DIR)/*.o
	rm -f $(BUILD_DIR)/*.mod
	rm -f $(DEBUG_BUILD_DIR)/*.o
	rm -f $(DEBUG_BUILD_DIR)/*.mod
	rm -f $(PROFILE_BUILD_DIR)/*.o
	rm -f $(PROFILE_BUILD_DIR)/*.mod
	test -n "$(LIB_DIR)" && rm -rf $(LIB_DIR)/HMx_cosmosis_interface.so.dSYM/
	test -n "$(BIN_DIR)" && rm -rf $(BIN_DIR)/HMx.dSYM/
	test -n "$(BIN_DIR)" && rm -rf $(BIN_DIR)/HMx_debug.dSYM/
