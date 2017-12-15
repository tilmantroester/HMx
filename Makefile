HMX_FFLAGS = -Warray-bounds -fmax-errors=4 -ffpe-trap=invalid,zero,overflow -fimplicit-none -O3 -fdefault-real-8 -fdefault-double-8
DEBUG_FLAGS = -Wall -fcheck=all -fbounds-check -fbacktrace -Og

ifeq ($(COSMOSIS_SRC_DIR),)
# No cosmosis
FC = gfortran
FFLAGS = $(HMX_FFLAGS) -std=gnu -ffree-line-length-none 
all: bin lib
else
# With cosmosis
include $(COSMOSIS_SRC_DIR)/config/compilers.mk
COSMOSIS_FFLAGS := $(FFLAGS)
FFLAGS = $(HMX_FFLAGS) $(COSMOSIS_FFLAGS)
all: bin lib cosmosis
endif

SRC_DIR = src
BUILD_DIR = build
DEBUG_BUILD_DIR = debug_build
LIB_DIR = lib
BIN_DIR = bin

_OBJ = constants.o \
	   random_numbers.o \
	   file_info.o \
	   logical_operations.o \
	   fix_polynomial.o \
	   array_operations.o \
	   table_integer.o \
	   special_functions.o \
	   interpolate.o \
	   solve_equations.o \
	   string_operations.o \
	   calculus.o \
	   calculus_table.o \
	   cosdef.o \
	   cosmology_functions.o \
	   HMx.o \
	   Limber.o

OBJ = $(addprefix $(BUILD_DIR)/,$(_OBJ))
DEBUG_OBJ = $(addprefix $(DEBUG_BUILD_DIR)/,$(_OBJ))

make_dirs = @mkdir -p $(@D)

lib: $(LIB_DIR)/libhmx.a
cosmosis: $(LIB_DIR)/HMx_cosmosis_interface.so
bin: $(BIN_DIR)/HMx

debug: FFLAGS += $(DEBUG_FLAGS)
debug: $(BIN_DIR)/HMx_debug

$(BUILD_DIR)/%.o: $(SRC_DIR)/%.f90
	$(make_dirs)
	$(FC) -c -o $@ $< -J$(BUILD_DIR) $(LDFLAGS) $(FFLAGS)

$(BIN_DIR)/HMx: $(OBJ) $(SRC_DIR)/HMx_driver.f90
	@echo "\nBuilding executable.\n"
	$(make_dirs)
	$(FC) -o $@ $^ -J$(BUILD_DIR) $(LDFLAGS) $(FFLAGS)

$(DEBUG_BUILD_DIR)/%.o: $(SRC_DIR)/%.f90
	$(make_dirs)
	$(FC) -c -o $@ $< -J$(DEBUG_BUILD_DIR) $(LDFLAGS) $(FFLAGS)

$(BIN_DIR)/HMx_debug: $(DEBUG_OBJ) $(SRC_DIR)/HMx_driver.f90
	@echo "\nBuilding debugging executable.\n"
	$(FC) -o $@ $^ -J$(DEBUG_BUILD_DIR) $(LDFLAGS) $(FFLAGS)

$(LIB_DIR)/libhmx.a: $(OBJ)
	@echo "\nBuilding static library.\n"
	$(make_dirs)
	$(AR) rc $@ $^

$(LIB_DIR)/HMx_cosmosis_interface.so: $(LIB_DIR)/libhmx.a $(SRC_DIR)/cosmosis_interface.f90
	@echo "\nBuilding cosmosis interface.\n"
	$(FC) $(FFLAGS) -shared -o $@ $^ -L$(LIB_DIR) -lhmx $(LDFLAGS) -lcosmosis_fortran -lcosmosis -I$(BUILD_DIR) -J$(BUILD_DIR)


.PHONY: clean
clean:
	rm -f $(BIN_DIR)/HMx
	rm -f $(BIN_DIR)/HMx_debug
	rm -f $(LIB_DIR)/libhmx.a
	rm -f $(LIB_DIR)/HMx_cosmosis_interface.so
	rm -f $(BUILD_DIR)/*.o
	rm -f $(BUILD_DIR)/*.mod
	rm -f $(SRC_DIR)/*.mod
	rm -f $(DEBUG_BUILD_DIR)/*.o
	rm -f $(DEBUG_BUILD_DIR)/*.mod
	test -n "$(LIB_DIR)" && rm -rf $(LIB_DIR)/HMx_cosmosis_interface.so.dSYM/
	test -n "$(BIN_DIR)" && rm -rf $(BIN_DIR)/HMx.dSYM/
	test -n "$(BIN_DIR)" && rm -rf $(BIN_DIR)/HMx_debug.dSYM/
