#HMx makefile

#MEAD: Added my favourite flags
FC = gfortran
FFLAGS = -std=gnu -Warray-bounds -ffree-line-length-none -fmax-errors=4 -ffpe-trap=invalid,zero,overflow -fimplicit-none -O3
DEBUGFLAGS = -Wall -fcheck=all -fbounds-check -fbacktrace -Og
#Mead: Done

SRC_DIR=src/
BUILD_DIR=build/

all: $(BUILD_DIR)HMx

#MEAD: Added debug options
debug: FFLAGS += $(DEBUGFLAGS)
debug: $(BUILD_DIR)HMx
#MEAD: Done

#MEAD: Changed compile command to add FFLAGS
$(BUILD_DIR)HMx: $(SRC_DIR)HMx.f90
	$(FC) $(FFLAGS) -o $@ $< -J$(BUILD_DIR)
#MEAD: Done

#$(FC) -std=gnu -ffree-line-length-none -o $@ $< -J$(BUILD_DIR)

clean:
	rm -f $(BUILD_DIR)HMx
	rm -f $(BUILD_DIR)hmx.o
	rm -f $(BUILD_DIR)cosdef.mod
	rm -f $(BUILD_DIR)HMx.mod
