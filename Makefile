include ${COSMOSIS_SRC_DIR}/config/compilers.mk

SRC_DIR=src/
BUILD_DIR=build/

interface: $(BUILD_DIR)libhmx.a $(BUILD_DIR)cosmosis_interface.so
bin: $(BUILD_DIR)HMx


$(BUILD_DIR)cosmosis_interface.so: $(BUILD_DIR)libhmx.a $(SRC_DIR)cosmosis_interface.f90
	$(FC) $(FFLAGS) -shared -o $@ $+ -L$(BUILD_DIR) -lhmx $(LDFLAGS) -lcosmosis_fortran -lcosmosis -I$(BUILD_DIR) -J$(BUILD_DIR)

$(BUILD_DIR)HMx: $(SRC_DIR)HMx.f90
	$(FC) -std=gnu -ffree-line-length-none -o $@ $< -J$(BUILD_DIR)
# test: test.f90 libhmcode.a
	# $(FC) $(FFLAGS) -o $@ $< -L. -lhmcode $(LDFLAGS)
	

clean:
	rm -f $(BUILD_DIR)HMx
	rm -f $(BUILD_DIR)libhmx.a
	rm -f $(BUILD_DIR)hmx.o
	rm -f $(BUILD_DIR)cosmosis_interface.so
	rm -f $(BUILD_DIR)cosdef.mod
	rm -f $(BUILD_DIR)HMx.mod
	rm -f $(BUILD_DIR)HMx_setup.mod
	rm -rf $(BUILD_DIR)*.dSYM/

$(BUILD_DIR)libhmx.a: $(SRC_DIR)HMx.f90
	$(FC) $(FFLAGS) -c  $+ $(LDFLAGS) -o $(BUILD_DIR)hmx.o -J$(BUILD_DIR)
	$(AR) rc $@ $(BUILD_DIR)hmx.o
