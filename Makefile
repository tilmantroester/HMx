#HMx makefile

#Set compiler and flags
FC = gfortran
FFLAGS = -std=gnu -Warray-bounds -ffree-line-length-none -fmax-errors=4 -ffpe-trap=invalid,zero,overflow -fimplicit-none -O3 -fdefault-real-8
DEBUGFLAGS = -Wall -fcheck=all -fbounds-check -fbacktrace -Og

#Directories
SRC=src
BUILD=build

#Execuatable
EXEC=$(BUILD)/HMx

#Objects to use for build
OBJECTS = $(BUILD)/constants.o $(BUILD)/random_numbers.o $(BUILD)/file_info.o $(BUILD)/fix_polynomial.o $(BUILD)/array_operations.o $(BUILD)/table_integer.o $(BUILD)/special_functions.o $(BUILD)/interpolate.o $(BUILD)/string_operations.o $(BUILD)/calculus.o $(BUILD)/calculus_table.o $(SRC)/HMx.f90

all: $(EXEC)

debug: FFLAGS += $(DEBUGFLAGS)
debug: $(EXEC)


$(EXEC): $(OBJECTS)
	$(FC) $(FFLAGS) -I$(BUILD) -L$(BUILD) $^ -o $@

$(BUILD)/constants.o: $(SRC)/constants.f90
	$(FC) $(FFLAGS) -I$(BUILD) -L$(BUILD) -c $^
	mv constants.o $(BUILD)/.
	mv constants.mod $(BUILD)/.

$(BUILD)/random_numbers.o: $(SRC)/random_numbers.f90
	$(FC) $(FFLAGS) -I$(BUILD) -L$(BUILD) -c $^
	mv random_numbers.o $(BUILD)/.
	mv random_numbers.mod $(BUILD)/.

$(BUILD)/file_info.o: $(SRC)/file_info.f90
	$(FC) $(FFLAGS) -I$(BUILD) -L$(BUILD) -c $^
	mv file_info.o $(BUILD)/.
	mv file_info.mod $(BUILD)/.

$(BUILD)/fix_polynomial.o: $(SRC)/fix_polynomial.f90
	$(FC) $(FFLAGS) -I$(BUILD) -L$(BUILD) -c $^
	mv fix_polynomial.o $(BUILD)/.
	mv fix_polynomial.mod $(BUILD)/.

$(BUILD)/array_operations.o: $(SRC)/array_operations.f90
	$(FC) $(FFLAGS) -I$(BUILD) -L$(BUILD) -c $^
	mv array_operations.o $(BUILD)/.
	mv array_operations.mod $(BUILD)/.

$(BUILD)/table_integer.o: $(SRC)/table_integer.f90
	$(FC) $(FFLAGS) -I$(BUILD) -L$(BUILD) -c $^
	mv table_integer.o $(BUILD)/.
	mv table_integer.mod $(BUILD)/.

$(BUILD)/special_functions.o: $(SRC)/special_functions.f90
	$(FC) $(FFLAGS) -I$(BUILD) -L$(BUILD) -c $^
	mv special_functions.o $(BUILD)/.
	mv special_functions.mod $(BUILD)/.

$(BUILD)/interpolate.o: $(SRC)/interpolate.f90
	$(FC) $(FFLAGS) -I$(BUILD) -L$(BUILD) -c $^
	mv interpolate.o $(BUILD)/.
	mv interpolate.mod $(BUILD)/.

$(BUILD)/string_operations.o: $(SRC)/string_operations.f90
	$(FC) $(FFLAGS) -I$(BUILD) -L$(BUILD) -c $^
	mv string_operations.o $(BUILD)/.
	mv string_operations.mod $(BUILD)/.

$(BUILD)/calculus.o: $(SRC)/calculus.f90
	$(FC) $(FFLAGS) -I$(BUILD) -L$(BUILD) -c $^
	mv calculus.o $(BUILD)/.
	mv calculus.mod $(BUILD)/.

$(BUILD)/calculus_table.o: $(SRC)/calculus_table.f90
	$(FC) $(FFLAGS) -I$(BUILD) -L$(BUILD) -c $^
	mv calculus_table.o $(BUILD)/.
	mv calculus_table.mod $(BUILD)/.

debug: FFLAGS += $(DEBUGFLAGS)
debug: $(EXEC)

clean:
	rm -f $(BUILD)/HMx
	rm -f $(BUILD)/*.mod
	rm -f $(BUILD)/*.o	

#CRAP

#$(FC) -std=gnu -ffree-line-length-none -o $@ $< -J$(BUILD)

#$(BUILD_DIR)constants.o: $(SRC_DIR)constants.f90
#	$(FC) $(FFLAGS) -c $@ $< -J$(BUILD_DIR)

#all: $(BUILD_DIR)HMx $(BUILD_DIR)constants.o $(BUILD_DIR)constants.mod

#$(BUILD_DIR)array_operations.o $(BUILD_DIR)file_info.o $(BUILD_DIR)special_functions.o $(BUILD_DIR)interpolate.o $(BUILD_DIR)string_operations.o $(BUILD_DIR)calculus.o $(BUILD_DIR)calculus_table.o
