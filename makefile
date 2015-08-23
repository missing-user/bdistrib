FC = ftn
LIBSTELL_DIR = /global/homes/l/landrema/20150410-02-stellinstall_245_edison/LIBSTELL/Release
EXTRA_COMPILE_FLAGS = -O3 -openmp
#EXTRA_COMPILE_FLAGS = -O0 -g -openmp
EXTRA_LINK_FLAGS =  -openmp -Wl,-ydgemm_

# Above, the link flag "-Wl,-ydgemm_" causes the linker to report which version of DGEMM (the BLAS3 matrix-matrix-multiplication subroutine) is used.

# End of system-dependent variable assignments

TARGET = bdistrib

.PHONY: all clean

all: $(TARGET)

include makefile.depend

%.o: %.f90
	$(FC) $(EXTRA_COMPILE_FLAGS) -I $(LIBSTELL_DIR) -c $<

%.o: %.f
	$(FC) $(EXTRA_COMPILE_FLAGS) -I $(LIBSTELL_DIR) -c $<

$(TARGET): $(OBJ_FILES)
	$(FC) -o $(TARGET) $(OBJ_FILES) $(LIBSTELL_DIR)/libstell.a $(EXTRA_LINK_FLAGS)

clean:
	rm -f *.o *.mod *.MOD *~ $(TARGET)

test: $(TARGET)
	@echo "Beginning functional tests." && cd examples && export BDISTRIB_RETEST=no && ./runExamples.py

retest: $(TARGET)
	@echo "Testing existing output files for examples without re-running then." && cd examples && export BDISTRIB_RETEST=yes && ./runExamples.py
