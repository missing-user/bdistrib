# makefile for NERSC Edison and Cori
# You must first load the cray-netcdf module:
#   module load cray-netcdf
# For Cori is is also necessary to run
#   module swap intel/16.0.0.109 intel/15.0.1.133
# to avoid a bug in the Intel MKL!!!
# It is convenient to run
#   module unload cray-libsci
# to avoid warning messages about libsci during compiling.

ifdef NERSC_HOST
        HOSTNAME = $(NERSC_HOST)
else
        HOSTNAME="laptop"
endif

ifeq ($(HOSTNAME),edison)
	FC = ftn
	## NERSC documentation recommends against specifying -O3
	## -mkl MUST APPEAR AT THE END!!
	EXTRA_COMPILE_FLAGS = -openmp -mkl
	EXTRA_LINK_FLAGS =  -openmp -mkl -Wl,-ydgemm_
	# Above, the link flag "-Wl,-ydgemm_" causes the linker to report which version of DGEMM (the BLAS3 matrix-matrix-multiplication subroutine) is used.
	# For batch systems, set the following variable to the command used to run jobs. This variable is used by 'make test'.
	BDISTRIB_COMMAND_TO_SUBMIT_JOB = srun -n 1 -c 24
else ifeq ($(HOSTNAME),cori)
	FC = ftn
	## NERSC documentation recommends against specifying -O3
	## -mkl MUST APPEAR AT THE END!!
	EXTRA_COMPILE_FLAGS = -qopenmp -mkl
	EXTRA_LINK_FLAGS =  -qopenmp -mkl -Wl,-ydgemm_
	# Above, the link flag "-Wl,-ydgemm_" causes the linker to report which version of DGEMM (the BLAS3 matrix-matrix-multiplication subroutine) is used.
	# For batch systems, set the following variable to the command used to run jobs. This variable is used by 'make test'.
	BDISTRIB_COMMAND_TO_SUBMIT_JOB = srun -n 1 -c 32
else
	FC = mpif90
	#EXTRA_COMPILE_FLAGS = -fopenmp -I/usr/include -ffree-line-length-none -cpp
	EXTRA_COMPILE_FLAGS = -fopenmp -I/usr/include -ffree-line-length-none
	EXTRA_LINK_FLAGS =  -fopenmp -L/usr/lib -lnetcdff  -lnetcdf -llapack -lblas
#-framework Accelerate

	# For batch systems, set the following variable to the command used to run jobs. This variable is used by 'make test'.
	BDISTRIB_COMMAND_TO_SUBMIT_JOB =
endif


# End of system-dependent variable assignments
LIBSTELL_DIR =/home/IPP-HGW/juph/LIBSTELL/build
LIBSTELL_LIB_BIN_NAME=libstell.a
#LIBSTELL_DIR = mini_libstell
#LIBSTELL_LIB_BIN_NAME=mini_libstell.a
TARGET = bdistrib

export

.PHONY: all clean

all: $(TARGET)

include makefile.depend

%.o: %.f90 $(LIBSTELL_DIR)/$(LIBSTELL_LIB_BIN_NAME)
	$(FC) $(EXTRA_COMPILE_FLAGS) -I $(LIBSTELL_DIR) -c $<

%.o: %.f $(LIBSTELL_DIR)/$(LIBSTELL_LIB_BIN_NAME)
	$(FC) $(EXTRA_COMPILE_FLAGS) -I $(LIBSTELL_DIR) -c $<

$(TARGET): $(OBJ_FILES) $(LIBSTELL_DIR)/$(LIBSTELL_LIB_BIN_NAME)
	$(FC) -o $(TARGET) $(OBJ_FILES) $(LIBSTELL_DIR)/$(LIBSTELL_LIB_BIN_NAME) $(EXTRA_LINK_FLAGS)
#	$(FC) -o $(TARGET) $(OBJ_FILES) $(LIBSTELL_DIR)/$(LIBSTELL_LIB_BIN_NAME) $(EXTRA_LINK_FLAGS)

$(LIBSTELL_DIR)/$(LIBSTELL_LIB_BIN_NAME):
	$(MAKE) -C $(LIBSTELL_DIR)

clean:
	rm -f *.o *.mod *.MOD *~ $(TARGET)
	cd $(LIBSTELL_DIR); rm -f *.o *.mod *.MOD *.a

test: $(TARGET)
	@echo "Beginning functional tests." && cd examples && export BDISTRIB_RETEST=no && python runExamples.py

retest: $(TARGET)
	@echo "Testing existing output files for examples without re-running then." && cd examples && export BDISTRIB_RETEST=yes && python runExamples.py
