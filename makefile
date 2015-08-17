FC = ftn
LIBSTELL_DIR = /global/homes/l/landrema/20150410-02-stellinstall_245_edison/LIBSTELL/Release
EXTRA_COMPILE_FLAGS = 
EXTRA_LINK_FLAGS = 

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
