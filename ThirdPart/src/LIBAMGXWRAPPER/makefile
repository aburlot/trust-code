AMGX_LIB=$(TRUST_ROOT)/lib/src/LIBAMGX/AmgX/lib/libamgxsh.so
AMGXWRAPPER_DIR=$(TRUST_ROOT)/lib/src/LIBAMGX/AmgXWrapper
AMGXWRAPPER_LIB=$(AMGXWRAPPER_DIR)/lib64/libAmgXWrapper.a
PATCH_AMGXWRAPPER := $(wildcard patch_AmgXWrapper/src/*.* patch_AmgXWrapper/example/*/*)
mpi=$(MPI_INCLUDE)/mpi.h

all: $(AMGXWRAPPER_LIB)

$(AMGXWRAPPER_LIB): install makefile $(mpi) $(PATCH_AMGXWRAPPER) $(AMGX_LIB)
	rm -r -f $(AMGXWRAPPER_DIR)
	./install $(AMGXWRAPPER_DIR)
	
clean:
	rm -r -f $(AMGXWRAPPER_DIR)
