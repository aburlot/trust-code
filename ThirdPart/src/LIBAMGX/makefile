AMGX_DIR=$(TRUST_ROOT)/lib/src/LIBAMGX/AmgX
AMGX_LIB=$(AMGX_DIR)/lib/libamgxsh.so
mpi=$(MPI_INCLUDE)/mpi.h

all:$(AMGX_LIB)

$(AMGX_LIB): install makefile $(mpi)
	rm -r -f $(AMGX_DIR)
	./install $(AMGX_DIR)

clean:
	rm -r -f $(AMGX_DIR)
