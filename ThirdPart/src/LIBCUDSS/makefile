# cuDSS
lib=$(TRUST_LIB)/src/LIBCUDSS/lib/libcudss_static.a

all: $(lib) 
$(lib): install.sh makefile
	@make clean
	./install.sh

clean:
	rm -rf $(TRUST_LIB)/src/LIBCUDSS
