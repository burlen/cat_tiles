
NETCDF_DIR=/opt/cray/pe/netcdf/4.8.1.1/GNU/8.2
CXX_FLAGS=-O3 -Wall -march=native -mtune=native

tmp := $(shell which icpc 2>/dev/null)
ifeq (1, $(.SHELLSTATUS))
COMP := g++
COMP_VER := $(shell gcc --version | grep ^gcc)
else
COMP := icpc
COMP_VER := $(shell icpc --version | grep ^icpc)
endif

all: cat_tiles.cpp
	@echo "Compiler : $(COMP_VER)"
	$(COMP) ${CXX_FLAGS} -I${NETCDF_DIR}/include -L${NETCDF_DIR}/lib -lnetcdf cat_tiles.cpp -o cat_tiles

clean:
	rm cat_tiles

