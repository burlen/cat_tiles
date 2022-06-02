

NETCDF_DIR=/opt/cray/pe/netcdf/4.8.1.1/GNU/8.2
CXX_FLAGS=-O3 -g -Wall


all: cat_tiles.cpp
	g++ ${CXX_FLAGS} -I${NETCDF_DIR}/include -L${NETCDF_DIR}/lib -lnetcdf cat_tiles.cpp -o cat_tiles

clean:
	rm cat_tiles

