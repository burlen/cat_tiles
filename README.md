# cat_tiles
some code to concatenate tiles of data in netcdf format


## Build on Cori
module swap PrgEnv-intel PrgEnv-gnu
module load cray-netcdf
make

## Usage
```
cat_tiles [input dir] [year] [day] [output dir]
```
See the script process_year.sh for an example.


