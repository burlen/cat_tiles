# cat_tiles
This is a code to concatenate tiles of data in netcdf format from the Western United States UCLA Daily Snow Reanalysis, Version 1 dataset. https://nsidc.org/data/WUS_UCLA_SR/versions/1

## Build on Cori
```
module swap PrgEnv-intel PrgEnv-gnu
module load cray-netcdf
make
```

## Usage
```
cat_tiles [input dir] [year] [day] [output dir]
```
See the script process_year.sh for an example.

## License
This utility was written by Burlen Loring <bloring@lbl.gov> and
Alan Rhoades <arhoades@lbl.gov> for the CASCADE project. The
authors retain all copyrights while making it freely available
under the terms of a BSD licence with the limitations required
by the US DOE and LBNL.
