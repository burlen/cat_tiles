set -x

max_jobs=32
year=1990

for day in `seq 0 364`
do
    n_jobs=`jobs | wc -l`
    if [[ ${n_jobs} -ge ${max_jobs} ]]
    then
        wait -n
    fi
    ./cat_tiles /global/cscratch1/sd/arhoades/TGW/data/Margulis/raw/080465779517797 ${year} ${day} ./ &
done
