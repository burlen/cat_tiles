#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <netcdf.h>

void check_error(int ierr, const char *call, int line)
{
    if (ierr != NC_NOERR)
    {
        fprintf(stderr, "\n\nERROR: In %s at line %d.\n %s\n\n", call, line, nc_strerror(ierr));
        abort();
    }
}

int main(int argc, char **argv)
{
    if (argc != 5)
    {
        fprintf(stderr, "ERROR: wrong number of arguments.\n"
            "usage:\n\n    cat_tiles [input dir] [year] [day] [output dir]\n\n"
            "note: days are numbered from zero and go to n_days - 1.");
        return 1;
    }

    const char *input_dir = argv[1];
    int year = atoi(argv[2]);
    size_t day = atoi(argv[3]);
    const char *output_dir = argv[4];

    fprintf(stderr, "Processing year %04d day %03zu\n", year, day);

    int ierr = 0;

    // hard coded values from the dataset documentation
    // https://nsidc.org/data/WUS_UCLA_SR/versions/1
    double bounds[] = {-125, -102, 31, 49};
    int n_tiles_x = 23;
    int n_tiles_y = 18;
    size_t tile_size = 225;
    double dx = 1. / tile_size;
    size_t nx =  tile_size*n_tiles_x;
    size_t ny =  tile_size*n_tiles_y;
    size_t ns = 5;
    size_t nd = 1; //365;
    float fill_val = -999.0;

    // create the output file
    char out_file[256] = {'\0'};
    sprintf(out_file, "%s/WUS_UCLA_SR_v01_ALL_0_agg_16_WY%04d_D%03zu_91_SWE_SCA_POST.nc",
        output_dir, year, day);

    int ofid = 0;
    ierr = nc_create(out_file, NC_CLOBBER|NC_NETCDF4, &ofid);
    check_error(ierr, "nc_create", __LINE__);

    // define dimensions
    int lon_did = 0;
    ierr = nc_def_dim(ofid, "Longitude", nx, &lon_did);
    check_error(ierr, "nc_def_dim", __LINE__);

    int lat_did = 0;
    ierr = nc_def_dim(ofid, "Latitude", ny, &lat_did);
    check_error(ierr, "nc_def_dim", __LINE__);

    int stat_did = 0;
    ierr = nc_def_dim(ofid, "Stats", ns, &stat_did);
    check_error(ierr, "nc_def_dim", __LINE__);

    int day_did = 0;
    ierr = nc_def_dim(ofid, "Day", nd, &day_did);
    check_error(ierr, "nc_def_dim", __LINE__);

    // define output variables
    // swep
    int swep_id = 0;
    int var_dims[4] = {day_did, stat_did, lon_did, lat_did};
    ierr = nc_def_var(ofid, "SWE_Post", NC_FLOAT, 4, var_dims, &swep_id);
    check_error(ierr, "nc_def_var", __LINE__);

    ierr = nc_def_var_fill(ofid, swep_id, NC_FILL, &fill_val);
    check_error(ierr, "nc_def_var_fill", __LINE__);

    ierr = nc_put_att(ofid, swep_id, "_FillValue", NC_FLOAT, 1, &fill_val);
    check_error(ierr, "nc_put_att", __LINE__);

    ierr = nc_put_att_text(ofid, swep_id, "units", strlen("meters"), "meters");
    check_error(ierr, "nc_put_att_text", __LINE__);

    // swac
    int scap_id = 0;
    ierr = nc_def_var(ofid, "SCA_Post", NC_FLOAT, 4, var_dims, &scap_id);
    check_error(ierr, "nc_def_var", __LINE__);

    ierr = nc_def_var_fill(ofid, scap_id, NC_FILL, &fill_val);
    check_error(ierr, "nc_def_var_fill", __LINE__);

    ierr = nc_put_att(ofid, scap_id, "_FillValue", NC_FLOAT, 1, &fill_val);
    check_error(ierr, "nc_put_att", __LINE__);

    // lat
    int lat_id = 0;
    ierr = nc_def_var(ofid, "Latitude", NC_FLOAT, 1, &lat_did, &lat_id);
    check_error(ierr, "nc_def_var", __LINE__);

    ierr = nc_put_att_text(ofid, lat_id, "units", strlen("degrees_north"), "degrees_north");
    check_error(ierr, "nc_put_att_text", __LINE__);

    // lon
    int lon_id = 0;
    ierr = nc_def_var(ofid, "Longitude", NC_FLOAT, 1, &lon_did, &lon_id);
    check_error(ierr, "nc_def_var", __LINE__);

    ierr = nc_put_att_text(ofid, lon_id, "units", strlen("degrees_east"), "degrees_east");
    check_error(ierr, "nc_put_att_text", __LINE__);

    // end definitions
    ierr = nc_enddef(ofid);
    check_error(ierr, "nc_enddef", __LINE__);

    // write coordinate arrays
    float *x = (float*)malloc(nx*sizeof(float));
    for (size_t i = 0; i < nx; ++i)
    {
        x[i] = bounds[0] + i * dx;
    }

    ierr = nc_put_var_float(ofid, lon_id, x);
    check_error(ierr, "nc_put_var_float", __LINE__);

    float *y = (float*)malloc(ny*sizeof(float));
    for (size_t i = 0; i < ny; ++i)
    {
        y[i] = bounds[2] + (ny - 1 - i) * dx; // descending
    }

    ierr = nc_put_var_float(ofid, lat_id, y);
    check_error(ierr, "nc_put_var_float", __LINE__);

    // read and then write tiles
    size_t in_size = tile_size*tile_size*ns*nd;
    float *swep = (float*)malloc(in_size*sizeof(float));
    float *scap = (float*)malloc(in_size*sizeof(float));

    for (int j = 0; j < n_tiles_y; ++j)
    {
        for (int i = 0; i < n_tiles_x; ++i)
        {
            int deg_w = -bounds[0] - i;
            int deg_n =  bounds[2] + j;

            // open tile
            char in_file[512] = {'\0'};
            sprintf(in_file, "%s/WUS_UCLA_SR_v01_N%02d_0W%03d_0_agg_16_WY%04d_91_SWE_SCA_POST.nc",
                input_dir, deg_n, deg_w, year);

            int ifid = 0;
            if ((ierr = nc_open(in_file, NC_NOWRITE, &ifid)))
            {
                // some tiles are missing
                fprintf(stderr, "WARNING: Missing tile %dN %dW. Initializing with _FillValue\n",
                    deg_n, deg_w);
                continue;
            }

            // read swe
            size_t istart[4] = {day, 0, 0, 0};
            size_t icount[4] = {nd, ns, tile_size, tile_size};

            int swep_id_in = 0;
            ierr = nc_inq_varid(ifid, "SWE_Post", &swep_id_in);
            check_error(ierr, "nc_inq_varid", __LINE__);

            ierr = nc_get_vara_float(ifid, swep_id_in, istart, icount, swep);
            check_error(ierr, "nc_get_var_float", __LINE__);

            // read sca
            int scap_id_in = 0;
            ierr = nc_inq_varid(ifid, "SCA_Post", &scap_id_in);
            check_error(ierr, "nc_inq_varid", __LINE__);

            ierr = nc_get_vara_float(ifid, scap_id_in, istart, icount, scap);
            check_error(ierr, "nc_get_var_float", __LINE__);

            // close the tile
            ierr = nc_close(ifid);
            check_error(ierr, "nc_close", __LINE__);

            // write the tile
            size_t ostart[4] = {0, 0, i*tile_size, (n_tiles_y - 1 - j)*tile_size};
            size_t ocount[4] = {nd, ns, tile_size, tile_size};

            ierr = nc_put_vara_float(ofid, swep_id,ostart, ocount, swep);
            check_error(ierr, "nc_put_vara_float", __LINE__);

            ierr = nc_put_vara_float(ofid, scap_id, ostart, ocount, scap);
            check_error(ierr, "nc_put_vara_float", __LINE__);

            fprintf(stderr, "STATUS: copied tile WY%04d_D%03zu_N%d_0_W%d_0"
                " \t tile loc = [%zu, %zu] \t corner = [%gN, %gW]\n",
                year, day, deg_n, deg_w, ostart[2] / tile_size, ostart[3] / tile_size,
                y[ostart[3]], x[ostart[2]]);
        }
    }

    free(swep);
    free(scap);
    free(x);
    free(y);

    ierr = nc_close(ofid);
    check_error(ierr, "nc_close", __LINE__);

    return 0;
}
