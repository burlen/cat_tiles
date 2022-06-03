#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <netcdf.h>

/******************************************************************
 * This utility was written by Burlen Loring <bloring@lbl.gov> and
 * Alan Rhoades <arhoades@lbl.gov> for the CASCADE project. The
 * authors retain all copyrights while making it freely available
 * under the terms of a BSD licence with the limitations required
 * by the US DOE and LBNL.
 ******************************************************************/

#define USE_MEMCPY

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
            "usage:\n\n    cat_tiles [input dir] [w_year_start] [day] [output dir]\n\n"
            "note: days are numbered from zero and go to n_days - 1.");
        return 1;
    }

    int verbose = 0;
    if (getenv("VERBOSE"))
    {
        verbose = 1;
    }

    const char *input_dir = argv[1];
    int w_year_start = atoi(argv[2]);
    int w_year_end = w_year_start - 1900 + 1;
    size_t day = atoi(argv[3]);
    const char *output_dir = argv[4];

    int ierr = 0;

    // hard coded values from the dataset documentation
    // https://nsidc.org/data/WUS_UCLA_SR/versions/1
    float bounds[] = {-125, -102, 31, 49};
    int n_tiles_x = 23;
    int n_tiles_y = 18;
    size_t tile_size = 225;
    float dx = 1. / tile_size;
    size_t nx =  tile_size*n_tiles_x;
    size_t ny =  tile_size*n_tiles_y;
    size_t ns = 5;
    size_t nd = 1; //365;
    float fill_val = -999.0;

    // generate coordinate arrays
    // lon
    size_t ni = nx > ny ? nx  : ny;
    int *idx = (int*)malloc(nx*sizeof(int));
    for (size_t i = 0; i < ni; ++i)
    {
        idx[i] = i;
    }

    float *lon = (float*)malloc(nx*sizeof(float));
    for (size_t i = 0; i < nx; ++i)
    {
        lon[i] = bounds[0] + idx[i] * dx;
    }

    // lat
    float *lat = (float*)malloc(ny*sizeof(float));
    for (size_t i = 0; i < ny; ++i)
    {
        lat[i] = bounds[2] + (int(ny) - 1 - idx[i]) * dx; // descending
    }

    // stats
    int stats[5] = {0,1,2,3,4};

    // time
    float *time = (float*)malloc(nd*sizeof(float));
    for (size_t i = 0; i < nd; ++i)
    {
        time[i] = day + idx[i];
    }

    // read the tiles into a buffer, copy into a larger buffer and write the
    // file in a single call
    size_t nxy = nx*ny;
    size_t nxys = nx*ny*ns;
    size_t out_size = nxys*nd;
    float *swep_out = (float*)malloc(out_size*sizeof(float));
    float *scap_out = (float*)malloc(out_size*sizeof(float));
    float *sdp_out = (float*)malloc(out_size*sizeof(float));

    size_t ta = tile_size*tile_size;
    size_t tas = ta*ns;
    size_t in_size = tas*nd;
    float *swep_in = (float*)malloc(in_size*sizeof(float));
    float *scap_in = (float*)malloc(in_size*sizeof(float));
    float *sdp_in = (float*)malloc(in_size*sizeof(float));

    // input tile data for nc get var
    size_t istart[4] = {day, 0, 0, 0};
    size_t icount[4] = {nd, ns, tile_size, tile_size};

    int n_tiles_swep = 0;
    int n_tiles_sdp = 0;

    for (int j = 0; j < n_tiles_y; ++j)
    {
        for (int i = 0; i < n_tiles_x; ++i)
        {
            int deg_w = -bounds[0] - i;
            int deg_n =  bounds[2] + j;

            // open SWE/SCA tile
            char in_file[512] = {'\0'};
            sprintf(in_file, "%s/WUS_UCLA_SR_v01_N%02d_0W%03d_0_agg_16_WY%04d_%02d_SWE_SCA_POST.nc",
                input_dir, deg_n, deg_w, w_year_start, w_year_end);

            int ifid = 0;
            if ((ierr = nc_open(in_file, NC_NOWRITE, &ifid)))
            {
                // some tiles are missing, try to catch other errors
                if (ierr != ENOENT)
                {
                    fprintf(stderr, "ERROR: Failed to open %s\n", in_file);
                    check_error(ierr, "nc_open", __LINE__);
                }

                if (verbose)
                {
                    fprintf(stderr, "WARNING: Missing SWE/SCA tile."
                        " WY%04d_%02d  D%03zu  N%d_0  W%d_0\n",
                        w_year_start, w_year_end, day, deg_n, deg_w);
                }

                 for (size_t i = 0; i < in_size; ++i)
                 {
                     swep_in[i] = fill_val;
                 }

                 for (size_t i = 0; i < in_size; ++i)
                 {
                     scap_in[i] = fill_val;
                 }
            }
            else
            {
                // read swe
                int swep_id_in = 0;
                ierr = nc_inq_varid(ifid, "SWE_Post", &swep_id_in);
                check_error(ierr, "nc_inq_varid", __LINE__);

                ierr = nc_get_vara_float(ifid, swep_id_in, istart, icount, swep_in);
                check_error(ierr, "nc_get_var_float", __LINE__);

                // read sca
                int scap_id_in = 0;
                ierr = nc_inq_varid(ifid, "SCA_Post", &scap_id_in);
                check_error(ierr, "nc_inq_varid", __LINE__);

                ierr = nc_get_vara_float(ifid, scap_id_in, istart, icount, scap_in);
                check_error(ierr, "nc_get_var_float", __LINE__);

                // close the tile
                ierr = nc_close(ifid);
                check_error(ierr, "nc_close", __LINE__);

                if (verbose)
                {
                    fprintf(stderr, "STATUS: Copied SWE/SCA tile."
                        " WY%04d_%02d  D%03zu  N%d_0  W%d_0\n",
                        w_year_start, w_year_end, day, deg_n, deg_w);
                }

                ++n_tiles_swep;
            }

            // open SD tile
            sprintf(in_file, "%s/WUS_UCLA_SR_v01_N%02d_0W%03d_0_agg_16_WY%04d_%02d_SD_POST.nc",
                input_dir, deg_n, deg_w, w_year_start, w_year_end);

            if ((ierr = nc_open(in_file, NC_NOWRITE, &ifid)))
            {
                // some tiles are missing, try to catch other errors
                if (ierr != ENOENT)
                {
                    fprintf(stderr, "ERROR: Failed to open %s\n", in_file);
                    check_error(ierr, "nc_open", __LINE__);
                }

                if (verbose)
                {
                    fprintf(stderr, "WARNING: Missing SD tile."
                        " WY%04d_%02d  D%03zu  N%d_0  W%d_0\n",
                        w_year_start, w_year_end, day, deg_n, deg_w);
                }

                 for (size_t i = 0; i < in_size; ++i)
                 {
                     sdp_in[i] = fill_val;
                 }
            }
            else
            {
                // read sd
                int sdp_id_in = 0;
                ierr = nc_inq_varid(ifid, "SD_Post", &sdp_id_in);
                check_error(ierr, "nc_inq_varid", __LINE__);

                ierr = nc_get_vara_float(ifid, sdp_id_in, istart, icount, sdp_in);
                check_error(ierr, "nc_get_var_float", __LINE__);

                // close the tile
                ierr = nc_close(ifid);
                check_error(ierr, "nc_close", __LINE__);

                if (verbose)
                {
                    fprintf(stderr, "STATUS: Copied SD tile."
                        " WY%04d_%02d  D%03zu  N%d_0  W%d_0\n",
                        w_year_start, w_year_end, day, deg_n, deg_w);
                }

                ++n_tiles_sdp;
            }

            // copy into the output buffer
            size_t ostart[4] = {0, 0, i*tile_size, (n_tiles_y - 1 - j)*tile_size};
            size_t ocount[4] = {nd, ns, tile_size, tile_size};

#ifdef DEBUG
            fprintf(stderr, "\ttile loc = [%zu, %zu] \t corner = [%gN, %gW]\n",
                ostart[2] / tile_size, ostart[3] / tile_size, lat[ostart[3]], lon[ostart[2]]);
#endif

            for (size_t p = 0; p < ocount[0]; ++p)
            {
                size_t pp = p + ostart[0];
                for (size_t o = 0; o < ocount[1]; ++o)
                {
                    size_t oo = o + ostart[1];
                    for (size_t n = 0; n < ocount[2]; ++n)
                    {
                        size_t nn = n + ostart[2];
#ifdef USE_MEMCPY
                        memcpy(swep_out + pp*nxys + oo*nxy + nn*ny + ostart[3],
                            swep_in + p*tas + o*ta + n*tile_size, ocount[3]*sizeof(float));

                        memcpy(scap_out + pp*nxys + oo*nxy + nn*ny + ostart[3],
                            scap_in + p*tas + o*ta + n*tile_size, ocount[3]*sizeof(float));

                        memcpy(sdp_out + pp*nxys + oo*nxy + nn*ny + ostart[3],
                            sdp_in + p*tas + o*ta + n*tile_size, ocount[3]*sizeof(float));
#else
                        for (size_t m = 0; m < ocount[3]; ++m)
                        {
                            size_t mm = m + ostart[3];

                            swep_out[pp*nxys + oo*nxy + nn*ny + mm]
                                = swep_in[p*tas + o*ta + n*tile_size + m];

                            scap_out[pp*nxys + oo*nxy + nn*ny + mm]
                                = scap_in[p*tas + o*ta + n*tile_size + m];

                            sdp_out[pp*nxys + oo*nxy + nn*ny + mm]
                                = sdp_in[p*tas + o*ta + n*tile_size + m];
                        }
#endif
                    }
                }
            }
        }
    }

    // check that we did something
    if (n_tiles_swep == 0)
    {
        fprintf(stderr, "ERROR: Failed to read any SWE_SCA_POST tiles!"
            "  WY%04d_%02d D%03zu\n", w_year_start, w_year_end, day);
        abort();
    }

    if (n_tiles_sdp == 0)
    {
        fprintf(stderr, "ERROR: Failed to read any SD_POST tiles!"
            "  WY%04d_%02d D%03zu\n", w_year_start, w_year_end, day);
    }

    if (n_tiles_swep != n_tiles_sdp)
    {
        fprintf(stderr, "ERROR: The number of SWE_SCA_POST tiles(%d) differs"
            " from the number of SD_POST tiles(%d).  WY%04d_%02d D%03zu\n",
            n_tiles_swep, n_tiles_sdp, w_year_start, w_year_end, day);
    }

    if (verbose)
    {
        fprintf(stderr, "STATUS: Writing concatenated data"
            " for WY%04d-%02d D%03zu.\n", w_year_start, w_year_end, day);
    }

    // create the output file
    char out_file[256] = {'\0'};
    sprintf(out_file, "%s/WUS_UCLA_SR_v01_ALL_0_agg_16_WY%04d_%02d_D%03zu_SD_SWE_SCA_POST.nc",
        output_dir, w_year_start, w_year_end, day);

    int ofid = 0;
    ierr = nc_create(out_file, NC_CLOBBER|NC_NETCDF4, &ofid);
    check_error(ierr, "nc_create", __LINE__);

    // turn off fill mode
    int old_fill = 0;
    ierr = nc_set_fill(ofid, NC_NOFILL, &old_fill);
    check_error(ierr, "nc_set_fill", __LINE__);

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

    int time_did = 0;
    ierr = nc_def_dim(ofid, "time", nd, &time_did);
    check_error(ierr, "nc_def_dim", __LINE__);

    // define output variables
    // swep
    int swep_id = 0;
    int var_dims[4] = {time_did, stat_did, lon_did, lat_did};
    ierr = nc_def_var(ofid, "SWE_Post", NC_FLOAT, 4, var_dims, &swep_id);
    check_error(ierr, "nc_def_var", __LINE__);

    ierr = nc_def_var_fill(ofid, swep_id, NC_FILL, &fill_val);
    check_error(ierr, "nc_def_var_fill", __LINE__);

    ierr = nc_put_att(ofid, swep_id, "_FillValue", NC_FLOAT, 1, &fill_val);
    check_error(ierr, "nc_put_att", __LINE__);

    ierr = nc_put_att_text(ofid, swep_id, "units", strlen("meters"), "meters");
    check_error(ierr, "nc_put_att_text", __LINE__);

    // scap
    int scap_id = 0;
    ierr = nc_def_var(ofid, "SCA_Post", NC_FLOAT, 4, var_dims, &scap_id);
    check_error(ierr, "nc_def_var", __LINE__);

    ierr = nc_def_var_fill(ofid, scap_id, NC_FILL, &fill_val);
    check_error(ierr, "nc_def_var_fill", __LINE__);

    ierr = nc_put_att(ofid, scap_id, "_FillValue", NC_FLOAT, 1, &fill_val);
    check_error(ierr, "nc_put_att", __LINE__);

    // sdp
    int sdp_id = 0;
    ierr = nc_def_var(ofid, "SD_Post", NC_FLOAT, 4, var_dims, &sdp_id);
    check_error(ierr, "nc_def_var", __LINE__);

    ierr = nc_def_var_fill(ofid, sdp_id, NC_FILL, &fill_val);
    check_error(ierr, "nc_def_var_fill", __LINE__);

    ierr = nc_put_att(ofid, sdp_id, "_FillValue", NC_FLOAT, 1, &fill_val);
    check_error(ierr, "nc_put_att", __LINE__);

    ierr = nc_put_att_text(ofid, sdp_id, "units", strlen("meters"), "meters");
    check_error(ierr, "nc_put_att_text", __LINE__);

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

    // stats
    int stats_id = 0;
    ierr = nc_def_var(ofid, "Stats", NC_INT, 1, &lon_did, &stats_id);
    check_error(ierr, "nc_def_var", __LINE__);

    ierr = nc_put_att_text(ofid, stats_id, "units", strlen("degrees_east"), "degrees_east");
    check_error(ierr, "nc_put_att_text", __LINE__);

    // time
    int time_id = 0;
    ierr = nc_def_var(ofid, "time", NC_FLOAT, 1, &time_did, &time_id);
    check_error(ierr, "nc_def_var", __LINE__);

    char time_units[64] = {'\0'};
    sprintf(time_units, "days since %04d-10-01", w_year_start);
    ierr = nc_put_att_text(ofid, time_id, "units", strlen(time_units), time_units);
    check_error(ierr, "nc_put_att_text", __LINE__);

    ierr = nc_put_att_text(ofid, time_id, "calendar", strlen("standard"), "standard");
    check_error(ierr, "nc_put_att_text", __LINE__);

    // end definitions
    ierr = nc_enddef(ofid);
    check_error(ierr, "nc_enddef", __LINE__);

    // write coordinate arrays

    ierr = nc_put_var_float(ofid, lon_id, lon);
    check_error(ierr, "nc_put_var_float", __LINE__);

    ierr = nc_put_var_float(ofid, lat_id, lat);
    check_error(ierr, "nc_put_var_float", __LINE__);

    ierr = nc_put_var_int(ofid, stats_id, stats);
    check_error(ierr, "nc_put_var_int", __LINE__);


    ierr = nc_put_var_float(ofid, time_id, time);
    check_error(ierr, "nc_put_var_float", __LINE__);

    // write the buffer
    ierr = nc_put_var_float(ofid, swep_id, swep_out);
    check_error(ierr, "nc_put_vara_float", __LINE__);

    ierr = nc_put_var_float(ofid, scap_id, scap_out);
    check_error(ierr, "nc_put_vara_float", __LINE__);

    ierr = nc_put_var_float(ofid, sdp_id, sdp_out);
    check_error(ierr, "nc_put_vara_float", __LINE__);

    free(swep_out);
    free(scap_out);
    free(sdp_out);
    free(swep_in);
    free(scap_in);
    free(sdp_in);
    free(lon);
    free(lat);
    free(time);

    ierr = nc_close(ofid);
    check_error(ierr, "nc_close", __LINE__);

    fprintf(stderr, "STATUS: Finished processing WY%04d-%02d D%03zu! %d SWE/SCA tiles. %d SD tiles\n",
        w_year_start, w_year_end, day, n_tiles_swep, n_tiles_sdp);

    return 0;
}
