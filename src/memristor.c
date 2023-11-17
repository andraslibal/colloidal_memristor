#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include "random_numrec.h"
#include "util.h"
#include "memristor.h"

#define PI 3.14159265358979323846

struct global_struct global;

/* IO operations an initialization */

void initialize_statistics()
{
    global.largest_cluster_size_sum = 0;
    global.avg_fx = 0;
    global.avg_fy = 0;
}

void initialize_output_files()
{
    char statistics_file_name[255];
    char movie_file_name[255];
    char histogram_file_name[255];

    // statistic file creation and opening
    strcpy(statistics_file_name, "./statistics/");

    strcat(statistics_file_name, global.baseName);
    strcat(statistics_file_name, ".txt");
    global.statistics_file = fopen(statistics_file_name, "wt");

    if (global.statistics_file == NULL)
    {
        printf("Could not open statistics file %s\n", statistics_file_name);
        exit(1);
    }
    printf("Successfully opened statistics file %s\n", statistics_file_name);

    // movie file creation and opening
    sprintf(movie_file_name, "./movies/%s.mvi", global.baseName);
    global.movie_file = fopen(movie_file_name, "wb");

    if (global.movie_file == NULL)
    {
        printf("Could not open movie file %s\n", movie_file_name);
        exit(1);
    }
    printf("Successfully opened movie file %s\n", movie_file_name);
}

void read_line_from_file(FILE *parameterfile, char descript[], char type[], void *value)
{
    double temp_double;
    int temp_int;
    uint64_t temp_uint64_t;
    char temp_descript[100];
    char temp_string[200];
    int read_items;

    if (strcmp(type, "double") == 0)
        read_items = fscanf(parameterfile, "%s %lf", temp_descript, &temp_double);
    else if (strcmp(type, "int") == 0)
        read_items = fscanf(parameterfile, "%s %d", temp_descript, &temp_int);
    else if (strcmp(type, "uint64_t") == 0)
        read_items = fscanf(parameterfile, "%s %ld", temp_descript, &temp_uint64_t);
    else if (strcmp(type, "string") == 0)
        read_items = fscanf(parameterfile, "%s %s", temp_descript, temp_string);

    if (read_items != 2)
    {
        printf("Error on parameter file line %s\n", descript);
        exit(1);
    }

    if (strcmp(temp_descript, descript) == 0)
    {

        if (strcmp(type, "double") == 0)
        {
            *(double *)value = temp_double;
            printf("Parameter : %s = %lf\n", descript, *(double *)value);
        }
        else if (strcmp(type, "int") == 0)
        {
            *(int *)value = temp_int;
            printf("Parameter : %s = %d\n", descript, *(int *)value);
        }
        else if (strcmp(type, "uint64_t") == 0)
        {
            *(uint64_t *)value = temp_uint64_t;
            printf("Parameter : %s = %ld\n", descript, *(uint64_t *)value);
        }
        else if (strcmp(type, "string") == 0)
        {
            strcpy(value, temp_string);
            // printf("%s\n",temp_string);
            printf("Parameter : %s = %s\n", descript, temp_string);
        }
    }
    else
    {
        printf("Expected %s in parameterfile\n", descript);
        exit(1);
    }
}

void initialize_from_file(int argc, char *argv[])
{
    FILE *parameterfile;
    char mvi_fname[300];

    if (argc < 2)
    {
        printf("Usage: ./memristor <parameter_file>\n");
        fflush(stdout);
        exit(1);
    }

    parameterfile = fopen(argv[1], "rt");
    if (parameterfile == NULL)
    {
        printf("Cannot open parameter file %s\n", argv[1]);
        fflush(stdout);
        exit(1);
    }
    printf("Opened parameter file %s\n", argv[1]);

    read_line_from_file(parameterfile, "basename", "string", &global.baseName);
    read_line_from_file(parameterfile, "SX", "double", &global.SX);
    read_line_from_file(parameterfile, "SY", "double", &global.SY);
    global.SX2 = global.SX / 2.0;
    global.SY2 = global.SY / 2.0;
    read_line_from_file(parameterfile, "N_particle", "int", &global.N_particle);
    read_line_from_file(parameterfile, "N_pinningsites", "int", &global.N_pinningsites);
    read_line_from_file(parameterfile, "dt", "double", &global.dt);
    read_line_from_file(parameterfile, "total_runtime", "int", &global.total_runtime);
    read_line_from_file(parameterfile, "echo_time", "int", &global.echo_time);
    read_line_from_file(parameterfile, "movie_time", "int", &global.movie_time);
    read_line_from_file(parameterfile, "statistics_time", "int", &global.statistics_time);
    read_line_from_file(parameterfile, "relax_time", "int", &global.relax_time);
    read_line_from_file(parameterfile, "relax_duration", "int", &global.relax_duration);

    read_line_from_file(parameterfile, "N_additional_particles", "int", &global.N_additional_particles);
    read_line_from_file(parameterfile, "N_steps_to_add_more_particles", "int", &global.N_steps_to_add_more_particles);
    read_line_from_file(parameterfile, "N_times_to_add_more_particles", "int", &global.N_times_to_add_more_particles);

    read_line_from_file(parameterfile, "generic_particle_R", "double", &global.generic_particle_R);
    read_line_from_file(parameterfile, "generic_particle_k_spring", "double", &global.generic_particle_k_spring);
    read_line_from_file(parameterfile, "generic_particle_motor_force", "double", &global.generic_particle_motor_force);
    read_line_from_file(parameterfile, "generic_particle_active_fraction", "double", &global.generic_particle_active_fraction);
    read_line_from_file(parameterfile, "generic_particle_motor_minimum_time", "int", &global.generic_particle_motor_minimum_time);
    global.generic_particle_motor_maximum_time = 2 * global.generic_particle_motor_minimum_time;
    printf("Calculated: generic_particle_motor_maximum_time = %d\n", global.generic_particle_motor_maximum_time);

    read_line_from_file(parameterfile, "generic_pinningsite_R", "double", &global.generic_pinningsite_R);
    read_line_from_file(parameterfile, "generic_pinningsite_fxmax", "double", &global.generic_pinningsite_fxmax);
    read_line_from_file(parameterfile, "dragging_force_x", "double", &global.dragging_force_x);
    read_line_from_file(parameterfile, "dragging_force_y", "double", &global.dragging_force_y);

    read_line_from_file(parameterfile, "eta_small", "double", &global.eta_small);
    read_line_from_file(parameterfile, "eta_large", "double", &global.eta_large);
    read_line_from_file(parameterfile, "obstacle_ratio", "double", &global.obstacle_ratio);
    read_line_from_file(parameterfile, "temperature", "double", &global.temperature);

    read_line_from_file(parameterfile, "seedToSet", "uint64_t", &global.seedToSet);
    if (global.seedToSet != 0)
        setseed(global.seedToSet);

    // Verlet list pre-init
    global.Verlet_distance_to_rebuild = global.generic_particle_R;
    printf("Calculated: Verlet_distance_to_rebuild = %lf\n", global.Verlet_distance_to_rebuild);
    global.Verlet_i = NULL;
    global.Verlet_j = NULL;
    global.N_Verlet = 0;
    global.Verlet_R = 4.0 * global.generic_particle_R;
    printf("Calculated: Verlet_R = %lf\n", global.Verlet_R);
    printf("Particle density is = %lf\n", global.N_particle * global.generic_particle_R * global.generic_particle_R * PI / global.SX / global.SY);
    printf("Base name: %s\n", global.baseName);
    printf("================================\n");
}

void write_cmovie_frame()
{
    int i;
    float floatholder;
    int intholder;

    intholder = global.N_particle + global.N_pinningsites;
    fwrite(&intholder, sizeof(int), 1, global.movie_file);

    intholder = global.time;
    fwrite(&intholder, sizeof(int), 1, global.movie_file);

    for (i = 0; i < global.N_particle; i++)
    {
        /* color */
        intholder = global.particle_color[i];
        fwrite(&intholder, sizeof(int), 1, global.movie_file);

        /* id */
        intholder = i; // ID
        fwrite(&intholder, sizeof(int), 1, global.movie_file);

        /* x, y */
        floatholder = (float)global.particle_x[i];
        fwrite(&floatholder, sizeof(float), 1, global.movie_file);
        floatholder = (float)global.particle_y[i];
        fwrite(&floatholder, sizeof(float), 1, global.movie_file);

        /* angle */
        floatholder = global.particle_angle_rad[i];
        fwrite(&floatholder, sizeof(float), 1, global.movie_file);
    }

    for (i = 0; i < global.N_pinningsites; i++)
    {
        /* color */
        intholder = 4;
        fwrite(&intholder, sizeof(int), 1, global.movie_file);

        /* id */
        intholder = i + global.N_particle;
        fwrite(&intholder, sizeof(int), 1, global.movie_file);

        /* x, y */
        floatholder = (float)global.pinningsite_x[i];
        fwrite(&floatholder, sizeof(float), 1, global.movie_file);
        floatholder = (float)global.pinningsite_y[i];
        fwrite(&floatholder, sizeof(float), 1, global.movie_file);

        /* angle */
        floatholder = 0.0;
        fwrite(&floatholder, sizeof(float), 1, global.movie_file);
    }
    fflush(global.statistics_file);
}

void copy_parameters_file(int argc, char *argv[])
{
    FILE *parameterfile;
    FILE *savedparameterfile;

    char parameterfilename[300];
    char savedparameterfilename[300];
    char chartoread;
    char numread;

    strcpy(parameterfilename, argv[1]);
    strcpy(savedparameterfilename, "savedparams/save_param_");
    strcat(savedparameterfilename, global.baseName);
    strcat(savedparameterfilename, ".txt");

    printf("Creating a saved copy for %s \n -> in %s\n", parameterfilename, savedparameterfilename);

    parameterfile = fopen(parameterfilename, "rb");
    if (parameterfile == NULL)
    {
        printf("Cannot open parameter file %s\n", parameterfilename);
        fflush(stdout);
        exit(1);
    }

    savedparameterfile = fopen(savedparameterfilename, "wb");
    if (savedparameterfile == NULL)
    {
        printf("Cannot open parameter file %s\n", savedparameterfilename);
        fflush(stdout);
        exit(1);
    }

    while (!feof(parameterfile))
    {
        numread = fread(&chartoread, sizeof(char), 1, parameterfile);
        // printf("%c\n",chartoread);
        fwrite(&chartoread, sizeof(char), 1, savedparameterfile);
    }

    fclose(parameterfile);
    fclose(savedparameterfile);
}

/* Calculations */

void calculate_PBC_folded_distance(double *dx, double *dy, double x1, double y1, double x2, double y2)
{
    double dxtemp, dytemp;

    dxtemp = x2 - x1;
    dytemp = y2 - y1;

#ifndef non_periodic_boundary
    // Periodic Boundary Conditions Fold-Back
    if (dxtemp > global.SX2)
        dxtemp -= global.SX;
    if (dxtemp < -global.SX2)
        dxtemp += global.SX;
    if (dytemp > global.SY2)
        dytemp -= global.SY;
    if (dytemp < -global.SY2)
        dytemp += global.SY;
#endif

    *dx = dxtemp;
    *dy = dytemp;
}

void initialize_particles()
{
    int j, N;
    double min_dist;
    double x_temp, y_temp, theta_temp;
    double dx, dy, dr;
    char overlap;

    global.particle_x = (double *)malloc((global.N_particle + global.N_additional_particles) * sizeof(double));
    global.particle_y = (double *)malloc((global.N_particle + global.N_additional_particles) * sizeof(double));
    global.particle_R = (double *)malloc((global.N_particle + global.N_additional_particles) * sizeof(double));
    global.particle_color = (int *)malloc((global.N_particle + global.N_additional_particles) * sizeof(int));
    global.particle_fx = (double *)malloc((global.N_particle + global.N_additional_particles) * sizeof(double));
    global.particle_fy = (double *)malloc((global.N_particle + global.N_additional_particles) * sizeof(double));
    global.particle_dx_so_far = (double *)malloc((global.N_particle + global.N_additional_particles) * sizeof(double));
    global.particle_dy_so_far = (double *)malloc((global.N_particle + global.N_additional_particles) * sizeof(double));
    global.particle_dx = (double *)malloc((global.N_particle + global.N_additional_particles) * sizeof(double));
    global.particle_dy = (double *)malloc((global.N_particle + global.N_additional_particles) * sizeof(double));
    global.particle_motor_force = (double *)malloc((global.N_particle + global.N_additional_particles) * sizeof(double));
    global.particle_angle_rad = (double *)malloc((global.N_particle + global.N_additional_particles) * sizeof(double));
    global.particle_cosfi = (double *)malloc((global.N_particle + global.N_additional_particles) * sizeof(double));
    global.particle_sinfi = (double *)malloc((global.N_particle + global.N_additional_particles) * sizeof(double));
    global.particle_motor_ellapsed_time = (int *)malloc((global.N_particle + global.N_additional_particles) * sizeof(int));
    global.particle_motor_total_time = (int *)malloc((global.N_particle + global.N_additional_particles) * sizeof(int));

    global.particle_is_tumbling = (int *)malloc((global.N_particle + global.N_additional_particles) * __SIZEOF_INT__);
    global.particle_is_active = (int *)malloc((global.N_particle + global.N_additional_particles) * __SIZEOF_INT__);
    global.particle_tumbling_torque = (double *)malloc((global.N_particle + global.N_additional_particles) * __SIZEOF_DOUBLE__);

    global.cluster = (struct cluster *)malloc((global.N_particle + global.N_additional_particles) * sizeof(struct cluster));
    global.particle = (struct particle *)malloc((global.N_particle + global.N_additional_particles) * sizeof(struct particle));
    global.particle_neighbor_nr = (int *)malloc((global.N_particle + global.N_additional_particles) * sizeof(int));

    global.particle_eta = (double *)malloc((global.N_particle + global.N_additional_particles) * sizeof(double));

    // this means the discs cover about half of the system (0.4)
    // with random deposition we can go up to around 0.5
    min_dist = 2.0 * global.generic_particle_R;
    printf("Min dist = %lf Max radius = %lf\n", min_dist, min_dist / 2.0);

    for (int i = 0; i < global.N_particle; i++)
    {
        overlap = 1;

        while (overlap)
        {
            x_temp = global.SX * Rand();
            y_temp = global.SY * Rand();
            overlap = 0;

            for (j = 0; j < i; j++)
            {
                calculate_PBC_folded_distance(&dx, &dy, x_temp, y_temp, global.particle_x[j], global.particle_y[j]);
                dr = sqrt(dx * dx + dy * dy);
                if (dr < min_dist)
                {
                    overlap = 1;
                    break;
                }
            }
        }

        global.particle_x[i] = x_temp;
        global.particle_y[i] = y_temp;
        global.particle_angle_rad[i] = 2.0 * PI * Rand();
        global.particle_R[i] = global.generic_particle_R;
        global.particle_cosfi[i] = cos(global.particle_angle_rad[i]);
        global.particle_sinfi[i] = sin(global.particle_angle_rad[i]);

        global.particle_neighbor_nr[i] = 0;

        if (Rand() < global.generic_particle_active_fraction)
            global.particle_is_active[i] = 1;
        else
            global.particle_is_active[i] = 0;

        if (Rand() < global.obstacle_ratio)
        {
            global.particle_eta[i] = global.eta_large;
            global.particle_color[i] = 8;
        }
        else
        {
            global.particle_eta[i] = global.eta_small;
            global.particle_color[i] = 6;
        }

        global.particle_motor_force[i] = global.particle_is_active[i] * global.generic_particle_motor_force;
        global.particle_motor_ellapsed_time[i] = 0;
        global.particle_motor_total_time[i] = global.generic_particle_motor_minimum_time;
        global.particle_motor_total_time[i] += (int)floor(Rand() * (global.generic_particle_motor_maximum_time - global.generic_particle_motor_minimum_time));

        global.particle_fx[i] = 0.0;
        global.particle_fy[i] = 0.0;

        global.particle_dx_so_far[i] = 0.0;
        global.particle_dy_so_far[i] = 0.0;
        global.particle_dx[i] = 0.0;
        global.particle_dy[i] = 0.0;
    }

    printf("%d Particles initialized successfully\n", global.N_particle);
}

void add_particles_at_later_time()
{
    int i, j, N;
    double min_dist;
    double x_temp, y_temp, theta_temp;
    double dx, dy, dr;
    char overlap;

    int exit_from_loop = 0;
    int N_particles_to_be_added = global.N_additional_particles / global.N_times_to_add_more_particles;
    int tries = 0;


    min_dist = 2.0 * global.generic_particle_R;
    printf("Min dist = %lf Max radius = %lf\n", min_dist, min_dist / 2.0);

    i = global.N_particle;
    while ((i < global.N_particle + N_particles_to_be_added) && (exit_from_loop == 0)) {

        overlap = 1;

        while (overlap && tries <= 10000)
        {
            x_temp = global.SX * Rand();
            y_temp = global.SY * Rand();
            overlap = 0;

            for (j = 0; j < i; j++)
            {
                calculate_PBC_folded_distance(&dx, &dy, x_temp, y_temp, global.particle_x[j], global.particle_y[j]);
                dr = sqrt(dx * dx + dy * dy);
                if (dr < min_dist)
                {
                    overlap = 1;
                    break;
                }
            }
        }
        if (tries > 10000)        
        {
            printf("System too dense cannot place extra additional particles\n");
        exit_from_loop = 1;
        }
        else {

            global.particle_x[i] = x_temp;
            global.particle_y[i] = y_temp;
            global.particle_angle_rad[i] = 2.0 * PI * Rand();
            global.particle_R[i] = global.generic_particle_R;
            global.particle_cosfi[i] = cos(global.particle_angle_rad[i]);
            global.particle_sinfi[i] = sin(global.particle_angle_rad[i]);

            global.particle_neighbor_nr[i] = 0;

            if (Rand() < global.generic_particle_active_fraction)
                global.particle_is_active[i] = 1;
            else
                global.particle_is_active[i] = 0;

            global.particle_eta[i] = global.eta_large;
            global.particle_color[i] = 8;

            global.particle_motor_force[i] = global.particle_is_active[i] * global.generic_particle_motor_force;
            global.particle_motor_ellapsed_time[i] = 0;
            global.particle_motor_total_time[i] = global.generic_particle_motor_minimum_time;
            global.particle_motor_total_time[i] += (int)floor(Rand() * (global.generic_particle_motor_maximum_time - global.generic_particle_motor_minimum_time));

            global.particle_fx[i] = 0.0;
            global.particle_fy[i] = 0.0;

            global.particle_dx_so_far[i] = 0.0;
            global.particle_dy_so_far[i] = 0.0;
            global.particle_dx[i] = 0.0;
            global.particle_dy[i] = 0.0;
        }

        i++;
    }

    printf("Added %d more particles from %d total amount. Remaining particles to be added: %d\n", i - global.N_particle, N_particles_to_be_added, global.N_additional_particles - i + global.N_particle);
    global.N_particle = i;
    printf("%d total particles \n", global.N_particle);
    fflush(stdout);
}

void calculate_internal_motor_forces()
{
    int i, j, i_particle;
    double theta;

    for (i = 0; i < global.N_particle; i++)
    {
        if (global.particle_is_active[i])
        {
            global.particle_fx[i] += global.particle_motor_force[i] * global.particle_cosfi[i];
            global.particle_fy[i] += global.particle_motor_force[i] * global.particle_sinfi[i];
            global.particle_motor_ellapsed_time[i]++;

            // it finished running
            if (global.particle_motor_ellapsed_time[i] >= global.particle_motor_total_time[i])
            {
                global.particle_motor_total_time[i] = global.generic_particle_motor_minimum_time;
                global.particle_motor_total_time[i] += (int)floor(Rand() * (global.generic_particle_motor_maximum_time - global.generic_particle_motor_minimum_time));
                global.particle_motor_ellapsed_time[i] = 0;

                global.particle_angle_rad[i] = 2.0 * PI * Rand();
            }
        }
    }
}

void calculate_interparticle_forces()
{
    int i, j, ii;
    double dx, dy, dr, delta_compressed;
    double fx, fy;

    for (ii = 0; ii < global.N_Verlet; ii++)
    {
        i = global.Verlet_i[ii];
        j = global.Verlet_j[ii];
        calculate_PBC_folded_distance(&dx, &dy, global.particle_x[i], global.particle_y[i], global.particle_x[j], global.particle_y[j]);
        dr = sqrt(dx * dx + dy * dy);

        delta_compressed = dr - global.particle_R[i] - global.particle_R[j];

        if (delta_compressed < 0)
        {
            fx = delta_compressed * global.generic_particle_k_spring * dx / dr;
            fy = delta_compressed * global.generic_particle_k_spring * dy / dr;

            global.particle_fx[i] += fx;
            global.particle_fy[i] += fy;
            global.particle_fx[j] -= fx;
            global.particle_fy[j] -= fy;
        }
    }
}

void calculate_dragging_force()
{
    int i, j, i_particle;
    double theta;

    for (i = 0; i < global.N_particle; i++)
    {
        if (global.particle_eta[i] == global.eta_small)
        {
            global.particle_fx[i] += global.dragging_force_x;
            global.particle_fy[i] += global.dragging_force_y;
        }
    }
}

void calculate_thermal_noise()
{
    for (int i = 0; i < global.N_particle; i++)
    {
        global.particle_fx[i] += global.temperature * gasdev();
        global.particle_fy[i] += global.temperature * gasdev();
    }
}

void check_Verlet_need_to_rebuild(int i)
{
    double dr2;

    dr2 = global.particle_dx_so_far[i] * global.particle_dx_so_far[i] + global.particle_dy_so_far[i] * global.particle_dy_so_far[i];
    if (dr2 > global.Verlet_distance_to_rebuild * global.Verlet_distance_to_rebuild)
        global.flag_Verlet_needs_rebuild = 1;
}

/* Pinningsites */

void initialize_pinningsites()
{
    int i, j;
    double min_dist;
    double density;
    double x_temp, y_temp;
    double dx, dy, dr;
    char overlap;

    int attempts, max_attempts;

    min_dist = 3.5;

    printf("Minimum distance for pinning sites: %lf\n", min_dist);
    density = (double)(global.N_pinningsites * PI * global.generic_pinningsite_R * global.generic_pinningsite_R) / (global.SX * global.SY);

    printf("Pinning site density is = %lf\n", density);

    global.pinningsite_x = (double *)malloc(global.N_pinningsites * sizeof(double));
    global.pinningsite_y = (double *)malloc(global.N_pinningsites * sizeof(double));
    global.pinningsite_R = (double *)malloc(global.N_pinningsites * sizeof(double));
    global.pinningsite_fmax = (double *)malloc(global.N_pinningsites * sizeof(double));

    max_attempts = 10000;

    for (i = 0; i < global.N_pinningsites; i++)
    {

        overlap = 1;
        attempts = 0;

        while ((overlap) && (attempts <= max_attempts))
        {

            x_temp = global.SX * Rand();
            y_temp = global.SY * Rand();

            overlap = 0;

            for (j = 0; j < i; j++)
            {
                calculate_PBC_folded_distance(&dx, &dy, x_temp, y_temp, global.pinningsite_x[j], global.pinningsite_y[j]);
                dr = sqrt(dx * dx + dy * dy);
                if (dr < min_dist)
                {
                    overlap = 1;
                    break;
                }
            }
        }

        if (attempts > max_attempts)
        {
            printf("Cannot place pinning sites, too dense\n");
            fflush(stdout);
            exit(1);
        }

        // found a good position that does not overlap
        global.pinningsite_x[i] = x_temp;
        global.pinningsite_y[i] = y_temp;

        // radius of the pinningsite
        global.pinningsite_R[i] = global.generic_pinningsite_R;
        global.pinningsite_fmax[i] = global.generic_pinningsite_fxmax;
    }

    printf("%d pinning sites initialized successfully\n", global.N_pinningsites);
}

void add_pinningsite_to_grid(int i, int j, int l)
{
    int k;

    global.pinningsite_grid[i][j][0]++;
    k = global.pinningsite_grid[i][j][0];
    global.pinningsite_grid[i][j] = (int *)realloc(global.pinningsite_grid[i][j], (k + 1) * sizeof(int));
    global.pinningsite_grid[i][j][k] = l;
}

void initialize_pinning_grid()
{
    int i, j, l;

    global.pinning_grid_cell_size = 2.0 * global.generic_pinningsite_R;
    global.pinningsite_grid_Nx = (int)floor(global.SX / global.pinning_grid_cell_size) + 1;
    global.pinningsite_grid_Ny = (int)floor(global.SY / global.pinning_grid_cell_size) + 1;

    printf("Pinning grid size = %lf\n", global.pinning_grid_cell_size);
    printf("Pinning grid size %d x %d\n", global.pinningsite_grid_Nx, global.pinningsite_grid_Ny);

    global.pinningsite_grid = (int ***)malloc(global.pinningsite_grid_Nx * sizeof(int **));
    for (i = 0; i < global.pinningsite_grid_Nx; i++)
    {
        global.pinningsite_grid[i] = (int **)malloc(global.pinningsite_grid_Ny * sizeof(int *));
        for (j = 0; j < global.pinningsite_grid_Ny; j++)
        {
            // first number determines the number of pinningsites in a grid cell
            global.pinningsite_grid[i][j] = (int *)malloc(sizeof(int));
            global.pinningsite_grid[i][j][0] = 0;
        }
    }

    for (l = 0; l < global.N_pinningsites; l++)
    {
        i = (int)floor(global.pinningsite_x[l] / global.pinning_grid_cell_size);
        j = (int)floor(global.pinningsite_y[l] / global.pinning_grid_cell_size);

        add_pinningsite_to_grid(i, j, l);
    }
}

void interaction_particle_pinningsiteGridCell(int k, int i, int j)
{
    // k - particle id
    // i,j - pinning grid cell coordinates
    double dx, dy;
    double pinning_force_x, pinning_force_y;
    int l, i_pinningsite;

    for (l = 1; l <= global.pinningsite_grid[i][j][0]; l++)
    {
        i_pinningsite = global.pinningsite_grid[i][j][l];
        calculate_PBC_folded_distance(&dx, &dy, global.particle_x[k], global.particle_y[k], global.pinningsite_x[i_pinningsite], global.pinningsite_y[i_pinningsite]);

        if ((dx * dx + dy * dy) <= global.pinningsite_R[i_pinningsite] * global.pinningsite_R[i_pinningsite])
        {
            pinning_force_x = global.pinningsite_fmax[i_pinningsite] * dx / global.pinningsite_R[i_pinningsite];
            pinning_force_y = global.pinningsite_fmax[i_pinningsite] * dy / global.pinningsite_R[i_pinningsite];
            global.particle_fx[k] += pinning_force_x;
            global.particle_fy[k] += pinning_force_y;
        }
    }
}

void calculate_pinning_forces()
{
    int i, j, i_left, i_right, j_down, j_up, k;

    for (k = 0; k < global.N_particle; k++)
    {

        i = (int)floor(global.particle_x[k] / global.pinning_grid_cell_size);
        j = (int)floor(global.particle_y[k] / global.pinning_grid_cell_size);

        i_left = i - 1;
        i_right = i + 1;
        j_down = j - 1;
        j_up = j + 1;

        // PBC foldback
        if (i_left < 0)
            i_left = global.pinningsite_grid_Nx - 1;
        else if (i_right >= global.pinningsite_grid_Nx)
            i_right = 0;
        if (j_down < 0)
            j_down = global.pinningsite_grid_Ny - 1;
        else if (j_up >= global.pinningsite_grid_Ny)
            j_up = 0;

        interaction_particle_pinningsiteGridCell(k, i, j);
        interaction_particle_pinningsiteGridCell(k, i, j_down);
        interaction_particle_pinningsiteGridCell(k, i, j_up);
        interaction_particle_pinningsiteGridCell(k, i_left, j);
        interaction_particle_pinningsiteGridCell(k, i_left, j_down);
        interaction_particle_pinningsiteGridCell(k, i_left, j_up);
        interaction_particle_pinningsiteGridCell(k, i_right, j);
        interaction_particle_pinningsiteGridCell(k, i_right, j_down);
        interaction_particle_pinningsiteGridCell(k, i_right, j_up);
    }
}

void write_pinningsites()
{
    FILE *testout;
    int i;
    char fname[300];

    sprintf(fname, "pinningsites/%s.txt", global.baseName);
    testout = fopen(fname, "wt");
    for (i = 0; i < global.N_pinningsites; i++)
        fprintf(testout, "%lf %lf\n", global.pinningsite_x[i], global.pinningsite_y[i]);
    fclose(testout);
    printf("%d Pinning sites written to testfile\n", global.N_pinningsites);
}

/* Updates */

void rebuild_Verlet()
{
    int i, j;
    double dx, dy;
    double dr2;
    double Verlet_R2;

    // discs are in contact at a distance less than 2R
    // the distance to which the Verlet list goes out is 4R (R usually is 1.0)
    Verlet_R2 = global.Verlet_R * global.Verlet_R; // square should be 16.0

    global.N_Verlet = 0;
    global.Verlet_i = (int *)realloc(global.Verlet_i, global.N_Verlet * sizeof(int));
    global.Verlet_j = (int *)realloc(global.Verlet_j, global.N_Verlet * sizeof(int));

    for (i = 0; i < global.N_particle; i++)
        for (j = i + 1; j < global.N_particle; j++)
        {
            calculate_PBC_folded_distance(&dx, &dy, global.particle_x[i], global.particle_y[i], global.particle_x[j], global.particle_y[j]);
            dr2 = dx * dx + dy * dy;
            if (dr2 < Verlet_R2)
            {
                global.N_Verlet++;
                global.Verlet_i = (int *)realloc(global.Verlet_i, global.N_Verlet * sizeof(int));
                global.Verlet_j = (int *)realloc(global.Verlet_j, global.N_Verlet * sizeof(int));
                global.Verlet_i[global.N_Verlet - 1] = i;
                global.Verlet_j[global.N_Verlet - 1] = j;
            }
        }

    // refresh the distance couters
    for (i = 0; i < global.N_particle; i++)
    {
        global.particle_dx_so_far[i] = 0.0;
        global.particle_dy_so_far[i] = 0.0;
    }
    global.flag_Verlet_needs_rebuild = 0;
}

void move_particles()
{

    int i, j;
    double dx, dy, dr;

    for (i = 0; i < global.N_particle; i++)
    {
        dx = 1 / global.particle_eta[i] * global.particle_fx[i] * global.dt;
        dy = 1 / global.particle_eta[i] * global.particle_fy[i] * global.dt;

        // gather the dx,dy so far
        global.particle_dx_so_far[i] += dx;
        global.particle_dy_so_far[i] += dy;

        check_Verlet_need_to_rebuild(i);

        global.particle_x[i] += dx;
        global.particle_y[i] += dy;

        // PBC put the long particle back into the box if it exited the box
        if (global.particle_x[i] > global.SX)
            global.particle_x[i] -= global.SX;
        if (global.particle_y[i] > global.SY)
            global.particle_y[i] -= global.SY;
        if (global.particle_x[i] <= 0.0)
            global.particle_x[i] += global.SX;
        if (global.particle_y[i] <= 0.0)
            global.particle_y[i] += global.SY;

        /* new motion angle for particle*/
        // global.particle_angle_rad[i] += global.particle_tumbling_torque[i] * global.dt;

        if (global.particle_angle_rad[i] >= 2 * PI)
        {
            global.particle_angle_rad[i] -= 2 * PI;
        }
        else if (global.particle_angle_rad[i] < 0)
        {
            global.particle_angle_rad[i] += 2 * PI;
        }

        global.particle_cosfi[i] = cos(global.particle_angle_rad[i]);
        global.particle_sinfi[i] = sin(global.particle_angle_rad[i]);

        if (global.particle_eta[i] == global.eta_small)
        {
            global.avg_fx += global.particle_fx[i];
            global.avg_fy += global.particle_fy[i];
        }
    }

    // zero all forces on particles
    for (i = 0; i < global.N_particle; i++)
    {
        global.particle_fx[i] = 0;
        global.particle_fy[i] = 0;
    }
}

/* Statistics */

void clusterize_particles()
{
    int i, j, k, cluster_i, cluster_j;
    double dx, dy;
    struct particle *small_cluster_particle;
    int particle_i;
    int first_particle_i;
    int last_particle_i;
    int next_particle_i;

    for (i = 0; i < global.N_particle; i++)
    {
        global.particle[i].clusterID = i;
        global.particle[i].next = NULL;

        global.cluster[i].size = 1;
        global.cluster[i].first = global.particle + i;
        global.cluster[i].last = global.particle + i;
    }

    // linking clusters together
    for (int ii = 0; ii < global.N_Verlet; ii++)
    {
        i = global.Verlet_i[ii];
        j = global.Verlet_j[ii];

        calculate_PBC_folded_distance(&dx, &dy, global.particle_x[i], global.particle_y[i], global.particle_x[j], global.particle_y[j]);
        // two partciles from different clusters touch
        if ((global.particle[i].clusterID != global.particle[j].clusterID) && (dx * dx + dy * dy <= (global.particle_R[i] + global.particle_R[j]) * (global.particle_R[i] + global.particle_R[j])))
        {

            if (global.cluster[global.particle[i].clusterID].size > global.cluster[global.particle[j].clusterID].size)
            {
                cluster_i = global.particle[i].clusterID; // cluster_i is the larger
                cluster_j = global.particle[j].clusterID;
            }
            else
            {
                cluster_i = global.particle[j].clusterID;
                cluster_j = global.particle[i].clusterID;
            }

            // update small cluster cluster_IDs
            small_cluster_particle = global.cluster[cluster_j].first;
            while (small_cluster_particle)
            {
                small_cluster_particle->clusterID = cluster_i;
                small_cluster_particle = small_cluster_particle->next;
            }

            // set large cluster
            global.cluster[cluster_i].size += global.cluster[cluster_j].size;
            global.cluster[cluster_i].last->next = global.cluster[cluster_j].first;
            global.cluster[cluster_i].last = global.cluster[cluster_j].last;

            // empty small cluster
            global.cluster[cluster_j].size = 0;
            global.cluster[cluster_j].first = NULL;
            global.cluster[cluster_j].last = NULL;
        }
    }
}

int largest_clusterID()
{
    int largest_clusterID = 0;

    for (int i = 1; i < global.N_particle; i++)
    {
        if (global.cluster[largest_clusterID].size < global.cluster[i].size)
        {
            largest_clusterID = i;
        }
    }

    return largest_clusterID;
}

void calculate_statistics()
{
    global.largest_cluster_size_sum += global.cluster[largest_clusterID()].size;
}

void write_statistics()
{
    double largest_cluster_size_avg;
    double avg_fx, avg_fy;

    avg_fx = global.avg_fx / (double)global.N_particle / (double)global.statistics_time;
    avg_fy = global.avg_fy / (double)global.N_particle / (double)global.statistics_time;
    largest_cluster_size_avg = global.largest_cluster_size_sum / (double)global.statistics_time;

    fprintf(global.statistics_file, "%d %f %f %f\n", global.time, largest_cluster_size_avg, avg_fx, avg_fy);
    // fprintf(global.statistics_file, "%d %f\n", global.time, largest_cluster_size_avg);
    fflush(global.statistics_file);
    global.largest_cluster_size_sum = 0;
    global.avg_fx = 0;
    global.avg_fy = 0;
}

/* Simulation */

void run_simulation()
{
    double time_ellapsed_percent;
    clock_t t_start = clock();

    int N_bands = 50;
    double *bands = (double *)calloc(N_bands, sizeof(double));

    for (global.time = 0; global.time < global.total_runtime; global.time++)
    {
        calculate_interparticle_forces();
        calculate_thermal_noise();

        if (global.time < global.relax_time || global.time >= global.relax_time + global.relax_duration)
        {
            calculate_dragging_force();
        }

        if (global.time % global.echo_time == 0)
        {
            time_ellapsed_percent = (double)global.time / global.total_runtime * 100.0;
            printf("%d/%d (%.2lf%%; %s Remaining: %.2lf s)\n", global.time, global.total_runtime,
                   time_ellapsed_percent, global.baseName, (100.0 - time_ellapsed_percent) * ((double)(clock() - t_start) / CLOCKS_PER_SEC) / time_ellapsed_percent);
            fflush(stdout);
        }

        if (global.time % global.movie_time == 0)
        {
            clusterize_particles();
            write_cmovie_frame();
        }

        calculate_statistics();
        move_particles();

        if (global.time % global.statistics_time == 0)
        {
            write_statistics();
        }

        if (global.flag_Verlet_needs_rebuild == 1)
        {
            rebuild_Verlet();
        }

        if ( global.time && (global.time % global.N_steps_to_add_more_particles == 0) && 
            global.time/global.N_steps_to_add_more_particles <= global.N_times_to_add_more_particles)
        {
            printf("Adding extra particles\n");
            fflush(stdout);
            add_particles_at_later_time();
            rebuild_Verlet();
            clusterize_particles();
            initialize_statistics();
        }
    }

    printf("\r%d/%d (100.00%%; ETA: 0.00 s)\n", global.time, global.total_runtime);
    printf("                     in %.2lf s\r", (double)(clock() - t_start) / CLOCKS_PER_SEC);
}

int main(int argc, char *argv[])
{
    initialize_from_file(argc, argv);
    copy_parameters_file(argc, argv);
    initialize_output_files();

    initialize_particles();
    initialize_statistics();

    rebuild_Verlet();

    run_simulation();

    fclose(global.movie_file);
    fclose(global.statistics_file);
}
