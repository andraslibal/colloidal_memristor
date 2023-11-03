#ifndef __MEMRISTOR_H__
#define __MEMRISTOR_H__

struct particle
{
    int i;                        
    int clusterID;                
    struct particle *next; 
};

struct cluster
{
    int size;                 
    struct particle *first; 
    struct particle *last;
};

struct global_struct {

    /*  General */

    char baseName[200]; // file basename
    uint64_t seedToSet; // seed modifier for multiple runs
    int time;

    FILE *movie_file;
    FILE *statistics_file;
    FILE *histogram_file;
    
    /* System */

    struct particle *particle;
    struct cluster *cluster;

    double SX, SY;   // system size X,Y direction
    double SX2, SY2; // half of the system size X,Y direction

    int N_particle;

    double dt;
    int total_runtime;
    int echo_time;
    int movie_time;
    int statistics_time;

    int *Verlet_i;
    int *Verlet_j;
    int N_Verlet;
    char flag_Verlet_needs_rebuild;
    double Verlet_distance_to_rebuild;
    double Verlet_R;

    /* Pinningsites */
    int N_pinningsites;

    double *pinningsite_x;
    double *pinningsite_y;
    double *pinningsite_R;
    double *pinningsite_fmax;

    int ***pinningsite_grid;
    double pinning_grid_cell_size;
    int pinningsite_grid_Nx, pinningsite_grid_Ny;

    double generic_pinningsite_R;
    double generic_pinningsite_fxmax;

    /* Particles */

    double *particle_x;
    double *particle_y;
    double *particle_R;
    double *particle_angle_rad;
    double *particle_cosfi;
    double *particle_sinfi;
    int *particle_color;

    double *particle_fx;
    double *particle_fy;
    double *particle_dx_so_far;
    double *particle_dy_so_far;
    double *particle_dx;
    double *particle_dy;

    double *particle_motor_force;
    int *particle_motor_ellapsed_time;
    int *particle_motor_total_time;

    double *particle_tumbling_torque;
    int *particle_is_tumbling;
    int *particle_is_active;

    double *particle_eta;

    /* Generic */

    double generic_particle_motor_force;
    int generic_particle_motor_minimum_time;
    int generic_particle_motor_maximum_time;

    double generic_particle_R;
    double generic_particle_k_spring;
    double generic_particle_active_fraction;

    double eta_small;
    double eta_large;
    double obstacle_ratio;
    double dragging_force_x;
    double dragging_force_y;

    /* Statistics */

    int largest_cluster_size_sum;
    int *particle_neighbor_nr;
    double avg_fx;
    double avg_fy;
};

#endif