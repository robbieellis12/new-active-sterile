#ifndef SN_ICS_H
#define SN_ICS_H

#include <stdio.h>
#include "sn_math.h"
#include "sn_cross_sections.h"
#include "math.h"
#include "sn_computations.h"
#include "diffeqs.h"
/*
struct massless_coll_table
{
    double T0;
    double T1;
    int NT;
    double hT;

    double c0;
    double c1;
    int Nc;
    double hc;

    double* table;
    double* table_dT;
    double* table_dc;
    double* table_dTc;

    double (*cs)(double,double);
    double m;
};
struct massless_coll_vals
{
    double f;
    double f_dT;
    double f_dc;
    double f_dTc;
};
void fill_massless_coll_table(struct massless_coll_table* mct,double y);
struct massless_coll_vals get_massless_coll_vals(double m, double T, double c,double y,double (*f)(double,double));
void free_massless_coll_table(struct massless_coll_table* mct);
double interp_massless_col(double m,double T,double c,struct massless_coll_table* mct);

double ts_quad_massless_coll(double m, double T, double c,double y,double (*f)(double,double));



*/


struct massless_coll_table
{
    double x0;
    double x1;
    int Nx;
    double hx;

    double c0;
    double c1;
    int Nc;
    double hc;

    double* table;
    double* table_dx;
    double* table_dc;
    double* table_dxc;

    double* table_moment;
    double* table_dx_moment;
    double* table_dc_moment;
    double* table_dxc_moment;

    double (*cs)(double,double);
    double m;
};
struct massless_coll_vals
{
    double f;
    double f_dx;
    double f_dc;
    double f_dxc;

    double f_moment;
    double f_dx_moment;
    double f_dc_moment;
    double f_dxc_moment;
};

struct massless_coll_return_vals
{
    double f;
    double f_moment;
};


struct massive_coll_vals
{
    double f;
    double f_db;
    double f_dc;
    double f_dbc;
    double f_moment;
    double f_db_moment;
    double f_dc_moment;
    double f_dbc_moment;
};

struct massive_coll_return_vals
{
    double f;
    double f_moment;
};


struct massive_coll_table
{
    double b0;
    double b1;
    int Nb;
    double hb;

    double c0;
    double c1;
    int Nc;
    double hc;

    double* table;
    double* table_db;
    double* table_dc;
    double* table_dbc;
    double* table_moment;
    double* table_db_moment;
    double* table_dc_moment;
    double* table_dbc_moment;

    double m;
};

struct interpolation_tables{

    struct massless_coll_table* aa_ss_table;
    struct massless_coll_table* ss_aa_table;
    struct massless_coll_table* aa_pp_table;
    struct massless_coll_table* ss_pp_table;
    struct massive_coll_table* pp_vv_table;
    struct polylogs* plgs;
};

struct dx_coll_params
{
    struct fixed_params* fixed_params;
};


void fill_massless_coll_table(struct massless_coll_table* mct,double y);
struct massless_coll_vals get_massless_coll_vals(double x, double c,double y,double (*f)(double,double));
void free_massless_coll_table(struct massless_coll_table* mct);
struct massless_coll_return_vals interp_massless_col(double x,double c,double y, struct massless_coll_table* mct);

struct massless_coll_return_vals ts_quad_massless_coll(double x, double c,double y,double (*f)(double,double));


struct massive_coll_vals get_collision_vals(double m, double b, double c);

struct massive_coll_vals get_collision_vals2(double m, double b, double c);

void fill_massive_coll_table(struct massive_coll_table* mct,double y);
void free_massive_coll_table(struct massive_coll_table* mct);
//struct collision_vals ts_quad_massive_coll(double m, double b, double c,double y);
struct collision_vals interp_massive_colls(double m,double b,double c,struct massive_coll_table* mct);
struct derivs integrated_de_omp_colls(double xi, double xf, struct fixed_params fixed_params, struct sp_params sp_params,struct interpolation_tables* tables);
struct final_derivs step_omp_coll(struct sp_params sp, struct fixed_params* fixed_params,struct interpolation_tables* tables);
#endif