#ifndef DIFFEQS_H
#define DIFFEQS_H
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "sn_computations.h"
#include "sn_math.h"
#include "sn_cross_sections.h"

struct dx_params
{
    struct polylogs* plgs;
    struct fixed_params* fixed_params;
    struct massive_interpolation_table* mip;
};

struct dx_eq_params
{
    double T;
    struct polylogs* plgs;
    struct fixed_params* fixed_params;
    struct massive_interpolation_table* mip;
};

struct final_derivs{
double dTs;
double dcs;
double dTp;
double dcp;
};
struct final_derivs_eq{
double dTs;
double dcs;
double dTp;
double dcp;
double n_s;
double p_s;
double n_p;
double p_p;
};
struct final_derivs step(struct sp_params sp, struct dx_params* dx_params);
struct solved_params
{
    double* T;
    double* Ts;
    double* Tphi;
    double* c_s;
    double* c_phi;
};

struct solved_params rk_solve(double T0, struct dx_params* dx_params,double h, double N);
struct solved_params rk_solve_omp(double T0, struct dx_params* dx_params,double h, double N);
struct final_derivs step_omp(struct sp_params sp, struct dx_params* dx_params);
struct final_derivs step_omp2(struct sp_params sp, struct dx_params* dx_params);


struct final_derivs_eq step_omp_eq(struct sp_params sp, struct dx_eq_params* dx_params);
#endif