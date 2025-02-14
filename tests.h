#ifndef TESTS_H
#define TESTS_H
#include <stdio.h>
#include <math.h>
#include <omp.h>
#include "sn_computations.h"
#include "derivs.h"
#include "sn_ics.h"
#include "sundials_solve.h"

struct collision_data{
    double C_aa_ss;
    double C_ss_aa;
    double C_aa_pp;
    double C_ss_pp;
    double C_aa_ss_moment;
    double C_ss_aa_moment;
    double C_aa_pp_moment;
    double C_ss_pp_moment;
    double C_as_p_num;
    double C_as_p;
    double C_p_as;
    double C_p_as_moment;
    double C_as_p_moment;
    double oscillation;
    double oscillation_moment;
};


void export_doubles(char* location, double* data, int size);

void generate_lin_spaced(double* x, double a, double b, int N);
double g_data(char* loc, int N,double a,double b, struct fixed_params fixed_params,struct sp_params sp_params, struct polylogs* plgs,struct massive_interpolation_table* mip);

void test_rates(double p, struct fixed_params fixed_params,struct sp_params sp_params, struct polylogs* plgs,struct massive_interpolation_table* mip);
void send_collision_data(char* loc,struct fixed_params fixed_params, struct sp_params sp_params, struct polylogs* plgs,struct massive_interpolation_table* mip);
void moment_test();
void moment_derivs_test();
void va_derivs_test(struct polylogs* plgs);
//void test_collision_interp(double m, double T, double c, struct massless_coll_table* mct);

void test_collision_interp();

void test_massive_colls();
void test_massive_coll_speed();
void test_collision_step();
#endif