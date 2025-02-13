#ifndef DERIVS_H
#define DERIVS_H
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "sn_computations.h"
#include "sn_math.h"
#include "sn_cross_sections.h"
struct massless_rates_derivs{
    double g_aa_ss;
    double g_aa_ss_T;
    double g_aa_ss_c;
    double g_ss_aa;
    double g_ss_aa_T;
    double g_ss_aa_c;
    double g_as_as;
    double g_as_as_T;
    double g_as_as_c;
    double g_sa_as;
    double g_sa_as_T;
    double g_sa_as_c;
    double g_aa_pp;
    double g_aa_pp_T;
    double g_aa_pp_c;
    double g_ss_pp;
    double g_ss_pp_T;
    double g_ss_pp_c;
};
struct collision_vals_derivs
{
    double coll;
    double moment;
    double coll_dT;
    double coll_dc;
    double moment_dT;
    double moment_dc;

};

struct vp_derivs
{
    double vp;
    double vp_dT;
    double vp_dc;
};

struct va_derivs
{
    double va;
    double va_dT;
    double va_dc;
};


struct collision_vals_derivs Collision_deriv(double m, double T, double c);
struct va_derivs va_deriv(double m, double T, double p, double c,struct polylogs* plgs);
struct vp_derivs vp_deriv(double m, double T, double p, double c);
#endif