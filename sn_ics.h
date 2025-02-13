#ifndef SN_ICS_H
#define SN_ICS_H

#include <stdio.h>
#include "sn_math.h"
#include "sn_cross_sections.h"
#include "math.h"

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
#endif