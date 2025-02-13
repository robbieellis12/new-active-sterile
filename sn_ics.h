#ifndef SN_ICS_H
#define SN_ICS_H

#include <stdio.h>
#include "sn_math.h"
#include "sn_cross_sections.h"
#include "math.h"
#include "sn_computations.h"

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

    double (*cs)(double,double);
    double m;
};
struct massless_coll_vals
{
    double f;
    double f_dx;
    double f_dc;
    double f_dxc;
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


void fill_massless_coll_table(struct massless_coll_table* mct,double y);
struct massless_coll_vals get_massless_coll_vals(double x, double c,double y,double (*f)(double,double));
void free_massless_coll_table(struct massless_coll_table* mct);
double interp_massless_col(double x,double c,struct massless_coll_table* mct);

double ts_quad_massless_coll(double x, double c,double y,double (*f)(double,double));


struct massive_coll_vals get_collision_vals(double m, double b, double c);

struct massive_coll_vals get_collision_vals2(double m, double b, double c);

void fill_massive_coll_table(struct massive_coll_table* mct,double y);
void free_massive_coll_table(struct massive_coll_table* mct);
//struct collision_vals ts_quad_massive_coll(double m, double b, double c,double y);
struct collision_vals interp_massive_colls(double m,double b,double c,struct massive_coll_table* mct);
#endif