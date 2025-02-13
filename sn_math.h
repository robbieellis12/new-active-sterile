#ifndef sn_math_h
#define sn_math_h

#define _USE_MATH_DEFINES
#define TS_STEP_SIZE 1/50.0//1/40.0 //40 50
#define TS_NUM_PTS 150 //120 //120 //120 200
#define TWO_PI_2 39.4784176044
#define TWO_PI_3 248.050213442
#define TWO_PI_4 1558.54545654
#define PI_4 97.409091034
#define SQRT2 1.41421356237
#define PI_4 97.409091034

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <omp.h>

extern double ts_quad_fixed_weights[2*TS_NUM_PTS+1];
extern double ts_quad_fixed_nodes[2*TS_NUM_PTS+1];
struct polylogs
{
    double* li2;
    double* li3;
    double* li4;
    double* dli2;
    double* dli3;
    double* dli4;
    
    double x0;
    double x1;
    double h;
    
    int N;
};

double coth(double x);

struct n_rho_explicit {double n; double nb;double nc; double p; double pb; double pc;};
struct n_rho_explicit2 {double n; double nT;double nc; double p; double pT; double pc; double other_fphi_deriv;};
struct n_p_vals{double n_s; double n_s_T; double n_s_c; double p_s; double p_s_T; double p_s_c; double In_p; double In_p_b; double In_p_c; double Ip_p; double Ip_p_b; double Ip_p_c;};
struct n_p_vals2{double n_s; double n_s_T; double n_s_c; double p_s; double p_s_T; double p_s_c; double In_p; double In_p_T; double In_p_c; double Ip_p; double Ip_p_T; double Ip_p_c; double other_fphi_deriv;};


void read_paired_list(char* filename,double** xdat, double** ydat,int* len);
double two_point_2d_hermite_interp(double x0, double y0, double x1, double x2, double y1, double y2, double* f, double* fx, double* fy, double* fxy);

double interpolate(double* x,double*y,int L,double x_0);

void ts_quad_info(double*, double*, double, int);

double ts_quad(double (*f)(double, void*), void*,  double, double, double, int);

double L1(double, double, double);

double L1p(double, double, double);

double L2(double, double, double);

double L2p(double, double, double);

double K1(double, double, double);

double K1p(double, double, double);

double K2(double, double, double);

double K2p(double, double, double);

double two_point_hermite_interp(double, double, double, double, double, double, double);

void fill_polylog_table(struct polylogs* plgs);
double interp_Li2(double x,struct polylogs* plgs);
double interp_Li3(double x,struct polylogs* plgs);
double interp_Li4(double x,struct polylogs* plgs);
void ts_quad_li(double x, double N, double* li,double* dli, double h, int n);
void free_polylog_table(struct polylogs* plgs);

struct n_rho_explicit ts_quad_n_rho_explicit(double b,double c);
struct n_rho_explicit2 ts_quad_n_rho_explicit2(double m,double T, double c);
double n_rho_fermi(double m, double T);
struct n_p_vals get_n_p_vals(double m, double Ts, double cs, double Tp, double cp,struct polylogs* plgs);
struct n_p_vals2 get_n_p_vals2(double m, double Ts, double cs, double Tp, double cp,struct polylogs* plgs);

#endif