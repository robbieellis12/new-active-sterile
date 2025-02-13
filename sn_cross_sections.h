#ifndef sn_cross_sections_h
#define sn_cross_sections_h

#include <stdio.h>
#include <math.h>
#include <omp.h>

extern const double zeta3_sn;
extern const double zeta5_sn;
extern const double zeta7_sn;
extern const double m_pl;
extern const double VI_2_LINEAR_CONST; //(1-gamma+ln(pi/2)) gamma=euler-marshoni (however its spelled) constant

double c_as_as(double, double);
double c_as_as_num(double u);
double c_vv_pp(double, double);
double c_vv_vv(double, double);
double c_pp_vv(double u,double y);
double c_vp_vp(double u,void* y);

#endif