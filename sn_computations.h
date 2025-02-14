#ifndef sn_computations_h
#define sn_computations_h

#include "sn_math.h"
#include "sn_cross_sections.h"

#define FERMI_CONST 1.1663787*1e-23
#define M_W_BOSON 80.377e3
#define M_Z_BOSON 91.1876e3
#define M_PL 1.221e22
#define RHO_CDM 9.337e-12

struct fixed_params
{
    double m_phi;
    double m_s;
    double y;
    double theta;
    double m_l_a;

    double* T_g;
    double* g_s;
    int len_g;
};

struct massless_rates{
    double g_aa_ss;
    double g_ss_aa;
    double g_as_as;
    double g_sa_as;
    double g_aa_pp;
    double g_ss_pp;
};

struct collision_vals
{
    double coll;
    double moment;
};

struct derivs{
double C_n_s;
double C_p_s;
double C_n_p;
double C_p_p;
};
struct sp_params
{
    double T;
    double Ts;
    double Tphi;
    double c_s;
    double c_phi;
};

struct massive_interpolation_table
{
    double* table;
    double* d_table;
    
    double x0;
    double x1;
    double h;
    
    //redundant
    int N;
    };

struct collision_vals Moment1(double m, double T, double c);
double Moment2(double m, double T, double c);
double va(double w, double c,struct polylogs* plgs);
double vp(double a,double b,double c,double w);
double v_tot(double p, double m_a, double m_phi, double y, struct sp_params* sp_params,struct polylogs* plgs);
double oscillation_peak_func(double p, double m_s, double cos_2theta, double m_a, double m_phi, double y, struct sp_params* sp_params,struct polylogs* plgs);
double find_oscillation_peak(double m_a, double m_s, double cos_2theta, double m_phi, double y, struct sp_params* sp_params,struct polylogs* plgs);

void fill_massive_tables(struct massive_interpolation_table* table, double y);
double interp_massive(double i0,double i1,double y, struct massive_interpolation_table* mip);
struct massless_rates compute_massless_rates(double T,double Ts,double w,double w_s,double c_s,double y,struct polylogs* plgs);
double compute_vp_vp_rate(double mphi,double p1,double cphi,double Tphi,double y,struct massive_interpolation_table* mip);
void free_massive_table(struct massive_interpolation_table* mip);
double C_as_as_narrow(double p1, double c, double m, double T, double y);

double oscillation_term(double p, struct fixed_params fixed_params,struct sp_params sp_params, struct polylogs* plgs,struct massive_interpolation_table* mip);
struct derivs integrated_de(double xi, double xf, struct fixed_params fixed_params, struct sp_params sp_params, struct polylogs* plgs,struct massive_interpolation_table* mip);

double H(double T, struct fixed_params* snp);
struct derivs integrated_de_omp(double xi, double xf, struct fixed_params fixed_params, struct sp_params sp_params, struct polylogs* plgs,struct massive_interpolation_table* mip);
struct derivs integrated_de_omp_eq(double xi, double xf, struct fixed_params fixed_params, struct sp_params sp_params, struct polylogs* plgs);
struct collision_vals Moment1_omp(double m, double T, double c);

#endif