#include "derivs.h"

struct massless_rates_derivs compute_massless_rates_derivs(double T,double Ts,double w,double w_s,double c_s,double y,struct polylogs* plgs)
{

    int n=TS_NUM_PTS;
    double int_val, weight, u;
    int_val = 0.0;
    double xf=100;
    double int_base_a;
    double int_base_s;

    double li_term_a;
    double li_term_s;

    double g_aa_ss,g_ss_aa, g_as_as, g_sa_as, g_aa_pp, g_ss_pp, g_pv_pv;
    g_aa_pp=0;
    g_ss_pp=0;
    g_as_as=0;
    g_sa_as=0;
    g_aa_ss=0;
    g_ss_aa=0;

    double g_aa_pp_c=0;
    double g_aa_pp_T=0;
    double g_ss_pp_c=0;
    double g_ss_pp_T=0;
    double g_aa_ss_c=0;
    double g_aa_ss_T=0;
    double g_ss_aa_c=0;
    double g_ss_aa_T=0;
    double g_as_as_c=0;
    double g_as_as_T=0;
    double g_sa_as_c=0;
    double g_sa_as_T=0;

    double c_vv_vv_T=0;
    double c_vv_vv_Ts=0;

    double c_as_as_T=0;
    double c_as_as_Ts=0;

    double c_vv_pp_T=0;
    double c_vv_pp_Ts=0;

    double f_a;
    double f_s;

double expterm=0;
double expterm_cs=0;
    for (int k = -n; k < n; ++k)
    {
        u = (xf * (ts_quad_fixed_nodes[k + n] + 1) ) / 2;
        weight=ts_quad_fixed_weights[k+n];
        expterm=exp(-u);
        expterm_cs=exp(c_s-u);

        f_a=u*expterm/(1-expterm); //I use the variables defined above because exp is more expensive to compute than using algebra
        f_s=u*expterm_cs/(1-expterm_cs);

        c_vv_vv_T=c_vv_vv(u/w,y);
        c_vv_vv_Ts=c_vv_vv(u/w_s,y);
        
        c_as_as_T=c_as_as(u/w,y);
        c_as_as_Ts=c_as_as(u/w_s,y);

        c_vv_pp_T=c_vv_pp(u/w,y);
        c_vv_pp_Ts=c_vv_pp(u/w_s,y);

        int_base_a=log(1+expterm)*u;
        int_base_s=log(1+expterm_cs)*u;
        //li_term_a=interp_Li2(-expterm, plgs);
        //li_term_s=interp_Li2(-expterm_cs, plgs);

        g_aa_ss+=weight*int_base_a*c_vv_vv_T;
        g_aa_ss_T+=weight*u*f_a*c_vv_vv_T;
        g_aa_ss_c+=weight*f_a*c_vv_vv_T;

        g_ss_aa+=weight*int_base_s*c_vv_vv_Ts;
        g_ss_aa_T+=weight*u*f_s*c_vv_vv_Ts;
        g_ss_aa_c+=weight*f_s*c_vv_vv_Ts;

        g_as_as+=weight*int_base_s*c_as_as_Ts;
        g_as_as_T+=weight*u*f_s*c_as_as_Ts;
        g_as_as_c+=weight*f_s*c_as_as_Ts;

        g_sa_as+=weight*int_base_a*c_as_as_T;
        g_sa_as_T+=weight*u*f_s*c_as_as_T;
        g_sa_as_c+=weight*f_s*c_as_as_T;

        g_aa_pp+=weight*int_base_a*c_vv_pp_T;
        g_aa_pp_T+=weight*u*f_a*c_vv_pp_T;
        g_aa_pp_c+=weight*f_a*c_vv_pp_T;

        g_ss_pp+=weight*int_base_s*c_vv_pp_Ts;
        g_ss_pp_T+=weight*u*f_s*c_vv_pp_Ts;
        g_ss_pp_c+=weight*f_s*c_vv_pp_Ts;

    }

    struct massless_rates_derivs rates;

    rates.g_aa_ss=g_aa_ss*xf/2;
    rates.g_aa_ss_T=g_aa_ss_T*xf/2;
    rates.g_aa_ss_c=g_aa_ss_c*xf/2;

    rates.g_ss_aa=g_ss_aa*xf/2;
    rates.g_ss_aa_T=g_ss_aa_T*xf/2;
    rates.g_ss_aa_c=g_ss_aa_c*xf/2;

    rates.g_as_as=g_as_as*xf/2;
    rates.g_as_as_T=g_as_as_T*xf/2;
    rates.g_as_as_c=g_as_as_c*xf/2;

    rates.g_sa_as=g_sa_as*xf/2;
    rates.g_sa_as_T=g_sa_as_T*xf/2;
    rates.g_sa_as_c=g_sa_as_c*xf/2;

    rates.g_aa_pp=g_aa_pp*xf/2;
    rates.g_aa_pp_T=g_aa_pp_T*xf/2;
    rates.g_aa_pp_c=g_aa_pp_c*xf/2;

    rates.g_ss_pp=g_ss_pp*xf/2;
    rates.g_ss_pp_T=g_ss_pp_T*xf/2;
    rates.g_ss_pp_c=g_ss_pp_c*xf/2;

    return rates;
}




struct vp_derivs vp_deriv(double m, double T, double p, double c)
{
    /*
    if(c>=b)
    {
        return -1;
    }
    if(a>90)
    {


    int n=TS_NUM_PTS;
    double h=TS_STEP_SIZE;

    double int_val1, u1, E1;
    int_val1=0.0;
    double xf=100+c;
    double threea3=3*a*a*a;
    for(int k=-n;k<n;++k)
    {
        u1=( xf*(ts_quad_fixed_nodes[k+n]+1))/2;
        E1=sqrt(u1*u1+b*b);
        int_val1+=ts_quad_fixed_weights[k+n]*(2*u1*(1/a- 1/w)+2*u1*u1*u1/threea3)/(exp(E1-c)-1)*u1/E1;
        
    }
   
    return xf/2*int_val1;



    }
*/
    int n=TS_NUM_PTS;
    double h=TS_STEP_SIZE;
    double w=m*m/(4*p*T);
    double b=m/T;
    double a=w-p/T;

    double int_val1, int_val2, u1, u2, E1, E2;
    int_val1=0.0;
    int_val2=0.0;
double intval_dT_1=0;
double intval_dT_2=0;
double intval_dc_1=0;
double intval_dc_2=0;


    double xf=100+c;

    double exp_term1=0;
    double exp_term2=0;
    double be_term1=0;
    double be_term2=0;
    double log_term1=0;
    double log_term2=0;

    for(int k=-n;k<n;++k)
    {
        u1=( (a-1e-8)*(ts_quad_fixed_nodes[k+n]+1))/2;
        u2=( xf*(ts_quad_fixed_nodes[k+n]+1) - (a+1e-8)*(ts_quad_fixed_nodes[k+n]-1) )/2;
        E1=sqrt(u1*u1+b*b);
        E2=sqrt(u2*u2+b*b);
        exp_term1=exp(E1-c);
        exp_term2=exp(E2-c);
        be_term1=1/(exp_term1-1);
        be_term2=1/(exp_term2-1);
        log_term1=log(fabs((u1+a)/(u1-a)));
        log_term2=log(fabs((u2+a)/(u2-a)));

        int_val1+=ts_quad_fixed_weights[k+n]*(log_term1 - 2*u1/w)*be_term1*u1/E1;
        int_val2+=ts_quad_fixed_weights[k+n]*(log_term2-2*u2/w)*be_term2*u2/E2;

        intval_dc_1+=ts_quad_fixed_weights[k+n]*(log_term1 - 2*u1/w)*be_term1*be_term1*exp_term1*u1/E1;
        intval_dc_2+=ts_quad_fixed_weights[k+n]*(log_term2-2*u2/w)*be_term2*be_term2*exp_term2*u2/E2;

        intval_dT_1+=ts_quad_fixed_weights[k+n]*(log_term1 - 2*u1/w)*be_term1*be_term1*exp_term1*u1;
        intval_dT_2+=ts_quad_fixed_weights[k+n]*(log_term2-2*u2/w)*be_term2*be_term2*exp_term2*u2;
    }
   struct vp_derivs vals;
    vals.vp=T*(int_val1*a/2+int_val2*(xf-a)/2);
    vals.vp_dc=T*(intval_dc_1*a/2+intval_dc_2*(xf-a)/2);
    vals.vp_dT=intval_dT_1*a/2+intval_dT_2*(xf-a)/2;
    return vals;
    //return a/2*int_val1+(xf-a)/2*int_val2;
}


//NOTE: the returns from this function include the factor of T coming from change of variables--this isnt included in the other function
struct va_derivs va_deriv(double m, double T, double p, double c,struct polylogs* plgs)
{
    double w=m*m/(4*p*T);
    /*
    if(w>30)
    {
        return -4*interp_Li4(-exp(c),plgs)/(w*w*w);
    }
    if(w<0.1)
    {
        return 2*interp_Li2(-exp(c),plgs)/w;
    }*/
    int n=TS_NUM_PTS;
    double int_val1,int_val2, u1,u2;
    double xm=w;
    double xf=100+c;
    int_val1=0;
    int_val2=0;

double int_val_dT1=0;
double int_val_dT2=0;
double int_val_dc1=0;
double int_val_dc2=0;

    double exp_term1;
    double exp_term2;
    double be_term1;
    double be_term2;
    double log_term1;
    double log_term2;
    for(int k=-n;k<n;++k)
    {
        u1=( (xm-1e-8)*(ts_quad_fixed_nodes[k+n]+1))/2;
        u2=( xf*(ts_quad_fixed_nodes[k+n]+1) - (xm+1e-8)*(ts_quad_fixed_nodes[k+n]-1) )/2;
        exp_term1=exp(u1-c);
        exp_term2=exp(u2-c);
        be_term1=1/(exp_term1+1);
        be_term2=1/(exp_term2+1);
        log_term1=log(fabs((u1+w)/(u1-w)));
        log_term2=log(fabs((u2+w)/(u2-w)));
        int_val1+=ts_quad_fixed_weights[k+n]*(log_term1-2*u1/w)*be_term1;
        int_val2+=ts_quad_fixed_weights[k+n]*(log_term2-2*u2/w)*be_term2;

        int_val_dc1+=ts_quad_fixed_weights[k+n]*(log_term1-2*u1/w)*be_term1*be_term1*exp_term1;
        int_val_dc2+=ts_quad_fixed_weights[k+n]*(log_term2-2*u2/w)*be_term2*be_term2*exp_term2;

        int_val_dT1+=ts_quad_fixed_weights[k+n]*(log_term1-2*u1/w)*be_term1*be_term1*exp_term1*u1;
        int_val_dT2+=ts_quad_fixed_weights[k+n]*(log_term2-2*u2/w)*be_term2*be_term2*exp_term2*u2;
    }

    struct va_derivs vals;
    vals.va=(int_val1*xm/2+int_val2*(xf-xm)/2)*T;
    vals.va_dc=(int_val_dc1*xm/2+int_val_dc2*(xf-xm)/2)*T;
    vals.va_dT=int_val_dT1*xm/2+int_val_dT2*(xf-xm)/2;
   
    return vals;
}



struct collision_vals_derivs Collision_deriv(double m, double T, double c)
{
    int n=TS_NUM_PTS;
    double int_val, moment_int_val, E;
    int_val=0.0;
    moment_int_val=0;

double b=m/T;

double intval_in=0;
double moment_intval_in=0;
double test_term=0;
double xf=50;
double xi=b;

double zf=70;
double zi=1;

double v;

double E_sqrt_term;
double v_sqrt_term;
double sinh_term_plus;
double sinh_term_minus;
double coth_term_plus;
double coth_term_minus;
double exp_term;
double cs_term;
double fd_term;
double lnsinh_term;

double dT_term=0;
double dc_term=0;

double dT_int_in=0;
double dT_moment_in=0;
double dc_int_in=0;
double dc_moment_in=0;

double intval_dT=0;
double intval_dc=0;
double moment_intval_dT=0;
double moment_intval_dc=0;
#pragma omp taskgroup task_reduction(+:int_val,moment_int_val,intval_dT, moment_intval_dT,intval_dc,moment_intval_dc)
{
#pragma omp taskloop grainsize(10) in_reduction(+:int_val,moment_int_val,intval_dT,moment_intval_dT,intval_dc,moment_intval_dc)
    for(int k=-n;k<n;++k)
    {
        E=( xf*(ts_quad_fixed_nodes[k+n]+1) - xi*(ts_quad_fixed_nodes[k+n]-1) )/2;
        E_sqrt_term=sqrt(E*E-b*b);
        cs_term=c_vv_pp(4*E*E*(T*T),1);
        zf=20/E;
        if(zf<1) continue;
         for(int j=-n;j<n;++j)
        {
        v=( zf*(ts_quad_fixed_nodes[j+n]+1) - zi*(ts_quad_fixed_nodes[j+n]-1) )/2;
        v_sqrt_term=sqrt(v*v-1);
        sinh_term_plus=sinh( (E*v+E_sqrt_term*v_sqrt_term-c)/2 );
        coth_term_plus=coth( (E*v+E_sqrt_term*v_sqrt_term-c)/2 );
        sinh_term_minus=sinh( (E*v-E_sqrt_term*v_sqrt_term-c)/2 );
        coth_term_minus=coth( (E*v-E_sqrt_term*v_sqrt_term-c)/2 );
        exp_term=exp(2*(E*v-c));
        fd_term=1/(exp_term-1);
        lnsinh_term=log(fabs(sinh_term_plus/sinh_term_minus));



        test_term=ts_quad_fixed_weights[j+n]*fd_term*lnsinh_term;
        dT_term=ts_quad_fixed_weights[j+n]*fd_term*(fd_term*exp_term*2*E*v*lnsinh_term-((E*v+E_sqrt_term*v_sqrt_term)*coth_term_plus-(E*v-E_sqrt_term*v_sqrt_term)*coth_term_minus)/2);
        dc_term=ts_quad_fixed_weights[j+n]*fd_term*(2*exp_term*fd_term*lnsinh_term-(coth_term_plus-coth_term_minus)/2);


         if(isnan(test_term) ||isnan(dT_term) || isnan(dc_term))
        {
            break;
        }
            intval_in+=test_term;
            moment_intval_in+=test_term*v;

            dT_int_in+=dT_term;
            dT_moment_in+=dT_term*v;

            dc_int_in+=dc_term;
            dc_moment_in+=dc_term*v;
       
        }
        int_val+=ts_quad_fixed_weights[k+n]*(zf-zi)/2*intval_in*E*E*cs_term;
        moment_int_val+=ts_quad_fixed_weights[k+n]*(zf-zi)/2*moment_intval_in*E*E*E*cs_term;

        intval_dT+=ts_quad_fixed_weights[k+n]*(zf-zi)/2*dT_int_in*E*E*cs_term;
        intval_dc+=ts_quad_fixed_weights[k+n]*(zf-zi)/2*dc_int_in*E*E*cs_term;

        moment_intval_dT+=ts_quad_fixed_weights[k+n]*(zf-zi)/2*dT_moment_in*E*E*E*cs_term;
        moment_intval_dc+=ts_quad_fixed_weights[k+n]*(zf-zi)/2*dc_moment_in*E*E*E*cs_term;

        intval_in=0;
        moment_intval_in=0;
        dT_int_in=0;
        dT_moment_in=0;
        dc_int_in=0;
        dc_moment_in=0;
    }
    }
  struct collision_vals_derivs vals;
  double T4=T*T*T*T;
  vals.coll=(xf-xi)/2*int_val*16*T4/(2*TWO_PI_4);
  vals.moment=(xf-xi)/2*moment_int_val*16*T4*T/(2*TWO_PI_4);
vals.coll_dT=(xf-xi)/2*intval_dT*16*T4/(2*TWO_PI_4*T)+vals.coll/T;
vals.moment_dT=(xf-xi)/2*moment_intval_dT*16*T4/(2*TWO_PI_4)+vals.moment/T;

vals.coll_dc=(xf-xi)/2*intval_dc*16*T*T*T*T/(2*TWO_PI_4);
vals.moment_dc=(xf-xi)/2*moment_intval_dc*16*T*T*T*T*T/(2*TWO_PI_4);

    return vals;
}
