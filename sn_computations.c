#include "sn_computations.h"

double ts_quad_fixed_weights[2*TS_NUM_PTS+1];
double ts_quad_fixed_nodes[2*TS_NUM_PTS+1];

struct collision_vals Moment1(double m, double T, double c)
{
    int n=TS_NUM_PTS;
    double int_val, dint_val, E;
    int_val=0.0;
    dint_val=0;

double b=m/T;

double intval_in=0;
double dintval_in=0;
double test_term=0;
double xf=50;
double xi=b;

double zf=70;
double zi=1;

double v;

    for(int k=-n;k<n;++k)
    {
        E=( xf*(ts_quad_fixed_nodes[k+n]+1) - xi*(ts_quad_fixed_nodes[k+n]-1) )/2;
        zf=20/E;
         for(int j=-n;j<n;++j)
        {
        v=( zf*(ts_quad_fixed_nodes[j+n]+1) - zi*(ts_quad_fixed_nodes[j+n]-1) )/2;
        test_term=ts_quad_fixed_weights[j+n]*E*E/(exp(2*(E*v-c))-1)*log(fabs(sinh( (E*v+sqrt((E*E-b*b)*(v*v-1))-c)/2 )/sinh( (E*v-sqrt((E*E-b*b)*(v*v-1))-c)/2 )))*c_vv_pp(4*E*E*(T*T),1);
         if(isnan(test_term))
        {
            break;
        }
            intval_in+=test_term;
            dintval_in+=test_term*v;
       
        }
        int_val+=ts_quad_fixed_weights[k+n]*(zf-zi)/2*intval_in;
        dint_val+=ts_quad_fixed_weights[k+n]*(zf-zi)/2*dintval_in*E;
        intval_in=0;
        dintval_in=0;
    }
  struct collision_vals vals;
  vals.coll=(xf-xi)/2*int_val*16*T*T*T*T/(2*TWO_PI_4);
  vals.moment=(xf-xi)/2*dint_val*16*T*T*T*T*T/(2*TWO_PI_4);
    return vals;
}

struct collision_vals Moment1_omp(double m, double T, double c)
{
    int n=TS_NUM_PTS;
    double int_val, dint_val, E;
    int_val=0.0;
    dint_val=0;

double b=m/T;

double intval_in=0;
double dintval_in=0;
double test_term=0;
double xf=50;
double xi=b;

double zf=70;
double zi=1;

double v;
#pragma omp taskgroup task_reduction(+:int_val,dint_val)
{
#pragma omp taskloop grainsize(10) in_reduction(+:int_val,dint_val)
    for(int k=-n;k<n;++k)
    {
        E=( xf*(ts_quad_fixed_nodes[k+n]+1) - xi*(ts_quad_fixed_nodes[k+n]-1) )/2;
        zf=20/E;
        if(zf<1) continue;
         for(int j=-n;j<n;++j)
        {
        v=( zf*(ts_quad_fixed_nodes[j+n]+1) - zi*(ts_quad_fixed_nodes[j+n]-1) )/2;
        test_term=ts_quad_fixed_weights[j+n]/(exp(2*(E*v-c))-1)*log(fabs(sinh( (E*v+sqrt((E*E-b*b)*(v*v-1))-c)/2 )/sinh( (E*v-sqrt((E*E-b*b)*(v*v-1))-c)/2 )))*c_pp_vv(4*E*E*(T*T),1);
         if(isnan(test_term))
        {
            break;
        }
            intval_in+=test_term;
            dintval_in+=test_term*v;
       
        }
        int_val+=ts_quad_fixed_weights[k+n]*(zf-zi)/2*intval_in*E*E;
        dint_val+=ts_quad_fixed_weights[k+n]*(zf-zi)/2*dintval_in*E*E*E;
        intval_in=0;
        dintval_in=0;
    }
    }
  struct collision_vals vals;
  vals.coll=(xf-xi)/2*int_val*16*T*T*T*T/(2*TWO_PI_4);
  vals.moment=(xf-xi)/2*dint_val*16*T*T*T*T*T/(2*TWO_PI_4);
    return vals;
}

double Moment2(double m, double T, double c)
{
    int n=TS_NUM_PTS;
    double  u1,u2,u3;
double intval_in=0;
double sigint = 0;
double int_val = 0;

double b = m / T;

double xf = 30;
double xi = 0;

double zf = 30;
double zi = 0;

double ti = -1;
double tf = 1;

double f1, f2;
double E1, E2;

for (int k = -n; k < n; ++k)
{
    intval_in = 0;

    u1 = (xf * (ts_quad_fixed_nodes[k + n] + 1) - xi * (ts_quad_fixed_nodes[k + n] - 1)) / 2;
    E1 = sqrt(u1 * u1 + b * b);
    f1 = 1 / (exp(E1-c) - 1);
    for (int j = -n; j < n; ++j)
    {
    sigint = 0;
        u2 = (zf * (ts_quad_fixed_nodes[j + n] + 1) - zi * (ts_quad_fixed_nodes[j + n] - 1)) / 2;
        E2 = sqrt(u2 * u2 + b * b);
        f2 = 1 / (exp(E2-c) - 1);
        //tf=2*m*m + T*T*2 * (E1*E2+u1 * u2 );
        //ti=2*m*m + T*T*2 * (E1*E2-u1 * u2 );
        for (int l = -n; l < n; ++l)
        {
            u3 = (tf * (ts_quad_fixed_nodes[l + n] + 1) - ti * (ts_quad_fixed_nodes[l + n] - 1)) / 2;
          // sigint += ts_quad_fixed_weights[l + n] * f(u3);
          sigint += ts_quad_fixed_weights[l + n] * c_vv_pp(2*m*m + T*T*2 * (E1 * E2-u1 * u2 * u3),1);
        }
        intval_in += ts_quad_fixed_weights[j + n] * u2  *u2* f2 /(E2) * sigint;
    }
    int_val += ts_quad_fixed_weights[k + n] * u1  *u1* f1 / E1 * intval_in;
}

return (tf - ti) / 2 * (zf - zi) / 2 * (xf - xi) / 2 * int_val*T*T*T*T/(2*TWO_PI_2);
}







double vp(double a,double b,double c,double w)
{
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

    int n=TS_NUM_PTS;
    double h=TS_STEP_SIZE;

    double int_val1, int_val2, u1, u2, E1, E2;
    int_val1=0.0;
    int_val2=0.0;
    double xf=100+c;
    for(int k=-n;k<n;++k)
    {
        u1=( (a-1e-8)*(ts_quad_fixed_nodes[k+n]+1))/2;
        u2=( xf*(ts_quad_fixed_nodes[k+n]+1) - (a+1e-8)*(ts_quad_fixed_nodes[k+n]-1) )/2;
        E1=sqrt(u1*u1+b*b);
        E2=sqrt(u2*u2+b*b);
        int_val1+=ts_quad_fixed_weights[k+n]*(log(fabs((u1+a)/(u1-a))) - 2*u1/w)/(exp(E1-c)-1)*u1/E1;
        int_val2+=ts_quad_fixed_weights[k+n]*(log(fabs((u2+a)/(u2-a)))-2*u2/w)/(exp(E2-c)-1)*u2/E2;
        
    }
   
    return a/2*int_val1+(xf-a)/2*int_val2;
}

double va(double w, double c,struct polylogs* plgs)
{
    if(w>30)
    {
        return -4*interp_Li4(-exp(c),plgs)/(w*w*w);
    }
    if(w<0.1)
    {
        return 2*interp_Li2(-exp(c),plgs)/w;
    }
    int n=TS_NUM_PTS;
    double int_val1,int_val2, u1,u2,div1,div2;
    double xm=w;
    double xf=100+c;
    int_val1=0;
    int_val2=0;
    for(int k=-n;k<n;++k)
    {
        u1=( (xm-1e-8)*(ts_quad_fixed_nodes[k+n]+1))/2;
        u2=( xf*(ts_quad_fixed_nodes[k+n]+1) - (xm+1e-8)*(ts_quad_fixed_nodes[k+n]-1) )/2;
        div1=(fabs((u1+w)/(u1-w)));
        div2=(fabs((u2+w)/(u2-w)));
        int_val1+=ts_quad_fixed_weights[k+n]*((log(div1)-2*u1/w)/(exp(u1-c)+1));
        int_val2+=ts_quad_fixed_weights[k+n]*(log(div2)-2*u2/w)/(exp(u2-c)+1);
    }
   
    return int_val1*xm/2+int_val2*(xf-xm)/2;
}


double v_tot(double p, double m_a, double m_phi, double y, struct sp_params* sp_params,struct polylogs* plgs)
{
    double T=sp_params->T;
    double Ts=sp_params->Ts;
    double Tphi=sp_params->Tphi;
    double c_s=sp_params->c_s;
    double c_phi=sp_params->c_phi;

    double p_l_a=n_rho_fermi(m_a,T); //energy density of the lepton corresponding to v_a
    double p_a=7/240*M_PI*T*T*T*T; //energy density of v_a

    double v_SM =-8*SQRT2*FERMI_CONST/3*p*(2*p_l_a/M_W_BOSON+2*p_a/M_Z_BOSON); //CHECK IF I NEED 2 HERE AND G=2 ABOVE
    double mphi2=m_phi*m_phi;
    double wbase=mphi2/(4*p);

    double v_a=va(wbase/T, 0,plgs)*mphi2/2*T;
    double v_s=va(wbase/Ts,c_s,plgs)*mphi2/2*Ts;
    double v_phi=vp(wbase/Tphi-p/Tphi,m_phi/Tphi,c_phi,wbase/Tphi)*mphi2/2*Tphi;

    return v_SM-y*y/(4*TWO_PI_2*p*p)*(v_a+v_s+v_phi);
}



//This is only the integral part--need to multiply these by y^4*m^2/((2pi)^2*4*p1^2*w^2)*T
struct massless_rates compute_massless_rates(double T,double Ts,double w,double w_s,double c_s,double y,struct polylogs* plgs)
{

    int n=TS_NUM_PTS;
    double int_val, weight, u;
    int_val = 0.0;
    double xf=100;
    double int_base_a;
    double int_base_s;

    double g_aa_ss,g_ss_aa, g_as_as, g_sa_as, g_aa_pp, g_ss_pp, g_pv_pv;
    g_aa_pp=0;
    g_ss_pp=0;
    g_as_as=0;
    g_sa_as=0;
    g_aa_ss=0;
    g_ss_aa=0;

double expterm=0;
double expterm_cs=0;
    for (int k = -n; k < n; ++k)
    {
        u = (xf * (ts_quad_fixed_nodes[k + n] + 1) ) / 2;
        weight=ts_quad_fixed_weights[k+n];
        expterm=exp(-u);
        expterm_cs=exp(c_s-u);
        int_base_a=log(1+expterm)*u;
        int_base_s=log(1+expterm_cs)*u;

        g_aa_ss+=weight*int_base_a*c_vv_vv(u/w,y);
        g_ss_aa+=weight*int_base_s*c_vv_vv(u/w_s,y);

        g_as_as+=weight*int_base_s*c_as_as(u/w_s,y);
        g_sa_as+=weight*int_base_a*c_as_as(u/w,y);

        g_aa_pp+=weight*int_base_a*c_vv_pp(u/w,y);
        g_ss_pp+=weight*int_base_s*c_vv_pp(u/w_s,y);

    }

    struct massless_rates rates;
    rates.g_aa_ss=g_aa_ss*xf/2;
    rates.g_ss_aa=g_ss_aa*xf/2;
    rates.g_as_as=g_as_as*xf/2;
    rates.g_sa_as=g_sa_as*xf/2;
    rates.g_aa_pp=g_aa_pp*xf/2;
    rates.g_ss_pp=g_ss_pp*xf/2;
    return rates;
}

void free_massive_table(struct massive_interpolation_table* mip)
{
    free(mip->table);
    free(mip->d_table);
}

double compute_vp_vp_rate(double mphi,double p1,double cphi,double Tphi,double y,struct massive_interpolation_table* mip)
{
    double b=mphi/Tphi;
    if(b>25)
    {
        return 0;
    }
    int n=TS_NUM_PTS;
    double int_val, weight, node, u2 , E2,frac1,frac2;
    double wphi=mphi*mphi/(4*p1*Tphi);
    int_val=0.0;
    double xi=0;
    double xf=100;
    double lb,ub;
    double mphi2=mphi*mphi;
    for(int k=-n;k<n;++k)
    {
        u2=( xf*(ts_quad_fixed_nodes[k+n]+1) - xi*(ts_quad_fixed_nodes[k+n]-1) )/2;
        E2=sqrt(u2*u2+b*b);
        lb=1+(E2-u2)/(2*wphi);
        ub=1+(E2+u2)/(2*wphi);
        int_val+=ts_quad_fixed_weights[k+n]*interp_massive(lb,ub,y,mip)*u2/(E2*(exp(E2-cphi)-1));
        
    }
    return (xf-xi)/2*int_val;
}


double C_as_as_narrow(double p1, double c, double m, double T, double y)
{
   return 0.159154943092*y*y*T*m*m/(p1*p1)*log(1+exp(c-m*m/(4*T*p1)));
}


double get_massive_rate(double b, double mphi, struct massive_interpolation_table* massive_table, double y, double a3, double b3, double a4, double b4)
{
    int n=TS_NUM_PTS;
    double int_val, dint_val, moment_int_val, d_moment_int_val, weight, node, u3, u4, lb, ub, E3, E4, be3, be4, tempintval,tempdintval;
    int_val=0.0;
    dint_val=0;
    moment_int_val=0;
    d_moment_int_val=0;
    lb=0;
    ub=0;
    E3=0;
    E4=0;
    be3=0;
    be4=0;
    //b=mphi/Tphi;
    for(int i=-n;i<n;++i)
    {
        for(int j=-n;j<n;++j)
        {
            u3=( b3*(ts_quad_fixed_nodes[i+n]+1) - a3*(ts_quad_fixed_nodes[i+n]-1) )/2;
            u4=( b4*(ts_quad_fixed_nodes[j+n]+1) - a4*(ts_quad_fixed_nodes[j+n]-1) )/2;
            E3=sqrt(u3*u3+b*b);
            E4=sqrt(u4*u4+b*b);
            be3=1/(exp(E3)-1);
            be4=1/(exp(E4)-1);
            lb=2*(1+(E3*E4 - u3*u4 )/(b*b));
            ub=2*(1+(E3*E4 + u3*u4 )/(b*b));
            
            tempintval=ts_quad_fixed_weights[i+n]*ts_quad_fixed_weights[j+n]*interp_massive(lb, ub, y, massive_table)*u3*u4/(E3*E4)*be3*be4;
            int_val+=tempintval;
        }
    }
   // double extfac=mphi*mphi/(4*M_PI*M_PI);
    return (b4-a4)/2*(b3-a3)/2*int_val;// *extfac*extfac*Tphi*Tphi;
}

double interp_massive(double i0,double i1,double y, struct massive_interpolation_table* mip)
{
    double x0=mip->x0;
    double x1=mip->x1;
    if(i1>1000)
    {
        i1=1000;
    }
    int idx0=floor((i0-x0)/mip->h);
    int idx1=floor((i1-x0)/mip->h);
    int N=mip->N;
    if(i0<0 || i1<0 ||idx0>=N || idx1>=N || idx1<0 ||idx0<0)
    {
        return ts_quad(&c_vp_vp, &y, i0, i1, TS_STEP_SIZE, TS_NUM_PTS);
    }
    
    return two_point_hermite_interp(i1, x0+mip->h*idx1, x0+mip->h*(idx1+1), mip->table[idx1], mip->d_table[idx1], mip->table[idx1+1], mip->d_table[idx1+1])
    - two_point_hermite_interp(i0, x0+mip->h*idx0, x0+mip->h*(idx0+1), mip->table[idx0], mip->d_table[idx0], mip->table[idx0+1], mip->d_table[idx0+1]);
}

void fill_massive_tables(struct massive_interpolation_table* table, double y)
{
    double x0=table->x0;
    double x1=table->x1;
    double h=table->h;
    int N=ceil((x1-x0)/h);
    table->N=N;
    
    table->x1=x0+h*N;
    x1=table->x1;
    
    table->table=(double*) malloc(N*sizeof(double));
    table->d_table=(double*) malloc(N*sizeof(double));
    

        #pragma omp parallel for
        for(int i=0;i<N;++i)
        {
            table->table[i]=ts_quad(&c_vp_vp, &y, 4, x0+i*h, TS_STEP_SIZE, TS_NUM_PTS);
            table->d_table[i]=c_vp_vp(x0+i*h,&y);
        }

    
    
}



double oscillation_peak_func(double p, double m_s, double cos_2theta, double m_a, double m_phi, double y, struct sp_params* sp_params,struct polylogs* plgs)
{
    return m_s*m_s/(2*p)*cos_2theta-v_tot(p,m_a,m_phi,y,sp_params,plgs);

}



double find_oscillation_peak(double m_a, double m_s, double cos_2theta, double m_phi, double y, struct sp_params* sp_params,struct polylogs* plgs)
{
    double T=sp_params->T;
    double Ts=sp_params->Ts;
    double Tphi=sp_params->Tphi;
    double c_s=sp_params->c_s;
    double c_phi=sp_params->c_phi;
    
    double x=0.1*m_phi*m_phi/(4*T);
    double func_eval=oscillation_peak_func(x,m_s,cos_2theta,m_a,m_phi,y,sp_params,plgs);

    int num_its=0;
    while (func_eval > 0)
    {
        ++num_its;
        x *= 10;
        func_eval = oscillation_peak_func(x, m_s, cos_2theta, m_a, m_phi, y, sp_params,plgs);
        if (num_its > 20)
        {
            printf("No oscillation peak found\n");
            return -1;
        }
    }

    double xL=0;
    double xM=0;
    if(num_its>1)
    {
        xL=x-1;
    }
    num_its=0;
    while(fabs(func_eval)>1e-10)
    {
        ++num_its;
        xM = (xL + x) / 2;
        func_eval = oscillation_peak_func(xM, m_s, cos_2theta, m_a, m_phi, y, sp_params, plgs);
        if (func_eval > 0)
        {
            xL = xM;
        }
        else
        {
            x = xM;
        }
        if(num_its>50)
        {
            printf("No oscillation peak found\n");
            return -1;
        }
    }
    return xM;
}



double oscillation_term(double p, struct fixed_params fixed_params,struct sp_params sp_params, struct polylogs* plgs,struct massive_interpolation_table* mip)
{
    double m_phi=fixed_params.m_phi;
    double m_s=fixed_params.m_s;
    double y=fixed_params.y;
    double theta=fixed_params.theta;
    double T=sp_params.T;
    double Ts=sp_params.Ts;
    double Tphi=sp_params.Tphi;
    double c_s=sp_params.c_s;
    double c_phi=sp_params.c_phi;
    double cos_2theta=cos(2*theta);
    double sin2_2theta=sin(2*theta)*sin(2*theta);
    double m_a=fixed_params.m_l_a;

    double w = m_phi * m_phi / (4 * p * T);
    double w_s = m_phi * m_phi / (4 * p * Ts);
    double w_phi = m_phi * m_phi / (4 * p * Tphi);
    double del = m_s * m_s / (2 * p);
    double y2 = y * y;
    double mphi2 = m_phi * m_phi;
    double y4 = y2 * y2;

    double v_tot_val = v_tot(p, m_a, m_phi, y, &sp_params, plgs);
    struct massless_rates rates = compute_massless_rates(T, Ts, w, w_s, c_s, y, plgs);

    double massless_coeff = y4 * mphi2 / (TWO_PI_2 * 4 * p * p);
    double g_aa_ss = rates.g_aa_ss*massless_coeff*T/(w*w);
    double g_ss_aa = rates.g_ss_aa*massless_coeff*Ts/(w_s*w_s);
    double g_as_as = C_as_as_narrow(p, c_s, m_phi, Ts, y*massless_coeff*Ts/(w_s*w_s));
    double g_sa_as = C_as_as_narrow(p, 0, m_phi, T, y)*massless_coeff*T/(w*w);
    double g_aa_pp = rates.g_aa_pp*massless_coeff*T/(w*w);
    double g_ss_pp = rates.g_ss_pp*massless_coeff*Ts/(w_s*w_s);

    double M_as_p = 8 * y2 * (mphi2 - m_s * m_s);
    double g_sa_p = M_as_p / (32 * M_PI * p * p) * (T * log(1 + exp(-w)));   // check
    double g_as_p = M_as_p / (32 * M_PI * p * p) * (Ts * log(1 + exp(-w_s))); // check
    double g_vp_vp = compute_vp_vp_rate(m_phi, p, c_phi, Tphi, y, mip) *Tphi*m_phi/(8*TWO_PI_2*p*p);

    double gamma_tot = g_aa_ss + g_aa_pp + g_as_p + g_as_as + g_ss_aa + g_ss_pp + g_sa_p + g_sa_as + 2 * g_vp_vp;
    double denom = del * del * sin2_2theta + gamma_tot * gamma_tot / 4 + (del * cos_2theta - v_tot_val) * (del * cos_2theta - v_tot_val);
    return gamma_tot / 4 * del * del * sin2_2theta / denom;
}



struct derivs integrated_de(double xi, double xf, struct fixed_params fixed_params, struct sp_params sp_params, struct polylogs* plgs,struct massive_interpolation_table* mip)
{
int n=TS_NUM_PTS;
double m_phi=fixed_params.m_phi;
double m_s=fixed_params.m_s;
double y=fixed_params.y;
double theta=fixed_params.theta;
double T=sp_params.T;
double Ts=sp_params.Ts;
double Tphi=sp_params.Tphi;
double c_s=sp_params.c_s;
double c_phi=sp_params.c_phi;
double cos_2theta=cos(2*theta);
double sin2_2theta=sin(2*theta)*sin(2*theta);
double m_a=fixed_params.m_l_a;

double T3=T*T*T;
double Ts3=Ts*Ts*Ts;
double Tphi3=Tphi*Tphi*Tphi;


double p1;

double vt_tot=0;
double g_aa_ss=0;
double g_ss_aa=0;
double g_as_as=0;
double g_sa_as=0;
double g_aa_pp=0;
double g_ss_pp=0;
double g_vp_vp=0;
double g_as_p=0;
double g_sa_p=0;
double g_tot=0;
double osc_term=0;
double massless_coeff=0;
double massless_moment_coeff=0;
double massless_coll_coeff=0;

double C_p_as_create=0;
double C_p_as_destroy=0;
double C_as_p_create=0;
double C_as_p_destroy=0;
double C_p_as_exp_term_c=0;
double C_p_as_exp_term_d=0;
double C_as_p_num=0;

double C_aa_ss=0;
double C_ss_aa=0;
double C_aa_pp=0;
double C_ss_pp=0;

double C_aa_ss_moment=0;
double C_ss_aa_moment=0;
double C_aa_pp_moment=0;
double C_ss_pp_moment=0;


struct massless_rates rates;

double m_phi2=m_phi*m_phi;
double y4=y*y*y*y;
double M_as_p = 8 * y*y * (m_phi2 - m_s * m_s);
double del;
double osc_denom;
double integrated_osc=0;
double moment_osc=0;
double weight;

double w;
double ws;
double wphi;
double f_s;
double f_fd;


for (int k = -n; k < n; ++k)
    {
        p1 = (xf * (ts_quad_fixed_nodes[k + n] + 1) - xi * (ts_quad_fixed_nodes[k + n] - 1)) / 2;
        weight = ts_quad_fixed_weights[k + n];
        w = m_phi2 / (4 * p1 * T);
        ws = m_phi2 / (4 * p1 * Ts);
        wphi = m_phi2 / (4 * p1 * Tphi);
        f_s=1/(exp(p1-c_s)+1);
        f_fd=1/(exp(p1)+1);

        //This is the integral of  f_s*f_a*p1--it is the CONSUMPTION of rho_s via the process a+s->p
        C_as_p_destroy+=weight * p1*f_s*log(1+exp(-w/Ts));//T*T_s*T_s
        //This is the integral of  f_p*p1--it is the PRODUCTION of rho_s via the process p->a+s
        C_as_p_create+=-weight *  p1*log(1-exp(c_phi-wphi/Tphi-p1));//T_phi^3

        //This is the integral of f_s*f_a*(p1+p2)--this is the PRODUCTION of rho_p via the process a+s->p
        C_p_as_exp_term_c=exp(c_s-ws/T);
        C_p_as_create+=weight*f_fd*( (T*p1+Ts/T *ws)*log(1+C_p_as_exp_term_c)-Ts*interp_Li2(-C_p_as_exp_term_c,plgs) ); //Multiply bs T_s*T
        //C_p_as_moment_R+=weight*interp_Li2(-exp(-p1),plgs);
        //This is the integral of f_phi*(p1+p2)--this is the CONSUMPTION of rho_p via the process p->a+s
        C_p_as_exp_term_d=exp(c_phi-wphi/Tphi-p1);
        C_p_as_destroy+=weight*(interp_Li2(C_p_as_exp_term_d,plgs)-(wphi/Tphi+p1)*log(1-C_p_as_exp_term_d)); //Multiply by T_phi^3

 //This is the collision term for the process a+s->p. It represents the production of s via this this process, not the consumption
        C_as_p_num-=weight*(Tphi*Tphi*log(1-exp(c_phi-wphi/Tphi-p1))+Ts*T*f_fd*log(1+exp(c_s-ws/T))); //Done!

    //now compute all other collision terms. We use the substitution p1->T*p1
        
        w=w/T;
        ws=ws/T;
        wphi=wphi/T;
        f_s=1/(exp(p1*T/Ts-c_s)+1);
        del = m_s * m_s / (2 * p1*T);

        rates = compute_massless_rates(T, Ts, w, ws, c_s, y, plgs);
        //massless_coeff = y4*m_phi2 / (TWO_PI_2 * 4 * p1 * p1);
        massless_coeff = 4*y4 / (TWO_PI_2 *m_phi2);
        massless_coll_coeff=weight*p1*p1*T3*2/TWO_PI_2;// This is meant to absorb the extra 4*pi*p1^2/(2*TWO_PI_3) term that comes from the outermost integral
        massless_moment_coeff =  massless_coll_coeff*p1*T; 
        g_aa_ss = rates.g_aa_ss * massless_coeff * T3;
        g_ss_aa = rates.g_ss_aa * massless_coeff * Ts3;

        g_as_as=C_as_as_narrow(p1,c_s,m_phi,Ts,y);
        g_sa_as=C_as_as_narrow(p1,0,m_phi,T,y);
        g_aa_pp = rates.g_aa_pp * massless_coeff *T3;
        g_ss_pp = rates.g_ss_pp * massless_coeff * Ts3;

        g_vp_vp=compute_vp_vp_rate(m_phi,p1,c_phi,Tphi,y,mip)*y4*Tphi*m_phi/(8*TWO_PI_2*p1*p1);


        g_sa_p = M_as_p / (32 * M_PI * p1 * p1*T*T) * (T * log(1 + exp(-w)));   // check
        g_as_p = M_as_p / (32 * M_PI * p1 * p1*T*T) * (Ts * log(1 + exp(-ws))); // check

        g_tot = g_aa_ss + g_aa_pp + g_as_p + g_as_as + g_ss_aa + g_ss_pp + g_sa_p + g_sa_as + 2 * g_vp_vp;


        vt_tot=v_tot(p1,m_a,m_phi,y,&sp_params,plgs);

        osc_denom=del*del*sin2_2theta+g_tot*g_tot/4+(del*cos_2theta-vt_tot)*(del*cos_2theta-vt_tot);
        osc_term=del*del*sin2_2theta/osc_denom*g_tot/4*(f_fd-f_s);
        integrated_osc+=weight*osc_term*p1*T*T/TWO_PI_2;
        moment_osc+=weight*osc_term*p1*p1*T3/TWO_PI_2;

        C_aa_pp+=g_aa_pp*massless_coll_coeff*f_fd;
        C_ss_pp+=g_ss_pp*massless_coll_coeff*f_s;
        C_aa_ss+=g_aa_ss*massless_coll_coeff*f_fd;
        C_ss_aa+=g_ss_aa*massless_coll_coeff*f_s;

        C_aa_pp_moment+=g_aa_pp*massless_moment_coeff*f_fd;
        C_ss_pp_moment+=g_ss_pp*massless_moment_coeff*f_s;
        C_aa_ss_moment+=g_aa_ss*massless_moment_coeff*f_fd;
        C_ss_aa_moment+=g_ss_aa*massless_moment_coeff*f_s;




    }
    struct collision_vals C_pp_aa_vals=Moment1(m_phi,T,0);
    struct collision_vals C_pp_ss_vals=Moment1(m_phi,Ts,c_s);
    double bound_weight = (xf - xi) / 2;

    double C_as_p_moment_tot=bound_weight*(C_as_p_create*Tphi*Tphi*Tphi-C_as_p_destroy*T*Ts*Ts)*M_as_p/(8*TWO_PI_3);
    double C_p_as_moment_tot=bound_weight*(C_p_as_create*Ts*T-C_p_as_destroy*Tphi*Tphi*Tphi)*M_as_p/(8*TWO_PI_3);

    C_as_p_num=bound_weight*C_as_p_num*M_as_p/(8*TWO_PI_3);

    struct derivs derivs;

    derivs.C_n_s=bound_weight*(integrated_osc+C_aa_ss-C_ss_aa)+C_as_p_num;
    derivs.C_p_s=bound_weight*(moment_osc+C_aa_ss_moment-C_ss_aa_moment)+C_as_p_moment_tot;

    derivs.C_n_p=bound_weight*(C_aa_pp+C_ss_pp)-y4*C_pp_aa_vals.coll-y4*C_pp_ss_vals.coll-C_as_p_num;
    derivs.C_p_p=bound_weight*(C_aa_pp_moment+C_ss_pp_moment)-y4*C_pp_aa_vals.moment-y4*C_pp_ss_vals.moment+C_p_as_moment_tot;
/*
printf("C_aa_ss: %.10e\n",C_aa_ss*bound_weight);
printf("C_ss_aa: %.10e\n",C_ss_aa*bound_weight);
printf("C_aa_pp: %.10e\n",C_aa_pp*bound_weight);
printf("C_ss_pp: %.10e\n",C_ss_pp*bound_weight);
printf("C_as_p_num: %.10e\n",C_as_p_num);
printf("C_as_p_moment_tot: %.10e\n",C_as_p_moment_tot);
printf("C_p_as_moment_tot: %.10e\n",C_p_as_moment_tot);
printf("C_pp_aa_vals.coll: %.10e\n",y4*C_pp_aa_vals.coll);
printf("C_pp_ss_vals.coll: %.10e\n",y4*C_pp_ss_vals.coll);
printf("C_pp_aa_vals.moment: %.10e\n",y4*C_pp_aa_vals.moment);
printf("C_pp_ss_vals.moment: %.10e\n",y4*C_pp_ss_vals.moment);
printf("C_aa_ss_moment: %.10e\n",C_aa_ss_moment*bound_weight);
printf("C_ss_aa_moment: %.10e\n",C_ss_aa_moment*bound_weight);
printf("C_aa_pp_moment: %.10e\n",C_aa_pp_moment*bound_weight);
printf("C_ss_pp_moment: %.10e\n",C_ss_pp_moment*bound_weight);
printf("integrated_osc: %.10e\n",integrated_osc);
printf("moment_osc: %.10e\n",moment_osc);

printf("C_n_s: %.10e\n",derivs.C_n_s);
printf("C_p_s: %.10e\n",derivs.C_p_s);
printf("C_n_p: %.10e\n",derivs.C_n_p);
printf("C_p_p: %.10e\n",derivs.C_p_p);*/
    return derivs;
}


double H(double T, struct fixed_params* fp)
{
    
    double geff=interpolate(fp->T_g, fp->g_s, fp->len_g, T);
    return sqrt((8*M_PI*M_PI*M_PI*geff/90))*(T*T)/M_PL;
}


struct derivs integrated_de_omp(double xi, double xf, struct fixed_params fixed_params, struct sp_params sp_params, struct polylogs* plgs,struct massive_interpolation_table* mip)
{

int n=TS_NUM_PTS;
double m_phi=fixed_params.m_phi;
double m_s=fixed_params.m_s;
double y=fixed_params.y;
double theta=fixed_params.theta;
double T=sp_params.T;
double Ts=sp_params.Ts;
double Tphi=sp_params.Tphi;
double c_s=sp_params.c_s;
double c_phi=sp_params.c_phi;
double cos_2theta=cos(2*theta);
double sin2_2theta=sin(2*theta)*sin(2*theta);
double m_a=fixed_params.m_l_a;

double T3=T*T*T;
double Ts3=Ts*Ts*Ts;
double Tphi3=Tphi*Tphi*Tphi;


double p1;

double vt_tot=0;
double g_aa_ss=0;
double g_ss_aa=0;
double g_as_as=0;
double g_sa_as=0;
double g_aa_pp=0;
double g_ss_pp=0;
double g_vp_vp=0;
double g_as_p=0;
double g_sa_p=0;
double g_tot=0;
double osc_term=0;
double massless_coeff=0;
double massless_moment_coeff=0;
double massless_coll_coeff=0;

double C_p_as_create=0;
double C_p_as_destroy=0;
double C_as_p_create=0;
double C_as_p_destroy=0;
double C_p_as_exp_term_c=0;
double C_p_as_exp_term_d=0;
double C_as_p_num=0;

double C_aa_ss=0;
double C_ss_aa=0;
double C_aa_pp=0;
double C_ss_pp=0;

double C_aa_ss_moment=0;
double C_ss_aa_moment=0;
double C_aa_pp_moment=0;
double C_ss_pp_moment=0;


struct massless_rates rates;

double m_phi2=m_phi*m_phi;
double y4=y*y*y*y;
double M_as_p = 8 * y*y * (m_phi2 - m_s * m_s);
double del;
double osc_denom;
double integrated_osc=0;
double moment_osc=0;
double weight;

double w;
double ws;
double wphi;
double f_s;
double f_fd;

#pragma omp taskgroup task_reduction(+:C_as_p_destroy) task_reduction(+:C_as_p_create) task_reduction(+:C_p_as_create) task_reduction(+:C_p_as_destroy) task_reduction(+:C_as_p_num) task_reduction(+:integrated_osc) task_reduction(+:moment_osc) task_reduction(+:C_aa_ss) task_reduction(+:C_ss_aa) task_reduction(+:C_aa_pp) task_reduction(+:C_ss_pp) task_reduction(+:C_aa_ss_moment) task_reduction(+:C_ss_aa_moment) task_reduction(+:C_aa_pp_moment) task_reduction(+:C_ss_pp_moment)
{ 
    #pragma omp taskloop grainsize(2) in_reduction(+:C_as_p_destroy) in_reduction(+:C_as_p_create) in_reduction(+:C_p_as_create) in_reduction(+:C_p_as_destroy) in_reduction(+:C_as_p_num) in_reduction(+:integrated_osc) in_reduction(+:moment_osc) in_reduction(+:C_aa_ss) in_reduction(+:C_ss_aa) in_reduction(+:C_aa_pp) in_reduction(+:C_ss_pp) in_reduction(+:C_aa_ss_moment) in_reduction(+:C_ss_aa_moment) in_reduction(+:C_aa_pp_moment) in_reduction(+:C_ss_pp_moment)
for (int k = -n; k < n; ++k)
    {
        p1 = (xf * (ts_quad_fixed_nodes[k + n] + 1) - xi * (ts_quad_fixed_nodes[k + n] - 1)) / 2;
        weight = ts_quad_fixed_weights[k + n];
        w = m_phi2 / (4 * p1 * T);
        ws = m_phi2 / (4 * p1 * Ts);
        wphi = m_phi2 / (4 * p1 * Tphi);
        f_s=1/(exp(p1-c_s)+1);
        f_fd=1/(exp(p1)+1);

        //This is the integral of  f_s*f_a*p1--it is the CONSUMPTION of rho_s via the process a+s->p
        C_as_p_destroy+=weight * p1*f_s*log(1+exp(-w/Ts));//T*T_s*T_s
        //This is the integral of  f_p*p1--it is the PRODUCTION of rho_s via the process p->a+s
        C_as_p_create+=-weight *  p1*log(1-exp(c_phi-wphi/Tphi-p1));//T_phi^3

        //This is the integral of f_s*f_a*(p1+p2)--this is the PRODUCTION of rho_p via the process a+s->p
        C_p_as_exp_term_c=exp(c_s-ws/T);
        C_p_as_create+=weight*f_fd*( (T*p1+Ts/T *ws)*log(1+C_p_as_exp_term_c)-Ts*interp_Li2(-C_p_as_exp_term_c,plgs) ); //Multiply bs T_s*T
        //C_p_as_moment_R+=weight*interp_Li2(-exp(-p1),plgs);
        //This is the integral of f_phi*(p1+p2)--this is the CONSUMPTION of rho_p via the process p->a+s
        C_p_as_exp_term_d=exp(c_phi-wphi/Tphi-p1);
        C_p_as_destroy+=weight*(interp_Li2(C_p_as_exp_term_d,plgs)-(wphi/Tphi+p1)*log(1-C_p_as_exp_term_d)); //Multiply by T_phi^3

 //This is the collision term for the process a+s->p. It represents the production of s via this this process, not the consumption
        C_as_p_num-=weight*(Tphi*Tphi*log(1-exp(c_phi-wphi/Tphi-p1))+Ts*T*f_fd*log(1+exp(c_s-ws/T))); //Done!

    //now compute all other collision terms. We use the substitution p1->T*p1
        
        w=w/T;
        ws=ws/T;
        wphi=wphi/T;
        f_s=1/(exp(p1*T/Ts-c_s)+1);
        del = m_s * m_s / (2 * p1*T);
        p1=p1*T;

        rates = compute_massless_rates(T, Ts, w, ws, c_s, y, plgs);
        //massless_coeff = y4*m_phi2 / (TWO_PI_2 * 4 * p1 * p1);
        massless_coeff = 4*y4 / (TWO_PI_2 *m_phi2);
        massless_coll_coeff=weight*p1*p1*T3*2/TWO_PI_2;// This is meant to absorb the extra 4*pi*p1^2/(2*TWO_PI_3) term that comes from the outermost integral
        massless_moment_coeff =  massless_coll_coeff*p1*T; 
        g_aa_ss = rates.g_aa_ss * massless_coeff * T3;
        g_ss_aa = rates.g_ss_aa * massless_coeff * Ts3;

        g_as_as=C_as_as_narrow(p1,c_s,m_phi,Ts,y);
        g_sa_as=C_as_as_narrow(p1,0,m_phi,T,y);
        g_aa_pp = rates.g_aa_pp * massless_coeff *T3;
        g_ss_pp = rates.g_ss_pp * massless_coeff * Ts3;

        g_vp_vp=compute_vp_vp_rate(m_phi,p1,c_phi,Tphi,y,mip)*y4*Tphi*m_phi/(8*TWO_PI_2*p1*p1);


        g_sa_p = M_as_p / (32 * M_PI * p1 * p1*T*T) * (T * log(1 + exp(-w)));   // check
        g_as_p = M_as_p / (32 * M_PI * p1 * p1*T*T) * (Ts * log(1 + exp(-ws))); // check

        g_tot = g_aa_ss + g_aa_pp + g_as_p + g_as_as + g_ss_aa + g_ss_pp + g_sa_p + g_sa_as + 2 * g_vp_vp;


        vt_tot=v_tot(p1,m_a,m_phi,y,&sp_params,plgs);

        osc_denom=del*del*sin2_2theta+g_tot*g_tot/4+(del*cos_2theta-vt_tot)*(del*cos_2theta-vt_tot);
        osc_term=del*del*sin2_2theta/osc_denom*g_tot/4*(f_fd-f_s);
        integrated_osc+=weight*osc_term*p1*p1/TWO_PI_2;
        moment_osc+=weight*osc_term*p1*p1*p1/TWO_PI_2;

        C_aa_pp+=g_aa_pp*massless_coll_coeff*f_fd;
        C_ss_pp+=g_ss_pp*massless_coll_coeff*f_s;
        C_aa_ss+=g_aa_ss*massless_coll_coeff*f_fd;
        C_ss_aa+=g_ss_aa*massless_coll_coeff*f_s;

        C_aa_pp_moment+=g_aa_pp*massless_moment_coeff*f_fd;
        C_ss_pp_moment+=g_ss_pp*massless_moment_coeff*f_s;
        C_aa_ss_moment+=g_aa_ss*massless_moment_coeff*f_fd;
        C_ss_aa_moment+=g_ss_aa*massless_moment_coeff*f_s;




    }
}

    struct collision_vals C_pp_aa_vals=Moment1_omp(m_phi,T,0);
    struct collision_vals C_pp_ss_vals=Moment1_omp(m_phi,Ts,c_s);
    double bound_weight = (xf - xi) / 2;
/*printf("C_as_p_destroy: %.10e\n",bound_weight*C_as_p_destroy*T*Ts*Ts);
printf("C_as_p_create: %.10e\n",bound_weight*C_as_p_create*Tphi*Tphi*Tphi);
printf("C_p_as_create: %.10e\n",bound_weight*C_p_as_create*Ts*T);
printf("C_p_as_destroy: %.10e\n",bound_weight*C_p_as_destroy*Tphi*Tphi*Tphi);
printf("C_as_p_num: %.10e\n",bound_weight*C_as_p_num);*/
    double C_as_p_moment_tot=bound_weight*(C_as_p_create*Tphi*Tphi*Tphi-C_as_p_destroy*T*Ts*Ts)*M_as_p/(8*TWO_PI_3);
    double C_p_as_moment_tot=bound_weight*(C_p_as_create*Ts*T-C_p_as_destroy*Tphi*Tphi*Tphi)*M_as_p/(8*TWO_PI_3);
//printf("C_as_p_moment_tot: %.10e\n",C_as_p_moment_tot);
//printf("C_p_as_moment_tot: %.10e\n",C_p_as_moment_tot);
//printf("as->p moment: %.10e\n",bound_weight*(C_as_p_create*Tphi*Tphi*Tphi-C_as_p_destroy*T*Ts*Ts));
//printf("p_>as: %.10e\n",bound_weight*(C_p_as_create*Ts*T-C_p_as_destroy*Tphi*Tphi*Tphi));
    C_as_p_num=bound_weight*C_as_p_num*M_as_p/(8*TWO_PI_3);
//printf("C_as_p_num: %.10e\n",C_as_p_num);
    struct derivs derivs;

    derivs.C_n_s=T*bound_weight*(integrated_osc+C_aa_ss-C_ss_aa)+C_as_p_num;
    derivs.C_p_s=T*bound_weight*(moment_osc+C_aa_ss_moment-C_ss_aa_moment)+C_as_p_moment_tot;

    derivs.C_n_p=T*bound_weight*(C_aa_pp+C_ss_pp)-y4*C_pp_aa_vals.coll-y4*C_pp_ss_vals.coll-C_as_p_num;
    derivs.C_p_p=T*bound_weight*(C_aa_pp_moment+C_ss_pp_moment)-y4*C_pp_aa_vals.moment-y4*C_pp_ss_vals.moment+C_p_as_moment_tot;

printf("C_aa_ss: %.10e\n",C_aa_ss*bound_weight);
printf("C_ss_aa: %.10e\n",C_ss_aa*bound_weight);
printf("C_aa_pp: %.10e\n",C_aa_pp*bound_weight);
printf("C_ss_pp: %.10e\n",C_ss_pp*bound_weight);
printf("C_as_p_num: %.10e\n",C_as_p_num);
printf("C_as_p_moment_tot: %.10e\n",C_as_p_moment_tot);
printf("C_p_as_moment_tot: %.10e\n",C_p_as_moment_tot);
printf("C_pp_aa_vals.coll: %.10e\n",y4*C_pp_aa_vals.coll);
printf("C_pp_ss_vals.coll: %.10e\n",y4*C_pp_ss_vals.coll);
printf("C_pp_aa_vals.moment: %.10e\n",y4*C_pp_aa_vals.moment);
printf("C_pp_ss_vals.moment: %.10e\n",y4*C_pp_ss_vals.moment);
printf("C_aa_ss_moment: %.10e\n",C_aa_ss_moment*bound_weight);
printf("C_ss_aa_moment: %.10e\n",C_ss_aa_moment*bound_weight);
printf("C_aa_pp_moment: %.10e\n",C_aa_pp_moment*bound_weight);
printf("C_ss_pp_moment: %.10e\n",C_ss_pp_moment*bound_weight);
printf("integrated_osc: %.10e\n",integrated_osc);
printf("moment_osc: %.10e\n",moment_osc);

printf("C_n_s: %.10e\n",derivs.C_n_s);
printf("C_p_s: %.10e\n",derivs.C_p_s);
printf("C_n_p: %.10e\n",derivs.C_n_p);
printf("C_p_p: %.10e\n",derivs.C_p_p);

    return derivs;
}



struct derivs integrated_de_omp_eq(double xi, double xf, struct fixed_params fixed_params, struct sp_params sp_params, struct polylogs* plgs,struct massive_interpolation_table* mip)
{

int n=TS_NUM_PTS;
double m_phi=fixed_params.m_phi;
double m_s=fixed_params.m_s;
double y=fixed_params.y;
double theta=fixed_params.theta;
double T=sp_params.T;
double Ts=sp_params.Ts;
double Tphi=sp_params.Tphi;
double c_s=sp_params.c_s;
double c_phi=sp_params.c_phi;
double m_a=fixed_params.m_l_a;

double T3=T*T*T;
double Ts3=Ts*Ts*Ts;
double Tphi3=Tphi*Tphi*Tphi;


double p1;


double C_p_as_create=0;
double C_p_as_destroy=0;
double C_as_p_create=0;
double C_as_p_destroy=0;
double C_p_as_exp_term_c=0;
double C_p_as_exp_term_d=0;
double C_as_p_num=0;




struct massless_rates rates;

double m_phi2=m_phi*m_phi;
double y4=y*y*y*y;
double M_as_p = 8 * y*y * (m_phi2 - m_s * m_s);
double del;
double osc_denom;
double integrated_osc=0;
double moment_osc=0;
double weight;

double w;
double ws;
double wphi;
double f_s;
double f_fd;

#pragma omp taskgroup task_reduction(+:C_as_p_destroy) task_reduction(+:C_as_p_create) task_reduction(+:C_p_as_create) task_reduction(+:C_p_as_destroy) task_reduction(+:C_as_p_num) 
{ 
    #pragma omp taskloop grainsize(50) untied in_reduction(+:C_as_p_destroy) in_reduction(+:C_as_p_create) in_reduction(+:C_p_as_create) in_reduction(+:C_p_as_destroy) in_reduction(+:C_as_p_num) 
    for (int k = -n; k < n; ++k)
    {
        p1 = (xf * (ts_quad_fixed_nodes[k + n] + 1) - xi * (ts_quad_fixed_nodes[k + n] - 1)) / 2;
        weight = ts_quad_fixed_weights[k + n];
        w = m_phi2 / (4 * p1 * T);
        ws = m_phi2 / (4 * p1 * Ts);
        wphi = m_phi2 / (4 * p1 * Tphi);
        f_s=1/(exp(p1-c_s)+1);
        f_fd=1/(exp(p1)+1);

        //This is the integral of  f_s*f_a*p1--it is the CONSUMPTION of rho_s via the process a+s->p
        C_as_p_destroy+=weight * p1*f_s*log(1+exp(-w/Ts));//T*T_s*T_s
        //This is the integral of  f_p*p1--it is the PRODUCTION of rho_s via the process p->a+s
        C_as_p_create+=-weight *  p1*log(1-exp(c_phi-wphi/Tphi-p1));//T_phi^3

        //This is the integral of f_s*f_a*(p1+p2)--this is the PRODUCTION of rho_p via the process a+s->p
        C_p_as_exp_term_c=exp(c_s-ws/T);
        C_p_as_create+=weight*f_fd*( (T*p1+Ts/T *ws)*log(1+C_p_as_exp_term_c)-Ts*interp_Li2(-C_p_as_exp_term_c,plgs) ); //Multiply bs T_s*T
        //C_p_as_moment_R+=weight*interp_Li2(-exp(-p1),plgs);
        //This is the integral of f_phi*(p1+p2)--this is the CONSUMPTION of rho_p via the process p->a+s
        C_p_as_exp_term_d=exp(c_phi-wphi/Tphi-p1);
        C_p_as_destroy+=weight*(interp_Li2(C_p_as_exp_term_d,plgs)-(wphi/Tphi+p1)*log(1-C_p_as_exp_term_d)); //Multiply by T_phi^3

 //This is the collision term for the process a+s->p. It represents the production of s via this this process, not the consumption
        C_as_p_num-=weight*(Tphi*Tphi*log(1-exp(c_phi-wphi/Tphi-p1))+Ts*T*f_fd*log(1+exp(c_s-ws/T))); //Done!
    }
}
    double bound_weight = (xf - xi) / 2;
/*printf("C_as_p_destroy: %.10e\n",bound_weight*C_as_p_destroy*T*Ts*Ts);
printf("C_as_p_create: %.10e\n",bound_weight*C_as_p_create*Tphi*Tphi*Tphi);
printf("C_p_as_create: %.10e\n",bound_weight*C_p_as_create*Ts*T);
printf("C_p_as_destroy: %.10e\n",bound_weight*C_p_as_destroy*Tphi*Tphi*Tphi);
printf("C_as_p_num: %.10e\n",bound_weight*C_as_p_num);*/
    double C_as_p_moment_tot=bound_weight*(C_as_p_create*Tphi*Tphi*Tphi-C_as_p_destroy*T*Ts*Ts)*M_as_p/(8*TWO_PI_3);
    double C_p_as_moment_tot=bound_weight*(C_p_as_create*Ts*T-C_p_as_destroy*Tphi*Tphi*Tphi)*M_as_p/(8*TWO_PI_3);
//printf("C_as_p_moment_tot: %.10e\n",C_as_p_moment_tot);
//printf("C_p_as_moment_tot: %.10e\n",C_p_as_moment_tot);
//printf("as->p moment: %.10e\n",bound_weight*(C_as_p_create*Tphi*Tphi*Tphi-C_as_p_destroy*T*Ts*Ts));
//printf("p_>as: %.10e\n",bound_weight*(C_p_as_create*Ts*T-C_p_as_destroy*Tphi*Tphi*Tphi));
    C_as_p_num=bound_weight*C_as_p_num*M_as_p/(8*TWO_PI_3);
//printf("C_as_p_num: %.10e\n",C_as_p_num);
    struct derivs derivs;

    derivs.C_n_s=C_as_p_num;
    derivs.C_p_s=C_as_p_moment_tot;

    derivs.C_n_p=-C_as_p_num;
    derivs.C_p_p=C_p_as_moment_tot;
/*
printf("C_aa_ss: %.10e\n",C_aa_ss*bound_weight);
printf("C_ss_aa: %.10e\n",C_ss_aa*bound_weight);
printf("C_aa_pp: %.10e\n",C_aa_pp*bound_weight);
printf("C_ss_pp: %.10e\n",C_ss_pp*bound_weight);
printf("C_as_p_num: %.10e\n",C_as_p_num);
printf("C_as_p_moment_tot: %.10e\n",C_as_p_moment_tot);
printf("C_p_as_moment_tot: %.10e\n",C_p_as_moment_tot);
printf("C_pp_aa_vals.coll: %.10e\n",y4*C_pp_aa_vals.coll);
printf("C_pp_ss_vals.coll: %.10e\n",y4*C_pp_ss_vals.coll);
printf("C_pp_aa_vals.moment: %.10e\n",y4*C_pp_aa_vals.moment);
printf("C_pp_ss_vals.moment: %.10e\n",y4*C_pp_ss_vals.moment);
printf("C_aa_ss_moment: %.10e\n",C_aa_ss_moment*bound_weight);
printf("C_ss_aa_moment: %.10e\n",C_ss_aa_moment*bound_weight);
printf("C_aa_pp_moment: %.10e\n",C_aa_pp_moment*bound_weight);
printf("C_ss_pp_moment: %.10e\n",C_ss_pp_moment*bound_weight);
printf("integrated_osc: %.10e\n",integrated_osc);
printf("moment_osc: %.10e\n",moment_osc);*/
/*
printf("C_n_s: %.10e\n",derivs.C_n_s);
printf("C_p_s: %.10e\n",derivs.C_p_s);
printf("C_n_p: %.10e\n",derivs.C_n_p);
printf("C_p_p: %.10e\n",derivs.C_p_p);
*/
    return derivs;
}