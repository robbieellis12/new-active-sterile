#include "tests.h"

void generate_lin_spaced(double* x, double a, double b, int N)
{
    double h=(b-a)/(N-1);
    for(int i=0;i<N;++i)
    {
        x[i]=a+i*h;
    }
}


void export_doubles(char* location, double* data, int size)
{
    FILE* f=fopen(location,"w");
    for(int i=0;i<size;++i)
    {
        fprintf(f,"%.10e\n",data[i]);
    }
    fclose(f);
}


double g_data(char* loc, int N,double a,double b, struct fixed_params fixed_params,struct sp_params sp_params, struct polylogs* plgs,struct massive_interpolation_table* mip)
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

  
    double y2 = y * y;
    double mphi2 = m_phi * m_phi;
    double y4 = y2 * y2;


    double w;
    double w_s;
    double w_phi;
    double v_tot_val;
    double del;
    double gamma_tot;
    double denom;
    double osc;
    double p;
    double g_aa_ss;
    double g_ss_aa;
    double g_as_as;
    double g_sa_as;
    double g_aa_pp;
    double g_ss_pp;
    double g_vp_vp;
    double M_as_p;
    double g_sa_p;
    double g_as_p;
    double massless_coeff;
    struct massless_rates rates;

    double vals[N];
    double vtot_vals[N];
    double g_aa_ss_vals[N];
    double g_ss_aa_vals[N];
    double g_as_as_vals[N];
    double g_sa_as_vals[N];
    double g_aa_pp_vals[N];
    double g_ss_pp_vals[N];
    double g_vp_vp_vals[N];
    double g_as_p_vals[N];
    double g_sa_p_vals[N];
    double osc_vals[N];

    double ps[N];
    generate_lin_spaced(ps,a,b,N);
    //#pragma omp parallel for
    for(int i=0;i<N;++i)
    {
    p=ps[i];
    del = m_s * m_s / (2 * p);
       w = m_phi * m_phi / (4 * p * T);
     w_s = m_phi * m_phi / (4 * p * Ts);
     w_phi = m_phi * m_phi / (4 * p * Tphi);

     v_tot_val = v_tot(p, m_a, m_phi, y, &sp_params, plgs);
    vtot_vals[i]=v_tot_val;
      rates = compute_massless_rates(T, Ts, w, w_s, c_s, y, plgs);

     massless_coeff = y4 * mphi2 / (TWO_PI_2 * 4 * p * p);
     g_aa_ss = rates.g_aa_ss*massless_coeff*T/(w*w);
     g_ss_aa = rates.g_ss_aa*massless_coeff*Ts/(w_s*w_s);
     g_as_as = C_as_as_narrow(p, c_s, m_phi, Ts, y*massless_coeff*Ts/(w_s*w_s));
     g_sa_as = C_as_as_narrow(p, 0, m_phi, T, y)*massless_coeff*T/(w*w);
     g_aa_pp = rates.g_aa_pp*massless_coeff*T/(w*w);
     g_ss_pp = rates.g_ss_pp*massless_coeff*Ts/(w_s*w_s);
     g_vp_vp = compute_vp_vp_rate(m_phi, p, c_phi, Tphi, y, mip) *Tphi*m_phi/(8*TWO_PI_2*p*p);

    g_aa_ss_vals[i]=g_aa_ss;
    g_ss_aa_vals[i]=g_ss_aa;
    g_as_as_vals[i]=g_as_as;
    g_sa_as_vals[i]=g_sa_as;
    g_aa_pp_vals[i]=g_aa_pp;
    g_ss_pp_vals[i]=g_ss_pp;
    g_vp_vp_vals[i]=g_vp_vp;



     M_as_p = 8 * y2 * (mphi2 - m_s * m_s);
     g_sa_p = M_as_p / (32 * M_PI * p * p) * (T * log(1 + exp(-w)));   // check
     g_as_p = M_as_p / (32 * M_PI * p * p) * (Ts * log(1 + exp(-w_s))); // check
    g_as_p_vals[i]=g_as_p;
    g_sa_p_vals[i]=g_sa_p;

     gamma_tot = g_aa_ss + g_aa_pp + g_as_p + g_as_as + g_ss_aa + g_ss_pp + g_sa_p + g_sa_as + 2 * g_vp_vp;
     denom = del * del * sin2_2theta + gamma_tot * gamma_tot / 4 + (del * cos_2theta - v_tot_val) * (del * cos_2theta - v_tot_val);
     osc= gamma_tot / 4 * del * del * sin2_2theta / denom;
    osc_vals[i]=osc;
}

 FILE* f=fopen(loc,"w");
 fprintf(f,"p,v_tot,g_aa_ss,g_ss_aa,g_as_as,g_sa_as,g_aa_pp,g_ss_pp,g_vp_vp,g_as_p,g_sa_p,osc\n");
    for(int i=0;i<N;++i)
    {
        fprintf(f, "%.10e, %.10e, %.10e, %.10e, %.10e, %.10e, %.10e, %.10e, %.10e, %.10e, %.10e, %.10e\n",ps[i],del*cos_2theta-vtot_vals[i],g_aa_ss_vals[i],g_ss_aa_vals[i],g_as_as_vals[i],g_sa_as_vals[i],g_aa_pp_vals[i],g_ss_pp_vals[i],g_vp_vp_vals[i],g_as_p_vals[i],g_sa_p_vals[i],osc_vals[i]);
    }
    fclose(f);
}

void test_rates(double p, struct fixed_params fixed_params,struct sp_params sp_params, struct polylogs* plgs,struct massive_interpolation_table* mip)
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

    double T3=T*T*T;
    double Ts3=Ts*Ts*Ts;
  
    double y2 = y * y;
    double mphi2 = m_phi * m_phi;
    double y4 = y2 * y2;
    double w = m_phi * m_phi / (4 * p * T);
    double w_s = m_phi * m_phi / (4 * p * Ts);

    struct massless_rates rates = compute_massless_rates(T, Ts,w, w_s, c_s, y, plgs);
    double massless_coeff = y4 /(M_PI*M_PI * mphi2);

    double g_aa_ss = rates.g_aa_ss*massless_coeff*T3;
    double g_ss_aa = rates.g_ss_aa*massless_coeff*Ts3;
    double g_as_as = C_as_as_narrow(p, c_s, m_phi, Ts, y);
    double g_sa_as = C_as_as_narrow(p, 0, m_phi, T, y);
    double g_aa_pp = rates.g_aa_pp*massless_coeff*T3;
    double g_ss_pp = rates.g_ss_pp*massless_coeff*Ts3;

    double M_as_p = 8 * y2 * (mphi2 - m_s * m_s);
    double g_sa_p = M_as_p / (32 * M_PI * p * p) * (T * log(1 + exp(-w)));   // check
    double g_as_p = M_as_p / (32 * M_PI * p * p) * (Ts * log(1 + exp(-w_s))); // check

    double g_vp_vp = y4*compute_vp_vp_rate(m_phi, p, c_phi, Tphi, y, mip)*Tphi*m_phi/(8*TWO_PI_2*p*p);

    printf("g_aa_ss: %.10e\n",g_aa_ss);
    printf("g_ss_aa: %.10e\n",g_ss_aa);
    printf("g_as_as: %.10e\n",g_as_as);
    printf("g_sa_as: %.10e\n",g_sa_as);
    printf("g_aa_pp: %.10e\n",g_aa_pp);
    printf("g_ss_pp: %.10e\n",g_ss_pp);
    printf("g_vp_vp: %.10e\n",g_vp_vp);
    printf("g_sa_p: %.10e\n",g_sa_p);
    printf("g_as_p: %.10e\n",g_as_p);

}



















void integrated_de_omp_save(struct collision_data* cd, double xi, double xf, struct fixed_params fixed_params, struct sp_params sp_params, struct polylogs* plgs,struct massive_interpolation_table* mip)
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
    #pragma omp taskloop grainsize(5) in_reduction(+:C_as_p_destroy) in_reduction(+:C_as_p_create) in_reduction(+:C_p_as_create) in_reduction(+:C_p_as_destroy) in_reduction(+:C_as_p_num) in_reduction(+:integrated_osc) in_reduction(+:moment_osc) in_reduction(+:C_aa_ss) in_reduction(+:C_ss_aa) in_reduction(+:C_aa_pp) in_reduction(+:C_ss_pp) in_reduction(+:C_aa_ss_moment) in_reduction(+:C_ss_aa_moment) in_reduction(+:C_aa_pp_moment) in_reduction(+:C_ss_pp_moment)
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
}
 
    double bound_weight = (xf - xi) / 2;

    double C_as_p_moment_tot=bound_weight*(C_as_p_create*Tphi*Tphi*Tphi-C_as_p_destroy*T*Ts*Ts)*M_as_p/(8*TWO_PI_3);
    double C_p_as_moment_tot=bound_weight*(C_p_as_create*Ts*T-C_p_as_destroy*Tphi*Tphi*Tphi)*M_as_p/(8*TWO_PI_3);

    C_as_p_num=bound_weight*C_as_p_num*M_as_p/(8*TWO_PI_3);

   // struct derivs derivs;

    //derivs.C_n_s=bound_weight*(integrated_osc+C_aa_ss-C_ss_aa)+C_as_p_num;
    //derivs.C_p_s=bound_weight*(moment_osc+C_aa_ss_moment-C_ss_aa_moment)+C_as_p_moment_tot;

    //derivs.C_n_p=bound_weight*(C_aa_pp+C_ss_pp)-y4*C_pp_aa_vals.coll-y4*C_pp_ss_vals.coll-C_as_p_num;
    //derivs.C_p_p=bound_weight*(C_aa_pp_moment+C_ss_pp_moment)-y4*C_pp_aa_vals.moment-y4*C_pp_ss_vals.moment+C_p_as_moment_tot;
cd->C_aa_ss=C_aa_ss*bound_weight;
cd->C_ss_aa=C_ss_aa*bound_weight;
cd->C_aa_pp=C_aa_pp*bound_weight;
cd->C_ss_pp=C_ss_pp*bound_weight;
cd->C_as_p_num=C_as_p_num;
cd->C_as_p_moment=C_as_p_moment_tot;
cd->C_p_as_moment=C_p_as_moment_tot;
cd->C_aa_ss_moment=C_aa_ss_moment*bound_weight;
cd->C_ss_aa_moment=C_ss_aa_moment*bound_weight;
cd->C_aa_pp_moment=C_aa_pp_moment*bound_weight;
cd->C_ss_pp_moment=C_ss_pp_moment*bound_weight;
cd->oscillation=integrated_osc;
cd->oscillation_moment=moment_osc;
}

void send_collision_data(char* loc,struct fixed_params fixed_params, struct sp_params sp_params, struct polylogs* plgs,struct massive_interpolation_table* mip)
{
    double m_phi=fixed_params.m_phi;
    double T=sp_params.T;
    double Ts=sp_params.Ts;
    double c_s=sp_params.c_s;
    double Tphi=sp_params.Tphi;
    double m_s=fixed_params.m_s;
    double y=fixed_params.y;
    double theta=fixed_params.theta;
    double m_a=fixed_params.m_l_a;
double cos_2theta=cos(2*theta);
    FILE* f=fopen(loc,"w");
    if(f==NULL)
    {
        printf("Error opening file\n");
        return;
    }
    int num_pts=100;
    double h=5;
    struct collision_data cd[num_pts];
    struct collision_data temp_cd[2];
    double x0=0;
    //label columns of output seperated by spaces:
    fprintf(f,"C_aa_ss,    C_ss_aa,    C_aa_pp,    C_ss_pp,    C_pp_aa,    C_pp_ss,    C_aa_ss_moment,    C_ss_aa_moment,    C_aa_pp_moment,    C_ss_pp_moment,    C_as_p_num,    C_as_p,    C_p_as,    C_p_as_moment,    C_as_p_moment,    oscillation,    oscillation_moment\n");
    //fprintf(f,"C_aa_ss,C_ss_aa,C_aa_pp,C_ss_pp,C_pp_aa,C_pp_ss,C_aa_ss_moment,C_ss_aa_moment,C_aa_pp_moment,C_ss_pp_moment,C_as_p_num,C_as_p,C_p_as,C_p_as_moment,C_as_p_moment,oscillation,oscillation_moment\n");
    #pragma omp parallel master
    {
    for(int i=0;i<num_pts;++i)
    {
        sp_params.T=T+h*i;
        //sp_params.c_phi=m_phi/sp_params.Tphi -1/(1+i);
        x0=find_oscillation_peak(m_a,m_s,cos_2theta,m_phi,y,&sp_params,plgs);
        if(x0!=-1)
        {
        integrated_de_omp_save(&temp_cd[1],0,x0,fixed_params,sp_params,plgs,mip);
        integrated_de_omp_save(&temp_cd[2],x0,100,fixed_params,sp_params,plgs,mip);
        cd[i].C_aa_pp=temp_cd[1].C_aa_pp+temp_cd[2].C_aa_pp;
        cd[i].C_ss_pp=temp_cd[1].C_ss_pp+temp_cd[2].C_ss_pp;
        cd[i].C_aa_ss=temp_cd[1].C_aa_ss+temp_cd[2].C_aa_ss;
        cd[i].C_ss_aa=temp_cd[1].C_ss_aa+temp_cd[2].C_ss_aa;
        cd[i].C_aa_pp_moment=temp_cd[1].C_aa_pp_moment+temp_cd[2].C_aa_pp_moment;
        cd[i].C_ss_pp_moment=temp_cd[1].C_ss_pp_moment+temp_cd[2].C_ss_pp_moment;
        cd[i].C_aa_ss_moment=temp_cd[1].C_aa_ss_moment+temp_cd[2].C_aa_ss_moment;
        cd[i].C_ss_aa_moment=temp_cd[1].C_ss_aa_moment+temp_cd[2].C_ss_aa_moment;
        cd[i].C_as_p_num=temp_cd[1].C_as_p_num+temp_cd[2].C_as_p_num;
        cd[i].C_as_p_moment=temp_cd[1].C_as_p_moment+temp_cd[2].C_as_p_moment;
        cd[i].C_p_as_moment=temp_cd[1].C_p_as_moment+temp_cd[2].C_p_as_moment;
        cd[i].oscillation=temp_cd[1].oscillation+temp_cd[2].oscillation;
        cd[i].oscillation_moment=temp_cd[1].oscillation_moment+temp_cd[2].oscillation_moment;
        }
        else
        {
            integrated_de_omp_save(&cd[i],0,100,fixed_params,sp_params,plgs,mip);
        }
       struct collision_vals C_pp_aa_vals=Moment1_omp(m_phi,sp_params.T,0);
    struct collision_vals C_pp_ss_vals=Moment1_omp(m_phi,sp_params.Ts,sp_params.c_s);
    
    fprintf(f,"%.10e, %.10e, %.10e, %.10e, %.10e, %.10e, %.10e, %.10e, %.10e, %.10e, %.10e, %.10e, %.10e, %.10e, %.10e, %.10e, %.10e\n",cd[i].C_aa_ss,cd[i].C_ss_aa,cd[i].C_aa_pp,cd[i].C_ss_pp,C_pp_aa_vals.coll,C_pp_ss_vals.coll,cd[i].C_aa_ss_moment,cd[i].C_ss_aa_moment,cd[i].C_aa_pp_moment,cd[i].C_ss_pp_moment,cd[i].C_as_p_num,cd[i].C_as_p_moment,cd[i].C_p_as_moment,cd[i].C_aa_ss_moment,cd[i].C_ss_aa_moment,cd[i].oscillation,cd[i].oscillation_moment);
    
    }
    }
}



struct collision_vals Moment2all(double m, double T, double c)
{
    int n=TS_NUM_PTS;
    double  u1,u2,u3;
double intval_in=0;
double sigint = 0;
double int_val = 0;
double momintval=0;

double b = m / T;

double xf = 120;
double xi = 0;

double zf = 120;
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
   // momintval+=ts_quad_fixed_weights[j + n] * ts_quad_fixed_weights[k + n] * u1  *u1* u2*u2*f2*f1 / (E1*E2) * sigint*(E1+E2)/2;
    }
    int_val += ts_quad_fixed_weights[k + n] * u1  *u1* f1 / E1 * intval_in;
    momintval+=ts_quad_fixed_weights[k + n] * u1  *u1* f1 * intval_in;
}
struct collision_vals vals;
vals.coll=(tf - ti) / 2 * (zf - zi) / 2 * (xf - xi) / 2 *int_val*T*T*T*T/(2*TWO_PI_2);
vals.moment=(tf - ti) / 2 * (zf - zi) / 2 * (xf - xi) / 2 * momintval*T*T*T*T*T/(2*TWO_PI_2);
    return vals;
}
struct collision_vals Moment1new(double m, double T, double c)
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
        zf=50/E;
        if(zf<1) continue;
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
  vals.coll=(xf-xi)/2*int_val*16*T*T*T*T/(2*TWO_PI_2);
  vals.moment=(xf-xi)/2*dint_val*16*T*T*T*T*T/(2*TWO_PI_2);
    return vals;
}

void moment_test()
{
    double m=100;
    double T=400;
    double c=0;
    struct collision_vals vals=Moment2all(m,T,c);
    printf("Collision: %.10e\n",vals.coll);
    printf("Moment: %.10e\n",vals.moment);
    struct collision_vals vals2=Moment1(m,T,c);
    printf("Collision2: %.10e\n",vals2.coll);
    printf("Moment2: %.10e\n",vals2.moment);
    struct collision_vals vals3=Moment1new(m,T,c);
    printf("Collision: %.10e\n",vals3.coll);
    printf("Moment: %.10e\n",vals3.moment);
    printf("Moment: %.10e\n",vals3.moment/vals.moment);
}

void moment_derivs_test()
{
#pragma omp parallel master
{
double m=2;
double T=5;
double c=-1;
double h=0.0001;

struct collision_vals_derivs cvd=Collision_deriv( m,  T,  c);

struct collision_vals vals_T_plus=Moment1_omp(m,T+h,c);
struct collision_vals vals_T_minus=Moment1_omp(m,T-h,c);

struct collision_vals vals_c_plus=Moment1_omp(m,T,c+h);
struct collision_vals vals_c_minus=Moment1_omp(m,T,c-h);

printf("Collision deriv: %.10e\n",cvd.coll_dT);
printf("Estimated collision deriv: %.10e\n",(vals_T_plus.coll-vals_T_minus.coll)/(2*h));

printf("Moment deriv: %.10e\n",cvd.moment_dT);
printf("Estimated moment deriv: %.10e\n",(vals_T_plus.moment-vals_T_minus.moment)/(2*h));

printf("Collision deriv (c): %.10e\n",cvd.coll_dc);
printf("Estimated collision deriv (c): %.10e\n",(vals_c_plus.coll-vals_c_minus.coll)/(2*h));

printf("Moment deriv (c): %.10e\n",cvd.moment_dc);
printf("Estimated moment deriv (c): %.10e\n",(vals_c_plus.moment-vals_c_minus.moment)/(2*h));

}
}

void va_derivs_test(struct polylogs* plgs)
{

    double m=4;
double T=5;
double p=1;
double c=-2;
double w=m*m/(4*p*T);
double a=w-p/T;
double b=m/T;
double h=0.0001;
double w_plus=m*m/(4*p*(T+h));
double w_minus=m*m/(4*p*(T-h));
double a_plus=w_plus-p/(T+h);
double a_minus=w_minus-p/(T-h);
double b_plus=m/(T+h);
double b_minus=m/(T-h);

double va_plus=(T+h)*va(w_plus,c,plgs);
double va_minus=(T-h)*va(w_minus,c,plgs);
double va_c_plus=T*va(w,c+h,plgs);
double va_c_minus=T*va(w,c-h,plgs);

struct va_derivs vd=va_deriv(m,T,p,c,plgs);

printf("va deriv: %.10e\n",vd.va_dT);
printf("Estimated va deriv: %.10e\n",(va_plus-va_minus)/(2*h));

printf("va deriv (c): %.10e\n",vd.va_dc);
printf("Estimated va deriv (c): %.10e\n",(va_c_plus-va_c_minus)/(2*h));

struct vp_derivs vpd=vp_deriv(m,T,p,c);
double vp_plus=(T+h)*vp(a_plus,b_plus,c,w_plus);
double vp_minus=(T-h)*vp(a_minus,b_minus,c,w_minus);
double vp_c_plus=T*vp(a,b,c+h,w);
double vp_c_minus=T*vp(a,b,c-h,w);

printf("vp deriv: %.10e\n",vpd.vp_dT);
printf("Estimated vp deriv: %.10e\n",(vp_plus-vp_minus)/(2*h));

printf("vp deriv (c): %.10e\n",vpd.vp_dc);
printf("Estimated vp deriv (c): %.10e\n",(vp_c_plus-vp_c_minus)/(2*h));

printf("vp val: %.10e\n",vpd.vp);



}








void test_collision_interp()
{
    struct massless_coll_table mct;
    mct.x0=0.1;
    mct.x1=20;
    mct.Nx=200;
    mct.c0=0;
    mct.c1=20;
    mct.Nc=50;
    double m=100;
    mct.m=m;
    mct.cs=&c_vv_vv;
    double y=1;
    double t1=omp_get_wtime();
    fill_massless_coll_table(&mct,y);
    double t2=omp_get_wtime();
    printf("Time to fill table: %.10e\n",t2-t1);
    double x=0.2;
    double c=1.45;
    double testval=interp_massless_col(x,c,&mct);
    printf("Interp val: %.10e\n",testval);
    double predictval=ts_quad_massless_coll(x,c,y,&c_vv_vv);
    printf("Predicted val: %.10e\n",predictval);
/*
    struct massless_coll_vals vals=get_massless_coll_vals(x,c,y,&c_vv_vv);
double h=0.001;
    struct massless_coll_vals valsT=get_massless_coll_vals(x+h,c,y,&c_vv_vv);
    struct massless_coll_vals valsc=get_massless_coll_vals(x,c+h,y,&c_vv_vv);
    struct massless_coll_vals valsTc=get_massless_coll_vals(x+h,c+h,y,&c_vv_vv);
    //printf("Predicted val f: %.10e\n",vals.f);
    printf("Predicted val f_dt: %.10e\n",vals.f_dx);
    printf("Test val f_dt: %.10e\n",(valsT.f-vals.f)/h);
    printf("Predicted val f_dc: %.10e\n",vals.f_dc);
    printf("Test val f_dc: %.10e\n",(valsc.f-vals.f)/h);
    printf("Predicted val f_dt_dc: %.10e\n",vals.f_dxc);
    printf("Test val f_dt_dc: %.10e\n",(valsTc.f-valsT.f-valsc.f+vals.f)/(h*h));*/
    free_massless_coll_table(&mct);
}


void test_massive_colls()
{
struct massive_coll_table mct;


    double m=1;
    double b=1;
    double T=m/b;
    double c=0;
//struct massive_coll_vals vals=get_collision_vals2(m,b,c);
/*
struct massive_coll_vals vals1=get_collision_vals(m,b,c);
    struct collision_vals vals2=Moment1_omp(m,T,c);
    printf("Interp val: %.10e\n",vals.f*8*T/(m*m*m*TWO_PI_4));
    //printf("Interp val: %.10e\n",vals1.f);

    printf("Predicted val: %.10e\n",vals2.coll);
    printf("ratio: %.10e\n",(vals2.coll/(vals.f)));*/




    mct.b0=0.1;
    mct.b1=5;
    mct.Nb=60;
    mct.c0=0;
    mct.c1=5;
    mct.Nc=50;
    mct.m=m;
    double y=1;
    fill_massive_coll_table(&mct,y);





    struct collision_vals vals=interp_massive_colls(m,b,c,&mct);
    struct collision_vals vals2=Moment1_omp(m,m/b,c);
    printf("Interp val: %.10e\n",vals.coll);
    printf("Predicted val: %.10e\n",vals2.coll);
    printf("Interp moment: %.10e\n",vals.moment);
    printf("Predicted moment: %.10e\n",vals2.moment);
    printf("ratio: %.10e\n",vals2.moment/vals.moment);
    printf("ratio: %.10e\n",vals2.coll/vals.coll);
    free_massive_coll_table(&mct);
    double h=1e-5;
    /*
    struct massive_coll_vals vals=get_collision_vals(m,b,c);
    struct massive_coll_vals valsb=get_collision_vals(m,b+h,c);
    struct massive_coll_vals valsc=get_collision_vals(m,b,c+h);
    printf("f: %.10e\n",vals.f);
    printf("f_db: %.10e\n",vals.f_db);
    printf("Estimated f_db: %.10e\n",(valsb.f-vals.f)/h);

    printf("f_dc: %.10e\n",vals.f_dc);
    printf("Estimated f_dc: %.10e\n",(valsc.f-vals.f)/h);

    printf("f_dbc: %.10e\n",vals.f_dbc);
    printf("Estimated f_dbc: %.10e\n",(get_collision_vals(m,b+h,c+h).f-valsb.f-valsc.f+vals.f)/(h*h));
    
    printf("f_moment: %.10e\n",vals.f_moment);

    printf("f_moment_db: %.10e\n",vals.f_db_moment);
     printf("Estimated f_moment_db: %.10e\n",(valsb.f_moment-vals.f_moment)/h);
    
    printf("f_moment_dc: %.10e\n",vals.f_dc_moment);
     printf("Estimated f_moment_dc: %.10e\n",(valsc.f_moment-vals.f_moment)/h);
   
    printf("f_moment_dbc: %.10e\n",vals.f_dbc_moment);
   printf("Estimated f_moment_dbc: %.10e\n",(get_collision_vals(m,b+h,c+h).f_moment-valsb.f_moment-valsc.f_moment+vals.f_moment)/(h*h));
*/
   //repeat with new function:
   struct massive_coll_vals valsnew=get_collision_vals2(m,b,c);
   struct massive_coll_vals valsnewb=get_collision_vals2(m,b+h,c);
   struct massive_coll_vals valsnewc=get_collision_vals2(m,b,c+h);
   struct massive_coll_vals valsnewbc=get_collision_vals2(m,b+h,c+h);

   printf("f: %.10e\n",valsnew.f);
   printf("f_db: %.10e\n",valsnew.f_db);
    printf("Estimated f_db: %.10e\n",(valsnewb.f-valsnew.f)/h);

    printf("f_dc: %.10e\n",valsnew.f_dc);
    printf("Estimated f_dc: %.10e\n",(valsnewc.f-valsnew.f)/h);

    printf("f_dbc: %.10e\n",valsnew.f_dbc);
    printf("Estimated f_dbc: %.10e\n",(valsnewbc.f-valsnewb.f-valsnewc.f+valsnew.f)/(h*h));

    printf("f_moment: %.10e\n",valsnew.f_moment);

    printf("f_moment_db: %.10e\n",valsnew.f_db_moment);
     printf("Estimated f_moment_db: %.10e\n",(valsnewb.f_moment-valsnew.f_moment)/h);

    printf("f_moment_dc: %.10e\n",valsnew.f_dc_moment);
        printf("Estimated f_moment_dc: %.10e\n",(valsnewc.f_moment-valsnew.f_moment)/h);

    printf("f_moment_dbc: %.10e\n",valsnew.f_dbc_moment);
    printf("Estimated f_moment_dbc: %.10e\n",(valsnewbc.f_moment-valsnewb.f_moment-valsnewc.f_moment+valsnew.f_moment)/(h*h));
   

}