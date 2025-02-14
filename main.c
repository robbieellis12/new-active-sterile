#include "sn_math.h"
#include "sn_cross_sections.h"
#include "sn_computations.h"
#include "tests.h"
#include "diffeqs.h"
#include "sundials_solve.h"
int main()
{
ts_quad_info(ts_quad_fixed_weights,ts_quad_fixed_nodes,TS_STEP_SIZE,TS_NUM_PTS);
//test_massive_colls();
test_collision_step();
return 0;
//test_collision_interp();
//return 0;

//double m1=Moment1(m,T,c);
//double m2=Moment2(m,T,c);
//printf("Moment1: %f\n",m1);
//printf("Moment2: %f\n",m2);

struct polylogs plgs;
    plgs.N=1000;
    plgs.x0=-2;
    plgs.x1=1;
    //omp_set_num_threads(5);
    fill_polylog_table(&plgs);
//double vav=va(1,0,&plgs);
//printf("va: %.10e\n",vav);
//vav=va(0.01,0,&plgs);
//printf("va: %.10e\n",vav);
//vav=va(100,0,&plgs);
//printf("va: %.10e\n",vav);


//double vap=vp(w-p/T,m_phi/T,c,w);
//printf("vp: %.10e\n",vap);


/*
double vt=v_tot(p,m_a,m_phi,y,&sp_params,&plgs);
printf("v_tot: %.10e\n",vt);
double osc_pk=oscillation_peak_func(p,m_s,cos_2theta,m_a,m_phi,y,&sp_params,&plgs);
printf("oscillation peak: %.10e\n",osc_pk);
double x0=find_oscillation_peak( m_a,  m_s,  cos_2theta,  m_phi,  y,  &sp_params, &plgs);
printf("oscillation peak: %.10e\n",x0);

int N=1000;
double vals[N];
generate_lin_spaced(vals,1,2*x0,N);
for(int i=0;i<N;++i)
{
    vals[i]=v_tot(vals[i],m_a,m_phi,y,&sp_params,&plgs);
}
char* loc="/home/robbie-ellis/VSCode/new active_sterile/exported/v_tot.txt";
export_doubles(loc,vals,N);
*/
//omp_set_num_threads(24);

double m_phi=100;
double theta=1e-8;
double cos_2theta=cos(2*theta);
double y=1;
double m_s=0.7;
double m_l_a=100;
struct sp_params sp_params;
double T=100;
sp_params.T=T;
sp_params.Ts=T;
sp_params.Tphi=T;
sp_params.c_s=0;
sp_params.c_phi=-2;
struct massive_interpolation_table mip;
mip.x0=0;
mip.x1=500;
mip.N=5000;
fill_massive_tables(&mip,y);

struct fixed_params fixed_params;
fixed_params.m_phi=m_phi;
fixed_params.m_s=m_s;
fixed_params.y=y;
fixed_params.theta=theta;
fixed_params.m_l_a=m_l_a;
char* gloc="/home/robbie-ellis/VSCode/new active_sterile/gstar.csv";
read_paired_list(gloc,&fixed_params.T_g,&fixed_params.g_s,&fixed_params.len_g);

//char* loc="/home/robbie-ellis/VSCode/new active_sterile/exported/interaction_data.txt";
//loc="/home/robbie-ellis/VSCode/new active_sterile/exported/collision_data.txt";
/*
struct dx_params dxp;
dxp.fixed_params=&fixed_params;
dxp.plgs=&plgs;
dxp.mip=&mip;*/
/*
struct final_derivs ds=step_omp(sp_params,&dxp);
printf("ts: %.10e\n",ds.dTs);
printf("cs: %.10e\n",ds.dcs);
printf("tp: %.10e\n",ds.dTp);
printf("cp: %.10e\n",ds.dcp);
double x=4*T*T/(m_phi*m_phi);
struct massless_coll_return_vals predictval=ts_quad_massless_coll(x,0,y,&c_vv_vv);
printf("coll: %.10e\n",predictval.f*8*pow(T,6)/(m_phi*m_phi*TWO_PI_4));*/


struct dx_eq_params dxp_eq;
dxp_eq.fixed_params=&fixed_params;
dxp_eq.plgs=&plgs;
dxp_eq.T=T;
//struct final_derivs ds=step_omp2(sp_params,&dxp);
//test_collision_interp();

/*
struct final_derivs ds_eq=step_omp_eq(sp_params,&dxp_eq);
printf("ts: %.10e\n",ds_eq.dTs);
printf("cs: %.10e\n",ds_eq.dcs);
printf("tp: %.10e\n",ds_eq.dTp);
printf("cp: %.10e\n",ds_eq.dcp);*/

double t1=omp_get_wtime();
omp_set_num_threads(24);
#pragma omp parallel master
{
solve_de_eq(T,sp_params, &dxp_eq);
}
double t2=omp_get_wtime();
printf("Time taken: %.10e\n",t2-t1);
/*
omp_set_num_threads(24);
#pragma omp parallel master
{
solve_de(T,sp_params, &dxp);
}*/
/*
struct n_rho_explicit2 nr=ts_quad_n_rho_explicit2(1,1,-1);
printf("n: %.10e\n",nr.n);
printf("p: %.10e\n",nr.p);
printf("nb: %.10e\n",nr.nT);
printf("pb: %.10e\n",nr.pT);
printf("nc: %.10e\n",nr.nc);
printf("pc: %.10e\n",nr.pc);*/
/*
struct final_derivs ds=step_omp2(sp_params,&dxp);
printf("ts: %.10e\n",ds.dTs);
printf("cs: %.10e\n",ds.dcs);
printf("tp: %.10e\n",ds.dTp);
printf("cp: %.10e\n",ds.dcp);
ds=step_omp(sp_params,&dxp);
printf("ts: %.10e\n",ds.dTs);
printf("cs: %.10e\n",ds.dcs);
printf("tp: %.10e\n",ds.dTp);
printf("cp: %.10e\n",ds.dcp);*/

//send_collision_data(loc,fixed_params,sp_params,&plgs,&mip);

//moment_test();
//moment_derivs_test();
//va_derivs_test(&plgs);
/*
struct n_p_vals npv=get_n_p_vals(100,10,3,10,-3,&plgs);
printf("In_p: %.10e\n",npv.In_p);
printf("In_p_b: %.10e\n",npv.In_p_b);
printf("In_p_c: %.10e\n",npv.In_p_c);
printf("Ip_p: %.10e\n",npv.Ip_p);
printf("Ip_p_b: %.10e\n",npv.Ip_p_b);
printf("Ip_p_c: %.10e\n",npv.Ip_p_c);
printf("n_s: %.10e\n",npv.n_s);
printf("p_s: %.10e\n",npv.p_s);
printf("n_s_T: %.10e\n",npv.n_s_T);
printf("n_s_c: %.10e\n",npv.n_s_c);
printf("p_s_T: %.10e\n",npv.p_s_T);
printf("p_s_c: %.10e\n",npv.p_s_c);
*/

/*
struct final_derivs ds=step_omp(sp_params,&dxp);
printf("ts: %.10e\n",ds.dTs);
printf("cs: %.10e\n",ds.dcs);
printf("tp: %.10e\n",ds.dTp);
printf("cp: %.10e\n",ds.dcp);*/


//double t1=omp_get_wtime();
//struct solved_params sp=rk_solve_omp(T,&dxp,-0.00001,50);
/*double t2=omp_get_wtime();
printf("Time taken: %.10e\n",t2-t1);*/

/*
struct derivs ds;
#pragma omp parallel master
{
ds=integrated_de_omp(0,100,fixed_params,sp_params,&plgs,&mip);
}
printf("single thread\n");
struct derivs ds2;
ds2=integrated_de(0,100,fixed_params,sp_params,&plgs,&mip);
printf("1st: %.10e, 2nd: %.10e\n",ds.C_n_s,ds2.C_n_s);
printf("1st: %.10e, 2nd: %.10e\n",ds.C_p_s,ds2.C_p_s);
printf("1st: %.10e, 2nd: %.10e\n",ds.C_n_p,ds2.C_n_p);
printf("1st: %.10e, 2nd: %.10e\n",ds.C_p_p,ds2.C_p_p);*/

/*
double x0=find_oscillation_peak( m_l_a,  m_s,  cos_2theta,  m_phi,  y,  &sp_params, &plgs);
if(x0==-1)
{
    x0=10;
}
printf("oscillation peak: %.10e\n",x0);
int N=1000;

test_rates(1,fixed_params,sp_params,&plgs,&mip);
g_data(loc, N,x0/8,20*x0,fixed_params,sp_params,&plgs,&mip);

struct derivs ds=integrated_de(0,100,fixed_params,sp_params,&plgs,&mip);
printf("C_n_s: %.10e\n",ds.C_n_s);
printf("C_p_s: %.10e\n",ds.C_p_s);
printf("C_n_p: %.10e\n",ds.C_n_p);
printf("C_p_p: %.10e\n",ds.C_p_p);
*/

/*
double vals[N];
double p[N];
generate_lin_spaced(p,0.1,2*x0,N);
#pragma omp parallel for
for(int i=0;i<N;++i)
{
    vals[i]=oscillation_term(p[i],fixed_params,sp_params,&plgs,&mip);
}
char* loc="/home/robbie-ellis/VSCode/new active_sterile/exported/oscillation.txt";
export_doubles(loc,vals,N);
free_polylog_table(&plgs);
*/
free_massive_table(&mip);
free_polylog_table(&plgs);
return 0;
}