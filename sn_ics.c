#include "sn_ics.h"

/*
struct massless_coll_vals get_massless_coll_vals(double m, double T, double c,double y,double (*f)(double,double))
{
    int n=TS_NUM_PTS;
 double  u,p1;

    double int_val=0.0;
    double int_val_dT=0.0;
    double int_val_dc=0.0;
    double int_val_dTc=0.0;

    double int_in=0;
    double int_in_dT=0;
    double int_in_dc=0;
    double int_in_dTc=0;


    double xf=100;
    double xi=0;

    double expterm_out,expterm_in,fd_in,fd_out,ln_term;
    double I_coeff=4*T*T/(m*m);
    double cs_term;
    for(int k=-n;k<n;++k)
    {
        p1=( xf*(ts_quad_fixed_nodes[k+n]+1) - xi*(ts_quad_fixed_nodes[k+n]-1) )/2;
        expterm_out=exp(p1-c);
        fd_out=1/(expterm_out+1);
        for(int l=-n;l<n;++l)
        {
            u=( xf*(ts_quad_fixed_nodes[l+n]+1) - xi*(ts_quad_fixed_nodes[l+n]-1) )/2;
            expterm_in=exp(u-c);
            ln_term=log(1+1/expterm_in);
            fd_in=1/(expterm_in+1);
            cs_term=f(I_coeff*p1*u,y);

            int_in+=ts_quad_fixed_weights[l+n]*u*ln_term*cs_term;
            int_in_dT+=ts_quad_fixed_weights[l+n]*u*(-ln_term+2*u*fd_in)*cs_term;
            int_in_dc+=ts_quad_fixed_weights[l+n]*u*(expterm_out*fd_out*ln_term+fd_in)*cs_term;
           // int_in_dTc+=ts_quad_fixed_weights[l+n]*u*(expterm_out*fd_out*fd_out*ln_term+fd_in*fd_out)*cs_term;
            int_in_dTc+=ts_quad_fixed_weights[l+n]*u*(fd_out*expterm_out*(-ln_term+2*u*fd_in)+fd_in*(-1+2*u*expterm_in*fd_in))*cs_term;
             }
        int_val+=ts_quad_fixed_weights[k+n]*p1*p1*fd_out*int_in;
        int_val_dT+=ts_quad_fixed_weights[k+n]*p1*p1*fd_out*int_in_dT;
        int_val_dc+=ts_quad_fixed_weights[k+n]*p1*p1*fd_out*int_in_dc;
        int_val_dTc+=ts_quad_fixed_weights[k+n]*p1*p1*fd_out*int_in_dTc;
        
        int_in=0;
        int_in_dT=0;
        int_in_dc=0;
        int_in_dTc=0;
    }
    double bound_weight=(xf-xi)*(xf-xi)/4;
   struct massless_coll_vals vals;
   vals.f=T*T*T*bound_weight*int_val;
    vals.f_dT=T*T*bound_weight*int_val_dT;
    vals.f_dc=T*T*T*bound_weight*int_val_dc;
    vals.f_dTc=T*T*bound_weight*int_val_dTc;
    return vals;
}

double ts_quad_massless_coll(double m, double T, double c,double y,double (*f)(double,double))
{
    int n=TS_NUM_PTS;
    double  u,p1;

    double int_val=0.0;
    double int_in=0;

    double xf=100;
    double xi=0;

    double I_coeff=4*T*T/(m*m);
    double cs_term;
    for(int k=-n;k<n;++k)
    {
        p1=( xf*(ts_quad_fixed_nodes[k+n]+1) - xi*(ts_quad_fixed_nodes[k+n]-1) )/2;
        for(int l=-n;l<n;++l)
        {
            u=( xf*(ts_quad_fixed_nodes[l+n]+1) - xi*(ts_quad_fixed_nodes[l+n]-1) )/2;

            int_in+=ts_quad_fixed_weights[l+n]*u*log(1+exp(c-u))*f(I_coeff*p1*u,y);
        }
        int_val+=ts_quad_fixed_weights[k+n]*p1*p1*1/(exp(p1-c)+1)*int_in;
        int_in=0;
    }
    return (xf-xi)*(xf-xi)/4*int_val*T*T*T*T*T*T/(m*m*M_PI*M_PI);
}

void fill_massless_coll_table(struct massless_coll_table* mct,double y)
{
    double T0=mct->T0;
    double T1=mct->T1;
    int NT=mct->NT;
    mct->hT=(T1-T0)/(NT-1.0);
    double hT=mct->hT;

    double c0=mct->c0;
    double c1=mct->c1;
    int Nc=mct->Nc;
    mct->hc=(c1-c0)/(Nc-1.0);
    double hc=mct->hc;

    mct->table=malloc(NT*Nc*sizeof(double));
    mct->table_dT=malloc(NT*Nc*sizeof(double));
    mct->table_dc=malloc(NT*Nc*sizeof(double));
    mct->table_dTc=malloc(NT*Nc*sizeof(double));

    double m=mct->m;
    double (*f)(double,double)=mct->cs;
    struct massless_coll_vals mcv;
    #pragma omp parallel for collapse(2) schedule(dynamic)
    for(int i=0;i<NT;++i)
    { 
        for(int j=0;j<Nc;++j)
        {
        mcv=get_massless_coll_vals(m,T0+i*hT,c0+j*hc,y,f);
        mct->table[j+i*Nc]=mcv.f;
        mct->table_dT[j+i*Nc]=mcv.f_dT;
        mct->table_dc[j+i*Nc]=mcv.f_dc;
        mct->table_dTc[j+i*Nc]=mcv.f_dTc;
        }
    }
}
double interp_massless_col(double m,double T,double c,struct massless_coll_table* mct)
{

double T0=mct->T0;
double hT=mct->hT;
double c0=mct->c0;
double hc=mct->hc;
int NT=mct->NT;
int Nc=mct->Nc;

int idx_T=floor((T-T0)/hT);
int idx_c=floor((c-c0)/hc);

if(idx_T<NT-1 && idx_T>=0 && idx_c<Nc-1 && idx_c>=0)
{
    printf("Interpolating\n");
double pts[4]={mct->table[idx_c+idx_T*Nc],mct->table[idx_c+1+idx_T*Nc],mct->table[idx_c+(idx_T+1)*Nc],mct->table[idx_c+1+(idx_T+1)*Nc]};
double dT_pts[4]={mct->table_dT[idx_c+idx_T*Nc],mct->table_dT[idx_c+1+idx_T*Nc],mct->table_dT[idx_c+(idx_T+1)*Nc],mct->table_dT[idx_c+1+(idx_T+1)*Nc]};
double dc_pts[4]={mct->table_dc[idx_c+idx_T*Nc],mct->table_dc[idx_c+1+idx_T*Nc],mct->table_dc[idx_c+(idx_T+1)*Nc],mct->table_dc[idx_c+1+(idx_T+1)*Nc]};
double dTc_pts[4]={mct->table_dTc[idx_c+idx_T*Nc],mct->table_dTc[idx_c+1+idx_T*Nc],mct->table_dTc[idx_c+(idx_T+1)*Nc],mct->table_dTc[idx_c+1+(idx_T+1)*Nc]};
return two_point_2d_hermite_interp(T,c,T0+idx_T*hT,c0+idx_c*hc,T0+(idx_T+1)*hT,c0+(idx_c+1)*hc,pts,dT_pts,dc_pts,dTc_pts);//*T*T*T/(m*m*M_PI*M_PI);
}
    return ts_quad_massless_coll(m,T,c,mct->m,mct->cs);//*T*T*T/(m*m*M_PI*M_PI);
}

void free_massless_coll_table(struct massless_coll_table* mct)
{
    free(mct->table);
    free(mct->table_dT);
    free(mct->table_dc);
    free(mct->table_dTc);
}
*/




struct massless_coll_vals get_massless_coll_vals(double x, double c,double y,double (*f)(double,double))
{
    int n=TS_NUM_PTS;
 double  u,p1;

    double int_val=0.0;
    double int_val_dx=0.0;
    double int_val_dc=0.0;
    double int_val_dxc=0.0;

    double int_val_moment=0.0;
    double int_val_dx_moment=0.0;
    double int_val_dc_moment=0.0;
    double int_val_dxc_moment=0.0;

    double int_in=0;
    double int_in_dx=0;
    double int_in_dc=0;
    double int_in_dxc=0;


    double xf=100;
    double xi=0;

    double expterm_out,expterm_in,fd_in,fd_out,ln_term;
    double cs_term;
    for(int k=-n;k<n;++k)
    {
        p1=( xf*(ts_quad_fixed_nodes[k+n]+1) - xi*(ts_quad_fixed_nodes[k+n]-1) )/2;
        expterm_out=exp(p1-c);
        fd_out=1/(expterm_out+1);
        for(int l=-n;l<n;++l)
        {
            u=( xf*(ts_quad_fixed_nodes[l+n]+1) - xi*(ts_quad_fixed_nodes[l+n]-1) )/2;
            expterm_in=exp(u-c);
            ln_term=log(1+1/expterm_in);
            fd_in=1/(expterm_in+1);
            cs_term=f(x*p1*u,y);

            int_in+=ts_quad_fixed_weights[l+n]*u*ln_term*cs_term;
            int_in_dx+=ts_quad_fixed_weights[l+n]*u*(-2*ln_term+u*fd_in)*cs_term;
            int_in_dc+=ts_quad_fixed_weights[l+n]*u*(expterm_out*fd_out*ln_term+fd_in)*cs_term;
           // int_in_dTc+=ts_quad_fixed_weights[l+n]*u*(expterm_out*fd_out*fd_out*ln_term+fd_in*fd_out)*cs_term;
            int_in_dxc+=ts_quad_fixed_weights[l+n]*u*(fd_out*expterm_out*(-2*ln_term+u*fd_in)+fd_in*(-2+u*expterm_in*fd_in))*cs_term;
             }
        int_val+=ts_quad_fixed_weights[k+n]*p1*p1*fd_out*int_in;
        int_val_dx+=ts_quad_fixed_weights[k+n]*p1*p1*fd_out*int_in_dx;
        int_val_dc+=ts_quad_fixed_weights[k+n]*p1*p1*fd_out*int_in_dc;
        int_val_dxc+=ts_quad_fixed_weights[k+n]*p1*p1*fd_out*int_in_dxc;

        int_val_moment+=ts_quad_fixed_weights[k+n]*p1*p1*p1*fd_out*int_in;
        int_val_dx_moment+=ts_quad_fixed_weights[k+n]*p1*p1*p1*fd_out*int_in_dx;
        int_val_dc_moment+=ts_quad_fixed_weights[k+n]*p1*p1*p1*fd_out*int_in_dc;
        int_val_dxc_moment+=ts_quad_fixed_weights[k+n]*p1*p1*p1*fd_out*int_in_dxc;
        
        int_in=0;
        int_in_dx=0;
        int_in_dc=0;
        int_in_dxc=0;
    }
    double bound_weight=(xf-xi)*(xf-xi)/4;
   struct massless_coll_vals vals;
   vals.f=bound_weight*int_val;
    vals.f_dx=bound_weight*int_val_dx/x;
    vals.f_dc=bound_weight*int_val_dc;
    vals.f_dxc=bound_weight*int_val_dxc/x;

    vals.f_moment=bound_weight*int_val_moment;
    vals.f_dx_moment=bound_weight*int_val_dx_moment/x;
    vals.f_dc_moment=bound_weight*int_val_dc_moment;
    vals.f_dxc_moment=bound_weight*int_val_dxc_moment/x;
    return vals;
}

struct massless_coll_return_vals ts_quad_massless_coll(double x, double c,double y,double (*f)(double,double))
{
    int n=TS_NUM_PTS;
    double  u,p1;

    double int_val=0.0;
    double int_val_moment=0.0;
    double int_in=0;

    double xf=100;
    double xi=0;

    double cs_term;

    double fd_pp_term=0;
    #pragma omp taskgroup task_reduction(+:int_val,int_val_moment)
    {
    #pragma omp taskloop grainsize(10) in_reduction(+:int_val,int_val_moment)
    for(int k=-n;k<n;++k)
    {
        p1=( xf*(ts_quad_fixed_nodes[k+n]+1) - xi*(ts_quad_fixed_nodes[k+n]-1) )/2;
        for(int l=-n;l<n;++l)
        {
            u=( xf*(ts_quad_fixed_nodes[l+n]+1) - xi*(ts_quad_fixed_nodes[l+n]-1) )/2;

            int_in+=ts_quad_fixed_weights[l+n]*u*log(1+exp(c-u))*f(x*p1*u,y);
        }
        fd_pp_term=p1*p1/(exp(p1-c)+1);
        int_val+=ts_quad_fixed_weights[k+n]*fd_pp_term*int_in;
        int_val_moment+=ts_quad_fixed_weights[k+n]*p1*fd_pp_term*int_in;
        int_in=0;
    }
    }
    struct massless_coll_return_vals mcrv;
    mcrv.f=(xf-xi)*(xf-xi)/4*int_val;
    mcrv.f_moment=(xf-xi)*(xf-xi)/4*int_val_moment;
    return mcrv;
}

void fill_massless_coll_table(struct massless_coll_table* mct,double y)
{
    double x0=mct->x0;
    double x1=mct->x1;
    int Nx=mct->Nx;
    mct->hx=(x1-x0)/(Nx-1.0);
    double hx=mct->hx;

    double c0=mct->c0;
    double c1=mct->c1;
    int Nc=mct->Nc;
    mct->hc=(c1-c0)/(Nc-1.0);
    double hc=mct->hc;

    mct->table=malloc(Nx*Nc*sizeof(double));
    mct->table_dx=malloc(Nx*Nc*sizeof(double));
    mct->table_dc=malloc(Nx*Nc*sizeof(double));
    mct->table_dxc=malloc(Nx*Nc*sizeof(double));

    mct->table_moment=malloc(Nx*Nc*sizeof(double));
    mct->table_dx_moment=malloc(Nx*Nc*sizeof(double));
    mct->table_dc_moment=malloc(Nx*Nc*sizeof(double));
    mct->table_dxc_moment=malloc(Nx*Nc*sizeof(double));

    double m=mct->m;
    double (*f)(double,double)=mct->cs;
    struct massless_coll_vals mcv;
    #pragma omp parallel for collapse(2) schedule(dynamic)
    for(int i=0;i<Nx;++i)
    { 
        for(int j=0;j<Nc;++j)
        {
        mcv=get_massless_coll_vals(x0+i*hx,c0+j*hc,y,f);
        mct->table[i+j*Nx]=mcv.f;
        mct->table_dx[i+j*Nx]=mcv.f_dx;
        mct->table_dc[i+j*Nx]=mcv.f_dc;
        mct->table_dxc[i+j*Nx]=mcv.f_dxc;

        mct->table_moment[i+j*Nx]=mcv.f_moment;
        mct->table_dx_moment[i+j*Nx]=mcv.f_dx_moment;
        mct->table_dc_moment[i+j*Nx]=mcv.f_dc_moment;
        mct->table_dxc_moment[i+j*Nx]=mcv.f_dxc_moment;
        }
    }
}
struct massless_coll_return_vals interp_massless_col(double x,double c,double y, struct massless_coll_table* mct)
{

double x0=mct->x0;
double hx=mct->hx;
double c0=mct->c0;
double hc=mct->hc;
int Nx=mct->Nx;
int Nc=mct->Nc;

int idx_x=floor((x-x0)/hx);
int idx_c=floor((c-c0)/hc);


if(idx_x<Nx-1 && idx_x>=0 && idx_c<Nc-1 && idx_c>=0)
{
    //printf("Interpolating\n");
    struct massless_coll_return_vals mcrv;
    double pts[4]={mct->table[idx_x+idx_c*Nx],mct->table[idx_x+1+idx_c*Nx],mct->table[idx_x+(idx_c+1)*Nx],mct->table[idx_x+1+(idx_c+1)*Nx]};
    double dx_pts[4]={mct->table_dx[idx_x+idx_c*Nx],mct->table_dx[idx_x+1+idx_c*Nx],mct->table_dx[idx_x+(idx_c+1)*Nx],mct->table_dx[idx_x+1+(idx_c+1)*Nx]};
    double dc_pts[4]={mct->table_dc[idx_x+idx_c*Nx],mct->table_dc[idx_x+1+idx_c*Nx],mct->table_dc[idx_x+(idx_c+1)*Nx],mct->table_dc[idx_x+1+(idx_c+1)*Nx]};
    double dxc_pts[4]={mct->table_dxc[idx_x+idx_c*Nx],mct->table_dxc[idx_x+1+idx_c*Nx],mct->table_dxc[idx_x+(idx_c+1)*Nx],mct->table_dxc[idx_x+1+(idx_c+1)*Nx]};
    mcrv.f=two_point_2d_hermite_interp(x,c,x0+idx_x*hx,x0+(idx_x+1)*hx,c0+idx_c*hc,c0+(idx_c+1)*hc,pts,dx_pts,dc_pts,dxc_pts);
    double moment_pts[4]={mct->table_moment[idx_x+idx_c*Nx],mct->table_moment[idx_x+1+idx_c*Nx],mct->table_moment[idx_x+(idx_c+1)*Nx],mct->table_moment[idx_x+1+(idx_c+1)*Nx]};
    double dx_moment_pts[4]={mct->table_dx_moment[idx_x+idx_c*Nx],mct->table_dx_moment[idx_x+1+idx_c*Nx],mct->table_dx_moment[idx_x+(idx_c+1)*Nx],mct->table_dx_moment[idx_x+1+(idx_c+1)*Nx]};
    double dc_moment_pts[4]={mct->table_dc_moment[idx_x+idx_c*Nx],mct->table_dc_moment[idx_x+1+idx_c*Nx],mct->table_dc_moment[idx_x+(idx_c+1)*Nx],mct->table_dc_moment[idx_x+1+(idx_c+1)*Nx]};
    double dxc_moment_pts[4]={mct->table_dxc_moment[idx_x+idx_c*Nx],mct->table_dxc_moment[idx_x+1+idx_c*Nx],mct->table_dxc_moment[idx_x+(idx_c+1)*Nx],mct->table_dxc_moment[idx_x+1+(idx_c+1)*Nx]};
    mcrv.f_moment=two_point_2d_hermite_interp(x,c,x0+idx_x*hx,x0+(idx_x+1)*hx,c0+idx_c*hc,c0+(idx_c+1)*hc,moment_pts,dx_moment_pts,dc_moment_pts,dxc_moment_pts);
    return mcrv;
}
    return ts_quad_massless_coll(x,c,y,mct->cs);//*T*T*T/(m*m*M_PI*M_PI);
}

void free_massless_coll_table(struct massless_coll_table* mct)
{
    free(mct->table);
    free(mct->table_dx);
    free(mct->table_dc);
    free(mct->table_dxc);

    free(mct->table_moment);
    free(mct->table_dx_moment);
    free(mct->table_dc_moment);
    free(mct->table_dxc_moment);
}








struct massive_coll_vals get_collision_vals(double m, double b, double c)
{
    int n=TS_NUM_PTS;
    double int_val, moment_int_val, E;
    int_val=0.0;
    moment_int_val=0;

double intval_in=0;
double moment_intval_in=0;
double test_term=0;
double xf=10*b;
double xi=1;

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

double db_term=0;
double dc_term=0;

double db_int_in=0;
double db_moment_in=0;
double dc_int_in=0;
double dc_moment_in=0;

double intval_db=0;
double intval_dc=0;
double moment_intval_db=0;
double moment_intval_dc=0;

double arg_plus=0;
double arg_minus=0;
double dbc_term1=0;
double dbc_term2=0;
double dbc_term3=0;
double dbcterm_tot=0;
double dbc_int_in=0;
double dbc_moment_in=0;
double moment_intval_dbc=0;
double intval_dbc=0;
    for(int k=-n;k<n;++k)
    {
        E=( xf*(ts_quad_fixed_nodes[k+n]+1) - xi*(ts_quad_fixed_nodes[k+n]-1) )/2;
        E_sqrt_term=sqrt(E*E-1);
        cs_term=c_vv_pp(4*E*E*(m*m),1);
        zf=40/E;
        if(zf<1) continue;
         for(int j=-n;j<n;++j)
        {
        v=( zf*(ts_quad_fixed_nodes[j+n]+1) - zi*(ts_quad_fixed_nodes[j+n]-1) )/2;
        v_sqrt_term=sqrt(v*v-1);
        arg_plus=E*v+E_sqrt_term*v_sqrt_term;
        arg_minus=E*v-E_sqrt_term*v_sqrt_term;
        sinh_term_plus=sinh( (b*arg_plus-c)/2 );
        coth_term_plus=coth( (b*arg_plus-c)/2 );
        sinh_term_minus=sinh( (b*arg_minus-c)/2 );
        coth_term_minus=coth( (b*arg_minus-c)/2 );
        exp_term=exp(2*(b*E*v-c));
        fd_term=1/(exp_term-1);
        lnsinh_term=log(fabs(sinh_term_plus/sinh_term_minus));



        test_term=ts_quad_fixed_weights[j+n]*fd_term*lnsinh_term;
        db_term=ts_quad_fixed_weights[j+n]*fd_term*(-2*E*v*exp_term*fd_term*lnsinh_term+0.5*(arg_plus*coth_term_plus-arg_minus*coth_term_minus));
        dc_term=ts_quad_fixed_weights[j+n]*fd_term*(2*exp_term*fd_term*lnsinh_term-0.5*(coth_term_plus-coth_term_minus));

        dbc_term1=-4*fd_term*fd_term*fd_term*exp_term*E*v*(exp_term+1)*lnsinh_term;
       dbc_term2=0.25*fd_term*(arg_plus/(sinh_term_plus*sinh_term_plus)-arg_minus/(sinh_term_minus*sinh_term_minus));
         dbc_term3=exp_term*fd_term*fd_term*((2*E*v+E_sqrt_term*v_sqrt_term)*coth_term_plus-(2*E*v-E_sqrt_term*v_sqrt_term)*coth_term_minus);
        dbcterm_tot=ts_quad_fixed_weights[j+n]*(dbc_term1+dbc_term2+dbc_term3);
         if(isnan(test_term) ||isnan(db_term) || isnan(dc_term) || isnan(dbcterm_tot))
        {
            break;
        }
            intval_in+=test_term;
            moment_intval_in+=test_term*v;

            db_int_in+=db_term;
            db_moment_in+=db_term*v;

            dc_int_in+=dc_term;
            dc_moment_in+=dc_term*v;

            dbc_int_in+=dbcterm_tot;
            dbc_moment_in+=dbcterm_tot*v;
       
        }
        int_val+=ts_quad_fixed_weights[k+n]*(zf-zi)/2*intval_in*E*E*cs_term;
        moment_int_val+=ts_quad_fixed_weights[k+n]*(zf-zi)/2*moment_intval_in*E*E*E*cs_term;

        intval_db+=ts_quad_fixed_weights[k+n]*(zf-zi)/2*db_int_in*E*E*cs_term;
        intval_dc+=ts_quad_fixed_weights[k+n]*(zf-zi)/2*dc_int_in*E*E*cs_term;
        intval_dbc+=ts_quad_fixed_weights[k+n]*(zf-zi)/2*dbc_int_in*E*E*cs_term;

        moment_intval_db+=ts_quad_fixed_weights[k+n]*(zf-zi)/2*db_moment_in*E*E*E*cs_term;
        moment_intval_dc+=ts_quad_fixed_weights[k+n]*(zf-zi)/2*dc_moment_in*E*E*E*cs_term;
        moment_intval_dbc+=ts_quad_fixed_weights[k+n]*(zf-zi)/2*dbc_moment_in*E*E*E*cs_term;


        intval_in=0;
        moment_intval_in=0;
        db_int_in=0;
        db_moment_in=0;
        dc_int_in=0;
        dc_moment_in=0;
        dbc_int_in=0;
        dbc_moment_in=0;
    }
  struct massive_coll_vals vals;
  vals.f=(xf-xi)/2*int_val;//*16*T4/(2*TWO_PI_4);
  vals.f_moment=(xf-xi)/2*moment_int_val;//*16*T4*T/(2*TWO_PI_4);
vals.f_db=(xf-xi)/2*intval_db;//*16*T4/(2*TWO_PI_4*T);
vals.f_db_moment=(xf-xi)/2*moment_intval_db;//*16*T4/(2*TWO_PI_4);

vals.f_dc=(xf-xi)/2*intval_dc;//*16*T*T*T*T/(2*TWO_PI_4);
vals.f_dc_moment=(xf-xi)/2*moment_intval_dc;//*16*T*T*T*T*T/(2*TWO_PI_4);
vals.f_dbc=(xf-xi)/2*intval_dbc;//*16*T*T*T*T/(2*TWO_PI_4);
vals.f_dbc_moment=(xf-xi)/2*moment_intval_dbc;//*16*T*T*T*T*T/(2*TWO_PI_4);
    return vals;
}


struct massive_coll_vals get_collision_vals2(double m, double b, double c)
{
    int n=TS_NUM_PTS;
    double int_val, moment_int_val, E;
    int_val=0.0;
    moment_int_val=0;

double intval_in=0;
double moment_intval_in=0;
double test_term=0;
double xf=70;
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

double db_term=0;
double dc_term=0;

double db_int_in=0;
double db_moment_in=0;
double dc_int_in=0;
double dc_moment_in=0;

double intval_db=0;
double intval_dc=0;
double moment_intval_db=0;
double moment_intval_dc=0;

double arg_plus=0;
double arg_minus=0;
double dbc_term1=0;
double dbc_term2=0;
double dbc_term3=0;
double dbcterm_tot=0;
double dbc_int_in=0;
double dbc_moment_in=0;
double moment_intval_dbc=0;
double intval_dbc=0;
    for(int k=-n;k<n;++k)
    {
        E=( xf*(ts_quad_fixed_nodes[k+n]+1) - xi*(ts_quad_fixed_nodes[k+n]-1) )/2;
        E_sqrt_term=sqrt(E*E-b*b);
        cs_term=c_pp_vv(4*E*E*m*m/(b*b),1);
        zf=20/E;
        if(zf<1) continue;
         for(int j=-n;j<n;++j)
        {
        v=( zf*(ts_quad_fixed_nodes[j+n]+1) - zi*(ts_quad_fixed_nodes[j+n]-1) )/2;
        v_sqrt_term=sqrt(v*v-1);
        arg_plus=E*v+E_sqrt_term*v_sqrt_term;
        arg_minus=E*v-E_sqrt_term*v_sqrt_term;
        sinh_term_plus=sinh( (arg_plus-c)/2 );
        coth_term_plus=coth( (arg_plus-c)/2 );
        sinh_term_minus=sinh( (arg_minus-c)/2 );
        coth_term_minus=coth( (arg_minus-c)/2 );
        exp_term=exp(2*(E*v-c));
        fd_term=1/(exp_term-1);
        lnsinh_term=log(fabs(sinh_term_plus/sinh_term_minus));



        test_term=ts_quad_fixed_weights[j+n]*fd_term*lnsinh_term;
        db_term=ts_quad_fixed_weights[j+n]*fd_term*(-2*E*v*exp_term*fd_term*lnsinh_term+0.5*(arg_plus*coth_term_plus-arg_minus*coth_term_minus));
        dc_term=ts_quad_fixed_weights[j+n]*fd_term*(2*exp_term*fd_term*lnsinh_term-0.5*(coth_term_plus-coth_term_minus));

        dbc_term1=-4*fd_term*fd_term*fd_term*exp_term*E*v*(exp_term+1)*lnsinh_term;
       dbc_term2=0.25*fd_term*(arg_plus/(sinh_term_plus*sinh_term_plus)-arg_minus/(sinh_term_minus*sinh_term_minus));
         dbc_term3=exp_term*fd_term*fd_term*((2*E*v+E_sqrt_term*v_sqrt_term)*coth_term_plus-(2*E*v-E_sqrt_term*v_sqrt_term)*coth_term_minus);
        dbcterm_tot=ts_quad_fixed_weights[j+n]*(dbc_term1+dbc_term2+dbc_term3);
         if(isnan(test_term) ||isnan(db_term) || isnan(dc_term) || isnan(dbcterm_tot))
        {
            break;
        }
            intval_in+=test_term;
            moment_intval_in+=test_term*v;

            db_int_in+=db_term;
            db_moment_in+=db_term*v;

            dc_int_in+=dc_term;
            dc_moment_in+=dc_term*v;

            dbc_int_in+=dbcterm_tot;
            dbc_moment_in+=dbcterm_tot*v;
       
        }
        int_val+=ts_quad_fixed_weights[k+n]*(zf-zi)/2*intval_in*E*E*cs_term;
        moment_int_val+=ts_quad_fixed_weights[k+n]*(zf-zi)/2*moment_intval_in*E*E*E*cs_term;

        intval_db+=ts_quad_fixed_weights[k+n]*(zf-zi)/2*db_int_in*E*E*cs_term;
        intval_dc+=ts_quad_fixed_weights[k+n]*(zf-zi)/2*dc_int_in*E*E*cs_term;
        intval_dbc+=ts_quad_fixed_weights[k+n]*(zf-zi)/2*dbc_int_in*E*E*cs_term;

        moment_intval_db+=ts_quad_fixed_weights[k+n]*(zf-zi)/2*db_moment_in*E*E*E*cs_term;
        moment_intval_dc+=ts_quad_fixed_weights[k+n]*(zf-zi)/2*dc_moment_in*E*E*E*cs_term;
        moment_intval_dbc+=ts_quad_fixed_weights[k+n]*(zf-zi)/2*dbc_moment_in*E*E*E*cs_term;


        intval_in=0;
        moment_intval_in=0;
        db_int_in=0;
        db_moment_in=0;
        dc_int_in=0;
        dc_moment_in=0;
        dbc_int_in=0;
        dbc_moment_in=0;
    }
  struct massive_coll_vals vals;
  double bm2=1/(b*b);
  double bm3=1/(b*b*b);
  double bm4=bm2*bm2;
  vals.f=(xf-xi)/2*int_val*bm3;//*16*T4/(2*TWO_PI_4);
  vals.f_moment=(xf-xi)/2*moment_int_val*bm4;//*16*T4*T/(2*TWO_PI_4);
vals.f_db=(xf-xi)/2*intval_db*bm4;//*16*T4/(2*TWO_PI_4*T);
vals.f_db_moment=(xf-xi)/2*moment_intval_db*bm3*bm2;//*16*T4/(2*TWO_PI_4);

vals.f_dc=(xf-xi)/2*intval_dc*bm3;//*16*T*T*T*T/(2*TWO_PI_4);
vals.f_dc_moment=(xf-xi)/2*moment_intval_dc*bm4;//*16*T*T*T*T*T/(2*TWO_PI_4);
vals.f_dbc=(xf-xi)/2*intval_dbc*bm4;//*16*T*T*T*T/(2*TWO_PI_4);
vals.f_dbc_moment=(xf-xi)/2*moment_intval_dbc*bm3*bm2;//*16*T*T*T*T*T/(2*TWO_PI_4);
    return vals;
}



void fill_massive_coll_table(struct massive_coll_table* mct,double y)
{
    double b0=mct->b0;
    double b1=mct->b1;
    int Nb=mct->Nb;
    mct->hb=(b1-b0)/(Nb-1.0);
    double hb=mct->hb;

    double c0=mct->c0;
    double c1=mct->c1;
    int Nc=mct->Nc;
    mct->hc=(c1-c0)/(Nc-1.0);
    double hc=mct->hc;

    mct->table=malloc(Nb*Nc*sizeof(double));
    mct->table_db=malloc(Nb*Nc*sizeof(double));
    mct->table_dc=malloc(Nb*Nc*sizeof(double));
    mct->table_dbc=malloc(Nb*Nc*sizeof(double));
    mct->table_moment=malloc(Nb*Nc*sizeof(double));
    mct->table_db_moment=malloc(Nb*Nc*sizeof(double));
    mct->table_dc_moment=malloc(Nb*Nc*sizeof(double));
    mct->table_dbc_moment=malloc(Nb*Nc*sizeof(double));

    double m=mct->m;
    struct massive_coll_vals mcv;
    #pragma omp parallel for collapse(2) schedule(dynamic)
    for(int i=0;i<Nb;++i)
    { 
        for(int j=0;j<Nc;++j)
        {
        mcv=get_collision_vals2(m,b0+i*hb,c0+j*hc);
        mct->table[i+j*Nb]=mcv.f;
        mct->table_db[i+j*Nb]=mcv.f_db;
        mct->table_dc[i+j*Nb]=mcv.f_dc;
        mct->table_dbc[i+j*Nb]=mcv.f_dbc;
        mct->table_moment[i+j*Nb]=mcv.f_moment;
        mct->table_db_moment[i+j*Nb]=mcv.f_db_moment;
        mct->table_dc_moment[i+j*Nb]=mcv.f_dc_moment;
        mct->table_dbc_moment[i+j*Nb]=mcv.f_dbc_moment;
        }
    }
}


/*struct massive_coll_return_vals ts_quad_massive_coll(double m, double b, double c,double y,double (*f)(double,double))
{
    int n=TS_NUM_PTS;
    double int_val, moment_int_val, E;
    int_val=0.0;
    moment_int_val=0;

double intval_in=0;
double moment_intval_in=0;
double test_term=0;
double xf=30;

}*/

void free_massive_coll_table(struct massive_coll_table* mct)
{
    free(mct->table);
    free(mct->table_db);
    free(mct->table_dc);
    free(mct->table_dbc);
    free(mct->table_moment);
    free(mct->table_db_moment);
    free(mct->table_dc_moment);
    free(mct->table_dbc_moment);
}


struct collision_vals interp_massive_colls(double m,double b,double c,struct massive_coll_table* mct)
{
double b0=mct->b0;
double hb=mct->hb;
double c0=mct->c0;
double hc=mct->hc;
int Nb=mct->Nb;
int Nc=mct->Nc;

int idx_b=floor((b-b0)/hb);
int idx_c=floor((c-c0)/hc);

struct collision_vals mcrv;

if(idx_b<Nb-1 && idx_b>=0 && idx_c<Nc-1 && idx_c>=0)
{
   // printf("Interpolating\n");
    double pts[4]={mct->table[idx_b+idx_c*Nb],mct->table[idx_b+1+idx_c*Nb],mct->table[idx_b+(idx_c+1)*Nb],mct->table[idx_b+1+(idx_c+1)*Nb]};
    double db_pts[4]={mct->table_db[idx_b+idx_c*Nb],mct->table_db[idx_b+1+idx_c*Nb],mct->table_db[idx_b+(idx_c+1)*Nb],mct->table_db[idx_b+1+(idx_c+1)*Nb]};
    double dc_pts[4]={mct->table_dc[idx_b+idx_c*Nb],mct->table_dc[idx_b+1+idx_c*Nb],mct->table_dc[idx_b+(idx_c+1)*Nb],mct->table_dc[idx_b+1+(idx_c+1)*Nb]};
    double dbc_pts[4]={mct->table_dbc[idx_b+idx_c*Nb],mct->table_dbc[idx_b+1+idx_c*Nb],mct->table_dbc[idx_b+(idx_c+1)*Nb],mct->table_dbc[idx_b+1+(idx_c+1)*Nb]};
    mcrv.coll=two_point_2d_hermite_interp(b,c,b0+idx_b*hb,b0+(idx_b+1)*hb,c0+idx_c*hc,c0+(idx_c+1)*hc,pts,db_pts,dc_pts,dbc_pts)*8*m*m*m*m/(TWO_PI_4*b);
    double mpts[4]={mct->table_moment[idx_b+idx_c*Nb],mct->table_moment[idx_b+1+idx_c*Nb],mct->table_moment[idx_b+(idx_c+1)*Nb],mct->table_moment[idx_b+1+(idx_c+1)*Nb]};
    double mdb_pts[4]={mct->table_db_moment[idx_b+idx_c*Nb],mct->table_db_moment[idx_b+1+idx_c*Nb],mct->table_db_moment[idx_b+(idx_c+1)*Nb],mct->table_db_moment[idx_b+1+(idx_c+1)*Nb]};
    double mdc_pts[4]={mct->table_dc_moment[idx_b+idx_c*Nb],mct->table_dc_moment[idx_b+1+idx_c*Nb],mct->table_dc_moment[idx_b+(idx_c+1)*Nb],mct->table_dc_moment[idx_b+1+(idx_c+1)*Nb]};
    double mdbc_pts[4]={mct->table_dbc_moment[idx_b+idx_c*Nb],mct->table_dbc_moment[idx_b+1+idx_c*Nb],mct->table_dbc_moment[idx_b+(idx_c+1)*Nb],mct->table_dbc_moment[idx_b+1+(idx_c+1)*Nb]};
    mcrv.moment=two_point_2d_hermite_interp(b,c,b0+idx_b*hb,b0+(idx_b+1)*hb,c0+idx_c*hc,c0+(idx_c+1)*hc,mpts,mdb_pts,mdc_pts,mdbc_pts)*8*m*m*m*m*m/(TWO_PI_4*b);
    return mcrv;
}
return Moment1_omp(m,m/b,c);
}















struct derivs integrated_de_omp_colls(double xi, double xf, struct fixed_params fixed_params, struct sp_params sp_params,struct interpolation_tables* tables)
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
struct polylogs* plgs=tables->plgs;

double T3=T*T*T;
double Ts3=Ts*Ts*Ts;
double Tphi3=Tphi*Tphi*Tphi;
double T6=T3*T3;
double Tphi6=Tphi3*Tphi3;
double Ts6=Tphi3*Tphi3;


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
double M_as_p = 8  *y* (m_phi2 - m_s * m_s); //reduce by y2
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

struct massless_coll_return_vals C_aa_ss;
struct massless_coll_return_vals C_aa_pp;
struct massless_coll_return_vals C_ss_aa;
struct massless_coll_return_vals C_ss_pp;
struct collision_vals C_pp_vv;
double collision_factor=y*y*y*8/(TWO_PI_4*m_phi2); //reduce by y2
double arg_factor=4/(m_phi2);
#pragma omp taskgroup task_reduction(+:C_as_p_destroy) task_reduction(+:C_as_p_create) task_reduction(+:C_p_as_create) task_reduction(+:C_p_as_destroy) task_reduction(+:C_as_p_num) 
{ 
    /*
        #pragma omp task
    {
        C_aa_ss=interp_massless_col(arg_factor*T*T,c_s,y,tables->aa_ss_table);
        C_aa_ss.f*=T6*collision_factor;
        C_aa_ss.f_moment*=T*T6*collision_factor;
    }
    #pragma omp task
    {
        C_aa_pp=interp_massless_col(arg_factor*T*T,c_phi,y,tables->aa_pp_table);
        C_aa_pp.f*=T6*collision_factor;
        C_aa_pp.f_moment*=T*T6*collision_factor;
    }
    #pragma omp task
    {
        C_ss_aa=interp_massless_col(arg_factor*Ts*Ts,c_s,y,tables->ss_aa_table);
        C_ss_aa.f*=Ts6*collision_factor;
        C_ss_aa.f_moment*=Ts*Ts6*collision_factor;
    }
    #pragma omp task
    {
        C_ss_pp=interp_massless_col(arg_factor*Ts*Ts,c_phi,y,tables->ss_pp_table);
        C_ss_pp.f*=Ts6*collision_factor;
        C_ss_pp.f_moment*=Ts*Ts6*collision_factor;
    }
    #pragma omp task
    {
        C_pp_vv=interp_massive_colls(m_phi,m_phi/Tphi,c_phi,tables->pp_vv_table);
        C_pp_vv.coll*=y*y*y*y;
        C_pp_vv.moment*=y*y*y*y;
    }*/
    #pragma omp taskloop grainsize(10) untied in_reduction(+:C_as_p_destroy) in_reduction(+:C_as_p_create) in_reduction(+:C_p_as_create) in_reduction(+:C_p_as_destroy) in_reduction(+:C_as_p_num) 
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
       C_aa_ss=interp_massless_col(arg_factor*T*T,c_s,y,tables->aa_ss_table);
        C_aa_ss.f*=T6*collision_factor;
        C_aa_ss.f_moment*=T*T6*collision_factor;

        C_aa_pp=interp_massless_col(arg_factor*T*T,c_phi,y,tables->aa_pp_table);
        C_aa_pp.f*=T6*collision_factor;
        C_aa_pp.f_moment*=T*T6*collision_factor;

        C_ss_aa=interp_massless_col(arg_factor*Ts*Ts,c_s,y,tables->ss_aa_table);
        C_ss_aa.f*=Ts6*collision_factor;
        C_ss_aa.f_moment*=Ts*Ts6*collision_factor;

 C_ss_pp=interp_massless_col(arg_factor*Ts*Ts,c_phi,y,tables->ss_pp_table);
        C_ss_pp.f*=Ts6*collision_factor;
        C_ss_pp.f_moment*=Ts*Ts6*collision_factor;

C_pp_vv=interp_massive_colls(m_phi,m_phi/Tphi,c_phi,tables->pp_vv_table);
        C_pp_vv.coll*=y*y*y;
        C_pp_vv.moment*=y*y*y;
/*
printf("c_aa_ss: %.10e\n",C_aa_ss);
printf("c_aa_pp: %.10e\n",C_aa_pp);
printf("c_ss_aa: %.10e\n",C_ss_aa);
printf("c_ss_pp: %.10e\n",C_ss_pp);
printf("c_pp_vv: %.10e\n",C_pp_vv);*/

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

    derivs.C_n_s=C_as_p_num+C_aa_ss.f+C_pp_vv.coll-C_ss_pp.f-C_ss_aa.f;
    derivs.C_p_s=C_as_p_moment_tot+C_aa_ss.f_moment+C_pp_vv.moment-C_ss_pp.f_moment-C_ss_aa.f_moment;

    derivs.C_n_p=-C_as_p_num +C_ss_pp.f+C_aa_pp.f-2*C_pp_vv.coll;
    derivs.C_p_p=C_p_as_moment_tot+C_ss_pp.f_moment+C_aa_pp.f_moment-2*C_pp_vv.moment;

    return derivs;
}





struct final_derivs step_omp_coll(struct sp_params sp, struct fixed_params* fixed_params,struct interpolation_tables* tables)
{
double mphi=fixed_params->m_phi;
double T=sp.T;
struct derivs new_derivs=integrated_de_omp_colls(0,100,*fixed_params,sp,tables);
/*printf("C_n_s: %.10e\n",new_derivs.C_n_s);
printf("C_p_s: %.10e\n",new_derivs.C_p_s);
printf("C_n_p: %.10e\n",new_derivs.C_n_p);
printf("C_p_p: %.10e\n",new_derivs.C_p_p);*/
double Tp=sp.Tphi;
double cp=sp.c_phi;
double Ts=sp.Ts;
double cs=sp.c_s;
struct polylogs* plgs=tables->plgs;


struct n_p_vals2 npv=get_n_p_vals2(mphi,Ts, cs, Tp,  cp, plgs);

double In_p=npv.In_p;
double In_p_T=npv.In_p_T;
double In_p_c=npv.In_p_c;
double n_p=In_p;

double Ip_p=npv.Ip_p;
double Ip_p_T=npv.Ip_p_T;
double Ip_p_c=npv.Ip_p_c;
double p_p=Ip_p;


//double ap=(3/Tp*In_p-mphi/(Tp*Tp)*In_p_b);
//double bp=In_p_c;
//double Cp=(4/Tp*Ip_p-mphi/(Tp*Tp)*Ip_p_b);
//double dp=Ip_p_c;

double ap=In_p_T;
double bp=In_p_c;
double Cp=Ip_p_T;
double dp=Ip_p_c;
double detp=ap*dp-bp*Cp;




double n_p_T=new_derivs.C_n_p;
double p_p_T=new_derivs.C_p_p;
double Ap=n_p_T;
double Bp=p_p_T;


double as=npv.n_s_T;
double bs=npv.n_s_c;
double Cs=npv.p_s_T;
double ds=npv.p_s_c;
double dets=as*ds-bs*Cs;

double n_s=npv.n_s;
double p_s=npv.p_s;

double n_s_T=new_derivs.C_n_s;
double p_s_T=new_derivs.C_p_s;

double dTs=1/dets*(ds*n_s_T-bs*p_s_T);
double dcs=1/dets*(-Cs*n_s_T+as*p_s_T);

double A_p=n_p_T;
double B_p=p_p_T;
double dTp=1/detp*(dp*A_p-bp*B_p);
double dcp=1/detp*(-Cp*A_p+ap*B_p);
/*
printf("T: %.10e, T_s: %.10e, T_p: %.10e, c_s: %.10e, c_p: %.10e\n",T,Ts,Tp,cs,cp);
printf("n_s: %.10e\n",n_s);
printf("p_s: %.10e\n",p_s);
printf("n_p: %.10e\n",In_p);
printf("p_p: %.10e\n",Ip_p);
printf("n_s_c: %.10e\n",npv.n_s_c);
printf("p_s_c: %.10e\n",npv.p_s_c);
printf("n_s_T: %.10e\n",npv.n_s_T);
printf("p_s_T: %.10e\n",npv.p_s_T);
printf("In_p: %.10e\n",In_p);
printf("In_p_b: %.10e\n",In_p_c);
printf("In_p_T: %.10e\n",In_p_T);
printf("Ip_p: %.10e\n",Ip_p);
printf("Ip_p_b: %.10e\n",Ip_p_c);
printf("Ip_p_T: %.10e\n",Ip_p_T);*/

//printf("a*dT_s+b*dc_s: %.10e, C_n_s: %.10e\n",as*dTs+bs*dcs,new_derivs.C_n_s);
//printf("c*dT_s+d*dc_s: %.10e, C_p_s: %.10e\n",Cs*dTs+ds*dcs,new_derivs.C_p_s);
//printf("a*dT_p+b*dc_p: %.10e, C_n_p: %.10e\n",ap*dTp+bp*dcp,new_derivs.C_n_p);
//printf("c*dT_p+d*dc_p: %.10e, C_p_p: %.10e\n",Cp*dTp+dp*dcp,new_derivs.C_p_p);

struct final_derivs final_derivs;

final_derivs.dTs=dTs;
final_derivs.dTp=dTp;
final_derivs.dcs=dcs;
final_derivs.dcp=dcp;


return final_derivs;

}