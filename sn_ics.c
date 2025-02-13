#include "sn_ics.h"


struct massless_coll_vals get_massless_coll_vals(double x, double c,double y,double (*f)(double,double))
{
    int n=TS_NUM_PTS;
 double  u,p1;

    double int_val=0.0;
    double int_val_dx=0.0;
    double int_val_dc=0.0;
    double int_val_dxc=0.0;

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
    return vals;
}

double ts_quad_massless_coll(double x, double c,double y,double (*f)(double,double))
{
    int n=TS_NUM_PTS;
    double  u,p1;

    double int_val=0.0;
    double int_in=0;

    double xf=100;
    double xi=0;

    double cs_term;
    for(int k=-n;k<n;++k)
    {
        p1=( xf*(ts_quad_fixed_nodes[k+n]+1) - xi*(ts_quad_fixed_nodes[k+n]-1) )/2;
        for(int l=-n;l<n;++l)
        {
            u=( xf*(ts_quad_fixed_nodes[l+n]+1) - xi*(ts_quad_fixed_nodes[l+n]-1) )/2;

            int_in+=ts_quad_fixed_weights[l+n]*u*log(1+exp(c-u))*f(x*p1*u,y);
        }
        int_val+=ts_quad_fixed_weights[k+n]*p1*p1*1/(exp(p1-c)+1)*int_in;
        int_in=0;
    }
    return (xf-xi)*(xf-xi)/4*int_val;
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
        }
    }
}
double interp_massless_col(double x,double c,struct massless_coll_table* mct)
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
    printf("Interpolating\n");
    double pts[4]={mct->table[idx_x+idx_c*Nx],mct->table[idx_x+1+idx_c*Nx],mct->table[idx_x+(idx_c+1)*Nx],mct->table[idx_x+1+(idx_c+1)*Nx]};
    double dx_pts[4]={mct->table_dx[idx_x+idx_c*Nx],mct->table_dx[idx_x+1+idx_c*Nx],mct->table_dx[idx_x+(idx_c+1)*Nx],mct->table_dx[idx_x+1+(idx_c+1)*Nx]};
    double dc_pts[4]={mct->table_dc[idx_x+idx_c*Nx],mct->table_dc[idx_x+1+idx_c*Nx],mct->table_dc[idx_x+(idx_c+1)*Nx],mct->table_dc[idx_x+1+(idx_c+1)*Nx]};
    double dxc_pts[4]={mct->table_dxc[idx_x+idx_c*Nx],mct->table_dxc[idx_x+1+idx_c*Nx],mct->table_dxc[idx_x+(idx_c+1)*Nx],mct->table_dxc[idx_x+1+(idx_c+1)*Nx]};
    return two_point_2d_hermite_interp(x,c,x0+idx_x*hx,x0+(idx_x+1)*hx,c0+idx_c*hc,c0+(idx_c+1)*hc,pts,dx_pts,dc_pts,dxc_pts);
}
    return ts_quad_massless_coll(x,c,mct->m,mct->cs);//*T*T*T/(m*m*M_PI*M_PI);
}

void free_massless_coll_table(struct massless_coll_table* mct)
{
    free(mct->table);
    free(mct->table_dx);
    free(mct->table_dc);
    free(mct->table_dxc);
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
    printf("Interpolating\n");
    double pts[4]={mct->table[idx_b+idx_c*Nb],mct->table[idx_b+1+idx_c*Nb],mct->table[idx_b+(idx_c+1)*Nb],mct->table[idx_b+1+(idx_c+1)*Nb]};
    double db_pts[4]={mct->table_db[idx_b+idx_c*Nb],mct->table_db[idx_b+1+idx_c*Nb],mct->table_db[idx_b+(idx_c+1)*Nb],mct->table_db[idx_b+1+(idx_c+1)*Nb]};
    double dc_pts[4]={mct->table_dc[idx_b+idx_c*Nb],mct->table_dc[idx_b+1+idx_c*Nb],mct->table_dc[idx_b+(idx_c+1)*Nb],mct->table_dc[idx_b+1+(idx_c+1)*Nb]};
    double dbc_pts[4]={mct->table_dbc[idx_b+idx_c*Nb],mct->table_dbc[idx_b+1+idx_c*Nb],mct->table_dbc[idx_b+(idx_c+1)*Nb],mct->table_dbc[idx_b+1+(idx_c+1)*Nb]};
    mcrv.coll=two_point_2d_hermite_interp(b,c,b0+idx_b*hb,b0+(idx_b+1)*hb,c0+idx_c*hc,c0+(idx_c+1)*hc,pts,db_pts,dc_pts,dbc_pts)*8/(TWO_PI_4*m*m*b);
    double mpts[4]={mct->table_moment[idx_b+idx_c*Nb],mct->table_moment[idx_b+1+idx_c*Nb],mct->table_moment[idx_b+(idx_c+1)*Nb],mct->table_moment[idx_b+1+(idx_c+1)*Nb]};
    double mdb_pts[4]={mct->table_db_moment[idx_b+idx_c*Nb],mct->table_db_moment[idx_b+1+idx_c*Nb],mct->table_db_moment[idx_b+(idx_c+1)*Nb],mct->table_db_moment[idx_b+1+(idx_c+1)*Nb]};
    double mdc_pts[4]={mct->table_dc_moment[idx_b+idx_c*Nb],mct->table_dc_moment[idx_b+1+idx_c*Nb],mct->table_dc_moment[idx_b+(idx_c+1)*Nb],mct->table_dc_moment[idx_b+1+(idx_c+1)*Nb]};
    double mdbc_pts[4]={mct->table_dbc_moment[idx_b+idx_c*Nb],mct->table_dbc_moment[idx_b+1+idx_c*Nb],mct->table_dbc_moment[idx_b+(idx_c+1)*Nb],mct->table_dbc_moment[idx_b+1+(idx_c+1)*Nb]};
    mcrv.moment=two_point_2d_hermite_interp(b,c,b0+idx_b*hb,b0+(idx_b+1)*hb,c0+idx_c*hc,c0+(idx_c+1)*hc,mpts,mdb_pts,mdc_pts,mdbc_pts)*8/(TWO_PI_4*m*m*m*b);
    return mcrv;
}
return Moment1_omp(m,m/b,c);
}
