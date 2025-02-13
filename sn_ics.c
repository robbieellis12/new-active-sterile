#include "sn_ics.h"


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