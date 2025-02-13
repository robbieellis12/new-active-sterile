
#include "sn_math.h"

double coth(double x)
{
    return 1/tanh(x);
}

void read_paired_list(char* filename,double** xdat, double** ydat,int* len)
{
    FILE* file;
    file=fopen(filename,"r");
    if(file==NULL)
    {
        printf("File not found\n");
        return;
    }
    int idx=0;
    double temp1;
    double temp2;
    while(fscanf(file, "%lf, %lf",&temp1,&temp2)==2)
    {
        ++idx;
    }
    fclose(file);
    printf("File length: %d\n",idx);
    *xdat=(double*) malloc(sizeof(double)*idx);
    *ydat=(double*) malloc(sizeof(double)*idx);
    *len=idx;
    idx=0;
    file=fopen(filename,"r");
    while(fscanf(file, "%lf, %lf\n",&temp1,&temp2)!=EOF)
    {
        ++idx;
        (*xdat)[idx]=temp1;
        (*ydat)[idx]=temp2;
    }
    printf("File length: %d\n",idx);
}
double interpolate(double* x,double*y,int L,double x_0)
{
    int firstGreater=-1;
    for(int i=0;i<L;++i)
    {
        if(x[i]>x_0)
        {
            firstGreater=i;
            break;
        }
    }
    if(firstGreater==0)
    {
        return y[0];
    }
    if(firstGreater==-1)
    {
        return y[L-1];
    }
    return (y[firstGreater]-y[firstGreater-1])/(x[firstGreater]-x[firstGreater-1])*(x_0-x[firstGreater])+y[firstGreater];
}

void ts_quad_info(double* weights_arr, double* nodes_arr,double h, int n)
{
    for(int k=-n;k<n;++k)
    {
        nodes_arr[k+n]=tanh(M_PI/2*sinh(k*h));
        weights_arr[k+n]=h/2*M_PI*cosh(k*h)/( cosh(M_PI/2*sinh(k*h)) * cosh(M_PI/2*sinh(k*h)) );
    }
}

double ts_quad(double (*f)(double, void*),void* params,double a, double b, double h, int n)
{
    double int_val, weight, node, u;
    int_val=0.0;
    for(int k=-n;k<n;++k)
    {
        u=( b*(ts_quad_fixed_nodes[k+n]+1) - a*(ts_quad_fixed_nodes[k+n]-1) )/2;
        int_val+=ts_quad_fixed_weights[k+n]*(*f)(u,params);
        
    }
   
    return (b-a)/2*int_val;
}


//Define functions used for Hermite interpolation

double L1(double x0, double x1, double x2)
{
    double dx21=x2-x1;
    double dx01=x0-x1;
    double dx02=x0-x2;
    
    return dx02*dx02/  (dx21*dx21) * ( 1 + 2 * dx01/dx21);
}

double L1p(double x0, double x1, double x2)
{
    double dx21=x2-x1;
    double dx01=x0-x1;
    double dx02=x0-x2;
    
    return 2*dx02*dx02/(dx21*dx21*dx21) + 2*dx02/(dx21*dx21)*(1+2*dx01/dx21);
}

double L2(double x0, double x1, double x2)
{
    double dx21=x2-x1;
    double dx01=x0-x1;
    double dx02=x0-x2;
    
    return dx01*dx01/  (dx21*dx21) * ( 1 - 2 * dx02/dx21);
}

double L2p(double x0, double x1, double x2)
{
    double dx21=x2-x1;
    double dx01=x0-x1;
    double dx02=x0-x2;
    
    return 2*dx01/(dx21*dx21)*( 1-2*dx02/dx21 ) - 2*dx01*dx01/(dx21*dx21*dx21);
}

double K1(double x0, double x1, double x2)
{
    double dx21=x2-x1;
    double dx01=x0-x1;
    double dx02=x0-x2;
    
    return dx02*dx02/(dx21*dx21)*(dx01);
}

double K1p(double x0, double x1, double x2)
{
    double dx21=x2-x1;
    double dx01=x0-x1;
    double dx02=x0-x2;
    
    return 2*dx02*dx01/(dx21*dx21)+dx02*dx02/(dx21*dx21);
}

double K2(double x0, double x1, double x2)
{
    double dx21=x2-x1;
    double dx01=x0-x1;
    double dx02=x0-x2;
    
    return dx01*dx01/(dx21*dx21)*(dx02);
}

double K2p(double x0, double x1, double x2)
{
    double dx21=x2-x1;
    double dx01=x0-x1;
    double dx02=x0-x2;
    
    return 2*dx02*dx01/(dx21*dx21)+dx01*dx01/(dx21*dx21);
}

double two_point_hermite_interp(double x0, double x1, double x2, double y1, double yp1, double y2, double yp2)
{
    double dx21=x2-x1;
    double dx01=x0-x1;
    double dx02=x0-x2;
    
    double L1= dx02*dx02/  (dx21*dx21) * ( 1 + 2 * dx01/dx21);
    double L2= dx01*dx01/  (dx21*dx21) * ( 1 - 2 * dx02/dx21);
    
    double K1= dx02*dx02/(dx21*dx21)*(dx01);
    double K2= dx01*dx01/(dx21*dx21)*(dx02);
    
    return L1*y1+L2*y2+K1*yp1+K2*yp2;
}
//Generalization of hermite interpolation using first order partials+mixed partials
//for interpolating f(x,y) inside of a rectangle. Same logic as 1d Hermite interp.
//Points for f,fx,fy must be ordered as follows (in the xy plane):
// 3 4
// 1 2
double two_point_2d_hermite_interp(double x0, double y0, double x1, double x2, double y1, double y2, double* f, double* fx, double* fy, double* fxy)
{
    //Sum over bottom row, then top row from left to right
    double L1x=L1(x0,x1,x2);
    double L2x=L2(x0,x1,x2);
    double L1y=L1(y0,y1,y2);
    double L2y=L2(y0,y1,y2);
    
    double K1x=K1(x0,x1,x2);
    double K2x=K2(x0,x1,x2);
    double K1y=K1(y0,y1,y2);
    double K2y=K2(y0,y1,y2);
    
    return
    L1x*L1y*f[0] + L2x*L1y*f[1]
    + L1x*L2y*f[2] + L2x*L2y*f[3]
    
    + K1x*L1y*fx[0] + K2x*L1y*fx[1]
    + K1x*L2y*fx[2] + K2x*L2y*fx[3]
    
    + L1x*K1y*fy[0] + L2x*K1y*fy[1]
    + L1x*K2y*fy[2] + L2x*K2y*fy[3]
    
    + K1x*K1y*fxy[0] + K2x*K1y*fxy[1]
    + K1x*K2y*fxy[2] + K2x*K2y*fxy[3];
}

void fill_polylog_table(struct polylogs* plgs)
{
    double x0=plgs->x0;
    double x1=plgs->x1;
    
    plgs->h=(x1-x0)/plgs->N;
    int N=plgs->N;
    double h=plgs->h;
    plgs->N=N;
    
    plgs->x1=x0+h*N;
    x1=plgs->x1;
    
    plgs->li2=(double*) malloc(N*sizeof(double));
    plgs->dli2=(double*) malloc(N*sizeof(double));
    plgs->li3=(double*) malloc(N*sizeof(double));
    plgs->dli3=(double*) malloc(N*sizeof(double));
    plgs->li4=(double*) malloc(N*sizeof(double));
    plgs->dli4=(double*) malloc(N*sizeof(double));
    
   
    
    double li,dli;
#pragma omp parallel for private(li,dli)
    for(int i=0;i<N;++i)
    {
        ts_quad_li(x0+i*h, 2, &li, &dli, TS_STEP_SIZE, TS_NUM_PTS);
        plgs->li2[i]=li;
        plgs->dli2[i]=dli;
        ts_quad_li(x0+i*h, 3, &li, &dli, TS_STEP_SIZE, TS_NUM_PTS);
        plgs->li3[i]=li;
        plgs->dli3[i]=dli;
        ts_quad_li(x0+i*h, 4, &li, &dli, TS_STEP_SIZE, TS_NUM_PTS);
        plgs->li4[i]=li;
        plgs->dli4[i]=dli;
    }

    
}



double interp_Li2(double x,struct polylogs* plgs)
{
    double x0=plgs->x0;
    double x1=plgs->x1;
    int N=plgs->N;
    int idx=floor((x-x0)/plgs->h);
    
    if(idx<0)
    {
        return -log(-x)*log(-x)/2-M_PI*M_PI/6-interp_Li2(1/x,plgs);
    }
    else if(fabs(x)<0.05)
    {
        return x+x*x/4;
    }
    else if(N-1==idx)
    {
        return M_PI*M_PI/6;
    }
    else if(N-1<=idx)
    {
        return 0;
    }
    
    return two_point_hermite_interp(x, x0+plgs->h*idx, x0+plgs->h*(idx+1), plgs->li2[idx], plgs->dli2[idx], plgs->li2[idx+1], plgs->dli2[idx+1]);
}

double interp_Li3(double x,struct polylogs* plgs)
{
    double x0=plgs->x0;
    double x1=plgs->x1;
    int N=plgs->N;
    int idx=floor((x-x0)/plgs->h);
    
    if (idx<0) {
        return -M_PI*M_PI/6*log(-x)-log(-x)*log(-x)*log(-x)/6+interp_Li3(1/x,plgs);
    }
    else if(fabs(x)<0.15)
    {
        return x+x*x/8+x*x*x/27;
    }
    else if(idx>=N)
    {
        return 0;
    }
    
    return two_point_hermite_interp(x, x0+plgs->h*idx, x0+plgs->h*(idx+1), plgs->li3[idx], plgs->dli3[idx], plgs->li3[idx+1], plgs->dli3[idx+1])/2;
}

double interp_Li4(double x,struct polylogs* plgs)
{
    double x0=plgs->x0;
    double x1=plgs->x1;
    int N=plgs->N;
    int idx=floor((x-x0)/plgs->h);
    
    if (idx<0) {
        return -M_PI*M_PI*M_PI*M_PI*7/360 - M_PI*M_PI*log(-x)*log(-x)/12 -log(-x)*log(-x)*log(-x)*log(-x)/24- interp_Li4(1/x,plgs);
    }
    else if(fabs(x)<0.3)
    {
        return x+x*x/16+x*x*x/81;
    }
    else if(idx>=N)
    {
        return 0;
    }
    
    return two_point_hermite_interp(x, x0+plgs->h*idx, x0+plgs->h*(idx+1), plgs->li4[idx], plgs->dli4[idx], plgs->li4[idx+1], plgs->dli4[idx+1])/6;
}


void ts_quad_li(double x, double N, double* li,double* dli, double h, int n)
{
    
    double int_val, d_int_val, weight, node, u;
    int_val=0.0;
    d_int_val=0.0;
    double xf=100;
    for(int k=-n;k<n;++k)
    {
        u=( xf*(ts_quad_fixed_nodes[k+n]+1) )/2;
        int_val+=ts_quad_fixed_weights[k+n]*pow(u,N-1)/(1-exp(u)/x);
        d_int_val+=-ts_quad_fixed_weights[k+n]*pow(u,N-1)*exp(u)*1/(x-exp(u))*1/(x-exp(u));
    }
    *li= -xf/2*int_val;
    *dli= -xf/2*d_int_val;
}

void free_polylog_table(struct polylogs* plgs)
{
    free(plgs->li2);
    free(plgs->li3);
    free(plgs->li4);
    free(plgs->dli2);
    free(plgs->dli3);
    free(plgs->dli4);
}


double n_rho_fermi(double m, double T)
{
    int n=TS_NUM_PTS;
     double int_val, weight, E,b, u;
    int_val=0.0;
    double xf=100;
    b=m/T;
    for(int k=-n;k<n;++k)
    {
        u=( xf*(ts_quad_fixed_nodes[k+n]+1) )/2;
        E=sqrt(u*u+b*b);
        int_val+=ts_quad_fixed_weights[k+n]*u*u*E/(exp(E)+1);
    }
    return xf/2*int_val*2*T*T*T*T/(2*M_PI*M_PI); //the 2 is g=2
}


struct n_rho_explicit ts_quad_n_rho_explicit(double b,double c)
{
    int n=TS_NUM_PTS;
    double int_val_n, int_val_nb,int_val_nc,int_val_p,int_val_pb,int_val_pc, u;
    int_val_n=0.0;
    int_val_nb=0.0;
    int_val_nc=0.0;
    int_val_p=0.0;
    int_val_pb=0.0;
    int_val_pc=0.0;
    double E,expterm,fdterm;
    double integration_endp=150;
    for(int k=-n;k<n;++k)
    {
        u=( integration_endp*(ts_quad_fixed_nodes[k+n]+1))/2;
        E=sqrt(u*u+b*b);
        expterm=exp(E-c);
        double fdterm=1/(expterm-1);
        int_val_n+=ts_quad_fixed_weights[k+n]*u*u*fdterm;
        int_val_nb+=-ts_quad_fixed_weights[k+n]*u*u*fdterm*fdterm*expterm*b/E;
        int_val_nc+=ts_quad_fixed_weights[k+n]*u*u*fdterm*fdterm*expterm;
        
        
        int_val_p+=ts_quad_fixed_weights[k+n]*u*u*fdterm*E;
        int_val_pb+=ts_quad_fixed_weights[k+n]*b*u*u*(1/E-expterm*fdterm)*fdterm;
        int_val_pc+=ts_quad_fixed_weights[k+n]*u*u*fdterm*fdterm*E*expterm;
        
    }
    struct n_rho_explicit nrv;
    nrv.n=int_val_n*integration_endp/2;
    nrv.nb=int_val_nb*integration_endp/2;
    nrv.nc=int_val_nc*integration_endp/2;
    
    nrv.p=int_val_p*integration_endp/2;
    nrv.pb=int_val_pb*integration_endp/2;
    nrv.pc=int_val_pc*integration_endp/2;
    return nrv;
}

struct n_p_vals get_n_p_vals(double mphi, double Ts, double cs, double Tp, double cp,struct polylogs* plgs)
{
    struct n_p_vals npv;

    double li2val=-interp_Li2(-exp(cs),plgs);
    double li3val=-2*interp_Li3(-exp(cs),plgs);
    double li4val=-6*interp_Li4(-exp(cs),plgs);

    double s_coeff=Ts*Ts*Ts/(2*M_PI*M_PI);
    double p_coeff=Tp*Tp*Tp/(2*M_PI*M_PI);
    
    npv.n_s=s_coeff*li3val;
    npv.n_s_c=s_coeff*2*li2val;
    npv.n_s_T=3*Ts*Ts/(2*M_PI*M_PI)*li3val;
    npv.p_s=Ts*s_coeff*li4val;
    npv.p_s_c=Ts*s_coeff*3*li3val;
    npv.p_s_T=2*Ts*Ts*Ts/(M_PI*M_PI)*li4val;


    struct n_rho_explicit nre;
    nre=ts_quad_n_rho_explicit(mphi/Tp,cp);

    npv.In_p=nre.n;
    npv.In_p_c=nre.nc;
    npv.In_p_b=nre.nb;
    npv.Ip_p=nre.p;
    npv.Ip_p_c=nre.pc;
    npv.Ip_p_b=nre.pb;

    return npv;
}


struct n_rho_explicit2 ts_quad_n_rho_explicit2(double m, double T,double c)
{
    int n=TS_NUM_PTS;
    double b=m/T;
    double int_val_n, int_val_nT,int_val_nc,int_val_p,int_val_pT,int_val_pc, u;
    int_val_n=0.0;
    int_val_nT=0.0;
    int_val_nc=0.0;
    int_val_p=0.0;
    int_val_pT=0.0;
    int_val_pc=0.0;
    double other_fphi_deriv=0;
    double E,expterm,fdterm;
    double integration_endp=150;
    for(int k=-n;k<n;++k)
    {
        u=( integration_endp*(ts_quad_fixed_nodes[k+n]+1))/2;
        E=sqrt(u*u+b*b);
        expterm=exp(E-c);
        double fdterm=1/(expterm-1);
        int_val_n+=ts_quad_fixed_weights[k+n]*u*u*fdterm;
        int_val_nT+=ts_quad_fixed_weights[k+n]*u*u*fdterm*fdterm*expterm*E;
        int_val_nc+=ts_quad_fixed_weights[k+n]*u*u*fdterm*fdterm*expterm;
        
        
        int_val_p+=ts_quad_fixed_weights[k+n]*u*u*fdterm*E;
        int_val_pT+=ts_quad_fixed_weights[k+n]*u*u*fdterm*fdterm*expterm*E*E;
        int_val_pc+=ts_quad_fixed_weights[k+n]*u*u*fdterm*fdterm*expterm*E;
        other_fphi_deriv+=ts_quad_fixed_weights[k+n]*(3*u*u*E+u*u*u*u/E)*fdterm*fdterm*expterm;
        
    }
    struct n_rho_explicit2 nrv;
    double const_term=1/(2*M_PI*M_PI);
    nrv.n=T*T*T*int_val_n*integration_endp/2*const_term;
    nrv.nT=T*T*int_val_nT*integration_endp/2*const_term;
    nrv.nc=T*T*T*int_val_nc*integration_endp/2*const_term;
    
    nrv.p=T*T*T*T*int_val_p*integration_endp/2*const_term;
    nrv.pT=T*T*T*int_val_pT*integration_endp/2*const_term;
    nrv.pc=T*T*T*T*int_val_pc*integration_endp/2*const_term;
    nrv.other_fphi_deriv=T*T*T*T*other_fphi_deriv*integration_endp/2*const_term;
    return nrv;
}

struct n_p_vals2 get_n_p_vals2(double mphi, double Ts, double cs, double Tp, double cp,struct polylogs* plgs)
{
    struct n_p_vals2 npv;

    double li2val=-interp_Li2(-exp(cs),plgs);
    double li3val=-2*interp_Li3(-exp(cs),plgs);
    double li4val=-6*interp_Li4(-exp(cs),plgs);

    double s_coeff=Ts*Ts*Ts/(2*M_PI*M_PI);
    double p_coeff=Tp*Tp*Tp/(2*M_PI*M_PI);
    
    npv.n_s=s_coeff*li3val;
    npv.n_s_c=s_coeff*2*li2val;
    npv.n_s_T=3*Ts*Ts/(2*M_PI*M_PI)*li3val;
    npv.p_s=Ts*s_coeff*li4val;
    npv.p_s_c=Ts*s_coeff*3*li3val;
    npv.p_s_T=2*Ts*Ts*Ts/(M_PI*M_PI)*li4val;


    struct n_rho_explicit2 nre;
    nre=ts_quad_n_rho_explicit2(mphi,Tp,cp);

    npv.In_p=nre.n;
    npv.In_p_c=nre.nc;
    npv.In_p_T=nre.nT;
    npv.Ip_p=nre.p;
    npv.Ip_p_c=nre.pc;
    npv.Ip_p_T=nre.pT;
    npv.other_fphi_deriv=nre.other_fphi_deriv;

    return npv;
}