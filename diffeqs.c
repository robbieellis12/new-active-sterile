#include "diffeqs.h"

struct final_derivs step(struct sp_params sp, struct dx_params* dx_params)
{
double mphi=dx_params->fixed_params->m_phi;
double T=sp.T;
struct derivs new_derivs=integrated_de(0,100,*dx_params->fixed_params,sp,dx_params->plgs,dx_params->mip);
double current_H=H(T,dx_params->fixed_params);
double Tp=sp.Tphi;
double cp=sp.c_phi;
double Ts=sp.Ts;
double cs=sp.c_s;
struct polylogs* plgs=dx_params->plgs;


struct n_p_vals npv=get_n_p_vals(mphi,Ts, cs, Tp,  cp, plgs);

double In_p=npv.In_p;
double In_p_b=npv.In_p_b;
double In_p_c=npv.In_p_c;
double n_p=Tp*Tp*Tp/(2*M_PI*M_PI)*In_p;

double Ip_p=npv.Ip_p;
double Ip_p_b=npv.Ip_p_b;
double Ip_p_c=npv.Ip_p_c;
double p_p=Tp*Tp*Tp*Tp/(2*M_PI*M_PI)*Ip_p;


//double ap=(3/Tp*In_p-mphi/(Tp*Tp)*In_p_b);
//double bp=In_p_c;
//double Cp=(4/Tp*Ip_p-mphi/(Tp*Tp)*Ip_p_b);
//double dp=Ip_p_c;

double ap=(3*Tp*Tp*In_p-mphi*Tp*In_p_b);
double bp=In_p_c*Tp*Tp*Tp;
double Cp=(4*Tp*Tp*Tp*Ip_p-mphi*Tp*Tp*Ip_p_b);
double dp=Ip_p_c*Tp*Tp*Tp*Tp;
double detp=ap*dp-bp*cp;




double n_p_T=new_derivs.C_n_p-3*current_H*n_p;
double p_p_T=new_derivs.C_p_p-4*current_H*p_p;
double Ap=n_p_T*2*M_PI*M_PI;
double Bp=p_p_T*2*M_PI*M_PI;


double as=npv.n_s_T;
double bs=npv.n_s_c;
double Cs=npv.p_s_T;
double ds=npv.p_s_c;
double dets=as*ds-bs*cs;

double n_s=npv.n_s;
double p_s=npv.p_s;

double n_s_T=new_derivs.C_n_s-3*current_H*n_s;
double p_s_T=new_derivs.C_p_s-4*current_H*p_s;

double dTs=1/dets*(ds*n_s_T-bs*p_s_T);
double dcs=1/dets*(-Cs*n_s_T+as*p_s_T);

double A_p=n_p_T*2*M_PI*M_PI;
double B_p=p_p_T*2*M_PI*M_PI;
double dTp=1/detp*(dp*A_p-bp*B_p);
double dcp=1/detp*(-Cp*A_p+ap*B_p);

struct final_derivs final_derivs;

final_derivs.dTs=-dTs/(T*current_H);
final_derivs.dTp=-dTp/(T*current_H);
final_derivs.dcs=-dcs/(T*current_H);
final_derivs.dcp=-dcp/(T*current_H);
return final_derivs;

}




struct final_derivs step_omp(struct sp_params sp, struct dx_params* dx_params)
{
double mphi=dx_params->fixed_params->m_phi;
double T=sp.T;
struct derivs new_derivs=integrated_de_omp(0,100,*dx_params->fixed_params,sp,dx_params->plgs,dx_params->mip);
double current_H=H(T,dx_params->fixed_params);
double Tp=sp.Tphi;
double cp=sp.c_phi;
double Ts=sp.Ts;
double cs=sp.c_s;
struct polylogs* plgs=dx_params->plgs;


struct n_p_vals npv=get_n_p_vals(mphi,Ts, cs, Tp,  cp, plgs);

double In_p=npv.In_p;
double In_p_b=npv.In_p_b;
double In_p_c=npv.In_p_c;
double n_p=Tp*Tp*Tp/(2*M_PI*M_PI)*In_p;

double Ip_p=npv.Ip_p;
double Ip_p_b=npv.Ip_p_b;
double Ip_p_c=npv.Ip_p_c;
double p_p=Tp*Tp*Tp*Tp/(2*M_PI*M_PI)*Ip_p;


//double ap=(3/Tp*In_p-mphi/(Tp*Tp)*In_p_b);
//double bp=In_p_c;
//double Cp=(4/Tp*Ip_p-mphi/(Tp*Tp)*Ip_p_b);
//double dp=Ip_p_c;

double ap=(3*Tp*Tp*In_p-mphi*Tp*In_p_b);
double bp=In_p_c*Tp*Tp*Tp;
double Cp=(4*Tp*Tp*Tp*Ip_p-mphi*Tp*Tp*Ip_p_b);
double dp=Ip_p_c*Tp*Tp*Tp*Tp;
double detp=ap*dp-bp*cp;




double n_p_T=new_derivs.C_n_p-3*current_H*n_p;
double p_p_T=new_derivs.C_p_p-4*current_H*p_p;
double Ap=n_p_T*2*M_PI*M_PI;
double Bp=p_p_T*2*M_PI*M_PI;


double as=npv.n_s_T;
double bs=npv.n_s_c;
double Cs=npv.p_s_T;
double ds=npv.p_s_c;
double dets=as*ds-bs*cs;

double n_s=npv.n_s;
double p_s=npv.p_s;

double n_s_T=new_derivs.C_n_s-3*current_H*n_s;
double p_s_T=new_derivs.C_p_s-4*current_H*p_s;

double dTs=1/dets*(ds*n_s_T-bs*p_s_T);
double dcs=1/dets*(-Cs*n_s_T+as*p_s_T);

double A_p=n_p_T*2*M_PI*M_PI;
double B_p=p_p_T*2*M_PI*M_PI;
double dTp=1/detp*(dp*A_p-bp*B_p);
double dcp=1/detp*(-Cp*A_p+ap*B_p);

struct final_derivs final_derivs;

final_derivs.dTs=-dTs/(T*current_H);
final_derivs.dTp=-dTp/(T*current_H);
final_derivs.dcs=-dcs/(T*current_H);
final_derivs.dcp=-dcp/(T*current_H);
return final_derivs;

}



struct final_derivs step_omp2(struct sp_params sp, struct dx_params* dx_params)
{
double mphi=dx_params->fixed_params->m_phi;
double T=sp.T;
struct derivs new_derivs=integrated_de_omp(0,100,*dx_params->fixed_params,sp,dx_params->plgs,dx_params->mip);
double current_H=H(T,dx_params->fixed_params);
double Tp=sp.Tphi;
double cp=sp.c_phi;
double Ts=sp.Ts;
double cs=sp.c_s;
struct polylogs* plgs=dx_params->plgs;


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




double n_p_T=new_derivs.C_n_p-3*current_H*n_p;
double p_p_T=new_derivs.C_p_p-4*current_H*p_p;
double Ap=n_p_T;
double Bp=p_p_T;


double as=npv.n_s_T;
double bs=npv.n_s_c;
double Cs=npv.p_s_T;
double ds=npv.p_s_c;
double dets=as*ds-bs*Cs;

double n_s=npv.n_s;
double p_s=npv.p_s;

double n_s_T=new_derivs.C_n_s-3*current_H*n_s;
double p_s_T=new_derivs.C_p_s-current_H*npv.other_fphi_deriv;

double dTs=1/dets*(ds*n_s_T-bs*p_s_T);
double dcs=1/dets*(-Cs*n_s_T+as*p_s_T);

double A_p=n_p_T;
double B_p=p_p_T;
double dTp=1/detp*(dp*A_p-bp*B_p);
double dcp=1/detp*(-Cp*A_p+ap*B_p);


//printf("n_s: %.10e\n",n_s);
//printf("p_s: %.10e\n",p_s);
//printf("n_p: %.10e\n",In_p);
//printf("p_p: %.10e\n",Ip_p);

struct final_derivs final_derivs;

final_derivs.dTs=-dTs/(T*current_H);
final_derivs.dTp=-dTp/(T*current_H);
final_derivs.dcs=-dcs/(T*current_H);
final_derivs.dcp=-dcp/(T*current_H);
return final_derivs;

}


struct final_derivs_eq step_omp_eq(struct sp_params sp, struct dx_eq_params* dx_params)
{
double mphi=dx_params->fixed_params->m_phi;
double T=dx_params->T;
struct derivs new_derivs=integrated_de_omp_eq(0,100,*dx_params->fixed_params,sp,dx_params->plgs);

double Tp=sp.Tphi;
double cp=sp.c_phi;
double Ts=sp.Ts;
double cs=sp.c_s;
struct polylogs* plgs=dx_params->plgs;


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

struct final_derivs_eq final_derivs;

final_derivs.dTs=dTs;
final_derivs.dTp=dTp;
final_derivs.dcs=dcs;
final_derivs.dcp=dcp;
final_derivs.n_s=n_s;
final_derivs.p_s=p_s;
final_derivs.n_p=n_p;
final_derivs.p_p=p_p;

return final_derivs;

}


struct solved_params rk_solve(double T0, struct dx_params* dx_params,double h, double N)
{
    double* T=(double*) malloc((N+1)*sizeof(double));
    
    double* Ts=(double*) malloc((N+1)*sizeof(double));
    double* c_s=(double*) malloc((N+1)*sizeof(double));
    double* Tphi=(double*) malloc((N+1)*sizeof(double));
    double* c_phi=(double*) malloc((N+1)*sizeof(double));
    
    //ret_val.x=x;
    //ret_val.y=y;
    
    //x[0]=x0;
struct sp_params sp_params;
    sp_params.T=T0;
    sp_params.Ts=49.188372;
    sp_params.Tphi=36.679233;
    sp_params.c_s=0.089456; //0
    sp_params.c_phi=-2.201947; //-2
    T[0]=sp_params.T;
    Ts[0]=sp_params.Ts;
    Tphi[0]=sp_params.Tphi;
    c_s[0]=sp_params.c_s;
    c_phi[0]=sp_params.c_phi;

    double df;
    
    double a[21]=
    {1/5.0,
    3/40.0, 9/40.0,
    44/45.0, -56/15.0, 32/9.0,
    19372/6561.0,    -25360/2187.0,    64448/6561.0,    -212/729.0,
    9017/3168.0,    -355/33.0,    46732/5247.0,    49/176.0,    -5103/18656.0,
    35/384.0,    0.0,    500/1113.0,    125/192.0,    -2187/6784.0,    11/84.0,
    };
    double cs[6]={1/5.0, 3/10.0, 4/5.0, 8/9.0, 1, 1};
    
    double bs[7]={35/384.0,    0.0,    500/1113.0,    125/192.0,    -2187/6784.0,    11/84.0,    0.0};
    //double ki[7];
    struct final_derivs sp_ki[7];
    
    //double yki=0;
double Ts_ki=0;
double Tphi_ki=0;
double cs_ki=0;
double cphi_ki=0;

    int aidx=0;
    //double ytemp=0;
    double Ts_temp=0;
    double Tphi_temp=0;
    double cs_temp=0;
    double cphi_temp=0;
    
    
    for(int i=1;i<N+1;++i)
    {
       //ki[0]=f(x[i-1],y[i-1]);
       sp_ki[0]=step(sp_params,dx_params);
       for(int j=1;j<7;++j)
       {
           // yki=y[i-1];
           Ts_ki = Ts[i - 1];
           Tphi_ki = Tphi[i - 1];
           cs_ki = c_s[i - 1];
           cphi_ki = c_phi[i - 1];
           aidx=j*(j-1)/2;
           for(int k=0;k<j;++k)
           {
               //yki+=a[aidx+k]*ki[k]*h;
                Ts_ki+=a[aidx+k]*sp_ki[k].dTs*h;
                Tphi_ki+=a[aidx+k]*sp_ki[k].dTp*h;
                cs_ki+=a[aidx+k]*sp_ki[k].dcs*h;
                cphi_ki+=a[aidx+k]*sp_ki[k].dcp*h;
           }
           // ki[j]=f(x[i-1]+cs[j]*h,yki);
           sp_params.Ts = Ts_ki;
           sp_params.Tphi = Tphi_ki;
           sp_params.c_s = cs_ki;
           sp_params.c_phi = cphi_ki;
           sp_params.T = T[i - 1] + cs[j-1] * h;
           sp_ki[j]=step(sp_params,dx_params);
       }
        //ytemp=y[i-1];
        Ts_temp=Ts[i-1];
        Tphi_temp=Tphi[i-1];
        cs_temp=c_s[i-1];
        cphi_temp=c_phi[i-1];
        for(int j=0;j<7;++j)
        {
            //ytemp=ytemp+bs[j]*ki[j]*h;
            Ts_temp+=bs[j]*sp_ki[j].dTs*h;
            Tphi_temp+=bs[j]*sp_ki[j].dTp*h;
            cs_temp+=bs[j]*sp_ki[j].dcs*h;
            cphi_temp+=bs[j]*sp_ki[j].dcp*h;
        }
       // y[i]=ytemp;
        Ts[i]=Ts_temp;
        Tphi[i]=Tphi_temp;
        c_s[i]=cs_temp;
        c_phi[i]=cphi_temp;
        //x[i]=x[i-1]+h;
        T[i]=T[i-1]+h;
        printf("Now at T: %f\n",T[i]);
        printf("Ts: %f\n",Ts[i]);
        printf("Tphi: %f\n",Tphi[i]);
        printf("cs: %f\n",c_s[i]);
        printf("cphi: %f\n",c_phi[i]);
        sp_params.T=T[i];
        sp_params.Ts=Ts[i];
        sp_params.Tphi=Tphi[i];
        sp_params.c_s=c_s[i];
        sp_params.c_phi=c_phi[i];
    }
    struct solved_params solved_params;
    solved_params.T=T;
    solved_params.Ts=Ts;
    solved_params.Tphi=Tphi;
    solved_params.c_s=c_s;
    solved_params.c_phi=c_phi;
    return solved_params;
}




struct solved_params rk_solve_omp(double T0, struct dx_params* dx_params,double h, double N)
{
    double* T=(double*) malloc((N+1)*sizeof(double));
    
    double* Ts=(double*) malloc((N+1)*sizeof(double));
    double* c_s=(double*) malloc((N+1)*sizeof(double));
    double* Tphi=(double*) malloc((N+1)*sizeof(double));
    double* c_phi=(double*) malloc((N+1)*sizeof(double));
    #pragma omp parallel master
    {
    
    //ret_val.x=x;
    //ret_val.y=y;
    
    //x[0]=x0;
struct sp_params sp_params;
    sp_params.T=T0;
    sp_params.Ts=49.188372;
    sp_params.Tphi=36.679233;
    sp_params.c_s=0.089456; //0
    sp_params.c_phi=-2.201947; //-2
    T[0]=sp_params.T;
    Ts[0]=sp_params.Ts;
    Tphi[0]=sp_params.Tphi;
    c_s[0]=sp_params.c_s;
    c_phi[0]=sp_params.c_phi;

    double df;
    
    double a[21]=
    {1/5.0,
    3/40.0, 9/40.0,
    44/45.0, -56/15.0, 32/9.0,
    19372/6561.0,    -25360/2187.0,    64448/6561.0,    -212/729.0,
    9017/3168.0,    -355/33.0,    46732/5247.0,    49/176.0,    -5103/18656.0,
    35/384.0,    0.0,    500/1113.0,    125/192.0,    -2187/6784.0,    11/84.0,
    };
    double cs[6]={1/5.0, 3/10.0, 4/5.0, 8/9.0, 1, 1};
    
    double bs[7]={35/384.0,    0.0,    500/1113.0,    125/192.0,    -2187/6784.0,    11/84.0,    0.0};
    //double ki[7];
    struct final_derivs sp_ki[7];
    
    //double yki=0;
double Ts_ki=0;
double Tphi_ki=0;
double cs_ki=0;
double cphi_ki=0;

    int aidx=0;
    //double ytemp=0;
    double Ts_temp=0;
    double Tphi_temp=0;
    double cs_temp=0;
    double cphi_temp=0;
    
    
    for(int i=1;i<N+1;++i)
    {
       //ki[0]=f(x[i-1],y[i-1]);
       sp_ki[0]=step_omp(sp_params,dx_params);
       for(int j=1;j<7;++j)
       {
           // yki=y[i-1];
           Ts_ki = Ts[i - 1];
           Tphi_ki = Tphi[i - 1];
           cs_ki = c_s[i - 1];
           cphi_ki = c_phi[i - 1];
           aidx=j*(j-1)/2;
           for(int k=0;k<j;++k)
           {
               //yki+=a[aidx+k]*ki[k]*h;
                Ts_ki+=a[aidx+k]*sp_ki[k].dTs*h;
                Tphi_ki+=a[aidx+k]*sp_ki[k].dTp*h;
                cs_ki+=a[aidx+k]*sp_ki[k].dcs*h;
                cphi_ki+=a[aidx+k]*sp_ki[k].dcp*h;
           }
           // ki[j]=f(x[i-1]+cs[j]*h,yki);
           sp_params.Ts = Ts_ki;
           sp_params.Tphi = Tphi_ki;
           sp_params.c_s = cs_ki;
           sp_params.c_phi = cphi_ki;
           sp_params.T = T[i - 1] + cs[j-1] * h;
           sp_ki[j]=step_omp(sp_params,dx_params);
       }
        //ytemp=y[i-1];
        Ts_temp=Ts[i-1];
        Tphi_temp=Tphi[i-1];
        cs_temp=c_s[i-1];
        cphi_temp=c_phi[i-1];
        for(int j=0;j<7;++j)
        {
            //ytemp=ytemp+bs[j]*ki[j]*h;
            Ts_temp+=bs[j]*sp_ki[j].dTs*h;
            Tphi_temp+=bs[j]*sp_ki[j].dTp*h;
            cs_temp+=bs[j]*sp_ki[j].dcs*h;
            cphi_temp+=bs[j]*sp_ki[j].dcp*h;
        }
       // y[i]=ytemp;
        Ts[i]=Ts_temp;
        Tphi[i]=Tphi_temp;
        c_s[i]=cs_temp;
        c_phi[i]=cphi_temp;
        //x[i]=x[i-1]+h;
        T[i]=T[i-1]+h;
        printf("Now at T: %f\n",T[i]);
        printf("Ts: %f\n",Ts[i]);
        printf("Tphi: %f\n",Tphi[i]);
        printf("cs: %f\n",c_s[i]);
        printf("cphi: %f\n",c_phi[i]);
        sp_params.T=T[i];
        sp_params.Ts=Ts[i];
        sp_params.Tphi=Tphi[i];
        sp_params.c_s=c_s[i];
        sp_params.c_phi=c_phi[i];
    }
    }
    struct solved_params solved_params;
    solved_params.T=T;
    solved_params.Ts=Ts;
    solved_params.Tphi=Tphi;
    solved_params.c_s=c_s;
    solved_params.c_phi=c_phi;
    return solved_params;
}