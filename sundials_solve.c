#include "sundials_solve.h"


#define Ith(v, i) NV_Ith_S(v, i - 1) 
#define IJth(A, i, j) \
  SM_ELEMENT_D(A, i - 1, j - 1)

int solve_de(double T0, struct sp_params sp_params_init,struct dx_params* dx_params)
{
    int NEQ=4;
SUNContext sunctx;
  sunrealtype t, tout;
  N_Vector sp_param_vec;
  N_Vector abstol;
  SUNMatrix A;
  SUNLinearSolver LS;
  void* cvode_mem;
  int retval, iout;
  int retvalr;
  int rootsfound[2];
  FILE* FID;
  retval=SUNContext_Create(SUN_COMM_NULL, &sunctx);
  if (check_retval(&retval, "SUNContext_Create", 1)) { return (1); }

  sp_param_vec = N_VNew_Serial(NEQ, sunctx);
  if (check_retval((void*)sp_param_vec, "N_VNew_Serial", 0)) { return (1); }

    Ith(sp_param_vec, 1) = sp_params_init.Ts;
    Ith(sp_param_vec, 2) = sp_params_init.c_s;
    Ith(sp_param_vec, 3) = sp_params_init.Tphi;
    Ith(sp_param_vec, 4) = sp_params_init.c_phi;

abstol = N_VNew_Serial(NEQ, sunctx);
  if (check_retval((void*)abstol, "N_VNew_Serial", 0)) { return (1); }

  Ith(abstol, 1) = 1e-3;
  Ith(abstol, 2) = 1e-5;
  Ith(abstol, 3) = 1e-3;
  Ith(abstol, 4) = 1e-5;


  cvode_mem = CVodeCreate(CV_BDF, sunctx);
  if (check_retval((void*)cvode_mem, "CVodeCreate", 0)) { return (1); }

    retval = CVodeInit(cvode_mem, f, T0, sp_param_vec);
  if (check_retval(&retval, "CVodeInit", 1)) { return (1); }

  retval = CVodeSVtolerances(cvode_mem, 1e-4, abstol);
  if (check_retval(&retval, "CVodeSVtolerances", 1)) { return (1); }

  A = SUNDenseMatrix(NEQ, NEQ, sunctx);
  if (check_retval((void*)A, "SUNDenseMatrix", 0)) { return (1); }

  LS = SUNLinSol_Dense(sp_param_vec, A, sunctx);
  if (check_retval((void*)LS, "SUNLinSol_Dense", 0)) { return (1); }

int userdata_retval=CVodeSetUserData(cvode_mem, dx_params);
if (check_retval(&userdata_retval, "CVodeSetUserData", 1)) { return (1); }


  retval = CVodeSetLinearSolver(cvode_mem, LS, A);
  if (check_retval(&retval, "CVodeSetLinearSolver", 1)) { return (1); }

  retval = CVodeSetJacFn(cvode_mem, Jac);
  if (check_retval(&retval, "CVodeSetJacFn", 1)) { return (1); }

  FID = fopen("/home/robbie-ellis/VSCode/new active_sterile/exported/sundialsdata.txt", "w");
  if(FID==NULL)
  {
      printf("Error opening file\n");
      return 1;
  }
  iout = 0;
  tout = T0-0.1;
  while (1)
  {
    retval = CVode(cvode_mem, tout, sp_param_vec, &t, CV_NORMAL);
    printf("Here\n");
    PrintOutput(t, Ith(sp_param_vec, 1), Ith(sp_param_vec, 2), Ith(sp_param_vec, 3), Ith(sp_param_vec, 4));

  //fprintf(FID,"%.10e, %.10e, %.10e, %.10e, %.10e\n", t, Ith(sp_param_vec, 1), Ith(sp_param_vec, 2), Ith(sp_param_vec, 3));
    if (check_retval(&retval, "CVode", 1)) { break; }
    if (retval == CV_SUCCESS)
    {
      iout++;
      tout -= 0.1;
    }

    //retval = CVodePrintAllStats(cvode_mem, FID, SUN_OUTPUTFORMAT_CSV);

    if (iout == 300) { break; }
  }
  fclose(FID);
}




static int f(sunrealtype T, N_Vector sp_params_vec, N_Vector d_sp_params, void* user_data)
{
struct dx_params* params=(struct dx_params*) user_data;
struct sp_params sp_params;
printf("T: %f\n",T);
printf("Tphi: %f\n",Ith(sp_params_vec,3));
printf("cphi: %f\n",Ith(sp_params_vec,4));
sp_params.T=T;
sp_params.Ts=Ith(sp_params_vec,1);
sp_params.c_s=Ith(sp_params_vec,2);
sp_params.Tphi=Ith(sp_params_vec,3);
sp_params.c_phi=Ith(sp_params_vec,4);
struct final_derivs derivs=step_omp2(sp_params,params);
Ith(d_sp_params,1)=derivs.dTs;
Ith(d_sp_params,2)=derivs.dcs;
Ith(d_sp_params,3)=derivs.dTp;
Ith(d_sp_params,4)=derivs.dcp;
//printf("dT: %f\n",derivs.dTs);
//printf("dTphi: %f\n",derivs.dTp);
return 0;
}




static int Jac(sunrealtype T, N_Vector sp_param_vec, N_Vector fy, SUNMatrix J,
               void* user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
{
struct dx_params* params=(struct dx_params*) user_data;
struct sp_params sp_params;
sp_params.T=T;
sp_params.Ts=Ith(sp_param_vec,1);
sp_params.c_s=Ith(sp_param_vec,2);
sp_params.Tphi=Ith(sp_param_vec,3);
sp_params.c_phi=Ith(sp_param_vec,4);
struct final_derivs derivs_x0=step_omp2(sp_params,params);

double h=1e-3;

struct sp_params sp_params_Ts;
sp_params_Ts.T=T;
sp_params_Ts.Ts=Ith(sp_param_vec,1)+h;
sp_params_Ts.c_s=Ith(sp_param_vec,2);
sp_params_Ts.Tphi=Ith(sp_param_vec,3);
sp_params_Ts.c_phi=Ith(sp_param_vec,4);
struct final_derivs derivs_Ts=step_omp2(sp_params_Ts,params);

struct sp_params sp_params_cs;
sp_params_cs.T=T;
sp_params_cs.Ts=Ith(sp_param_vec,1);
sp_params_cs.c_s=Ith(sp_param_vec,2)+h;
sp_params_cs.Tphi=Ith(sp_param_vec,3);
sp_params_cs.c_phi=Ith(sp_param_vec,4);
struct final_derivs derivs_cs=step_omp2(sp_params_cs,params);

struct sp_params sp_params_Tphi;
sp_params_Tphi.T=T;
sp_params_Tphi.Ts=Ith(sp_param_vec,1);
sp_params_Tphi.c_s=Ith(sp_param_vec,2);
sp_params_Tphi.Tphi=Ith(sp_param_vec,3)+h;
sp_params_Tphi.c_phi=Ith(sp_param_vec,4);
struct final_derivs derivs_Tphi=step_omp2(sp_params_Tphi,params);

struct sp_params sp_params_cphi;
sp_params_cphi.T=T;
sp_params_cphi.Ts=Ith(sp_param_vec,1);
sp_params_cphi.c_s=Ith(sp_param_vec,2);
sp_params_cphi.Tphi=Ith(sp_param_vec,3);
sp_params_cphi.c_phi=Ith(sp_param_vec,4)+h;
struct final_derivs derivs_cphi=step_omp2(sp_params_cphi,params);



  IJth(J, 1, 1) = (derivs_Tphi.dTs-derivs_x0.dTs)/h;
  IJth(J, 1, 2) = (derivs_Tphi.dcs-derivs_x0.dcs)/h;
  IJth(J, 1, 3) = (derivs_Tphi.dTp-derivs_x0.dTp)/h;
  IJth(J, 1, 4) = (derivs_Tphi.dcp-derivs_x0.dcp)/h;

  IJth(J, 2, 1) = (derivs_cs.dTs-derivs_x0.dTs)/h;
    IJth(J, 2, 2) = (derivs_cs.dcs-derivs_x0.dcs)/h;
    IJth(J, 2, 3) = (derivs_cs.dTp-derivs_x0.dTp)/h;
    IJth(J, 2, 4) = (derivs_cs.dcp-derivs_x0.dcp)/h;

    IJth(J, 3, 1) = (derivs_Tphi.dTs-derivs_x0.dTs)/h;
    IJth(J, 3, 2) = (derivs_Tphi.dcs-derivs_x0.dcs)/h;
    IJth(J, 3, 3) = (derivs_Tphi.dTp-derivs_x0.dTp)/h;
    IJth(J, 3, 4) = (derivs_Tphi.dcp-derivs_x0.dcp)/h;

    IJth(J, 4, 1) = (derivs_cphi.dTs-derivs_x0.dTs)/h;
    IJth(J, 4, 2) = (derivs_cphi.dcs-derivs_x0.dcs)/h;
    IJth(J, 4, 3) = (derivs_cphi.dTp-derivs_x0.dTp)/h;
    IJth(J, 4, 4) = (derivs_cphi.dcp-derivs_x0.dcp)/h;

  return (0);
}

 void PrintOutput(sunrealtype t, sunrealtype y1, sunrealtype y2,
                        sunrealtype y3, sunrealtype y4)
{
  //printf("At t = %0.4e      y =%14.6e  %14.6e  %14.6e  %14.6e\n", t, y1, y2, y3, y4);
  
#if defined(SUNDIALS_EXTENDED_PRECISION)
  printf("At t = %0.4Le      y =%14.6Le  %14.6Le  %14.6Le\n", t, y1, y2, y3);
#elif defined(SUNDIALS_DOUBLE_PRECISION)
  printf("At t = %0.4e      y =%14.6e  %14.6e  %14.6e  %14.6e\n", t, y1, y2, y3, y4);
#else
  printf("At t = %0.4e      y =%14.6e  %14.6e  %14.6e\n", t, y1, y2, y3);
#endif

  return;
}


 int check_retval(void* returnvalue, const char* funcname, int opt)
{
  int* retval;

  /* Check if SUNDIALS function returned NULL pointer - no memory allocated */
  if (opt == 0 && returnvalue == NULL)
  {
    fprintf(stderr, "\nSUNDIALS_ERROR: %s() failed - returned NULL pointer\n\n",
            funcname);
    return (1);
  }

  /* Check if retval < 0 */
  else if (opt == 1)
  {
    retval = (int*)returnvalue;
    if (*retval < 0)
    {
      fprintf(stderr, "\nSUNDIALS_ERROR: %s() failed with retval = %d\n\n",
              funcname, *retval);
      return (1);
    }
  }

  /* Check if function returned NULL pointer - no memory allocated */
  else if (opt == 2 && returnvalue == NULL)
  {
    fprintf(stderr, "\nMEMORY_ERROR: %s() failed - returned NULL pointer\n\n",
            funcname);
    return (1);
  }

  return (0);
}



























struct sp_params solve_de_eq(double T0, struct sp_params sp_params_init,struct dx_eq_params* dx_params)
{

  struct sp_params first_eq;
  first_eq.Ts=-1;
  first_eq.c_s=-1;
  first_eq.Tphi=-1;
  first_eq.c_phi=-1;

    int NEQ=4;
SUNContext sunctx;
  sunrealtype t, tout;
  N_Vector sp_param_vec;
  N_Vector abstol;
  SUNMatrix A;
  SUNLinearSolver LS;
  void* cvode_mem;
  int retval, iout;
  int retvalr;
  int rootsfound[2];
  FILE* FID;
  retval=SUNContext_Create(SUN_COMM_NULL, &sunctx);
  if (check_retval(&retval, "SUNContext_Create", 1)) { return (first_eq); }

  sp_param_vec = N_VNew_Serial(NEQ, sunctx);
  if (check_retval((void*)sp_param_vec, "N_VNew_Serial", 0)) { return (first_eq); }

    Ith(sp_param_vec, 1) = sp_params_init.Ts;
    Ith(sp_param_vec, 2) = sp_params_init.c_s;
    Ith(sp_param_vec, 3) = sp_params_init.Tphi;
    Ith(sp_param_vec, 4) = sp_params_init.c_phi;

abstol = N_VNew_Serial(NEQ, sunctx);
  if (check_retval((void*)abstol, "N_VNew_Serial", 0)) { return (first_eq); }

  Ith(abstol, 1) = 1e-3;
  Ith(abstol, 2) = 1e-8;
  Ith(abstol, 3) = 1e-3;
  Ith(abstol, 4) = 1e-8;


  cvode_mem = CVodeCreate(CV_ADAMS, sunctx);
  if (check_retval((void*)cvode_mem, "CVodeCreate", 0)) { return (first_eq); }

 CVodeSetMaxNumSteps(cvode_mem, 2000);

    retval = CVodeInit(cvode_mem, f_eq, 1, sp_param_vec);
  if (check_retval(&retval, "CVodeInit", 1)) { return (first_eq); }

  retval = CVodeSVtolerances(cvode_mem, 1e-9, abstol);
  if (check_retval(&retval, "CVodeSVtolerances", 1)) { return (first_eq); }

  A = SUNDenseMatrix(NEQ, NEQ, sunctx);
  if (check_retval((void*)A, "SUNDenseMatrix", 0)) { return (first_eq); }

  LS = SUNLinSol_Dense(sp_param_vec, A, sunctx);
  if (check_retval((void*)LS, "SUNLinSol_Dense", 0)) { return (first_eq); }

int userdata_retval=CVodeSetUserData(cvode_mem, dx_params);
if (check_retval(&userdata_retval, "CVodeSetUserData", 1)) { return (first_eq); }


  retval = CVodeSetLinearSolver(cvode_mem, LS, A);
  if (check_retval(&retval, "CVodeSetLinearSolver", 1)) { return (first_eq); }

  retval = CVodeSetJacFn(cvode_mem, Jac_eq);
  if (check_retval(&retval, "CVodeSetJacFn", 1)) { return (first_eq); }

  FID = fopen("/home/robbie-ellis/VSCode/new active_sterile/exported/sundialsdata.txt", "w");
  if(FID==NULL)
  {
      printf("Error opening file\n");
      return first_eq;
  }


  double T_prev=0;
  double T_cur=sp_params_init.Ts;
  iout = 0;
  tout = 1.001;
  double tfin=10000;
  /*
  double T_phi_prev=0;
  double T_phi_cur=sp_params_init.Tphi;
  double c_phi_prev=0;
  double c_phi_cur=sp_params_init.c_phi;
  double c_s_prev=0;
  double c_s_cur=sp_params_init.c_s;*/



  while (tout<tfin)
  {
    //t_ask=tfin;
    retval = CVode(cvode_mem, tout, sp_param_vec, &t, CV_NORMAL);
    printf("iout: %d\n",iout);
    //PrintOutput(t, Ith(sp_param_vec, 1), Ith(sp_param_vec, 2), Ith(sp_param_vec, 3), Ith(sp_param_vec, 4));
    
  fprintf(FID,"%.10e, %.10e, %.10e, %.10e, %.10e\n", t, Ith(sp_param_vec, 1), Ith(sp_param_vec, 2), Ith(sp_param_vec, 3), Ith(sp_param_vec, 4));
  //printf("t-tout: %.10e\n",t-tout);
    if (check_retval(&retval, "CVode", 1)) { break; }
    if (retval == CV_SUCCESS)
    {
      iout++;
      tout+=0.2;
    }

    retval = CVodePrintAllStats(cvode_mem, FID, SUN_OUTPUTFORMAT_CSV);

    if (iout == 2000) { break; }
    /*if(change<max_change/100)
    {
      break;
    }*/
  }
  fclose(FID);
  first_eq.Ts=Ith(sp_param_vec, 1);
  first_eq.c_s=Ith(sp_param_vec, 2);
  first_eq.Tphi=Ith(sp_param_vec, 3);
  first_eq.c_phi=Ith(sp_param_vec, 4);
  return first_eq;

}




static int f_eq(sunrealtype t, N_Vector sp_params_vec, N_Vector d_sp_params, void* user_data)
{
struct dx_eq_params* params=(struct dx_eq_params*) user_data;
struct sp_params sp_params;
//printf("T: %f\n",params->T);
//printf("Tphi: %f\n",Ith(sp_params_vec,3));
//printf("cphi: %f\n",Ith(sp_params_vec,4));
sp_params.Ts=Ith(sp_params_vec,1);
sp_params.c_s=Ith(sp_params_vec,2);
sp_params.Tphi=Ith(sp_params_vec,3);
sp_params.c_phi=Ith(sp_params_vec,4);
sp_params.T=params->T;
struct final_derivs_eq derivs=step_omp_eq(sp_params,params);
Ith(d_sp_params,1)=derivs.dTs;
Ith(d_sp_params,2)=derivs.dcs;
Ith(d_sp_params,3)=derivs.dTp;
Ith(d_sp_params,4)=derivs.dcp;
//printf("dT: %f\n",derivs.dTs);
//printf("dTphi: %f\n",derivs.dTp);
return 0;
}




static int Jac_eq(sunrealtype T, N_Vector sp_param_vec, N_Vector fy, SUNMatrix J,
               void* user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
{
struct dx_eq_params* params=(struct dx_eq_params*) user_data;
struct sp_params sp_params;
sp_params.T=T;
sp_params.Ts=Ith(sp_param_vec,1);
sp_params.c_s=Ith(sp_param_vec,2);
sp_params.Tphi=Ith(sp_param_vec,3);
sp_params.c_phi=Ith(sp_param_vec,4);
struct final_derivs_eq derivs_x0=step_omp_eq(sp_params,params);

double h=1e-4;

struct sp_params sp_params_Ts;
sp_params_Ts.T=T;
sp_params_Ts.Ts=Ith(sp_param_vec,1)+h;
sp_params_Ts.c_s=Ith(sp_param_vec,2);
sp_params_Ts.Tphi=Ith(sp_param_vec,3);
sp_params_Ts.c_phi=Ith(sp_param_vec,4);
struct final_derivs_eq derivs_Ts=step_omp_eq(sp_params_Ts,params);

struct sp_params sp_params_cs;
sp_params_cs.T=T;
sp_params_cs.Ts=Ith(sp_param_vec,1);
sp_params_cs.c_s=Ith(sp_param_vec,2)+h;
sp_params_cs.Tphi=Ith(sp_param_vec,3);
sp_params_cs.c_phi=Ith(sp_param_vec,4);
struct final_derivs_eq derivs_cs=step_omp_eq(sp_params_cs,params);

struct sp_params sp_params_Tphi;
sp_params_Tphi.T=T;
sp_params_Tphi.Ts=Ith(sp_param_vec,1);
sp_params_Tphi.c_s=Ith(sp_param_vec,2);
sp_params_Tphi.Tphi=Ith(sp_param_vec,3)+h;
sp_params_Tphi.c_phi=Ith(sp_param_vec,4);
struct final_derivs_eq derivs_Tphi=step_omp_eq(sp_params_Tphi,params);

struct sp_params sp_params_cphi;
sp_params_cphi.T=T;
sp_params_cphi.Ts=Ith(sp_param_vec,1);
sp_params_cphi.c_s=Ith(sp_param_vec,2);
sp_params_cphi.Tphi=Ith(sp_param_vec,3);
sp_params_cphi.c_phi=Ith(sp_param_vec,4)+h;
struct final_derivs_eq derivs_cphi=step_omp_eq(sp_params_cphi,params);

struct sp_params sp_params_Tsm;
sp_params_Ts.T=T;
sp_params_Ts.Ts=Ith(sp_param_vec,1)-h;
sp_params_Ts.c_s=Ith(sp_param_vec,2);
sp_params_Ts.Tphi=Ith(sp_param_vec,3);
sp_params_Ts.c_phi=Ith(sp_param_vec,4);
struct final_derivs_eq derivs_Tsm=step_omp_eq(sp_params_Ts,params);

struct sp_params sp_params_csm;
sp_params_cs.T=T;
sp_params_cs.Ts=Ith(sp_param_vec,1);
sp_params_cs.c_s=Ith(sp_param_vec,2)-h;
sp_params_cs.Tphi=Ith(sp_param_vec,3);
sp_params_cs.c_phi=Ith(sp_param_vec,4);
struct final_derivs_eq derivs_csm=step_omp_eq(sp_params_cs,params);

struct sp_params sp_params_Tphim;
sp_params_Tphi.T=T;
sp_params_Tphi.Ts=Ith(sp_param_vec,1);
sp_params_Tphi.c_s=Ith(sp_param_vec,2);
sp_params_Tphi.Tphi=Ith(sp_param_vec,3)-h;
sp_params_Tphi.c_phi=Ith(sp_param_vec,4);
struct final_derivs_eq derivs_Tphim=step_omp_eq(sp_params_Tphi,params);

struct sp_params sp_params_cphim;
sp_params_cphi.T=T;
sp_params_cphi.Ts=Ith(sp_param_vec,1);
sp_params_cphi.c_s=Ith(sp_param_vec,2);
sp_params_cphi.Tphi=Ith(sp_param_vec,3);
sp_params_cphi.c_phi=Ith(sp_param_vec,4)-h;
struct final_derivs_eq derivs_cphim=step_omp_eq(sp_params_cphi,params);
  IJth(J, 1, 1) = (derivs_Tphi.dTs-derivs_Tphim.dTs)/(2*h);
  IJth(J, 1, 2) = (derivs_Tphi.dcs-derivs_Tphim.dcs)/(2*h);
  IJth(J, 1, 3) = (derivs_Tphi.dTp-derivs_Tphim.dTp)/(2*h);
  IJth(J, 1, 4) = (derivs_Tphi.dcp-derivs_Tphim.dcp)/(2*h);

  IJth(J, 2, 1) = (derivs_cs.dTs-derivs_csm.dTs)/(2*h);
    IJth(J, 2, 2) = (derivs_cs.dcs-derivs_csm.dcs)/(2*h);
    IJth(J, 2, 3) = (derivs_cs.dTp-derivs_csm.dTp)/(2*h);
    IJth(J, 2, 4) = (derivs_cs.dcp-derivs_csm.dcp)/(2*h);

    IJth(J, 3, 1) = (derivs_Tphi.dTs-derivs_Tphim.dTs)/(2*h);
    IJth(J, 3, 2) = (derivs_Tphi.dcs-derivs_Tphim.dcs)/(2*h);
    IJth(J, 3, 3) = (derivs_Tphi.dTp-derivs_Tphim.dTp)/(2*h);
    IJth(J, 3, 4) = (derivs_Tphi.dcp-derivs_Tphim.dcp)/(2*h);

    IJth(J, 4, 1) = (derivs_cphi.dTs-derivs_cphim.dTs)/(2*h);
    IJth(J, 4, 2) = (derivs_cphi.dcs-derivs_cphim.dcs)/(2*h);
    IJth(J, 4, 3) = (derivs_cphi.dTp-derivs_cphim.dTp)/(2*h);
    IJth(J, 4, 4) = (derivs_cphi.dcp-derivs_cphim.dcp)/(2*h);

/*
  IJth(J, 1, 1) = (derivs_Tphi.dTs-derivs_x0.dTs)/h;
  IJth(J, 1, 2) = (derivs_Tphi.dcs-derivs_x0.dcs)/h;
  IJth(J, 1, 3) = (derivs_Tphi.dTp-derivs_x0.dTp)/h;
  IJth(J, 1, 4) = (derivs_Tphi.dcp-derivs_x0.dcp)/h;

  IJth(J, 2, 1) = (derivs_cs.dTs-derivs_x0.dTs)/h;
    IJth(J, 2, 2) = (derivs_cs.dcs-derivs_x0.dcs)/h;
    IJth(J, 2, 3) = (derivs_cs.dTp-derivs_x0.dTp)/h;
    IJth(J, 2, 4) = (derivs_cs.dcp-derivs_x0.dcp)/h;

    IJth(J, 3, 1) = (derivs_Tphi.dTs-derivs_x0.dTs)/h;
    IJth(J, 3, 2) = (derivs_Tphi.dcs-derivs_x0.dcs)/h;
    IJth(J, 3, 3) = (derivs_Tphi.dTp-derivs_x0.dTp)/h;
    IJth(J, 3, 4) = (derivs_Tphi.dcp-derivs_x0.dcp)/h;

    IJth(J, 4, 1) = (derivs_cphi.dTs-derivs_x0.dTs)/h;
    IJth(J, 4, 2) = (derivs_cphi.dcs-derivs_x0.dcs)/h;
    IJth(J, 4, 3) = (derivs_cphi.dTp-derivs_x0.dTp)/h;
    IJth(J, 4, 4) = (derivs_cphi.dcp-derivs_x0.dcp)/h;
*/
  return (0);
}












int solve_de_coll(double T0, struct sp_params sp_params_init,struct fixed_params* fixed_params,struct interpolation_tables* all_tables)
{

  struct coll_parameters coll_params;
  coll_params.all_tables=all_tables;
  coll_params.fixed_params=fixed_params;
  coll_params.T=T0;
    int NEQ=4;
SUNContext sunctx;
  sunrealtype t, tout;
  N_Vector sp_param_vec;
  N_Vector abstol;
  SUNMatrix A;
  SUNLinearSolver LS;
  void* cvode_mem;
  int retval, iout;
  int retvalr;
  int rootsfound[2];
  FILE* FID;
  retval=SUNContext_Create(SUN_COMM_NULL, &sunctx);
  if (check_retval(&retval, "SUNContext_Create", 1)) { return (1); }

  sp_param_vec = N_VNew_Serial(NEQ, sunctx);
  if (check_retval((void*)sp_param_vec, "N_VNew_Serial", 0)) { return (1); }

    Ith(sp_param_vec, 1) = sp_params_init.Ts;
    Ith(sp_param_vec, 2) = sp_params_init.c_s;
    Ith(sp_param_vec, 3) = sp_params_init.Tphi;
    Ith(sp_param_vec, 4) = sp_params_init.c_phi;

abstol = N_VNew_Serial(NEQ, sunctx);
  if (check_retval((void*)abstol, "N_VNew_Serial", 0)) { return (1); }

  Ith(abstol, 1) = 1e-3;
  Ith(abstol, 2) = 1e-5;
  Ith(abstol, 3) = 1e-3;
  Ith(abstol, 4) = 1e-5;


  cvode_mem = CVodeCreate(CV_BDF, sunctx);
  if (check_retval((void*)cvode_mem, "CVodeCreate", 0)) { return (1); }

 CVodeSetMaxNumSteps(cvode_mem, 2000);
    retval = CVodeInit(cvode_mem, f_coll, T0, sp_param_vec);
  if (check_retval(&retval, "CVodeInit", 1)) { return (1); }

  retval = CVodeSVtolerances(cvode_mem, 1e-5, abstol);
  if (check_retval(&retval, "CVodeSVtolerances", 1)) { return (1); }

  A = SUNDenseMatrix(NEQ, NEQ, sunctx);
  if (check_retval((void*)A, "SUNDenseMatrix", 0)) { return (1); }

  LS = SUNLinSol_Dense(sp_param_vec, A, sunctx);
  if (check_retval((void*)LS, "SUNLinSol_Dense", 0)) { return (1); }

int userdata_retval=CVodeSetUserData(cvode_mem, &coll_params);
if (check_retval(&userdata_retval, "CVodeSetUserData", 1)) { return (1); }


  retval = CVodeSetLinearSolver(cvode_mem, LS, A);
  if (check_retval(&retval, "CVodeSetLinearSolver", 1)) { return (1); }

  retval = CVodeSetJacFn(cvode_mem, Jac_coll);
  if (check_retval(&retval, "CVodeSetJacFn", 1)) { return (1); }

  FID = fopen("/home/robbie-ellis/VSCode/new active_sterile/exported/sundials_coll.txt", "w");
  if(FID==NULL)
  {
      printf("Error opening file\n");
      return 1;
  }
  iout = 0;
  tout = T0+0.1;
  while (1)
  {
    retval = CVode(cvode_mem, tout, sp_param_vec, &t, CV_NORMAL);
    printf("Here\n");
    PrintOutput(t, Ith(sp_param_vec, 1), Ith(sp_param_vec, 2), Ith(sp_param_vec, 3), Ith(sp_param_vec, 4));

  fprintf(FID,"%.10e, %.10e, %.10e, %.10e, %.10e\n", t, Ith(sp_param_vec, 1), Ith(sp_param_vec, 2), Ith(sp_param_vec, 3), Ith(sp_param_vec, 4));
    if (check_retval(&retval, "CVode", 1)) { break; }
    if (retval == CV_SUCCESS)
    {
      iout++;
      tout += 50;
    }

    retval = CVodePrintAllStats(cvode_mem, FID, SUN_OUTPUTFORMAT_CSV);

    if (iout == 1000) { break; }
  }
  fclose(FID);
}




static int f_coll(sunrealtype t, N_Vector sp_params_vec, N_Vector d_sp_params, void* user_data)
{
struct coll_parameters* params=(struct coll_parameters*) user_data;
struct sp_params sp_params;
//printf("T: %f\n",params->T);
//printf("Tphi: %f\n",Ith(sp_params_vec,3));
//printf("cphi: %f\n",Ith(sp_params_vec,4));
sp_params.Ts=Ith(sp_params_vec,1);
sp_params.c_s=Ith(sp_params_vec,2);
sp_params.Tphi=Ith(sp_params_vec,3);
sp_params.c_phi=Ith(sp_params_vec,4);
sp_params.T=params->T;
struct final_derivs derivs=step_omp_coll(sp_params,params->fixed_params,params->all_tables);
Ith(d_sp_params,1)=derivs.dTs;
Ith(d_sp_params,2)=derivs.dcs;
Ith(d_sp_params,3)=derivs.dTp;
Ith(d_sp_params,4)=derivs.dcp;
//printf("dT: %f\n",derivs.dTs);
//printf("dTphi: %f\n",derivs.dTp);
return 0;
}




static int Jac_coll(sunrealtype t, N_Vector sp_param_vec, N_Vector fy, SUNMatrix J,
               void* user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
{
struct coll_parameters* params=(struct coll_parameters*) user_data;
struct fixed_params* fixed_params=params->fixed_params;
struct interpolation_tables* all_tables=params->all_tables;
struct sp_params sp_params;
double T=params->T;
sp_params.T=T;
sp_params.Ts=Ith(sp_param_vec,1);
sp_params.c_s=Ith(sp_param_vec,2);
sp_params.Tphi=Ith(sp_param_vec,3);
sp_params.c_phi=Ith(sp_param_vec,4);
struct final_derivs derivs_x0=step_omp_coll(sp_params,fixed_params,all_tables);

double h=1e-3;

struct sp_params sp_params_Ts;
sp_params_Ts.T=T;
sp_params_Ts.Ts=Ith(sp_param_vec,1)+h;
sp_params_Ts.c_s=Ith(sp_param_vec,2);
sp_params_Ts.Tphi=Ith(sp_param_vec,3);
sp_params_Ts.c_phi=Ith(sp_param_vec,4);
struct final_derivs derivs_Ts=step_omp_coll(sp_params_Ts,fixed_params,all_tables);

struct sp_params sp_params_cs;
sp_params_cs.T=T;
sp_params_cs.Ts=Ith(sp_param_vec,1);
sp_params_cs.c_s=Ith(sp_param_vec,2)+h;
sp_params_cs.Tphi=Ith(sp_param_vec,3);
sp_params_cs.c_phi=Ith(sp_param_vec,4);
struct final_derivs derivs_cs=step_omp_coll(sp_params_cs,fixed_params,all_tables);

struct sp_params sp_params_Tphi;
sp_params_Tphi.T=T;
sp_params_Tphi.Ts=Ith(sp_param_vec,1);
sp_params_Tphi.c_s=Ith(sp_param_vec,2);
sp_params_Tphi.Tphi=Ith(sp_param_vec,3)+h;
sp_params_Tphi.c_phi=Ith(sp_param_vec,4);
struct final_derivs derivs_Tphi=step_omp_coll(sp_params_Tphi,fixed_params,all_tables);

struct sp_params sp_params_cphi;
sp_params_cphi.T=T;
sp_params_cphi.Ts=Ith(sp_param_vec,1);
sp_params_cphi.c_s=Ith(sp_param_vec,2);
sp_params_cphi.Tphi=Ith(sp_param_vec,3);
sp_params_cphi.c_phi=Ith(sp_param_vec,4)+h;
struct final_derivs derivs_cphi=step_omp_coll(sp_params_cphi,fixed_params,all_tables);



  IJth(J, 1, 1) = (derivs_Tphi.dTs-derivs_x0.dTs)/h;
  IJth(J, 1, 2) = (derivs_Tphi.dcs-derivs_x0.dcs)/h;
  IJth(J, 1, 3) = (derivs_Tphi.dTp-derivs_x0.dTp)/h;
  IJth(J, 1, 4) = (derivs_Tphi.dcp-derivs_x0.dcp)/h;

  IJth(J, 2, 1) = (derivs_cs.dTs-derivs_x0.dTs)/h;
    IJth(J, 2, 2) = (derivs_cs.dcs-derivs_x0.dcs)/h;
    IJth(J, 2, 3) = (derivs_cs.dTp-derivs_x0.dTp)/h;
    IJth(J, 2, 4) = (derivs_cs.dcp-derivs_x0.dcp)/h;

    IJth(J, 3, 1) = (derivs_Tphi.dTs-derivs_x0.dTs)/h;
    IJth(J, 3, 2) = (derivs_Tphi.dcs-derivs_x0.dcs)/h;
    IJth(J, 3, 3) = (derivs_Tphi.dTp-derivs_x0.dTp)/h;
    IJth(J, 3, 4) = (derivs_Tphi.dcp-derivs_x0.dcp)/h;

    IJth(J, 4, 1) = (derivs_cphi.dTs-derivs_x0.dTs)/h;
    IJth(J, 4, 2) = (derivs_cphi.dcs-derivs_x0.dcs)/h;
    IJth(J, 4, 3) = (derivs_cphi.dTp-derivs_x0.dTp)/h;
    IJth(J, 4, 4) = (derivs_cphi.dcp-derivs_x0.dcp)/h;

  return (0);
}
