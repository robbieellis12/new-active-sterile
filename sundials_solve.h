#ifndef SUNDIALS_SOLVE_H
#define SUNDIALS_SOLVE_H
#include "sn_math.h"
#include "sn_cross_sections.h"
#include "sn_computations.h"
#include "tests.h"
#include "diffeqs.h"
#include <sundials/sundials_core.h>
#include <cvode/cvode.h>
#include <nvector/nvector_serial.h>
#include <sunlinsol/sunlinsol_dense.h>       
#include <sunmatrix/sunmatrix_dense.h>

int solve_de(double T0, struct sp_params, struct dx_params* dx_params);
static int f(sunrealtype t, N_Vector y, N_Vector ydot, void* user_data);
 int check_retval(void* returnvalue, const char* funcname, int opt);
 void PrintOutput(sunrealtype t, sunrealtype y1, sunrealtype y2,
                        sunrealtype y3, sunrealtype y4);
static int Jac(sunrealtype T, N_Vector sp_param_vec, N_Vector fy, SUNMatrix J,
               void* user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);
static int f_eq(sunrealtype t, N_Vector sp_params_vec, N_Vector d_sp_params, void* user_data);
static int Jac_eq(sunrealtype T, N_Vector sp_param_vec, N_Vector fy, SUNMatrix J,
               void* user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);
int solve_de_eq(double T0, struct sp_params sp_params_init,struct dx_eq_params* dx_params);

#endif