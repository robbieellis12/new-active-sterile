#include "sn_cross_sections.h"

const double zeta3_sn=1.2020569031595942853;
const double zeta5_sn=1.0369277551433699263;
const double zeta7_sn=1.0083492773819228268;
const double m_pl=1.221e28;
const double VI_2_LINEAR_CONST=0.874367040388; //(1-gamma+ln(pi/2))

// u=s/m_phi^2 -- cross section of va+vs<->va+vs. No constant term of y^4m_phi^2/pi
double c_as_as(double u,double y)
{
    double width=y*y/(8*M_PI);
    double inside_term=-(2-3*u)*(1-u)*log(1+u)+u*(2-4*u+u*u+3*u*u*u)/(1+u);
    return inside_term/((width*width+(1-u)*(1-u))*u*u*M_PI);
}

double c_as_as_num(double u)
{
    double inside_term=-(2-3*u)*(1-u)*log(1+u)+u*(2-4*u+u*u+3*u*u*u)/(1+u);
    return inside_term/(u*u*M_PI);
}

// u=s/m_phi^2 -- cross section of va+va<->va+va and vs+vs<->va+va. Constant term of y^4/2pi
double c_vv_pp(double u,double y)
{
    if(u<4)
    {
        return 0;
    }
    double sqrt_term=sqrt(u*(u-4));
    double arg=sqrt_term/(2-u);
    if(fabs(arg)>1-1e-8)
    {
        return 0;
    }
    return (2*(6-4*u+u*u)/(2-u)*atanh(arg) +3*sqrt_term) /(4*M_PI*u*u);
   // return (2*(6-4*u+u*u)/(2-u)*log((sqrt_term+u-2)/(sqrt_term-u+2))-6*sqrt_term) * y*y*y*y/(4*M_PI*u*u);
}
double c_pp_vv(double u,double y)
{
    if(u<=4)
    {
        return 0;
    }
    return c_vv_pp(u, y)*(u)*u/(u-4);
}
double c_vp_vp(double u,void* y)
{
    if(u<=1)
    {
        return 0;
    }
    return 2*(1/(u*u)+4*u*u*atanh((1-u)*(1-u)/(1+u*u))/((1-u)*(1-u)) -3)*(u-1)/(sqrt(u)*M_PI);
}


//u=s/m_phi^2 -- constant term of y^4/(pi*m_phi^2) omitted
double c_vv_vv(double u, double y)
{
    return 4*(5+3*u)*( u/(1+u) - 2/(2+u)*log(1+u) )/(M_PI*u*u);
    
}
