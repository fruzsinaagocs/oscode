#include <solver.hpp>
//#include "catch.hpp"
#include <math.h>
#include <nag.h>
#include <nag_stdlib.h>
#include <nagd02.h>

using namespace RKWKB;

extern int n;
extern double m, k, phi_p, mp;

Vector background(Scalar);

// z = [phi, dphi, a, H]
//TODO: function prototypes here?
// avoid using globals?

Vector background(Scalar t){
    // Gives background at a given time t
    Vector result(4);
    result << phi_p - std::sqrt(2.0/3.0)*mp*std::log(t), -std::sqrt(2.0/3.0)*mp/t, std::pow(t, 1.0/3.0), 1.0/(3.0*t); 
    return result;
};

Scalar V(Scalar phi){
    return std::pow(m,2)*std::pow(phi,n);
};

Scalar dV(Scalar phi){
    return std::pow(m,2)*n*std::pow(phi,n-1);
};

Scalar ddV(Scalar phi){
    if(n<2)
        return 0.0;
    else
        return std::pow(m,2)*n*(n-1)*std::pow(phi, n-2);
};

Scalar dddV(Scalar phi){
    if(n<3)
        return 0.0;
    else
        return std::pow(m,2)*n*(n-1)*(n-2)*std::pow(phi, n-3);
};

Vector F(Vector z){
    Vector result(z.size());
    result << z(1), -3.0*z(1)*std::sqrt(1.0/(3*std::pow(mp, 2))*(0.5*std::pow(z(1), 2) + V(z(0)))) - dV(z(0)), z(2)*z(3), -1.0/(3*std::pow(mp, 2))*(std::pow(z(1), 2) - V(z(0))) - std::pow(z(3), 2);
    return result;
};

Matrix DF(Vector z){
    Matrix result(z.size(),z.size());
    result << 0.0, 1.0, 0.0, 0.0,
              -3.0/2.0*z(1)*std::pow(1.0/(3.0*std::pow(mp,2))*(0.5*std::pow(z(1),2)+V(z(0))),-0.5)*1.0/(3.0*std::pow(mp,2))*dV(z(0))-ddV(z(0)), -3.0*std::pow(1.0/(3.0*std::pow(mp,2))*(0.5*std::pow(z(1),2)+V(z(0))),0.5)-3.0/2.0*z(1)*std::pow(1.0/(3.0*std::pow(mp,2))*(0.5*std::pow(z(1),2)+V(z(0))),-0.5)*1.0/(3.0*std::pow(mp,2))*z(1), 0.0, 0.0,
              0.0, 0.0, z(3), z(2),
              1.0/(3.0*std::pow(mp,2))*dV(z(0)), -2.0/(3.0*std::pow(mp,2))*z(1), 0.0, -2.0*z(3);
    return result;
};

Scalar w(Vector z){
    return k/z(2);
};

Vector Dw(Vector z){
    Vector result(z.size());
    result << 0.0, 0.0, -k/std::pow(z(2),2), 0.0;
    return result;
}

Matrix DDw(Vector z){
    Matrix result(z.size(), z.size());
    result << 0.0, 0.0, 0.0, 0.0,
              0.0, 0.0, 0.0, 0.0,
              0.0, 0.0, 2.0*k/std::pow(z(2),3), 0.0,
              0.0, 0.0, 0.0, 0.0;
    return result;
}

Scalar g(Vector z){
    Vector dz=F(z);
    return 1.5*z(3) + dz(1)/z(1) - dz(3)/z(3);;
};

Vector Dg(Vector z){
    Vector result(z.size());
    result << -3.0/2.0*std::pow(1.0/(3.0*std::pow(mp,2))*(0.5*std::pow(z(1),2) + V(z(0))),-0.5)*1.0/(3.0*std::pow(mp,2))*dV(z(0)) - ddV(z(0))/z(1) - 1.0/(3.0*std::pow(mp,2)*z(3))*dV(z(0)), -3.0/2.0*std::pow(1.0/(3.0*std::pow(mp,2))*(0.5*std::pow(z(1),2)+V(z(0))),-0.5)*1.0/(3.0*std::pow(mp,2))*z(1) + dV(z(0))/std::pow(z(1),2) + 2.0/(3.0*std::pow(mp,2)*z(3))*z(1), 0.0, 5.0/2.0 ;
    return result;
};

Matrix DDg(Vector z){
    Matrix result(z.size(), z.size());
    result << -3.0/2.0*std::pow(1.0/(3.0*std::pow(mp,2))*(0.5*std::pow(z(1),2) + V(z(0))),-0.5)*1.0/(3*std::pow(mp,2))*ddV(z(0)) + 3.0/4.0*std::pow(1.0/(3.0*std::pow(mp,2))*(0.5*std::pow(z(1),2) + V(z(0))),-1.5)*std::pow(1.0/(3.0*std::pow(mp,2))*dV(z(0)),2) -dddV(z(0))/z(1) - ddV(z(0))/(3.0*std::pow(mp,2)*z(3)), 3.0/4.0*std::pow(1.0/(3.0*std::pow(mp,2))*(0.5*std::pow(z(1),2) + V(z(0))),-1.5)*std::pow(1.0/(3.0*std::pow(mp,2)),2)*dV(z(0))*z(1) + ddV(z(0))/std::pow(z(1),2), 0.0, dV(z(0))/(3.0*std::pow(mp,2)*std::pow(z(3),2)),
              3.0/4.0*std::pow(1.0/(3.0*std::pow(mp,2))*(0.5*std::pow(z(1),2) + V(z(0))),-1.5)*std::pow(1.0/(3.0*std::pow(mp,2)),2)*dV(z(0))*z(1) + ddV(z(0))/std::pow(z(1),2), -3.0/2.0*std::pow(1.0/(3.0*std::pow(mp,2))*(0.5*std::pow(z(1),2) + V(z(0))),-0.5)*1.0/(3.0*std::pow(mp,2)) + 3.0/4.0*std::pow(1.0/(3.0*std::pow(mp,2))*(0.5*std::pow(z(1),2) + V(z(0))),-1.5)*std::pow(1.0/(3.0*std::pow(mp,2))*z(1),2) - 2.0*dV(z(0))/std::pow(z(1),3) + 2.0/(3.0*std::pow(mp,2)*z(3)), 0.0, -2.0/(3.0*std::pow(mp,2)*std::pow(z(3),2)),
             0.0, 0.0, 0.0, 0.0,
             dV(z(0))/(3.0*std::pow(mp,2)*std::pow(z(3),2)), -2.0/(3.0*std::pow(mp,2)*std::pow(z(3),2)), 0.0, 0.0;
    return result;
};

Scalar f_end(Vector y, Scalar t){
    // A function to end integration of the ODE when f_end=0.
    double t_f = 100000.0;
    return (t - t_f);
};

static void NAG_CALL f(double t, Integer n, const double *y, double *yp, Nag_Comm *comm){
    // system of differential equations to integrate, compatible with NAG types.
    
    
    yp[0] = y[1];
    yp[1] = (-2.0*(1.5*y[5] + (-3.0*y[3]*std::sqrt(1.0/(3*std::pow(mp, 2))*(0.5*std::pow(y[3], 2) + V(y[2]))) - dV(y[2]))/y[3] - (-1.0/(3*std::pow(mp, 2))*(std::pow(y[3], 2) - V(y[2])) - std::pow(y[5], 2))/y[5])*y[1] - std::pow(k/y[4],2)*y[0]).real();
    yp[2] = y[3];
    yp[3] = (-3.0*y[3]*std::sqrt(1.0/(3*std::pow(mp, 2))*(0.5*std::pow(y[3], 2) + V(y[2]))) - dV(y[2])).real();
    yp[4] = y[4]*y[5];
    yp[5] = (-1.0/(3*std::pow(mp, 2))*(std::pow(y[3], 2) - V(y[2])) - std::pow(y[5], 2)).real();

};
