#include <solver.hpp>
#include <math.h>

using namespace RKWKB;

extern int n;
extern double m, k, phi_p, mp;

Vector background(Scalar);

// z = [phi, dphi, a, H]

Vector background(Scalar t){
    // Gives background at a given time t in kinetic dominance.
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

Scalar inflation_boundaries(Vector y, Scalar t){
    return (-1.0/(3.0*std::pow(mp,2))*(std::pow(y(3),2)-V(y(2)))*y(4));
};

Scalar outside_horizon(Vector y, Scalar t){
    return y(4)*y(5) - 100*k;
};

Scalar until_start(Vector y, Scalar t){
    double tend = 1e4;
    return tend - t; 
};

double HD(double k, double Rk1, double Rk2, Vector ybg, Vector dybg){
    double z = std::real(ybg(1)*ybg(2)/ybg(3));
    std::complex<double> a(-1.0/(z*std::sqrt(2.0*k)*100*k), 0.0);
    std::complex<double> b(0.0, -std::real(a)*k*10/(std::real(ybg(2))*k));
    return std::norm(a*Rk1 + b*Rk2)*std::pow(k, 3)/(2*M_PI*M_PI); 
};

double RST(double k, double Rk1, double Rk2, Vector ybg, Vector dybg){
    double z = std::real(ybg(1)*ybg(2)/ybg(3));
    double dz_z = std::real(dybg(1)/dybg(0) + ybg(3) - dybg(3)/ybg(3));
    std::complex<double> a(-1.0/(z*std::sqrt(2.0*k)*100*k), 0.0);
    std::complex<double> b(-std::real(a)*dz_z*10/k, -std::real(a)*k*10/(std::real(ybg(2))*k));
    return std::norm(a*Rk1 + b*Rk2)*std::pow(k, 3)/(2*M_PI*M_PI); 
};
