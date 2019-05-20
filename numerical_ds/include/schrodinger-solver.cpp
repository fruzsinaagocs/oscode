#include "interpolator.hpp"
#include "solver.hpp"
//#include <vector>
#include <cmath>
#include <chrono>
#include <boost/math/special_functions/hermite.hpp>
#include <fstream>
#include <string>
#include <stdlib.h>

int n=1;
double m=1, K=2, E=std::sqrt(K/m)*(n-0.5)*(1.0);

std::complex<double> g(double t){
    return 0.0;
};

double V(double t){
// quantum mechanical potential well
    return 0.5*K*t*t;
};

std::complex<double> w(double t){
    std::complex<double> result = std::sqrt(std::complex<double>(2*m*(E-V(t))));
    return result;
};

Eigen::Vector2cd ic(double t){

    double gam = std::sqrt(m*K);
    Eigen::Vector2cd result;
    double A = std::pow(gam/M_PI,1/4.0)*1.0/std::sqrt(std::pow(2.0,n-1)*tgamma(n));
    result << boost::math::hermite(n-1,std::sqrt(gam)*t)*std::exp(-0.5*gam*t*t), -gam*t*boost::math::hermite(n-1,std::sqrt(gam)*t)*std::exp(-0.5*gam*t*t) + 2*(n-1)*boost::math::hermite(n-2,std::sqrt(gam)*t)*std::sqrt(gam)*std::exp(-0.5*gam*t*t);
    return A*result;
};

int main(){
     
    // Example with w(t), g(t) analytically given
    std::ofstream f;
    std::complex<double> x0, dx0;
    double ti, tf, rtol, atol, h0;
    bool full_output = false;//true;
    int order = 3;   
    de_system sys(&w, &g);
    ti = (-std::sqrt((n-0.5)*2/std::sqrt(K*m))-2.0);
    tf = 0.5;
    Eigen::Vector2cd bcL=ic(ti), bcR=ic(-ti);
    rtol = 1e-4;
    atol = 0.0;
    h0 = 1.0;
    x0 = ic(ti)(0);
    dx0 = ic(ti)(1);
    // Left
    Solution solutionL(sys,x0,dx0,ti,tf,order,rtol,atol,h0,full_output); 
    solutionL.solve();
    // Right
    x0 = ic(-ti)(0);
    dx0 = ic(-ti)(1);
    Solution solutionR(sys,x0,dx0,-ti,tf,order,rtol,atol,-h0,full_output); 
    solutionR.solve();

    std::complex<double> xL,xR,dxL,dxR;
    xL = solutionL.sol.back();
    dxL = solutionL.dsol.back();
    xR = solutionR.sol.back();
    dxR = solutionR.dsol.back();
    std::cout << std::setprecision(7) << "xL, dxL, xR, dxR: " << xL << ", " << dxL << ", " << xR << ", " << dxR << std::endl;
    std::cout << "Difference at the middle is: " << dxL/xL - dxR/xR << std::endl;
//    a = (bcf(0)*x2b - bcb(0)*x2f)/(x1f*x2b - x1b*x2f);
//    b = (bcf(0)*x1b - bcb(0)*x1f)/(x2f*x1b - x2b*x1f);
    
};


