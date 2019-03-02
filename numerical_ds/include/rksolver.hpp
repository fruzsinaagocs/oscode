#pragma once
#include "system.hpp"

class RKSolver
{
    private: 
    // Frequency and friction term
    std::function<std::complex<double>(double)> w;
    std::function<std::complex<double>(double)> g;
    

    // Butcher tablaus
    Eigen::Matrix<double,5,5> butcher_a5;
    Eigen::Matrix<double,3,3> butcher_a4;
    Eigen::Matrix<double,6,1> butcher_b5, butcher_c5;
    Eigen::Matrix<double,4,1> butcher_b4, butcher_c4;
   
    public:
    // constructors
    RKSolver();
    RKSolver(de_system &de_sys);
    // grid of ws, gs
    Eigen::Matrix<std::complex<double>,6,1> ws, gs;
    Eigen::Matrix<std::complex<double>,5,1> ws5, gs5;

    // single step function 
    Eigen::Matrix<std::complex<double>,2,4> step(std::complex<double>, std::complex<double>, double, double);
    // ODE to solve
    Eigen::Matrix<std::complex<double>,1,4> f(double t, const Eigen::Matrix<std::complex<double>,1,4> &y);  

};

RKSolver::RKSolver(){
};

RKSolver::RKSolver(de_system &de_sys){
    
    // Set frequency and friction terms
    RKSolver::w = de_sys.w;
    RKSolver::g = de_sys.g;
    
    // Set Butcher tableaus
    RKSolver::butcher_a5 << 0.1174723380352676535740,0,0,0,0,
                 -0.1862479800651504276304,0.5436322218248278794734,0,0,0,
                 -0.6064303885508280518989,1,0.2490461467911506000559,0,0,
                 2.899356540015731406420,-4.368525611566240669139,2.133806714786316899991,0.2178900187289247091542,0,
                 18.67996349995727204273,-28.85057783973131956546,10.72053408420926869789,1.414741756508049078612,-0.9646615009432702537787;
	RKSolver::butcher_a4 << 0.172673164646011428100,0,0,
                -1.568317088384971429762,2.395643923738960001662,0,
                -8.769507466172720011410,10.97821961869480000808,-1.208712152522079996671;
	RKSolver::butcher_b5 << 0.1127557227351729739820,0,0.5065579732655351768382,0.04830040376995117617928,0.3784749562978469803166,-0.04608905606850630731611;
	RKSolver::butcher_c5 << 0,0.117472338035267653574,0.357384241759677451843,0.642615758240322548157,0.882527661964732346426,1;
	RKSolver::butcher_b4 << -0.08333333333333333333558,0.5833333333333333333357,0.5833333333333333333356,-0.08333333333333333333558;
	RKSolver::butcher_c4 << 0,0.172673164646011428100,0.827326835353988571900,1;
};

Eigen::Matrix<std::complex<double>,1,4> RKSolver::f(double t, const Eigen::Matrix<std::complex<double>,1,4> &y){
    
    std::complex<double> wi = w(t);
    std::complex<double> gi = g(t);
    Eigen::Matrix<std::complex<double>,1,4> result;
    result << y[1], -wi*wi*y[0]-2.0*gi*y[1], wi, gi;
    return result;
};

Eigen::Matrix<std::complex<double>,2,4> RKSolver::step(std::complex<double> x0, std::complex<double> dx0, double t0, double h){

    Eigen::Matrix<std::complex<double>,1,4> y0, y, y4, y5, delta, k5_i, k4_i;
    y4 = Eigen::Matrix<std::complex<double>,1,4>::Zero();
    y0 << x0, dx0, 0.0, 0.0;
    y5 = y4;
    // TODO: resizing of ws5, gs5, insertion
    Eigen::Matrix<std::complex<double>,6,4> k5;
    Eigen::Matrix<std::complex<double>,4,4> k4;
    Eigen::Matrix<std::complex<double>,2,4> result;
    k5.row(0) = h*f(t0, y0);
    ws(0) = k5.row(0)(2)/h;
    gs(0) = k5.row(0)(3)/h;
    for(int s=1; s<=5; s++){
        y = y0;
        for(int i=0; i<=(s-1); i++)
            y += butcher_a5(s-1,i)*k5.row(i);
        k5_i = h*f(t0 + butcher_c5(s)*h, y);
        k5.row(s) = k5_i;
        ws(s) = k5_i(2)/h;
        gs(s) = k5_i(3)/h;
    }
    k4.row(0) = k5.row(0);
    ws5(0) = ws(0);
    gs5(0) = gs(0);
    for(int s=1; s<=3; s++){
       y = y0;
       for(int i=0; i<=(s-1); i++)
           y += butcher_a4(s-1,i)*k4.row(i);
        k4_i = h*f(t0 + butcher_c4(s)*h, y);
        k4.row(s) = k4_i;
        ws5(s) = k4_i(2)/h;
        gs5(s) = k4_i(3)/h;
    }
    for(int j=0; j<=5; j++)
        y5 += butcher_b5(j)*k5.row(j);
    for(int j=0; j<=3; j++)
        y4 += butcher_b4(j)*k4.row(j);
    delta = y5 - y4;
    y5 += y0;
    result << y5, delta;
    // Add in missing w, g at t+h/2
    ws5(4) = ws5(3);
    ws5(3) = ws5(2);
    ws5(2) = w(t0 + h/2);
    gs5(4) = gs5(3);
    gs5(3) = gs5(2);
    gs5(2) = g(t0 + h/2);
    return result;
};

