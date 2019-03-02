#include <math.h>
#include <nag.h>
#include <nag_stdlib.h>
#include <nagd02.h>
#include <iostream>
#include <list>
#include <fstream>
#include <Eigen/Dense>
#include <chrono>
#include "interpolator.hpp"
#include "system.hpp"
#include "solver.hpp"
#include "mpi.h"

int nv=2;
double m=4.51e-6;
double k;
double phi_p=23.293;
double mp=1;
LinearInterpolator<double, std::complex<double>> winterp, ginterp;

double hd(double ki, const double *ybg, const double *dybg, std::complex<double> rk1, std::complex<double> rk2){
    // Initial conditions for the perturbations
    
    double z = ybg[1]*ybg[2]/ybg[3];
    double dz_z = ybg[3] + dybg[1]/ybg[1] - dybg[3]/ybg[3];
    double a = -1.0/(z*std::pow(2.0*ki,0.5)*100*ki);
    std::complex<double> b = std::complex<double>(-dz_z*10.0/ki*a, -10.0/ybg[2]*a);
    //std::cout << "ybg: " << ybg[0] << " " << ybg[1] << " " << ybg[2] << " " << ybg[3] << std::endl;
    //std::cout << a << ", "<< b << std::endl;
    return std::pow(std::abs(a*rk1 + b*rk2),2)*std::pow(ki,3)/(2.0*M_PI*M_PI);
};

double rst(double ki, const double *ybg, const double *dybg, std::complex<double> rk1, std::complex<double> rk2){
    // Initial conditions for the perturbations
    double z = ybg[1]*ybg[2]/ybg[3];
    double a = -1.0/(z*std::pow(2.0*ki,0.5)*100*ki);
    std::complex<double> b = std::complex<double>(0.0, -10.0/ybg[2]*a);
    return std::pow(std::abs(a*rk1 + b*rk2),2)*std::pow(ki,3)/(2.0*M_PI*M_PI);
};

std::complex<double> w(double t){
    return k*std::exp(winterp(t));
};

std::complex<double> g(double t){
    return ginterp(t); 
};

Eigen::Matrix<double,1,4> background(double t){
    // Gives background at a given time t
    Eigen::Matrix<double,1,4> result;
    result << phi_p - std::sqrt(2.0/3.0)*mp*std::log(t),
    -std::sqrt(2.0/3.0)*mp/t, std::pow(t, 1.0/3.0), 1.0/(3.0*t); 
    return result;
};

double V(double phi){
    return std::pow(m,2)*std::pow(phi,nv);
};

double dV(double phi){
    return std::pow(m,2)*nv*std::pow(phi,nv-1);
};

static void NAG_CALL f(double t, Integer n, const double *y, double *yp,
Nag_Comm *comm){
    yp[0] = y[1];
    yp[1] = (-3.0*y[1]*std::sqrt(1.0/(3*std::pow(mp, 2))*(0.5*std::pow(y[1], 2)
    + V(y[0]))) - dV(y[0]));
    yp[2] = y[2]*y[3];
    yp[3] = (-1.0/(3*std::pow(mp, 2))*(std::pow(y[1], 2) - V(y[0])) -
    std::pow(y[3], 2));
};

int main(){


    // Range of wavenumbers
    int nk=3000;
    Eigen::VectorXd ks=Eigen::VectorXd::LinSpaced(nk,-3,5);
    for(int i=0; i<nk; i++)
        ks(i) = std::pow(10,ks(i));
    Eigen::VectorXd exits=Eigen::VectorXd::Zero(nk);
    
    // These will contain t,w,g
    int N = round(1e6);
    Eigen::VectorXd ts=Eigen::VectorXd::Zero(N);
    Eigen::VectorXcd logws=Eigen::VectorXcd::Zero(N), gs=Eigen::VectorXcd::Zero(N);
    
    std::cout << "Starting background calculation" << std::endl;
    // Use the NAG library to solve for the scale factor a(t) (and phi(t),
    // dphi(t), H(t))
    Integer liwsav, lrwsav, exit_status=0, npts=round(N-1), n=4;
    double tgot, tol, tnext, twant, tstart=1.0, ti=1e4, tend=1e6, tinc=(tend-tstart)/npts;
    double *rwsav=0, *thresh=0, *ygot=0, *yinit=0, *dyinit=0, *ybg=0, *dybg=0, *ymax=0;
    double *ypgot=0;
    Integer *iwsav=0;
    NagError fail;
    Nag_RK_method method;
    Nag_ErrorAssess errass;
    Nag_Comm comm;
    INIT_FAIL(fail);
    liwsav = 130;
    lrwsav = 350 + 32 * n;
    thresh = NAG_ALLOC(n, double);
    ygot = NAG_ALLOC(n, double);
    yinit = NAG_ALLOC(n, double);
    dyinit = NAG_ALLOC(n, double);
    ybg = NAG_ALLOC(n, double);
    dybg = NAG_ALLOC(n, double);
    ypgot = NAG_ALLOC(n, double);
    ymax = NAG_ALLOC(n, double);
    iwsav = NAG_ALLOC(liwsav, Integer);
    rwsav = NAG_ALLOC(lrwsav, double);
    method = (Nag_RK_method) nag_enum_name_to_value("Nag_RK_7_8");
    errass = (Nag_ErrorAssess) nag_enum_name_to_value("Nag_ErrorAssess_off");
    
    Eigen::Matrix<double,1,4> y0 = background(tstart);
    yinit[0] = y0(0);
    yinit[1] = y0(1);
    yinit[2] = y0(2);
    yinit[3] = y0(3);
    f(tstart,n,yinit,dyinit,&comm);
    thresh[0] = 1.0e-14;
    thresh[1] = 1.0e-14;
    thresh[2] = 1.0e-14;
    thresh[3] = 1.0e-14;
    tol = 1.0e-14;
    ts(0) = tstart;
    logws(0) = -std::log(yinit[2]);
    gs(0) = 1.5*yinit[3] + dyinit[1]/yinit[1] - dyinit[3]/yinit[3];
    nag_ode_ivp_rkts_setup(n, tstart, tend, yinit, tol, thresh, method, errass, 0.0, iwsav, rwsav, &fail);
    twant=tstart;
    int start=0;
    double prev_horizon = 100.0/dyinit[2];
    double next_horizon = prev_horizon; 
    for(int i=0; i<(npts-2); i++){
        tnext = twant;
        if(twant+tinc > ti and twant < ti)
            twant = ti;
        else
            twant+=tinc; 
        while(tnext < twant){
            nag_ode_ivp_rkts_range(f,4,twant,&tgot,ygot,ypgot,ymax,&comm,iwsav,rwsav,&fail);
            tnext = tgot; 
            };
        ts(i) = twant;
        logws(i) = -std::log(ygot[2]);
        gs(i) = 1.5*ygot[3] + ypgot[1]/ygot[1] - ypgot[3]/ygot[3];
        next_horizon = 100.0/ypgot[2];
        if(std::abs(twant-ti)<1e-6){
            ybg[0] = ygot[0];
            ybg[1] = ygot[1];
            ybg[2] = ygot[2];
            ybg[3] = ygot[3];
            f(ti,n,ybg,dybg,&comm);
        };
        if(next_horizon < prev_horizon){
            for(int j=start; j<nk; j++){
                if(100.0/ypgot[2] < 1.0/ks(j)){
                    exits(j) = ts(i);
                    start=j+1;
                    break;
                };
            };
        };
        prev_horizon = next_horizon;

    };
        
    NAG_FREE(thresh);
    NAG_FREE(yinit);
    NAG_FREE(ygot);
    NAG_FREE(ypgot);
    NAG_FREE(ymax);
    NAG_FREE(rwsav);
    NAG_FREE(iwsav);
    NAG_FREE(dyinit);
    std::cout << "Done solving background" << std::endl;

    // Construct interpolating functions w,g from grids
    std::cout << "Making w(t), g(t)" << std::endl;
    for(int i=0; i<N; i++){
        winterp.insert(ts(i),logws(i));
        ginterp.insert(ts(i),gs(i));
    }; 
    de_system system(&w,&g);
    std::cout << "Done making w(t), g(t)" << std::endl;
    
    // Solve the evolution of each perturbation
    double tf, rtol, atol, h0;
    std::complex<double> x0, dx0;
    int order=3;
    bool full_output=false;
    rtol=1e-4;
    atol=0.0;
    h0=1.0;
    
    std::list<std::complex<double>> rk1, rk2;
    std::list<double> times;
    double t1,t2;   

    for(int i=0; i<nk; i++){
        k = ks(i);
        tf = exits(i); 
        dx0 = 0.0;
        x0 = 100.0*k;
        
        Solution solution1(system, x0, dx0, ti, tf, order, rtol, atol, h0, full_output);
        t1 = MPI_Wtime();
        solution1.solve();
        rk1.push_back(solution1.sol.back());
        x0 = 0.0;
        dx0 = 10.0*k*k;
        Solution solution2(system, x0, dx0, ti, tf, order, rtol, atol, h0, full_output);
        solution2.solve();
        t2 = MPI_Wtime();
        rk2.push_back(solution2.sol.back());
        times.push_back(t2-t1);

    };
    
    auto it_1 = rk1.begin();
    auto it_2 = rk2.begin();
    auto it_3 = times.begin();
    while(it_1!=rk1.end() && it_2!=rk2.end() && it_3!=times.end()){
        //std::cout << "rk1: " << *it_1 << ", rk2: " << *it_2 << std::endl;
        ++it_1;
        ++it_2;
        ++it_3;
    };

    // Write PPS to file
    std::ofstream f;
    f.open("test/ms/pps-timed.txt");
    it_1 = rk1.begin();
    it_2 = rk2.begin();
    it_3 = times.begin();
    for(int k=0; k<nk; k++){
        f << ks(k) << ", " << hd(ks(k),ybg,dybg,*it_1,*it_2) << ", " << rst(ks(k),ybg,dybg,*it_1,*it_2) << ", " << *it_3 << std::endl;
        ++it_1;
        ++it_2;
        ++it_3;
    };
    f.close();

    return 0;
}
