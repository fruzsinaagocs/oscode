#include <math.h>
#include <nag.h>
#include <nag_stdlib.h>
#include <nagd02.h>
#include <iostream>
#include <list>
#include <fstream>
#include <Eigen/Dense>
#include <chrono>
#include <valgrind/callgrind.h>
#include <boost/math/special_functions/bessel_prime.hpp>
#include <boost/math/special_functions/bessel.hpp>
#include "interpolator.hpp"
#include "system.hpp"
#include "solver.hpp"
#include "mpi.h"
int nv=2;
double m=4.51e-6;
double k;
double phi_p=23.293;
double mp=1;
double alpha=0.005;
double phi_s=22.50;
double delta_phi=0.03;
int N=round(1e6), npts=round(N-1);
double tstart=1.0, tend=1e6, tinc=(tend-tstart)/npts;
Eigen::VectorXd logws, listgs;

std::complex<double> hankel1_0(double x){
    return std::complex<double>(boost::math::cyl_bessel_j(0,x), boost::math::cyl_neumann(0,x));
};

std::complex<double> hankel1_0_prime(double x){
    return std::complex<double>(boost::math::cyl_bessel_j_prime(0,x), boost::math::cyl_neumann_prime(0,x));
};

double hd(double ki, const double *ybg, const double *dybg, std::complex<double> rk1, std::complex<double> rk2){
    // Initial conditions for the perturbations, Bunch--Davies
    double z = ybg[1]*ybg[2]/ybg[3];
    double dz_z = ybg[3] + dybg[1]/ybg[1] - dybg[3]/ybg[3];
    double a = -1.0/(z*std::pow(2.0*ki,0.5)*100*ki);
    std::complex<double> b = std::complex<double>(-dz_z*10.0/ki*a, -10.0/ybg[2]*a);
    return std::pow(std::abs(a*rk1 + b*rk2),2)*std::pow(ki,3)/(2.0*M_PI*M_PI);
};

double rst(double ki, const double *ybg, const double *dybg, std::complex<double> rk1, std::complex<double> rk2){
    // Initial conditions for the perturbations with renormalised SET
    double z = ybg[1]*ybg[2]/ybg[3];
    double a = -1.0/(z*std::pow(2.0*ki,0.5)*100*ki);
    std::complex<double> b = std::complex<double>(0.0, -10.0/ybg[2]*a);
    return std::pow(std::abs(a*rk1 + b*rk2),2)*std::pow(ki,3)/(2.0*M_PI*M_PI);
};

double kd(double t0, double ki, const double *ybg, const double *dybg, std::complex<double> rk1, std::complex<double> rk2){
    // Initial conditions for the perturbations in Kinetic Dominance at t0
    double z = ybg[1]*ybg[2]/ybg[3];
    double dz_z = ybg[3] + dybg[1]/ybg[1] - dybg[3]/ybg[3];
    std::complex<double> a = 1.0/(100.0*ki*z)*std::pow(3.0*M_PI/8.0,0.5)*std::pow(t0,1.0/3.0)*hankel1_0(3.0/2.0*ki*std::pow(t0,2.0/3.0));
    std::complex<double> b = 1.0/(10.0*ki)*(-dz_z + 1.0/(3.0*t0) + ki*std::pow(t0,-1.0/3.0)*hankel1_0_prime(3.0/2.0*ki*std::pow(t0,2.0/3.0))/hankel1_0(3.0/2.0*ki*std::pow(t0,2.0/3.0)))*a;
    return std::pow(std::abs(a*rk1 + b*rk2),2)*std::pow(ki,3)/(2.0*M_PI*M_PI);
};

double RKSolver::w(double t){
    int i;
    i=int((t-tstart)/tinc);
    
    double logw0 = logws(i);
    double logw1 = logws(i+1);
    return k*std::exp(logw0+(logw1-logw0)*(t-tstart-tinc*i)/tinc);
};

double RKSolver::g(double t){
    int i;
    i=int((t-tstart)/tinc);
    
    double g0 = listgs(i);
    double g1 = listgs(i+1);
    return (g0+(g1-g0)*(t-tstart-tinc*i)/tinc);
};

double win(double){
    return 0.0;
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
    return
    std::pow(m,2)*nv*std::pow(phi,nv-1);
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
    Eigen::VectorXd ks=Eigen::VectorXd::LinSpaced(nk,-3,3);
    for(int i=0; i<nk; i++)
        ks(i) = std::pow(10,ks(i));
    Eigen::VectorXd exits=Eigen::VectorXd::Zero(nk);
    
    // These will contain t,w,g
    
    Eigen::VectorXd ts=Eigen::VectorXd::Zero(N);
    logws=Eigen::VectorXd::Zero(N);
    listgs=Eigen::VectorXd::Zero(N);
    // To log background evolution
    Eigen::VectorXd aHs=Eigen::VectorXd::Zero(N), phis=Eigen::VectorXd::Zero(N);
    std::cout << "Starting background calculation" << std::endl;
    // Use the NAG library to solve for the scale factor a(t) (and phi(t),
    // dphi(t), H(t))
    Integer liwsav, lrwsav, n=4;
    double tgot, tol, tnext, twant, ti=1e4;
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
    // Set initial value of ybg as y0
    ybg[0] = y0(0);
    ybg[1] = y0(1);
    ybg[2] = y0(2);
    ybg[3] = y0(3);
    f(tstart,n,yinit,dyinit,&comm);
    thresh[0] = 1.0e-14;
    thresh[1] = 1.0e-14;
    thresh[2] = 1.0e-14;
    thresh[3] = 1.0e-14;
    tol = 1.0e-14;
    ts(0) = tstart;
    aHs(0)=y0(2)*y0(3);
    phis(0)=y0(0);
    logws(0) = -std::log(yinit[2]);
    listgs(0) = 1.5*yinit[3] + dyinit[1]/yinit[1] - dyinit[3]/yinit[3];
    nag_ode_ivp_rkts_setup(n, tstart, tend, yinit, tol, thresh, method, errass, 0.0, iwsav, rwsav, &fail);
    twant=tstart;
    int start=0;
    double prev_horizon = 100.0/dyinit[2];
    double next_horizon = prev_horizon; 
    for(int i=1; i<(npts-1); i++){
        tnext = twant;
        if(twant+tinc > ti and twant < ti)
            twant = ti;
        twant+=tinc; 
        while(tnext < twant){
            nag_ode_ivp_rkts_range(f,4,twant,&tgot,ygot,ypgot,ymax,&comm,iwsav,rwsav,&fail);
            tnext = tgot; 
            };
        ts(i) = twant;
        aHs(i) = ygot[2]*ygot[3];
        phis(i) = ygot[0];
        logws(i) = -std::log(ygot[2]);
        listgs(i) = 1.5*ygot[3] + ypgot[1]/ygot[1] - ypgot[3]/ygot[3];
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
                    exits(j) = tstart+tinc*i;
                    start=j+1;
                    //std::cout << exits(j) << std::endl;
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

    // Write phi(t) to file:
    std::ofstream fbg;
    fbg.open("test/ms/pps-bg.txt");
    for(int i=0; i<N; i++){
        fbg << ts(i) << ", " << aHs(i) << ", " << phis(i) << std::endl;
    };
    fbg.close();
    std::cout << "Done solving background" << std::endl;
    
    // Construct system to solve
    de_system system(&win,&win);
    
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
    double t1,t2,timing0,timing1;   
    
    timing0 = MPI_Wtime();
    CALLGRIND_START_INSTRUMENTATION;
    CALLGRIND_TOGGLE_COLLECT;
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
    CALLGRIND_TOGGLE_COLLECT;
    CALLGRIND_STOP_INSTRUMENTATION;
    timing1 = MPI_Wtime();
    std::cout << "total time: " << timing1-timing0 << " s." << std::endl;
    
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
    f.open("test/ms/pps-recap");
    it_1 = rk1.begin();
    it_2 = rk2.begin();
    it_3 = times.begin();
    for(int k=0; k<nk; k++){
        f << ks(k) << ", "  << *it_1 << ", " << *it_2 << ", " << hd(ks(k),ybg,dybg,*it_1,*it_2) << ", " << rst(ks(k),ybg,dybg,*it_1,*it_2) << ", " << kd(ti, ks(k),ybg,dybg,*it_1,*it_2) << ", " << *it_3 << std::endl;
        ++it_1;
        ++it_2;
        ++it_3;
    };
    f.close();

    return 0;
}
