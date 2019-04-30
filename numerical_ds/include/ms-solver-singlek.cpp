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

int nv=2;
double m=4.51e-6;
double k=0.001;
double phi_p=23.293;
double mp=1;
LinearInterpolator<double, double> winterp, ginterp;

double w(double t){
    return k*winterp(t);
};

double g(double t){
    return ginterp(t); 
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

static void NAG_CALL g(double t, Integer n, const double *y, double *yp,
Nag_Comm *comm){
    yp[0] = y[1];
    yp[1] = (-2.0*(1.5*y[5] + (-3.0*y[3]*std::sqrt(1.0/(3*std::pow(mp,
    2))*(0.5*std::pow(y[3], 2) + V(y[2]))) - dV(y[2]))/y[3] -
    (-1.0/(3*std::pow(mp, 2))*(std::pow(y[3], 2) - V(y[2])) - std::pow(y[5],
    2))/y[5])*y[1] - std::pow(k/y[4],2)*y[0]);
    yp[2] = y[3];
    yp[3] = (-3.0*y[3]*std::sqrt(1.0/(3*std::pow(mp, 2))*(0.5*std::pow(y[3], 2)
    + V(y[2]))) - dV(y[2]));
    yp[4] = y[4]*y[5];
    yp[5] = (-1.0/(3*std::pow(mp, 2))*(std::pow(y[3], 2) - V(y[2])) -
    std::pow(y[5], 2));
};

int main(){

    // These will contain t,w,g
    int N = round(1e6);
    Eigen::VectorXd ts=Eigen::VectorXd::Zero(N);
    Eigen::VectorXd logws=Eigen::VectorXd::Zero(N), gs=Eigen::VectorXd::Zero(N);
    
    std::cout << "Starting background calculation" << std::endl;
    // Use the NAG library to solve for the scale factor a(t) (and phi(t),
    // dphi(t), H(t))
    Integer liwsav, lrwsav, exit_status=0, npts1=1e3, npts=round(N-1-npts1), n=4;
    double tgot, tol, tnext, twant, tstart=1.0, tend=1e6, tend1=1e2,
    tinc1=(tend1-tstart)/npts1, tinc=(tend-tend1)/npts;
    double *rwsav=0, *thresh=0, *ygot=0, *yinit=0, *dyinit=0, *yend=0, *ymax=0;
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
    yend = NAG_ALLOC(n, double);
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
    logws(0) = 1.0/yinit[2];
    gs(0) = 1.5*yinit[3] + dyinit[1]/yinit[1] - dyinit[3]/yinit[3];
    nag_ode_ivp_rkts_setup(n, tstart, tend, yinit, tol, thresh, method, errass, 0.0, iwsav, rwsav, &fail);
    twant=tstart;
    for(int i=1; i<(npts-1); i++){
        tnext = twant;
        if(i < npts1)
            twant+=tinc1; 
        else
            twant+=tinc;
        while(tnext < twant){
            nag_ode_ivp_rkts_range(f,4,twant,&tgot,ygot,ypgot,ymax,&comm,iwsav,rwsav,&fail);
            tnext = tgot; 
            };
        ts(i) = twant;
        logws(i) = -std::log(ygot[2]);
        gs(i) = 1.5*ygot[3] + ypgot[1]/ygot[1] - ypgot[3]/ygot[3];
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
    //std::cout << "Done making w(t), g(t)" << std::endl;
    
    // Solve the evolution of the perturbation
    auto t1 = std::chrono::high_resolution_clock::now();
    std::cout << "Solving the perturbation" << std::endl;
    double ti, tf, rtol, atol, h0;
    std::complex<double> x0, dx0;
    int order=3;
    bool full_output=true;
    rtol=1e-4;
    atol=0.0;
    h0=1.0;
    ti=1.6;
    tf=5e5;
    x0=0.0;
    dx0=10.0*k*k;
    Solution solution(system, x0, dx0, ti, tf, order, rtol, atol, h0, full_output);
    solution.solve();
    auto t2 = std::chrono::high_resolution_clock::now();
    std::cout << "Done solving the perturbation" << std::endl;
    std::chrono::duration<double,std::milli> t12 = t2-t1;
    std::cout << "Time: " << t12.count() << " ms" << std::endl;

    // Use the NAG library to solve for the perturbation separately
    std::cout << "Starting solving for perturbation with NAG" << std::endl;
    auto t3 = std::chrono::high_resolution_clock::now();
    // First solve up to the start 
    exit_status=0;
    tstart=1.0; tend=ti;   
    INIT_FAIL(fail);
    n = 4;
    liwsav = 130;
    lrwsav = 350 + 32 * n;
    thresh = NAG_ALLOC(n, double);
    ygot = NAG_ALLOC(n, double);
    yinit = NAG_ALLOC(n, double);
    dyinit = NAG_ALLOC(n, double);
    ypgot = NAG_ALLOC(n, double);
    ymax = NAG_ALLOC(n, double);
    iwsav = NAG_ALLOC(liwsav, Integer);
    rwsav = NAG_ALLOC(lrwsav, double);
    method = (Nag_RK_method) nag_enum_name_to_value("Nag_RK_7_8");
    errass = (Nag_ErrorAssess) nag_enum_name_to_value("Nag_ErrorAssess_off");
    
    yinit[0] = y0(0);
    yinit[1] = y0(1);
    yinit[2] = y0(2);
    yinit[3] = y0(3);
    thresh[0] = 1.0e-7;
    thresh[1] = 1.0e-7;
    thresh[2] = 1.0e-7;
    thresh[3] = 1.0e-7;
    tol = 1.0e-6;
    nag_ode_ivp_rkts_setup(n, tstart, tend, yinit, tol, thresh, method, errass, 0.0, iwsav, rwsav, &fail);
    twant = tend;
    tnext = tstart;
    while(tnext < twant){
        nag_ode_ivp_rkts_range(f,4,twant,&tgot,ygot,ypgot,ymax,&comm,iwsav,rwsav,&fail);
        tnext = tgot; 
    };
    yend[0] = ygot[0];
    yend[1] = ygot[1];
    yend[2] = ygot[2];
    yend[3] = ygot[3];
    // For the KD specific case
    if(tstart == ti){
        yend[0] = y0(0);
        yend[1] = y0(1);
        yend[2] = y0(2);
        yend[3] = y0(3);
    };

    NAG_FREE(thresh);
    NAG_FREE(yinit);
    NAG_FREE(ygot);
    NAG_FREE(ypgot);
    NAG_FREE(ymax);
    NAG_FREE(rwsav);
    NAG_FREE(iwsav);
    
    std::cout << "background done with nag" << std::endl;
    // Then from tstart with the perturbation
    exit_status=0; npts=round(1e5);
    tstart=ti; tend=3e5; tinc=(tend-tstart)/npts;
    INIT_FAIL(fail);
    n = 6;
    liwsav = 130;
    lrwsav = 350 + 32 * n;
    thresh = NAG_ALLOC(n, double);
    ygot = NAG_ALLOC(n, double);
    yinit = NAG_ALLOC(n, double);
    dyinit = NAG_ALLOC(n, double);
    ypgot = NAG_ALLOC(n, double);
    ymax = NAG_ALLOC(n, double);
    iwsav = NAG_ALLOC(liwsav, Integer);
    rwsav = NAG_ALLOC(lrwsav, double);
    method = (Nag_RK_method) nag_enum_name_to_value("Nag_RK_7_8");
    errass = (Nag_ErrorAssess) nag_enum_name_to_value("Nag_ErrorAssess_off");
    
    std::list<double> tnag;
    std::list<double> xnag;
    yinit[0] = 0.0;
    yinit[1] = 10.0*k*k;
    yinit[2] = yend[0];
    yinit[3] = yend[1];
    yinit[4] = yend[2]; 
    yinit[5] = yend[3]; 
    thresh[0] = 1.0e-7;
    thresh[1] = 1.0e-7;
    thresh[2] = 1.0e-7;
    thresh[3] = 1.0e-7;
    tnag.push_back(tstart);
    xnag.push_back(yinit[0]);
    std::cout << "y0: " << yinit[0] << " " << yinit[1] << " " << yinit[2] <<
    " "<<  yinit[3] << " " << yinit[4] << " " << yinit[5] << std::endl;
    tol = 1.0e-6;
    nag_ode_ivp_rkts_setup(n, tstart, tend, yinit, tol, thresh, method, errass, 0.0, iwsav, rwsav, &fail);
    twant=tstart;
    for(int i=1; i<(npts-1); i++){
        tnext = twant;
        twant+=tinc; 
        while(tnext < twant){
            nag_ode_ivp_rkts_range(g,6,twant,&tgot,ygot,ypgot,ymax,&comm,iwsav,rwsav,&fail);
            tnext = tgot; 
            };
        tnag.push_back(tgot);
        xnag.push_back(ygot[0]);
    };

    NAG_FREE(thresh);
    NAG_FREE(yinit);
    NAG_FREE(ygot);
    NAG_FREE(ypgot);
    NAG_FREE(ymax);
    NAG_FREE(rwsav);
    NAG_FREE(iwsav);
    auto t4 = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double,std::milli> t34 = t4-t3;
    std::cout << "Done solving for perturbation with NAG" << std::endl;
    std::cout << "Time: " << t34.count() << " ms" << std::endl; 
    // Write results to file
    std::ofstream f;
    f.open("test/ms/nag-ms-kd1.txt");
    auto it_t = tnag.begin();
    auto it_x = xnag.begin();
    while(it_t != tnag.end() || it_x != xnag.end()){
        f << *it_t << ", " << *it_x << std::endl; 
        ++it_t;
        ++it_x;
    };

    return 0;
}
