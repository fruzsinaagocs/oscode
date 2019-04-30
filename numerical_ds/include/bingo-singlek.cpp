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
double m=7.147378e-6;
double k;
double phi_p=16.5;
double mp=1;
double kpivot=0.05;
double Npivot=50.0;
int Nbg=round(5e5), npts=round(Nbg-1);
double Nstart=0.0, Nend=75.0, Ninc=(Nend-Nstart)/npts;
Eigen::VectorXd logws, listgs;
double ai;

double hd(double ki, double phi, double dphi, double ddphi, double H, double N,
double ai, std::complex<double>
rk1, std::complex<double> rk2){
    // Initial conditions for the perturbations, Bunch--Davies
    double z = ai*std::exp(N)*dphi;
    double dz_z = ddphi/dphi + 1.0;
    double a = 1.0/(z*std::pow(2.0*ki,0.5)*10.0);
    std::complex<double> b = std::complex<double>(-dz_z*a/1e3,
    -ki/(H*ai*std::exp(N)*1e3)*a);
    std::cout << 10*a << ", " << 1e4*b << std::endl;
    return std::pow(std::abs(a*rk1 + b*rk2),2)*std::pow(ki,3)/(2.0*M_PI*M_PI);
};

std::complex<double> w(double N){
    int i;
    i=int((N-Nstart)/Ninc);
    
    double logw0 = logws(i);
    double logw1 = logws(i+1);
    return k*std::exp(logw0+(logw1-logw0)*(N-Nstart-Ninc*i)/Ninc);
};

std::complex<double> g(double N){
    int i;
    i=int((N-Nstart)/Ninc);
    
    double g0 = listgs(i);
    double g1 = listgs(i+1);
    return (g0+(g1-g0)*(N-Nstart-Ninc*i)/Ninc);
};

double V(double phi){
    return 0.5*std::pow(m,2)*std::pow(phi,nv);
};

double dV(double phi){
    return 0.5*std::pow(m,2)*nv*std::pow(phi,nv-1);
};

static void NAG_CALL f(double t, Integer n, const double *y, double *yp,
Nag_Comm *comm){
    yp[0] = y[1];
    yp[1] = -(3.0-(0.5*y[1]*y[1]))*y[1] - (6.0-(y[1]*y[1]))*dV(y[0])/(2.0*V(y[0]));
};

static void NAG_CALL ftot(double t, Integer n, const double *y, double *yp,
Nag_Comm *comm){
// phi, dphi, Re(Rk), Im(Rk), Re(dRk), Im(dRk)
    
    yp[0] = y[1];
    yp[1] = -(3.0-(0.5*y[1]*y[1]))*y[1] - (6.0-(y[1]*y[1]))*dV(y[0])/(2.0*V(y[0]));

    double H = std::pow(V(y[0])/(3.0-0.5*y[1]*y[1]),0.5); 
    double a = ai*std::exp(t);
    double z = a*y[1];
    double dz = a*(y[1]+yp[1]);
    double gt = -0.25*y[1]*y[1] + yp[1]/y[1] + 1.5; 
    double wt = k/(a*H);
    
    yp[2] = y[4]; 
    yp[3] = y[5]; 
    yp[4] = -wt*wt*y[2]-2.0*gt*y[4];
    yp[5] = -wt*wt*y[3]-2.0*gt*y[5];
};

int main(){

    // Range of wavenumbers
    int nk=2;
    Eigen::VectorXd ks=Eigen::VectorXd::LinSpaced(nk,-5,8);
    for(int i=0; i<nk; i++)
        ks(i) = std::pow(10,ks(i));
    Eigen::VectorXd exits=Eigen::VectorXd::Zero(nk);
    Eigen::VectorXd starts=Eigen::VectorXd::Zero(nk);
    Eigen::VectorXd startindices=Eigen::VectorXd::Zero(nk);
    
    // These will contain N,w,g
    Eigen::VectorXd Ns=Eigen::VectorXd::Zero(Nbg);
    //logws=Eigen::VectorXd::Zero(Nbg);
    //listgs=Eigen::VectorXd::Zero(Nbg);
    // To log background evolution
    Eigen::VectorXd Hs=Eigen::VectorXd::Zero(Nbg), phis=Eigen::VectorXd::Zero(Nbg),
    dphis=Eigen::VectorXd::Zero(Nbg), ddphis=Eigen::VectorXd::Zero(Nbg);
    std::cout << "Starting background calculation" << std::endl;
    // Use the NAG library to solve for phi(N), dphi(N).
    Integer liwsav, lrwsav, n=2;
    double tgot, tol, tnext, twant;
    double *rwsav=0, *thresh=0, *ygot=0, *yinit=0, *ymax=0;
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
    ypgot = NAG_ALLOC(n, double);
    ymax = NAG_ALLOC(n, double);
    iwsav = NAG_ALLOC(liwsav, Integer);
    rwsav = NAG_ALLOC(lrwsav, double);
    method = (Nag_RK_method) nag_enum_name_to_value("Nag_RK_7_8");
    errass = (Nag_ErrorAssess) nag_enum_name_to_value("Nag_ErrorAssess_off");
    
    // Define slow-roll solution
    double Vinit = V(phi_p), dVinit = dV(phi_p);
    double Hinit = std::sqrt(Vinit/3.0);
    double phidot_init = -dVinit/(3.0*Hinit);
    Hinit = std::sqrt(1.0/3.0*(0.5*phidot_init*phidot_init + Vinit));
    double dphi_init = phidot_init/Hinit;
    // Hubble slow-roll parameter, end of inflation, scale factor
    double Nendinfl;
    double epsilon_h=0.0;
    // Set initial values of phi, dphi
    yinit[0] = phi_p;
    yinit[1] = dphi_init;
    // Set zeroth values for background vectors: N, phi, dphi, ddphi, H
    Ns(0) = Nstart; 
    phis(0) = phi_p;
    dphis(0) = dphi_init;
    ddphis(0) = -(3.0-0.5*dphi_init*dphi_init)*dphi_init - dVinit/(Hinit*Hinit);
    Hs(0) = std::sqrt(Vinit/(3.0-0.5*dphi_init*dphi_init));
    // Solver settings 
    thresh[0] = 1.0e-14;
    thresh[1] = 1.0e-14;

    tol = 1.0e-14;
    nag_ode_ivp_rkts_setup(n, Nstart, Nend, yinit, tol, thresh, method, errass, 0.0, iwsav, rwsav, &fail);
    twant=Nstart;
    
    for(int i=1; i<=(npts); i++){
        tnext = twant;
        twant+=Ninc;
        while(tnext < twant){
            nag_ode_ivp_rkts_range(f,2,twant,&tgot,ygot,ypgot,ymax,&comm,iwsav,rwsav,&fail);
            tnext = tgot; 
            };
        Ns(i) = twant;
        phis(i) = ygot[0];
        dphis(i) = ygot[1];
        ddphis(i) = ypgot[1];
        Hs(i) = std::pow(V(ygot[0])/(3.0-0.5*ygot[1]*ygot[1]),0.5); 
        epsilon_h = 0.5*dphis(i)*dphis(i);
        if(epsilon_h>1.0){
            Nendinfl = twant;
            break;
        };
    };

    // Find the scale factor (for rescaling)
    double Nexit = Nendinfl-Npivot, Hexit;
    for(int i=0; i<npts; i++){
        if(Ns(i)>Nexit){
            Hexit = Hs(i);
            break;
        };
    };
    ai = kpivot/(Hexit*std::exp(Nexit));
    std::cout << "scale factor is: " << ai << std::endl;
        
    NAG_FREE(thresh);
    NAG_FREE(yinit);
    NAG_FREE(ygot);
    NAG_FREE(ypgot);
    NAG_FREE(ymax);
    NAG_FREE(rwsav);
    NAG_FREE(iwsav);

    // Finding beginning and end of integration for each mode
    double epsilon0, epsilon1;
    for(int i=0; i<nk; i++){
        for(int j=0; j<Nbg; j++){
            epsilon0 = ks(i)/(100.0*ai*std::exp(Ns(j))) - Hs(j);
            if(epsilon0 <= 0.0){
               starts(i) = Ns(j);
               startindices(i) = j;
               break;
            };
        };
    };

    for(int i=0; i<nk; i++){
        for(int j=0; j<Nbg; j++){
            epsilon1 = ks(i)/(1e-2*ai*std::exp(Ns(j))) - Hs(j);
            if(epsilon1 <= 0.0){
               exits(i) = Ns(j);
               break;
            };
        };
    };

    // Solve perturbation from sub-horizon to super-horizon scales with NAG
    // 1. Get continuous reference line
    int pert_npts = 5e3, startsi;
    double z, dz_z, a0, ti, tf;
    Eigen::VectorXd rks = Eigen::VectorXd::Zero(pert_npts), rksim = Eigen::VectorXd::Zero(pert_npts);
    Eigen::VectorXd ts = Eigen::VectorXd::Zero(pert_npts);
    double tinc;
    for(int i=0; i<nk; i++){
        k = ks(i);
        startsi = startindices(i);
        ti = starts(i);
        tf = exits(i); 
        a0 = ai*std::exp(Ns(startsi));
        z = a0*dphis(startsi);
        dz_z = ddphis(startsi)/dphis(startsi) + 1.0;
        tinc = (tf-ti)/pert_npts;

        Integer liwsav, lrwsav; 
        n=6;
        double *rwsav=0, *thresh=0, *ygot=0, *yinit=0, *ymax=0;
        double *ypgot=0;
        Integer *iwsav=0;
        liwsav = 130;
        lrwsav = 350 + 32 * n;
        thresh = NAG_ALLOC(n, double);
        ygot = NAG_ALLOC(n, double);
        yinit = NAG_ALLOC(n, double);
        ypgot = NAG_ALLOC(n, double);
        ymax = NAG_ALLOC(n, double);
        iwsav = NAG_ALLOC(liwsav, Integer);
        rwsav = NAG_ALLOC(lrwsav, double);
        method = (Nag_RK_method) nag_enum_name_to_value("Nag_RK_7_8");
        errass = (Nag_ErrorAssess) nag_enum_name_to_value("Nag_ErrorAssess_off");
        
        // Set initial values of phi, dphi
        yinit[0] = phis(startsi);
        yinit[1] = dphis(startsi);
        yinit[2] = 1.0/(std::sqrt(2.0*k))/z;
        yinit[3] = 0.0;
        yinit[4] = std::real(-yinit[2]*dz_z);
        yinit[5] = -std::sqrt(k/2.0)/(a0*z*Hs(startsi));
        std::cout << "initial conditions" << std::endl;
        std::cout << yinit[0] << std::endl;
        std::cout << yinit[1] << std::endl;
        std::cout << yinit[2] << std::endl;
        std::cout << yinit[4] << std::endl;
        std::cout << yinit[3] << std::endl;
        std::cout << yinit[5] << std::endl;

        // Set zeroth value of solution vector
        rks(0) = yinit[2];
        rksim(0) = yinit[3];
        ts(0) = ti;
        // Solver settings 
        thresh[0] = 1.0e-14;
        thresh[1] = 1.0e-14;
        thresh[2] = 1.0e-14;
        thresh[3] = 1.0e-14;
        thresh[4] = 1.0e-14;
        thresh[5] = 1.0e-14;

        tol = 1e-4;
        nag_ode_ivp_rkts_setup(n, ti, tf, yinit, tol, thresh, method, errass, 0.0, iwsav, rwsav, &fail);
        twant=ti;
        
        for(int i=1; i<(pert_npts); i++){
            tnext = twant;
            twant+=tinc;
            while(tnext < twant){
                nag_ode_ivp_rkts_range(ftot,n,twant,&tgot,ygot,ypgot,ymax,&comm,iwsav,rwsav,&fail);
                tnext = tgot; 
                };
            ts(i) = tgot;
            rks(i) = ygot[2]; 
            rksim(i) = ygot[3];
        };
        
        NAG_FREE(thresh);
        NAG_FREE(yinit);
        NAG_FREE(ygot);
        NAG_FREE(ypgot);
        NAG_FREE(ymax);
        NAG_FREE(rwsav);
        NAG_FREE(iwsav);

        std::ofstream f; std::string filename;
        std::cout << "name of outputfile: " << std::endl;
        std::cin >> filename;
        f.open(filename);
        for(int i=0; i<pert_npts; i++){
            f << ts(i) << " " << rks(i) << " " << rksim(i) << std::endl;  
        };
        f.close();
        
        // 2. Get individual RK steps:

        //Integer liwsav, lrwsav; 
        //n=6;
        //double *rwsav=0, *thresh=0, *ygot=0, *yinit=0, *ymax=0;
        //double *ypgot=0;
        //Integer *iwsav=0;
        liwsav = 130;
        lrwsav = 350 + 32 * n;
        thresh = NAG_ALLOC(n, double);
        ygot = NAG_ALLOC(n, double);
        yinit = NAG_ALLOC(n, double);
        ypgot = NAG_ALLOC(n, double);
        ymax = NAG_ALLOC(n, double);
        iwsav = NAG_ALLOC(liwsav, Integer);
        rwsav = NAG_ALLOC(lrwsav, double);
        method = (Nag_RK_method) nag_enum_name_to_value("Nag_RK_7_8");
        errass = (Nag_ErrorAssess) nag_enum_name_to_value("Nag_ErrorAssess_off");
        // Define solution vectors
        std::vector<double> rks2, ts2;

        // Set initial values of phi, dphi
        yinit[0] = phis(startsi);
        yinit[1] = dphis(startsi);
        yinit[2] = 1.0/(std::sqrt(2.0*k))/z;
        yinit[3] = 0.0;
        yinit[4] = std::real(-yinit[2]*dz_z);
        yinit[5] = -std::sqrt(k/2.0)/(a0*z*Hs(startsi));
        std::cout << "initial conditions" << std::endl;
        std::cout << yinit[0] << std::endl;
        std::cout << yinit[1] << std::endl;
        std::cout << yinit[2] << std::endl;
        std::cout << yinit[4] << std::endl;
        std::cout << yinit[3] << std::endl;
        std::cout << yinit[5] << std::endl;

        // Set zeroth value of solution vector
        rks2.push_back(yinit[2]);
        ts2.push_back(ti);
        // Solver settings 
        thresh[0] = 1.0e-14;
        thresh[1] = 1.0e-14;
        thresh[2] = 1.0e-14;
        thresh[3] = 1.0e-14;
        thresh[4] = 1.0e-14;
        thresh[5] = 1.0e-14;
        tol = 1e-4;
        nag_ode_ivp_rkts_setup(n, ti, tf, yinit, tol, thresh, method, errass, 0.0, iwsav, rwsav, &fail);
        
        tgot = ti;
        while(tgot < tf){ 
            nag_ode_ivp_rkts_onestep(ftot,n,&tgot,ygot,ypgot,&comm,iwsav,rwsav,&fail);
            ts2.push_back(tgot);
            rks2.push_back(ygot[2]); 
        };

        std::ofstream f2; std::string filename2;
        std::cout << "name of outputfile with individual RK steps: " << std::endl;
        std::cin >> filename2;
        f2.open(filename2);
        auto it1 = ts2.begin(), it2 = rks2.begin();
        while(it1!=ts2.end()){
            f2 << *it1 << " " << *it2 << std::endl;  
            ++it1;
            ++it2;
        };
        f2.close();
    };




//    // Construct list of logws, gs
//    for(int i=0; i<Nbg; i++){
//        logws(i) = -std::log(ai*Hs(i)) - Ns(i);
//        listgs(i) = -0.25*dphis(i)*dphis(i) + ddphis(i)/dphis(i) + 1.5;
//    };
//
//    // Construct system to solve
//    de_system system(&w,&g);
//    
//    // Solve the evolution of each perturbation
//    double ti, tf, rtol, atol, h0;
//    std::complex<double> x0, dx0;
//    int order=3;
//    bool full_output=false;
//    rtol=1e-4;
//    atol=0.0;
//    h0=1.0;
//    
//    std::list<std::complex<double>> rk1, rk2;
//    std::list<double> times;
//    double t1,t2,timing0,timing1,startsi;   
//    double z, dz_z, a0;
//    
//    timing0 = MPI_Wtime();
//    CALLGRIND_START_INSTRUMENTATION;
//    CALLGRIND_TOGGLE_COLLECT;
//    for(int i=0; i<nk; i++){
//        k = ks(i);
//        startsi = startindices(i);
//        ti = starts(i);
//        tf = exits(i); 
//        a0 = ai*std::exp(Ns(startsi));
//        z = a0*dphis(startsi);
//        dz_z = ddphis(startsi)/dphis(startsi) + 1.0;
//        x0 = 1.0/(std::sqrt(2.0*k))/z;
//        dx0 = std::complex<double>(std::real(-x0*dz_z),-std::sqrt(k/2.0)/(a0*z*Hs(startsi)));
//       
//        Solution solution1(system, x0, dx0, ti, tf, order, rtol, atol, h0, full_output);
//        t1 = MPI_Wtime();
//        solution1.solve();
//        rk1.push_back(solution1.sol.back());
//        //x0 = 0.0;
//        //dx0 = 1e4;
//        //Solution solution2(system, x0, dx0, ti, tf, order, rtol, atol, h0, full_output);
//        //solution2.solve();
//        t2 = MPI_Wtime();
//        //rk2.push_back(solution2.sol.back());
//        times.push_back(t2-t1);
//        
//    };
//    CALLGRIND_TOGGLE_COLLECT;
//    CALLGRIND_STOP_INSTRUMENTATION;
//    timing1 = MPI_Wtime();
//    std::cout << "total time: " << timing1-timing0 << " s." << std::endl;
//    
//    auto it_1 = rk1.begin();
//    auto it_2 = rk2.begin();
//    auto it_3 = times.begin();
//    while(it_1!=rk1.end() && it_2!=rk2.end() && it_3!=times.end()){
//        //std::cout << "rk1: " << *it_1 << ", rk2: " << *it_2 << std::endl;
//        ++it_1;
//        ++it_2;
//        ++it_3;
//    };
//
//    // Write PPS to file
//    std::ofstream f;
//    f.open("test/ms/pps-bingo-start2e2-highk-lowk.txt");
//    it_1 = rk1.begin();
//    it_2 = rk2.begin();
//    it_3 = times.begin();
//    double starti;
//    for(int k=0; k<nk; k++){
//        starti = startindices(k);
//        f << ks(k) << ", "  << *it_1 << ", " << *it_2 << ", " <<
////        hd(ks(k),phis(starti),dphis(starti),ddphis(starti),Hs(starti),Ns(starti),ai,*it_1,*it_2) << ", " << *it_3 << std::endl;
//        std::abs(*it_1)*std::abs(*it_1)*std::pow(ks(k),3)/(2.0*M_PI*M_PI) << ", " << *it_3 << std::endl;
//        ++it_1;
//        ++it_2;
//        ++it_3;
//    };
//    f.close();

    return 0;
}
