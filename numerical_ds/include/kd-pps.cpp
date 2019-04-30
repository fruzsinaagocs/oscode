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
#include <boost/math/special_functions/hankel.hpp>
#include "interpolator.hpp"
#include "system.hpp"
#include "solver.hpp"
#include "mpi.h"

int nv=2;
double m=5.0e-6;
double k;
double phi_p=23.1;
double mp=1;
double kpivot=0.05;
double Npivot=54.0;
int Nbg=round(2e6), npts=round(Nbg-1);
double Nstart=0.0, Nend=68.0, Ninc=(Nend-Nstart)/npts;
Eigen::VectorXd logws, listgs;
double Ak=0.0, Bk=1.0;

// For setting Kinetically Dominated initial conditions

std::complex<double> H1_0(double x){
    return boost::math::cyl_hankel_1(0,x);
};

std::complex<double> H2_0(double x){
    return boost::math::cyl_hankel_2(0,x);
};

std::complex<double> H1_0_prime(double x){
    return std::complex<double>(boost::math::cyl_bessel_j_prime(0,x),boost::math::cyl_neumann_prime(0,x));
};

std::complex<double> H2_0_prime(double x){
    return std::complex<double>(boost::math::cyl_bessel_j_prime(0,x),-boost::math::cyl_neumann_prime(0,x));
};

Eigen::Vector2cd kd(double ki, double phi, double dphi, double ddphi, double N, double ai){
    // Initial conditions for the perturbations in Kinetic Dominance at t0
    Eigen::Vector2cd result;
    N = std::log(ai) + N;
    double z = std::exp(N)*dphi;
    double dz_z = ddphi/dphi + 1.0;
    std::complex<double> R =
    1.0/(2.0*z)*std::sqrt(1.5*M_PI)*std::exp(N)*(Ak*H1_0(std::pow(ai,-3)*1.5*ki*std::exp(2.0*N))
    + Bk*H2_0(std::pow(ai,-3)*1.5*ki*std::exp(2.0*N)));
    std::complex<double> dR = (-dz_z+1.0)*R +
    3.0*ki/(2.0*z)*std::sqrt(1.5*M_PI)*std::exp(3.0*N)*std::pow(ai,-3)*(Ak*H1_0_prime(std::pow(ai,-3)*1.5*ki*std::exp(2.0*N))
    + Bk*H2_0_prime(1.5*ki*std::exp(std::pow(ai,-3)*2.0*N)));
    result << R, dR; 
    return result;
};

double pps(double ki, std::complex<double> rk1, std::complex<double> rk2,
std::complex<double> x01, std::complex<double> dx01, std::complex<double> x02,
std::complex<double> dx02, std::complex<double> x0, std::complex<double> dx0){ 
    
    // Find linear combination of arbitrary, independent i.c. that give the
    // required i.c.

    std::complex<double> a = (x0*dx02 - dx0*x02)/(x01*dx02 - dx01*x02);
    std::complex<double> b = (x0*dx01 - dx0*x01)/(x02*dx01 - dx02*x01);
    //std::cout << "a: "<< a << ", b: " << b << std::endl;
    //std::cout << "rk1: " << rk1 << ", rk2: " << rk2 << std::endl;
    double power = std::pow(std::abs(a*rk1 + b*rk2),2)*std::pow(ki,3)/(2.0*M_PI*M_PI);
    return power;

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

// Dummy frequency after having moved away from std::function, TODO:remove and
// move back to std::function
std::complex<double> win(double){
    return 0.0;
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

int main(){

    // Range of wavenumbers
    int nk=2000;
    Eigen::VectorXd ks=Eigen::VectorXd::LinSpaced(nk,-5,1);
    for(int i=0; i<nk; i++)
        ks(i) = std::pow(10,ks(i));
    Eigen::VectorXd exits=Eigen::VectorXd::Zero(nk);
    Eigen::VectorXd starts=Eigen::VectorXd::Zero(nk);
    Eigen::VectorXd startindices=Eigen::VectorXd::Zero(nk);
    
    // These will contain N,w,g
    Eigen::VectorXd Ns=Eigen::VectorXd::Zero(Nbg);
    logws=Eigen::VectorXd::Zero(Nbg);
    listgs=Eigen::VectorXd::Zero(Nbg);
    // To log background evolution
    Eigen::VectorXd Hs=Eigen::VectorXd::Zero(Nbg), phis=Eigen::VectorXd::Zero(Nbg),
    dphis=Eigen::VectorXd::Zero(Nbg), ddphis=Eigen::VectorXd::Zero(Nbg);
    std::cout << "Starting background calculation" << std::endl;
    // Use the NAG library to solve for phi(N), dphi(N), H(N).
    Integer liwsav, lrwsav, n=2;
    double tgot, tol, tnext, twant;
    double *rwsav=0, *thresh=0, *ygot=0, *yinit=0, *ybg=0, *dybg=0, *ymax=0;
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
    ybg = NAG_ALLOC(n, double);
    dybg = NAG_ALLOC(n, double);
    ypgot = NAG_ALLOC(n, double);
    ymax = NAG_ALLOC(n, double);
    iwsav = NAG_ALLOC(liwsav, Integer);
    rwsav = NAG_ALLOC(lrwsav, double);
    method = (Nag_RK_method) nag_enum_name_to_value("Nag_RK_7_8");
    errass = (Nag_ErrorAssess) nag_enum_name_to_value("Nag_ErrorAssess_off");
    
    // Define KD solution
    double Vinit = V(phi_p), dVinit = dV(phi_p);
    double Hinit = 1.0/3.0*std::exp(-3.0*Nstart);
    double dphi_init = -std::sqrt(6.0-2.0*Vinit/(Hinit*Hinit));
    // Hubble slow-roll parameter, end of inflation, scale factor
    double Nendinfl;
    double epsilon_h=0.0, epsilon_prev=1.2;
    double ai;
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
        Hs(i) = std::pow(V(ygot[0])/(3.0-0.5*ygot[1]*ygot[1]),0.5); 
        ddphis(i) = -(3.0-0.5*dphis(i)*dphis(i))*dphis(i) - dVinit/(Hs(i)*Hs(i));;
        epsilon_h = 0.5*dphis(i)*dphis(i);
        if(epsilon_h>1.0 and epsilon_prev<1.0){
            Nendinfl = twant;
            std::cout << "inflation ends at N=" << Nendinfl << std::endl;
            break;
        };
        epsilon_prev = epsilon_h;
    };

    // Find the scale factor (for rescaling)
    double Nexit = Nendinfl-Npivot, Hexit;
    for(int i=0; i<npts; i++){
        if(Ns(i)>Nexit){
            Hexit = Hs(i);
            break;
        };
    };

    kpivot = Hexit*std::exp(Nexit);
    std::cout << "kpivot: " << kpivot << std::endl;
    for(int i=0; i<nk; i++)
        ks(i) = ks(i)*kpivot/0.05;
    //std::cout << "new scales: " << ks(0) << ".." << ks(nk-1) << std::endl;
    ai = 1.0;
    //std::cout << "scale factor is: " << ai << std::endl;
        
    NAG_FREE(thresh);
    NAG_FREE(yinit);
    NAG_FREE(ygot);
    NAG_FREE(ypgot);
    NAG_FREE(ymax);
    NAG_FREE(rwsav);
    NAG_FREE(iwsav);

    // Write phi(N) to file:
    std::ofstream fbg;
    fbg.open("test/ms/kd-bg.txt");
    for(int i=0; i<Nbg; i++){
        fbg << Ns(i) << ", " << phis(i) << ", " << dphis(i) << ", " << ddphis(i) << ", " << Hs(i) << std::endl;
    };
    fbg.close();
    std::cout << "Done solving background" << std::endl;

    // Finding beginning and end of integration for each mode
    for(int i=0; i<nk; i++){
        for(int j=0; j<Nbg; j++){
            if(Ns(j) > 1.1){
               starts(i) = Ns(j);
               startindices(i) = j;
               break;
            };
        };
    };

    double epsilon1;
    for(int i=0; i<nk; i++){
        for(int j=startindices(i); j<Nbg; j++){
            epsilon1 = ks(i)/(1e-2*ai*std::exp(Ns(j))) - Hs(j);
            if(epsilon1 <= 0.0){
               exits(i) = Ns(j);
               break;
            };
        };
    };
    
    //std::cout << "starts: " << starts(0) << ", " << starts(nk-1) << std::endl;
    //std::cout << "start indices: " << startindices(0) << ", " << startindices(nk-1) << std::endl;
    //std::cout << "exits: " << exits(0) << ", " << exits(nk-1) << std::endl;

    // Construct list of logws, gs
    for(int i=0; i<Nbg; i++){
        logws(i) = -std::log(ai*Hs(i)) - Ns(i);
        listgs(i) = -0.25*dphis(i)*dphis(i) + -(3.0-(0.5*dphis(i)*dphis(i))) - (6.0-(dphis(i)*dphis(i)))*dV(phis(i))/(2.0*V(phis(i))*dphis(i)) + 1.5;
    };

    // Construct system to solve
    de_system system(&w,&g);
    
    // Solve the evolution of each perturbation
    double ti, tf, rtol, atol, h0;
    std::complex<double> x0, dx0;
    int order=3;
    bool full_output=false;
    rtol=1e-4;
    atol=0.0;
    h0=1.0;
    
    std::list<std::complex<double>> rk1, rk2;
    std::list<double> times;
    double t1,t2,timing0,timing1,startsi;   
    double phi,dphi,ddphi,N,z,dz_z,a0;
    Eigen::Vector2cd ic;
    Eigen::VectorXcd x01s(nk), x02s(nk), dx01s(nk), dx02s(nk), x0s(nk), dx0s(nk);
    
    timing0 = MPI_Wtime();
    CALLGRIND_START_INSTRUMENTATION;
    CALLGRIND_TOGGLE_COLLECT;
    for(int i=0; i<nk; i++){
        k = ks(i);
        startsi = startindices(i);
        ti = starts(i);
        tf = exits(i); 
        phi = phis(startsi);
        dphi = dphis(startsi);
        ddphi = ddphis(startsi);
        N = Ns(startsi);
        ic = kd(k,phi,dphi,ddphi,N,ai);
        x0s(i) = ic(0);
        dx0s(i) = ic(1);
        a0 = ai*std::exp(N);
        z = a0*dphi;
        dz_z = ddphi/dphi + 1.0;
        x0 = 1.0/(std::sqrt(2.0*k))/z;
        dx0 =
        std::complex<double>(std::real(-x0*dz_z),-std::sqrt(k/2.0)/(a0*z*Hs(startsi)));
        //std::cout << "ic: " << x0 << "\\n" << dx0 << std::endl;
        x01s(i) = x0;
        dx01s(i) = dx0;
        Solution solution1(system, x0, dx0, ti, tf, order, rtol, atol, h0,
        full_output);
        t1 = MPI_Wtime();
        solution1.solve();
        rk1.push_back(solution1.sol.back());
        
        x0 = 1.0/(std::sqrt(2.0*k))/z;
        dx0 =
        -std::complex<double>(std::real(-x0*dz_z),-std::sqrt(k/2.0)/(a0*z*Hs(startsi)));
        x02s(i) = x0;
        dx02s(i) = dx0;
        Solution solution2(system, x0, dx0, ti, tf, order, rtol, atol, h0, full_output);
        solution2.solve();
        t2 = MPI_Wtime();
        rk2.push_back(solution2.sol.back());
        times.push_back(t2-t1);
        std::cout << "done k=" << k << " in " << t2-t1 << " s." << std::endl;
        
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
    f.open("test/ms/pps-testcomplexfn.txt");
    it_1 = rk1.begin();
    it_2 = rk2.begin();
    it_3 = times.begin();
    double starti;
    for(int u=0; u<nk; u++){
        starti = startindices(u);
        f << ks(u)*0.05/kpivot << ", "  << *it_1 << ", " << *it_2 << ", " <<
        pps(ks(u),*it_1,*it_2,x01s(u),dx01s(u),x02s(u),dx02s(u),x01s(u),dx02s(u))
        << ", " <<
        pps(ks(u),*it_1,*it_2,x01s(u),dx01s(u),x02s(u),dx02s(u),x0s(u),dx0s(u))
        << ", " << *it_3 << std::endl;
        ++it_1;
        ++it_2;
        ++it_3;
    };
    f.close();

    return 0;
}
