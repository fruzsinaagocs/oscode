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
double phi_p=23.7;
double mp=1;
double kpivot=0.05;
double Npivot=50.0;
int Nbg=round(5e5), npts=round(Nbg-1);
double Nstart=0.0, Nend=68.0, Ninc=(Nend-Nstart)/npts;
Eigen::VectorXd logws, listgs;
double Ak=0.0, Bk=1.0;
double ai;

std::complex<double> H1_0(double x){
    return std::complex<double>(boost::math::cyl_bessel_j(0,x),boost::math::cyl_neumann(0,x));
};

std::complex<double> H2_0(double x){
    return std::complex<double>(boost::math::cyl_bessel_j(0,x),-boost::math::cyl_neumann(0,x));
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
    double z = ai*std::exp(N)*dphi;
    double dz_z = ddphi/dphi + 1.0;
    std::complex<double> R =
    1.0/(2.0*z)*std::sqrt(1.5*M_PI)*std::exp(N)*(Ak*H1_0(1.5*ki*std::exp(2.0*N))
    + Bk*H2_0(1.5*ki*std::exp(2.0*N)));
    std::complex<double> dR = (-dz_z+1.0)*R +
    3.0*ki/(2.0*z)*std::sqrt(1.5*M_PI)*std::exp(N)*(Ak*H1_0_prime(1.5*ki*std::exp(2.0*N))
    + Bk*H2_0_prime(1.5*ki*std::exp(2.0*N)));
    result << R, dR; 
    return result;
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
    int nk=1500;
    Eigen::VectorXd ks=Eigen::VectorXd::LinSpaced(nk,-5,3);
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
    
    // Define KD solution
    double Vinit = V(phi_p), dVinit = dV(phi_p);
    double Hinit = 1.0/3.0*std::exp(-3.0*Nstart);
    double dphi_init = -std::sqrt(6.0-2.0*Vinit/(Hinit*Hinit));
    // Hubble slow-roll parameter, end of inflation, scale factor
    double Nendinfl;
    double epsilon_h=0.0, epsilon_prev=1.2;
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
        if(epsilon_h>1.0 and epsilon_prev<1.0){
            Nendinfl = twant;
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
    for(int i=0; i<nk; i++){
        for(int j=0; j<Nbg; j++){
            if(Ns(j) > 1.0){
               starts(i) = Ns(j);
               startindices(i) = j;
               break;
            };
        };
    };

    double epsilon1;
    for(int i=0; i<nk; i++){
        for(int j=0; j<Nbg; j++){
            epsilon1 = ks(i)/(1e-2*ai*std::exp(Ns(j))) - Hs(j);
            if(epsilon1 <= 0.0){
               exits(i) = Ns(j);
               break;
            };
        };
    };

    std::cout << "starts: " << starts(0) << ", " << starts(nk-1) << std::endl;
    std::cout << "start indices: " << startindices(0) << ", " <<
    startindices(nk-1) << std::endl;
    std::cout << "exits: " << exits(0) << ", " << exits(nk-1) << std::endl;


    // Solve perturbation from sub-horizon to super-horizon scales with NAG
    // 1. Get continuous reference line
    int pert_npts = 5e3, startsi;
    double ti,tf;
    double phi,dphi,ddphi,N,z,dz_z,a0;
    Eigen::Vector2cd ic;
    Eigen::VectorXd rks = Eigen::VectorXd::Zero(pert_npts), rksim = Eigen::VectorXd::Zero(pert_npts);
    Eigen::VectorXd ts = Eigen::VectorXd::Zero(pert_npts);
    double tinc;
    std::complex<double> x0,dx0;
    for(int i=0; i<2; i++){
        k = ks(i);
        startsi = startindices(i);
        ti = starts(i);
        tf = exits(i); 
        phi = phis(startsi);
        dphi = dphis(startsi);
        ddphi = ddphis(startsi);
        N = Ns(startsi);
        a0 = ai*std::exp(N);
        z = a0*dphi;
        dz_z = ddphi/dphi + 1.0;
        ic = kd(k,phi,dphi,ddphi,N,ai);
        x0 = 1.0/(std::sqrt(2.0*k))/z;
        dx0 = std::complex<double>(std::real(-x0*dz_z),-std::sqrt(k/2.0)/(a0*z*Hs(startsi)));
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
        yinit[0] = phi;
        yinit[1] = dphi;
        yinit[2] = std::real(x0);
        yinit[3] = std::imag(x0);
        yinit[4] = std::real(dx0);
        yinit[5] = std::imag(dx0);
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
            std::cout << "t: " << ts(i) << ", rk: " << rks(i) << std::endl;
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
        
    };

    return 0;
}
