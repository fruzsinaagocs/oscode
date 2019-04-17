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
double m=5e-6;
double k; // wavevector of perturbation
double K=1; // curvature
double mp=1;
double kpivot=0.05;
double Npivot=54.0;
int Nbg=round(2e6), npts=round(Nbg-1);
double Nstart=0.0, Nend=75.0, Ninc=(Nend-Nstart)/npts;
Eigen::VectorXd logws, listgs;
double Ak=0.0, Bk=1.0;
// Initial conditions at horizon turnover
double Ni = 10.16;
double o_ki = 1e2;
// For constructing w,g
Eigen::VectorXd dE_E=Eigen::VectorXd::Zero(Nbg),
E=Eigen::VectorXd::Zero(Nbg), o_ks=Eigen::VectorXd::Zero(Nbg),
logo_ks=Eigen::VectorXd::Zero(Nbg);


// For setting Kinetically Dominated initial conditions for the perturbation
// modes

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

Eigen::Vector2cd kd(double ki, double phi, double dphi, double ddphi, double N){
    // Initial conditions for the perturbations in Kinetic Dominance at t0
    Eigen::Vector2cd result;
    double z = std::exp(N)*dphi;
    double dz_z = ddphi/dphi + 1.0;
    std::complex<double> R =
    1.0/(2.0*z)*std::sqrt(1.5*M_PI)*std::exp(N)*(Ak*H1_0(1.5*ki*std::exp(2.0*N))
    + Bk*H2_0(1.5*ki*std::exp(2.0*N)));
    std::complex<double> dR = (-dz_z+1.0)*R +
    3.0*ki/(2.0*z)*std::sqrt(1.5*M_PI)*std::exp(3.0*N)*(Ak*H1_0_prime(1.5*ki*std::exp(2.0*N))
    + Bk*H2_0_prime(1.5*ki*std::exp(2.0*N)));
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
    double power = std::pow(std::abs(a*rk1 + b*rk2),2)*std::pow(ki,3)/(2.0*M_PI*M_PI);
    return power;

};

//std::sqrt(o_ks(i)*(ksq - K -2.0*K*ksq*dE_E(i)/(E(i)*K+ksq)))

double RKSolver::w(double N){
    int i;
    i=int((N-Nstart)/Ninc);
    double k2 = k*(k+2.0)-3.0*K;
    double logw0 = 0.5*(logo_ks(i) +  std::log(k2 - K -2.0*K*k2*dE_E(i)/(E(i)*K+k2)));
    double logw1 = 0.5*(logo_ks(i) + std::log(k2 - K -2.0*K*k2*dE_E(i+1)/(E(i+1)*K+k2)));
    return std::exp(logw0+(logw1-logw0)*(N-Nstart-Ninc*i)/Ninc);
};

double RKSolver::g(double N){
    int i;
    i=int((N-Nstart)/Ninc);
    double k2 = k*(k+2.0)-3.0*K; 
    double g0 = 0.5*(K*o_ks(i) + 3.0 - E(i) + dE_E(i)*k2/(E(i)*K+k2));
    double g1 = 0.5*(K*o_ks(i+1) + 3.0 - E(i+1) + dE_E(i+1)*k2/(E(i+1)*K+k2));
    return (g0+(g1-g0)*(N-Nstart-Ninc*i)/Ninc);
};

// Dummy frequency after having moved away from std::function, TODO:remove and
// move back to std::function
double win(double){
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
    // equation giving dy with y = [logomega_k, phi/mp], 
    // dy = [ dlogomega_k, dphi/mp]
    //std::cout << "N: " << t << ", dlogok: " << yp[0] << std::endl;
    yp[0] = 4.0 + std::exp(y[0])*(4.0*K - 2.0*std::exp(2.0*t)*V(y[1]*mp)/(mp*mp));
    yp[1] = -std::sqrt(6.0 + std::exp(y[0])*(6.0*K -
    2.0*std::exp(2.0*t)*V(y[1]*mp)/(mp*mp)));
};

int main(){

    // Range of wavenumbers
    int nk=1000;
    Eigen::VectorXd ks=Eigen::VectorXd::LinSpaced(nk,-4,1);
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
    Eigen::VectorXd phis=Eigen::VectorXd::Zero(Nbg),
    dphis=Eigen::VectorXd::Zero(Nbg), 
    dlogo_ks=Eigen::VectorXd::Zero(Nbg);
    std::cout << "Starting background calculation" << std::endl;
    // Set up the NAG library to solve for the background
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
   
    // 1. Integrate forward towards SR

    // Determine closest point in grid to set i.c.
    int mid = int(npts*(Ni-Nstart)/(Nend-Nstart));
    Ni = Nstart + Ninc*mid;
    // Define i.c. at top of Hubble horizon
    double phi_i = std::sqrt(4*(1.0/o_ki + K)*std::exp(-2.0*Ni)/(m*m));
    std::cout << "middle phi: " << phi_i << std::endl;
    // Hubble slow-roll parameter, end of inflation, scale factor
    double Nendinfl;
    double epsilon_h=0.0, epsilon_prev=1.2;
    double ai=1.0;
    // Set initial values of phi, dphi
    yinit[0] = std::log(o_ki);
    yinit[1] = phi_i;
    std::cout << "i.c.: " << yinit[0] << ", " << yinit[1] << std::endl;
    // Set middle values for background vectors: N, phi, dphi, ddphi, H
    Ns(mid) = Ni;
    phis(mid) = phi_i;
    dphis(mid) = -std::sqrt(6.0 + o_ki*(6.0*K -
    2.0*std::exp(2.0*Ni)*V(phi_i*mp)/(mp*mp)));
    logo_ks(mid) = std::log(o_ki);
    dlogo_ks(mid) = 4.0 + o_ki*(4.0*K - 2.0*std::exp(2.0*Ni)*V(phi_i*mp)/(mp*mp));
    // Solver settings 
    thresh[0] = 1.0e-14;
    thresh[1] = 1.0e-14;
    tol = 1.0e-14;
    nag_ode_ivp_rkts_setup(n, Ni, Nend, yinit, tol, thresh, method, errass, 0.0, iwsav, rwsav, &fail);
    twant=Ni;
    
    for(int i=mid+1; i<=(npts); i++){
        tnext = twant;
        twant+=Ninc;
        //std::cout << tnext << std::endl;
        while(tnext < twant){
            nag_ode_ivp_rkts_range(f,2,twant,&tgot,ygot,ypgot,ymax,&comm,iwsav,rwsav,&fail);
            //if(fail.code != NE_NOERROR){
            //    printf("error from NAG.\n%s\n",fail.message);
            //}
            tnext = tgot; 
            };
        Ns(i) = twant;
        phis(i) = ygot[1];
        dphis(i) = ypgot[1];
        logo_ks(i) = ygot[0];
        dlogo_ks(i) = ypgot[0];
        epsilon_h = 0.5*dphis(i)*dphis(i)*mp*mp;
        if(epsilon_h>1.0 and epsilon_prev<1.0){
            Nendinfl = twant;
            std::cout << "inflation ends at N=" << Nendinfl << std::endl;
            break;
        };
        epsilon_prev = epsilon_h;
    };
       
    NAG_FREE(thresh);
    NAG_FREE(yinit);
    NAG_FREE(ygot);
    NAG_FREE(ypgot);
    NAG_FREE(ymax);
    NAG_FREE(rwsav);
    NAG_FREE(iwsav);
    
    // 2. Integrate backwards towards KD
    
    thresh = NAG_ALLOC(n, double);
    ygot = NAG_ALLOC(n, double);
    yinit = NAG_ALLOC(n, double);
    ypgot = NAG_ALLOC(n, double);
    ymax = NAG_ALLOC(n, double);
    iwsav = NAG_ALLOC(liwsav, Integer);
    rwsav = NAG_ALLOC(lrwsav, double);
    yinit[0] = std::log(o_ki);
    yinit[1] = phi_i;
    // Solver settings 
    thresh[0] = 1.0e-14;
    thresh[1] = 1.0e-14;
    tol = 1.0e-14;

    std::cout << "starting backwards integration" << std::endl;
    nag_ode_ivp_rkts_setup(n, Ni, Nstart, yinit, tol, thresh, method, errass, 0.0, iwsav, rwsav, &fail);
    twant=Ni;
    
    for(int i=mid-1; i>=0; i--){
        tnext = twant;
        twant-=Ninc;
        if(twant < Nstart)
            twant = Nstart;
        //std::cout << twant << std::endl;
        while(tnext > twant){
            nag_ode_ivp_rkts_range(f,2,twant,&tgot,ygot,ypgot,ymax,&comm,iwsav,rwsav,&fail);
            tnext = tgot; 
            };
        Ns(i) = twant;
        phis(i) = ygot[1];
        dphis(i) = ypgot[1];
        logo_ks(i) = ygot[0];
        dlogo_ks(i) = ypgot[0];
    };
    
    NAG_FREE(thresh);
    NAG_FREE(yinit);
    NAG_FREE(ygot);
    NAG_FREE(ypgot);
    NAG_FREE(ymax);
    NAG_FREE(rwsav);
    NAG_FREE(iwsav);

    // Construct variables used in logws, gs
    for(int i=0; i<Nbg; i++){
       dE_E(i) =
       dlogo_ks(i)-4.0-2.0*dV(phis(i)*mp)*std::exp(logo_ks(i))*std::exp(2.0*Ns(i))/dphis(i)/mp;
       E(i) = 0.5*dphis(i)*dphis(i)*mp*mp;
       o_ks(i) = std::exp(logo_ks(i));
    };

    double atoday = 4.3e4; 
    for(int i=0; i<nk; i++)
        ks(i) = int(ks(i)*atoday);
    std::cout << "new scales: " << ks(0) << ".." << ks(nk-1) << std::endl;
 
    // Write background to file
    std::ofstream fbg;
    fbg.open("test/ms/kd-closed-bg-tip.txt");
    double ksq = 4000*(4000+2.0)-3*K;
    for(int i=0; i<Nbg; i++){
        fbg << Ns(i) << ", " << phis(i) << ", " << dphis(i) << ", " <<
        logo_ks(i) << ", " << dlogo_ks(i) << ", " << std::sqrt(o_ks(i)*(ksq - K
        -2.0*K*ksq*dE_E(i)/(E(i)*K+ksq))) << ", " << 0.5*(K*o_ks(i) + 3.0 - E(i) +
        dE_E(i)*ksq/(E(i)*K+ksq)) << ", " << dE_E(i) << std::endl;
    };
    fbg.close();
    std::cout << "Done solving background" << std::endl;

    std::string break_msg;
    std::cout << "press anything to continue" << std::endl;
    std::cin >> break_msg;

    // Finding beginning and end of integration for each mode
    for(int i=0; i<nk; i++){
        for(int j=0; j<Nbg; j++){
            if(Ns(j) > Ni){
               starts(i) = Ns(j);
               startindices(i) = j;
               break;
            };
        };
    };

    double epsilon1;
    for(int i=0; i<nk; i++){
        for(int j=startindices(i); j<Nbg; j++){
            epsilon1 = 1e2*ks(i)*std::sqrt(o_ks(j)) - 1.0;
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

    // Construct system to solve
    de_system system(&win,&win);
    
    // Solve the evolution of each perturbation
    double ti, tf, rtol, atol, h0;
    std::complex<double> x0, dx0;
    int order=3;
    bool full_output=false;//true;
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
        phi = phis(startsi)*mp;
        dphi = dphis(startsi)*mp;
        ddphi = 0.5*dE_E(startsi)*dphi*mp;
        N = Ns(startsi);
        ic = kd(k,phi,dphi,ddphi,N);
        x0s(i) = ic(0);
        dx0s(i) = ic(1);
        a0 = std::exp(N);
        z = a0*dphi;
        dz_z = ddphi/dphi + 1.0;
        //std::cout << "starting vars. phi: " << phi << ", dphi: " << dphi << ", ddphi: " << ddphi << ", N: " << N << ", ic: " << ic << ", z: " << z << ", dz_z: " << dz_z << std::endl;
        //std::cout << "dz_z: " << dz_z << ", ddphi: " << ddphi << ", dphi: " << dphi << std::endl;
        x0 = 1.0/(std::sqrt(2.0*k))/z;
        dx0 = std::complex<double>(std::real(-x0*dz_z),-std::sqrt(k/2.0)*std::sqrt(o_ks(startsi))/z);
        x01s(i) = x0;//1.0/(std::sqrt(2.0*k))/z;
        dx01s(i) = dx0;//std::complex<double>(std::real(-x01s(i)*dz_z),-std::sqrt(k/2.0)*std::sqrt(o_ks(startsi))/z);
        Solution solution1(system, x0, dx0, ti, tf, order, rtol, atol, h0,
        full_output);
        t1 = MPI_Wtime();
        solution1.solve();
        rk1.push_back(solution1.sol.back());
        
        x0 = 1.0/(std::sqrt(2.0*k))/z;
        dx0 = -std::complex<double>(std::real(-x0*dz_z),-std::sqrt(k/2.0)*std::sqrt(o_ks(startsi))/z);
        x02s(i) = x0;//1.0/(std::sqrt(2.0*k))/z;
        dx02s(i) = dx0;//-std::complex<double>(std::real(-x02s(i)*dz_z),-std::sqrt(k/2.0)*std::sqrt(o_ks(startsi))/z);
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
    f.open("test/ms/pps-kd-closed-test7.txt");
    it_1 = rk1.begin();
    it_2 = rk2.begin();
    it_3 = times.begin();
    double starti;
    for(int u=0; u<nk; u++){
        starti = startindices(u);
        f << ks(u)/atoday << ", "  << *it_1 << ", " << *it_2 << ", " <<
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
