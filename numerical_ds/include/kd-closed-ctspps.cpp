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
double m=1;
//double m=5e-6;
double k; // wavevector of perturbation
double K=1; // curvature
double mp=1;
double kpivot=0.05;
double Npivot=54.0;
int Nbg=round(1e5), npts=round(Nbg-1);
double Nstart=-3.0, Nend=85.0, Ninc=(Nend-Nstart)/npts;
Eigen::VectorXd logws, listgs;
double Ak=0.0, Bk=1.0;
// Initial conditions at horizon turnover
double Ni = -2.08;//10.16;
double o_ki = 1e1; 
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
    // Initial conditions for the perturbations in Kinetic Dominance
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

Eigen::Vector2cd bd(double ki, double phi, double dphi, double ddphi, double N){
    // Initial conditions for the perturbations using Bunch--Davies vacuum
    Eigen::Vector2cd result;
    double dz_z = ddphi/dphi + 1.0;
    double a0 = std::exp(N);
    double z = a0*dphi;
    std::complex<double> R = 1.0/(std::sqrt(2.0*ki))/z;
    std::complex<double> dR =
    -std::complex<double>(-std::real(R*dz_z),-std::sqrt(ki/2.0*o_ki)/z);
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


std::complex<double> w(double N){
    int i;
    i=int((N-Nstart)/Ninc);
    double k2;
    if(K>0)
        k2 = k*(k+2.0)-3.0*K;
    else
        k2 = k*k-3.0*K;
    std::complex<double> logw0 = 0.5*(std::log(o_ks(i)) + std::log(std::complex<double>(k2 - K -2.0*K*k2*dE_E(i)/(E(i)*K+k2))));
    std::complex<double> logw1 = 0.5*(std::log(o_ks(i+1)) + std::log(std::complex<double>(k2 - K -2.0*K*k2*dE_E(i+1)/(E(i+1)*K+k2))));
    return std::exp(logw0+(logw1-logw0)*(N-Nstart-Ninc*i)/Ninc);
};

std::complex<double> g(double N){
    int i;
    i=int((N-Nstart)/Ninc);
    double k2;
    if(K>0)
        k2 = k*(k+2.0)-3.0*K;
    else
        k2 = k*k-3.0*K;
    double g0 = 0.5*(K*o_ks(i) + 3.0 - E(i) + dE_E(i)*k2/(E(i)*K+k2));
    double g1 = 0.5*(K*o_ks(i+1) + 3.0 - E(i+1) + dE_E(i+1)*k2/(E(i+1)*K+k2));
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
    // equation giving dy with y = [logomega_k, phi/mp], 
    // dy = [ dlogomega_k, dphi/mp]
    yp[0] = 4.0 + std::exp(y[0])*(4.0*K - 2.0*std::exp(2.0*t)*V(y[1]*mp)/(mp*mp));
    yp[1] = -std::sqrt(6.0 + std::exp(y[0])*(6.0*K -
    2.0*std::exp(2.0*t)*V(y[1]*mp)/(mp*mp)));
};

int main(){

    int nospectra=300;
    Eigen::VectorXd Nis=Eigen::VectorXd::Zero(nospectra);
    Eigen::VectorXd okis=Eigen::VectorXd::Zero(nospectra);
    Nis << 1.39,        1.37124214,  1.35247388,  1.33369522,  1.31490616,  1.29610671,
 1.27729686,  1.25847662,  1.23964597,  1.22080493,  1.20195349,  1.18309165, 
 1.16421942,  1.14533679,  1.12644376,  1.10754034,  1.08862651,  1.06970229,
 1.05076767,  1.03182266,  1.01286725,  0.99390144,  0.97492523,  0.95593862,
 0.93694162,  0.91793422,  0.89891643,  0.87988823,  0.86084964,  0.84180065,
 0.82274127,  0.80367148,  0.7845913 ,  0.76550072,  0.74639975,  0.72728837,
 0.7081666 ,  0.68903444,  0.66989187,  0.65073891,  0.63157555,  0.61240179,
 0.59321764,  0.57402308,  0.55481814,  0.53560279,  0.51637704,  0.4971409 ,
 0.47789436,  0.45863743,  0.4393701 ,  0.42009236,  0.40080424,  0.38150571,
 0.36219679,  0.34287747,  0.32354775,  0.30420764,  0.28485712,  0.26549621,
 0.24612491,  0.2267432 ,  0.2073511 ,  0.1879486 ,  0.16853571,  0.14911241,
 0.12967872,  0.11023463,  0.09078015,  0.07131526,  0.05183998,  0.0323543 ,
 0.01285823, -0.00664824, -0.02616511, -0.04569238, -0.06523005, -0.08477811,
-0.10433657, -0.12390542, -0.14348468, -0.16307433, -0.18267438, -0.20228482,
-0.22190567, -0.24153691, -0.26117855, -0.28083058, -0.30049302, -0.32016585,
-0.33984533, -0.35947241, -0.37902671, -0.39850823, -0.41791697, -0.43725293,
-0.4565161 , -0.47570649, -0.49482411, -0.51386894, -0.53284098, -0.55174025,
-0.57056674, -0.58932044, -0.60800136, -0.6266095 , -0.64514486, -0.66360744,
-0.68199724, -0.70031425, -0.71855848, -0.73672994, -0.75482861, -0.77285449,
-0.7908076 , -0.80868793, -0.82649547, -0.84423023, -0.86189222, -0.87948142,
-0.89699783, -0.91444147, -0.93181233, -0.9491104 , -0.96633569, -0.9834882 ,
-1.00056793, -1.01757488, -1.03450904, -1.05137043, -1.06815903, -1.08487485,
-1.10151789, -1.11808815, -1.13458563, -1.15101032, -1.16736224, -1.18364137,
-1.19984772, -1.21598129, -1.23204208, -1.24803008, -1.26394531, -1.27978775,
-1.29555741, -1.31125429, -1.32687839, -1.34242971, -1.35790825, -1.373314  ,
-1.38863385, -1.4037891 , -1.41876661, -1.4335664 , -1.44818846, -1.46263279,
-1.47689939, -1.49098827, -1.50489941, -1.51863283, -1.53218852, -1.54556648,
-1.55876671, -1.57178922, -1.58463399, -1.59730104, -1.60979036, -1.62210195,
-1.63423581, -1.64619195, -1.65797035, -1.66957103, -1.68099398, -1.6922392 ,
-1.70330669, -1.71419645, -1.72490849, -1.7354428 , -1.74579937, -1.75597822,
-1.76597935, -1.77580274, -1.7854484 , -1.79491634, -1.80420655, -1.81331903,
-1.82225378, -1.8310108 , -1.8395901 , -1.84799167, -1.8562155 , -1.86426161,
-1.87213   , -1.87982065, -1.88733357, -1.89466877, -1.90182624, -1.90880598,
-1.91560799, -1.92223227, -1.92867883, -1.93494765, -1.94103875, -1.94695212,
-1.95268776, -1.95824567, -1.96362586, -1.96882832, -1.97385304, -1.97870004,
-1.98340146, -1.98805046, -1.99265293, -1.9972089 , -2.00171834, -2.00618127,
-2.01059767, -2.01496757, -2.01929094, -2.0235678 , -2.02779814, -2.03198196,
-2.03611926, -2.04021005, -2.04425432, -2.04825207, -2.05220331, -2.05610802,
-2.05996623, -2.06377791, -2.06754307, -2.07126172, -2.07493385, -2.07855947,
-2.08213856, -2.08567114, -2.0891572 , -2.09259674, -2.09598977, -2.09933628,
-2.10263627, -2.10588975, -2.1090967 , -2.11225714, -2.11537106, -2.11843847,
-2.12145936, -2.12443373, -2.12736158, -2.13024291, -2.13307773, -2.13586603,
-2.13860781, -2.14130308, -2.14395183, -2.14655406, -2.14910977, -2.15161897,
-2.15408165, -2.15649781, -2.15886745, -2.16119058, -2.16346719, -2.16569728,
-2.16788085, -2.17001791, -2.17210845, -2.17415247, -2.17614998, -2.17810096,
-2.18000543, -2.18186339, -2.18367482, -2.18543974, -2.18715814, -2.18883002,
-2.19045539, -2.19203424, -2.19356657, -2.19505238, -2.19649168, -2.19788446,
-2.19923072, -2.20053046, -2.20178369, -2.2029904 , -2.20415059, -2.20526426,
-2.20633142, -2.20735206, -2.20832618, -2.20925379, -2.21013487, -2.21096944,
-2.2117575 , -2.21249903, -2.21319405, -2.21384255, -2.21444453, -2.215;

    for(int nsp=0;nsp<nospectra;nsp++){
//        Nis(nsp)=1.39-3.06/nospectra*nsp;
        okis(nsp)=1e-3*std::pow(1e6,nsp/double(nospectra));
    };

    for(int nsp=0; nsp<nospectra; nsp++){
        Ni=Nis(nsp);
        o_ki=okis(nsp);

        // Range of wavenumbers
        int nk=1500;
        Eigen::VectorXd ks=Eigen::VectorXd::LinSpaced(nk,-1,4);
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
        thresh[0] = 1.0e-10;
        thresh[1] = 1.0e-10;
        tol = 1.0e-8;
        nag_ode_ivp_rkts_setup(n, Ni, Nend, yinit, tol, thresh, method, errass, 0.0, iwsav, rwsav, &fail);
        twant=Ni;
        
        for(int i=mid+1; i<=(npts); i++){
            tnext = twant;
            twant+=Ninc;
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
        thresh[0] = 1.0e-10;
        thresh[1] = 1.0e-10;
        tol = 1.0e-8;
    
        std::cout << "starting backwards integration" << std::endl;
        nag_ode_ivp_rkts_setup(n, Ni, Nstart, yinit, tol, thresh, method, errass, 0.0, iwsav, rwsav, &fail);
        twant=Ni;
        
        for(int i=mid-1; i>=0; i--){
            tnext = twant;
            twant-=Ninc;
            if(twant < Nstart)
                twant = Nstart;
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
        //for(int i=0; i<nk; i++)
        //    ks(i) = int(ks(i)*atoday);
        //std::cout << "new scales: " << ks(0) << ".." << ks(nk-1) << std::endl;
     
        // Write background to file
        std::ofstream fbg;
        fbg.open("test/ms/kd-closed-bg-tip.txt");
        double ksq,kexample;
        kexample=ks(0);
        if(K>0)
            ksq = kexample*(kexample+2.0)-3*K;
        else
            ksq = kexample*kexample-3*K;
        for(int i=0; i<Nbg; i++){
            fbg << Ns(i) << ", " << phis(i) << ", " << dphis(i) << ", " <<
            logo_ks(i) << ", " << dlogo_ks(i) << ", " << (ksq - K
            -2.0*K*ksq*dE_E(i)/(E(i)*K+ksq)) << ", " << 0.5*(K*o_ks(i) + 3.0 - E(i) +
            dE_E(i)*ksq/(E(i)*K+ksq)) << ", " << dE_E(i) << std::endl;
        };
        fbg.close();
//        std::cout << "Done solving background" << std::endl;
//        std::string brkmsg;
//        std::cout << "press anything" << std::endl;
//        std::cin >> brkmsg;
    
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
        std::cout << "exits: " << exits(0) << ", " << exits(nk-1) << std::endl;
    
        // Construct system to solve
        de_system system(&w,&g);
        
        // Solve the evolution of each perturbation
        double ti, tf, rtol, atol, h0;
        std::complex<double> x0, dx0;
        int order=3;
        bool full_output=false;//true;
        rtol=1e-3;
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
            ic = bd(k,phi,dphi,ddphi,N);
            x0s(i) = ic(0);
            dx0s(i) = ic(1);
            a0 = std::exp(N);
            z = a0*dphi;
            dz_z = ddphi/dphi + 1.0;
            x0 = ic(0);
            dx0 = ic(1);
            x01s(i) = x0;
            dx01s(i) = dx0;
            Solution solution1(system, x0, dx0, ti, tf, order, rtol, atol, h0,
            full_output);
            t1 = MPI_Wtime();
            solution1.solve();
            rk1.push_back(solution1.sol.back());
            
            x0 = 0.0; 
            dx0 = k; 
            x02s(i) = x0;
            dx02s(i) = dx0;
            Solution solution2(system, x0, dx0, ti, tf, order, rtol, atol, h0, full_output);
            solution2.solve();
            t2 = MPI_Wtime();
            rk2.push_back(solution2.sol.back());
            times.push_back(t2-t1);
            //std::cout << "done k=" << k << " in " << t2-t1 << " s." << std::endl;
            
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
        f.open("test/ms/pps-kd-closed-kcts-corr"+std::to_string(nsp)+".txt");
        it_1 = rk1.begin();
        it_2 = rk2.begin();
        it_3 = times.begin();
        double starti;
        for(int u=0; u<nk; u++){
            starti = startindices(u);
            f << ks(u) << ", "  << *it_1 << ", " << *it_2 << ", " <<
            pps(ks(u),*it_1,*it_2,x01s(u),dx01s(u),x02s(u),dx02s(u),x0s(u),dx0s(u))
            << ", " << *it_3 << std::endl;
            ++it_1;
            ++it_2;
            ++it_3;
        };
        f.close();
    };

    return 0;
}
