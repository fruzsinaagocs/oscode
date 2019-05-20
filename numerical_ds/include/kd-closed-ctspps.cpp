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
double Nstart=0.0, Nend=65.0, Ninc=(Nend-Nstart)/npts;
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
    std::complex<double> w0 = std::sqrt(o_ks(i)*std::log(std::complex<double>(k2 - K -2.0*K*k2*dE_E(i)/(E(i)*K+k2))));
    std::complex<double> w1 = std::sqrt(o_ks(i+1)*std::log(std::complex<double>(k2 - K -2.0*K*k2*dE_E(i+1)/(E(i+1)*K+k2))));
    return (w0+(w1-w0)*(N-Nstart-Ninc*i)/Ninc); 
//    std::complex<double> logw0 = 0.5*(std::log(o_ks(i)) + std::log(std::complex<double>(k2 - K -2.0*K*k2*dE_E(i)/(E(i)*K+k2))));
//    std::complex<double> logw1 = 0.5*(std::log(o_ks(i+1)) + std::log(std::complex<double>(k2 - K -2.0*K*k2*dE_E(i+1)/(E(i+1)*K+k2))));
//    return std::exp(logw0+(logw1-logw0)*(N-Nstart-Ninc*i)/Ninc);
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
    Nis << 1.38031136428, 1.36042399313, 1.34053741753, 1.32065166917,
    1.30076678155,
    1.28088278954, 1.26099972911, 1.24111763858, 1.22123655698, 1.20135652516,
    1.18147758563, 1.16159978264, 1.14172316223, 1.12184777224, 1.1019736625,
    1.08210088492, 1.06222949341, 1.04235954399, 1.02249109499, 1.00262420711,
    0.982758943476, 0.962895369764, 0.943033554205, 0.923173568017,
    0.903315485301, 0.883459381721, 0.863605339029, 0.843753439753,
    0.823903770272, 0.804056420725, 0.784211484739, 0.764369059757,
    0.74452924695, 0.724692151546, 0.704857883495, 0.685026556178,
    0.665198288699, 0.645373204096, 0.625551430167, 0.605733100908,
    0.585918353701, 0.566107332956, 0.5463001892, 0.52649707789, 0.506698161214,
    0.486903607084, 0.467113590748, 0.447328294022, 0.427547906285,
    0.407772622063, 0.388002651138, 0.368238202214, 0.348479497168,
    0.328726765812, 0.308980245091, 0.289240184521, 0.269506842375,
    0.249780484807, 0.230061390336, 0.210349848657, 0.190646160142,
    0.170950635048, 0.151263601713, 0.131585392547, 0.111916359475,
    0.0922568654613, 0.0726072873725, 0.0529680169669, 0.0333394595625,
    0.0137220408642, -0.00588380296113, -0.0254776177663, -0.0450589304482,
    -0.0646272504491, -0.0841820695199, -0.103722859874, -0.123249074692,
    -0.142760143522, -0.162255480346, -0.181734472198, -0.201196486369,
    -0.220640865577, -0.240066930222, -0.259473975794, -0.278861265443,
    -0.298228046246, -0.317573530003, -0.33689690814, -0.356197334834,
    -0.375473939111, -0.394725818831, -0.413952040792, -0.43315163789,
    -0.452323609754, -0.471466922548, -0.490580507737, -0.509663259472,
    -0.528714035989, -0.547731658124, -0.566714907264, -0.585662523961,
    -0.604573207256, -0.62344561829, -0.642278374093, -0.661070049531,
    -0.679819174832, -0.698524234228, -0.71718366776, -0.735795870808,
    -0.754359183588, -0.772871909715, -0.791332300959, -0.809738557805,
    -0.828088833906, -0.846381232804, -0.864613810133, -0.882784568511,
    -0.900891465021, -0.91893239438, -0.936905215596, -0.954807735176,
    -0.972637705163, -0.990392831287, -1.00807077069, -1.02566913439,
    -1.04318548324, -1.06061733481, -1.07796216758, -1.09521739719,
    -1.11238042223, -1.12944858908, -1.14641921024, -1.1632895629,
    -1.18005688872, -1.19671840137, -1.21327128917, -1.22971271195,
    -1.24603981068, -1.26224971312, -1.27833952686, -1.29430635452,
    -1.31014729834, -1.32585944544, -1.34143989978, -1.35688577296,
    -1.3721941886, -1.38736228948, -1.40238724719, -1.41726625822,
    -1.43199656048, -1.44657543227, -1.46100019917, -1.47526824585,
    -1.48937701133, -1.50332400738, -1.51710681433, -1.53072309467,
    -1.54417059518, -1.55744715357, -1.57055070634, -1.58347929094,
    -1.59623105614, -1.60880426187, -1.62119729037, -1.63340864335,
    -1.64543695783, -1.65728099897, -1.66893967055, -1.68041201519,
    -1.6916972209, -1.70279462098, -1.71370369898, -1.7244240864,
    -1.73495556692, -1.7452980776, -1.75545171058, -1.76541670638,
    -1.77519346136, -1.78478252235, -1.79418458547, -1.80340049575,
    -1.81243124211, -1.82127795986, -1.82994191798, -1.83842452647,
    -1.84672732483, -1.85485197836, -1.86280027712, -1.87057412872,
    -1.87817555305, -1.88560667611, -1.89286972575, -1.8999670275,
    -1.90690099338, -1.91367412, -1.92028898158, -1.92674822571, -1.93305456048,
    -1.93921075925, -1.94521964154, -1.95108408029, -1.956806984, -1.962391299,
    -1.96784000258, -1.9731560923, -1.97834258731, -1.98340251907,
    -1.98833892765, -1.99315485767, -1.99785335177, -2.0024374483,
    -2.00691017624, -2.0112745515, -2.01553357537, -2.01969022546,
    -2.02374745958, -2.02770820909, -2.0315753752, -2.03535183147,
    -2.03904041543, -2.04264393095, -2.04616514576, -2.04960679009,
    -2.0529715524, -2.05626208289, -2.05948099095, -2.06263084162,
    -2.06571416028, -2.06873342963, -2.07169108677, -2.07458952794,
    -2.07743110893, -2.08021813955, -2.08295288958, -2.08563759007,
    -2.08827442789, -2.09086555295, -2.09341307688, -2.09591907349,
    -2.09838558238, -2.10081460603, -2.10320811647, -2.10556805285,
    -2.10789632887, -2.11019482465, -2.112465399, -2.11470988849, -2.1169301044,
    -2.11912784382, -2.12130488528, -2.12346299531, -2.12560393163,
    -2.12772944185, -2.1298412738, -2.13194117157, -2.13403088627,
    -2.13611217553, -2.1381868085, -2.14025657625, -2.1423232859,
    -2.14438877771, -2.14645492481, -2.14852364105, -2.15059689046,
    -2.15267669232, -2.15476513228, -2.15686437469, -2.15897666552,
    -2.16110435689, -2.1632499089, -2.16541591416, -2.16760511229,
    -2.16982040874, -2.17206490293, -2.17434190611, -2.1766549796,
    -2.17900796571, -2.18140502842, -2.18385070059, -2.18634993799,
    -2.18890818649, -2.19153145071, -2.1942263944, -2.19700043574,
    -2.19986188495, -2.20282009354, -2.2058856417, -2.20907056827,
    -2.21238864616, -2.21585572502, -2.21949015819, -2.22331333002,
    -2.22735030858, -2.23163067704, -2.23618955886, -2.24106890214,
    -2.2463190626, -2.2520007219, -2.25818710394, -2.26496635029; 

    for(int nsp=0;nsp<nospectra;nsp++){
//        Nis(nsp)=1.39-3.06/nospectra*nsp;
        okis(nsp)=1e-3*std::pow(147910.8388,nsp/double(nospectra));
    };

    for(int nsp=0; nsp<1; nsp++){
//      for(int nsp=0; nsp<nospectra; nsp++){
        Ni=Nis(nsp);
        o_ki=okis(nsp);

        // Range of wavenumbers
        int nk=1000;
        
        // CTS k
        Eigen::VectorXd ks=Eigen::VectorXd::LinSpaced(nk,-2,2);
        for(int i=0; i<nk; i++)
            ks(i) = std::pow(10,ks(i));

        // INT k
//        Eigen::VectorXd ks=Eigen::VectorXd::Zero(nk);
//        Eigen::VectorXd kscts=Eigen::VectorXd::LinSpaced(nk-100,2,4);
//        for(int i=1; i<=100; i++)
//            ks(i-1) = i;
//        for(int i=100; i<nk; i++)
//            ks(i) = std::pow(10,kscts(i-100)); 
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
        kexample=ks(500);
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
                if(Ns(j) >= Ni){
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
        rtol=1e-6;
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
            // Avoid WKB steps at k=1
            if(k < 1.00001)
                order=3;
            else
                order=3;
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
            x0 = 1.0;
            dx0 = 0.0;
            x01s(i) = x0;
            dx01s(i) = dx0;
            Solution solution1(system, x0, dx0, ti, tf, order, rtol, atol, h0,
            full_output);
            t1 = MPI_Wtime();
            solution1.solve();
            rk1.push_back(solution1.sol.back());
            
            x0 = 0.0; 
            dx0 = 1.0; 
            x02s(i) = x0;
            dx02s(i) = dx0;
            Solution solution2(system, x0, dx0, ti, tf, order, rtol, atol, h0, full_output);
            solution2.solve();
            t2 = MPI_Wtime();
            rk2.push_back(solution2.sol.back());
            times.push_back(t2-t1);
            //if(k==1.0)
            std::cout << nsp <<" done k=" << k << " with " << solution1.wkbsteps << "," << solution2.wkbsteps << " wkb steps, from " << ti << " to " << tf << std::endl;
            
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
        f.open("test/ms/pps-test.txt"); 
//        f.open("test/ms/pps-kd-closed-kcts-corr"+std::to_string(nsp)+".txt");
//        f.open("test/ms/pps-kd-closed-kint-corr"+std::to_string(nsp)+".txt");
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
