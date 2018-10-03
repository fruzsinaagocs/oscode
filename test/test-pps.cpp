#include "test-pps.hpp"
#include <chrono>

// Settings for the Mukhanov--Sasaki equation.
int n=2;
double m=4.5e-6;
double k;
double phi_p=23.293;
double mp=1.0;
double tstart=1e4;


void test_pps_nag(){

    std::string outputfile = "outputs/test-pps-nag-lowtol3.txt";
    std::ofstream fout;
    fout.open(outputfile, std::ios_base::app);
    fout << "# PPS with NAG" << std::endl;
    fout << "# parameters: n=" << n << ", m=" << m << ", phi_p=" << phi_p << ", mp=" << mp << ", ic set at tstart=" << tstart << ", rtol=1e-4" << ", atol=1e-7" << ", order=1" << std::endl;
    fout << "# k, R_k1, R_k2, P_R(HD), P_R(RST), dt" << std::endl;
 
    // NAG integrator setup
    Integer liwsav, lrwsav, n, j;
    double hnext, tend, tgot, tol, twant, waste, tnext;
    Integer fevals, stepcost, stepsok;
    double *rwsav = 0, *thresh = 0, *ygot = 0, *yinit = 0, *ymax = 0;
    double *ypgot = 0;
    Integer *iwsav = 0;
    NagError fail;
    Nag_RK_method method;
    Nag_ErrorAssess errass;
    Nag_Comm comm;
    INIT_FAIL(fail);
    n = 6;
    liwsav = 130;
    lrwsav = 350 + 32 * n;
    thresh = NAG_ALLOC(n, double);
    ygot = NAG_ALLOC(n, double);
    yinit = NAG_ALLOC(n, double);
    ypgot = NAG_ALLOC(n, double);
    ymax = NAG_ALLOC(n, double);
    iwsav = NAG_ALLOC(liwsav, Integer);
    rwsav = NAG_ALLOC(lrwsav, double);
    method = (Nag_RK_method) nag_enum_name_to_value("Nag_RK_4_5");
    errass = (Nag_ErrorAssess) nag_enum_name_to_value("Nag_ErrorAssess_off");
    
    // First integrate up to tstart, where we set initial conditions 
    tend = 1.0e4;
    k=0.1; 
    yinit[0] = 100.0*k;
    yinit[1] = 0.0;
    yinit[2] = background(1.0)(0).real();
    yinit[3] = background(1.0)(1).real();
    yinit[4] = background(1.0)(2).real();
    yinit[5] = background(1.0)(3).real();
    thresh[0] = 1.0e-7;
    thresh[1] = 1.0e-7;
    thresh[2] = 1.0e-7;
    thresh[3] = 1.0e-7;
    thresh[4] = 1.0e-7;
    thresh[5] = 1.0e-7;
    tol = 1.0e-7;
    nag_ode_ivp_rkts_setup(n, 1.0, tend, yinit, 1.0e-7, thresh, method, errass, 0.0, iwsav, rwsav, &fail);
    tnext=1.0;
    while(tnext < tend){
        nag_ode_ivp_rkts_range(f,n,tend,&tgot,ygot,ypgot,ymax,&comm,iwsav,rwsav,&fail);
        tnext = tgot; 
    };
    Vector ybg(4);
    Vector dybg(4);
    ybg << ygot[2], ygot[3], ygot[4], ygot[5];
    dybg << ypgot[2], ypgot[3], ypgot[4], ypgot[5];
 
    NAG_FREE(thresh);
    NAG_FREE(yinit);
    NAG_FREE(ygot);
    NAG_FREE(ypgot);
    NAG_FREE(ymax);
    NAG_FREE(rwsav);
    NAG_FREE(iwsav);

        
    // Then for each k mode, solve the MS.
    int npts = 6000;
    double kmin = 1e-3;
    double kmax = 1e+3;
    double kinc = std::exp((std::log(kmax) - std::log(kmin))/(double) npts);
    double Rk1, Rk2, power1, power2;
    tend = 1.0e+10;
    k = kmin;
    
    for(int i=0; i<npts; i++){
        thresh = NAG_ALLOC(n, double);
        ygot = NAG_ALLOC(n, double);
        yinit = NAG_ALLOC(n, double);
        ypgot = NAG_ALLOC(n, double);
        ymax = NAG_ALLOC(n, double);
        iwsav = NAG_ALLOC(liwsav, Integer);
        rwsav = NAG_ALLOC(lrwsav, double);
        INIT_FAIL(fail);
 
        yinit[0] = 100.0*k;
        yinit[1] = 0.0;
        yinit[2] = ybg[0].real(); 
        yinit[3] = ybg[1].real();
        yinit[4] = ybg[2].real();
        yinit[5] = ybg[3].real();
        thresh[0] = 1.0e-7;
        thresh[1] = 1.0e-7;
        thresh[2] = 1.0e-7;
        thresh[3] = 1.0e-7;
        thresh[4] = 1.0e-7;
        thresh[5] = 1.0e-7;
        tol = 1.0e-7;

        nag_ode_ivp_rkts_setup(n, tstart, tend, yinit, tol, thresh, method, errass, 0.0, iwsav, rwsav, &fail);
        auto start1 = std::chrono::system_clock::now(); 
        tnext=tstart;
        while(tnext < tend){
            nag_ode_ivp_rkts_range(f, n, tend, &tgot, ygot, ypgot, ymax, &comm, iwsav, rwsav, &fail);
            tnext = tgot;
            if(ygot[4]*ygot[5] > 100.0*k){
                break;
            };
        };
        auto end1 = std::chrono::system_clock::now();
        std::chrono::duration<double> elapsed1 = (end1-start1);
        Rk1 = ygot[0];

        NAG_FREE(thresh);
        NAG_FREE(yinit);
        NAG_FREE(ygot);
        NAG_FREE(ypgot);
        NAG_FREE(ymax);
        NAG_FREE(rwsav);
        NAG_FREE(iwsav);

        thresh = NAG_ALLOC(n, double);
        ygot = NAG_ALLOC(n, double);
        yinit = NAG_ALLOC(n, double);
        ypgot = NAG_ALLOC(n, double);
        ymax = NAG_ALLOC(n, double);
        iwsav = NAG_ALLOC(liwsav, Integer);
        rwsav = NAG_ALLOC(lrwsav, double);

        yinit[0] = 0.0;
        yinit[1] = 10.0*std::pow(k,2);
        yinit[2] = ybg[0].real(); 
        yinit[3] = ybg[1].real();
        yinit[4] = ybg[2].real();
        yinit[5] = ybg[3].real();
        thresh[0] = 1.0e-7;
        thresh[1] = 1.0e-7;
        thresh[2] = 1.0e-7;
        thresh[3] = 1.0e-7;
        thresh[4] = 1.0e-7;
        thresh[5] = 1.0e-7;
        tol = 1.0e-7;
    
        nag_ode_ivp_rkts_setup(n, tstart, tend, yinit, tol, thresh, method, errass, 0.0, iwsav, rwsav, &fail);
        auto start2 = std::chrono::system_clock::now(); 
        tnext=tstart;
        while(tnext < tend){
            nag_ode_ivp_rkts_range(f, n, tend, &tgot, ygot, ypgot, ymax, &comm, iwsav, rwsav, &fail);
            tnext = tgot;
            if(ygot[4]*ygot[5] > 100.0*k)
                break;
        };
        auto end2 = std::chrono::system_clock::now();
        std::chrono::duration<double> elapsed2 = (end2-start2);
        Rk2 = ygot[0]; 
 
        NAG_FREE(thresh);
        NAG_FREE(yinit);
        NAG_FREE(ygot);
        NAG_FREE(ypgot);
        NAG_FREE(ymax);
        NAG_FREE(rwsav);
        NAG_FREE(iwsav);

        power1 = HD(k, Rk1, Rk2, ybg, dybg);
        power2 = RST(k, Rk1, Rk2, ybg, dybg);
        std::cout << k << " " << Rk1 << " " << Rk2 << " " << power1 << " " << power2 << " " << elapsed1.count()+elapsed2.count() << std::endl;
        fout << k << " " << Rk1 << " " << Rk2 << " " << power1 << " " << power2  << " " << elapsed1.count()+elapsed2.count() << std::endl;
        k *= kinc; 
    };

};

void test_pps_rkwkb(){

    std::string outputfile = "outputs/pps-newerr3.txt";
    std::ofstream fout;
    int npts = 20;
    double kmin = 4e+7;
    double kmax = 1e+8;
    double kinc = std::exp((std::log(kmax) - std::log(kmin))/(double) npts);
    double Rk1, Rk2, power1, power2;
    double t = 1.0;
    double rtol=1e-3;
    double atol=1e-7;
    int order=1;
    
    fout.open(outputfile, std::ios_base::app);
    fout << "# PPS" << std::endl;
    fout << "# parameters: n=" << n << ", m=" << m << ", phi_p=" << phi_p << ", mp=" << mp << ", ic set at tstart=" << tstart << ", rtol=" << rtol << ", atol=" << atol << ", order=" << order << std::endl;
    fout << "# k, R_k1, R_k2, P_R(HD), P_R(RST), dt" << std::endl;
    
    // First solve until t=tstart to see what background is there.
    Vector ic(6), ic1(6), ic2(6);                                      
    k = 0.1;
    ic << 1000.0*k, 0.0, background(t);
    de_system MSSystem(F, DF, w, Dw, DDw, g, Dg, DDg);   
    Solution BGSolution(MSSystem, ic, t, until_start, order, 1e-7, atol, 1.0);
    BGSolution.evolve();
    Vector ybg = BGSolution.y;
    Vector dybg = BGSolution.f_tot(ybg); 
    k = kmin;
    
    // Then solve the Mukhanov--Sasaki equation for each k-mode, obtaining two
    // linearly indepdendent solutions.
    for(int i=0; i<npts; i++){
        ic1 << 100.0*k, 0.0, ybg.segment(2,4);
        ic2 << 0.0, 10*std::pow(k,2), ybg.segment(2,4);
        
        auto start1 = std::chrono::system_clock::now();
        Solution MSSolution1(MSSystem, ic1, tstart, outside_horizon, order, rtol, atol, 1.0);
        MSSolution1.evolve();
        auto end1 = std::chrono::system_clock::now();
        std::chrono::duration<double> elapsed1 = (end1-start1);
        
        auto start2 = std::chrono::system_clock::now();
        Solution MSSolution2(MSSystem, ic2, tstart, outside_horizon, order, rtol, atol, 1.0);
        MSSolution2.evolve();
        auto end2 = std::chrono::system_clock::now();
        std::chrono::duration<double> elapsed2 = (end2-start2);

        Rk1 = std::real(MSSolution1.y(0));
        Rk2 = std::real(MSSolution2.y(0)); 
        power1 = HD(k, Rk1, Rk2, ybg.segment(2,4), dybg.segment(2,4));
        power2 = RST(k, Rk1, Rk2, ybg.segment(2,4), dybg.segment(2,4));
        std::cout << k << " " << Rk1 << " " << Rk2 << " " << power1 << " " << power2 << " " << elapsed1.count()+elapsed2.count() << std::endl;
        fout << k << " " << Rk1 << " " << Rk2 << " " << power1 << " " << power2  << " " << elapsed1.count()+elapsed2.count() << std::endl;
        k *= kinc;  
    };
    fout.close();
     
};

int main(){

    test_pps_nag();
    return 0;
};

