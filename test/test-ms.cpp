#include "test-pps.hpp" 
#include "catch.hpp"
// Settings for the Mukhanov--Sasaki equation.
int n=2;
double m=4.51e-6;
double k=0.1;
double phi_p=23.293;
double mp=1;
double tstart = 1e4;
double twant = 5e4;


TEST_CASE("background", "[ms]"){

    std::string outputfile = "outputs/ms-bg-midt.txt";
    double rtol=1e-4;
    double atol=1e-7;
    int order=1;
    double t=1.0;
    Vector ic(6);
    ic << 100*k, 0.0, background(tstart);
    de_system MSSystem(F, DF, w, Dw, DDw, g, Dg, DDg);
    Solution BGSolution(MSSystem,ic,t,until_start,order,rtol,atol,1.0);
    BGSolution.evolve();
    Vector ybg = BGSolution.y.segment(2,4);
    ic << 100*k, 0.0, ybg;
    Solution MSSolution(MSSystem,ic,tstart,until_twant,order,rtol,atol,1.0); 
    MSSolution.evolve();
    Vector ybg2 = MSSolution.y.head(6);
    ic << ybg2;
    double hstart = 1e-2;
    double hend = 1e4;
    int npts = 6000;
    double hinc = std::exp((std::log(hend)-std::log(hstart))/(double) npts);
    double h = hstart;
    std::cout << h << std::endl;
    Step rkf_step, wkb_step;
    double a, H, phi, dS0;
    std::ofstream fout;
    fout.open(outputfile, std::ios_base::app);
    fout << "# Background in MS as function of h at fixed t" << std::endl;
    fout << "# h, phi, a, H, dS0" << std::endl;
    
    for(int i=0; i<npts; i++){
        rkf_step = MSSolution.step(&MSSolution.rkfsolver, MSSolution.f_tot, MSSolution.y, h);
        MSSolution.wkbsolver->y0=MSSolution.y;
        MSSolution.wkbsolver->y1=rkf_step.y;
        MSSolution.wkbsolver->error0=Vector::Zero(9);
        MSSolution.wkbsolver->error1=rkf_step.error;
        wkb_step = MSSolution.step(MSSolution.wkbsolver,MSSolution.f_tot,MSSolution.y,h);
        phi = std::real(rkf_step.y(2));
        a = std::real(rkf_step.y(4));
        H = std::real(rkf_step.y(5));
        dS0 = std::imag(MSSolution.wkbsolver->ds_all(order-1));

        fout << h << " " << phi << " " << a << " " << H << " " << dS0 << std::endl;
        h *= hinc;
    };
    fout.close();

};

TEST_CASE("error","[ms]"){

    std::string outputfile = "outputs/ms-errors-midt.txt";
    double rtol=1e-4;
    double atol=1e-7;
    int order=1;
    double t=1.0;
    Vector ic(6);
    ic << 100*k, 0.0, background(tstart);
    de_system MSSystem(F, DF, w, Dw, DDw, g, Dg, DDg);
    Solution BGSolution(MSSystem,ic,t,until_start,order,rtol,atol,1.0);
    BGSolution.evolve();
    Vector ybg = BGSolution.y.segment(2,4);
    ic << 100*k, 0.0, ybg;
    Solution MSSolution(MSSystem,ic,tstart,until_twant,order,rtol,atol,1.0); 
    MSSolution.evolve();
    Vector ybg2 = MSSolution.y.head(6);
    ic << ybg2;
    double hstart = 1e-2;
    double hend = 1e4;
    int npts = 6000;
    double hinc = std::exp((std::log(hend)-std::log(hstart))/(double) npts);
    double h = hstart;
    std::cout << h << std::endl;
    Step rkf_step, wkb_step;
    Scalar truncerror, truncerror_crude, dS2, S1, S0;
    std::ofstream fout;
    fout.open(outputfile, std::ios_base::app);
    fout << "# Error progression in MS as function of h at fixed t" << std::endl;
    fout << "# h, RKF abs error, RKF rel error, WKB abs error, WKB rel error, WKB abs truncation error, WKB rel truncation error, WKB rel truncation error (crude), dS2, S1, S0" << std::endl;
    
    for(int i=0; i<npts; i++){
        rkf_step = MSSolution.step(&MSSolution.rkfsolver, MSSolution.f_tot, MSSolution.y, h);
        //for(int j=0; j<ic.size(); j++)
        //    CHECK(MSSolution.y(j) == ic(j));
        MSSolution.wkbsolver->y0=MSSolution.y;
        MSSolution.wkbsolver->y1=rkf_step.y;
        MSSolution.wkbsolver->error0=Vector::Zero(9);
        MSSolution.wkbsolver->error1=rkf_step.error;
        wkb_step = MSSolution.step(MSSolution.wkbsolver,MSSolution.f_tot,MSSolution.y,h);
        dS2 = MSSolution.wkbsolver->ds_all(order+1);
        S1 = MSSolution.wkbsolver->s_all(order);
        S0 = MSSolution.wkbsolver->s_all(order-1);
        truncerror = MSSolution.wkbsolver->trunc_error(0);
        truncerror_crude = MSSolution.wkbsolver->s_all(order+1);
        fout << h << " " << std::abs(rkf_step.error(0)) << " " << std::abs(rkf_step.error(0)/rkf_step.y(0)) << " " << std::abs(wkb_step.error(0)) << " " << std::abs(wkb_step.error(0)/wkb_step.y(0)) <<  " " << std::abs(truncerror) << " " << std::abs(truncerror/wkb_step.y(0)) << " " << std::abs(truncerror_crude) << " " << std::abs(dS2) << " " << std::abs(S1) << " " << std::abs(S0) << std::endl;
        h *= hinc;
    };
    fout.close();
};


TEST_CASE("rkwkb","[ms]"){
    // Solving the MS equation with RKWKB. 
    std::string outputfile = "outputs/ms-rkwkb-mixederr-trybg.txt";
    double t=1.0;
    double rtol=1e-6;
    double atol=1e-7;
    int order=1;
    Vector ic(6), ic1(6);
    std::cout << tstart << std::endl;
    // First solve background to t=tstart where we set i.c. for R_K
    ic << 100.0*k, 0.0, background(t);
    de_system MSSystem(F, DF, w, Dw, DDw, g, Dg, DDg);   
    Solution BGSolution(MSSystem,ic,t,until_start,order,rtol,atol,1.0);
    BGSolution.evolve();
    Vector ybg = BGSolution.y;
    
    // Then solve the MS up to after horizon exit. 
    ic1 << 100.0*k, 0.0, ybg.segment(2,4);
    Solution MSSolution(MSSystem,ic1,tstart,outside_horizon,order,rtol,atol,1.0,outputfile);
    MSSolution.evolve();
    //std::cout << "solution at " << my_solution.t << ": " << my_solution.y(0) << "," << my_solution.y(1) << " " << my_solution.y(2) << " " << my_solution.y(3) << " " << my_solution.y(4) << " " << my_solution.y(5) <<  std::endl;
    std::cout << "total steps: " << MSSolution.stepsall << std::endl;
    std::cout << "waste ratio: " << MSSolution.waste << std::endl;
    std::cout << MSSolution.stepsstr << std::endl;
};

TEST_CASE("nag", "[ms]"){

    std::string outputfile="outputs/ms-nag-truncerr.txt";
    std::ofstream fout;
    fout.open(outputfile, std::ios_base::app);
    fout << "# MS with NAG" << std::endl;
    fout << "# parameters: n=" << n << ", m=" << m << ", phi_p=" << phi_p << ", mp=" << mp << ", ic set at tstart=" << tstart << ", rtol=1e-4" << ", atol=1e-7" << ", order=1" << std::endl;
    fout << "# t, R_k" << std::endl;

    // NAG integrator setup
    Integer liwsav, lrwsav, n, j;
    double hnext, hstart, tend, tgot, tol, tnext, twant, tinc, waste;
    Integer fevals, stepcost, stepsok, npts;
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
    
    // First integrate up to tstart
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
    nag_ode_ivp_rkts_setup(n, 1.0, tstart, yinit, 1.0e-7, thresh, method, errass, 0.0, iwsav, rwsav, &fail);
    tnext=1.0;
    while(tnext < tstart){
        nag_ode_ivp_rkts_range(f,n,tstart,&tgot,ygot,ypgot,ymax,&comm,iwsav,rwsav,&fail);
        tnext = tgot; 
    };
    Vector ybg(4);
    ybg << ygot[2], ygot[3], ygot[4], ygot[5];
 
    NAG_FREE(thresh);
    NAG_FREE(yinit);
    NAG_FREE(ygot);
    NAG_FREE(ypgot);
    NAG_FREE(ymax);
    NAG_FREE(rwsav);
    NAG_FREE(iwsav);

    // Then solve the MS from tstart until horizon exit.
     
    thresh = NAG_ALLOC(n, double);
    ygot = NAG_ALLOC(n, double);
    yinit = NAG_ALLOC(n, double);
    ypgot = NAG_ALLOC(n, double);
    ymax = NAG_ALLOC(n, double);
    iwsav = NAG_ALLOC(liwsav, Integer);
    rwsav = NAG_ALLOC(lrwsav, double);
    INIT_FAIL(fail); 
    
    tend = 1e+10;
    npts = 1e5;
    tinc = std::exp((std::log(tend) - std::log(tstart))/(double) (npts));
    yinit[0] = 100.0*k;
    yinit[1] = 0.0;
    yinit[2] = background(tstart)(0).real();
    yinit[3] = background(tstart)(1).real();
    yinit[4] = background(tstart)(2).real();
    yinit[5] = background(tstart)(3).real();
    hstart = 0.0;
    thresh[0] = 1.0e-7;
    thresh[1] = 1.0e-7;
    thresh[2] = 1.0e-7;
    thresh[3] = 1.0e-7;
    thresh[4] = 1.0e-7;
    thresh[5] = 1.0e-7;
    tol = 1.0e-7;
    nag_ode_ivp_rkts_setup(n, tstart, tend, yinit, tol, thresh, method, errass, hstart, iwsav, rwsav, &fail);
    
    std::cout << tstart << "," << tend << "," << ybg << std::endl;
    twant = tstart;
    for(j=0; j<npts; j++){
        tnext = twant;
        twant *= tinc;
        while(tnext < twant){
            nag_ode_ivp_rkts_range(f, n, twant, &tgot, ygot, ypgot, ymax, &comm, iwsav, rwsav, &fail);
            tnext = tgot;
        };
                                                                                  
        // write current step to file
        std::ofstream fout;
        fout.open(outputfile, std::ios_base::app);
        fout << tgot << " " << ygot[0] << std::endl;
        fout.close();
        if(ygot[4]*ygot[5] > 100.0*k)
            break;

       // Get diagnostics on whole integration using nag_ode_ivp_rkts_diag (d02ptc).
        //nag_ode_ivp_rkts_diag(&fevals, &stepcost, &waste, &stepsok, &hnext, iwsav, rwsav, &fail);
    };
    //nag_ode_ivp_rkts_diag(&fevals, &stepcost, &waste, &stepsok, &hnext, iwsav, rwsav, &fail);
    //std::cout << "solution: " << ygot[0] << ", " << ygot[1] << ", " << ygot[2] << ", " << ygot[3] << ", " <<  ygot[4] << ", "<<  ygot[5] << ", " << std::endl;   
    //std::cout << "total steps: " << stepsok*(1.0/(1.0 - waste)) << std::endl;
    //std::cout << "waste ratio: " << waste << std::endl;

    NAG_FREE(thresh);
    NAG_FREE(yinit);
    NAG_FREE(ygot);
    NAG_FREE(ypgot);
    NAG_FREE(ymax);
    NAG_FREE(rwsav);
    NAG_FREE(iwsav);
}; 
