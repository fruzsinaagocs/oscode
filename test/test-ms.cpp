#include "test-pps.hpp" 
#include "catch.hpp"
// Settings for the Mukhanov--Sasaki equation.
int n=2;
double m=4.51e-6;
double k=0.1;
double phi_p=23.293;
double mp=1;
double tstart = 1e4;

TEST_CASE("rkwkb","[ms]"){
    // Solving the MS equation with RKWKB. 
    std::string outputfile = "outputs/ms-rkwkb-sameN3.txt";
    double t=1.0;
    double rtol=1e-4;
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
    //std::cout << "total steps: " << my_solution.stepsall << std::endl;
    //std::cout << "waste ratio: " << my_solution.waste << std::endl;
};

TEST_CASE("nag", "[ms]"){

    std::string outputfile="outputs/ms-nag4.txt";
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
    yinit[2] = ybg[0].real();
    yinit[3] = ybg[1].real();
    yinit[4] = ybg[2].real();
    yinit[5] = ybg[3].real();
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
