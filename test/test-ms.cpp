#include "test-ms.hpp"
// Settings for the Mukhanov--Sasaki equation.
int n=2;
double m=1e-5;
double k=1000;
double phi_p=23.0;
double mp=1.0;

TEST_CASE("rkwkb","[ms]"){
    // Solving the MS equation with RKWKB. 
    double t = 1.0;
    Vector ic(6);                                      
    ic << 1.0, 0.0, background(t);
    de_system my_system(F, DF, w, Dw, DDw, g, Dg, DDg);   
    Solution my_solution(my_system, ic, t, f_end, 1, 1e-4);
    my_solution.evolve();
//    std::cout << "solution at " << my_solution.t << ": " << my_solution.y(0) << "," << my_solution.y(1) << " " << my_solution.y(2) << " " << my_solution.y(3) << " " << my_solution.y(4) << " " << my_solution.y(5) <<  std::endl;
//    std::cout << "total steps: " << my_solution.stepsall << std::endl;
//    std::cout << "waste ratio: " << my_solution.waste << std::endl;
};

TEST_CASE("nag", "[ms]"){
    // Solving the MS equation with a basic RK solver from the NAG C library.
    Integer liwsav, lrwsav, n, j;
    double hnext, hstart, tend, tgot, tol, tstart, tnext, twant, tinc, waste;
    Integer fevals, stepcost, stepsok, npts;
    double *rwsav = 0, *thresh = 0, *ygot = 0, *yinit = 0, *ymax = 0;
    double *ypgot = 0;
    Integer *iwsav = 0;
    NagError fail;
    Nag_RK_method method;
    Nag_ErrorAssess errass;
    Nag_Comm comm;
    std::string nag_output = "nag_output.txt";
                                                                                  
    INIT_FAIL(fail);

    n = 6; // dimensions of solution vector
    liwsav = 130;
    lrwsav = 350 + 32 * n;
                                                                                  
    // These are all arrays, NAG_ALLOC allocates memory to them.
    thresh = NAG_ALLOC(n, double);
    ygot = NAG_ALLOC(n, double);
    yinit = NAG_ALLOC(n, double);
    ypgot = NAG_ALLOC(n, double);
    ymax = NAG_ALLOC(n, double);
    iwsav = NAG_ALLOC(liwsav, Integer);
    rwsav = NAG_ALLOC(lrwsav, double);
    
    // Set initial conditions for ODE and parameters for the integrator. 
    // nag_enum_name_to_value (x04nac) Converts NAG enum member name to value. 
    method = (Nag_RK_method) nag_enum_name_to_value("Nag_RK_4_5");
    errass = (Nag_ErrorAssess) nag_enum_name_to_value("Nag_ErrorAssess_off");
    tstart = 1.0;
    tend = 100000;
    npts = 1;
    //npts = 100000;
    tinc = (tend - tstart)/(double) (npts);
    yinit[0] = 1.0;
    yinit[1] = 0.0;
    // This routine cannot handle complex numbers, but the solution is real
    // anyway, so take real part.
    yinit[2] = background(tstart)(0).real();
    yinit[3] = background(tstart)(1).real();
    yinit[4] = background(tstart)(2).real();
    yinit[5] = background(tstart)(3).real();
    hstart = 0.0;
    thresh[0] = 1.0e-8;
    thresh[1] = 1.0e-8;
    thresh[2] = 1.0e-8;
    thresh[3] = 1.0e-8;
    tol = 1.0e-8;
                                                                                  
    // Initialize Runge-Kutta method for integrating ODE using nag_ode_ivp_rkts_setup (d02pqc).
    nag_ode_ivp_rkts_setup(n, tstart, tend, yinit, tol, thresh, method,
                             errass, hstart, iwsav, rwsav, &fail);
                                                                                  
    twant = tstart;
    for(j=0; j<npts; j++){
        tnext = twant;
        twant += tinc;
        while(tnext < twant){
            nag_ode_ivp_rkts_range(f, n, twant, &tgot, ygot, ypgot, ymax, &comm, iwsav, rwsav, &fail);
            tnext = tgot;
        };
                                                                                  
        // write current step to file
        std::ofstream fout;
        fout.open(nag_output, std::ios_base::app);
        fout << tgot << " ";
        for (int k = 0; k < n; k++)
            fout <<  ygot[k] << " ";
        fout << std::endl;
        fout.close();

       // Get diagnostics on whole integration using nag_ode_ivp_rkts_diag (d02ptc).
        //nag_ode_ivp_rkts_diag(&fevals, &stepcost, &waste, &stepsok, &hnext, iwsav, rwsav, &fail);
    };
    nag_ode_ivp_rkts_diag(&fevals, &stepcost, &waste, &stepsok, &hnext, iwsav, rwsav, &fail);
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
