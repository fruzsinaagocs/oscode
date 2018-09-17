#include "test-airy.hpp"

Vector background(Scalar t, int d){
    // Gives background at a given time t
    Vector result = Vector::Ones(d);
    return result*std::pow(t, 1.0/4.0);
};

Vector ana_s(Scalar t, int order){
    Vector result(4);
    result << std::complex<double>(0.0, 1.0)*2.0/3.0*std::pow(t, 3.0/2.0), -1.0/4.0*std::log(t), std::complex<double>(0.0, 1.0)*(-5.0)/48.0*std::pow(t, -3.0/2.0), -5.0/64.0*std::pow(t, -3.0);
    return result.head(order+2);
}

Vector ana_ds(Scalar t, int order){
    Vector result(4);
    result << std::complex<double>(0.0, 1.0)*std::pow(t, 1/2.0), -1.0/(4.0*t), std::complex<double>(0.0, 1.0)*5.0/32.0*std::pow(t, -5.0/2.0), 5.0/192.0*std::pow(t, -4);
    return result.head(order+2);
}

Vector ana_dds(Scalar t, int order){
    Vector result(3);
    result << 0.5*std::complex<double>(0.0, 1.0)*std::pow(t, -0.5), 0.25*std::pow(t, -2), -25.0/64.0*std::complex<double>(0.0, 1.0)*std::pow(t, -7.0/2.0);
    return result.head(order+2); 
}

Vector airy(double t){
    Vector result(2);
//    result << std::complex<double>(boost::math::airy_ai(-t), 0.0), std::complex<double>(-boost::math::airy_ai_prime(-t), 0.0);
    result << std::complex<double>(boost::math::airy_ai(-t), boost::math::airy_bi(-t)), std::complex<double>(-boost::math::airy_ai_prime(-t), -boost::math::airy_bi_prime(-t));
    return result; 
}


TEST_CASE("setting up a system of ODEs"){
    
    double t = 1.0;
    Vector ic(4);                                        
    ic << 1.0, 1.0, background(t, 2);
    Vector y = ic.tail(2);
    de_system my_system(F, DF, w, Dw, DDw, g, Dg, DDg);
     
    REQUIRE( std::real(my_system.w(y)) == Approx(std::pow(t, 1/2.0)) );
    REQUIRE( std::real(my_system.dw(y)) == Approx(1/2.0*std::pow(t, -1/2.0)) );
    REQUIRE( std::real(my_system.ddw(y)) == Approx(-1/4.0*std::pow(t, -3/2.0)) );
    REQUIRE( std::real(my_system.dg(y)) == Approx(0.0) );
    REQUIRE( std::real(my_system.ddg(y)) == Approx(0.0) );

};

TEST_CASE("s_ds_dds"){

    int order = 1;
    double t = 50000;
    Vector ic(4);                                      
    int d = ic.size()+(order+2);
    ic << airy(t), background(t, 2);
    de_system my_system(F, DF, w, Dw, DDw, g, Dg, DDg);   
    Solution my_solution(my_system, ic, 1000.0, f_end);
    REQUIRE(my_solution.y.size() == d);
    REQUIRE(my_solution.y.segment(2, d-(order+2)-2).size() == 2);
    Vector y_bg = my_solution.y.segment(2, d-(order+2)-2);
    REQUIRE(my_solution.wkbsolver->dS(y_bg).size() == order+2);
    for(int i=0; i<(order+2); i++){
        CHECK(std::real(my_solution.wkbsolver->dS(y_bg)(i)) == Approx(std::real(ana_ds(t,order)(i))) );
        CHECK(std::imag(my_solution.wkbsolver->dS(y_bg)(i)) == Approx(std::imag(ana_ds(t,order)(i))) );
        };
    for(int i=0; i<(order+1); i++){
        CHECK(std::real(my_solution.wkbsolver->ddS(y_bg)(i)) == Approx(std::real(ana_dds(t,order)(i))) );
        CHECK(std::imag(my_solution.wkbsolver->ddS(y_bg)(i)) == Approx(std::imag(ana_dds(t,order)(i))) );
    };
};

TEST_CASE("Single RKF step"){

    int order = 1;
    double t = 1.0;
    double step = 0.1;
    Vector ic(4);                                      
    int d = ic.size()+(order+2);
    ic << airy(t), background(t, 2);
    de_system my_system(F, DF, w, Dw, DDw, g, Dg, DDg);   
    Solution my_solution(my_system, ic, 1.0, f_end);
    Step rkf_step = my_solution.step(&my_solution.rkfsolver, my_solution.f_tot, my_solution.y, step); 
    t+=step;
    for(int i=0; i<2; i++){
        CHECK(std::real(rkf_step.y(i)) == Approx(std::real(airy(t)(i))));
        CHECK(std::imag(rkf_step.y(i)) == Approx(std::imag(airy(t)(i))));
    };
    for(int i=2; i<(d-(order+2)); i++){
        CHECK(std::real(rkf_step.y(i)) == Approx(std::real(background(t,2)(0))));
        CHECK(std::imag(rkf_step.y(i)) == Approx(std::imag(background(t,2)(0))));
    };
    for(int i=(d-(order+2)); i<d; i++){
        CHECK(std::real(rkf_step.y(i)) == Approx(std::real(ana_s(t,order)(i-d+(order+2)) - ana_s(t-step,order)(i-d+(order+2)))));
        CHECK(std::imag(rkf_step.y(i)) == Approx(std::imag(ana_s(t,order)(i-d+(order+2)) - ana_s(t-step,order)(i-d+(order+2)))));
    };
};

TEST_CASE("single_wkb_step"){

    int order = 1;
    double t = 101.0;
    double step = 899.0;
    Vector ic(4);                                      
    int d = ic.size()+(order+2);
    ic << airy(t), background(t, 2);
    de_system my_system(F, DF, w, Dw, DDw, g, Dg, DDg);   
    Solution my_solution(my_system, ic, t, f_end);
    Vector y0 = my_solution.y;
    REQUIRE(my_solution.f_tot(my_solution.y).size() == d);
    Step rkf_step = my_solution.step(&my_solution.rkfsolver, my_solution.f_tot, my_solution.y, step); 
    REQUIRE(rkf_step.y.size() == d);
    my_solution.wkbsolver->y0 = y0;
    my_solution.wkbsolver->y1 = rkf_step.y;
    my_solution.wkbsolver->error1 = rkf_step.error;
    my_solution.wkbsolver->error0 = Vector::Zero(d);
    SECTION("wkb step"){
        Step wkb_step = my_solution.step(my_solution.wkbsolver, my_solution.f_tot, my_solution.y, step);
        t+=step;
        for(int i=0; i<2; i++){
            CHECK(std::real(wkb_step.y(i)) == Approx(std::real(airy(t)(i))));
            CHECK(std::imag(wkb_step.y(i)) == Approx(std::imag(airy(t)(i))));
        };
        for(int i=2; i<(d-(order+2)); i++){
            CHECK(std::real(wkb_step.y(i)) == Approx(std::real(background(t,2)(0))));
            CHECK(std::imag(wkb_step.y(i)) == Approx(std::imag(background(t,2)(0))));
        };
        for(int i=(d-(order+2)); i<d; i++){
            CHECK(std::real(wkb_step.y(i)) == Approx(std::real(ana_s(t,order)(i-d+(order+2)) - ana_s(t-step,order)(i-d+(order+2)))));
            CHECK(std::imag(wkb_step.y(i)) == Approx(std::imag(ana_s(t,order)(i-d+(order+2)) - ana_s(t-step,order)(i-d+(order+2)))));
        };
        CHECK(my_solution.wkbsolver->trunc_error.size() == 2);
    };

};

TEST_CASE("rkwkb"){
    double t = 1.0;
    Vector ic(4);                                      
    ic << airy(t), background(t, 2);
    de_system my_system(F, DF, w, Dw, DDw, g, Dg, DDg);   
    Solution my_solution(my_system, ic, t, f_end, 1, 1e-4, 0.0, 1, "outputs/airy-rkwkb.txt");
    my_solution.evolve();
};

TEST_CASE("nag"){

    Integer liwsav, lrwsav, n;
    double hnext, hstart, tend, tgot, tol, tstart, tnext, waste;
    double *rwsav = 0, *thresh = 0, *ygot = 0, *yinit = 0, *ymax = 0;
    double *ypgot = 0;
    Integer *iwsav = 0;
    NagError fail;
    Nag_RK_method method;
    Nag_ErrorAssess errass;
    Nag_Comm comm;
    std::string nag_output = "outputs/airy-nag.txt";
                                                                                  
    INIT_FAIL(fail);

    // dimensions of solution vector
    n = 4;
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
                                                                                  
    
    /* Set initial conditions for ODE and parameters for the integrator. */
    /* nag_enum_name_to_value (x04nac) Converts NAG enum member name to value. */
    method = (Nag_RK_method) nag_enum_name_to_value("Nag_RK_4_5");
    errass = (Nag_ErrorAssess) nag_enum_name_to_value("Nag_ErrorAssess_off");
    tstart = 1.0;
    tend = 100000;
    yinit[0] = 1.0;
    yinit[1] = 1.0;
    yinit[2] = boost::math::airy_ai(-tstart);
    yinit[3] = -boost::math::airy_ai_prime(-tstart);
    hstart = 0.0;
    thresh[0] = 1.0e-8;
    thresh[1] = 1.0e-8;
    thresh[2] = 1.0e-8;
    thresh[3] = 1.0e-8;
    tol = 1.0e-5;
                                                                                  
    nag_ode_ivp_rkts_setup(n, tstart, tend, yinit, tol, thresh, method,
                             errass, hstart, iwsav, rwsav, &fail);
          
    tnext = tstart;
    while (tnext < tend){
        nag_ode_ivp_rkts_range(f, n, tend, &tgot, ygot, ypgot, ymax, &comm,
                               iwsav, rwsav, &fail);
                                                                                  
        tnext = tgot;
        // write current step to file
        std::ofstream fout;
        fout.open(nag_output, std::ios_base::app);
        fout << tgot << " ";
        for (int k = 0; k < n; k++)
            fout <<  ygot[k] << " ";
        fout << std::endl;
        fout.close();
    };
    
    NAG_FREE(thresh);
    NAG_FREE(yinit);
    NAG_FREE(ygot);
    NAG_FREE(ypgot);
    NAG_FREE(ymax);
    NAG_FREE(rwsav);
    NAG_FREE(iwsav);
};

