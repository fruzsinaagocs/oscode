#include "test-airy.hpp"

TEST_CASE("setting up a system of ODEs"){
    
    Scalar t = 6.0;
    Vector y(2);                                        
    y << std::pow(t, 1/4.0), std::pow(t, 1/4.0);
    de_system my_system(F, DF, w, Dw, DDw, g, Dg, DDg);
     
    REQUIRE( std::abs(my_system.w(y) - std::pow(t, 1/2.0)) <= 1e-14 );
    REQUIRE( std::abs(my_system.dw(y) - 1/2.0*std::pow(t, -1/2.0)) <= 1e-14 );
    REQUIRE( std::abs(my_system.ddw(y) - -1/4.0*std::pow(t, -3/2.0)) <= 1e-14 );
    REQUIRE( std::abs(my_system.dg(y) - 0.0) <= 1e-14 );
    REQUIRE( std::abs(my_system.ddg(y)- 0.0) <= 1e-14 );

    Solution my_solution(my_system, y, 1.0, f_end);
    Vector y_tot0(5), y_tot1(5);
    y_tot0 << std::complex<double>(0.2, 1.0), std::complex<double>(3e-5, -2.0), y, 0.0;
    y_tot1 << 3.0, 4.0, 4.07, 4.07, std::complex<double>(0,3e-4);
    my_solution.wkbsolver1.y0 = y_tot0;
    my_solution.wkbsolver1.y1 = y_tot1;
    Step wkb_step = my_solution.wkbsolver1.step(F, y, 1.0);
    //std::cout << "wkb step: " << wkb_step.y << " " << wkb_step.error << " " << wkb_step.wkb << std::endl;
    Vector my_dS = my_solution.wkbsolver1.dS(y);
    Vector my_S_odd = my_solution.wkbsolver1.S_odd(y);
    
    Vector analytic_dS(2);
    analytic_dS << std::complex<double>(0.0, 1.0)*std::pow(t, 1/2.0), -1.0/(4.0*t);
    REQUIRE( analytic_dS.size() == my_dS.size());
    for(int i=0; i<analytic_dS.size(); i++)
            REQUIRE( std::abs(my_dS(i) - analytic_dS(i)) <= 1e-14 );

}
