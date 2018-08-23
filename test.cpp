#include "test.hpp" 
using namespace RKWKB;

Vector F(Vector z){
    Vector result(z.size());
    result << 0.25*std::pow(z(0),-3), 0.25*std::pow(z(1),-3);
    return result;
};

Matrix DF(Vector z){
    Matrix result(z.size(),z.size());
    result << -0.75*std::pow(z(0), -4), 0.0,
               0.0, -0.75*std::pow(z(1), -4);
    return result;

};

Scalar w(Vector z){
    return z(0)*z(1);
};

Vector Dw(Vector z){
    Vector result(z.size());
    result << z(1), z(0);
    return result;
}

Matrix DDw(Vector z){
    Matrix result(z.size(), z.size());
    result << 0.0, 1.0,
              1.0, 0.0;
    return result;
}

Scalar g(Vector z){
    return 0.0;
};

Vector Dg(Vector z){
    Vector result(z.size());
    result << 0.0, 0.0;
    return result;
};

Matrix DDg(Vector z){
    Matrix result(z.size(), z.size());
    result << 0.0, 0.0,
              0.0, 0.0;
    return result;
};

Scalar f_end(Vector y, Scalar t){
    // A function to end integration of the ODE when f_end=0.
    double t_f = 10.3;
    return (t - t_f);
};

int main(){

    Vector y(2);
    y << 1.0, 1.0;
    de_system my_system(F, DF, w, Dw, DDw, g, Dg, DDg);
    std::cout << "DE system: " << std::endl;
    std::cout << "dw: " << my_system.dw(y) << std::endl;
    std::cout << "dg: " << my_system.dg(y) << std::endl;
    std::cout << "ddw: " << my_system.ddw(y) << std::endl;
    std::cout << "ddg: " << my_system.ddg(y) << std::endl;
    Solution my_solution(my_system, y, 1.0, f_end);
    Step wkb_step = my_solution.wkbsolver1.step(F, y, 1.0);
    std::cout << "wkb step: " << wkb_step.y << " " << wkb_step.error << " " << wkb_step.wkb << std::endl;
    //Vector my_dS = my_solution.wkbsolver1.dS(y);
    //std::cout << "dS: " << my_dS << std::endl;
    return 0;
};

 
