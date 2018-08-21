#include "test.hpp"
using namespace RKWKB;

Vector F(Vector z){
    // A trial function to describe the background evolution. 

    Vector result(z.size());
    result << 0.25*std::pow(z(0),-3), 0.25*std::pow(z(1),-3);
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
    de_system my_system(F);
    Solution my_solution(my_system, y, 1.0, f_end); 
    my_solution.evolve();
    return 0;
};

 
