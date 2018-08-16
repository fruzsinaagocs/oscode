#include "test.hpp"
using namespace RKWKB;

Vector F(Vector y){
    // A trial function to describe the background evolution. 

    Vector result(y.size());
    result << 0.25*std::pow(y(0),-3), 0.25*std::pow(y(1),-3);
    return result;

}


int main(){

    Vector y(2);
    y << 2.0, 3.0;
    de_system my_system(F);
    RKFsolver my_solver(my_system);
    my_solver.step(my_solver.sys.F, y, 1.0);
    //std::cout << "F(y): " << my_solver.sys.F(y) << std::endl;

}


