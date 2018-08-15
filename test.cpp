#include "test.hpp"
//#include <cstdio> 

V F(V y){
    // A trial function to describe the background evolution. 

    V result(y.size());
    result << 0.25*std::pow(y(0),-3), 0.25*std::pow(y(1),-3);
    return result;

}


int main(){

    V y(2);
    y << 2.0, 3.0;
    de_system my_system(F);
    WKBsolver my_solver(my_system);
    std::cout << "F(y): " << my_solver.sys.F(y) << std::endl;

}


