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
    Solution my_solution(my_system, y, 1.0, 10.0);
    Step my_step = my_solution.step(&my_solution.rkfsolver, my_solution.y, my_solution.h);
    std::cout << "result of first step: " <<  my_step.y << ", and its error: " << my_step.error << ". Was this a WKB step? " << std::boolalpha << my_step.wkb << std::endl;

}


