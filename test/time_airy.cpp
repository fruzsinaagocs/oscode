#include "test-airy.hpp"
#include <time.h>

Vector background(Scalar t, int d){
    // Gives background at a given time t
    Vector result = Vector::Ones(d);
    return result*std::pow(t, 1.0/4.0);
};

Vector airy(double t){
    Vector result(2);
    result << std::complex<double>(boost::math::airy_ai(-t), 0.0), std::complex<double>(-boost::math::airy_ai_prime(-t), 0.0);
    //result << std::complex<double>(boost::math::airy_ai(-t), boost::math::airy_bi(-t)), std::complex<double>(-boost::math::airy_ai_prime(-t), -boost::math::airy_bi_prime(-t));
    return result; 
}


int main(){

    double t = 1.0;
    Vector ic(4);                                      
    ic << airy(t), background(t, 2);
    clock_t begin = clock();
    for(int i=0; i<1000; i++){
        de_system my_system(F, DF, w, Dw, DDw, g, Dg, DDg);   
        Solution my_solution(my_system, ic, t, f_end);
        my_solution.evolve();
        //std::cout << "solution at " << my_solution.t << ": " << my_solution.y(0) << "," << my_solution.y(1) << std::endl;
        //std::cout << "total steps: " << my_solution.stepsall << std::endl;
        //std::cout << "waste ratio: " << my_solution.waste << std::endl;
        //std::cout << "accuracy at end of integration: " << abs((my_solution.y(0) - airy(my_solution.t)(0))/airy(my_solution.t)(0)) << "," << abs((my_solution.y(1) - airy(my_solution.t)(1))/airy(my_solution.t)(1)) << std::endl;
    }
    clock_t end = clock();
    double time_spent = (double)(end-begin)/(CLOCKS_PER_SEC*1e3);
    std::cout << "elapsed time per loop: " << time_spent << std::endl;
}
