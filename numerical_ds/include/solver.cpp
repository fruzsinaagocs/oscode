#include "solver.hpp"
    
std::complex<double> g(double t){
    return t*t*t;
};
    
std::complex<double> w(double t){
    std::cout << "w called" << std::endl;
    return t*t;
};

int main(){
    
    de_system sys(&w, &g);
    Solution solution(sys, 1.0, 3.0, 3.0, 100.0); 
    solution.solve();

};


