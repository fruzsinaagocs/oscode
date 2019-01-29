#include "solver.hpp"
    
std::complex<double> g(double t){
    return 0.0;//std::pow(t, 0.2);
};
    
std::complex<double> w(double t){
    return std::pow(t, 0.5);
};

int main(){
    
    de_system sys(&w, &g);
    // solution takes: system, x0, dx0, ti, tf, (order, rtol, atol, h, full_output)
    //for(int i=0; i<6000; i++){
    std::complex<double> x, dx;
    x = std::complex<double>(0.53556089, 0.10399739);
    dx = std::complex<double>(0.01016057, -0.59237562);
    Solution solution(sys, x, dx, 1.0, 50000.0); 
    solution.solve();
    //};


};


