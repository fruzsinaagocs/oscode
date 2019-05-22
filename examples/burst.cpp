#include "solver.hpp"
#include <cmath>
#include <fstream>
#include <string>
#include <stdlib.h>

double n = 100.0;

// Define the gamma term
std::complex<double> g(double t){
    return 0.0;
};

// Define the frequency
std::complex<double> w(double t){
    return std::pow(n*n - 1.0,0.5)/(1.0 + t*t);
};

// Initial conditions x, dx
std::complex<double> xburst(double t){
    return 100*std::pow(1.0 + t*t,
    0.5)/n*std::complex<double>(std::cos(n*std::atan(t)),std::sin(n*std::atan(t))); 
};

std::complex<double> dxburst(double t){
    return 100/std::pow(1.0 + t*t,
    0.5)/n*(std::complex<double>(t,n)*std::cos(n*std::atan(t)) +
    std::complex<double>(-n,t)*std::sin(n*std::atan(t))); 
};

int main(){

    std::ofstream f;
    std::string output = "output.txt";
    std::complex<double> x0, dx0;
    double ti, tf;
    // Create differential equation 'system'
    de_system sys(&w, &g);
    // Define integration range
    ti = -2*n;
    tf = 2*n;
    // Define initial conditions
    x0 = xburst(ti); 
    dx0 = dxburst(ti); 
    // Solve the equation
    Solution solution(sys, x0, dx0, ti, tf); 
    solution.solve();
    // The solution is stored in lists, copy the solution
    std::list<std::complex<double>> xs = solution.sol;
    std::list<double> ts = solution.times;
    std::list<bool> types = solution.wkbs;
    int steps = solution.ssteps;
    // Write result in file
    f.open(output);
    auto it_t = ts.begin();
    auto it_x = xs.begin();
    auto it_ty = types.begin();
    for(int i=0; i<=steps; i++){
        f << *it_t << ", " << std::real(*it_x) << ", " << *it_ty << std::endl; 
        ++it_t;
        ++it_x;
        ++it_ty;
    };
    f.close();

};


