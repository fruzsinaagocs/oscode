#include "interpolator.hpp"
#include "solver.hpp"
//#include <vector>
#include <cmath>
#include <chrono>
#include <boost/math/special_functions/hermite.hpp>
#include <fstream>
#include <string>
#include <stdlib.h>

int n=2;
double m=1, K=2, E=std::sqrt(K/m)*(n-0.5);

std::complex<double> RKSolver::g(double t){
    return 0.0;
};

std::complex<double> gdummy(double t){
    return t;
};

std::complex<double> wdummy(double t){
    return t;
};

double V(double t){
// quantum mechanical potential well
    return 0.5*K*t*t;
};

std::complex<double> RKSolver::w(double t){
    std::complex<double> result = std::sqrt(std::complex<double>(2*m*(E-V(t))));
    return result;
};

Eigen::Vector2cd ic(double t){

    double gam = std::sqrt(m*K);
    Eigen::Vector2cd result;
    double A = std::pow(gam/M_PI,1/4.0)*1.0/std::sqrt(std::pow(2.0,n-1)*tgamma(n));
    result << boost::math::hermite(n-1,std::sqrt(gam)*t)*std::exp(-0.5*gam*t*t), -gam*t*boost::math::hermite(n-1,std::sqrt(gam)*t)*std::exp(-0.5*gam*t*t) + 2*(n-1)*boost::math::hermite(n-2,std::sqrt(gam)*t)*std::sqrt(gam)*std::exp(-0.5*gam*t*t);
    return A*result;
};

int main(){
     
    // Example with w(t), g(t) analytically given
    std::ofstream f;
    int no = 1;
    std::vector<int> steps,totsteps,wkbsteps;
    std::vector<double> runtimes;
    std::complex<double> x0, dx0;
    double ti, tf, rtol, atol, h0, runtime;
    bool full_output = true;
    int order = 3;   
    for(int i=0; i<no; i++){
        de_system sys(&wdummy, &gdummy);
        ti = (-std::sqrt((n-0.5)*2/std::sqrt(K*m))-0.5);
        tf = -ti;
        Eigen::Vector2cd ics=ic(ti);
        x0 = ics(0);
        dx0 = ics(1); 
        rtol = 1e-4;
        atol = 0.0;
        h0 = 1.0;
        runtime = 0.0;
        for(int j=0; j<1; j++){
            Solution solution(sys, x0, dx0, ti, tf, order, rtol, atol, h0, full_output); 
            auto t1 = std::chrono::high_resolution_clock::now();
            solution.solve();
            auto t2 = std::chrono::high_resolution_clock::now();
            std::chrono::duration<double,std::milli> t12 = t2-t1;
            runtime += t12.count();
        };
        std::cout << "time: " << runtime/100.0 << " ms." << std::endl;
        
        Solution solution(sys, x0, dx0, ti, tf, order, rtol, atol, h0, full_output);
        solution.solve();
        steps.emplace_back(solution.ssteps);
        wkbsteps.emplace_back(solution.wkbsteps);
        totsteps.emplace_back(solution.totsteps);
        runtimes.emplace_back(runtime/100.0);
        std::cout << "steps: " << solution.totsteps << std::endl;
    };
    
};


