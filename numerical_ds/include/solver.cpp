#include "interpolator.hpp"
#include "solver.hpp"
//#include <vector>
#include <cmath>
#include <chrono>
#include <boost/math/special_functions/airy.hpp>
#include <fstream>
#include <string>

double n = 1e5;

std::complex<double> g0(double t){
    return 0.0;
};
    
std::complex<double> wairy(double t){
    return std::pow(t, 0.5);
};

std::complex<double> wburst(double t){
    return std::pow(n*n - 1.0,0.5)/(1.0 + t*t);
};

std::complex<double> xburst(double t){
    return 100*std::pow(1.0 + t*t,
    0.5)/n*std::complex<double>(std::cos(n*std::atan(t)),std::sin(n*std::atan(t))); 
};

std::complex<double> dxburst(double t){
    return 100/std::pow(1.0 + t*t,
    0.5)/n*(std::complex<double>(t,n)*std::cos(n*std::atan(t)) +
    std::complex<double>(-n,t)*std::sin(n*std::atan(t))); 
};

std::complex<double> xairy(double t){
    return std::complex<double>(boost::math::airy_ai(-t), boost::math::airy_bi(-t));
};

std::complex<double> dxairy(double t){
    return std::complex<double>(-boost::math::airy_ai_prime(-t), -boost::math::airy_bi_prime(-t));
};

int main(){
   
    // Example with w(t), g(t) analytically given
    std::ofstream f;
    int no = 400;
    Eigen::VectorXd rtols = Eigen::VectorXd::LinSpaced(no,-3.0,-7.0);
    std::vector<int> steps,totsteps,wkbsteps;
    std::vector<double> runtimes;
    std::complex<double> x0, dx0;
    double ti, tf, rtol, atol, h0, runtime;
    bool full_output = false;
    int order = 3;   
    for(int i=0; i<no; i++){
        std::cout << "n: " << n << std::endl;
        de_system sys(&wburst, &g0);
        ti = -2*n;
        tf = 2*n;
        x0 = xburst(ti); 
        dx0 = dxburst(ti); 
        rtol = std::pow(10, rtols(i));
        atol = 0.0;
        h0 = 1.0;
        runtime = 0.0;
        for(int j=0; j<100; j++){
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
    
    f.open("plots/bursttimingn1e5.txt");
    f << "# Testing how tolerance affects runtime in the burst equation\n" <<
    "# n = " << n << "\n" << 
    "# rtol, runtime/ms, total steps, successful steps, wkb steps" << std::endl;
    for(int i=0; i<no; i++){
        f << std::setprecision(20) << rtols(i) << ", " << runtimes[i] << ", " << totsteps[i] << ", " << steps[i] << ", " << wkbsteps[i] << std::endl;  
    };
    f.close();
    
    // Example of the same solution, but with t,w,g supplied as a grid
//    int n = 100000;
//    Eigen::VectorXd logts = Eigen::VectorXd::LinSpaced(n, 0.0, 5.2);
//    Eigen::VectorXd ts = logts;
//    Eigen::VectorXcd ws = Eigen::VectorXcd::Zero(n);
//    Eigen::VectorXcd gs = Eigen::VectorXcd::Zero(n); 
//    for(int i=0; i<n; i++){
//        ts(i) = std::pow(10,ts(i));
//        ws(i) = std::pow(ts(i), 0.5);
//    };
//    de_system sys2(ts, ws, gs);
//    Solution solution2(sys2, x0, dx0, ti, tf);
//    auto t3 = std::chrono::high_resolution_clock::now();
//    solution2.solve();
//    auto t4 = std::chrono::high_resolution_clock::now();
//    std::chrono::duration<double,std::milli> t34 = t4-t3;
//    std::cout << "time: " << t34.count() << " ms." << std::endl;

};


