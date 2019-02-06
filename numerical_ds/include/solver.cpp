#include "interpolator.hpp"
#include "solver.hpp"
//#include <vector>
#include <cmath>
#include <chrono>
#include <boost/math/special_functions/airy.hpp>

double n = 1e9;

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
    // Solving the Airy equation from t=1 to t=10^6
    de_system sys(&wburst, &g0);
    // solution takes: system, x0, dx0, ti, tf, (order, rtol, atol, h, full_output)
    std::complex<double> x0, dx0;
    double ti, tf;
    ti = -2*n;
    tf = 2*n;
    x0 = xairy(ti); 
    dx0 = dxairy(ti); 
    Solution solution(sys, x0, dx0, ti, tf, 3, 1e-4, 0.0, 1.0, true); 
    auto t1 = std::chrono::high_resolution_clock::now();
    solution.solve();
    auto t2 = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double,std::milli> t12 = t2-t1;
    std::cout << "time: " << t12.count() << " ms." << std::endl;
    std::cout << "\n\n\n" << std::endl;

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


