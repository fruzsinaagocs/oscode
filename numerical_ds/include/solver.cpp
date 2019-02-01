#include "interpolator.hpp"
#include "solver.hpp"
//#include <vector>
#include <chrono>
#include <boost/math/special_functions/airy.hpp>

std::complex<double> g(double t){
    return 0.0;
};
    
std::complex<double> w(double t){
    return std::pow(t, 0.5);
};

int main(){
   
    // Example with w(t), g(t) analytically given
    // Solving the Airy equation from t=1 to t=10^6
    de_system sys(&w, &g);
    // solution takes: system, x0, dx0, ti, tf, (order, rtol, atol, h, full_output)
    std::complex<double> x0, dx0;
    double ti, tf;
    ti = 1.0;
    tf = 1e7;
    x0 = std::complex<double>(boost::math::airy_ai(-ti), boost::math::airy_bi(-ti));
    dx0 = std::complex<double>(-boost::math::airy_ai_prime(-ti), -boost::math::airy_bi_prime(-ti));
    Solution solution(sys, x0, dx0, ti, tf); 
    auto t1 = std::chrono::high_resolution_clock::now();
    solution.solve();
    auto t2 = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double,std::milli> t12 = t2-t1;
    std::cout << "time: " << t12.count() << " ms." << std::endl;
    std::cout << "\n\n\n" << std::endl;

    // Example of the same solution, but with t,w,g supplied as a grid
    int n = 100000;
    Eigen::VectorXd logts = Eigen::VectorXd::LinSpaced(n, 0.0, 3.0);
    Eigen::VectorXd ts = logts;
    Eigen::VectorXcd ws = Eigen::VectorXcd::Zero(n);
    Eigen::VectorXcd gs = Eigen::VectorXcd::Zero(n); 
    for(int i=0; i<n; i++){
        ts(i) = std::pow(10,ts(i));
        ws(i) = std::pow(ts(i), 0.5);
    };
    auto t3 = std::chrono::high_resolution_clock::now();
    de_system sys2(ts, ws, gs);
    Solution solution2(sys2, x0, dx0, ti, 500.0);
    solution2.solve();
    auto t4 = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double,std::milli> t34 = t4-t3;
    std::cout << "time: " << t34.count() << " ms." << std::endl;

};


