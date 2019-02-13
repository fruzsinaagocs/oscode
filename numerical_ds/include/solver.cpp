#include "interpolator.hpp"
#include "solver.hpp"
//#include <vector>
#include <cmath>
#include <chrono>
#include <boost/math/special_functions/airy.hpp>
#include <fstream>
#include <string>

double n = 40.0;

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
     
    // Example of the same solution, but with t,w,g supplied as a grid
    const int no = round(9e5);
    double ti, tf, rtol, atol, h0;
    std::complex<double> x0, dx0;
    int order = 3;
    bool full_output = true;
    rtol = 1e-4;
    atol = 0.0;
    h0 = 1.0;
    ti = -2*n;
    tf = 2*n;
    x0 = xburst(ti);
    dx0 = dxburst(ti);
    Eigen::VectorXd logts = Eigen::VectorXd::LinSpaced(no, -6, std::log10(2*no));
    Eigen::VectorXd Ts(2*no+2);
    Ts << logts.colwise().reverse(), -16.0, -16.0, logts; 
    Eigen::VectorXd ts = Ts;
    Eigen::VectorXcd logws = Eigen::VectorXcd::Zero(2*no+2);
    Eigen::VectorXcd gs = Eigen::VectorXcd::Zero(2*no+2); 
    for(int i=0; i<(2*no+2); i++){
        if(i<no+1)
            ts(i) = -std::pow(10,ts(i));
        else
            ts(i) = std::pow(10,ts(i));
        logws(i) = std::log(std::pow(n*n - 1.0,0.5)/(1.0 + ts(i)*ts(i)));
    };
    de_system sys2(ts, logws, gs, true);
    Solution solution2(sys2, x0, dx0, ti, tf, order, rtol, atol, h0, full_output);
    auto t3 = std::chrono::high_resolution_clock::now();
    solution2.solve();
    auto t4 = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double,std::milli> t34 = t4-t3;
    std::cout << "time: " << t34.count() << " ms." << std::endl;
  
    // Example with w(t), g(t) analytically given
//    std::ofstream f;
//    int no = 500;
//    Eigen::VectorXd ns = Eigen::VectorXd::LinSpaced(no,1.0,10.0);
//    std::vector<int> steps,totsteps,wkbsteps;
//    std::vector<double> runtimes;
//    std::complex<double> x0, dx0;
//    double ti, tf, rtol, atol, h0, runtime;
//    bool full_output = false;
//    int order = 3;   
//    for(int i=0; i<no; i++){
//        n = round(std::pow(10,ns(i)));
//        std::cout << "n: " << n << std::endl;
//        de_system sys(&wburst, &g0);
//        ti = -2*n;
//        tf = 2*n;
//        x0 = xburst(ti); 
//        dx0 = dxburst(ti); 
//        rtol = 1e-6;
//        atol = 0.0;
//        h0 = 1.0;
//        runtime = 0.0;
//        for(int j=0; j<1; j++){
//            Solution solution(sys, x0, dx0, ti, tf, order, rtol, atol, h0, full_output); 
//            auto t1 = std::chrono::high_resolution_clock::now();
//            solution.solve();
//            auto t2 = std::chrono::high_resolution_clock::now();
//            std::chrono::duration<double,std::milli> t12 = t2-t1;
//            runtime += t12.count();
//        };
//        std::cout << "time: " << runtime/100.0 << " ms." << std::endl;
//        
//        Solution solution(sys, x0, dx0, ti, tf, order, rtol, atol, h0, full_output);
//        solution.solve();
//        steps.emplace_back(solution.ssteps);
//        wkbsteps.emplace_back(solution.wkbsteps);
//        totsteps.emplace_back(solution.totsteps);
//        runtimes.emplace_back(runtime/100.0);
//        std::cout << "steps: " << solution.totsteps << std::endl;
//    };
//    
//    f.open("plots/bursttimingtol-6_stepscorr.txt");
//    f << "# Testing how tolerance affects runtime in the burst equation\n" <<
//    "# tolerance rtol = " << rtol << "\n" << 
//    "# log10(n), total steps, successful steps, wkb steps" << std::endl;
//    for(int i=0; i<no; i++){
//        f << std::setprecision(20) << ns(i) << ", " << totsteps[i] << ", " << steps[i] << ", " << wkbsteps[i] << std::endl;  
//    };
//    f.close();
    

};


