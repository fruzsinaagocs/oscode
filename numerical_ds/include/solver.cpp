#include "interpolator.hpp"
#include "solver.hpp"
//#include <vector>
#include <cmath>
#include <chrono>
#include <boost/math/special_functions/airy.hpp>
#include <fstream>
#include <string>
#include <stdlib.h>

double n = 40.0;

std::complex<double> g(double t){
    return 0.0;
};

//std::complex<double> w(double t){
//    return std::pow(t, 0.5);
//};

std::complex<double> w(double t){
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

//de_system create_system(){
//    
//    const int no = round(1e6);
//    double ti = -2*n;
//    double tf = 2*n;
//    Eigen::VectorXd logts = Eigen::VectorXd::LinSpaced(no, -6, std::log10(tf));
//    Eigen::VectorXd Ts(2*no+2);
//    Ts << logts.colwise().reverse(), -360.0, -360.0, logts; 
//    Eigen::VectorXd ts = Ts;
//    Eigen::VectorXcd logws = Eigen::VectorXcd::Zero(2*no+2);
//    Eigen::VectorXcd gs = Eigen::VectorXcd::Zero(2*no+2); 
//    for(int i=0; i<(2*no+2); i++){
//        if(i<(no+1))
//            ts(i) = -std::pow(10,ts(i));
//        else
//            ts(i) = std::pow(10,ts(i));
//        logws(i) = std::log(std::pow(n*n - 1.0,0.5)/(1.0 + ts(i)*ts(i)));
//    };
//    de_system sys2(ts, logws, gs, true);
//    Ts.resize(0); ts.resize(0); logts.resize(0); logws.resize(0); gs.resize(0);
//    return sys2;
//};

int main(){
     
    // Example of the same solution, but with t,w,g supplied as a grid
    //const int no = round(1e6);
//    double ti, tf, rtol, atol, h0;
//    std::complex<double> x0, dx0;
//    int order = 3;
//    bool full_output = true;
//    rtol = 1e-4;
//    atol = 0.0;
//    h0 = 1.0;
//    ti = -2*n;
//    tf = 2*n;
//    x0 = xburst(ti);
//    dx0 = dxburst(ti);
//    auto t1 = std::chrono::high_resolution_clock::now();
//    de_system sys2 = create_system();
//    auto t2 = std::chrono::high_resolution_clock::now(); 
//    Solution solution2(sys2, x0, dx0, ti, tf, order, rtol, atol, h0);//, full_output);
//    auto t3 = std::chrono::high_resolution_clock::now();
//    solution2.solve();
//    auto t4 = std::chrono::high_resolution_clock::now();
//    std::chrono::duration<double,std::milli> t12 = t2-t1;
//    std::chrono::duration<double,std::milli> t23 = t3-t2;
//    std::chrono::duration<double,std::milli> t34 = t4-t3;
//    std::cout << "\nsetting up system: " << t12.count() << "\nsetting up solver: " <<
//    t23.count() << "\nsolving: " << t34.count() << std::endl;
//    auto t6 = std::chrono::high_resolution_clock::now();
//    srand (time(NULL));
//    std::complex<double> result;
//    double ttry;
//    for(int i=0; i<11*100; i++){
//        ttry = (rand() % 100)/10.0;
//        result = sys2.w(ttry);
//        result = sys2.g(ttry);
//    };
//    auto t7 = std::chrono::high_resolution_clock::now();
//    std::chrono::duration<double,std::milli> t67 = t7-t6;
//    std::cout << "Calling w, g takes: " << t67.count() << std::endl;
//    std::cout << "total steps: " << solution2.totsteps << std::endl;
//    //for(auto it=solution2.wkbs.begin(); it!=solution2.wkbs.end(); ++it)
//    //    std::cout << *it << " ";

    // Example with w(t), g(t) analytically given
    std::ofstream f;
    int no = 1;
//    Eigen::VectorXd ns = Eigen::VectorXd::LinSpaced(no,8.0,8.0);
    std::vector<int> steps,totsteps,wkbsteps;
    std::vector<double> runtimes;
    std::complex<double> x0, dx0;
    double ti, tf, rtol, atol, h0, runtime;
    bool full_output = true;
    int order = 3;   
    for(int i=0; i<no; i++){
//        n = 40;
        std::cout << "n: " << n << std::endl;
        de_system sys(&w, &g);
        ti = -2*n;
        tf = 2*n;
        x0 = xburst(ti); 
        dx0 = dxburst(ti); 
        rtol = 1e-3;
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
        
        Solution solution(sys, x0, dx0, ti, tf, order, rtol, atol, h0, false);
        solution.solve();
        steps.emplace_back(solution.ssteps);
        wkbsteps.emplace_back(solution.wkbsteps);
        totsteps.emplace_back(solution.totsteps);
        runtimes.emplace_back(runtime/100.0);
        std::cout << "steps: " << solution.totsteps << std::endl;
    };
    
    //f.open("plots/bursttimingtol-6_stepscorr.txt");
    //f << "# Testing how tolerance affects runtime in the burst equation\n" <<
    //"# tolerance rtol = " << rtol << "\n" << 
    //"# log10(n), total steps, successful steps, wkb steps" << std::endl;
    //for(int i=0; i<no; i++){
    //    f << std::setprecision(20) << ns(i) << ", " << totsteps[i] << ", " << steps[i] << ", " << wkbsteps[i] << std::endl;  
    //};
    //f.close();
    

};


