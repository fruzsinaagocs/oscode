#include "interpolator.hpp"
#include "solver.hpp"
//#include <vector>
#include <cmath>
#include <chrono>
#include <boost/math/special_functions/airy.hpp>
#include <fstream>
#include <string>
#include <stdlib.h>

double m=1, E=1e2, K=1;

double RKSolver::g(double t){
    return 0.0;
};

double gdummy(double t){
    return t;
};

double wdummy(double t){
    return t;
};

double V(double t){
// quantum mechanical potential well
    return 0.5*K*t*t;
};

double RKSolver::w(double t){
    return std::sqrt(2*m*(E-V(t)));
};

int main(){
     
    // Example with w(t), g(t) analytically given
    std::ofstream f;
    int no = 1;
    Eigen::VectorXd ns = Eigen::VectorXd::LinSpaced(no,8.0,8.0);
    std::vector<int> steps,totsteps,wkbsteps;
    std::vector<double> runtimes;
    std::complex<double> x0, dx0;
    double ti, tf, rtol, atol, h0, runtime;
    bool full_output = true;
    int order = 3;   
    for(int i=0; i<no; i++){
        //n = round(std::pow(10,ns(i)));
        //std::cout << "n: " << n << std::endl;
        de_system sys(&wdummy, &gdummy);
        ti = -10.0;
        tf = 10.0;
        x0 = 0.0; 
        dx0 = 1.0; 
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
    
    //f.open("plots/bursttimingtol-6_stepscorr.txt");
    //f << "# Testing how tolerance affects runtime in the burst equation\n" <<
    //"# tolerance rtol = " << rtol << "\n" << 
    //"# log10(n), total steps, successful steps, wkb steps" << std::endl;
    //for(int i=0; i<no; i++){
    //    f << std::setprecision(20) << ns(i) << ", " << totsteps[i] << ", " << steps[i] << ", " << wkbsteps[i] << std::endl;  
    //};
    //f.close();
    

};


