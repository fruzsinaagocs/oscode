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
    
    //de_system sys(&w, &g);
    
    //std::vector<double> times(10000,0);
    //std::vector<std::complex<double>> ws(10000,0), gs(10000,0);
    //for(int j=0; j<10000; j++){
    //    times[j] = j/100.0;
    //    ws[j] = w(j/100.0); 
    //}
    //de_system sys(times, ws, gs);
    
    //int N = 10000; 
    //Eigen::VectorXd times=Eigen::VectorXd::Zero(N);
    //Eigen::VectorXcd ws = Eigen::VectorXcd::Zero(N),  gs = Eigen::VectorXcd::Zero(N);
    //for(int j=0; j<N; j++){
    //    times(j) = j/100.0;
    //    ws(j) = w(j/100.0);
    //}

    int N = 10001;
    double times_arr[10001];
    std::complex<double> ws_arr[10001], gs_arr[10001];
    
    for(int j=0; j<N; j++){
        times_arr[j] = j/100.0;
        ws_arr[j] = w(j/100.0);
        gs_arr[j] = 0.0;
    }
    std::cout << times_arr[N-1] << std::endl;
    double * times;
    std::complex<double> * ws, * gs;
    times = times_arr;
    ws = ws_arr;
    gs = gs_arr;

    de_system sys(times, ws, gs, times, N);
    std::cout << "Built system" << std::endl;
    
    // Define integration range
    ti = 0.0;
    tf = n;
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
    std::cout << "Done! " << std::endl;

};


