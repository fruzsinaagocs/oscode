#pragma once
#include <Eigen/Dense>
#include <fstream>
#include <iostream>
#include "system.hpp"
#include "rksolver.hpp"
#include "wkbsolver.hpp"

class Solution
{
    private:
    // Parameters for solver
    std::complex<double> (*w) (double);
    std::complex<double> (*g) (double);
    double t, tf, rtol, atol, h0;
    std::complex<double> x, dx;
    int order;
    bool fo;
    RKSolver rksolver;


    public:
    // constructor
    Solution(de_system de_sys, std::complex<double> x0, std::complex<double> dx0, double t_i, double t_f, int o=4, double r_tol=1e-4, double a_tol=0.0, double h_0=1, bool full_output=false);
    void solve();
};


Solution::Solution(de_system de_sys, std::complex<double> x0, std::complex<double> dx0, double t_i, double t_f, int o, double r_tol, double a_tol, double h_0, bool full_output){
    
    // Set parameters for solver
    w = de_sys.w;
    g = de_sys.g;
    x = x0;
    dx = dx0;
    t = t_i;
    tf = t_f;
    order = o;
    rtol = r_tol;
    atol = a_tol;
    h0 = h_0;
    fo = full_output;
    rksolver = RKSolver(de_sys);

};

void Solution::solve(){ 
    
    std::cout << "Function w at t=" << t << " is: " << w(t) << ", the order is " << order << std::endl;
    Eigen::Matrix<std::complex<double>,2,1> y;
    y << x, dx;
    std::cout << "Function f at t=" << t << " is: " << rksolver.f(t,y) << std::endl; 


};
