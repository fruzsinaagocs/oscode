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
    // TODO: any way not to create all objects here?
    //WKBSolver wkbsolver
    WKBSolver * wkbsolver;
    WKBSolver1 wkbsolver1;
    WKBSolver2 wkbsolver2;
    WKBSolver3 wkbsolver3;

    public:
    // constructor
    Solution(de_system de_sys, std::complex<double> x0, std::complex<double> dx0, double t_i, double t_f, int o=3, double r_tol=1e-4, double a_tol=0.0, double h_0=1, bool full_output=false);
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
    switch(order){
        case 1: wkbsolver1 = WKBSolver1(de_sys);
                wkbsolver = &wkbsolver1;
                break;
        case 2: wkbsolver2 = WKBSolver2(de_sys);
                wkbsolver = &wkbsolver2;
                break;
        case 3: wkbsolver3 = WKBSolver3(de_sys);
                wkbsolver = &wkbsolver3;
                break;
    };
};

void Solution::solve(){ 
    
    std::cout << "Function w at t=" << t << " is: " << w(t) << ", the order is " << order << std::endl;
    Eigen::Matrix<std::complex<double>,1,4> y;
    y << x, dx, 0.0, 0.0;
    std::cout << "Function f at t=" << t << " is: " << rksolver.f(t,y) << std::endl; 
    std::cout << "Result of a RK step: " << rksolver.step(x, dx, t, h0) << std::endl;
    Eigen::Matrix<std::complex<double>,2,2> y1;
    //std::cout << "ws rk: " << rksolver.ws << std::endl;
    y1 = wkbsolver->step(x, dx, t, h0, rksolver.ws, rksolver.gs, rksolver.ws5, rksolver.gs5);
    std::cout << "WKB step: " << y1 << std::endl;


};
