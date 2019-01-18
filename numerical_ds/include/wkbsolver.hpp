#pragma once
#include "system.hpp"

class WKBSolver
{
    private: 
    // frequency and friction terms 
    std::complex<double> (*w)(double);
    std::complex<double> (*g)(double);   
    // declarations of derivative terms 
    std::complex<double> d1w1(double);
    std::complex<double> d1w2(double);
    std::complex<double> d1w3(double);
    std::complex<double> d1w4(double);
    std::complex<double> d1w5(double);
    std::complex<double> d1w6(double);
    std::complex<double> d2w1(double);
    std::complex<double> d2w6(double);
    std::complex<double> d3w1(double);
    std::complex<double> d3w6(double);
    std::complex<double> d4w1(double);
    std::complex<double> d1g1(double);
    std::complex<double> d1g6(double);
    std::complex<double> d2g1(double);
    std::complex<double> d2g6(double);
    std::complex<double> d3g1(double);
    std::complex<double> d1w2_5(double);
    std::complex<double> d1w3_5(double);
    std::complex<double> d1w4_5(double);

    public:
    // constructor
    WKBSolver();
    WKBSolver(de_system);
    //void step(); 

};

WKBSolver::WKBSolver(){
};

WKBSolver::WKBSolver(de_system de_sys){
     
    // Set frequency and friction terms
    w = de_sys.w;
    g = de_sys.g;     


};
