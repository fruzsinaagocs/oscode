#pragma once
#include <iterator>
#include <cmath>

// TODO: size checking
template<typename X, typename Y>
struct LinearInterpolator{

    double x0, dx;
//    int Nx, Ny;
    X x_; // array of indep. variable
    Y y_; // array of dep. variable

    LinearInterpolator(const X &x, const Y &y){
        // Constructor of struct, sets struct members
        x_ = x;
        y_ = y;
        x0 = x[0];
        dx = x[1]-x[0];
//        Nx = x_.size();
//        Ny = y_.size();
//        if(Nx==0 or Ny==0){
//            std::cout << "One or more of the supplied arrays has size 0" <<
//            std::endl;
//        }
//        if(Nx!=Ny){
//            std::cout << "Sizes of supplied arrays do not match" << std::endl;
//        }
//        std::cout << "size: " << Nx << ", " << Ny << std::endl;
    }
    
    std::complex<double> operator() (double x) const {
        // Does linear interpolation 
        int i=int((x-x0)/dx);
        std::complex<double> y0=y_[i];
        std::complex<double> y1=y_[i+1];
        return (y0+(y1-y0)*(x-x0-dx*i)/dx);
    }

    std::complex<double> expit (double x) const {
        // Does linear interpolation when the input is ln()-d
        int i=int((x-x0)/dx);
        std::complex<double> y0=y_[i];
        std::complex<double> y1=y_[i+1];
        return std::exp(y0+(y1-y0)*(x-x0-dx*i)/dx);
    }
   
};

