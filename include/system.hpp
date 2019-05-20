#pragma once
#include "interpolator.hpp"

// TODO: size checking
class de_system
    {
    private:
        // interpolation checker

    public:
        // constructors
        template<typename X, typename Y, typename Z> de_system(const X &ts, const Y &ws, const Z &gs, bool isglogw=false, bool islogg=false);
        de_system(std::complex<double> (*w)(double), std::complex<double> (*g)(double));
        std::function<std::complex<double>(double)> w;
        std::function<std::complex<double>(double)> g;
        bool interp;
};

template<typename X, typename Y, typename Z> de_system::de_system(const X &ts, const Y &ws, const Z &gs, bool islogw, bool islogg){
    
    // 3. User supplied t, w/ln(w), g/ln(g) as array-like objects
    // (Eigen::Vectors, std::vectors, or arrays)
    
    interp = true;
    LinearInterpolator<X, Y> winterp(ts,ws);
    LinearInterpolator<X, Z> ginterp(ts,gs);
    if(islogw)
        w = std::bind(&LinearInterpolator<X,Y>::expit, winterp,
        std::placeholders::_1);
    else
        w = std::bind(&LinearInterpolator<X,Y>::operator(), winterp,
        std::placeholders::_1);
    if(islogg)
        g = std::bind(&LinearInterpolator<X,Z>::expit, ginterp,
        std::placeholders::_1);
    else
        g = std::bind(&LinearInterpolator<X,Z>::operator(), ginterp,
        std::placeholders::_1);
}

de_system::de_system(std::complex<double> (*W)(double), std::complex<double> (*G)(double)){
    // Default constructor for a system of differential equations
    // 1. User (in c++) defined the functions w,g themselves and supplied
    // function pointers.

    interp = false;
    w = W;
    g = G;
};

