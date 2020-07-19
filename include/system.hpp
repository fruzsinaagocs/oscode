#pragma once
#include "interpolator.hpp"

class de_system
    {
    private:
        // interpolation checker

    public:
        // constructors
        template<typename X, typename Y, typename Z, typename X_it> de_system(X &ts, Y &ws, Z &gs, X_it x_it, int size, bool isglogw=false, bool islogg=false);
        de_system(std::complex<double> (*w)(double), std::complex<double> (*g)(double));
        std::function<std::complex<double>(double)> w;
        std::function<std::complex<double>(double)> g;
};

//void de_system::set_array_types(Eigen::VectorXd::InnerIterator it){
//}
//
//void de_system::set_array_types(std::vector<double>::iterator it){
//}
//
//void de_system::set_array_types(double * it){
//}



template<typename X, typename Y, typename Z, typename X_it> de_system::de_system(X &ts, Y &ws, Z &gs, X_it x_it, int size, bool islogw, bool islogg){
    
    // 3. User supplied t, w/ln(w), g/ln(g) as array-like objects
    // (Eigen::Vectors, std::vectors, or arrays)
    
    int even = 0;

    LinearInterpolator<X, Y, X_it> winterp(ts,ws,even);
    LinearInterpolator<X, Z, X_it> ginterp(ts,gs,even);
    
    if(even==0){
        winterp.set_interp_start(x_it);
        ginterp.set_interp_start(x_it);
        winterp.set_interp_bounds(ts,ts+size-1);
        ginterp.set_interp_bounds(ts,ts+size-1);
    }
    
    if(islogw)
        w = std::bind(&LinearInterpolator<X,Y,X>::expit, winterp,
        std::placeholders::_1);
    else
        w = std::bind(&LinearInterpolator<X,Y,X>::operator(), winterp,
        std::placeholders::_1);
    if(islogg)
        g = std::bind(&LinearInterpolator<X,Z,X>::expit, ginterp,
        std::placeholders::_1);
    else
        g = std::bind(&LinearInterpolator<X,Z,X>::operator(), ginterp,
        std::placeholders::_1);
}

de_system::de_system(std::complex<double> (*W)(double), std::complex<double> (*G)(double)){
    // Default constructor for a system of differential equations
    // 1. User (in c++) defined the functions w,g themselves and supplied
    // function pointers.

    w = W;
    g = G;
};

