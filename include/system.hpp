#pragma once
#include "interpolator.hpp"

class de_system
    {
    private:
        // interpolation checker
        int even_;

    public:
        // constructors
        template<typename X, typename Y, typename Z, typename X_it>de_system(X &ts, Y &ws, Z &gs, X_it x_it, int size, bool isglogw=false, bool islogg=false, int even=0, int check_grid=0);
        de_system(std::complex<double> (*w)(double), std::complex<double> (*g)(double));
        de_system();
        std::function<std::complex<double>(double)> w;
        std::function<std::complex<double>(double)> g;
        LinearInterpolator<> Winterp;
        LinearInterpolator<> Ginterp; 
        bool islogg_, islogw_;
        // Grid fineness check
        bool grid_fine_enough = 1;
};

de_system::de_system(){
    // Default constructor
}

template<typename X, typename Y, typename Z, typename X_it> de_system::de_system(X &ts, Y &ws, Z &gs, X_it x_it, int size, bool islogw, bool islogg, int even, int check_grid){
    
    // 3. User supplied t, w/ln(w), g/ln(g) as array-like objects
    // (Eigen::Vectors, std::vectors, or arrays)
    
    even_ = even;
    islogg_ = islogg;
    islogw_ = islogw;

    LinearInterpolator<X,Y,X_it> winterp(ts,ws,even_);
    LinearInterpolator<X,Z,X_it> ginterp(ts,gs,even_);
    
    Winterp = winterp;
    Ginterp = ginterp;

    if(even_==0){
        Winterp.set_interp_start(x_it);
        Ginterp.set_interp_start(x_it);
        Winterp.set_interp_bounds(ts,ts+size-1);
        Ginterp.set_interp_bounds(ts,ts+size-1);
    }
    
    // Grid fineness check
    if(check_grid == 1){
        int w_is_fine = Winterp.check_grid_fineness(size);
        int g_is_fine = Ginterp.check_grid_fineness(size);
        if(w_is_fine==1 && g_is_fine==1)
            grid_fine_enough = 1; 
        else
            grid_fine_enough = 0;
    }

    if(islogw)
        w = std::bind(&LinearInterpolator<X,Y,X_it>::expit, Winterp,
        std::placeholders::_1);
    else
        w = std::bind(&LinearInterpolator<X,Y,X_it>::operator(), Winterp,
        std::placeholders::_1);
    if(islogg)
        g = std::bind(&LinearInterpolator<X,Z,X_it>::expit, Ginterp,
        std::placeholders::_1);
    else
        g = std::bind(&LinearInterpolator<X,Z,X_it>::operator(), Ginterp,
        std::placeholders::_1);
     
      
}

de_system::de_system(std::complex<double> (*W)(double), std::complex<double> (*G)(double)){
    // Default constructor for a system of differential equations
    // 1. User (in c++) defined the functions w,g themselves and supplied
    // function pointers.

    w = W;
    g = G;
};

