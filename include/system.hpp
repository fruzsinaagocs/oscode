#pragma once
#include "interpolator.hpp"

/** */
class de_system
    {
    private:
        int even_;

    public:
        template<typename X, typename Y, typename Z, typename X_it>de_system(X
        &ts, Y &ws, Z &gs, X_it x_it, int size, bool isglogw=false, bool
        islogg=false, int even=0, int check_grid=0);
        de_system(std::complex<double> (*w)(double), std::complex<double> (*g)(double));
        de_system();
        std::function<std::complex<double>(double)> w;
        std::function<std::complex<double>(double)> g;
        LinearInterpolator<> Winterp;
        LinearInterpolator<> Ginterp; 
        bool islogg_, islogw_;
        bool grid_fine_enough = 1;
};

/** Default contructor */
de_system::de_system(){
}

/** Constructor for the case of the user having defined the frequency and
 * damping terms as sequences
 */
template<typename X, typename Y, typename Z, typename X_it>
de_system::de_system(X &ts, Y &ws, Z &gs, X_it x_it, int size, bool islogw, bool
islogg, int even, int check_grid){
    
    even_ = even;
    islogg_ = islogg;
    islogw_ = islogw;

    /** Set up interpolation on the supplied frequency and damping arrays */
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
    
    /** Check if supplied grids are sampled finely enough for the purposes of
     * linear interpolation
     */
    if(check_grid == 1){
        int w_is_fine = Winterp.check_grid_fineness(size);
        int g_is_fine = Ginterp.check_grid_fineness(size);
        if(w_is_fine==1 && g_is_fine==1)
            grid_fine_enough = 1; 
        else
            grid_fine_enough = 0;
    }

    /** Bind result of interpolation to a function, this will be called by the
     * routines taking RK and WKB steps
     */
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

/** Constructor for the case when the frequency and damping terms have been
 * defined as functions
 */
de_system::de_system(std::complex<double> (*W)(double), std::complex<double> (*G)(double)){

    w = W;
    g = G;
};

