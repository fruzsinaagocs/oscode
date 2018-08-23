#include <iostream>
#include "typedefs.hpp"

namespace RKWKB
{
    class de_system
    {
        private:
    
        public:
            // default constructor
            de_system();        
            // constructor overloads
            de_system(Vectorfn, Matrixfn, Scalarfn, Vectorfn, Matrixfn, Scalarfn, Vectorfn, Matrixfn, Rank3Tfn=std::function<Rank3T(Vector)>{}, Rank3Tfn=std::function<Rank3T(Vector)>{}, Rank3Tfn=std::function<Rank3T(Vector)>{}, Rank4Tfn=std::function<Rank4T(Vector)>{}, Rank4Tfn=std::function<Rank4T(Vector)>{}, Rank4Tfn=std::function<Rank4T(Vector)>{});
            
            // class data
            
            // class functions
            Vectorfn F, Dw, Dg;
            Matrixfn DF, DDw, DDg;
            Scalarfn w, g, dw, dg, ddw, ddg;
            // optional for higher order solvers (if not given but called, will
            // give bad_function_call
            Rank3Tfn D3w, D3g, DDF;
            Rank4Tfn D4w, D4g, D3F;
            Scalarfn d3w, d4w, d3g, d4g;
    };
    
    de_system::de_system(){
        // Default constructor of de_system (does nothing)
    };
    
    de_system::de_system(Vectorfn f, Matrixfn Df, Scalarfn W, Vectorfn DW, Matrixfn DDW, Scalarfn gam, Vectorfn Dgam, Matrixfn DDgam, Rank3Tfn D3W, Rank3Tfn D3gam, Rank3Tfn DDf, Rank4Tfn D4W, Rank4Tfn D4gam, Rank4Tfn D3f){
        // default constructor for a system of differential equations
    
        F = f;
        DF = Df;
        w = W;
        Dw = DW;
        DDw = DDW;
        g = gam;
        Dg = Dgam;
        DDg = DDgam;
        D3w = D3W;
        D3g = D3gam;
        DDF = DDf;
        D4w = D4W;
        D4g = D4gam;
        D3F = D3f;
    };
    
    Scalar dw(Vector y){
        return 1.0;
    };

    Scalar dg(Vector y){
        return 1.0;
    };

}




