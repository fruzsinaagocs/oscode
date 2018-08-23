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
            Vectorfn F, Dw, Dg;
            Matrixfn DF, DDw, DDg;
            Scalarfn w, g;
            // optional for higher order solvers (if not given but called, will
            // give bad_function_call
            Rank3Tfn D3w, D3g, DDF;
            Rank4Tfn D4w, D4g, D3F;
            
            // class functions
            Scalar dw(Vector);
            Scalar dg(Vector);
            Scalar ddw(Vector);
            Scalar ddg(Vector);
            Scalar d3w(Vector);
            Scalar d4w(Vector);
            Scalar d3g(Vector);
            Scalar d4g(Vector);
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
    
    Scalar de_system::dw(Vector y){
        return Dw(y).transpose()*(F(y));
    };

    Scalar de_system::dg(Vector y){
        return Dg(y).transpose()*(F(y));
    };

    Scalar de_system::ddw(Vector y){
        Matrix m = DDw(y)*(F(y)*F(y).transpose());
        Scalar num = Dw(y).transpose()*(DF(y).transpose()*F(y));
        return m.trace() + num; 
    };

    Scalar de_system::ddg(Vector y){
        Matrix m = DDg(y)*(F(y)*F(y).transpose());
        Scalar num = Dg(y).transpose()*(DF(y).transpose()*F(y));
        return m.trace() + num; 
    };

//    Scalar de_system::d3w(Vector y){};
//    Scalar de_system::d4w(Vector y){};
//    Scalar de_system::d3g(Vector y){};
//    Scalar de_system::d4g(Vector y){};
}




