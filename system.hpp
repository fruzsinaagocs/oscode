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
            de_system(Vectorfn);//, Mfn, Sfn, Vfn, Mfn, Sfn, Vfn, Mfn);
            
            // class data
            
            // class functions
            Vectorfn F;
            Matrixfn DF;
            Scalarfn w;
            Vectorfn Dw;
            Matrixfn DDw;
            Scalarfn g;
            Vectorfn Dg;
            Matrixfn DDg;
    
    };
    
    de_system::de_system(){
        // Default constructor of de_system (does nothing)
    }
    
    de_system::de_system(Vectorfn f){//, Mfn Df, Sfn freq, Vfn Dfreq, Mfn DDfreq, Sfn gam, Vfn Dgam, Mfn DDgam){
        // default constructor for a system of differential equations
    
        F = f;
        //DF = Df;
        //w = freq;
        //Dw = Dfreq;
        //DDw = DDfreq;
        //g = gam;
        //Dg = Dgam;
        //DDg = DDgam;
        std::cout << "Constructed a de_system object" << std::endl;
    }
}




