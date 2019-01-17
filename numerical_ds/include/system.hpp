#pragma once

class de_system
    {
    private:
                // interpolator
        // interpolation checker

    public:
        // constructor
        // TODO: will need arrays as input later
        de_system(std::complex<double> (*w)(double), std::complex<double> (*g)(double));
        std::complex<double> (*w) (double);
        std::complex<double> (*g) (double);
   
};
    
de_system::de_system(std::complex<double> (*W)(double), std::complex<double> (*G)(double)){
        // default constructor for a system of differential equations
    
        w = W;
        g = G;
//        std::cout << w(1.0) << std::endl;
    };
