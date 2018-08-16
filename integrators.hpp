# include "system.hpp"

struct Step
{
    // data structure to store the result of a step and its error.
    Vector y, error;
    // constructor takes in size and resizes y, error
    Step(int d){
        y.resize(d);
        error.resize(d);
    };
};

////////////////////////////////////////////////////////////////////////////

class RKFsolver
{
    private:

    public:
    // default constructor
    RKFsolver(); 
    // constructor overloads
    RKFsolver(de_system);

    // class data
    de_system sys;
    Matrix butcher_a;
    Vector butcher_c, butcher_b4, butcher_b5, butcher_r;
    
    // class functions
    Step step(Vectorfn F, Vector y, double h);

};

///////////////////////////////////////////////////////////////////////////

RKFsolver::RKFsolver(){
    // Default constructor of RKFsolver (does nothing)
};

RKFsolver::RKFsolver(de_system system){
    // Constructor for an RKFsolver from a system of differential equations
    sys = system;
    
    butcher_a.resize(5,5);
    butcher_c.resize(6);
    butcher_b4.resize(6);
    butcher_b5.resize(6);
    butcher_r.resize(6);
    
    butcher_a << 1.0/4.0, 0.0, 0.0, 0.0, 0.0,
    			3.0/32.0, 9.0/32.0, 0.0, 0.0, 0.0,
    			1932.0/2197.0, -7200.0/2197.0, 7296.0/2197.0, 0.0, 0.0,
    			439.0/216.0, -8.0, 3680.0/513.0, -845.0/4104.0, 0.0,
    			-8.0/27.0, 2.0, -3544.0/2565.0, 1859.0/4104.0, -11.0/40.0;
    butcher_c << 0.0, 1.0/4.0, 3.0/8.0, 12.0/13.0, 1.0, 1.0/2.0;
    butcher_b4 << 25.0/216.0 , 0.0, 1408.0/2565.0, 2197.0/4104.0, -1.0/5.0 , 0.0;
    butcher_b5 << 16.0/135.0 , 0.0 , 6656.0/12825.0 , 28561.0/56430.0, -9.0/50.0, 2.0/55.0;
    butcher_r << 1.0/360.0, 0.0, -128/4275.0, -2197/75240.0, 1/50.0, 2/55.0; 

    std::cout << "Constructed an RKFsolver object" << std::endl;
};

Step RKFsolver::step(Vectorfn F, Vector y, double h){
    // Stepper function using the RKF method
    
    int d = y.size();
    Step result(d);
    std::cout << "Took a step!" << std::endl;
    return result;
};

///////////////////////////////////////////////////////////////////////////

class WKBsolver
{
    private:

    public:
    // default constructor
    WKBsolver();    
    // constructor overloads
    WKBsolver(de_system);

    // class data
    de_system sys;
    
    // class methods
    void step(); 

};

WKBsolver::WKBsolver(){
    // Default constructor for a WKBsolver (does nothing)
};

WKBsolver::WKBsolver(de_system system){
    // Constructor for a WKBsolver from a system of differential equations
    sys = system;
    std::cout << "Constructed a WKBsolver object" << std::endl;
};

void WKBsolver::step(){
    // Stepper function using the RKF method   
};

