#include "integrators.hpp"

class solution
{
    private:


    public:
    // default constructor
    solution();
    // constructor overloads
    solution(de_system, Vector, double, double);

    // class data 
    double t_i, t_f, rtol, atol, h_i;
    Vector y_i;

    // class functions
    void evolve();
    void step();

};

solution::solution(){
    // default constructor for a solution object (does nothing)
}

solution::solution(de_system de_sys, Vector ic, double t_ini, double t_fin){
    // constructor for solution of a system of differential equations
    
    y_i = ic;
    t_i = t_ini;
    t_f = t_fin;
    rtol = 1e-4;
    atol = 0.0;
    h_i = 1.0;
    std::cout << "Constructed a solution object with properties: ic " << y_i << ", ti " << t_i << ", atol " << atol;

};

void solution::evolve(){
    // function to compute full numerical solution of the de_system.
};

void solution::step(){
    // function to take a single step in the numerical solution of the
    // de_system.
};
