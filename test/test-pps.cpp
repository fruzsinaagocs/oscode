#include "test-pps.hpp"
// Settings for the Mukhanov--Sasaki equation.
int n=2;
double m=1e-5;
double k=10000.0;
double phi_p=23.0;
double mp=1.0;

int main(){

    // First solve for the background only with NAG to determine k_pivot.

    // Solving the MS equation with RKWKB. 
    double t = 1.0;
    Vector ic(6);                                      
    ic << 1000.0*k, 0.0, background(t);
    de_system my_system(F, DF, w, Dw, DDw, g, Dg, DDg);   
    for(int i=0; i<100; i++){
    Solution my_solution(my_system, ic, t, f_end, 1, 1e-5, 1e-7, 1.0, "outputs/ms-rkwkb.txt");
    my_solution.evolve();
    //std::cout << "solution at " << my_solution.t << ": " << my_solution.y(0) << "," << my_solution.y(1) << " " << my_solution.y(2) << " " << my_solution.y(3) << " " << my_solution.y(4) << " " << my_solution.y(5) <<  std::endl;
    //std::cout << "total steps: " << my_solution.stepsall << std::endl;
    //std::cout << "waste ratio: " << my_solution.waste << std::endl;
    };
};

};
