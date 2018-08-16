#include "integrators.hpp"

namespace RKWKB{
    class Solution
    {
        private:
    
    
        public:
        // default constructor
        Solution();
        // constructor overloads
        Solution(de_system de_sys, Vector ic, double t_ini, double t_fin, double r_tol=1e-4, double a_tol=0.0, double stepsize=1.0);
    
        // class data 
        double t_i, t_f, rtol, atol, h; // current stepsize
        Vector y; // current solution vector
        Vector error; // current error on solution vector
        de_system de_sys;
        RKFsolver rkfsolver;
        WKBsolver wkbsolver;
    
        // class functions
        void evolve();
        Step step(Method * method, Vector Y, double H);
        
    };
    
    Solution::Solution(){
        // default constructor for a solution object (does nothing)
    }
    
    Solution::Solution(de_system system, Vector ic, double t_ini, double t_fin, double r_tol, double a_tol, double stepsize){
        // constructor for solution of a system of differential equations
        
        y = ic;
        t_i = t_ini;
        t_f = t_fin;
        rtol = r_tol;
        atol = a_tol;
        h = stepsize;

        de_sys = system;
        rkfsolver = RKFsolver(de_sys.F);
        wkbsolver = WKBsolver(de_sys.F);
    };
    
    void Solution::evolve(){
        // function to compute full numerical solution of the de_system.
    };
    
    Step Solution::step(Method * method, Vector Y, double H){
        // function to take a single RKWKB step in the numerical solution of the
        // de_system, from y with stepsize h.
    
        Step s = method->step(de_sys.F, Y, H);
        return s;
    };
}
