#include "integrators.hpp"
#include <fstream>

namespace RKWKB{
    
    class Solution
    {
        private:
    
    
        public:
        // default constructor
        Solution();
        // constructor overloads
        Solution(de_system de_sys, Vector ic, double t_ini, event F_end, int o=1, double r_tol=1e-4, double a_tol=0.0, double stepsize=1.0, std::string output="output.txt");
    
        // class data 
        int d, order; // size of solution vector, order of WKB solution
        double t_i, rtol, h; // current stepsize
        double t; // current time
        Vector y; // current solution vector
        Vector error; // current error on solution vector
        Vector atol;
        de_system de_sys;
        RKFsolver rkfsolver;
        WKBsolver * wkbsolver;
        WKBsolver1 wkbsolver1;
        WKBsolver2 wkbsolver2;
        WKBsolver3 wkbsolver3;
        event f_end; // function s.t. solution terminates at f(y,t)=0
        std::string outputfile;
    
        // class functions
        Vector F_tot(Vector); // gives time-derivative of all variables propagated
        void evolve();
        Step step(Integrator * method, Vectorfn F, Vector Y, double H); 
        void write(std::string output);
        void update_h(Vector err, Vector scale, bool wkb, bool success);
        
    };
    
    Solution::Solution(){
        // default constructor for a solution object (does nothing)
    };
    
    Solution::Solution(de_system system, Vector ic, double t_ini, event F_end, int o, double r_tol, double a_tol, double stepsize, std::string output){
        // constructor for solution of a system of differential equations
        
        order = o;
        d = ic.size() + (order+2)/2;
        y.resize(d);
        y << ic, Vector::Zero((order+2)/2);
        t_i = t_ini;
        rtol = r_tol;
        atol = a_tol*Vector::Ones(y.size());
        h = stepsize;
        t = t_i;
        f_end = F_end;
        error = Vector::Zero(y.size());
        outputfile = output;


        de_sys = system;
        rkfsolver = RKFsolver(de_sys);
        switch(order){
            case 1: wkbsolver1 = WKBsolver1(de_sys,1);
                    wkbsolver = &wkbsolver1;
                    break;
            case 2: wkbsolver2 = WKBsolver2(de_sys,2);
                    wkbsolver = &wkbsolver2;
                    break;
            case 3: wkbsolver3 = WKBsolver3(de_sys,3);
                    wkbsolver = &wkbsolver3;
                    break;
        };
        //write(outputfile);
    };
                
    void Solution::evolve(){
        // function to compute full numerical solution of the de_system.
    
        bool sign = (f_end(y, t_i).real() <= 0.0);
        bool next_sign = sign;
        Vector y_next = y;
        Vector error_next = error;
        Vector scale = rtol*y + atol;
        double t_next = t;
        Step step_rkf, step_wkb;
        double end_error;
        bool wkb=false;
                                
        while(true){
            // keep updating stepsize until step is successful
            while(true){
                step_rkf = step(&rkfsolver, de_sys.F, y, h);
                step_wkb = step(&rkfsolver, de_sys.F, y, h); // for testing, only take RKF steps.
                wkb = step_rkf.error.norm() > step_wkb.error.norm();
                if(wkb){
                    y_next = step_rkf.y;
                    error_next = step_rkf.error;
                }
                else{
                    y_next = step_wkb.y;
                    error_next = step_wkb.error;
                }
                scale = rtol*y_next + atol; 
                if(error_next.norm() <= scale.norm()){
                    t_next += h;   
                    break;
                }
                else{
                    update_h(error_next, scale, wkb, false);
                };
            };

            // check if f_end has changed sign
            next_sign = (f_end(y_next, t_next).real() <= 0.0);
            end_error = std::abs(f_end(y_next, t_next));
            if(next_sign != sign){
                h = (t_next - t)/2.0;
                t_next = t;
            }
            else{
                y = y_next;
                t = t_next;
                error = error_next;
                // update stepsize
                update_h(error_next, scale, wkb, true);
                std::cout << "time: " << t << ", solution: " << y << std::endl;
                write(outputfile);
                if(std::abs(end_error) < 1e-4)
                    break;
            };
        };
    };

    Step Solution::step(Integrator * integrator, Vectorfn F, Vector Y, double H){
        // function to take a single step with a given method in the numerical solution of the
        // de_system, from y with stepsize h.
    
        Step s = integrator->step(F, Y, H);
        return s;
    };
    
    void Solution::update_h(Vector err, Vector scale, bool wkb, bool success){
        // updates stepsize.

        if(success){
            if(wkb)
                h*=1;
            else
                h*=pow(err.norm()/scale.norm(), -1/5.0);
        }
        else{
            if(wkb)
                h*=1;
            else
                h*=pow(err.norm()/scale.norm(), -1/4.0);
        };
    };

    void Solution::write(std::string output){
        // write current step to file

        std::ofstream fout;
        fout.open(output, std::ios_base::app); // appends at the end of outputfile
        fout << t << " ";
        for(int i=0; i<y.size(); i++)
            fout << y(i) << " ";
        for(int i=0; i<error.size(); i++) 
            fout << error(i) << " ";
        fout << std::endl;
        fout.close();
    };

    Vector Solution::F_tot(Vector z){
        // time-derivative of all variables propagated. z = [x, x', y, S_i]
        Vector result(d);
        Vector z_bg = z.segment(2, d-2-(order+2)/2);
        result << z(1), -z(0)*std::pow(de_sys.w(z_bg),2) -2.0*z(1)*de_sys.g(z_bg), de_sys.F(z_bg), wkbsolver->dS_even(z_bg);
        return result;
    };

}






